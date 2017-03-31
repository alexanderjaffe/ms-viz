
CountVis = function(_parentElement, _data, _eventHandler){
    
    this.parentElement = _parentElement;
    this.data = _data;
    this.eventHandler = _eventHandler;
    this.displayData = [];
    this.metaData;

    this.margin = {top: 50, right: 20, bottom:20 , left: 100},
    this.width = window.innerWidth/2 - this.margin.left - this.margin.right,
    this.height = window.innerHeight/4; - this.margin.top - this.margin.bottom;

    this.initVis();
}

CountVis.prototype.initVis = function(){

    // get mode and save out
    themode = this.mode_id;
    var that = this;

    // from here down, adapted from CS171 section 6 
    this.svg = this.parentElement.append("svg")
        .attr("width", this.width + this.margin.left + this.margin.right)
        .attr("height", this.height + this.margin.top + this.margin.bottom)
        .append("g")
        .attr("transform", "translate(" + this.margin.left + "," + this.margin.top + ")");

    // creates axis and scales
    this.x = d3.scale.linear()
      .range([0, this.width]);

    this.y = d3.scale.linear()
      .range([this.height, 0]);

    this.xAxis = d3.svg.axis()
      .scale(this.x)
      .orient("bottom");

    this.yAxis = d3.svg.axis()
      .scale(this.y)
      .orient("left");

    this.format = d3.format(".1f")

    // Add axes visual elements
    this.svg.append("g")
        .attr("class", "x axis2")
        .attr("transform", "translate(0," + this.height + ")")
        .append("text")
        .attr("x", that.width)
        .attr("y", -10)
        .attr("dy", ".71em")
        .style("text-anchor", "end")
        .text("m/z");

    this.svg.append("g")
        .attr("class", "y axis2")
        .append("text")
        .attr("transform", "rotate(-90)")
        .attr("y", 6)
        .attr("dy", ".71em")
        .style("text-anchor", "end")
        .text("Raw Intensity");

    // filter, aggregate, modify data
    this.wrangleData(null);

    // call the up method
    this.updateVis();

}

CountVis.prototype.wrangleData= function(pass){

    // filter for passed cmpd and sample
    if (pass){

        cmpd = pass.cmpd.replace("#", "compound_")
        sample = pass.sample

        this.intData = this.data.filter(function(d){

            if (d.compound == cmpd && d.sample == sample){
                return true
            }
            else {return false}
        })

        this.displayData = this.intData[0]["spectrum"]
    }

    else {this.displayData = this.data[20]["spectrum"]}
    // else display data stays empty - no spectrum
}

CountVis.prototype.updateVis = function(){

    var that = this;

    // set scale domains using data
    this.xmax = d3.max(this.displayData, function(d){return d.mz})
    this.xmin = d3.min(this.displayData, function(d){return d.mz})
    this.ymax = d3.max(this.displayData, function(d){return d.i})
    this.ymin = d3.min(this.displayData, function(d){return d.i})
    // add padding to make room for text
    this.x.domain([(this.xmin - this.xmax/10),(this.xmax+this.xmax/10)])
    this.y.domain([this.ymin,(this.ymax+this.ymax/8)])

    // calculate 50th percentile value
    this.p = this.percentile(this.displayData.map(function(d){return d.i}),0.90)

    // updates axis
    this.svg.select(".x.axis2")
        .transition().duration(750)
        .call(this.xAxis);

    this.svg.select(".y.axis2")
        .transition().duration(750)
        .call(this.yAxis)

    // bind new data 
    var peaks = this.svg.selectAll(".peaks")
        .data(this.displayData)

    // EXIT - remove extra elements
    peaks.exit().remove();

    // update the line
    peaks.select("line")
        .transition().duration(750)
        .attr("x1", function(d){return that.x(d.mz)})
        .attr("y2", function(d){return that.y(that.ymin)})
        .attr("x2", function(d){return that.x(d.mz)})
        .attr("y1", function(d){return that.y(d.i)})

    // then update the labels
    peaks.select("text")
        .transition().duration(750)
        .attr("x", function(d){return that.x(d.mz)})
        .attr("y", function(d){return that.y(d.i)})
        //.attr("class", function(d,i){return "label" + i})
        .text(function(d){return that.format(d.mz)})
        .style("visibility", function(d){
            if (d.i >= that.p){return "visible"}
                else {return "hidden"}
        })

    // ENTER - add any new groups and set class attribute 
    var peaks = peaks.enter().append("g").attr("class", "peaks")

    // add elements to these groups
    peaks.append("svg:line")
    peaks.append("text")

    // set main attributes for lines
    peaks.select("line")
        .attr("x1", function(d){return that.x(d.mz)})
        .attr("y2", function(d){return that.y(that.ymin)})
        .attr("x2", function(d){return that.x(d.mz)})
        .attr("y1", function(d){return that.y(d.i)})
        .style("stroke","black")
        .on("mouseover", function(d,i){

            d3.select(".label" + i).style("fill", "red")
        })
        .on("mouseout", function(d,i){

            d3.select(".label" + i).style("fill", "black")
        })

    // set main attributes for labels
    peaks.select("text")
        .attr("class", "peaks-label")
        .attr("x", function(d){return that.x(d.mz)})
        .attr("y", function(d){return that.y(d.i)})
        .attr("dy", "-.4em")
        .attr("class", function(d,i){return "label" + i})
        .text(function(d){return that.format(d.mz)})
        .style("visibility", function(d){
            if (d.i >= that.p){return "visible"}
                else {return "hidden"}
        })
        .style("font-size", "10px")
        .style("text-anchor", "middle")
}

CountVis.prototype.oncellMouse= function (pass){

    this.wrangleData(pass);

    this.updateVis();
}

// from https://gist.github.com/IceCreamYou/6ffa1b18c4c8f6aeaad2
CountVis.prototype.percentile = function (arr, p) {
    
    if (arr.length === 0) return 0;
    if (typeof p !== 'number') throw new TypeError('p must be a number');
    if (p <= 0) return arr[0];
    if (p >= 1) return arr[arr.length - 1];

    var index = arr.length * p,
        lower = Math.floor(index),
        upper = lower + 1,
        weight = index % 1;

    if (upper >= arr.length) return arr[lower];
    
    return arr[lower] * (1 - weight) + arr[upper] * weight;
}
