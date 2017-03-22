/* Responsive radar chart to visualize MS compounds table. Adapted from
http://bl.ocks.org/nbremer/6506614. */

RadarChart = function(_parentElement, _data, _eventHandler){
    
    this.parentElement = _parentElement;
    this.data = _data;
    this.displayData = [];
    this.eventHandler = _eventHandler;
    this.margin = {top: 75, right: 20, bottom: 10, left: 100},

    // boot up the viz
    this.initVis();

}

/* Method that sets up the SVG and the variables */
RadarChart.prototype.initVis = function(){
    
    var that = this;

    this.height = this.width = window.innerHeight/2.5;
    this.cradius = 5;
    this.factor = 1;
    this.factorLegend = .85;
    this.levels = 6;
    this.maxValue = 0.6;
    this.radians = 2 * Math.PI;
    this.opacityArea = 0.5;
    this.ToRight = 5;
    this.color = d3.scale.category10();

    // add plotting space
    this.svg = this.parentElement.append("svg")
      .attr("width", this.width + this.margin.left + this.margin.right)
      .attr("height", this.height + this.margin.top + this.margin.bottom)
      .append("g")
      .attr("transform", "translate(" + this.margin.left + "," + this.margin.top + ")");

    // filter, aggregate, modify data
    this.wrangleData();

    // call the update method
    this.updateVis();
}

/* Wrassle the data.*/
RadarChart.prototype.wrangleData = function(){

    var that = this;

    this.displayData = [];
    var uniq_cmpds = _.uniq(this.data.map(function(d){return d.cmpd})).sort()

    uniq_cmpds.forEach(function(d){

        //filter the data for this compound
        var temp_data = that.data.filter(function(k){
            if(k.cmpd == d){return true}
                else {return false}
        })

        // generate an object for each sample
        temp_arr = []
        temp_data.forEach(function(l){
            temp_arr.push({axis:l.var,value:parseInt(l.value)})
        })

        that.displayData.push(temp_arr)
    })

}

/** the drawing function - should use the D3 selection, enter, exit*/
RadarChart.prototype.updateVis = function(){

    var that = this;

    this.maxValue = Math.max(this.maxValue, d3.max(this.displayData, function(i){return d3.max(i.map(function(o){return o.value;}))}));
    this.allAxis = (this.displayData[0].map(function(i, j){return i.axis}));
    this.total = this.allAxis.length;
    this.radius = this.factor*Math.min(this.width/2, this.height/2);
    this.Format = d3.format('%');
    this.tooltip;

    //Circular segments
    for(var j=0; j<this.levels-1; j++){
      var levelFactor = this.factor*this.radius*((j+1)/this.levels);
      this.svg.selectAll(".levels")
       .data(this.allAxis)
       .enter()
       .append("svg:line")
       .attr("x1", function(d, i){return levelFactor*(1-that.factor*Math.sin(i*that.radians/that.total));})
       .attr("y1", function(d, i){return levelFactor*(1-that.factor*Math.cos(i*that.radians/that.total));})
       .attr("x2", function(d, i){return levelFactor*(1-that.factor*Math.sin((i+1)*that.radians/that.total));})
       .attr("y2", function(d, i){return levelFactor*(1-that.factor*Math.cos((i+1)*that.radians/that.total));})
       .attr("class", "line")
       .style("stroke", "grey")
       .style("stroke-opacity", "0.75")
       .style("stroke-width", "0.3px")
       .attr("transform", "translate(" + (that.width/2-levelFactor) + ", " + (that.height/2-levelFactor) + ")");
    }

    //Text indicating at what % each level is
    for(var j=0; j<this.levels; j++){
      var levelFactor = this.factor*this.radius*((j+1)/this.levels);
      this.svg.selectAll(".levels")
       .data([1]) //dummy data
       .enter()
       .append("svg:text")
       .attr("x", function(d){return levelFactor*(1-that.factor*Math.sin(0));})
       .attr("y", function(d){return levelFactor*(1-that.factor*Math.cos(0));})
       .attr("class", "legend")
       .style("font-family", "sans-serif")
       .style("font-size", "10px")
       .attr("transform", "translate(" + (that.width/2-levelFactor + that.ToRight) + ", " + (that.height/2-levelFactor) + ")")
       .attr("fill", "#737373")
       .text(that.Format((j+1)*that.maxValue/that.levels));
    }

    var series = 0;

    var axis = this.svg.selectAll(".axis")
            .data(this.allAxis)
            .enter()
            .append("g")
            .attr("class", "axis");

    axis.append("line")
        .attr("x1", that.width/2)
        .attr("y1", that.height/2)
        .attr("x2", function(d, i){return that.width/2*(1-that.factor*Math.sin(i*that.radians/that.total));})
        .attr("y2", function(d, i){return that.height/2*(1-that.factor*Math.cos(i*that.radians/that.total));})
        .attr("class", "line")
        .style("stroke", "grey")
        .style("stroke-width", "1px");

    axis.append("text")
        .attr("class", "legend")
        .text(function(d){return d})
        .style("font-family", "sans-serif")
        .style("font-size", "11px")
        .attr("text-anchor", "middle")
        .attr("dy", "1.5em")
        .attr("transform", function(d, i){return "translate(0, -10)"})
        .attr("x", function(d, i){return that.width/2*(1-that.factorLegend*Math.sin(i*that.radians/that.total))-60*Math.sin(i*that.radians/that.total);})
        .attr("y", function(d, i){return that.height/2*(1-Math.cos(i*that.radians/that.total))-20*Math.cos(i*that.radians/that.total);});

    this.displayData.forEach(function(y, x){
      
      var dataValues = [];
      
      that.svg.selectAll(".nodes")
        .data(y, function(j, i){
          dataValues.push([
            that.width/2*(1-(parseFloat(Math.max(j.value, 0))/that.maxValue)*that.factor*Math.sin(i*that.radians/that.total)), 
            that.height/2*(1-(parseFloat(Math.max(j.value, 0))/that.maxValue)*that.factor*Math.cos(i*that.radians/that.total))
          ]);
        });

      dataValues.push(dataValues[0]);
      
      that.svg.selectAll(".area")
                     .data([dataValues])
                     .enter()
                     .append("polygon")
                     .attr("class", "radar-chart-serie"+series)
                     .style("stroke-width", "2px")
                     .style("stroke", that.color(series))
                     .attr("points",function(d) {
                         var str="";
                         for(var pti=0;pti<d.length;pti++){
                             str=str+d[pti][0]+","+d[pti][1]+" ";
                         }
                         return str;
                      })
                     .style("fill", function(j, i){return that.color(series)})
                     .style("fill-opacity", that.opacityArea)
                     .on('mouseover', function (d){
                                        z = "polygon."+d3.select(this).attr("class");
                                        that.svg.selectAll("polygon")
                                         .transition(200)
                                         .style("fill-opacity", 0.1); 
                                        that.svg.selectAll(z)
                                         .transition(200)
                                         .style("fill-opacity", .7);
                                      })
                     .on('mouseout', function(){
                                        that.svg.selectAll("polygon")
                                         .transition(200)
                                         .style("fill-opacity", that.opacityArea);
                     });
      series++;
    });
    
    series=0;

    this.displayData.forEach(function(y, x){
      dataValues = [];
      that.svg.selectAll(".nodes")
        .data(y).enter()
        .append("svg:circle")
        .attr("class", "radar-chart-serie"+series)
        .attr('r', that.cradius)
        .attr("alt", function(j){return Math.max(j.value, 0)})
        .attr("cx", function(j, i){
          dataValues.push([
            that.width/2*(1-(parseFloat(Math.max(j.value, 0))/that.maxValue)*that.factor*Math.sin(i*that.radians/that.total)), 
            that.height/2*(1-(parseFloat(Math.max(j.value, 0))/that.maxValue)*that.factor*Math.cos(i*that.radians/that.total))
        ]);
        return that.width/2*(1-(Math.max(j.value, 0)/that.maxValue)*that.factor*Math.sin(i*that.radians/that.total));
        })
        .attr("cy", function(j, i){
          return that.height/2*(1-(Math.max(j.value, 0)/that.maxValue)*that.factor*Math.cos(i*that.radians/that.total));
        })
        .attr("data-id", function(j){return j.axis})
        .style("fill", that.color(series)).style("fill-opacity", .9)
        /*.on('mouseover', function (d){
                    newX =  parseFloat(d3.select(this).attr('cx')) - 10;
                    newY =  parseFloat(d3.select(this).attr('cy')) - 5;
                    
                    that.tooltip
                        .attr('x', newX)
                        .attr('y', newY)
                        .text(that.Format(d.value))
                        .transition(200)
                        .style('opacity', 1);
                        
                    z = "polygon."+d3.select(this).attr("class");
                    that.svg.selectAll("polygon")
                        .transition(200)
                        .style("fill-opacity", 0.1); 
                    that.svg.selectAll(z)
                        .transition(200)
                        .style("fill-opacity", .7);
                  })
        .on('mouseout', function(){
                    that.tooltip
                        .transition(200)
                        .style('opacity', 0);
                    that.svg.selectAll("polygon")
                        .transition(200)
                        .style("fill-opacity", that.opacityArea);
                  })
        .append("svg:title")
        .text(function(j){return Math.max(j.value, 0)}); */

      series++;
    });

}

/* Define behavior on user input.*/
RadarChart.prototype.onSelectionChange= function(pass){

    // unpack passed object
    type = pass["type"]

    // call relevant function

}