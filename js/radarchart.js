/* Responsive radar chart to visualize MS compounds table. Adapted from
http://bl.ocks.org/nbremer/6506614. */

RadarChart = function(_parentElement, _data, _eventHandler){

    this.parentElement = _parentElement;
    this.data = _data;
    this.displayData = [];
    this.eventHandler = _eventHandler;
    this.margin = {top: 50, right: 50, bottom: 40, left: 200},

    // boot up the viz
    this.initVis();

}

/* Method that sets up the SVG and the variables */
RadarChart.prototype.initVis = function(){

    var that = this;

    this.height = window.innerHeight/4;
    this.width = window.innerWidth/2;
    this.cradius = 5;
    this.factor = 1;
    this.factorLegend = .85;
    this.levels = 6;
    this.maxValue = 1.25;
    this.radians = 2 * Math.PI;
    this.opacityArea = 0.5;
    this.ToRight = 5;
    this.color = d3.scale.category10();
    this.mode;
    // keep track of computed datasets
    this.quantdata = null;
    this.bindata = null;
    this.all_vals={quant:[], bin:[]};

    // add plotting space
    this.svg = this.parentElement.append("svg")
      .attr("width", this.width + this.margin.left + this.margin.right)
      .attr("height", this.height + this.margin.top + this.margin.bottom)
      .append("g")
      .attr("transform", "translate(" + this.margin.left + "," + this.margin.top + ")");

    // filter, aggregate, modify data
    this.wrangleData("quant");

    // call the update method
    this.updateVis();
}

/* Wrassle the data.*/
RadarChart.prototype.wrangleData = function(_param){

    var that = this;
    // update mode global
    this.mode = _param;
    this.tempData = [];
    // get unique compounds
    var uniq_cmpds = _.uniq(this.data.map(function(d){return d.cmpd})).sort()

    // only compute if not previously computed
    if (this.quantdata == null | this.bindata == null){
        uniq_cmpds.forEach(function(d){

            //filter the data for this compound
            var temp_data = that.data.filter(function(k){
                if(k.cmpd == d){return true}
                    else {return false}
            })

            // generate a data object for each sample
            temp_arr = []
            temp_data.forEach(function(l){

                // if computing as binary 1/0
                if (that.mode == "bin"){
                    temp = (parseFloat(l.log_value) > 0); disp_value = temp ? 1 : 0;
                }
                    // if quantitative
                    else {disp_value = parseFloat(l.log_value)}
                // add as data object
                temp_arr.push({axis:that.cleanSampleName(l.var),value: disp_value,cmpd:l.cmpd});
                // get flat array of vals for later use
                that.all_vals[_param].push(disp_value)
            })

            // add for each compound
            that.tempData.push(temp_arr)
        })

        // update data globals
        if (_param=="quant"){this.displayData = this.quantdata = this.tempData}
            else {this.displayData = this.bindata = this.tempData}
    }
    // if previously computed
    else {

        // update data globals
        if (_param == "quant"){this.displayData = this.quantdata}
            else {this.displayData = this.bindata}
    }

    this.updateVis();
}

/** the drawing function - should use the D3 selection, enter, exit*/
RadarChart.prototype.updateVis = function(){

    var that = this;

    // compute data-dependant variables
    this.maxValue = d3.max([d3.max(this.all_vals[this.mode]), 1.25])
    this.allAxis = (this.displayData[0].map(function(i, j){return i.axis}));
    this.total = this.allAxis.length;
    this.radius = this.factor*Math.min(this.height/2, this.height/2);
    this.Format = d3.format('.1f');

    // remove old elements
    // could implement select-enter-exit
    this.svg.selectAll("line").remove()
    this.svg.selectAll(".legend").remove()
    this.svg.selectAll(".label").remove()
    this.svg.selectAll("polygon").remove()
    this.svg.selectAll("circle").remove()

    //Circular segments
    for(var j=0; j<this.levels-1; j++){

      var levelFactor = this.factor*this.radius*((j+1)/this.levels);

      var seg = this.svg.selectAll(".levels")
        .data(this.allAxis)

      seg.enter()
        .append("svg:line")
        .attr("x1", function(d, i){return levelFactor*(1-that.factor*Math.sin(i*that.radians/that.total));})
        .attr("y1", function(d, i){return levelFactor*(1-that.factor*Math.cos(i*that.radians/that.total));})
        .attr("x2", function(d, i){return levelFactor*(1-that.factor*Math.sin((i+1)*that.radians/that.total));})
        .attr("y2", function(d, i){return levelFactor*(1-that.factor*Math.cos((i+1)*that.radians/that.total));})
        .attr("class", "line")
        .style("stroke", "grey")
        .style("stroke-opacity", "0.75")
        .style("stroke-width", "0.3px")
        .attr("transform", "translate(" + (that.height/2-levelFactor) + ", " + (that.height/2-levelFactor) + ")")
        .style("opacity", function(){
          if (that.mode=="quant"){return 1 }
              else {b=(((j+1)*that.maxValue/that.levels)==1); return b ? 1 : 0}
        })
    }

    //Text indicating at what % each level is
    for(var j=0; j<this.levels; j++){

      var levelFactor = this.factor*this.radius*((j+1)/this.levels);

      var labs = this.svg.selectAll(".levels")
        .data([1]) //dummy data

      labs.enter()
        .append("svg:text")
        .attr("x", function(d){return levelFactor*(1-that.factor*Math.sin(0));})
        .attr("y", function(d){return levelFactor*(1-that.factor*Math.cos(0));})
        .attr("class", "legend")
        .style("font-family", "sans-serif")
        .style("font-size", "8px")
        .attr("transform", "translate(" + (that.height/2-levelFactor + that.ToRight) + ", " + (that.height/2-levelFactor) + ")")
        .attr("fill", "#737373")
        .text(function(){return that.Format((j+1)*that.maxValue/that.levels)})
        .style("opacity", function(){
          if (that.mode=="quant"){return 1 }
              else {b=(((j+1)*that.maxValue/that.levels)==1); return b ? 1 : 0}
        })

    }

    // add axes
    var series = 0;
    var axis = this.svg.selectAll(".axis")
      .data(this.allAxis);

    axis.enter()
      .append("g")
      .attr("class", "axis");

    axis.append("line")
      .attr("x1", that.height/2)
      .attr("y1", that.height/2)
      .attr("x2", function(d, i){return that.height/2*(1-that.factor*Math.sin(i*that.radians/that.total));})
      .attr("y2", function(d, i){return that.height/2*(1-that.factor*Math.cos(i*that.radians/that.total));})
      .attr("class", "line")
      .style("stroke", "grey")
      .style("stroke-width", "1px");

    axis.append("text")
      .attr("class", function(d){return "label label" + d})
      .text(function(d){return d})
      .style("font-family", "sans-serif")
      .style("font-size", "10px")
      .attr("text-anchor", "middle")
      .attr("dy", "1.5em")
      .attr("transform", function(d, i){return "translate(0, -10)"})
      .attr("x", function(d, i){return that.height/2*(1-that.factorLegend*Math.sin(i*that.radians/that.total))-60*Math.sin(i*that.radians/that.total);})
      .attr("y", function(d, i){return that.height/2*(1-Math.cos(i*that.radians/that.total))-20*Math.cos(i*that.radians/that.total);});

    // create polygons for each compound
    this.displayData.forEach(function(y, x){

      var dataValues = [];

      // calculate nodes for polygon
      var nodes = that.svg.selectAll(".nodes")
        .data(y, function(j, i){
          dataValues.push([
            that.height/2*(1-(parseFloat(Math.max(j.value, 0))/that.maxValue)*that.factor*Math.sin(i*that.radians/that.total)),
            that.height/2*(1-(parseFloat(Math.max(j.value, 0))/that.maxValue)*that.factor*Math.cos(i*that.radians/that.total))
          ]);
        });

      // add in first one twice so polygon connects back
      dataValues.push(dataValues[0]);

      var polys = that.svg.selectAll(".area")
        .data([dataValues])

      polys.enter()
        .append("polygon")
        .attr("class", function(d){return "cmpd" + y[0]["cmpd"].replace(/#/g , "")})
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
        // initialize as hidden
        .style("fill-opacity", 0)
        .style("stroke-opacity", 0)

      series++;

    });

    series=0;

    // add in data points for each compound
    this.displayData.forEach(function(y, x){

      dataValues = [];

      var circs = that.svg.selectAll(".nodes")
        .data(y)

      circs.enter()
        .append("svg:circle")
        .attr("class", function(d){return "cmpd" + y[0]["cmpd"].replace(/#/g , "")})
        .attr('r', that.cradius)
        .attr("alt", function(j){return Math.max(j.value, 0)})
        .attr("cx", function(j, i){
          dataValues.push([
            that.height/2*(1-(parseFloat(Math.max(j.value, 0))/that.maxValue)*that.factor*Math.sin(i*that.radians/that.total)),
            that.height/2*(1-(parseFloat(Math.max(j.value, 0))/that.maxValue)*that.factor*Math.cos(i*that.radians/that.total))
          ]);
          return that.height/2*(1-(Math.max(j.value, 0)/that.maxValue)*that.factor*Math.sin(i*that.radians/that.total));
        })
        .attr("cy", function(j, i){
          return that.height/2*(1-(Math.max(j.value, 0)/that.maxValue)*that.factor*Math.cos(i*that.radians/that.total));
        })
        .attr("data-id", function(j){return j.axis})
        .style("fill", that.color(series))
        // initialize as hidden
        .style("fill-opacity", 0)

      series++;

    });
}

/* Define behavior on user input.*/
RadarChart.prototype.oncellMouse= function(pass){

    // unpack passed object
    cmpd = pass.cmpd
    cmpd_class = ".cmpd" + cmpd.replace(/#/g, "")
    // construct and execute class call
    // make visible on mouse over
    if (pass.dtype=="on"){
        this.svg.selectAll(cmpd_class).style("stroke-opacity",1).style("fill-opacity", 0.7)
    }
    // hide again on mouse off
    else if (pass.dtype=="off"){
        this.svg.selectAll(cmpd_class).style("stroke-opacity",0).style("fill-opacity", 0)
    }
}

/* Define behavior on user input.*/
RadarChart.prototype.onSelectionChange= function(pass){

    var that = this;

    // only interested in data type changes
    if (pass["type"] == "data_type"){

      // change some params
      if (pass.value=="bin"){
          that.levels = 5;
      }
      else {that.levels=6}

      // re-construct the data
      this.wrangleData(pass.value);
    }
}

/** cleanse sample names to avoid invalid class selectors */
RadarChart.prototype.cleanSampleName = function(name){

    // starts with numeric?
    if (!isNaN(name[0] - parseFloat(name[0]))){
      temp = "_" + name
    }
    else {temp = name}
    // remove invalid chars
    return temp.replace(/[~!@$%^&*()+=,./';:"?><\[\]\{\}|`#]/g, "_")

}
