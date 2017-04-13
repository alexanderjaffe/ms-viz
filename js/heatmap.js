/* Interactive heatmap to visualize MS compounds table. Adapted from
http://bl.ocks.org/ianyfchang/8119685. */

HeatMap = function(_parentElement, _data, _eventHandler){
    
    this.parentElement = _parentElement;
    this.data = _data;
    this.displayData = [];
    this.eventHandler = _eventHandler;
    this.margin = {top: 50, right: 20, bottom: 10, left: 100},

    // boot up the viz
    this.initVis();

}

/* Method that sets up the SVG and the variables */
HeatMap.prototype.initVis = function(){
    
    var that = this;

    this.width = window.innerWidth - this.margin.right - 120 ;
    //this.width = this.cellSize*this.col_number, // - margin.left - margin.right,
    this.height = window.innerHeight/2.5 // - margin.top - margin.bottom,
    //gridSize = Math.floor(width / 24),
    this.colorBuckets = 21,
    this.rowSortOrder=false;
    this.colSortOrder=false;
    this.colorScale = d3.scale.linear()
    // initalize data mode
    this.data_mode = "quant";
    // initialize GNPS hits
    this.hits = false
    
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
HeatMap.prototype.wrangleData = function(){

    // displayData should hold the data which is visualized
   //this.displayData = this.data.map(function(d) {
    this.intData = this.data.map(function(d) {
        return {

            row:   +parseInt(d.row),
            col:   +parseInt(d.col),
            value: +parseFloat(d.log_value),
            hrow: +parseInt(d.hrow),
            hcol: +parseInt(d.hcol),
            ccol: +parseInt(d.ccol),
            sample: d.var,
            cmpd: d.cmpd,
            mass: d.Mass,
            hit: d.hit,
            score: d.score
        };
    })

    /*this.displayData = this.intData.filter(function(d){
        if (d.col > 400){return false}
            else {return true}
    })*/
    
    this.displayData = this.intData;

}

/** the drawing function - should use the D3 selection, enter, exit*/
HeatMap.prototype.updateVis = function(){

    var that = this;

    // calculate min and max vals in data
    this.dmax = d3.max(this.displayData.map(function(d){return d.value}))
    // only search vals greater than zero
    this.dmin = d3.min(this.displayData.filter(function(d){
        if (d.value==0){return false}
            else {return true}
        })
        .map(function(d){return d.value}))
    // define color scale params using displaydata
    this.colorScale.domain([this.dmin,this.dmax])
      .range(["white","blue"]);

    // calculate min and max table indices
    this.col_number = d3.max(this.displayData.map(function(d){return d.col}))
    this.row_number = d3.max(this.displayData.map(function(d){return d.row}))
    // use this to determine cell width/height
    this.cellWidth = this.width/this.col_number;
    this.cellSize = this.height/this.row_number;
    // this.legendElementWidth = this.cellSize*2.5;

    // get row and col labels
    this.rowLabel = _.uniq(this.displayData.map(function(d){return d.sample})).sort()
    this.colLabel = _.uniq(this.displayData.map(function(d){return d.cmpd})).sort()
    // get position info for rows/cols
    this.rowinfo = _.uniq(this.displayData, function(d){return d.sample})
    this.colinfo = _.uniq(this.displayData, function(d){return d.cmpd})

    // col label tweak factor
    this.col_adjust = 2;

    // add in rowlabels
    var rowLabels = this.svg.append("g")
        .selectAll(".rowLabelg")
        .data(this.rowLabel)
        .enter()
        .append("text")
        .text(function (d) { return d.slice(0,15); })
        .attr("x", 0)
        .attr("y", function(d,i){return i*that.cellSize;})
        .style("text-anchor", "end")
        .style("font-size", "10px")
        .attr("transform", "translate(-6," + that.cellSize/1.5 + ")")
        .attr("class", function (d,i) { return "rowLabel mono r"+i;} ) 
        .on("mouseover", function(d) {d3.select(this).classed("text-hover",true);})
        .on("mouseout" , function(d) {d3.select(this).classed("text-hover",false);})
        .on("click", function(d,i) {that.rowSortOrder=!that.rowSortOrder; that.sortbylabel("r",i,that.rowSortOrder);d3.select("#order").property("selectedIndex", 4).node().focus();;});

    // add in collabels
    var colLabels = this.svg.append("g")
        .selectAll(".colLabelg")
        .data(this.colLabel)
        .enter()
        .append("text")
        .text(function (d) { return d; })
        .attr("x", 0)
        .attr("y", function(d,i){return i*that.cellWidth})
        .style("text-anchor", "left")
        .style("opacity","0")
        .attr("transform", "translate("+ that.cellWidth/1.5 + ",-6) rotate (-90)")
        .attr("class",  function (d,i) { return "colLabel mono c"+i;} )
        .on("mouseover", function(d) {d3.select(this).classed("text-hover",true);})
        .on("mouseout" , function(d) {d3.select(this).classed("text-hover",false);})
        //.on("click", function(d,i) {console.log(i); that.colSortOrder=!that.colSortOrder; that.sortbylabel("c",i,that.colSortOrder);d3.select("#order").property("selectedIndex", 4).node().focus();;});

    // main heatmap
    var heatMap = this.svg.append("g").attr("class","g3")
        .selectAll(".cellg")
        .data(this.displayData,function(d){return d.row+":"+d.col;})
        .enter()
        .append("rect")
        .attr("x", function(d) { return (d.col-1) * that.cellWidth; })
        .attr("y", function(d) { return (d.row-1) * that.cellSize; })
        .attr("class", function(d){return "cell cell-border cr"+(d.row-1)+" cc"+(d.col-1);})
        .attr("width", that.cellWidth)
        .attr("height", that.cellSize)
        .style("fill", function(d) { return that.getColor(d.value); })
        .on("click", function(d,i) {that.colSortOrder=!that.colSortOrder; that.sortbylabel("c",(d.col-1),that.colSortOrder);d3.select("#order").property("selectedIndex", 4).node().focus();;})
        .on("mouseover", function(d){
            
            //highlight text
            d3.select(this).classed("cell-hover",true);
            d3.selectAll(".rowLabel").classed("text-highlight",function(r,ri){return ri==(d.row-1)});
            // make column visible when hovered
            d3.selectAll(".colLabel")
                .style("opacity",function(c,ci){
                    if (ci==(d.col-1)) {return 1} 
                    else {return 0}})
                .classed("text-highlight",function(c,ci){return ci==(d.col-1);});

            // highlight radar label
            d3.select(".label" + d.sample).style("fill","red")
            // trigger event for spectrum viz
            $(that.eventHandler).trigger("cellMouseover", {cmpd:d.cmpd, mass:d.mass, sample:d.sample, dtype:"on", value: (d.value > 0), hit: d.hit, score: d.score})
            
        })
        .on("mouseout", function(d){
            d3.select(this).classed("cell-hover",false);
            d3.selectAll(".rowLabel").classed("text-highlight",false);
            //d3.selectAll(".colLabel").classed("text-highlight",false);
            d3.selectAll(".colLabel").style("opacity", 0);
            d3.select("#tooltip").classed("hidden", true);

            // unhighlight radar label
            d3.select(".label" + d.sample).style("fill","black")

            // trigger event for spectrum viz
            $(that.eventHandler).trigger("cellMouseover", {cmpd:d.cmpd, mass:d.mass, sample:d.sample, dtype:"off"})
        });

    // add in legend

    var lh = 50;
    var lw = 400;

    var key = d3.select("#legend").append("svg").attr("width", lw).attr("height", lh);

    var legend = key.append("defs").append("svg:linearGradient")
        .attr("id", "gradient").attr("x1", "0%").attr("y1", "100%")
        .attr("x2", "100%").attr("y2", "100%")
        //.attr("spreadMethod", "pad");

    legend.append("stop").attr("offset", "0%")
        .attr("stop-color", "white")
        .attr("stop-opacity", 1);

    legend.append("stop").attr("offset", "100%")
        .attr("class","gradienthi")
        .attr("stop-color", "blue")
        .attr("stop-opacity", 1);

    key.append("rect").attr("width", lw-100)
        .attr("height", 10)
        .style("fill", "url(#gradient)")
        .attr("transform", "translate(10,0)");

    this.y = d3.scale.linear().range([lw-100, 0]).domain([this.dmax,this.dmin])

    this.legend_axis = d3.svg.axis().scale(that.y).orient("bottom");

    key.append("g").attr("class", "axis3")
        .attr("transform", "translate(10,12)")
        .call(that.legend_axis).append("text")
        .attr("y", 30).attr("dy", ".71em")
        .style("text-anchor", "end")
        .text("log intensity")
        .attr("class","legendtext")
        .attr("transform", "translate("+(lw-30)+"," + (-22)+")")

}

/* Define behavior on user input.*/
HeatMap.prototype.onSelectionChange= function(pass){

    // unpack passed object
    type = pass["type"]

    // call relevant function
    if (type=="order"){this.order(pass["value"])}
    else if (type == "data_type"){this.heatmapUpdate(pass["value"])}

}

HeatMap.prototype.getColor = function(x){
  if (x == 0){
    return "white";
  }else{
    return this.colorScale(x)
  }
}

// method to return property of a row or column
HeatMap.prototype.getProp = function(array, field, val, out){

    var temp = array.filter(function(d){
        if (d[field] == val){
            return true
        }
        else {return false}
    })

    return temp[0][out]
}

// update heatmap for quantitative/binary modes
HeatMap.prototype.heatmapUpdate = function(value){

    var that = this;
    // update data mode
    this.data_mode = value;

    var t = this.svg.transition().duration(1500);

    if (value == "quant"){
        // cells
        t.selectAll(".cell")
            .style("fill", function(d) { 
                if (d.hit){ return "red"}
                else {return that.colorScale(d.value); }});

        //legend
        d3.selectAll(".gradienthi").transition().duration(1500).attr("stop-color","blue")
        that.y.domain([this.dmax,this.dmin])
        d3.select(".axis3")
            .transition().duration(1500)
            .call(that.legend_axis)
        d3.selectAll(".legendtext").style("visibility","visible")
    }
    
    else if (value == "bin"){
        // cells
        t.selectAll(".cell")
            .style("fill", function(d){
                if (d.hit) {return "red"}
                else if (d.value > 0){ return "darkgrey";}
                else {return "white";}
            });
        // legend
        d3.selectAll(".gradienthi").transition().duration(1500).attr("stop-color","darkgrey")
        that.y.domain([1,0])
        d3.select(".axis3")
            .transition().duration(1500)
            .call(that.legend_axis)
        d3.selectAll(".legendtext").style("visibility","hidden")
    }
}

// handle show gnps hits button push
HeatMap.prototype.show_hits = function(value){

    var that = this;

    if (this.hits){
        this.hits = false;
        this.svg.selectAll("rect").transition().duration(750).style("fill", function(d){ 
            if (that.data_mode == "quant"){
                return that.colorScale(d.value)
            }
            else {
                if (d.value > 0){ return "darkgrey";}
                else {return "white";}
            }
        })
    }
    
    else {
        this.hits = true;
        this.svg.selectAll("rect").transition().duration(750).style("fill", function(d){ 
            if (that.data_mode == "quant"){
                if (d.hit){ return "red"}
                    else {return that.colorScale(d.value)}
            }
            else {
                if (d.hit){return "red"}
                else if (d.value > 0){ return "darkgrey";}
                else {return "white";}
            }
        })
    }

}

// order cells on user input
HeatMap.prototype.order = function(value){
   
    var that = this;

    // make transitions for sorting by hierarchical clustering
    if(value=="hclust"){
        
        var t = this.svg.transition().duration(1200);
        
        t.selectAll(".cell")
            .attr("x", function(d) { return (d.hcol -1) * that.cellWidth; })
            .attr("y", function(d) { return (d.hrow -1) * that.cellSize; });

        t.selectAll(".rowLabel")
            .attr("y", function (d, i) { return that.cellSize * (that.getProp(that.rowinfo, "sample", d, "hrow") -1); });

        t.selectAll(".colLabel")
            .attr("y", function (d, i) { return that.cellWidth * (that.getProp(that.colinfo, "cmpd", d, "hcol") - 1); });
    }
    // make transitions for sorting by name
    else if (value=="sort_name"){

        var t = this.svg.transition().duration(1200);
        
        t.selectAll(".cell")
            .attr("x", function(d) { return (d.col -1) * that.cellWidth; })
            .attr("y", function(d) { return (d.row -1) * that.cellSize; });

        t.selectAll(".rowLabel")
            .attr("y", function (d, i) { return i * that.cellSize; });

        t.selectAll(".colLabel")
            .attr("y", function (d, i) { return i * that.cellWidth});
   }
   // make transitions for sorting by count across samples
    else if (value=="sort_count"){

        var t = this.svg.transition().duration(1200);
        
        t.selectAll(".cell")
            .attr("x", function(d) { return (d.ccol -1) * that.cellWidth; })

        t.selectAll(".colLabel")
            .attr("y", function (d, i) { return that.cellWidth * (that.getProp(that.colinfo, "cmpd", d, "ccol") - 1); });
   }
}

// order cells on row/column click
HeatMap.prototype.sortbylabel =   function(rORc,i,sortOrder){
       
       var that = this;
       var t = this.svg.transition().duration(1200);
       var all=[];
       var sorted; // sorted is zero-based index
       
       d3.selectAll(".c"+rORc+i) 
         .filter(function(d){
            all.push(d.value);
          });
       
        // sort compounds in a sample
        if(rORc=="r"){ 
         sorted=d3.range(this.col_number).sort(function(a,b){ if(sortOrder){ return all[b]-all[a];}else{ return all[a]-all[b];}});
         t.selectAll(".cell")
           .attr("x", function(d) { return sorted.indexOf(d.col-1) * that.cellWidth; });
         t.selectAll(".colLabel")
          .attr("y", function (d, i) { return sorted.indexOf(i) * that.cellWidth});
       }
       // sort samples for a compound
       else{ 
         sorted=d3.range(this.row_number).sort(function(a,b){if(sortOrder){ return all[b]-all[a];}else{ return all[a]-all[b];}});
         t.selectAll(".cell")
           .attr("y", function(d) { return sorted.indexOf(d.row-1) * that.cellSize; });
         t.selectAll(".rowLabel")
          .attr("y", function (d, i) { return sorted.indexOf(i) * that.cellSize; });
       }
  }