<head>
  <script src="https://unpkg.com/ideogram@1.3.0/dist/js/ideogram.min.js"></script>
  <script src="https://cdnjs.cloudflare.com/ajax/libs/d3/5.7.0/d3.js"></script>
  <style>
  .annot {
    stroke: black;
    #stroke-opacity:0.9;
    stroke-width:0.3;
  }
  svg rect {
   fill: gray;
 }
svg text {
   fill: yellow;
   font: 12px sans-serif;
   text-anchor: end;
}
</style>
</head>
<body>
  <h2>Impute report: ${run_name}</h2>
  <div class="ideogram-container" style="display: inline;"></div>
  <p>
    <div class="progress" style="display: inline;"></div>
  </p>
  <p>
    <div>Timestamp: ${timestamp}</div>
    <div>Chrms: ${chromosomes}</div>
    <div>Progress: ${progress}</div>
    <div>Finished: ${finished_chunks}</div>
    <div>Total: ${total_chunks}</div>
  </p>
  <script>
      var rundetails = ${rundetails};
      var config = {
        container: '.ideogram-container',
        organism: 'human',
        //chromosomes: ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22'],
        chrHeight: 500,
        chrMargin: 2,
        annotations: ${chunkstring},
        annotationsLayout: 'overlay',
        orientation: 'vertical',
        assembly : "GRCh37",
        //fullChromosomeLabels: true,
        //showBandLabels: true,
        //resolution: 550,
        //color : "F333"
      };
    var ideogram = new Ideogram(config);
  </script>
  <script>
     var data = [10, 5, 12, 15];
     var data2 =  {'22': [36,36], '3': [197,198], '5': [164,181], '4': [186,191], '7': [4,160], '6': [138,171], '9': [23,141], '8': [75,147]};
     var width = 800,
        scaleFactor = 20,
        barHeight = 8;
        var data = Object.keys(data2).map(function(key){
             return data2[key][1];
         });
        var graph = d3.select(".progress")
                     .append("svg")
                     .attr("width", width)
                     .attr("height", barHeight * data.length);
        var bar = graph.selectAll("g")
                 .data(data)
                 .enter()
                 .append("g")
                 .attr("transform", function(d, i) {
                    return "translate(0," + i * barHeight + ")";
                 });
       bar.append("rect")
                     .attr("width", function(d) {
                        return d * scaleFactor;
                     })
                     .attr("height", barHeight - 1);

      bar.append("text")
                     .attr("x", function(d) { return (d*scaleFactor); })
                     .attr("y", barHeight / 2)
                     .attr("dy", ".35em")
                     .text(function(d) { return d; });
               var bar2 = graph.selectAll("g")
                            .data(data2)
                            .enter()
                            .append("g")
                            .attr("transform", function(d, i) {
                               return "translate(0," + i[1] * barHeight + ")";
                            });
               bar2.append("rect")
                             .attr("width", function(d) {
                                return d * scaleFactor;
                             })
                             .attr("height", barHeight - 1);


</script>
</body>
