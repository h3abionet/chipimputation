<head>
  <script src="https://unpkg.com/ideogram@1.3.0/dist/js/ideogram.min.js"></script>
  <script src="https://cdnjs.cloudflare.com/ajax/libs/d3/5.7.0/d3.js"></script>
  <style>
  .annot {
    stroke: black;
    #stroke-opacity:0.9;
    stroke-width:0.3;
  }
</style>
</head>
<body>
  <h2>Impute report: ${run_name}</h2>
  <div class="ideogram-container" style="display: inline;"></div>
  <p>
    <div>Timestamp: ${timestamp}</div>
    <div>Chrms: ${chromosomes}</div>
    <div>Progress: ${progress}</div>
  </p>
  <script>
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
</body>
