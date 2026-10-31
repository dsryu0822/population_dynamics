
js(v::AbstractVector{<:Real}) = "[" * join(v, ",") * "]"
js(v::AbstractVector) = "[" * join(("\"$x\"" for x in v), ",") * "]"

function plotly_html(traces::String, layout::String = "{}"; file = "plot.html")
    write(file, """
    <!DOCTYPE html>
    <html>
    <head>
      <meta charset="utf-8">
      <script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
    </head>
    <body>
      <div id="plot" style="width:100%;height:95vh"></div>
      <script>Plotly.newPlot("plot", $traces, $layout);</script>
    </body>
    </html>
    """)
end