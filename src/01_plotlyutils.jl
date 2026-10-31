models = unique(result_fps.model)

traces = join((
    begin
        sub = result_fps[result_fps.model .== m, :]
        """{type: "scatter3d", mode: "markers", name: "$m",
            x: $(js(sub.x)), y: $(js(sub.y)), z: $(js(sub.z)),
            marker: {size: 3}}"""
    end
    for m in models
), ",")

plotly_html(
    "[" * traces * "]",
    """{
      scene: {
        xaxis: {range: [5, 9]},
        yaxis: {range: [5, 9]},
        zaxis: {range: [5, 9]}
      }
    }""",
)