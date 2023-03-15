using Plots, LaTeXStrings; mm = Plots.mm
# default(fontfamily = "Computer Modern")
default(fontfamily = "맑은 고딕")
using Colors, ColorSchemes

red = colorant"#C00000"
orange = colorant"#ED7D31"
yellow = colorant"#FFC000"
애플망고 = cgrad([red, orange, yellow])
무지개 = ColorSchemes.rainbow