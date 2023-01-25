coordinate = DataFrame([
"서울" 126.9783882  37.5666103  5   5
"부산" 129.0750223  35.1798160  9   3
"대구" 128.6017630  35.8713900  8   8
"인천" 126.7051505  37.4559418  3   6
"광주" 126.8513380  35.1600320  6   2
"대전" 127.3849508  36.3504396  2   8
"울산" 129.3112994  35.5394773  9   5
"세종" 127.2894325  36.4803512  3   9
"경기" 127.5508020  37.4363177  6   8
"강원" 128.3115261  37.8603672  4   9
"충북" 127.6551404  36.7853718  1   9
"충남" 126.8453965  36.6173379  1   7
"전북" 127.2368291  35.6910153  3   4
"전남" 126.9571667  34.9007274  3   3
"경북" 128.9625780  36.6308397  9   7
"경남" 128.2417453  35.4414209  8   1
"제주" 126.5758344  33.4273366  1   3
], ["locaiton", "lon", "lat", "x", "y"])
x = coordinate.x .- 5 + 0.1cos.(2 * (1:17))
y = coordinate.y .- 5 + 0.1cos.(2 * (1:17))

pp02_1_ = []
for t = 2012:2070
    mtrx_mgrn_ = Dict()
    vctr_mgrn = mgrn[:, "y$t"]
    mtrx_mgrn = reshape(sum(reshape(vctr_mgrn, 17,2,17,17), dims = 1:2), 17, 17)
    mtrx_mgrn .*= (1 .- I(17))

    vrtx_wght = sqrt.(combine(groupby(rslt, :location), "y$t" => sum => "y$t")[:,"y$t"])/80
    edge_wght = mtrx_mgrn' - mtrx_mgrn
    edge_wght[edge_wght .< 0] .= 0
    edge_wght = edge_wght ./ maximum(edge_wght)
    pp02_1 = scatter(x,y,text = text.(1:17, 8, :black, "Computer Modern"),
        ms = vrtx_wght, color = :white,
        title = "$t year",
        grid = false, lims = (-6, 6), legend = :bottomleft, label = :none, size = (500, 500), axis = false, dpi = 200)
    plot!([-100,-100],[-100,-100],arrow = true, color = :red, label = ">0.2")
    plot!([-100,-100],[-100,-100],arrow = true, color = :black, label = ">0.1")
    plot!([-100,-100],[-100,-100],arrow = true, color = :black, ls = :dash, label = ">.0.01")
    for i = 1:17
        for j = 1:17
            if edge_wght[i,j] > 0.01
                Δy = y[j] - y[i]
                Δx = x[j] - x[i]
                slope = Δy / Δx
                θ = atan(slope)
                if Δy > 0
                    x1 = x[i] + sign(slope) * cos(θ) * vrtx_wght[i] / 40
                    y1 = y[i] + sign(slope) * sin(θ) * vrtx_wght[i] / 40
                    x2 = x[j] - sign(slope) * cos(θ) * vrtx_wght[j] / 40
                    y2 = y[j] - sign(slope) * sin(θ) * vrtx_wght[j] / 40
                else
                    x1 = x[i] - sign(slope) * cos(θ) * vrtx_wght[i] / 40
                    y1 = y[i] - sign(slope) * sin(θ) * vrtx_wght[i] / 40
                    x2 = x[j] + sign(slope) * cos(θ) * vrtx_wght[j] / 40
                    y2 = y[j] + sign(slope) * sin(θ) * vrtx_wght[j] / 40
                end
                
                if edge_wght[i,j] > 0.2
                    ls_ = :solid
                    color_ = :red
                elseif edge_wght[i,j] > 0.1
                    ls_ = :solid
                    color_ = :black
                else
                    ls_ = :dash
                    color_ = :black
                end
                    
                plot!(pp02_1, [x1,x2],[y1,y2], arrow = true, lw = 1, color = color_, ls = ls_, label = :none)
            end
        end
    end
    push!(pp02_1_, pp02_1)
    png(pp02_1, "G:/figure/snapshot/02 00 $t.png")
end
png(plot(pp02_1_[1], pp02_1_[17], pp02_1_[33], pp02_1_[50], layout = (2,2), size = (1000, 1000), dpi = 200),
    "G:/figure/02 00.png")