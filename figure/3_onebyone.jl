# t_ = 1950:2022

# temp = data_["FRA"]; x_ = temp.X; y_ = temp.Y;
# grid11 = plot(t_, x_, color = 1); scatter!(grid11, t_[[end]], x_[[end]], msw = 0, color = 1, xticks = [1950, 2030], xlims = [1950, 2030], yticks = [])
# grid21 = plot(t_, y_, color = 1); scatter!(grid21, t_[[end]], y_[[end]], msw = 0, color = 1, xticks = [1950, 2030], xlims = [1950, 2030], yticks = [])
# grid31 = plot(x_, y_, color = 1); scatter!(grid31, x_[[end]], y_[[end]], msw = 0, color = 1, lims = [1e+3, 1e+8], ticks = [1e+3, 1e+8], scale = :log10, aspect_ratio = 1)

# temp = data_["BRA"]; x_ = temp.X; y_ = temp.Y;
# grid12 = plot(t_, x_, color = 2); scatter!(grid12, t_[[end]], x_[[end]], msw = 0, color = 2, xticks = [1950, 2030], xlims = [1950, 2030], yticks = [])
# grid22 = plot(t_, y_, color = 2); scatter!(grid22, t_[[end]], y_[[end]], msw = 0, color = 2, xticks = [1950, 2030], xlims = [1950, 2030], yticks = [])
# grid32 = plot(x_, y_, color = 2); scatter!(grid32, x_[[end]], y_[[end]], msw = 0, color = 2, lims = [1e+3, 1e+8], ticks = [1e+3, 1e+8], scale = :log10, aspect_ratio = 1)

# temp = data_["GMB"]; x_ = temp.X; y_ = temp.Y;
# grid13 = plot(t_, x_, color = 3); scatter!(grid13, t_[[end]], x_[[end]], msw = 0, color = 3, xticks = [1950, 2030], xlims = [1950, 2030], yticks = [])
# grid23 = plot(t_, y_, color = 3); scatter!(grid23, t_[[end]], y_[[end]], msw = 0, color = 3, xticks = [1950, 2030], xlims = [1950, 2030], yticks = [])
# grid33 = plot(x_, y_, color = 3); scatter!(grid33, x_[[end]], y_[[end]], msw = 0, color = 3, lims = [1e+3, 1e+8], ticks = [1e+3, 1e+8], scale = :log10, aspect_ratio = 1)

# temp = data_["MCO"]; x_ = temp.X; y_ = temp.Y;
# grid14 = plot(t_, x_, color = 4); scatter!(grid14, t_[[end]], x_[[end]], msw = 0, color = 4, xticks = [1950, 2030], xlims = [1950, 2030], yticks = [])
# grid24 = plot(t_, y_, color = 4); scatter!(grid24, t_[[end]], y_[[end]], msw = 0, color = 4, xticks = [1950, 2030], xlims = [1950, 2030], yticks = [])
# grid34 = plot(x_, y_, color = 4); scatter!(grid34, x_[[end]], y_[[end]], msw = 0, color = 4, lims = [1e+3, 1e+8], ticks = [1e+3, 1e+8], scale = :log10, aspect_ratio = 1)

# plot(grid11, grid12, grid13, grid14,
#      grid21, grid22, grid23, grid24,
#      grid31, grid32, grid33, grid34,
#      layout = (3, 4), size = (1200, 600), right_margin = 5mm)
# png("temp")