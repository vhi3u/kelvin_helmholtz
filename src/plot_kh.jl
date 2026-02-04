using CairoMakie
using Oceananigans

# ============ Plotting ============
u_timeseries = FieldTimeSeries("kelvin_helmholtz_1F.jld2", "u")
Ω_timeseries = FieldTimeSeries("kelvin_helmholtz_1F.jld2", "Ω")

x, y, z = nodes(u_timeseries)

# Compute max values across all timesteps for consistent colorbar
u_max = maximum(abs, interior(u_timeseries))
Ω_max = maximum(abs, interior(Ω_timeseries))

fig = Figure(size=(900, 800))
ax_u = Axis(fig[1, 1], title="u velocity", xlabel="x", ylabel="z")
ax_Ω = Axis(fig[2, 1], title="Vorticity", xlabel="x", ylabel="z")

n = Observable(1)
u_n = @lift interior(u_timeseries[$n], :, 1, :)
Ω_n = @lift interior(Ω_timeseries[$n], :, 1, :)

hm_u = heatmap!(ax_u, x, z, u_n, colormap=:balance, colorrange=(-u_max, u_max))
Colorbar(fig[1, 2], hm_u)

hm_Ω = heatmap!(ax_Ω, x, z, Ω_n, colormap=:balance, colorrange=(-Ω_max, Ω_max))
Colorbar(fig[2, 2], hm_Ω)

record(fig, "kelvin_helmholtz_1F.mp4", 1:length(u_timeseries.times), framerate=15) do i
    n[] = i
end
