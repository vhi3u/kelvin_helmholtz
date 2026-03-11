"""
Vorticity analysis script.
Plots the time series of vorticity at a specific point (x, z).
"""

using CairoMakie
using Oceananigans
using JLD2
using Printf

function plot_vorticity_at_point(filename="KH1F.jld2"; x_target=5.0, z_target=-2.5)
    filepath = joinpath("output", filename)
    if !isfile(filepath)
        error("File $filepath not found. Run the simulation first.")
    end

    println("Loading vorticity time series from $filepath...")
    ζ_timeseries = FieldTimeSeries(filepath, "ζ")
    times = ζ_timeseries.times
    
    # Get grid nodes
    x_nodes, y_nodes, z_nodes = nodes(ζ_timeseries)
    
    # Find indices closest to target coordinates
    i_idx = argmin(abs.(x_nodes .- x_target))
    k_idx = argmin(abs.(z_nodes .- z_target))
    
    actual_x = x_nodes[i_idx]
    actual_z = z_nodes[k_idx]
    
    println("Target: ($x_target, $z_target)")
    println("Actual: ($actual_x, $actual_z) at indices ($i_idx, $k_idx)")

    # Extract time series at that point
    # interior(ζ_timeseries) has shape (Nx, Ny, Nz, Nt)
    # ζ_at_point = [ζ_timeseries[n][i_idx, 1, k_idx] for n in 1:length(times)]
    # A more efficient way to extract a single point from a FieldTimeSeries:
    ζ_at_point = Float64[]
    for n in 1:length(times)
        push!(ζ_at_point, ζ_timeseries[n][i_idx, 1, k_idx])
    end

    # --- 1. Time Series Plot ---
    fig_ts = Figure(size=(800, 500))
    ax_ts = Axis(fig_ts[1, 1],
        xlabel = "Time (s)",
        ylabel = "Vorticity ζ (1/s)",
        title = @sprintf("Vorticity Time Series at (x=%.2f, z=%.2f)", actual_x, actual_z))
    
    lines!(ax_ts, times, ζ_at_point, linewidth=3, color=:crimson)
    save(joinpath("output", "vorticity_time_series.png"), fig_ts)

    # --- 2. Spatial Slices at Final Time ---
    n_final = length(times)
    ζ_final = ζ_timeseries[n_final]
    
    # Horizontal slice at z = -Lz/2
    fig_h = Figure(size=(800, 500))
    ax_h = Axis(fig_h[1, 1],
        xlabel = "x (m)",
        ylabel = "Vorticity ζ (1/s)",
        title = @sprintf("Horizontal Vorticity Slice at z=%.2f (t=%.2f s)", actual_z, times[n_final]))
    
    lines!(ax_h, x_nodes, interior(ζ_final, :, 1, k_idx), linewidth=3, color=:blue)
    save(joinpath("output", "vorticity_slice_horizontal.png"), fig_h)

    # Vertical slice at x = Lx/2
    fig_v = Figure(size=(800, 500))
    ax_v = Axis(fig_v[1, 1],
        xlabel = "z (m)",
        ylabel = "Vorticity ζ (1/s)",
        title = @sprintf("Vertical Vorticity Slice at x=%.2f (t=%.2f s)", actual_x, times[n_final]))
    
    lines!(ax_v, z_nodes, interior(ζ_final, i_idx, 1, :), linewidth=3, color=:green)
    save(joinpath("output", "vorticity_slice_vertical.png"), fig_v)

    println("All plots saved to output/ directory.")
    
    return fig_ts
end
