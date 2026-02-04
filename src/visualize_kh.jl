"""
Kelvin-Helmholtz visualization functions.
Usage:
    include("visualize_kh.jl")
    animate_kh("KH1F.jld2"; format=:mp4)
    animate_kh("KH1F.jld2"; format=:gif)
"""

using CairoMakie
using Oceananigans

"""
    animate_kh(filename; format=:mp4, framerate=8, time_unit="s")

Create an animation of u velocity and vorticity from a KH simulation.

Input is read from `output/` directory.
Output animation is saved to `animation/` directory.

# Arguments
- `filename`: JLD2 filename (e.g., "KH1F.jld2")
- `format`: Output format, either `:mp4` or `:gif` (default: :mp4)
- `framerate`: Animation framerate, lower = slower (default: 8)
- `time_unit`: Unit string to display for time (default: "s")

# Example
    animate_kh("KH1F.jld2"; format=:gif, framerate=5, time_unit="seconds")
    # Reads: output/KH1F.jld2
    # Saves: animation/KH1F.gif
"""
function animate_kh(filename::String; format=:mp4, framerate=8, time_unit="s")
    # Build input and output paths
    input_file = joinpath("output", filename)

    # Get base name without extension and add new extension
    base_name = replace(filename, r"\.jld2$" => "")
    ext = format == :gif ? ".gif" : ".mp4"
    output_file = joinpath("animation", base_name * ext)

    # Create animation directory if it doesn't exist
    mkpath("animation")

    println("Loading data from $input_file...")
    u_timeseries = FieldTimeSeries(input_file, "u")
    ζ_timeseries = FieldTimeSeries(input_file, "ζ")

    x, y, z = nodes(u_timeseries)

    # Compute max values across all timesteps for consistent colorbar
    u_max = maximum(abs, interior(u_timeseries))
    ζ_max = maximum(abs, interior(ζ_timeseries))

    println("u_max = $u_max, ζ_max = $ζ_max")

    fig = Figure(size=(900, 850))
    ax_u = Axis(fig[1, 1], title="u velocity", xlabel="x (m)", ylabel="z (m)")
    ax_ζ = Axis(fig[2, 1], title="Relative Vorticity", xlabel="x (m)", ylabel="z (m)")

    # Time display
    times = ζ_timeseries.times
    n = Observable(1)
    time_str = @lift @sprintf("t = %.2f %s", times[$n], time_unit)
    Label(fig[0, :], time_str, fontsize=20, tellwidth=false)

    u_n = @lift interior(u_timeseries[$n], :, 1, :)
    ζ_n = @lift interior(ζ_timeseries[$n], :, 1, :)

    hm_u = heatmap!(ax_u, x, z, u_n, colormap=:balance, colorrange=(-u_max, u_max))
    Colorbar(fig[1, 2], hm_u, label="u (m/s)")

    hm_ζ = heatmap!(ax_ζ, x, z, ζ_n, colormap=:balance, colorrange=(-ζ_max, ζ_max))
    Colorbar(fig[2, 2], hm_ζ, label="ζ (1/s)")

    println("Recording $format animation to $output_file...")
    record(fig, output_file, 1:length(ζ_timeseries.times), framerate=framerate) do i
        n[] = i
    end

    println("Animation saved to $output_file")
    return fig
end
