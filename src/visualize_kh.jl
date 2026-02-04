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
    animate_kh(filename; format=:mp4, framerate=15)

Create an animation of u velocity and vorticity from a KH simulation.

Input is read from `output/` directory.
Output animation is saved to `animation/` directory.

# Arguments
- `filename`: JLD2 filename (e.g., "KH1F.jld2")
- `format`: Output format, either `:mp4` or `:gif` (default: :mp4)
- `framerate`: Animation framerate (default: 15)

# Example
    animate_kh("KH1F.jld2"; format=:gif)
    # Reads: output/KH1F.jld2
    # Saves: animation/KH1F.gif
"""
function animate_kh(filename::String; format=:mp4, framerate=15)
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
    Ω_timeseries = FieldTimeSeries(input_file, "Ω")

    x, y, z = nodes(u_timeseries)

    # Compute max values across all timesteps for consistent colorbar
    u_max = maximum(abs, interior(u_timeseries))
    Ω_max = maximum(abs, interior(Ω_timeseries))

    println("u_max = $u_max, Ω_max = $Ω_max")

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

    println("Recording $format animation to $output_file...")
    record(fig, output_file, 1:length(u_timeseries.times), framerate=framerate) do i
        n[] = i
    end

    println("Animation saved to $output_file")
    return fig
end
