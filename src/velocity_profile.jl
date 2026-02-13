"""
Velocity profile plotting functions.
Usage:
    include("velocity_profile.jl")
    plot_velocity_profile(u_func, Lz, δ, u₀; output_file="output/velocity_profile.png")
"""

using CairoMakie
using Printf

"""
    plot_velocity_profile(u_func, Lz, δ, u₀; output_file="output/velocity_profile.png", show_table=true)

Plot and print the initial velocity profile U(z).

# Arguments
- `u_func`: Velocity function u(x, z)
- `Lz`: Vertical domain extent
- `δ`: Shear layer thickness
- `u₀`: Velocity magnitude
- `output_file`: Path for output PNG file
- `show_table`: Print velocity table to console (default: true)

# Returns
- `fig`: The CairoMakie Figure object
"""
function plot_velocity_profile(u_func, Lz, δ, u₀;
    output_file="output/velocity_profile.png",
    show_table=true)

    z_vals = range(-Lz, 0, length=50)
    u_profile = [u_func(0, z) for z in z_vals]

    if show_table
        println("\n=== Initial Velocity Profile ===")
        println("z (m)      | u (m/s)")
        println("-----------|---------")
        for z in [-Lz, -3 * Lz / 4, -Lz / 2, -Lz / 4, 0]
            @printf("%.2f      | %.3f\n", z, u_func(0, z))
        end
        println("\nShear layer thickness δ = $δ m")
        println("Max shear ∂u/∂z ≈ $(u₀/δ) s⁻¹ at z = $(-Lz/2) m")
    end

    # Plot velocity profile
    fig = Figure(size=(400, 500), backgroundcolor=:transparent)
    # ax = Axis(fig[1, 1],
    #     xlabel="u velocity (m/s)",
    #     ylabel="z (m)",
    #     title="Initial Velocity Profile U(z)")
    ax = Axis(fig[1, 1], backgroundcolor=:transparent)
    hidedecorations!(ax)
    hidespines!(ax)
    # lines!(ax, u_profile, collect(z_vals), linewidth=3, color=:blue)
    lines!(ax, u_profile, collect(z_vals), linewidth=6, color=:black, linestyle=:dash)
    # hlines!(ax, [-Lz / 2], color=:red, linestyle=:dash, label="Shear layer center")
    # vlines!(ax, [0], color=:gray, linestyle=:dot)
    # axislegend(ax, position=:lt)

    # # Draw horizontal arrows that scale with local velocity
    # n_arrows = 10
    # z_arrow = range(-Lz, 0, length=n_arrows)
    # x_start = zeros(n_arrows)           # arrows start at x=0
    # u_arrows = [u_func(0, z) for z in z_arrow]  # arrow length = local velocity
    # arrows!(ax, x_start, collect(z_arrow), u_arrows, zeros(n_arrows);
    #     tipwidth=15, tiplength=5, shaftwidth=6, color=:black)

    save(output_file, fig)
    println("\nVelocity profile saved to $output_file")

    return fig
end
