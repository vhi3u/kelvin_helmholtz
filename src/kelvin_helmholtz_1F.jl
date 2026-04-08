# download packages if not installed already
# using Pkg
# Pkg.add("Oceananigans")
# Pkg.add("CairoMakie")
# Pkg.add("JLD2")
# Pkg.add("Statistics")
# Pkg.add("Oceanostics")


# import packages
using Oceananigans
using Oceanostics
using Oceanostics.ProgressMessengers: SingleLineMessenger
using Printf
using CairoMakie
using JLD2

# Load utility scripts if they exist
include("visualize_kh.jl")
include("velocity_profile.jl")
include("vorticity_at_center.jl")

"""
    run_kh_simulation(u₀=0.1, δ=0.1; output_name="KH1F", stop_time=60, Lx=10.0, Lz=5.0, resolution=(128, 128))

Runs a 2D Kelvin-Helmholtz simulation with specified velocity magnitude and shear layer thickness.
"""
function run_kh_simulation(u₀=0.1, δ=0.1; output_name="KH1F", stop_time=60, Lx=10.0, Lz=5.0, resolution=(128, 128), show_progress=true)
    # grid setup
    grid = RectilinearGrid(size=resolution, x=(0, Lx), z=(-Lz, 0),
        topology=(Periodic, Flat, Bounded))

    # velocity function setup
    u_init(x, z) = (u₀ / 2) * (1 + tanh((z + Lz / 2) / δ))
    w_pert(x, z) = 1e-4 * randn()

    # model setup
    model = NonhydrostaticModel(grid;
        advection=WENO(),
        closure=nothing
    )

    set!(model, u=u_init, w=w_pert)

    # simulation setup
    simulation = Simulation(model, Δt=0.001, stop_time=stop_time)
    conjure_time_step_wizard!(simulation, cfl=0.8)

    if show_progress
        progress = SingleLineMessenger()
        simulation.callbacks[:progress] = Callback(progress, TimeInterval(5))
    end

    # output setup
    # Manual vorticity calculation for robustness
    u, v, w = model.velocities
    ζ = Field(@at (Face, Center, Face) ∂z(u) - ∂x(w))

    out_file = "output/$(output_name).jld2"
    mkpath("output")

    simulation.output_writers[:fields] = JLD2Writer(model, (; u, w, ζ),
        schedule=TimeInterval(1),
        filename=out_file,
        overwrite_existing=true)

    # Variables to store previous values for dΓ/dt calculation
    prev_Γ = Ref(0.0)
    prev_t = Ref(0.0)
    ΔA = (Lx / resolution[1]) * (Lz / resolution[2])

    # Create a function to print it during the simulation
    function print_circulation(sim)
        compute!(ζ)
        # Total circulation Γ = ∫ ζ dA
        Γ = sum(interior(ζ)) * ΔA
        t = sim.model.clock.time

        dΓdt = t > 0 ? (Γ - prev_Γ[]) / (t - prev_t[]) : 0.0

        if show_progress
            @info @sprintf("Time: %.2f s, Γ: %.6e, dΓ/dt: %.6e", t, Γ, dΓdt)
        end

        prev_Γ[] = Γ
        prev_t[] = t
    end
    simulation.callbacks[:circulation] = Callback(print_circulation, TimeInterval(5))

    # run simulation    
    @info "Starting simulation: u₀=$u₀, δ=$δ, output=$output_name"
    run!(simulation)

    # Calculate max vorticity
    compute!(ζ)
    max_ζ = maximum(abs, interior(ζ))
    @info "Simulation completed. Max vorticity: $max_ζ. Output saved to $out_file"

    # ============ Visualization ============
    try
        # Only save GIF version
        animate_kh("$(output_name).jld2"; format=:gif, time_unit="s", u₀=u₀, δ=δ)
    catch e
        @warn "Visualization failed: $e"
    end

    return (simulation=simulation, max_ζ=max_ζ)
end

# If this file is run directly, execute the default simulation
if abspath(PROGRAM_FILE) == @__FILE__
    run_kh_simulation(0.1, 0.1)
end


