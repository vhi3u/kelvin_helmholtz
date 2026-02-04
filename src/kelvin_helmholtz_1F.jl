# # download packages if not installed already
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

# grid setup
Lx, Lz = 10.0, 5.0
grid = RectilinearGrid(size=(256, 256), x=(-2 * Lx, 2 * Lx), z=(-Lz, 0),
    topology=(Periodic, Flat, Bounded)) # we use a periodic boundary condition along the x direction to allow flow to appear infinite

# velocity function setup
u₀ = 1.0   # velocity magnitude (m/s)
δ = 0.1  # shear layer thickness (m)
u_init(x, z) = u₀ * tanh((z + Lz / 2) / δ) # tanh function looks like positive u₀ at top, negative u₀ at bottom, and quick transition layer in between.
w_pert(x, z) = 0.01 * randn() # small perturbation

# model setup
model = NonhydrostaticModel(; grid,
    advection=WENO(),
    closure=ScalarDiffusivity(ν=1e-5, κ=1e-5), # ν can be tuned (kinematic viscosity), κ can be tuned (thermal diffusivity)
)

set!(model, u=u_init, w=w_pert)

# simulation setup
simulation = Simulation(model, Δt=0.001, stop_time=120) # simulation time in seconds
conjure_time_step_wizard!(simulation, cfl=0.8)

progress = SingleLineMessenger()
simulation.callbacks[:progress] = Callback(progress, TimeInterval(5)) # frequency of progress updates in seconds

# output setup
u, v, w = model.velocities
ζ = Field(∂z(u) - ∂x(w))  # relative vorticity

simulation.output_writers[:fields] = JLD2Writer(model, (; u, w, ζ),
    schedule=TimeInterval(1),
    filename="output/KH1F.jld2",
    overwrite_existing=true)


# run simulation    
run!(simulation)

# ============ Visualization (optional) ============
include("visualize_kh.jl")
include("velocity_profile.jl")
plot_velocity_profile(u_init, Lz, δ, u₀)
animate_kh("KH1F.jld2"; format=:mp4, time_unit="s")
animate_kh("KH1F.jld2"; format=:gif, time_unit="s")


