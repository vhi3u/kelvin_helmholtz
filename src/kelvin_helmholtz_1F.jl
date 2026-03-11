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

# grid setup
Lx, Lz = 10.0, 5.0
grid = RectilinearGrid(size=(128, 128), x=(0, Lx), z=(-Lz, 0),
    topology=(Periodic, Flat, Bounded)) # we use a periodic boundary condition along the x direction to allow flow to appear infinite

# velocity function setup
u₀ = 0.1   # velocity magnitude (m/s)
δ = 0.1  # shear layer thickness (m)
u_init(x, z) = (u₀ / 2) * (1 + tanh((z + Lz / 2) / δ))  # tanh from 0 at bottom to u₀ at top
w_pert(x, z) = 0.01 * randn() # small perturbation

# model setup
# For a purely inviscid fluid, we set closure=nothing. 
# This removes all physical viscosity (ν) and diffusivity (κ).
# Note: Numerical diffusion from the advection scheme (WENO) will still exist 
# at the grid scale to maintain stability.
model = NonhydrostaticModel(grid;
    advection=WENO(),
    closure=nothing
)

set!(model, u=u_init, w=w_pert)

# simulation setup
simulation = Simulation(model, Δt=0.001, stop_time=60) # simulation time in seconds
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

total_vorticity_integral = Field(Integral(ζ))

# Variables to store previous values for dΓ/dt calculation
prev_Γ = Ref(0.0)
prev_t = Ref(0.0)

# 2. Create a function to print it during the simulation
function print_circulation(sim)
    compute!(total_vorticity_integral)
    Γ = total_vorticity_integral[1, 1, 1]
    t = sim.model.clock.time

    # Calculate dΓ/dt via finite difference
    dΓdt = t > 0 ? (Γ - prev_Γ[]) / (t - prev_t[]) : 0.0

    @info @sprintf("Time: %.2f s, Γ: %.6e, dΓ/dt: %.6e", t, Γ, dΓdt)

    prev_Γ[] = Γ
    prev_t[] = t
end
simulation.callbacks[:circulation] = Callback(print_circulation, TimeInterval(5))
# 3. Add it as a callback


# run simulation    
run!(simulation)

# ============ Visualization (optional) ============
include("visualize_kh.jl")
include("velocity_profile.jl")
plot_velocity_profile(u_init, Lz, δ, u₀)
animate_kh("KH1F.jld2"; format=:mp4, time_unit="s")
animate_kh("KH1F.jld2"; format=:gif, time_unit="s")

include("vorticity_at_center.jl")
plot_vorticity_at_point("KH1F.jld2"; x_target=Lx/2, z_target=-Lz/2)


