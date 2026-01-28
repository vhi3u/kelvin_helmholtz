using Oceananigans
using Oceananigans.Units
using Printf
using CairoMakie

# ## 1. Grid and Domain
# Increased resolution to 128x128 to capture smaller turbulent eddies (high Re).
grid = RectilinearGrid(size=(128, 128), x=(-5, 5), z=(-5, 5),
    topology=(Periodic, Flat, Bounded))

# ## 2. Parameters and Background State
# Lower viscosity for higher Reynolds number.
Ri = 0.20
h = 0.25
U_scale = 1.0

shear_flow(x, z, t) = U_scale * tanh(z)
stratification(x, z, t) = h * Ri * tanh(z / h)

U = BackgroundField(shear_flow)
B = BackgroundField(stratification)

# ## 3. Model Setup
# Viscosity reduced to 1e-5.
model = NonhydrostaticModel(; grid,
    advection=UpwindBiased(order=5),
    background_fields=(u=U, b=B),
    closure=ScalarDiffusivity(ν=1e-5, κ=1e-5),
    buoyancy=BuoyancyTracer(),
    tracers=:b)

# ## 4. Initialization
u_noise(x, z) = 1e-4 * randn()
w_noise(x, z) = 1e-4 * randn()

set!(model, u=u_noise, w=w_noise)

# ## 5. Simulation Setup
# Using TimeStepWizard for adaptive time-stepping to maintain stability with high resolution.
wizard = TimeStepWizard(cfl=0.5, max_change=1.1, max_Δt=0.1)

simulation = Simulation(model, Δt=0.01, stop_time=150)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

# Progress message
function print_progress(sim)
    u, v, w = sim.model.velocities
    @info @sprintf("Iter: %d, time: %.1f, Δt: %.2e, max|u|: %.3f",
        sim.model.clock.iteration, sim.model.clock.time,
        sim.Δt, maximum(abs, u))
end

simulation.callbacks[:progress] = Callback(print_progress, IterationInterval(100))

# ## 6. Output Writer
u, v, w = model.velocities
b = model.tracers.b

Ω = Field(∂z(u) + ∂z(model.background_fields.velocities.u) - ∂x(w))
B_total = Field(b + model.background_fields.tracers.b)

simulation.output_writers[:fields] = JLD2Writer(model, (; Ω, B=B_total),
    schedule=TimeInterval(0.5),
    filename="kelvin_helmholtz_simple.jld2",
    overwrite_existing=true)

# ## 7. Run Simulation
@info "Starting simulation..."
run!(simulation)
@info "Simulation complete."

# ## 8. Visualization
@info "Visualizing results..."

filepath = "kelvin_helmholtz_simple.jld2"
Ω_timeseries = FieldTimeSeries(filepath, "Ω")
B_timeseries = FieldTimeSeries(filepath, "B")
times = Ω_timeseries.times

n = Observable(1)
Ωₙ = @lift Ω_timeseries[$n]
Bₙ = @lift B_timeseries[$n]
title = @lift @sprintf("t = %.2f", times[$n])

fig = Figure(size=(800, 400))
ax_Ω = Axis(fig[1, 1], title="Total Vorticity")
ax_B = Axis(fig[1, 2], title="Total Buoyancy")

hm_Ω = heatmap!(ax_Ω, Ωₙ, colormap=:balance, colorrange=(-1, 1))
hm_B = heatmap!(ax_B, Bₙ, colormap=:balance, colorrange=(-0.05, 0.05))

Label(fig[0, :], title, fontsize=20)

frames = 1:length(times)
record(fig, "kh_instability_simple.mp4", frames, framerate=12) do i
    n[] = i
end

@info "Visualization saved to kh_instability_simple.mp4"