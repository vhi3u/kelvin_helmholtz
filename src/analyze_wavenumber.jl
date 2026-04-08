"""
Calculates the spatial wavenumber and wavelength of Kelvin-Helmholtz eddies.
This is done by taking a horizontal slice of vorticity at the shear layer center 
and performing a Fourier Transform to find the dominant spatial frequency.
"""

using Oceananigans
using JLD2
using FFTW
using CairoMakie
using Printf
using Statistics

function analyze_wavenumber(filename="KH1F.jld2"; time_index=nothing, plot_name="wavenumber_analysis", u₀=nothing, δ=nothing)
    filepath = joinpath("output", filename)
    if !isfile(filepath)
        @warn "File $filepath not found. Skipping analysis."
        return nothing
    end

    @printf("\n--- Starting Physics-Based Analysis for %s ---\n", filename)
    ζ_timeseries = FieldTimeSeries(filepath, "ζ")
    w_timeseries = FieldTimeSeries(filepath, "w")

    times = ζ_timeseries.times
    grid = ζ_timeseries.grid
    Lx, Lz = grid.Lx, grid.Lz
    Nx, Nz = grid.Nx, grid.Nz
    ΔA = (Lx / Nx) * (Lz / Nz)

    # 1. Theoretical Values (Reference: Drazin & Reid, 2004)
    # Theory for tanh profile: σ ≈ 0.2 * (u₀/δ), Γ = -u₀ * Lx, k ≈ 1/(2δ)
    if u₀ !== nothing && δ !== nothing
        σ_theory = 0.2 * (u₀ / δ)
        Γ_theory = -u₀ * Lx
        k_theory = 1 / (2 * δ)
        ω₀_z_center = -(u₀ / (2 * δ)) # at z = -Lz/2, sech(0) = 1
    else
        σ_theory, Γ_theory, k_theory, ω₀_z_center = 0.0, 0.0, 0.0, 0.0
    end

    # Use the final time if not specified
    n_final = (time_index === nothing) ? length(times) : time_index
    t_curr = times[n_final]

    # Indices for diagnostics
    z_nodes = nodes(ζ_timeseries)[3]
    k_center = argmin(abs.(z_nodes .+ Lz / 2))
    x_nodes = nodes(ζ_timeseries)[1]
    i_center = argmin(abs.(x_nodes .- Lx / 2))

    # 2. Time-series diagnostics
    @info "Calculating integrated diagnostics (VKE, Circulation, Amplitude)..."
    ρ₀ = 1.0 # Reference density
    vke = Float64[]
    circulation = Float64[]
    amplitude = Float64[]
    ζ_center_t = Float64[]

    for i in 1:length(times)
        ζ_field = interior(ζ_timeseries[i])
        w_field = interior(w_timeseries[i])

        # VKE = sum(0.5 * ρ₀ * w² * ΔA)
        push!(vke, 0.5 * ρ₀ * sum(w_field .^ 2) * ΔA)
        push!(circulation, sum(ζ_field) * ΔA)
        push!(amplitude, maximum(abs, w_field))
        push!(ζ_center_t, Float64(ζ_timeseries[i][i_center, 1, k_center]))
    end

    # 3. Growth Rate Analysis (Linear fit to log(VKE))
    # Detect the "True" Growth Regime:
    # 1. First, find where energy is actually growing (ignore initial noise decay)
    valid_idx = findall(x -> x > 1e-25, vke)
    if isempty(valid_idx)
        σ_measured = 0.0
        fit_idx = []
        β = [0.0, 0.0]
    else
        # Find the minimum VKE (trough) before growth starts
        trough_idx_in_valid = argmin(vke[valid_idx])
        start_fit_idx = valid_idx[trough_idx_in_valid]
        
        # Fit from the trough to either the end or where saturation starts (approx)
        # We take a window of samples after the trough
        end_fit_idx = min(length(times), start_fit_idx + 30) 
        fit_idx = start_fit_idx:end_fit_idx
        
        if length(fit_idx) > 3
            X = [ones(length(fit_idx)) times[fit_idx]]
            Y = log.(vke[fit_idx])
            β = X \ Y
            σ_measured = β[2] / 2 # Slope of log(E) is 2σ
        else
            σ_measured = 0.0
            β = [0.0, 0.0]
        end
    end

    # 4. Spatial Analysis (Final time)
    ζ_slice = interior(ζ_timeseries[n_final], :, 1, k_center)
    ζ_prime = ζ_slice .- sum(ζ_slice) / length(ζ_slice)
    ζ_ft = rfft(ζ_prime)
    power_spectrum = abs2.(ζ_ft)

    modes = 0:length(power_spectrum)-1
    k_vals = 2π .* modes ./ Lx
    idx_max = length(power_spectrum) > 1 ? argmax(power_spectrum[2:end]) + 1 : 1
    m_dominant = modes[idx_max]
    k_dominant = k_vals[idx_max]
    λ_dominant = (m_dominant > 0) ? 2π / k_dominant : Inf
    max_ζ = maximum(abs, interior(ζ_timeseries))
    @printf("Analysis Metrics:\n")
    @printf("  Sim σ: %.4f | Theory σ: %.4f (diff: %.1f%%)\n", σ_measured, σ_theory, 100 * abs(σ_measured - σ_theory) / σ_theory)
    @printf("  Sim k: %.4f | Theory k: %.4f\n", k_dominant, k_theory)
    @printf("  Circulation: %.4f (Theory: %.4f)\n", mean(circulation), Γ_theory)

    # 5. Visualization (6 Panels Modeling reference figure)
    fig = Figure(size=(1400, 950), backgroundcolor=:white)

    # Use a GridLayout for robust column management
    gl = fig[1, 1] = GridLayout()

    title_text = "Kelvin-Helmholtz Stability Analysis: Simulation vs Linear Theory"
    if u₀ !== nothing && δ !== nothing
        title_text *= " (u₀ = $u₀ m/s, δ = $δ m)"
    end
    Label(fig[0, :], title_text, fontsize=30, font=:bold, padding=(0, 0, 15, 15))

    # (1, 1): Vorticity Magnitude at Center
    ax1 = Axis(gl[1, 1], xlabel="Time (s)", ylabel="|ζ| (1/s)", title="Vorticity magnitude at Center", alignmode=Outside())
    lines!(ax1, times, abs.(ζ_center_t), color=:blue, linewidth=3, label="Simulation")
    axislegend(ax1, position=:lt)

    # (1, 2): Perturbation Amplitude Growth
    ax2 = Axis(gl[1, 2], xlabel="Time (s)", ylabel="Amplitude A(t)", title="Perturbation Amplitude Growth", yscale=log10, alignmode=Outside())
    scatterlines!(ax2, times, amplitude, color=:red, markersize=8)

    # (1, 3): Vertical Kinetic Energy
    ax3 = Axis(gl[1, 3], xlabel="Time (s)", ylabel="E_z (J/m)", title="VKE: ∫ 0.5 ρ w² dA", yscale=log10, alignmode=Outside())
    scatterlines!(ax3, times, vke, color=:green, markersize=8)

    # (2, 1): Energy Growth & Fit (log scale)
    ax4 = Axis(gl[2, 1], xlabel="Time (s)", ylabel="log(E_y)", 
               title=@sprintf("Growth Rate Extraction (σ ≈ %.4f)", σ_measured), alignmode=Outside())
    lines!(ax4, times, log.(vke), color=:black, linestyle=:dash, label="log(VKE)")
    if length(fit_idx) > 1
        lines!(ax4, times[fit_idx], β[1] .+ (2 * σ_measured) .* times[fit_idx], color=:magenta, linewidth=4, label="Linear Growth Fit")
    end
    axislegend(ax4, position=:rb)

    # (2, 2): Spatial Power Spectrum
    max_mode_plot = max(10, Int(round(m_dominant * 2.5)))
    ax5 = Axis(gl[2, 2], xlabel="Mode number (m)", ylabel="Power |ζ̂|²", title="Dominant Mode Analysis",
        xticks=0:max(1, Int(round(max_mode_plot / 5))):max_mode_plot, alignmode=Outside())
    stem!(ax5, modes[1:max_mode_plot+1], power_spectrum[1:max_mode_plot+1], color=:crimson)
    if m_dominant > 0
        vlines!(ax5, [m_dominant], color=:black, linestyle=:dash)
        text!(ax5, m_dominant, maximum(power_spectrum) * 0.8, text="Mode $m_dominant", align=(:left, :top))
    end

    # (2, 3): Perturbation Circulation (Conserved)
    ax6 = Axis(gl[2, 3], xlabel="Time (s)", ylabel="Γ' (m²/s)", title="Perturbation Circulation Γ' (Conserved)", alignmode=Outside())
    
    # Calculate perturbation circulation relative to the first time step
    Γ_prime = circulation .- circulation[1]
    
    lines!(ax6, times, Γ_prime, color=:black, linewidth=4, label="Sim Γ'")
    
    # Center y-axis around zero to show how small the variation is
    max_dev = max(1e-12, maximum(abs, Γ_prime))
    ylims!(ax6, -1.2*max_dev, 1.2*max_dev)
    
    axislegend(ax6, position=:rt)

    # Fix the layout squashing
    for i in 1:3
        colsize!(gl, i, Relative(1 / 3))
    end
    rowgap!(gl, 15)
    colgap!(gl, 15)

    plot_path = joinpath("output", "$(plot_name).png")
    save(plot_path, fig)
    @printf("Comprehensive summary plot saved to %s\n", plot_path)

    λ_theory = (k_theory > 0) ? 2π / k_theory : Inf
    return (
        m=m_dominant, k=k_dominant, λ=λ_dominant, t=t_curr, max_ζ=max_ζ, σ=σ_measured,
        σ_theory=σ_theory, k_theory=k_theory, m_theory=Lx/λ_theory
    )
end

# Example usage when script is run directly
if abspath(PROGRAM_FILE) == @__FILE__
    analyze_wavenumber()
end
