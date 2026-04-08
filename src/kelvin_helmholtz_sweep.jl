using Printf
# Kelvin-Helmholtz Parameter Sweep Script
# This script performs a sweep over velocity (u₀) and shear layer thickness (δ)

include("kelvin_helmholtz_1F.jl")
include("analyze_wavenumber.jl")

# Parameter ranges
u0_values = 0.1:0.1:0.5
delta_values = [0.1, 0.2, 0.4, 0.6, 0.8, 1.0]

@info "Starting parameter sweep..."
@info "u₀ values: $u0_values"
@info "δ values: $delta_values"

# Optional: Store results in a list to print a summary at the end
results = []

for u₀ in u0_values
    for δ in delta_values
        # Create a unique name for the output file
        # Format: u0.1 -> u01, delta 0.1 -> d01
        u_str = @sprintf("%02d", Int(round(u₀ * 10)))
        d_str = @sprintf("%02d", Int(round(δ * 10)))
        output_name = "KH_u$(u_str)_d$(d_str)"

        @info ">>> Starting sweep iteration: u₀=$u₀, δ=$δ ($output_name) <<<"

        try
            # 1. Run the simulation (progress hidden)
            sim_res = run_kh_simulation(u₀, δ; output_name=output_name, stop_time=60, show_progress=false)

            # 2. Run the wavenumber analysis
            @info "Analyzing results for $output_name..."
            analysis = analyze_wavenumber("$(output_name).jld2";
                plot_name="$(output_name)_analysis",
                u₀=u₀, δ=δ)

            if analysis !== nothing
                push!(results, (
                    u0=u₀,
                    delta=δ,
                    m_sim=analysis.m,
                    m_theory=analysis.m_theory,
                    k_sim=analysis.k,
                    k_theory=analysis.k_theory,
                    max_ζ=analysis.max_ζ,
                    sigma_sim=analysis.σ,
                    sigma_theory=analysis.σ_theory
                ))
            end

        catch e
            @error "Failed for u₀=$u₀, δ=$δ: $e"
            Base.show_backtrace(stdout, catch_backtrace())
        end
        println() # space between iterations
    end
end

@info "Parameter sweep and analysis completed!"

# Print a comprehensive summary table
if !isempty(results)
    println("\n" * "="^85)
    println(@sprintf("%-8s %-8s | %-8s | %-12s | %-12s | %-12s",
        "u₀", "δ", "Mode", "k (rad/m)", "Max ζ", "σ (s⁻¹)"))
    println("-"^85)
    for r in results
        println(@sprintf("%-8.1f %-8.1f | %-8d | %-12.4f | %-12.4f | %-12.4f",
            r.u0, r.delta, r.m_sim, r.k_sim, r.max_ζ, r.sigma_sim))
    end
    println("="^85)

    # Export to CSV for Excel comparison
    using DelimitedFiles
    csv_file = "output/kh_sweep_results.csv"
    open(csv_file, "w") do io
        println(io, "u0,delta,m_sim,m_theory,k_sim,k_theory,max_vorticity_sim,sigma_measured,sigma_theory")
        for r in results
            println(io, "$(r.u0),$(r.delta),$(r.m_sim),$(r.m_theory),$(r.k_sim),$(r.k_theory),$(r.max_ζ),$(r.sigma_sim),$(r.sigma_theory)")
        end
    end
    @info "Sweep results exported to $csv_file"
end
