using CSV
using DataFrames
using CairoMakie
using Printf

function plot_sweep_comparison(csv_path="output/kh_sweep_results.csv")
    if !isfile(csv_path)
        @error "CSV file not found at $csv_path. Please run the sweep script first."
        return
    end

    df = CSV.read(csv_path, DataFrame)
    
    # Filter for a specific delta to show u0 dependence, 
    # and vice versa, or just plot everything with different colors
    
    u0s = unique(df.u0)
    deltas = unique(df.delta)
    
    fig = Figure(size=(1200, 800))
    
    # 1. Growth Rate (sigma) comparison
    ax_sigma = Axis(fig[1, 1], xlabel="Velocity u₀ (m/s)", ylabel="Growth Rate σ (s⁻¹)", 
                    title="Growth Rate: Simulation vs Theory")
    
    # 2. Number of Billows (m) comparison
    ax_m = Axis(fig[1, 2], xlabel="Velocity u₀ (m/s)", ylabel="Number of Billows (m)", 
                title="Dominant Mode: Simulation vs Theory")
    
    # 3. Max Vorticity comparison
    ax_zeta = Axis(fig[2, 1], xlabel="Velocity u₀ (m/s)", ylabel="Max Vorticity ζ (s⁻¹)", 
                   title="Max Vorticity Magnitude")

    colors = [:blue, :red, :green, :orange, :purple, :cyan]
    
    for (i, d) in enumerate(deltas)
        sub = df[df.delta .== d, :]
        sort!(sub, :u0)
        
        color = colors[mod1(i, length(colors))]
        
        # Plot sigma
        scatterlines!(ax_sigma, sub.u0, sub.sigma_measured, color=color, marker=:circle, label="Sim (δ=$d)")
        lines!(ax_sigma, sub.u0, sub.sigma_theory, color=color, linestyle=:dash, label="Theory (δ=$d)")
        
        # Plot m (Billows)
        scatterlines!(ax_m, sub.u0, sub.m_sim, color=color, marker=:circle)
        lines!(ax_m, sub.u0, sub.m_theory, color=color, linestyle=:dash)
        
        # Plot Max Zeta
        scatterlines!(ax_zeta, sub.u0, sub.max_vorticity_sim, color=color, marker=:circle)
    end
    
    # Add a global legend
    Legend(fig[2, 2], ax_sigma, "Parameters", framevisible=false)
    
    save("output/sweep_comparison.png", fig)
    @info "Sweep comparison plot saved to output/sweep_comparison.png"
end

if abspath(PROGRAM_FILE) == @__FILE__
    plot_sweep_comparison()
end
