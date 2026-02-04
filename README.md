# Kelvin-Helmholtz Instability Simulations

**MATH 563 - Fluid Mechanics**

Niyati Goswami, Matthew Mackowsky, William Mahoney, Madison McCaskill, and Victor Hieu Nguyen

University of North Carolina at Chapel Hill

---

## Requirements

- **Julia 1.12.0** (or later)

## Installation

### 1. Download Julia

Download and install Julia from the official website: [https://julialang.org/downloads/](https://julialang.org/downloads/)

Or use [juliaup](https://github.com/JuliaLang/juliaup) for version management:
```bash
curl -fsSL https://install.julialang.org | sh
```

### 2. Clone the Repository

```bash
git clone https://github.com/vhi3u/kelvin_helmholtz.git
```

### 3. Navigate to the Project Directory

```bash
cd kelvin_helmholtz
```

### 4. Launch Julia with the Project Environment

```bash
julia --project=.
```

This will automatically load the packages specified in `Project.toml`. On first run, you may need to instantiate the environment:

```julia
using Pkg
Pkg.instantiate()
```

## Usage

### 1. Modify Parameters

Open `src/kelvin_helmholtz_1F.jl` and adjust parameters as needed:

- **Grid resolution**: Modify `size=(128, 128)` in the `RectilinearGrid` constructor
- **Domain size**: Adjust `Lx` and `Lz` values
- **Velocity profile**: Modify `U_0` and `δ` for the shear layer

### 2. Run the Simulation

In the Julia REPL (with project activated):

```julia
include("src/kelvin_helmholtz_1F.jl")
```

### 3. Visualization

After the simulation completes, you can:

- **Generate animation**: Use `src/visualize_kh.jl` to create MP4/GIF animations
- **Plot velocity profile**: Use `src/velocity_profile.jl` to visualize the velocity field
- **Note**: These options are already checked in the code and will run automatically after the simulation. You can uncomment them if you want to generate animations or velocity profiles separately.

## Output

- Simulation data is saved to `output/KH1F.jld2`
- Animations are saved to the `animation/` directory
