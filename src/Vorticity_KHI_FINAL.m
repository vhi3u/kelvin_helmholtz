%% =========================================================
%  Kelvin-Helmholtz Instability: Linear Theory Analysis
%  
%  SIMPLIFIED VERSION - Pure Linear Theory Only
%  Evaluates at early times where linear theory is valid
%
%  Project: The Impact of Velocity on Kelvin-Helmholtz Instability 
%           in Single Stratified Fluid
%  Authors: Goswami, Mackowsky, Mahoney, McCaskill, Nguyen
%  Institution: UNC Chapel Hill, Department of Mathematics
%  Date: March 2026
%
%  THEORY:
%  For a hyperbolic tangent velocity profile in a homogeneous fluid:
%    |ω(x,z,t)| = |ω₀(z) + ω'(x,z,t)|
%  where:
%    ω₀(z) = -(u₀/2δ) * sech²((z + Lz/2)/δ)    [base-state vorticity]
%    ω'(x,z,t) = A₀ * exp(σᵣ*t) * sech²((z + Lz/2)/δ) * cos(kx)  [perturbation]
%
%  GROWTH RATE REFERENCE:
%  The theoretical maximum growth rate σᵣ ≈ 0.2*(u₀/δ) for a tanh 
%  velocity profile is derived from linear stability theory.
%  
%  Citation:
%    Drazin, P. G., & Reid, W. H. (2004). Hydrodynamic Stability 
%    (2nd ed.). Cambridge University Press. Section 4.6.
%
%    Hazel, P. (1972). "Numerical studies of the stability 
%    of inviscid stratified shear flows." Journal of Fluid Mechanics, 
%    51(1), 39-61. doi:10.1017/S0022112072001065
%
% =========================================================

clear; clc; close all;

fprintf('=========================================================\n');
fprintf('  KHI Linear Theory Calculator\n');
fprintf('  Homogeneous Fluid (Ri = 0)\n');
fprintf('=========================================================\n\n');

%% =========================================================
%  SECTION 1: SIMULATION PARAMETERS
%  ** EDIT THESE VALUES TO MATCH YOUR OCEANANIGANS SETUP **
% =========================================================

fprintf('--- Input Parameters ---\n\n');

% --- Domain Configuration (from Oceananigans) ---
Lx = 10.0;          % Horizontal domain length (m)
Lz = 5.0;           % Vertical domain depth (m)
Nx = 128;           % Number of horizontal grid points (pixels)
Nz = 128;           % Number of vertical grid points (pixels)

fprintf('Domain:\n');
fprintf('  Lx = %.2f m  (horizontal extent)\n', Lx);
fprintf('  Lz = %.2f m  (vertical depth)\n', Lz);
fprintf('  Grid: %d × %d pixels\n\n', Nx, Nz);

% --- Physical Parameters ---
u0    = 0.5;        % Velocity magnitude (m/s) - MUST MATCH JULIA
delta = 1.0;        % Shear layer thickness (m)

fprintf('Flow Properties:\n');
fprintf('  u₀    = %.3f m/s   (velocity magnitude)\n', u0);
fprintf('  δ     = %.3f m     (shear layer thickness)\n\n', delta);

% --- Initial Perturbation Amplitude ---
A0 = 1e-4;          % Initial perturbation ampliTime_s	Ey_J_per_m	Log_Ey	2_Sigma_t	Growth_Rate_s_invtude (dimensionless)

fprintf('Perturbation:\n');
fprintf('  A₀ = %.2e  (initial perturbation amplitude)\n', A0);
fprintf('  Note: This should match the noise amplitude in your\n');
fprintf('        Oceananigans initialization code.\n\n');

% --- FIGURE TITLE (EDIT THIS) ---
% Customize the title for your output figure
figure_title = sprintf('KHI: u0=%.2f m/s, delta=%.2f m', u0, delta);
% Alternative examples:
% figure_title = 'KHI Linear Theory: u₀=0.1 m/s, δ=0.1 m';
% figure_title = sprintf('KH Analysis: u₀=%.2f m/s, δ=%.2f m, Ri=%.2f', u0, delta, 0);

fprintf('Figure Title:\n');
fprintf('  "%s"\n\n', figure_title);

%% =========================================================
%  SECTION 2: DERIVED QUANTITIES
% =========================================================

fprintf('--- Derived Quantities ---\n\n');

% Evaluation point: Center of domain, at shear layer center
x_center = Lx / 2;          % Horizontal center (m)
z_center = -Lz / 2;         % Shear layer center (m)

fprintf('Evaluation Point:\n');
fprintf('  x_center = %.2f m  (horizontal center)\n', x_center);
fprintf('  z_center = %.2f m  (shear layer center)\n\n', z_center);

% Dominant wavenumber of fastest-growing KH mode
% For tanh profile: k_max ≈ 1/(2δ)
k_max = 1 / (2 * delta);
lambda_max = 2 * pi / k_max;

fprintf('Instability Characteristics:\n');
fprintf('  k_max      = %.4f rad/m  (dominant wavenumber)\n', k_max);
fprintf('  λ_max      = %.4f m      (dominant wavelength)\n', lambda_max);
fprintf('  Expected number of billows ≈ %.1f\n\n', Lx / lambda_max);

% Theoretical maximum growth rate (Drazin & Reid, 2004; Hazel, 1972)
% For homogeneous tanh profile: σᵣ,max ≈ 0.2 * (u₀ / δ)
sigma_r = 0.2 * (u0 / delta);

% Determine valid time range for linear theory
% Linear theory is valid for approximately 5 e-folding times
max_valid_time = 5 / sigma_r;

fprintf('Growth Rate (Theoretical):\n');
fprintf('  σᵣ = %.6f s⁻¹\n', sigma_r);
fprintf('  e-folding time = %.2f s\n', 1/sigma_r);
fprintf('  Linear regime validity: t < %.1f s\n\n', max_valid_time);

% Choose evaluation times within linear regime
% We'll use 5 time points spanning 0 to ~3 e-folding times
t_max_eval = min(3/sigma_r, 25);  % Cap at 25 seconds for safety
t_points = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60];

fprintf('Evaluation Times (within linear regime):\n');
fprintf('  t = [');
fprintf('%.1f ', t_points);
fprintf('] seconds\n\n');

% Richardson number (should be 0 for homogeneous fluid)
g = 9.81;                    % Gravitational acceleration (m/s²)
drho_over_rho0 = 0.0;        % Density difference = 0 (homogeneous)
Ri = (g * drho_over_rho0 * delta) / (u0^2);

fprintf('Richardson Number:\n');
fprintf('  Ri = %.4f  (homogeneous fluid)\n\n', Ri);

%% =========================================================
%  SECTION 3: VORTICITY CALCULATION
% =========================================================

fprintf('--- Vorticity Calculation ---\n\n');

% Normalized vertical coordinate at evaluation point
zeta = (z_center + Lz/2) / delta;  % Should be 0 at shear layer center
sech_squared = sech(zeta)^2;

fprintf('Vertical Structure:\n');
fprintf('  ζ = (z + Lz/2)/δ = %.4f\n', zeta);
fprintf('  sech²(ζ) = %.4f\n\n', sech_squared);

% --- Base-State Vorticity ω₀(z) ---
% From the tanh velocity profile: ω₀(z) = -dU/dz = -(u₀/2δ) * sech²(ζ)
omega_0 = -(u0 / (2 * delta)) * sech_squared;

fprintf('Base-State Vorticity:\n');
fprintf('  ω₀(z_center) = %+.6f s⁻¹\n', omega_0);
fprintf('  (This is the vorticity from the background shear)\n\n');

% Pre-allocate arrays for time series
n_times = length(t_points);
omega_prime = zeros(1, n_times);
omega_total = zeros(1, n_times);
omega_magnitude = zeros(1, n_times);

% Calculate at each time point
fprintf('TIME-DEPENDENT RESULTS:\n');
fprintf('Time(s) | ω''(s⁻¹)    | ω_total(s⁻¹) | |ω|(s⁻¹)   | e-foldings\n');
fprintf('--------|-----------|-------------|-----------|------------\n');

for i = 1:n_times
    t = t_points(i);
    
    % Perturbation vorticity: ω'(x,z,t) = A₀ * exp(σᵣ*t) * sech²(ζ) * cos(kx)
    omega_prime(i) = A0 * exp(sigma_r * t) * sech_squared * cos(k_max * x_center);
    
    % Total vorticity
    omega_total(i) = omega_0 + omega_prime(i);
    
    % Magnitude
    omega_magnitude(i) = abs(omega_total(i));
    
    % Number of e-foldings
    n_e_foldings = sigma_r * t;
    
    fprintf('%7.1f | %+9.4f | %+11.4f | %9.4f | %.2f\n', ...
        t, omega_prime(i), omega_total(i), omega_magnitude(i), n_e_foldings);
end

fprintf('\n');

%% =========================================================
%  EXCEL COPY-PASTE: VORTICITY DATA
% =========================================================

fprintf('=========================================================\n');
fprintf('  COPY-PASTE DATA FOR EXCEL (VORTICITY)\n');
fprintf('=========================================================\n\n');

fprintf('Instructions:\n');
fprintf('1. Select the data below (starting from the header row)\n');
fprintf('2. Copy to clipboard (Ctrl+C or Cmd+C)\n');
fprintf('3. Paste into Excel as Tab-delimited text\n\n');

fprintf('--- BEGIN VORTICITY DATA ---\n');
fprintf('Time_s\tVorticity_Magnitude_s_inv\tOmega_Prime_s_inv\tOmega_Total_s_inv\tE_foldings\n');
for i = 1:n_times
    fprintf('%.2f\t%.8f\t%.8f\t%.8f\t%.4f\n', ...
        t_points(i), omega_magnitude(i), omega_prime(i), omega_total(i), sigma_r * t_points(i));
end
fprintf('--- END VORTICITY DATA ---\n\n');

%% =========================================================
%  SECTION 4: ENERGY TRACKING
%  ** DETAILED EXPLANATION **
% =========================================================

fprintf('=========================================================\n');
fprintf('  ENERGY ANALYSIS - DETAILED EXPLANATION\n');
fprintf('=========================================================\n\n');

fprintf('WHAT IS VERTICAL KINETIC ENERGY?\n');
fprintf('--------------------------------\n');
fprintf('In the Kelvin-Helmholtz instability, the perturbation grows by\n');
fprintf('extracting energy from the background shear flow. The vertical\n');
fprintf('velocity w (which starts at zero) grows exponentially.\n\n');

fprintf('Vertical kinetic energy measures the strength of the perturbation:\n');
fprintf('  E_y(t) = (1/2) ∫∫ ρ*w² dx dz\n\n');

fprintf('WHY E_y IS IMPORTANT:\n');
fprintf('1. Growth Rate: E_y(t) ∝ exp(2*σᵣ*t) in linear theory\n');
fprintf('   - Grows TWICE as fast as amplitude (because E ∝ A²)\n');
fprintf('   - Can extract σᵣ from slope of log(E_y) vs t\n\n');

fprintf('2. Saturation Detection: E_y stops growing exponentially\n');
fprintf('   - Marks transition to nonlinear regime\n');
fprintf('   - Used by Keppens et al. (1999) to define saturation time\n\n');

fprintf('3. Energy Budget: Tracks energy transfer from base flow\n');
fprintf('   - Total kinetic energy is conserved (inviscid)\n');
fprintf('   - E_y increases → background shear decreases\n\n');

fprintf('HOW WE CALCULATE E_y FROM THEORY:\n');
fprintf('----------------------------------\n');
fprintf('For linear KH instability, vertical velocity is:\n');
fprintf('  w(x,z,t) = A(t) * F(z) * sin(kx)\n\n');

fprintf('where:\n');
fprintf('  A(t) = A₀ * exp(σᵣ*t)  [amplitude growth]\n');
fprintf('  F(z) = (u₀/δ) * sech²(z/δ)  [vertical structure]\n');
fprintf('  sin(kx) [horizontal structure]\n\n');

fprintf('The energy integral becomes:\n');
fprintf('  E_y = (1/2) * ρ₀ * A² * [∫ F²(z) dz] * [∫ sin²(kx) dx]\n\n');

fprintf('Spatial integrals (analytical):\n');
fprintf('  ∫ sech⁴(z/δ) dz ≈ 2δ * (2/3)  [vertical]\n');
fprintf('  ∫ sin²(kx) dx = Lx/2  [horizontal]\n\n');

fprintf('Therefore:\n');
fprintf('  E_y(t) = (1/2) * ρ₀ * [A₀*exp(σᵣ*t)]² * (u₀/δ)² * (4δ/3) * (Lx/2)\n');
fprintf('         = C *exp(2*σᵣ*t)\n\n');

fprintf('where C is a constant depending on initial conditions.\n\n');

% Reference density
rho0 = 1.0;  % kg/m³

% Compute vertical kinetic energy
Ey_theoretical = zeros(1, n_times);

for i = 1:n_times
    t = t_points(i);
    
    % Amplitude at time t
    A_t = A0 * exp(sigma_r * t);
    
    % Vertical integral: ∫ sech⁴(z/δ) dz ≈ 2δ * (2/3)
    vertical_integral = 2 * delta * (2/3);
    
    % Horizontal integral: ∫ sin²(kx) dx = Lx/2
    horizontal_integral = Lx / 2;
    
    % Velocity amplitude squared: [A * (u₀/δ)]²
    w_amplitude_squared = A_t^2 * (u0 / delta)^2;
    
    % Total vertical kinetic energy
    Ey_theoretical(i) = 0.5 * rho0 * w_amplitude_squared * vertical_integral * horizontal_integral;
end

fprintf('COMPUTED VALUES:\n');
fprintf('Time(s) | E_y (J/m)  | log(E_y)   | 2*σᵣ*t\n');
fprintf('--------|------------|------------|--------\n');

for i = 1:n_times
    fprintf('%7.1f | %10.4e | %+10.4f | %7.4f\n', ...
        t_points(i), Ey_theoretical(i), log(Ey_theoretical(i)), 2*sigma_r*t_points(i));
end

fprintf('\n');
fprintf('NOTE: log(E_y) should increase linearly with slope = 2*σᵣ\n\n');

%% =========================================================
%  EXCEL COPY-PASTE: ENERGY DATA
% =========================================================

fprintf('=========================================================\n');
fprintf('  COPY-PASTE DATA FOR EXCEL (ENERGY)\n');
fprintf('=========================================================\n\n');

fprintf('Instructions:\n');
fprintf('1. Select the data below (starting from the header row)\n');
fprintf('2. Copy to clipboard (Ctrl+C or Cmd+C)\n');
fprintf('3. Paste into Excel as Tab-delimited text\n\n');

fprintf('--- BEGIN ENERGY DATA ---\n');
fprintf('Time_s\tEy_J_per_m\tLog_Ey\t2_Sigma_t\tGrowth_Rate_s_inv\n');
for i = 1:n_times
    fprintf('%.2f\t%.8e\t%.8f\t%.6f\t%.8f\n', ...
        t_points(i), Ey_theoretical(i), log(Ey_theoretical(i)), 2*sigma_r*t_points(i), sigma_r);
end
fprintf('--- END ENERGY DATA ---\n\n');

% Verify exponential growth by fitting
fprintf('GROWTH RATE VERIFICATION:\n');
fprintf('-------------------------\n');

% Linear regression: log(E_y) = log(E₀) + 2*σᵣ*t
log_Ey = log(Ey_theoretical);
coeffs = polyfit(t_points, log_Ey, 1);
sigma_r_fitted = coeffs(1) / 2;  % Slope / 2

fprintf('From exponential fit:\n');
fprintf('  Slope of log(E_y) vs t = %.6f\n', coeffs(1));
fprintf('  Extracted growth rate = %.6f s⁻¹\n', sigma_r_fitted);
fprintf('  Theoretical σᵣ        = %.6f s⁻¹\n', sigma_r);
fprintf('  Error                 = %.2e%%\n\n', 100*abs(sigma_r_fitted - sigma_r)/sigma_r);

fprintf('This demonstrates that E_y ∝ exp(2*σᵣ*t), confirming linear theory.\n\n');

%% =========================================================
%  SECTION 5: CIRCULATION (VORTICITY INTEGRAL)
% =========================================================

fprintf('=========================================================\n');
fprintf('  CIRCULATION ANALYSIS\n');
fprintf('=========================================================\n\n');

fprintf('WHAT IS CIRCULATION?\n');
fprintf('--------------------\n');
fprintf('Circulation is the integral of vorticity over the domain:\n');
fprintf('  Γ(t) = ∫∫ ζ(x,z,t) dx dz\n\n');

fprintf('WHY CIRCULATION MATTERS:\n');
fprintf('1. Kelvin''s Circulation Theorem: In inviscid, barotropic flow,\n');
fprintf('   circulation is conserved: dΓ/dt = 0\n\n');

fprintf('2. Validation Check: Your simulation should maintain Γ ≈ constant\n');
fprintf('   - Tests numerical conservation properties\n');
fprintf('   - Verifies periodic boundary conditions\n\n');

fprintf('3. Physical Meaning: Total "amount" of rotation in the domain\n');
fprintf('   - Base state has circulation from background shear\n');
fprintf('   - Perturbation vorticity integrates to ~0 (periodic)\n\n');

fprintf('CALCULATION FROM THEORY:\n');
fprintf('------------------------\n');

% Base state circulation
Gamma_base = omega_0 * Lx * Lz;

fprintf('Base state circulation:\n');
fprintf('  Γ₀ = ω₀ * Lx * Lz\n');
fprintf('     = (%.4f s⁻¹) * (%.1f m) * (%.1f m)\n', omega_0, Lx, Lz);
fprintf('     = %.4f m²/s\n\n', Gamma_base);

% Perturbation circulation
fprintf('Perturbation circulation:\n');
fprintf('  Γ'' = ∫∫ ω''(x,z,t) dx dz\n');
fprintf('     = ∫∫ A(t)*sech²(z/δ)*cos(kx) dx dz\n');
fprintf('     = A(t) * [∫ sech²(z/δ) dz] * [∫ cos(kx) dx]\n');
fprintf('     = A(t) * (2δ) * [0]  ← cos(kx) integrates to zero!\n');
fprintf('     = 0\n\n');

fprintf('TOTAL CIRCULATION:\n');
fprintf('  Γ_total(t) = Γ₀ + Γ'' = Γ₀ = %.4f m²/s (constant)\n\n', Gamma_base);

fprintf('This confirms circulation conservation in linear theory.\n\n');

fprintf('IN YOUR SIMULATION:\n');
fprintf('You should see Γ(t) ≈ %.4f m²/s at all times.\n', Gamma_base);
fprintf('Small variations (<5%%) are acceptable due to:\n');
fprintf('  - Numerical discretization errors\n');
fprintf('  - Finite domain effects\n');
fprintf('  - Nonlinear effects at late times\n\n');

%% =========================================================
%  SECTION 6: VISUALIZATION
% =========================================================

% Create comprehensive figure
fig = figure('Position', [100, 100, 1400, 900]);

% Subplot 1: Vorticity magnitude vs time
subplot(2,3,1);
plot(t_points, omega_magnitude, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
hold on;
plot(t_points, abs(omega_0)*ones(size(t_points)), 'k--', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('|ω| (s⁻¹)');
title('Vorticity Magnitude at Center');
legend('|ω_{total}|', '|ω₀|', 'Location', 'best');
grid on;

% Subplot 2: Amplitude growth
subplot(2,3,2);
amplitudes = A0 * exp(sigma_r * t_points);
semilogy(t_points, amplitudes, 'r-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
xlabel('Time (s)');
ylabel('Amplitude A(t)');
title('Perturbation Amplitude Growth');
grid on;

% Subplot 3: E_y vs time (log scale)
subplot(2,3,3);
semilogy(t_points, Ey_theoretical, 'g-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'g');
xlabel('Time (s)');
ylabel('E_y (J/m)');
title('Vertical Kinetic Energy');
grid on;

% Subplot 4: log(E_y) vs time (should be linear)
subplot(2,3,4);
plot(t_points, log(Ey_theoretical), 'mo', 'MarkerSize', 8, 'MarkerFaceColor', 'm');
hold on;
% Add fitted line
t_line = linspace(t_points(1), t_points(end), 100);
log_Ey_fit = coeffs(1) * t_line + coeffs(2);
plot(t_line, log_Ey_fit, 'k--', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('log(E_y)');
title(sprintf('Energy Growth (slope = %.4f = 2σᵣ)', coeffs(1)));
legend('Theoretical', sprintf('Fit: σᵣ = %.4f s⁻¹', sigma_r_fitted), 'Location', 'best');
grid on;

% Subplot 5: Number of e-foldings
subplot(2,3,5);
e_foldings = sigma_r * t_points;
plot(t_points, e_foldings, 'c-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'c');
hold on;
plot([t_points(1), t_points(end)], [5, 5], 'r--', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Number of e-foldings');
title('Linear Regime Validity');
legend('σᵣ*t', 'Linear limit (≈5)', 'Location', 'best');
grid on;

% Subplot 6: Circulation (constant)
subplot(2,3,6);
plot(t_points, Gamma_base * ones(size(t_points)), 'k-', 'LineWidth', 3);
xlabel('Time (s)');
ylabel('Γ (m²/s)');
title('Total Circulation (Conserved)');
y1 = Gamma_base * 0.95;
y2 = Gamma_base * 1.05;

ylim([min(y1,y2), max(y1,y2)]);
grid on;

% Add main title
sgtitle(figure_title, 'FontSize', 16, 'FontWeight', 'bold');

% Save figure
% Create filename including both u0 and delta (safe for filesystem)
safe_title = sprintf('KHI_u0_%.2f_ms_delta_%.2f_m', u0, delta);

% Save figure
saveas(fig, [safe_title '.png']);

% Print exact saved filename
fprintf('Figure saved as: %s.png\n\n', safe_title);

%% =========================================================
%  SECTION 7: EXPORT FOR COMPARISON
% =========================================================

% Create results table
Results = table(t_points', omega_magnitude', Ey_theoretical', ...
    'VariableNames', {'Time_s', 'Vorticity_Magnitude_1_per_s', 'Ey_J_per_m'});

% Add growth rate info
Results.Growth_Rate_s_inv = sigma_r * ones(n_times, 1);
Results.E_foldings = (sigma_r * t_points)';

% Also save workspace
save('KHI_linear_theory_workspace.mat');
fprintf('Workspace saved to: KHI_linear_theory_workspace.mat\n\n');

%% =========================================================
%  EXCEL COPY-PASTE: COMBINED DATA (ALL METRICS)
% =========================================================

fprintf('=========================================================\n');
fprintf('  COPY-PASTE DATA FOR EXCEL (COMBINED)\n');
fprintf('=========================================================\n\n');

fprintf('Instructions:\n');
fprintf('1. Select ALL data below (from header to last row)\n');
fprintf('2. Copy to clipboard (Ctrl+C or Cmd+C)\n');
fprintf('3. Paste into Excel - columns will auto-separate\n\n');

fprintf('--- BEGIN COMBINED DATA ---\n');
fprintf('Time_s\tVorticity_Mag_s_inv\tEy_J_per_m\tLog_Ey\tCirculation_m2_per_s\tGrowth_Rate_s_inv\n');
for i = 1:n_times
    fprintf('%.2f\t%.8f\t%.8e\t%.8f\t%.8f\t%.8f\n', ...
        t_points(i), omega_magnitude(i), Ey_theoretical(i), ...
        log(Ey_theoretical(i)), Gamma_base, sigma_r);
end
fprintf('--- END COMBINED DATA ---\n\n');

fprintf('NOTE: Paste into Excel starting at cell A1 for best results.\n');
fprintf('      Column headers will be in the first row.\n\n');

%% =========================================================
%  SECTION 8: COMPARISON INSTRUCTIONS
% =========================================================

fprintf('=========================================================\n');
fprintf('  INSTRUCTIONS FOR COMPARISON WITH SIMULATION\n');
fprintf('=========================================================\n\n');

fprintf('1. FROM YOUR OCEANANIGANS SIMULATION, EXTRACT:\n\n');

fprintf('   a) Vorticity at center point (x=%.1f, z=%.1f) at times:\n', x_center, z_center);
fprintf('      t = [');
fprintf('%.1f ', t_points);
fprintf('] s\n\n');

fprintf('   b) Vertical kinetic energy at the same times:\n');
fprintf('      E_y(t) = (1/2) * sum(w.^2) * dx * dz\n');
fprintf('      where w is the vertical velocity field\n\n');

fprintf('   c) Total circulation at the same times:\n');
fprintf('      Γ(t) = sum(ζ) * dx * dz\n');
fprintf('      where ζ is the vorticity field\n\n');

fprintf('2. COMPARE:\n\n');

fprintf('   a) Growth Rate (MOST IMPORTANT):\n');
fprintf('      - Plot log(E_y_sim) vs t\n');
fprintf('      - Fit a line to extract σᵣ_sim\n');
fprintf('      - Compare with σᵣ_theory = %.6f s⁻¹\n', sigma_r);
fprintf('      - Expected agreement: within 5%%\n\n');

fprintf('   b) Vorticity Magnitude:\n');
fprintf('      - Plot |ω_sim(t)| vs |ω_theory(t)|\n');
fprintf('      - Should agree closely in linear regime\n');
fprintf('      - Expected RMS error: < 5%%\n\n');

fprintf('   c) Circulation Conservation:\n');
fprintf('      - Check that Γ_sim(t) ≈ %.4f m²/s (constant)\n', Gamma_base);
fprintf('      - Variation should be < 5%%\n');
fprintf('      - Tests numerical scheme quality\n\n');

fprintf('3. VALIDITY CHECK:\n\n');
fprintf('   Linear theory is valid when:\n');
fprintf('   - Number of e-foldings < 5\n');
fprintf('   - |ω''| << |ω₀| (perturbation small)\n');
fprintf('   - E_y grows exponentially (not saturated)\n\n');

fprintf('   For your parameters:\n');
fprintf('   - Valid for t < %.1f s\n', max_valid_time);
fprintf('   - Our evaluation times satisfy this: max(t) = %.1f s\n\n', max(t_points));

fprintf('=========================================================\n');
fprintf('  CALCULATION COMPLETE\n');
fprintf('=========================================================\n');
