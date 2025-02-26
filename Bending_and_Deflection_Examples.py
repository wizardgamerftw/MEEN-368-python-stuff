"""
These are all the Bending and Deflection Functions with example calculations

Each example prints the input values along with the computed result

Don't run them in here!! Copy and paste them into the calculations python file, or you'll run this whole thing!!
"""

from Bending_and_Deflection_Functions import *
from sympy import symbols

x = symbols('x')

# ----------------------------------------------------
# Example 1: Normal Stress
# ----------------------------------------------------
# Compute normal bending stress using:  σ = -M * y / I
M = 500         # Bending moment in Nm
y = 0.05        # Distance from the neutral axis in m
I = 8e-6        # Moment of inertia in m^4
sigma = normal_stress(M, y, I)
print("Example 1: Normal Stress")
print(f"Inputs: M = {M} Nm, y = {y} m, I = {I} m^4")
print(f"Normal stress, σ = {sigma:.3f} Pa\n")

# ----------------------------------------------------
# Example 2: Shear Stress
# ----------------------------------------------------
# Compute shear stress using:  τ = V * Q / (I * b)
V = 1000        # Shear force in N
Q = 2e-4        # Statical moment in m^3
b = 0.1         # Width in m
tau = shear_stress(V, Q, I, b)
print("Example 2: Shear Stress")
print(f"Inputs: V = {V} N, Q = {Q} m^3, I = {I} m^4, b = {b} m")
print(f"Shear stress, τ = {tau:.3f} Pa\n")

# ----------------------------------------------------
# Example 3: Statical Moment
# ----------------------------------------------------
# Compute statical moment using:  Q = (X * A_i * y_bar) / aa
X = 10          # Multiplying factor
A_i = 0.005     # Area in m^2
y_bar = 0.02    # Centroid location in m
aa = 0.03       # Reference distance in m
Q_val = statical_moment(X, A_i, y_bar, aa)
print("Example 3: Statical Moment")
print(f"Inputs: X = {X}, A_i = {A_i} m^2, y_bar = {y_bar} m, aa = {aa} m")
print(f"Statical moment, Q = {Q_val:.5f} m^3\n")

# ----------------------------------------------------
# Example 4: Parallel Axis Inertia
# ----------------------------------------------------
# Compute moment of inertia about a parallel axis using: I = I_na + A*d²
I_na = 1e-6     # Centroidal moment of inertia in m^4
A = 0.01        # Area in m^2
d = 0.05        # Distance in m
I_parallel = parallel_axis_inertia(I_na, A, d)
print("Example 4: Parallel Axis Inertia")
print(f"Inputs: I_na = {I_na} m^4, A = {A} m^2, d = {d} m")
print(f"Parallel axis inertia, I = {I_parallel:.6e} m^4\n")

# ----------------------------------------------------
# Example 5: Shear Force from Moment
# ----------------------------------------------------
# Let M_expr = 100 * x, so that V = dM/dx = 100
M_expr = 100 * x
V_expr = shear_force_from_moment(M_expr, x)
print("Example 5: Shear Force from Moment")
print(f"Input moment expression: M(x) = 100*x")
print(f"Computed shear force, V = {V_expr}\n")

# ----------------------------------------------------
# Example 6: Loading from Shear
# ----------------------------------------------------
# Let V_expr = 200 - 50*x, so that q = -dV/dx = 50
V_expr = 200 - 50 * x
q_expr = loading_from_shear(V_expr, x)
print("Example 6: Loading from Shear")
print(f"Input shear force expression: V(x) = 200 - 50*x")
print(f"Computed load, q = {q_expr}\n")

# ----------------------------------------------------
# Example 7: Strain Energy in Bending
# ----------------------------------------------------
# Using M_expr = 100 * x, E = 210e9, I = 8e-6, L = 2 m.
E = 210e9     # Modulus of elasticity in Pa
L = 2         # Beam length in m
U_bending = strain_energy_bending(M_expr, E, I, L, x)
print("Example 7: Strain Energy in Bending")
print(f"Inputs: M(x) = 100*x, E = {E} Pa, I = {I} m^4, L = {L} m")
print(f"Strain energy, U = {U_bending}\n")

# ----------------------------------------------------
# Example 8: Deflection using Castigliano's Theorem
# ----------------------------------------------------
# For demonstration, let dM/dP = 1 (constant)
dM_dP_expr = 1
delta_cast = deflection_castigliano(M_expr, E, I, L, x, dM_dP_expr)
print("Example 8: Deflection using Castigliano's Theorem")
print(f"Computed deflection, δ = {delta_cast}\n")

# ----------------------------------------------------
# Example 9: Rotation using Castigliano's Theorem
# ----------------------------------------------------
# For demonstration, let dM/dM = 1 (constant)
theta_cast = rotation_castigliano(M_expr, E, I, L, x, 1)
print("Example 9: Rotation using Castigliano's Theorem")
print(f"Computed rotation, θ = {theta_cast}\n")

# ----------------------------------------------------
# Example 10: Deflection from Loading (Symbolic)
# ----------------------------------------------------
# Let q(x) = 1000 (constant load)
q_expr = 1000
defl_results = deflection_from_loading(q_expr, E, I, x)
print("Example 10: Deflection from Loading")
print("Symbolic integration results:")
for i, expr in enumerate(defl_results, start=1):
    print(f"  Integration level {i}: {expr}")
print()

# ----------------------------------------------------
# Example 11: Discontinuity Function
# ----------------------------------------------------
# For x < 2, function = 0; for x >= 2, function = (x-2)^3.
disc_expr = discontinuity(x, 2, 3)
print("Example 11: Discontinuity Function")
print(f"Discontinuity expression: {disc_expr}\n")

# ----------------------------------------------------
# Example 12: Axial Deflection (Uniform)
# ----------------------------------------------------
N = 2000      # Axial load in N
axial_defl_uniform_val = axial_deflection_uniform(N, L, E, A)
print("Example 12: Axial Deflection (Uniform)")
print(f"Inputs: N = {N} N, L = {L} m, E = {E} Pa, A = {A} m^2")
print(f"Axial deflection, δ = {axial_defl_uniform_val:.6f} m\n")

# ----------------------------------------------------
# Example 13: Axial Deflection (Nonuniform)
# ----------------------------------------------------
# Let N(x) = 1000 * x
N_expr = 1000 * x
axial_defl_nonuni = axial_deflection_nonuniform(N_expr, E, A, L, x)
print("Example 13: Axial Deflection (Nonuniform)")
print(f"Computed axial deflection, δ = {axial_defl_nonuni}\n")

# ----------------------------------------------------
# Example 14: Axial Deflection (Multiple Segments)
# ----------------------------------------------------
# Define two segments: (N=2000, L=1) and (N=1500, L=1)
segments = [(2000, 1, E, A), (1500, 1, E, A)]
axial_defl_multiple_val = axial_deflection_multiple(segments)
print("Example 14: Axial Deflection (Multiple Segments)")
print(f"Computed axial deflection, δ = {axial_defl_multiple_val:.6f} m\n")

# ----------------------------------------------------
# Example 15: Axial Deflection due to Thermal Expansion
# ----------------------------------------------------
alpha = 1.2e-5  # Coefficient of thermal expansion (1/K)
deltaT = 30     # Temperature change (K)
axial_defl_thermal_val = axial_deflection_thermal(alpha, deltaT, L)
print("Example 15: Axial Deflection due to Thermal Expansion")
print(f"Inputs: α = {alpha}, ΔT = {deltaT} K, L = {L} m")
print(f"Thermal axial deflection, δ_thermal = {axial_defl_thermal_val:.6f} m\n")

# ----------------------------------------------------
# Example 16: Torsion Angle (Uniform)
# ----------------------------------------------------
T = 150       # Torque in Nm
G = 80e9      # Shear modulus in Pa
d_val = 0.1   # Diameter in m
J = polar_inertia_circular(d_val)
torsion_angle_uni = torsion_angle_uniform(T, L, G, J)
print("Example 16: Torsion Angle (Uniform)")
print(f"Inputs: T = {T} Nm, L = {L} m, G = {G} Pa, J = {J:.6e} m^4")
print(f"Torsion angle, ϕ = {torsion_angle_uni:.6f} radians\n")

# ----------------------------------------------------
# Example 17: Torsion Angle (Multiple Segments)
# ----------------------------------------------------
segments_torsion = [(150, 1, G, J), (150, 1, G, J)]
torsion_angle_mult = torsion_angle_multiple(segments_torsion)
print("Example 17: Torsion Angle (Multiple Segments)")
print(f"Computed torsion angle, ϕ = {torsion_angle_mult:.6f} radians\n")

# ----------------------------------------------------
# Example 18: Torsion Angle (Nonuniform)
# ----------------------------------------------------
# Let T(x) = 100 * x
T_expr = 100 * x
torsion_angle_nonuni = torsion_angle_nonuniform(T_expr, G, J, L, x)
print("Example 18: Torsion Angle (Nonuniform)")
print(f"Computed torsion angle, ϕ = {torsion_angle_nonuni}\n")

# ----------------------------------------------------
# Example 19: Axial Strain
# ----------------------------------------------------
delta = 0.001   # Change in length in m
L_orig = 2      # Original length in m
axial_str = axial_strain(delta, L_orig)
print("Example 19: Axial Strain")
print(f"Inputs: Δ = {delta} m, L = {L_orig} m")
print(f"Axial strain, ε = {axial_str:.6f}\n")

# ----------------------------------------------------
# Example 20: Thermal Strain
# ----------------------------------------------------
thermal_str = thermal_strain(alpha, deltaT)
print("Example 20: Thermal Strain")
print(f"Inputs: α = {alpha}, ΔT = {deltaT}")
print(f"Thermal strain, ε_thermal = {thermal_str:.6e}\n")

# ----------------------------------------------------
# Example 21: Poisson's Ratio
# ----------------------------------------------------
lat_strain = -0.0002
axial_strain_val = 0.001
nu_val = poisson_ratio(lat_strain, axial_strain_val)
print("Example 21: Poisson's Ratio")
print(f"Inputs: Lateral strain = {lat_strain}, Axial strain = {axial_strain_val}")
print(f"Computed Poisson's ratio, ν = {nu_val:.3f}\n")

# ----------------------------------------------------
# Example 22: Hooke's Law (Normal)
# ----------------------------------------------------
normal_stress_calc = hookes_law_normal(E, axial_strain_val)
print("Example 22: Hooke's Law (Normal)")
print(f"Inputs: E = {E} Pa, strain = {axial_strain_val}")
print(f"Normal stress via Hooke's law, σ = {normal_stress_calc:.3e} Pa\n")

# ----------------------------------------------------
# Example 23: Hooke's Law (Shear)
# ----------------------------------------------------
shear_strain_val = 0.005
shear_stress_calc = hookes_law_shear(G, shear_strain_val)
print("Example 23: Hooke's Law (Shear)")
print(f"Inputs: G = {G} Pa, shear strain = {shear_strain_val}")
print(f"Shear stress via Hooke's law, τ = {shear_stress_calc:.3e} Pa\n")

# ----------------------------------------------------
# Example 24: Axial Strain Energy (Uniform)
# ----------------------------------------------------
axial_strain_energy_uni = axial_strain_energy_uniform(N, L, E, A)
print("Example 24: Axial Strain Energy (Uniform)")
print(f"Inputs: N = {N} N, L = {L} m, E = {E} Pa, A = {A} m^2")
print(f"Axial strain energy, U = {axial_strain_energy_uni:.3f} J\n")

# ----------------------------------------------------
# Example 25: Axial Strain Energy (Multiple Segments)
# ----------------------------------------------------
axial_strain_energy_mult = axial_strain_energy_multiple(segments)
print("Example 25: Axial Strain Energy (Multiple Segments)")
print(f"Computed axial strain energy, U = {axial_strain_energy_mult:.3f} J\n")

# ----------------------------------------------------
# Example 26: Axial Strain Energy (Nonuniform)
# ----------------------------------------------------
axial_strain_energy_nonuni = axial_strain_energy_nonuniform(N_expr, E, A, L, x)
print("Example 26: Axial Strain Energy (Nonuniform)")
print(f"Computed axial strain energy, U = {axial_strain_energy_nonuni}\n")

# ----------------------------------------------------
# Example 27: Strain Energy Density
# ----------------------------------------------------
u_density = strain_energy_density(E, axial_strain_val)
print("Example 27: Strain Energy Density")
print(f"Inputs: E = {E} Pa, strain = {axial_strain_val}")
print(f"Strain energy density, u = {u_density:.3e} J/m³\n")

# ----------------------------------------------------
# Example 28: Torsion Shear Stress
# ----------------------------------------------------
rho = d_val / 2  # Evaluate at the surface (half the diameter)
torsion_shear = torsion_shear_stress(T, rho, J)
print("Example 28: Torsion Shear Stress")
print(f"Inputs: T = {T} Nm, ρ = {rho} m, J = {J:.6e} m^4")
print(f"Torsion shear stress, τ = {torsion_shear:.3e} Pa\n")

# ----------------------------------------------------
# Example 29: Rectangular Inertia
# ----------------------------------------------------
b_rect = 0.1   # width in m
h_rect = 0.2   # height in m
I_rect = rectangular_inertia(b_rect, h_rect, axis='x')
print("Example 29: Rectangular Inertia (x-axis)")
print(f"Inputs: b = {b_rect} m, h = {h_rect} m")
print(f"Rectangular inertia, I = {I_rect:.6e} m^4\n")

# ----------------------------------------------------
# Example 30: Circular Inertia
# ----------------------------------------------------
d_circ = 0.1  # Diameter in m
I_circ = circular_inertia(d_circ)
print("Example 30: Circular Inertia")
print(f"Input: d = {d_circ} m")
print(f"Circular inertia, I = {I_circ:.6e} m^4\n")

# ----------------------------------------------------
# Example 31: Polar Inertia Circular
# ----------------------------------------------------
J_circ = polar_inertia_circular(d_circ)
print("Example 31: Polar Inertia (Circular)")
print(f"Input: d = {d_circ} m")
print(f"Polar inertia, J = {J_circ:.6e} m^4\n")

# ----------------------------------------------------
# Example 32: Power Transmission
# ----------------------------------------------------
P_from_rpm = power_transmission(T, rpm=1800)
P_from_omega = power_transmission(T, omega_rad=10)
print("Example 32: Power Transmission")
print(f"Inputs: T = {T} Nm, rpm = 1800  -> P = {P_from_rpm:.3f} W")
print(f"Inputs: T = {T} Nm, ω = 10 rad/s -> P = {P_from_omega:.3f} W\n")

# ----------------------------------------------------
# Example 33: Torsion Strain Energy (Uniform)
# ----------------------------------------------------
torsion_energy_uni = torsion_strain_energy_uniform(T, L, G, J)
print("Example 33: Torsion Strain Energy (Uniform)")
print(f"Inputs: T = {T} Nm, L = {L} m, G = {G} Pa, J = {J:.6e} m^4")
print(f"Torsion strain energy, U = {torsion_energy_uni:.3f} J\n")

# ----------------------------------------------------
# Example 34: Torsion Strain Energy (Multiple Segments)
# ----------------------------------------------------
torsion_energy_mult = torsion_strain_energy_multiple(segments_torsion)
print("Example 34: Torsion Strain Energy (Multiple Segments)")
print(f"Computed torsion strain energy, U = {torsion_energy_mult:.3f} J\n")

# ----------------------------------------------------
# Example 35: Torsion Strain Energy (Nonuniform)
# ----------------------------------------------------
torsion_energy_nonuni = torsion_strain_energy_nonuniform(T_expr, G, J, L, x)
print("Example 35: Torsion Strain Energy (Nonuniform)")
print(f"Computed torsion strain energy, U = {torsion_energy_nonuni}\n")

# ----------------------------------------------------
# Example 36: Shear Strain Energy Density
# ----------------------------------------------------
u_shear_from_tau = shear_strain_energy_density(G, tau=shear_stress_calc)
u_shear_from_strain = shear_strain_energy_density(G, shear_strain=shear_strain_val)
print("Example 36: Shear Strain Energy Density")
print(f"From shear stress: u = {u_shear_from_tau:.3e} J/m³")
print(f"From shear strain: u = {u_shear_from_strain:.3e} J/m³\n")

# ----------------------------------------------------
# Example 37: Moment of Force Calculation
# ----------------------------------------------------
position = [1, 0, 0]
force = [0, 10, 0]
moment_vec = moment_of_force(position, force)
print("Example 37: Moment of Force Calculation")
print(f"Inputs: position = {position}, force = {force}")
print(f"Computed moment (torque): {moment_vec}")
