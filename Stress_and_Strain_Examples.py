"""
These are all the Stress and Strain Functions with example calculations

Each example prints the input values along with the computed result

Don't run them in here!! Copy and paste them into the calculations python file, or you'll run this whole thing!!
"""

from Stress_and_Strain_Functions import *
import math
import numpy as np
from sympy import symbols

x = symbols('x')

# ----------------------------------------------------
# Example 1: Stress Transformation - σθ
# ----------------------------------------------------
theta = 0.3         # radians
sigma_x = 100.0     # Pa
sigma_y = 50.0      # Pa
tau_xy = 25.0       # Pa
sigma_theta = stress_transformation_sigma(theta, sigma_x, sigma_y, tau_xy)
print("Example 1: Stress Transformation - σθ")
print(f"Inputs: θ = {theta} rad, σx = {sigma_x} Pa, σy = {sigma_y} Pa, τxy = {tau_xy} Pa")
print(f"Transformed normal stress, σθ = {sigma_theta:.3f} Pa\n")

# ----------------------------------------------------
# Example 2: Stress Transformation - σθ+π/2
# ----------------------------------------------------
sigma_theta_perp = stress_transformation_sigma_perp(theta, sigma_x, sigma_y, tau_xy)
print("Example 2: Stress Transformation - σθ+π/2")
print(f"Inputs: θ = {theta} rad, σx = {sigma_x} Pa, σy = {sigma_y} Pa, τxy = {tau_xy} Pa")
print(f"Transformed normal stress (perpendicular), σθ+π/2 = {sigma_theta_perp:.3f} Pa\n")

# ----------------------------------------------------
# Example 3: Shear Stress Transformation
# ----------------------------------------------------
tau_theta = shear_stress_transformation(theta, sigma_x, sigma_y, tau_xy)
print("Example 3: Shear Stress Transformation")
print(f"Inputs: θ = {theta} rad, σx = {sigma_x} Pa, σy = {sigma_y} Pa, τxy = {tau_xy} Pa")
print(f"Transformed shear stress, τθ = {tau_theta:.3f} Pa\n")

# ----------------------------------------------------
# Example 4: Principal Stresses
# ----------------------------------------------------
sigma1, sigma2 = principal_stresses(sigma_x, sigma_y, tau_xy)
print("Example 4: Principal Stresses")
print(f"Inputs: σx = {sigma_x} Pa, σy = {sigma_y} Pa, τxy = {tau_xy} Pa")
print(f"Principal stresses: σ1 = {sigma1:.3f} Pa, σ2 = {sigma2:.3f} Pa\n")

# ----------------------------------------------------
# Example 5: Maximum In-Plane Shear Stress
# ----------------------------------------------------
tau_max = max_inplane_shear_stress(sigma_x, sigma_y, tau_xy)
print("Example 5: Maximum In-Plane Shear Stress")
print(f"Inputs: σx = {sigma_x} Pa, σy = {sigma_y} Pa, τxy = {tau_xy} Pa")
print(f"Maximum in-plane shear stress, τmax = {tau_max:.3f} Pa\n")

# ----------------------------------------------------
# Example 6: Principal Stress Angle
# ----------------------------------------------------
theta_p = principal_stress_angle(sigma_x, sigma_y, tau_xy)
print("Example 6: Principal Stress Angle")
print(f"Inputs: σx = {sigma_x} Pa, σy = {sigma_y} Pa, τxy = {tau_xy} Pa")
print(f"Principal stress angle, θp = {theta_p:.3f} rad\n")

# ----------------------------------------------------
# Example 7: Maximum Shear Stress Angles
# ----------------------------------------------------
theta_s1, theta_s2 = max_shear_stress_angles(sigma_x, sigma_y, tau_xy)
print("Example 7: Maximum Shear Stress Angles")
print(f"Inputs: σx = {sigma_x} Pa, σy = {sigma_y} Pa, τxy = {tau_xy} Pa")
print(f"Maximum shear stress angles: θs1 = {theta_s1:.3f} rad, θs2 = {theta_s2:.3f} rad\n")

# ----------------------------------------------------
# Example 8: Strain Transformation - εx'
# ----------------------------------------------------
epsilon_x = 0.001
epsilon_y = 0.0005
gamma_xy = 0.002
epsilon_theta = strain_transformation_epsilon(theta, epsilon_x, epsilon_y, gamma_xy)
print("Example 8: Strain Transformation - εx'")
print(f"Inputs: θ = {theta} rad, εx = {epsilon_x}, εy = {epsilon_y}, γxy = {gamma_xy}")
print(f"Transformed normal strain, εx' = {epsilon_theta:.6f}\n")

# ----------------------------------------------------
# Example 9: Strain Transformation - εy'
# ----------------------------------------------------
epsilon_theta_perp = strain_transformation_epsilon_perp(theta, epsilon_x, epsilon_y, gamma_xy)
print("Example 9: Strain Transformation - εy'")
print(f"Inputs: θ = {theta} rad, εx = {epsilon_x}, εy = {epsilon_y}, γxy = {gamma_xy}")
print(f"Transformed normal strain (perpendicular), εy' = {epsilon_theta_perp:.6f}\n")

# ----------------------------------------------------
# Example 10: Shear Strain Transformation
# ----------------------------------------------------
gamma_trans = shear_strain_transformation(theta, epsilon_x, epsilon_y, gamma_xy)
print("Example 10: Shear Strain Transformation")
print(f"Inputs: θ = {theta} rad, εx = {epsilon_x}, εy = {epsilon_y}, γxy = {gamma_xy}")
# Returns half the shear strain; multiply by 2 for the full shear strain.
print(f"Full transformed shear strain, γx'y' = {2 * gamma_trans:.6f}\n")

# ----------------------------------------------------
# Example 11: Principal Strains
# ----------------------------------------------------
eps1, eps2 = principal_strains(epsilon_x, epsilon_y, gamma_xy)
print("Example 11: Principal Strains")
print(f"Inputs: εx = {epsilon_x}, εy = {epsilon_y}, γxy = {gamma_xy}")
print(f"Principal strains: ε1 = {eps1:.6f}, ε2 = {eps2:.6f}\n")

# ----------------------------------------------------
# Example 12: Maximum In-Plane Shear Strain
# ----------------------------------------------------
gamma_max = max_inplane_shear_strain(epsilon_x, epsilon_y, gamma_xy)
print("Example 12: Maximum In-Plane Shear Strain")
print(f"Inputs: εx = {epsilon_x}, εy = {epsilon_y}, γxy = {gamma_xy}")
print(f"Maximum in-plane shear strain, γmax = {gamma_max:.6f}\n")

# ----------------------------------------------------
# Example 13: Principal Strain Angle
# ----------------------------------------------------
theta_p_strain = principal_strain_angle(epsilon_x, epsilon_y, gamma_xy)
print("Example 13: Principal Strain Angle")
print(f"Inputs: εx = {epsilon_x}, εy = {epsilon_y}, γxy = {gamma_xy}")
print(f"Principal strain angle, θp = {theta_p_strain:.6f} rad\n")

# ----------------------------------------------------
# Example 14: Maximum Shear Strain Angles
# ----------------------------------------------------
theta_s1_strain, theta_s2_strain = max_shear_strain_angles(epsilon_x, epsilon_y, gamma_xy)
print("Example 14: Maximum Shear Strain Angles")
print(f"Inputs: εx = {epsilon_x}, εy = {epsilon_y}, γxy = {gamma_xy}")
print(f"Maximum shear strain angles: θs1 = {theta_s1_strain:.6f} rad, θs2 = {theta_s2_strain:.6f} rad\n")

# ----------------------------------------------------
# Example 15: Strain Rosette Conversion
# ----------------------------------------------------
angles = [0, math.pi / 4, math.pi / 2]  # Angles in radians
measured_strains = [0.001, 0.0012, 0.0008]
eps_x_rosette, eps_y_rosette, gamma_xy_rosette = strain_rosette(angles, measured_strains)
print("Example 15: Strain Rosette Conversion")
print(f"Inputs: angles = {angles}, measured strains = {measured_strains}")
print(f"Computed strains: εx = {eps_x_rosette:.6f}, εy = {eps_y_rosette:.6f}, γxy = {gamma_xy_rosette:.6f}\n")

# ----------------------------------------------------
# Example 16: Shear Modulus
# ----------------------------------------------------
E = 210e9  # Pa
nu = 0.3
G = shear_modulus(E, nu)
print("Example 16: Shear Modulus")
print(f"Inputs: E = {E} Pa, ν = {nu}")
print(f"Shear modulus, G = {G:.3e} Pa\n")

# ----------------------------------------------------
# Example 17: Strain from Stress
# ----------------------------------------------------
sigma_x_val = 100e6  # Pa
sigma_y_val = 50e6   # Pa
sigma_z_val = 20e6   # Pa
eps_from_stress = strain_from_stress(sigma_x_val, sigma_y_val, sigma_z_val, E, nu)
print("Example 17: Strain from Stress")
print(f"Inputs: σx = {sigma_x_val:.3e} Pa, σy = {sigma_y_val:.3e} Pa, σz = {sigma_z_val:.3e} Pa, E = {E} Pa, ν = {nu}")
print(f"Computed strains: εx = {eps_from_stress[0]:.6f}, εy = {eps_from_stress[1]:.6f}, εz = {eps_from_stress[2]:.6f}\n")

# ----------------------------------------------------
# Example 18: Stress from Strain
# ----------------------------------------------------
epsilon_x_val = 0.001
epsilon_y_val = 0.0005
epsilon_z_val = 0.0002
stress_from_eps = stress_from_strain(epsilon_x_val, epsilon_y_val, epsilon_z_val, E, nu)
print("Example 18: Stress from Strain")
print(f"Inputs: εx = {epsilon_x_val}, εy = {epsilon_y_val}, εz = {epsilon_z_val}, E = {E} Pa, ν = {nu}")
print(f"Computed stresses: σx = {stress_from_eps[0]:.3e} Pa, σy = {stress_from_eps[1]:.3e} Pa, σz = {stress_from_eps[2]:.3e} Pa\n")

# ----------------------------------------------------
# Example 19: Plane Stress - Stress Components
# ----------------------------------------------------
epsilon_x_plane = 0.001
epsilon_y_plane = 0.0005
gamma_xy_plane = 0.002
stress_plane = plane_stress_stress(epsilon_x_plane, epsilon_y_plane, gamma_xy_plane, E, nu)
print("Example 19: Plane Stress - Stress Components")
print(f"Inputs: εx = {epsilon_x_plane}, εy = {epsilon_y_plane}, γxy = {gamma_xy_plane}, E = {E} Pa, ν = {nu}")
print(f"Computed stresses: σx = {stress_plane[0]:.3e} Pa, σy = {stress_plane[1]:.3e} Pa, τxy = {stress_plane[2]:.3e} Pa, σz = {stress_plane[3]:.3e} Pa\n")

# ----------------------------------------------------
# Example 20: Plane Strain - Strain Components
# ----------------------------------------------------
sigma_x_plane = 100e6  # Pa
sigma_y_plane = 50e6   # Pa
tau_xy_plane = 20e6    # Pa
strain_plane = plane_strain_strain(sigma_x_plane, sigma_y_plane, tau_xy_plane, E, nu)
print("Example 20: Plane Strain - Strain Components")
print(f"Inputs: σx = {sigma_x_plane:.3e} Pa, σy = {sigma_y_plane:.3e} Pa, τxy = {tau_xy_plane:.3e} Pa, E = {E} Pa, ν = {nu}")
print(f"Computed strains: εx = {strain_plane[0]:.6f}, εy = {strain_plane[1]:.6f}, γxy = {strain_plane[2]:.6f}, εz = {strain_plane[3]:.6f}\n")
