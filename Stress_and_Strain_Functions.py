"""
Page 6-8 of ListOfEquations.pdf

Equations

  1. Stress Transformation Functions (Equations 1 - 9)
     stress_transformation_sigma(theta, sigma_x, sigma_y, tau_xy)
         Computes σθ, the normal stress on a plane rotated by angle θ.
         Equation: σθ = (σx + σy)/2 + ((σx - σy)/2)*cos(2θ) + τxy*sin(2θ)

     stress_transformation_sigma_perp(theta, sigma_x, sigma_y, tau_xy)
         Computes σθ+π/2, the normal stress on the plane perpendicular to the one at θ.
         Equation: σθ+π/2 = (σx + σy)/2 - ((σx - σy)/2)*cos(2θ) - τxy*sin(2θ)

     shear_stress_transformation(theta, sigma_x, sigma_y, tau_xy)
         Computes τθ, the shear stress on a plane rotated by angle θ.
         Equation: τθ = -((σx - σy)/2)*sin(2θ) + τxy*cos(2θ)

     principal_stresses(sigma_x, sigma_y, tau_xy)
         Computes the principal stresses σ1 and σ2.
         Equation: σ1,2 = σ_avg ± R, where σ_avg = (σx + σy)/2 and R = sqrt(((σx - σy)/2)^2 + τxy^2)

     max_inplane_shear_stress(sigma_x, sigma_y, tau_xy)
         Computes the maximum in-plane shear stress.
         Equation: τmax = R

     principal_stress_angle(sigma_x, sigma_y, tau_xy)
         Computes the principal stress angle.
         Equation: θp = 1/2 * arctan(2τxy/(σx - σy))

     max_shear_stress_angles(sigma_x, sigma_y, tau_xy)
         Computes the maximum shear stress angles.
         Equation: θs1 = θp - 45°, θs2 = θp + 45° (θp is the principal stress angle)

  2. Strain Transformation Functions (Equations 10 - 19)
     strain_transformation_epsilon(theta, epsilon_x, epsilon_y, gamma_xy)
         Computes εx' (or εθ), the normal strain on a plane rotated by angle θ.
         Equation: εx' = (εx + εy)/2 + ((εx - εy)/2)*cos(2θ) + (γxy/2)*sin(2θ)

     strain_transformation_epsilon_perp(theta, epsilon_x, epsilon_y, gamma_xy)
         Computes εy' (or εθ+π/2), the normal strain on the plane perpendicular to the one at θ.
         Equation: εy' = (εx + εy)/2 - ((εx - εy)/2)*cos(2θ) - (γxy/2)*sin(2θ)

     shear_strain_transformation(theta, epsilon_x, epsilon_y, gamma_xy)
         Computes the transformed shear strain.
         Equation: (γx'y'/2) = -((εx - εy)/2)*sin(2θ) + (γxy/2)*cos(2θ)
         (Returns the full shear strain, i.e. 2× the computed value)

     principal_strains(epsilon_x, epsilon_y, gamma_xy)
         Computes the principal strains ε1 and ε2.
         Equation: ε1,2 = εavg ± R, where εavg = (εx + εy)/2 and R = sqrt(((εx - εy)/2)^2 + (γxy/2)^2)

     max_inplane_shear_strain(epsilon_x, epsilon_y, gamma_xy)
         Computes the maximum in-plane shear strain.
         Equation: γmax = 2R

     principal_strain_angle(epsilon_x, epsilon_y, gamma_xy)
         Computes the principal strain angle.
         Equation: θp = 1/2 * arctan(γxy/(εx - εy))

     max_shear_strain_angles(epsilon_x, epsilon_y, gamma_xy)
         Computes the maximum shear strain angles.
         Equation: θs1 = θp - 45°, θs2 = θp + 45°

     strain_rosette(angles, strains)
         Converts strain rosette measurements (εa, εb, εc) to (εx, εy, γxy)
         (angles is a list/array of three angles in radians and strains is the corresponding measured values)

  3. Constitutive Relations (Equations 20 - 26)
     shear_modulus(E, nu)
         Computes the shear modulus G.
         Equation: G = E / [2*(1 + ν)]

     strain_from_stress(sigma_x, sigma_y, sigma_z, E, nu)
         Computes strain components (εx, εy, εz) from stress components.
         Equation: εx = (σx - ν*(σy + σz)) / E, etc.

     stress_from_strain(epsilon_x, epsilon_y, epsilon_z, E, nu)
         Computes stress components (σx, σy, σz) from strain components.
         Equation: σx = [E/((1+ν)(1-2ν))]*[(1-ν)εx + ν(εy+εz)], etc.

  4. Plane Stress & Plane Strain Functions (Equations 27 - 34)
     plane_stress_stress(epsilon_x, epsilon_y, gamma_xy, E, nu)
         Computes stress components under plane stress conditions.
         Equations:
           σx = E/(1-ν²)[εx + ν εy], σy = E/(1-ν²)[εy + ν εx],
           τxy = G*γxy (with G = E/(2*(1+ν))) and σz = 0.

     plane_strain_strain(sigma_x, sigma_y, tau_xy, E, nu)
         Computes strain components under plane strain conditions.
         Equations:
           εx = (1+ν)/E [σx(1-ν) - νσy], εy = (1+ν)/E [σy(1-ν) - νσx],
           γxy = τxy/G (with G = E/(2*(1+ν))) and εz = 0.

"""

import math
import numpy as np


# -------------------------------------------------
# 1. Stress Transformation Functions
# -------------------------------------------------
def stress_transformation_sigma(theta, sigma_x, sigma_y, tau_xy):
    """
    Computes σθ, the normal stress on a plane rotated by angle θ.
    Equation (1):
      σθ = (σx + σy)/2 + ((σx - σy)/2)*cos(2θ) + τxy*sin(2θ)
    """
    return (sigma_x + sigma_y) / 2 + ((sigma_x - sigma_y) / 2) * math.cos(2 * theta) + tau_xy * math.sin(2 * theta)


def stress_transformation_sigma_perp(theta, sigma_x, sigma_y, tau_xy):
    """
    Computes σθ+π/2, the normal stress on the plane perpendicular to the one at θ.
    Equation (2):
      σθ+π/2 = (σx + σy)/2 - ((σx - σy)/2)*cos(2θ) - τxy*sin(2θ)
    """
    return (sigma_x + sigma_y) / 2 - ((sigma_x - sigma_y) / 2) * math.cos(2 * theta) - tau_xy * math.sin(2 * theta)


def shear_stress_transformation(theta, sigma_x, sigma_y, tau_xy):
    """
    Computes τθ, the shear stress on a plane rotated by angle θ.
    Equation (3):
      τθ = -((σx - σy)/2)*sin(2θ) + τxy*cos(2θ)
    """
    return -((sigma_x - sigma_y) / 2) * math.sin(2 * theta) + tau_xy * math.cos(2 * theta)


def principal_stresses(sigma_x, sigma_y, tau_xy):
    """
    Computes the principal stresses σ1 and σ2.
    Equation (4):
      σ1,2 = σ_avg ± R, where σ_avg = (σx + σy)/2 and R = sqrt(((σx - σy)/2)^2 + τxy^2)
    """
    sigma_avg = (sigma_x + sigma_y) / 2
    R = math.sqrt(((sigma_x - sigma_y) / 2) ** 2 + tau_xy ** 2)
    return sigma_avg + R, sigma_avg - R


def max_inplane_shear_stress(sigma_x, sigma_y, tau_xy):
    """
    Computes the maximum in-plane shear stress.
    Equation (5):
      τmax = R, with R as defined above.
    """
    R = math.sqrt(((sigma_x - sigma_y) / 2) ** 2 + tau_xy ** 2)
    return R


def principal_stress_angle(sigma_x, sigma_y, tau_xy):
    """
    Computes the principal stress angle.
    Equation (6):
      θp = 1/2 * arctan(2τxy/(σx - σy))

    Note: Additional adjustments may be needed depending on whether σx > σy.
    """
    return 0.5 * math.atan2(2 * tau_xy, sigma_x - sigma_y)


def max_shear_stress_angles(sigma_x, sigma_y, tau_xy):
    """
    Computes the maximum shear stress angles.
    Equation (9):
      If θp is the principal stress angle, then:
         θs1 = θp - 45°   and   θs2 = θp + 45°.
      (Angles are in radians.)
    """
    theta_p = principal_stress_angle(sigma_x, sigma_y, tau_xy)
    return theta_p - math.pi / 4, theta_p + math.pi / 4


# -------------------------------------------------
# 2. Strain Transformation Functions
# -------------------------------------------------
def strain_transformation_epsilon(theta, epsilon_x, epsilon_y, gamma_xy):
    """
    Computes εx' (or εθ), the normal strain on a plane rotated by angle θ.
    Equation (10):
      εx' = (εx + εy)/2 + ((εx - εy)/2)*cos(2θ) + (γxy/2)*sin(2θ)
    """
    return (epsilon_x + epsilon_y) / 2 + ((epsilon_x - epsilon_y) / 2) * math.cos(2 * theta) + (
                gamma_xy / 2) * math.sin(2 * theta)


def strain_transformation_epsilon_perp(theta, epsilon_x, epsilon_y, gamma_xy):
    """
    Computes εy' (or εθ+π/2), the normal strain on the plane perpendicular to the one at θ.
    Equation (11):
      εy' = (εx + εy)/2 - ((εx - εy)/2)*cos(2θ) - (γxy/2)*sin(2θ)
    """
    return (epsilon_x + epsilon_y) / 2 - ((epsilon_x - epsilon_y) / 2) * math.cos(2 * theta) - (
                gamma_xy / 2) * math.sin(2 * theta)


def shear_strain_transformation(theta, epsilon_x, epsilon_y, gamma_xy):
    """
    Computes the transformed shear strain.
    Equation (12):
      (γx'y'/2) = -((εx - εy)/2)*sin(2θ) + (γxy/2)*cos(2θ)

    Returns the full shear strain γx'y' (i.e. 2*(value returned)).
    """
    return -((epsilon_x - epsilon_y) / 2) * math.sin(2 * theta) + (gamma_xy / 2) * math.cos(2 * theta)


def principal_strains(epsilon_x, epsilon_y, gamma_xy):
    """
    Computes the principal strains.
    Equation (13):
      ε1,2 = εavg ± R, where εavg = (εx + εy)/2 and R = sqrt(((εx - εy)/2)^2 + (γxy/2)^2)
    """
    eps_avg = (epsilon_x + epsilon_y) / 2
    R = math.sqrt(((epsilon_x - epsilon_y) / 2) ** 2 + (gamma_xy / 2) ** 2)
    return eps_avg + R, eps_avg - R


def max_inplane_shear_strain(epsilon_x, epsilon_y, gamma_xy):
    """
    Computes the maximum in-plane shear strain.
    Equation (14):
      γmax = 2R.
    """
    R = math.sqrt(((epsilon_x - epsilon_y) / 2) ** 2 + (gamma_xy / 2) ** 2)
    return 2 * R


def principal_strain_angle(epsilon_x, epsilon_y, gamma_xy):
    """
    Computes the principal strain angle.
    Equation (15):
      θp = 1/2 * arctan(γxy/(εx - εy))

    Note: Adjustments similar to the stress case might be needed.
    """
    return 0.5 * math.atan2(gamma_xy, epsilon_x - epsilon_y)


def max_shear_strain_angles(epsilon_x, epsilon_y, gamma_xy):
    """
    Computes the maximum shear strain angles.
    Equation (18):
      θs1 = θp - 45°   and   θs2 = θp + 45°.
    """
    theta_p = principal_strain_angle(epsilon_x, epsilon_y, gamma_xy)
    return theta_p - math.pi / 4, theta_p + math.pi / 4


def strain_rosette(angles, strains):
    """
    Solves for εx, εy, and γxy from strain rosette measurements.

    The strain rosette relation is given by:
        [1 + cos(2θ1)   1 - cos(2θ1)   2 sin(2θ1)]
        [1 + cos(2θ2)   1 - cos(2θ2)   2 sin(2θ2)]   * [εx, εy, γxy/2]^T = [2εa, 2εb, 2εc]^T
        [1 + cos(2θ3)   1 - cos(2θ3)   2 sin(2θ3)]

    Parameters:
      angles : list/array of three angles [θ1, θ2, θ3] in radians.
      strains: list/array of three measured strains [εa, εb, εc].

    Returns:
      (εx, εy, γxy) where γxy is doubled from the solved value.
    """
    theta1, theta2, theta3 = angles
    A = np.array([
        [1 + math.cos(2 * theta1), 1 - math.cos(2 * theta1), 2 * math.sin(2 * theta1)],
        [1 + math.cos(2 * theta2), 1 - math.cos(2 * theta2), 2 * math.sin(2 * theta2)],
        [1 + math.cos(2 * theta3), 1 - math.cos(2 * theta3), 2 * math.sin(2 * theta3)]
    ])
    b = np.array([2 * strains[0], 2 * strains[1], 2 * strains[2]])
    sol = np.linalg.solve(A, b)
    epsilon_x = sol[0]
    epsilon_y = sol[1]
    gamma_xy = 2 * sol[2]  # because we solved for γxy/2
    return epsilon_x, epsilon_y, gamma_xy


# -------------------------------------------------
# 3. Constitutive Relations
# -------------------------------------------------
def shear_modulus(E, nu):
    """
    Computes the shear modulus G from Young's modulus E and Poisson's ratio ν.
    G = E / [2*(1 + ν)]
    """
    return E / (2 * (1 + nu))


def strain_from_stress(sigma_x, sigma_y, sigma_z, E, nu):
    """
    Computes strain components (εx, εy, εz) from stress components using:
      εx = (σx - ν*(σy + σz)) / E, etc.
    """
    epsilon_x = (sigma_x - nu * (sigma_y + sigma_z)) / E
    epsilon_y = (sigma_y - nu * (sigma_x + sigma_z)) / E
    epsilon_z = (sigma_z - nu * (sigma_x + sigma_y)) / E
    return epsilon_x, epsilon_y, epsilon_z


def stress_from_strain(epsilon_x, epsilon_y, epsilon_z, E, nu):
    """
    Computes stress components (σx, σy, σz) from strain components using:
      σx = [E/((1+ν)(1-2ν))]*[(1-ν)εx + ν(εy+εz)], etc.
    """
    factor = E / ((1 + nu) * (1 - 2 * nu))
    sigma_x = factor * ((1 - nu) * epsilon_x + nu * (epsilon_y + epsilon_z))
    sigma_y = factor * ((1 - nu) * epsilon_y + nu * (epsilon_x + epsilon_z))
    sigma_z = factor * ((1 - nu) * epsilon_z + nu * (epsilon_x + epsilon_y))
    return sigma_x, sigma_y, sigma_z


# -------------------------------------------------
# 4. Plane Stress & Plane Strain Functions
# -------------------------------------------------
def plane_stress_stress(epsilon_x, epsilon_y, gamma_xy, E, nu):
    """
    Computes stress components for a plane stress condition.

    Equations:
      σx = E/(1-ν^2)[εx + ν εy]
      σy = E/(1-ν^2)[εy + ν εx]
      τxy = G * γxy,  with G = E/(2*(1+ν))
      σz = 0
    """
    sigma_x = E / (1 - nu ** 2) * (epsilon_x + nu * epsilon_y)
    sigma_y = E / (1 - nu ** 2) * (epsilon_y + nu * epsilon_x)
    G = shear_modulus(E, nu)
    tau_xy = G * gamma_xy
    sigma_z = 0
    return sigma_x, sigma_y, tau_xy, sigma_z


def plane_strain_strain(sigma_x, sigma_y, tau_xy, E, nu):
    """
    Computes strain components for a plane strain condition.

    Equations:
      εx = (1+ν)/E [σx(1-ν) - νσy]
      εy = (1+ν)/E [σy(1-ν) - νσx]
      γxy = τxy / G, with G = E/(2*(1+ν))
      εz = 0
    """
    epsilon_x = (1 + nu) / E * (sigma_x * (1 - nu) - nu * sigma_y)
    epsilon_y = (1 + nu) / E * (sigma_y * (1 - nu) - nu * sigma_x)
    G = shear_modulus(E, nu)
    gamma_xy = tau_xy / G
    epsilon_z = 0
    return epsilon_x, epsilon_y, gamma_xy, epsilon_z

