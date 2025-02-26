"""
Page 1-5 of ListOfEquations.pdf

Equations

  1. Bending Equations Functions
     normal_stress(M, y, I)
         Computes normal stress in bending.
         Equation: σ = -M*y / I

     shear_stress(V, Q, I, b)
         Computes shear stress.
         Equation: τ = V * Q / (I * b)

     statical_moment(X, A_i, y_bar, aa)
         Computes the statical moment Q.
         Equation: Q = (X * A_i * y_bar) / aa

     parallel_axis_inertia(I_na, A, d)
         Computes moment of inertia about a parallel axis.
         Equation: I = I_na + A*d²

     shear_force_from_moment(M_expr, x)
         Computes shear force V from the bending moment expression.
         Equation: V = dM/dx

     loading_from_shear(V_expr, x)
         Computes load q from the shear force expression.
         Equation: -q = dV/dx   (i.e. q = -dV/dx)

     strain_energy_bending(M_expr, E, I, L, x)
         Computes the strain energy in bending.
         Equation: U = ∫₀ᴸ [M(x)]²/(2*E*I) dx

     deflection_castigliano(M_expr, E, I, L, x, dM_dP_expr)
         Computes deflection using Castigliano's theorem.
         Equation: δ = ∫₀ᴸ (M/(E*I)) * (∂M/∂P) dx

     rotation_castigliano(M_expr, E, I, L, x, dM_dM_expr)
         Computes rotation using Castigliano's theorem.
         Equation: θ = ∫₀ᴸ (M/(E*I)) * (∂M/∂M) dx

  2. Deflection in Bending Functions
     deflection_from_loading(q_expr, E, I, x)
         Solves for deflection, slope, moment, and shear given the loading function.
         Governing equation: E*I*y''''(x) = -q(x)

  3. Discontinuity Function
     discontinuity(x_val, a, n)
         Returns the discontinuity function ⟨x - a⟩ⁿ:
           0, if x < a, and (x - a)ⁿ if x ≥ a.

  4. Deflection Under Axial Load Functions
     axial_deflection_uniform(N, L, E, A)
         Computes axial deflection for a uniform segment.
         Equation: δ = N*L / (E*A)

     axial_deflection_nonuniform(N_expr, E_expr, A_expr, L, x)
         Computes axial deflection for a non-uniform segment.
         Equation: δ = ∫₀ᴸ [N(x)]/(E(x)*A(x)) dx

     axial_deflection_multiple(segments)
         Computes axial deflection for multiple uniform segments.

     axial_deflection_thermal(alpha, deltaT, L)
         Computes thermal expansion deflection.
         Equation: δ_thermal = α * ΔT * L

  5. Deflection Under Torsion Functions
     torsion_angle_uniform(T, L, G, J)
         Computes the angle of twist for a uniform segment.
         Equation: ϕ = T * L / (G * J)

     torsion_angle_multiple(segments)
         Computes the angle of twist for multiple uniform segments.

     torsion_angle_nonuniform(T_expr, G_expr, J_expr, L, x)
         Computes the angle of twist for a non-uniform segment.
         Equation: ϕ = ∫₀ᴸ [T(x)]/(G(x)*J(x)) dx

  6. Useful Relations Functions
     axial_strain(delta, L)
         Computes axial strain.
         Equation: ε = δ / L

     thermal_strain(alpha, deltaT)
         Computes thermal strain.
         Equation: ε_thermal = α * ΔT

     poisson_ratio(lat_strain, axial_strain)
         Computes Poisson's ratio.
         Equation: ν = - (lateral strain / axial strain)

     hookes_law_normal(E, strain)
         Computes normal stress using Hooke's law.
         Equation: σ = E * strain

     hookes_law_shear(G, shear_strain)
         Computes shear stress using Hooke's law.
         Equation: τ = G * shear_strain

     axial_strain_energy_uniform(N, L, E, A)
         Computes axial strain energy for a uniform segment.
         Equation: U = N²*L/(2*E*A)

     axial_strain_energy_multiple(segments)
         Computes axial strain energy for multiple segments.

     axial_strain_energy_nonuniform(N_expr, E_expr, A_expr, L, x)
         Computes axial strain energy for a non-uniform segment.
         Equation: U = ∫₀ᴸ [N(x)]²/(2*E(x)*A(x)) dx

     strain_energy_density(E, strain, stress=None)
         Computes the strain energy density.
         (Uses stress²/(2*E) if stress is provided, otherwise E*strain²/2)

     torsion_shear_stress(T, rho, J)
         Computes shear stress in torsion.
         Equation: τ = T * ρ / J

     rectangular_inertia(b, h, axis='x')
         Computes the area moment of inertia for a rectangular cross-section.
         For axis 'x': I = b * h³/12; for axis 'y': I = h * b³/12

     circular_inertia(d)
         Computes the area moment of inertia for a circular cross-section.
         Equation: I = π * d⁴/64

     polar_inertia_circular(d)
         Computes the polar moment of inertia for a circular cross-section.
         Equation: J = π * d⁴/32

     power_transmission(T, rpm=None, omega_rad=None)
         Computes power transmission.
         Equation: P = T * (π*rpm/30) if rpm provided, or P = T * omega_rad

     torsion_strain_energy_uniform(T, L, G, J)
         Computes strain energy in torsion for a uniform segment.
         Equation: U = T²*L/(2*G*J)

     torsion_strain_energy_multiple(segments)
         Computes strain energy in torsion for multiple segments.

     torsion_strain_energy_nonuniform(T_expr, G_expr, J_expr, L, x)
         Computes strain energy in torsion for a non-uniform segment.
         Equation: U = ∫₀ᴸ [T(x)]²/(2*G(x)*J(x)) dx

     shear_strain_energy_density(G, shear_strain=None, tau=None)
         Computes shear strain energy density.
         (Uses τ²/(2*G) if shear stress is provided or G*(shear_strain)²/2 if shear strain is provided)

  7. Moment of Force Calculation
     moment_of_force(position, force)
         Computes the moment (torque) of a force about a point using the cross product.
"""

import math
import numpy as np
from sympy import symbols, diff, integrate, Piecewise

# ---------------------------
# 1. Bending Equations Functions
# ---------------------------

def normal_stress(M, y, I):
    """
    Computes the normal stress in bending.
    Equation: σ = -M*y / I.

    Parameters:
      M : Bending moment.
      y : Distance from the neutral axis.
      I : Moment of inertia.
    """
    return -M * y / I


def shear_stress(V, Q, I, b):
    """
    Computes the shear stress.
    Equation: τ = V * Q / (I * b).

    Parameters:
      V : Shear force.
      Q : Statical moment of area.
      I : Moment of inertia.
      b : Width of the section.
    """
    return V * Q / (I * b)


def statical_moment(X, A_i, y_bar, aa):
    """
    Computes the statical moment Q above/below a line.
    Equation: Q = (X * A_i * y_bar) / aa.

    Parameters:
      X    : A multiplying factor or moment arm.
      A_i  : Area of the segment.
      y_bar: Centroid location of the segment.
      aa   : Reference distance.
    """
    return X * A_i * y_bar / aa


def parallel_axis_inertia(I_na, A, d):
    """
    Computes the moment of inertia about a parallel axis.
    Equation: I = I_na + A*d^2.

    Parameters:
      I_na: Moment of inertia about the centroidal axis.
      A   : Area.
      d   : Distance between the centroidal axis and the new axis.
    """
    return I_na + A * d ** 2


def shear_force_from_moment(M_expr, x):
    """
    Computes the shear force V from the bending moment expression.
    Equation: V = dM/dx.

    Parameters:
      M_expr: A sympy expression for the bending moment as a function of x.
      x     : A sympy symbol.
    """
    return diff(M_expr, x)


def loading_from_shear(V_expr, x):
    """
    Computes the load q from the shear force expression.
    Equation: -q = dV/dx  →  q = -dV/dx.

    Parameters:
      V_expr: A sympy expression for the shear force as a function of x.
      x     : A sympy symbol.
    """
    return -diff(V_expr, x)


def strain_energy_bending(M_expr, E, I, L, x):
    """
    Computes the strain energy in bending.
    Equation: U = ∫₀ᴸ [M(x)]²/(2*E*I) dx.

    Parameters:
      M_expr: A sympy expression for the bending moment as a function of x.
      E     : Modulus of elasticity.
      I     : Moment of inertia.
      L     : Length of the beam.
      x     : A sympy symbol.
    """
    return integrate(M_expr ** 2 / (2 * E * I), (x, 0, L))


def deflection_castigliano(M_expr, E, I, L, x, dM_dP_expr):
    """
    Computes deflection using Castigliano's theorem.
    Equation: δ = ∫₀ᴸ (M/(E*I)) * (∂M/∂P) dx.

    Parameters:
      M_expr     : sympy expression for the bending moment.
      E          : Modulus of elasticity.
      I          : Moment of inertia.
      L          : Length of the beam.
      x          : A sympy symbol.
      dM_dP_expr : sympy expression for the derivative of M with respect to the applied load P.
    """
    return integrate(M_expr / (E * I) * dM_dP_expr, (x, 0, L))


def rotation_castigliano(M_expr, E, I, L, x, dM_dM_expr):
    """
    Computes rotation using Castigliano's theorem.
    Equation: θ = ∫₀ᴸ (M/(E*I)) * (∂M/∂M) dx.

    Parameters:
      M_expr     : sympy expression for the bending moment.
      E          : Modulus of elasticity.
      I          : Moment of inertia.
      L          : Length of the beam.
      x          : A sympy symbol.
      dM_dM_expr : sympy expression for the derivative of M with respect to a moment load.
    """
    return integrate(M_expr / (E * I) * dM_dM_expr, (x, 0, L))


# ---------------------------
# 2. Deflection in Bending Functions
# ---------------------------

def deflection_from_loading(q_expr, E, I, x):
    """
    Solves for deflection, slope, moment, and shear given the loading function.
    The governing differential equation is: E*I*y''''(x) = -q(x)

    Returns a tuple (y, y', y'', y''') where:
      y  : Deflection (with integration constant C4)
      y' : Angle of rotation (with integration constant C3)
      y'': Related to moment (with constant C2)
      y''': Related to shear (with constant C1)

    Note: Integration constants (C1, C2, C3, C4) appear in the expressions.

    Parameters:
      q_expr: A sympy expression for the load q(x).
      E     : Modulus of elasticity.
      I     : Moment of inertia.
      x     : A sympy symbol.
    """
    C1, C2, C3, C4 = symbols('C1 C2 C3 C4')
    y4 = -q_expr / (E * I)
    y3 = integrate(y4, x) + C1
    y2 = integrate(y3, x) + C2
    y1 = integrate(y2, x) + C3
    y = integrate(y1, x) + C4
    return y, y1, y2, y3


# ---------------------------
# 3. Discontinuity Function
# ---------------------------

def discontinuity(x_val, a, n):
    """
    Returns the discontinuity function ⟨x - a⟩^n defined as:
      0,         if x < a
      (x - a)^n, if x ≥ a

    Parameters:
      x_val: A sympy symbol or numeric value.
      a    : The threshold value.
      n    : The exponent.
    """
    return Piecewise((0, x_val < a), ((x_val - a) ** n, True))


# ---------------------------
# 4. Deflection Under Axial Load Functions
# ---------------------------

def axial_deflection_uniform(N, L, E, A):
    """
    Computes axial deflection for a uniform segment.
    Equation: δ = N*L / (E*A).

    Parameters:
      N : Axial load.
      L : Length.
      E : Modulus of elasticity.
      A : Cross-sectional area.
    """
    return N * L / (E * A)


def axial_deflection_nonuniform(N_expr, E_expr, A_expr, L, x):
    """
    Computes axial deflection for a non-uniform segment.
    Equation: δ = ∫₀ᴸ [N(x)]/(E(x)*A(x)) dx.

    Parameters:
      N_expr: sympy expression for the axial load as a function of x.
      E_expr: sympy expression or constant for the modulus of elasticity.
      A_expr: sympy expression or constant for the cross-sectional area.
      L     : Length.
      x     : A sympy symbol.
    """
    return integrate(N_expr / (E_expr * A_expr), (x, 0, L))


def axial_deflection_multiple(segments):
    """
    Computes axial deflection for multiple uniform segments.

    Parameters:
      segments: A list of tuples [(N, L, E, A), ...] for each segment.

    Returns:
      Sum of the deflections from each segment.
    """
    return sum(N * L / (E * A) for (N, L, E, A) in segments)


def axial_deflection_thermal(alpha, deltaT, L):
    """
    Computes the thermal expansion deflection.
    Equation: δ_thermal = α * ΔT * L.

    Parameters:
      alpha : Coefficient of thermal expansion.
      deltaT: Temperature change.
      L     : Length.
    """
    return alpha * deltaT * L


# ---------------------------
# 5. Deflection Under Torsion Functions
# ---------------------------

def torsion_angle_uniform(T, L, G, J):
    """
    Computes the angle of twist for a uniform segment.
    Equation: ϕ = T * L / (G * J).

    Parameters:
      T : Applied torque.
      L : Length.
      G : Shear modulus.
      J : Polar moment of inertia.
    """
    return T * L / (G * J)


def torsion_angle_multiple(segments):
    """
    Computes the angle of twist for multiple uniform segments.

    Parameters:
      segments: A list of tuples [(T, L, G, J), ...] for each segment.

    Returns:
      Sum of the twist angles from each segment.
    """
    return sum(T * L / (G * J) for (T, L, G, J) in segments)


def torsion_angle_nonuniform(T_expr, G_expr, J_expr, L, x):
    """
    Computes the angle of twist for a non-uniform segment.
    Equation: ϕ = ∫₀ᴸ [T(x)]/(G(x)*J(x)) dx.

    Parameters:
      T_expr: sympy expression for torque as a function of x.
      G_expr: sympy expression or constant for shear modulus.
      J_expr: sympy expression or constant for polar moment of inertia.
      L     : Length.
      x     : A sympy symbol.
    """
    return integrate(T_expr / (G_expr * J_expr), (x, 0, L))


# ---------------------------
# 6. Useful Relations Functions
# ---------------------------

def axial_strain(delta, L):
    """
    Computes axial strain.
    Equation: ε = δ / L.

    Parameters:
      delta: Change in length.
      L    : Original length.
    """
    return delta / L


def thermal_strain(alpha, deltaT):
    """
    Computes thermal strain.
    Equation: ε_thermal = α * ΔT.

    Parameters:
      alpha : Coefficient of thermal expansion.
      deltaT: Temperature change.
    """
    return alpha * deltaT


def poisson_ratio(lat_strain, axial_strain):
    """
    Computes Poisson's ratio.
    Equation: ν = - (lateral strain / axial strain).

    Parameters:
      lat_strain : Lateral strain.
      axial_strain: Axial strain.
    """
    return -lat_strain / axial_strain


def hookes_law_normal(E, strain):
    """
    Computes normal stress using Hooke's law.
    Equation: σ = E * strain.

    Parameters:
      E     : Modulus of elasticity.
      strain: Normal strain.
    """
    return E * strain


def hookes_law_shear(G, shear_strain):
    """
    Computes shear stress using Hooke's law.
    Equation: τ = G * shear_strain.

    Parameters:
      G          : Shear modulus.
      shear_strain: Shear strain.
    """
    return G * shear_strain


def axial_strain_energy_uniform(N, L, E, A):
    """
    Computes axial strain energy for a uniform segment.
    Equation: U = N² * L / (2 * E * A).

    Parameters:
      N : Axial load.
      L : Length.
      E : Modulus of elasticity.
      A : Cross-sectional area.
    """
    return N ** 2 * L / (2 * E * A)


def axial_strain_energy_multiple(segments):
    """
    Computes axial strain energy for multiple segments.

    Parameters:
      segments: A list of tuples [(N, L, E, A), ...] for each segment.

    Returns:
      Sum of the strain energies.
    """
    return sum(N ** 2 * L / (2 * E * A) for (N, L, E, A) in segments)


def axial_strain_energy_nonuniform(N_expr, E_expr, A_expr, L, x):
    """
    Computes axial strain energy for a non-uniform segment.
    Equation: U = ∫₀ᴸ [N(x)]²/(2*E(x)*A(x)) dx.

    Parameters:
      N_expr: sympy expression for axial load.
      E_expr: sympy expression or constant for modulus of elasticity.
      A_expr: sympy expression or constant for cross-sectional area.
      L     : Length.
      x     : A sympy symbol.
    """
    return integrate(N_expr ** 2 / (2 * E_expr * A_expr), (x, 0, L))


def strain_energy_density(E, strain, stress=None):
    """
    Computes the strain energy density.
    If stress is provided, uses: u = stress²/(2*E);
    otherwise, uses: u = E * strain² / 2.

    Parameters:
      E     : Modulus of elasticity.
      strain: Strain.
      stress: (Optional) Stress.
    """
    if stress is not None:
        return stress ** 2 / (2 * E)
    else:
        return E * strain ** 2 / 2


def torsion_shear_stress(T, rho, J):
    """
    Computes shear stress in torsion.
    Equation: τ = T * ρ / J.

    Parameters:
      T   : Applied torque.
      rho : Distance from the center (radius at which stress is evaluated).
      J   : Polar moment of inertia.
    """
    return T * rho / J


def rectangular_inertia(b, h, axis='x'):
    """
    Computes the area moment of inertia for a rectangular cross-section.

    For axis 'x' (base along x): I = b * h³/12.
    For axis 'y' (height along y): I = h * b³/12.

    Parameters:
      b   : Base width.
      h   : Height.
      axis: 'x' or 'y' (default is 'x').
    """
    if axis.lower() == 'x':
        return b * h ** 3 / 12
    elif axis.lower() == 'y':
        return h * b ** 3 / 12
    else:
        raise ValueError("Axis must be 'x' or 'y'.")


def circular_inertia(d):
    """
    Computes the area moment of inertia for a circular cross-section.
    Equation: I = π * d⁴/64.

    Parameters:
      d: Diameter.
    """
    return math.pi * d ** 4 / 64


def polar_inertia_circular(d):
    """
    Computes the polar moment of inertia for a circular cross-section.
    Equation: J = π * d⁴/32.

    Parameters:
      d: Diameter.
    """
    return math.pi * d ** 4 / 32


def power_transmission(T, rpm=None, omega_rad=None):
    """
    Computes power transmission.

    If rpm is provided: P = T * (π * rpm/30).
    If angular velocity (omega_rad) is provided: P = T * omega_rad.

    Parameters:
      T        : Torque.
      rpm      : Rotational speed in revolutions per minute.
      omega_rad: Angular velocity in radians per second.
    """
    if rpm is not None:
        return T * (math.pi * rpm / 30)
    elif omega_rad is not None:
        return T * omega_rad
    else:
        raise ValueError("Either rpm or omega_rad must be provided.")


def torsion_strain_energy_uniform(T, L, G, J):
    """
    Computes strain energy in torsion for a uniform segment.
    Equation: U = T² * L/(2*G*J).

    Parameters:
      T : Applied torque.
      L : Length.
      G : Shear modulus.
      J : Polar moment of inertia.
    """
    return T ** 2 * L / (2 * G * J)


def torsion_strain_energy_multiple(segments):
    """
    Computes strain energy in torsion for multiple segments.

    Parameters:
      segments: A list of tuples [(T, L, G, J), ...] for each segment.

    Returns:
      Sum of the strain energies.
    """
    return sum(T ** 2 * L / (2 * G * J) for (T, L, G, J) in segments)


def torsion_strain_energy_nonuniform(T_expr, G_expr, J_expr, L, x):
    """
    Computes strain energy in torsion for a non-uniform segment.
    Equation: U = ∫₀ᴸ [T(x)]²/(2*G(x)*J(x)) dx.

    Parameters:
      T_expr: sympy expression for torque.
      G_expr: sympy expression or constant for shear modulus.
      J_expr: sympy expression or constant for polar moment of inertia.
      L     : Length.
      x     : A sympy symbol.
    """
    return integrate(T_expr ** 2 / (2 * G_expr * J_expr), (x, 0, L))


def shear_strain_energy_density(G, shear_strain=None, tau=None):
    """
    Computes shear strain energy density.

    If shear stress τ is provided, uses: u = τ²/(2*G).
    Otherwise, if shear strain is provided, uses: u = G*(shear_strain)²/2.

    Parameters:
      G           : Shear modulus.
      shear_strain: (Optional) Shear strain.
      tau         : (Optional) Shear stress.
    """
    if tau is not None:
        return tau ** 2 / (2 * G)
    elif shear_strain is not None:
        return G * shear_strain ** 2 / 2
    else:
        raise ValueError("Either shear_strain or tau must be provided.")


# ---------------------------
# 7. Moment of Force Calculation
# ---------------------------

def moment_of_force(position, force):
    """
    Computes the moment (torque) of a force about a point using the cross product.

    Parameters:
      position: Array-like, position vector [x, y, z] from the reference point to the point of application.
      force   : Array-like, force vector [Fx, Fy, Fz].

    Returns:
      A NumPy array representing the moment vector.
    """
    position = np.array(position)
    force = np.array(force)
    return np.cross(position, force)

