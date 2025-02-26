import math

"""
This is Page 10-11 of ListOfEquations.pdf

Inputs:
    stress_unit: 'kpsi' or 'MPa', the units of stress
    loading_type: "bending", "axial", or "torsion", the type of loading
    sigma_max: the maximum stress (σmax)
    sigma_min: the minimum stress (σmin)
    Kt: the stress concentration factor
    r: the notch radius
    Sut: the ultimate tensile strength
    Sy: the yield strength
    bar_type: "Rotating", "Non-Rotating Circular", or "Non-Rotating Rectangular", the type of bar
    diameter: the diameter of the bar (if circular)
    h: the height of the bar (if rectangular)
    b: the width of the bar (if rectangular)
    surface_type: "Ground", "Machined/CD", "Hot-Rolled", or "As-forged", the surface condition
    
    *temp: the temperature value.
    *temp_unit: "F" or "C", the temperature unit.
    *reliability_percent: the desired reliability percentage (50, 90, 95, 99, 99.9, 99.99).
    
Note: professor has said sometimes temp and reliability won't be given, and to set that factor as 1
    
Outputs:
    a_notch: the calculated notch parameter 'a'
    Kf: the calculated fatigue stress concentration factor
    sigma_a: the calculated alternating stress (σa).
    sigma_m: the calculated mean stress (σm)
    Ka: surface factor
    Kb: size factor
    Kc: load factor
    Kd: temperature factor
    Ke: reliability factor
    Se: the calculated endurance limit
    ny: the calculated yield factor of safety 
    nf: the calculated fatigue factor of safety
    N: the calculated finite life (N) if nf < 1.
"""
import sys

# =======================
# Inputs
# =======================
print("This file is an input based one, type 'quit' if you ever want out immediately\n------ \n")

def input_or_quit(prompt):
    """
    Gets user input and exits if the user types 'quit'
    """
    user_input = input(prompt).strip()
    if user_input.lower() == "quit":
        sys.exit("User requested exit. Quitting")
    return user_input

# Stress and material units: 'kpsi' or 'MPa'
stress_unit = "kpsi"

# Loading type: "bending", "axial", or "torsion"
loading_type = "bending"

# Applies stress
sigma_max = 150.0  # Maximum stress (σmax)
sigma_min = 50.0  # Minimum stress (σmin)

# Stress concentration factors
Kt = 2.5  # Stress concentration factor
r = 0.1  # Radius (in the same length units as 'a'; inches for kpsi, mm for MPa)

# Material strengths
Sut = 200.0  # Ultimate tensile strength (Sut)
Sy = 120.0  # Yield strength (Sy)

# Distance
# For kpsi: d in inches
# For MPa: d in mm
# Round Rotating Bar: d = diameter
# Non-Rotating Bar:
#   d_e = 0.37d for circular section
#   d_e = 0.808 * sqrt(h*b) for rectangular
# Choose: "Rotating", "Not-Rotating Circular", or "Non-Rotating Rectangular"
bar_type = "Rotating"
diameter = 1.0
h = 1.0 # Height, you won't have this if round rotating bar
b = 1.0 # Diameter, you won't have this if round rotating bar

# Surface condition
# Choose: "Ground", "Machined/CD", "Hot-Rolled", or "As-forged"
surface_type = "Machined/CD"

# Temperature settings
# If you don't have a temperature specified, set these to None
# temp = None  # Not given
# temp_unit = None  # Not given
temp = 70          # Example: 70.0 for a specific temperature, or None if not given
temp_unit = "F"     # Example: "F" or "C", or None if not given


# Reliability setting (choose from 50, 90, 95, 99, 99.9, 99.99)
# If you don't have a temperature specified, set to None
reliability_percent = 95
# reliability_percent = None

# =======================
# End of Configuration
# =======================


def compute_a_notch(Sut, loading_type, stress_unit):
    """
    Compute the notch parameter 'a' based on loading type and stress unit.
    The formulas given are for √a; we square the result to obtain a.
    """
    loading_type = loading_type.lower()
    if loading_type in ["bending", "axial"]:
        if stress_unit.lower() == "kpsi":
            if 50 <= Sut <= 250:
                sqrt_a = 0.246 - 3.08e-3 * Sut + 1.51e-5 * (Sut ** 2) - 2.67e-8 * (Sut ** 3)
            else:
                raise ValueError("For bending/axial in kpsi, Sut must be between 50 and 250 kpsi.")
        elif stress_unit.lower() == "mpa":
            if 340 <= Sut <= 1700:
                sqrt_a = 1.24 - 2.25e-3 * Sut + 1.6e-6 * (Sut ** 2) - 4.11e-10 * (Sut ** 3)
            else:
                raise ValueError("For bending/axial in MPa, Sut must be between 340 and 1700 MPa.")
        else:
            raise ValueError("Unknown stress unit. Use 'kpsi' or 'MPa'.")
    elif loading_type == "torsion":
        if stress_unit.lower() == "kpsi":
            if 50 <= Sut <= 220:
                sqrt_a = 0.19 - 2.51e-3 * Sut + 1.35e-5 * (Sut ** 2) - 2.67e-8 * (Sut ** 3)
            else:
                raise ValueError("For torsion in kpsi, Sut must be between 50 and 220 kpsi.")
        elif stress_unit.lower() == "mpa":
            if 340 <= Sut <= 1500:
                sqrt_a = 0.958 - 1.83e-3 * Sut + 1.43e-6 * (Sut ** 2) - 4.11e-10 * (Sut ** 3)
            else:
                raise ValueError("For torsion in MPa, Sut must be between 340 and 1500 MPa.")
        else:
            raise ValueError("Unknown stress unit. Use 'kpsi' or 'MPa'.")
    else:
        raise ValueError("Unknown loading type. Use 'bending', 'axial', or 'torsion'.")

    a_val = sqrt_a ** 2
    return a_val


def compute_Kf(Kt, a_notch, r):
    """
    Calculate the fatigue stress concentration factor, Kf.
    Kf = 1 + (Kt - 1) / (1 +(a/r))
    """
    return 1 + ((Kt - 1) / (1 * (a_notch / r)))


def compute_surface_factor(surface_type, Sut, stress_unit):
    """
    Compute the surface factor, k_a.
    The coefficients depend on the chosen surface type and the stress unit.
    """
    surface_type = surface_type.lower()
    if stress_unit.lower() == "kpsi":
        if surface_type == "ground":
            a_coeff, b_exp = 1.21, -0.067
        elif surface_type in ["machined", "machined/cd", "cd"]:
            a_coeff, b_exp = 2.00, -0.217
        elif surface_type == "hot-rolled":
            a_coeff, b_exp = 11.0, -0.650
        elif surface_type == "as-forged":
            a_coeff, b_exp = 12.7, -0.758
        else:
            raise ValueError("Unknown surface type for kpsi. Choose Ground, Machined/CD, Hot-Rolled, or As-forged.")
    elif stress_unit.lower() == "mpa":
        if surface_type == "ground":
            a_coeff, b_exp = 1.38, -0.067
        elif surface_type in ["machined", "machined/cd", "cd"]:
            a_coeff, b_exp = 3.04, -0.217
        elif surface_type == "hot-rolled":
            a_coeff, b_exp = 38.6, -0.650
        elif surface_type == "as-forged":
            a_coeff, b_exp = 54.9, -0.758
        else:
            raise ValueError("Unknown surface type for MPa. Choose Ground, Machined/CD, Hot-Rolled, or As-forged.")
    else:
        raise ValueError("Unknown stress unit. Use 'kpsi' or 'MPa'.")

    k_a = a_coeff * (Sut ** b_exp)
    return k_a

if bar_type == "Rotating":
    d = diameter
elif bar_type == "Non-Rotating Circular":
    d = 0.37 * diameter
elif bar_type == "Non-Rotating Rectangular":
    d = 0.808 * math.sqrt(h * b)
else:
    raise ValueError("Unknown bar type, check spelling")

def compute_size_factor(d, stress_unit, loading_type):
    """
    Compute the size factor, k_b.
    For axial loading, k_b = 1.
    Otherwise, the formula depends on the diameter (d) and the units.
    """

    if loading_type.lower() == "axial":
        return 1.0
    if stress_unit.lower() == "kpsi":
        # d is in inches
        if 0.3 <= d <= 2:
            return 0.879 * (d ** -0.107)
        elif 2 < d <= 10:
            return 0.91 * (d ** -0.157)
        else:
            raise ValueError("For kpsi, d must be between 0.3 and 10 inches.")
    elif stress_unit.lower() == "mpa":
        # d is in mm
        if 7.62 <= d <= 51:
            return 1.24 * (d ** -0.107)
        elif 51 < d <= 254:
            return 1.51 * (d ** -0.157)
        else:
            raise ValueError("For MPa, d must be between 7.62 and 254 mm.")
    else:
        raise ValueError("Unknown stress unit. Use 'kpsi' or 'MPa'.")


def compute_load_factor(loading_type):
    """
    Compute the load factor, k_c.
      - k_c = 1.0 for bending,
      - k_c = 0.85 for axial,
      - k_c = 0.59 for torsion.
    """
    lt = loading_type.lower()
    if lt == "bending":
        return 1.0
    elif lt == "axial":
        return 0.85
    elif lt == "torsion":
        return 0.59
    else:
        raise ValueError("Unknown loading type for load factor.")


def compute_temperature_factor(temp=None, temp_unit=None):
    """
    Compute the temperature factor, k_d.
      - For Fahrenheit: k_d = 0.98 + 3.5e-4*T - 6.3e-7*T^2.
      - For Celsius: k_d = 0.99 + 5.9e-4*T - 2.1e-6*T^2.
    If temp or temp_unit is not provided (i.e. is None), return 1.
    """
    if temp is None or temp_unit is None:
        return 1.0  # No temperature effect
    if temp_unit.lower() == "f":
        return 0.98 + 3.5e-4 * temp - 6.3e-7 * (temp**2)
    elif temp_unit.lower() == "c":
        return 0.99 + 5.9e-4 * temp - 2.1e-6 * (temp**2)
    else:
        raise ValueError("Unknown temperature unit. Use 'F' or 'C'.")

def compute_reliability_factor(reliability_percent=None):
    """
    Compute the reliability factor, k_e.
    k_e = 1 - 0.08*z, where z depends on the desired percent reliability.
    The mapping is:
      50%   -> z = 0
      90%   -> z = 1.288
      95%   -> z = 1.645
      99%   -> z = 2.326
      99.9% -> z = 3.091
      99.99%-> z = 3.719
    If reliability_percent is not provided (None), return 1.
    """
    if reliability_percent is None:
        return 1.0  # No reliability effect
    mapping = {
        50: 0,
        90: 1.288,
        95: 1.645,
        99: 2.326,
        99.9: 3.091,
        99.99: 3.719
    }
    if reliability_percent in mapping:
        z_val = mapping[reliability_percent]
    else:
        # Choose the nearest key if an exact match isn't provided
        closest = min(mapping.keys(), key=lambda x: abs(x - reliability_percent))
        z_val = mapping[closest]
    return 1 - 0.08 * z_val



def compute_S_prime_e(Sut, stress_unit):
    """
    Compute S'_e.
      - For kpsi: S'_e = 0.5*Sut if Sut ≤ 200, otherwise 100.
      - For MPa:  S'_e = 0.5*Sut if Sut ≤ 1400, otherwise 700.
    """
    if stress_unit.lower() == "kpsi":
        return 0.5 * Sut if Sut <= 200 else 100.0
    elif stress_unit.lower() == "mpa":
        return 0.5 * Sut if Sut <= 1400 else 700.0
    else:
        raise ValueError("Unknown stress unit. Use 'kpsi' or 'MPa'.")


def compute_endurance_limit(Sut, d, surface_type, loading_type, temp, temp_unit, reliability_percent, stress_unit):
    """
    Calculate the endurance limit, Se.
    Se = k_a * k_b * k_c * k_d * k_e * S'_e
    """
    k_a = compute_surface_factor(surface_type, Sut, stress_unit)
    k_b = compute_size_factor(d, stress_unit, loading_type)
    k_c = compute_load_factor(loading_type)
    k_d = compute_temperature_factor(temp, temp_unit)  # If temp is None, k_d = 1
    k_e = compute_reliability_factor(reliability_percent)  # If reliability_percent is None, k_e = 1
    S_prime_e = compute_S_prime_e(Sut, stress_unit)
    Se = k_a * k_b * k_c * k_d * k_e * S_prime_e
    return Se


def compute_yield_FOS(Sy, sigma_a, sigma_m):
    """
    Calculate the yield factor of safety (ny).
      - If σa > 0 and σm ≤ 0 (compressive): ny = Sy / (σa - σm)
      - If σa > 0 and σm > 0 (tensile): ny = Sy / (σa + σm)
      - If σa = 0 and σm ≠ 0: ny = Sy / |σm|
    """
    if sigma_a > 0 and sigma_m <= 0:
        return Sy / (sigma_a - sigma_m)
    elif sigma_a > 0 and sigma_m > 0:
        return Sy / (sigma_a + sigma_m)
    elif sigma_a == 0 and sigma_m != 0:
        return Sy / abs(sigma_m)
    else:
        return None


def compute_fatigue_FOS(Se, sigma_a, sigma_m, Sut):
    """
    Calculate the fatigue factor of safety (nf).
      - For σa > 0, σm ≤ 0: nf = Se / σa
      - For σa > 0, σm > 0: nf = [ (σa/Se) + (σm/Sut) ]⁻¹
    """
    if sigma_a > 0 and sigma_m <= 0:
        return Se / sigma_a
    elif sigma_a > 0 and sigma_m > 0:
        return 1 / ((sigma_a / Se) + (sigma_m / Sut))
    else:
        return None


def compute_finite_life(sigma_a, sigma_m, Sut, Se, stress_unit):
    """
    Calculate the finite life (N) when nf < 1.
    First compute factor f based on Sut:
      - For kpsi: if 70 < Sut < 200, f = 1.06 - 2.8e-3*Sut + 6.9e-6*Sut².
      - For MPa:  if 500 < Sut < 1400, f = 1.06 - 4.1e-4*Sut + 1.5e-7*Sut².
    Then compute:
      a_fatigue = (f * Sut)² / Se
      b_fatigue = -1/3 * log(f * Sut / Se)
    Finally,
      - If σm = 0: N = (σa / a_fatigue)^(1/b_fatigue)
      - Otherwise: σ_ar = σa / (1 - σm/Sut) and N = (σ_ar / a_fatigue)^(1/b_fatigue)
    """
    if stress_unit.lower() == "kpsi":
        if 70 < Sut < 200:
            f = 1.06 - 2.8e-3 * Sut + 6.9e-6 * (Sut ** 2)
        else:
            raise ValueError("For finite life in kpsi, Sut must be between 70 and 200 kpsi.")
    elif stress_unit.lower() == "mpa":
        if 500 < Sut < 1400:
            f = 1.06 - 4.1e-4 * Sut + 1.5e-7 * (Sut ** 2)
        else:
            raise ValueError("For finite life in MPa, Sut must be between 500 and 1400 MPa.")
    else:
        raise ValueError("Unknown stress unit. Use 'kpsi' or 'MPa'.")

    a_fatigue = (f * Sut) ** 2 / Se
    ratio = f * Sut / Se
    if ratio <= 0:
        raise ValueError("Invalid ratio encountered in finite life calculation.")
    b_fatigue = -1 / 3 * math.log(ratio)

    if sigma_m == 0:
        N = (sigma_a / a_fatigue) ** (1 / b_fatigue)
    else:
        sigma_ar = sigma_a / (1 - sigma_m / Sut)
        N = (sigma_ar / a_fatigue) ** (1 / b_fatigue)
    return N


def main():
    # --- Step 1: Stress State Calculation ---
    a_notch = compute_a_notch(Sut, loading_type, stress_unit)
    print(f"Calculated notch parameter (a): {a_notch:.5f}")

    Kf = compute_Kf(Kt, a_notch, r)
    print(f"Calculated fatigue stress concentration factor (Kf): {Kf:.5f}")

    sigma_a = Kf * (sigma_max - sigma_min) / 2
    sigma_m = Kf * (sigma_max + sigma_min) / 2
    print(f"Calculated alternating stress (σa): {sigma_a:.5f}")
    print(f"Calculated mean stress (σm): {sigma_m:.5f}")

    # --- Step 2: Endurance Limit Calculation ---
    Se = compute_endurance_limit(Sut, d, surface_type, loading_type, temp, temp_unit, reliability_percent, stress_unit)
    print(f"Calculated endurance limit (Se): {Se:.5f}")

    # --- Step 3: Factor of Safety ---
    ny = compute_yield_FOS(Sy, sigma_a, sigma_m)
    nf = compute_fatigue_FOS(Se, sigma_a, sigma_m, Sut)
    if ny is not None:
        print(f"Yield factor of safety (ny): {ny:.5f}")
    else:
        print("Yield factor of safety (ny) could not be calculated with the given inputs.")
    if nf is not None:
        print(f"Fatigue factor of safety (nf): {nf:.5f}")
    else:
        print("Fatigue factor of safety (nf) could not be calculated with the given inputs.")

    # --- Step 4: Finite Life Calculation (if nf < 1) ---
    if nf is not None and nf < 1:
        N = compute_finite_life(sigma_a, sigma_m, Sut, Se, stress_unit)
        print(f"Finite life (N): {N:.5f}")
    else:
        print("Finite life calculation not required (nf ≥ 1).")


if __name__ == "__main__":
    main()
