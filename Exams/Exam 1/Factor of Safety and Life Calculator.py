import math

# -------------------------------
# Step 0: Define helper functions
# -------------------------------

def calculate_Kf(Kt, p, a_notch, r):
    """
    Notch stress concentration factor corrected for notch sensitivity.
    Equation: Kf = 1 + (Kt - 1) / (1 + p*(a_notch/r))

    Parameters:
      Kt      : The theoretical stress concentration factor.
      p       : Notch sensitivity parameter.
      a_notch : Characteristic flaw size (in same units as r).
      r       : Notch radius.

    Returns:
      Kf      : Effective stress concentration factor.
    """
    return 1 + (Kt - 1) / (1 + p * (a_notch / r))


def calculate_stress_state(sigma_max, sigma_min, Kf):
    """
    Computes the alternating and mean stresses.

    Equations:
       σa = Kf*(σmax - σmin)/2
       σm = Kf*(σmax + σmin)/2

    Parameters:
      sigma_max : Maximum stress.
      sigma_min : Minimum stress.
      Kf        : Effective stress concentration factor.

    Returns:
      sigma_a : Alternating stress.
      sigma_m : Mean stress.
    """
    sigma_a = Kf * (sigma_max - sigma_min) / 2
    sigma_m = Kf * (sigma_max + sigma_min) / 2
    return sigma_a, sigma_m


def compute_crack_size(Sut, loading_type='bending', units='kpsi'):
    """
    Computes the characteristic crack (or flaw) size using the empirical formula
    for √a. (Remember: a = (sqrt_a)^2.)

    For bending/axial loading:
      - kpsi:  √a = 0.246 - 3.08e-3 * Sut + 1.51e-5 * Sut**2 - 2.67e-8 * Sut**3,
               valid for 50 ≤ Sut ≤ 250 kpsi.
      - MPa:   √a = 1.24 - 2.25e-3 * Sut + 1.6e-6 * Sut**2 - 4.11e-10 * Sut**3,
               valid for 340 ≤ Sut ≤ 1700 MPa.

    For torsion:
      - kpsi:  √a = 0.19 - 2.51e-3 * Sut + 1.35e-5 * Sut**2 - 2.67e-8 * Sut**3,
               valid for 50 ≤ Sut ≤ 220 kpsi.
      - MPa:   √a = 0.958 - 1.83e-3 * Sut + 1.43e-6 * Sut**2 - 4.11e-10 * Sut**3,
               valid for 340 ≤ Sut ≤ 1500 MPa.

    Parameters:
      Sut         : Ultimate tensile strength.
      loading_type: 'bending' (or 'axial') or 'torsion'
      units       : 'kpsi' or 'MPa'

    Returns:
      a_crack: The computed flaw size (a) (inches or mm depending on units).
    """
    if units == 'kpsi':
        if loading_type in ['bending', 'axial']:
            sqrt_a = 0.246 - 3.08e-3 * Sut + 1.51e-5 * Sut ** 2 - 2.67e-8 * Sut ** 3
        elif loading_type == 'torsion':
            sqrt_a = 0.19 - 2.51e-3 * Sut + 1.35e-5 * Sut ** 2 - 2.67e-8 * Sut ** 3
        else:
            raise ValueError("Invalid loading type.")
    elif units == 'MPa':
        if loading_type in ['bending', 'axial']:
            sqrt_a = 1.24 - 2.25e-3 * Sut + 1.6e-6 * Sut ** 2 - 4.11e-10 * Sut ** 3
        elif loading_type == 'torsion':
            sqrt_a = 0.958 - 1.83e-3 * Sut + 1.43e-6 * Sut ** 2 - 4.11e-10 * Sut ** 3
        else:
            raise ValueError("Invalid loading type.")
    else:
        raise ValueError("Invalid unit system.")

    a_crack = sqrt_a ** 2
    return a_crack


def calculate_endurance_limit(Sut, d, loading_type='bending', surface_type='Machined/CD',
                              temperature=70, temp_unit='F', reliability=0.90, units='kpsi'):
    """
    Calculates the endurance limit (Se) based on surface, size, load, temperature, and reliability factors.

    The endurance limit is given by:
      Se = ka * kb * kc * kd * ke * S'_e

    where:
      • ka = constant * Sut^(exponent), with values from the surface type table.
      • kb = size factor (for non-axial loading, depends on diameter d).
      • kc = load factor (1 for bending, 0.85 for axial, 0.59 for torsion).
      • kd = temperature factor.
      • ke = reliability factor = 1 - 0.08*z_a (z_a based on reliability level).
      • S'_e = 0.5*Sut if Sut ≤ limit, otherwise a fixed value (100 kpsi or 700 MPa).

    Parameters:
      Sut         : Ultimate tensile strength.
      d           : Diameter (in inches for kpsi or mm for MPa).
      loading_type: 'bending', 'axial', or 'torsion'.
      surface_type: One of 'Ground', 'Machined/CD', 'Hot-Rolled', or 'As-forged'.
      temperature : Temperature value.
      temp_unit   : 'F' for Fahrenheit or 'C' for Celsius.
      reliability : Reliability as a fraction (e.g., 0.90).
      units       : 'kpsi' or 'MPa'.

    Returns:
      Se : The endurance limit.
    """
    # Surface factor (ka) from table:
    surface_factors = {
        'Ground': {'k': 1.21 if units == 'kpsi' else 1.38, 'b': -0.067},
        'Machined/CD': {'k': 2.00 if units == 'kpsi' else 3.04, 'b': -0.217},
        'Hot-Rolled': {'k': 11.0 if units == 'kpsi' else 38.6, 'b': -0.650},
        'As-forged': {'k': 12.7 if units == 'kpsi' else 54.9, 'b': -0.758}
    }
    if surface_type not in surface_factors:
        raise ValueError("Invalid surface type.")
    k_val = surface_factors[surface_type]['k']
    b_exp = surface_factors[surface_type]['b']
    ka = k_val * (Sut ** b_exp)

    # Size factor (kb)
    # For axial loading, kb = 1.
    if loading_type == 'axial':
        kb = 1.0
    else:
        if units == 'kpsi':
            if 0.3 <= d <= 2:
                kb = 0.879 * (d ** -0.107)
            elif 2 < d <= 10:
                kb = 0.91 * (d ** -0.157)
            else:
                raise ValueError("Diameter out of range for kpsi size factor.")
        elif units == 'MPa':
            if 7.62 <= d <= 51:
                kb = 1.24 * (d ** -0.107)
            elif 51 < d <= 254:
                kb = 1.51 * (d ** -0.157)
            else:
                raise ValueError("Diameter out of range for MPa size factor.")
        else:
            raise ValueError("Invalid units for size factor.")

    # Load factor (kc)
    if loading_type == 'bending':
        kc = 1.0
    elif loading_type == 'axial':
        kc = 0.85
    elif loading_type == 'torsion':
        kc = 0.59
    else:
        raise ValueError("Invalid loading type for load factor.")

    # Temperature factor (kd)
    if temp_unit.upper() == 'F':
        kd = 0.98 + 3.5e-4 * temperature - 6.3e-7 * (temperature ** 2)
    elif temp_unit.upper() == 'C':
        kd = 0.99 + 5.9e-4 * temperature - 2.1e-6 * (temperature ** 2)
    else:
        raise ValueError("Invalid temperature unit.")

    # Reliability factor (ke)
    reliability_map = {
        0.50: 0.0,
        0.90: 1.288,
        0.95: 1.645,
        0.99: 2.326,
        0.999: 3.091,
        0.9999: 3.719
    }
    if reliability in reliability_map:
        z_a = reliability_map[reliability]
    else:
        raise ValueError("Unsupported reliability value. Use one of: 0.50, 0.90, 0.95, 0.99, 0.999, 0.9999.")
    ke = 1 - 0.08 * z_a

    # Unmodified endurance limit, S'_e
    if units == 'kpsi':
        if Sut <= 200:
            Se_prime = 0.5 * Sut
        else:
            Se_prime = 100.0
    elif units == 'MPa':
        if Sut <= 1400:
            Se_prime = 0.5 * Sut
        else:
            Se_prime = 700.0
    else:
        raise ValueError("Invalid units for endurance limit.")

    Se = ka * kb * kc * kd * ke * Se_prime
    return Se


def calculate_f(Sut, units='kpsi'):
    """
    Calculates the fatigue factor f.

    For 70 < Sut < 200 kpsi:
         f = 1.06 - 2.8e-3 * Sut + 6.9e-6 * Sut**2
    For 500 < Sut < 1400 MPa:
         f = 1.06 - 4.1e-4 * Sut + 1.5e-7 * Sut**2

    Parameters:
      Sut   : Ultimate tensile strength.
      units : 'kpsi' or 'MPa'

    Returns:
      f_val: Fatigue factor.
    """
    if units == 'kpsi':
        f_val = 1.06 - 2.8e-3 * Sut + 6.9e-6 * Sut ** 2
    elif units == 'MPa':
        f_val = 1.06 - 4.1e-4 * Sut + 1.5e-7 * Sut ** 2
    else:
        raise ValueError("Invalid unit system.")
    return f_val


def calculate_factor_of_safety(sigma_a, sigma_m, Sy, Se, Sut):
    """
    Calculates the factors of safety. Different equations apply depending on the stress state:

    Case 1: σa > 0, σm ≤ 0 (compressive)
         Yield FOS: ny = Sy / (σa - σm)
         Fatigue FOS: nf = Se / σa

    Case 2: σa > 0, σm > 0 (tensile)
         Yield FOS: ny = Sy / (σa + σm)
         Fatigue FOS: nf = (σa / Se) + (σm / (Sut - 1))   [*Check unit consistency*]

    Case 3: σa = 0, σm ≠ 0 (static)
         Yield FOS: ny = Sy / |σm|
         (Fatigue factor may not be applicable)

    Parameters:
      sigma_a : Alternating stress.
      sigma_m : Mean stress.
      Sy      : Yield strength.
      Se      : Endurance limit.
      Sut     : Ultimate tensile strength.

    Returns:
      ny : Yielding factor of safety.
      nf : Fatigue factor of safety (if applicable; else None).
    """
    if sigma_a > 0 and sigma_m <= 0:
        ny = Sy / (sigma_a - sigma_m)
        nf = Se / sigma_a
    elif sigma_a > 0 and sigma_m > 0:
        ny = Sy / (sigma_a + sigma_m)
        # The following formulation is based on the provided process.
        nf = (sigma_a / Se) + (sigma_m / (Sut - 1))
    elif sigma_a == 0 and sigma_m != 0:
        ny = Sy / abs(sigma_m)
        nf = None
    else:
        ny = None
        nf = None
    return ny, nf


def calculate_finite_life(sigma_a, sigma_m, Sut, Se, f):
    """
    Calculates the finite fatigue life (N). Two cases are considered:

    (A) When σm = 0:
         a_fatigue = (f * Sut)**2 / Se
         b = -1/3 * log10( f*Sut / Se )
         N = (σa / a_fatigue)^(1/b)

    (B) When σm ≠ 0:
         σa,r = σa / (1 - σm/Sut)
         N = (σa,r / a_fatigue)^(1/b)

    Parameters:
      sigma_a : Alternating stress.
      sigma_m : Mean stress.
      Sut     : Ultimate tensile strength.
      Se      : Endurance limit.
      f       : Fatigue factor (from calculate_f).

    Returns:
      N       : Number of cycles to failure.
      a_fatigue: Parameter 'a' used in the life calculation.
      b       : Exponent parameter.
    """
    # Compute a_fatigue using the provided formula:
    a_fatigue = (f * Sut) ** 2 / Se
    ratio = (f * Sut) / Se
    if ratio <= 0:
        raise ValueError("Invalid ratio for life calculation.")
    # Using log base 10 as indicated in the process:
    b = - (1 / 3) * math.log10(ratio)

    if sigma_m == 0:
        N = (sigma_a / a_fatigue) ** (1 / b)
    else:
        sigma_ar = sigma_a / (1 - sigma_m / Sut)
        N = (sigma_ar / a_fatigue) ** (1 / b)
    return N, a_fatigue, b


# -------------------------------
# Main: Set initial conditions and run calculations
# -------------------------------
if __name__ == "__main__":
    # --- User-defined initial conditions ---
    # Units can be 'kpsi' or 'MPa'. (This example uses kpsi/in.)
    units = 'kpsi'

    # Material properties (example values)
    Sut = 150.0  # Ultimate tensile strength (kpsi)
    Sy = 90.0  # Yield strength (kpsi)

    # Loading conditions
    sigma_max = 50.0  # Maximum applied stress (kpsi)
    sigma_min = 10.0  # Minimum applied stress (kpsi)

    # Notch and geometric parameters for Kf:
    Kt = 2.5  # Theoretical stress concentration factor
    p = 0.8  # Notch sensitivity parameter
    # Compute a_notch from material properties and loading type (bending or torsion)
    loading_type = 'bending'
    a_notch = compute_crack_size(Sut, loading_type, units)  # characteristic flaw size
    r = 0.1  # Notch radius (inches)

    # Endurance limit factors:
    d = 1.0  # Diameter (inches for kpsi system)
    surface_type = 'Machined/CD'
    temperature = 70  # Temperature in Fahrenheit
    temp_unit = 'F'
    reliability = 0.90

    # -------------------------------
    # Step 1: Stress state
    # -------------------------------
    Kf = calculate_Kf(Kt, p, a_notch, r)
    sigma_a, sigma_m = calculate_stress_state(sigma_max, sigma_min, Kf)
    print("Effective Kf =", Kf)
    print("Alternating stress σa =", sigma_a)
    print("Mean stress σm =", sigma_m)

    # -------------------------------
    # Step 2: Endurance Limit
    # -------------------------------
    Se = calculate_endurance_limit(Sut, d, loading_type, surface_type, temperature, temp_unit, reliability, units)
    print("Endurance Limit Se =", Se)

    # -------------------------------
    # Step 3: Factor of Safety
    # -------------------------------
    ny, nf = calculate_factor_of_safety(sigma_a, sigma_m, Sy, Se, Sut)
    print("Yield Factor of Safety (ny) =", ny)
    if nf is not None:
        print("Fatigue Factor of Safety (nf) =", nf)

    # -------------------------------
    # Step 4: Finite Life Calculation
    # Only applicable if nf < 1 (finite life)
    # -------------------------------
    # Calculate fatigue factor f:
    f = calculate_f(Sut, units)
    # Compute life N:
    try:
        N, a_fatigue, b = calculate_finite_life(sigma_a, sigma_m, Sut, Se, f)
        print("Finite life (N cycles) =", N)
        print("Life calculation parameters: a =", a_fatigue, ", b =", b)
    except Exception as e:
        print("Error in life calculation:", e)

    # -------------------------------
    # Yielding Check (if desired)
    # -------------------------------
    # One might also check for yielding using similar equations.
    # For example, using: a_yield = (f * Sy)**2 / Se and b_yield = -1/3 * log10(f * Sy / Se)
    a_yield = (f * Sy) ** 2 / Se
    ratio_yield = (f * Sy) / Se
    if ratio_yield > 0:
        b_yield = - (1 / 3) * math.log10(ratio_yield)
        N_yield = (sigma_a / a_yield) ** (1 / b_yield)
        print("Yielding check - N =", N_yield)
    else:
        print("Invalid ratio for yielding check.")
