# i got an hour this isn't going to be pretty
import numpy as np

F = 6000  # Force in Newtons (N)
Sy = 390 * 10**6  # Yield strength in Pascals (Pa)
Sut = 470e6  # Ultimate strength in Pascals (Pa)
E = 186e9  # Young's modulus in Pascals (Pa)
D = 50 * 10**-3 # Larger shaft diameter in meters, 50 mm
d = 25 * 10**-3  # Shaft diameter in meters, 25 mm
L = 0.5  # Shaft length in meters, 500 mm
r = 3 * 10**-3 # i think in mm?

print(f"D_d = {D/d}")
print(f"r_d = {r/d}")





# # import numpy as np
# #
# # # Converstions TODO: is this a bad idea?
# # MPA_to_Pa = 10**6
# # MM_to_M = 10**-3
# #
# # # Givens #########
# # # Forces
# # F = 6000  # Force in Newtons (N)
# #
# # # Material Properties (?)
# # Sy = 390 * 10**6  # Yield strength in Pascals (Pa)
# # Sut = 470e6  # Ultimate strength in Pascals (Pa)
# # E = 186e9  # Young's modulus in Pascals (Pa)
# #
# # # Geometry
# # D = 50 * 10**-3 # Larger shaft diameter in meters, 50 mm
# # d = 25 * 10**-3  # Shaft diameter in meters, 25 mm
# # L = 0.5  # Shaft length in meters, 500 mm
# # r = 3 * 10**-3 # i think in mm?
# # # Kf = 1.5  # Fatigue stress concentration factor (assumed)
# # # reliability = 0.99  # 99% reliability
# #
# # # TODO: for later
# # # # Step 1: Calculate the stress state (sigma_a and sigma_m)
# # # load_type = bending, torsion
# # # if load_type = bending:
# # #     if 50 <= Sut <= 250:
# # #         blah
# # #     if
# # #
# # # kf = 1 + (kt - 1)/(1 + np.sqrt(a/r))
# # #
# # # sigma_min = 1
# # # sigma_max = 1
# # #
# # # sigma_a = Kf * (sigma_max - sigma_min) / 2
# # # sigma_m = Kf * (sigma_max + sigma_min) / 2
# # #
# # #
# # #
# # # # M = F * L  # Moment (Nm)
# # # # I = (np.pi * d**4) / 64  # Moment of inertia
# # # # sigma_max = (M * (d / 2)) / I  # Max bending stress
# # # # sigma_min = 0  # Assuming fully reversed loading
# # # #
# # # # # Alternating and Mean Stresses
# # # # sigma_a = Kf * (sigma_max - sigma_min) / 2
# # # # sigma_m = Kf * (sigma_max + sigma_min) / 2
# # # #
# # # # # Step 2: Compute Endurance Limit
# # # # Se_prime = 0.5 * Sut  # Theoretical endurance limit
# # # # ka = 3.04 * (Sut / 1e6) ** -0.217  # Surface factor for machined/CD
# # # # kb = 1.24 * (d * 1e3) ** -0.107  # Size factor (converted to mm)
# # # # kc = 1.0  # Loading factor (bending)
# # # # kd = 1.0  # Temperature factor (assumed room temp)
# # # # ke = 1 - 0.08 * 2.326  # Reliability factor for 99% reliability
# # # #
# # # # Se = ka * kb * kc * kd * ke * Se_prime  # Adjusted endurance limit
# # # #
# # # # # Step 3: Compute Factor of Safety
# # # # nf = 1 / ((sigma_a / Se) + (sigma_m / Sut))
# # # # ny = Sy / (sigma_a + sigma_m)
# # # #
# # # # # Step 4: Compute Maximum Deflection
# # # # delta_max = (F * L**3) / (3 * E * I)  # Deflection formula for simply supported beam
# # # #
# # # # # Display results
# # # # print(f"Alternating Stress: {sigma_a / 1e6:.2f} MPa")
# # # # print(f"Mean Stress: {sigma_m / 1e6:.2f} MPa")
# # # # print(f"Endurance Limit: {Se / 1e6:.2f} MPa")
# # # # print(f"Fatigue Factor of Safety: {nf:.2f}")
# # # # print(f"Yielding Factor of Safety: {ny:.2f}")
# # # # print(f"Maximum Deflection: {delta_max * 1e3:.4f} mm")
# #
# # # TODO: ok stop here
# #
# # # import numpy as np
# # #
# # # # Given values
# # # F = 6000  # Force in Newtons (N)
# # # Sy = 390e6  # Yield strength in Pascals (Pa)
# # # Sut = 470e6  # Ultimate strength in Pascals (Pa)
# # # E = 186e9  # Young's modulus in Pascals (Pa)
# # # d = 25e-3  # Shaft diameter in meters (assumed 25 mm)
# # # L = 0.5  # Shaft length in meters (assumed 500 mm)
# # # Kf = 1.5  # Fatigue stress concentration factor (assumed)
# # # reliability = 0.99  # 99% reliability
# # #
# # # # Step 1: Compute maximum and minimum stress (assuming bending)
# # # M = F * L  # Moment (Nm)
# # # I = (np.pi * d**4) / 64  # Moment of inertia
# # # sigma_max = (M * (d / 2)) / I  # Max bending stress
# # # sigma_min = 0  # Assuming fully reversed loading
# # #
# # # # Alternating and Mean Stresses
# # # sigma_a = Kf * (sigma_max - sigma_min) / 2
# # # sigma_m = Kf * (sigma_max + sigma_min) / 2
# # #
# # # # Step 2: Compute Endurance Limit
# # # Se_prime = 0.5 * Sut  # Theoretical endurance limit
# # # ka = 3.04 * (Sut / 1e6) ** -0.217  # Surface factor for machined/CD
# # # kb = 1.24 * (d * 1e3) ** -0.107  # Size factor (converted to mm)
# # # kc = 1.0  # Loading factor (bending)
# # # kd = 1.0  # Temperature factor (assumed room temp)
# # # ke = 1 - 0.08 * 2.326  # Reliability factor for 99% reliability
# # #
# # # Se = ka * kb * kc * kd * ke * Se_prime  # Adjusted endurance limit
# # #
# # # # Step 3: Compute Factor of Safety
# # # nf = 1 / ((sigma_a / Se) + (sigma_m / Sut))
# # # ny = Sy / (sigma_a + sigma_m)
# # #
# # # # Step 4: Compute Maximum Deflection
# # # delta_max = (F * L**3) / (3 * E * I)  # Deflection formula for simply supported beam
# # #
# # # # Display results
# # # print(f"Alternating Stress: {sigma_a / 1e6:.2f} MPa")
# # # print(f"Mean Stress: {sigma_m / 1e6:.2f} MPa")
# # # print(f"Endurance Limit: {Se / 1e6:.2f} MPa")
# # # print(f"Fatigue Factor of Safety: {nf:.2f}")
# # # print(f"Yielding Factor of Safety: {ny:.2f}")
# # # print(f"Maximum Deflection: {delta_max * 1e3:.4f} mm")
#
# # TODO: next script
#
# # import math
# #
# # def calculate_stress_state(sigma_max, sigma_min, Kt, p, a, r):
# #     Kf = 1 + (Kt - 1) / (1 + p * (a / r))
# #     sigma_a = Kf * (sigma_max - sigma_min) / 2
# #     sigma_m = Kf * (sigma_max + sigma_min) / 2
# #     return sigma_a, sigma_m
# #
# # def determine_a(Sut, loading_type):
# #     if loading_type == "bending" or loading_type == "axial":
# #         if 50 <= Sut <= 250:
# #             return 0.246 - 3.08e-3 * Sut + 1.51e-5 * Sut**2 - 2.67e-8 * Sut**3
# #         elif 340 <= Sut <= 1700:
# #             return 1.24 - 2.25e-3 * Sut + 1.6e-6 * Sut**2 - 4.11e-10 * Sut**3
# #     elif loading_type == "torsion":
# #         if 50 <= Sut <= 220:
# #             return 0.19 - 2.51e-3 * Sut + 1.35e-5 * Sut**2 - 2.67e-8 * Sut**3
# #         elif 340 <= Sut <= 1500:
# #             return 0.958 - 1.83e-3 * Sut + 1.43e-6 * Sut**2 - 4.11e-10 * Sut**3
# #     return None
# #
# # def calculate_endurance_limit(Sut, ka, kb, kc, kd, ke):
# #     Se_prime = 0.5 * Sut if Sut <= 200 else 100 if Sut > 200 else 700
# #     return ka * kb * kc * kd * ke * Se_prime
# #
# # def calculate_surface_factor(Sut, surface_type):
# #     surface_factors = {
# #         "Ground": (1.21, 1.38, -0.067),
# #         "Machined/CD": (2.00, 3.04, -0.217),
# #         "Hot-Rolled": (11.0, 38.6, -0.650),
# #         "As-forged": (12.7, 54.9, -0.758)
# #     }
# #     a, _, b = surface_factors[surface_type]
# #     return a * Sut**b
# #
# # def calculate_size_factor(d, rotating=True):
# #     if rotating:
# #         if 0.3 <= d <= 2:
# #             return 0.879 * d**-0.107
# #         elif 2 <= d <= 10:
# #             return 0.91 * d**-0.157
# #     else:
# #         return 1  # For axial loading
# #     return None
# #
# # def calculate_load_factor(loading_type):
# #     return {"bending": 1, "axial": 0.85, "torsion": 0.59}.get(loading_type, None)
# #
# # def calculate_temperature_factor(T, unit="F"):
# #     if unit == "F":
# #         return 0.98 + 3.5e-4 * T - 6.3e-7 * T**2
# #     elif unit == "C":
# #         return 0.99 + 5.9e-4 * T - 2.1e-6 * T**2
# #     return None
# #
# # def calculate_reliability_factor(za):
# #     return 1 - 0.08 * za
# #
# # def calculate_factor_of_safety(sigma_a, sigma_m, Se, Sut, Sy):
# #     if sigma_a > 0 and sigma_m <= 0:
# #         ny = Sy / (sigma_a - sigma_m)
# #         nf = Se / sigma_a
# #     elif sigma_a > 0 and sigma_m > 0:
# #         ny = Sy / (sigma_a + sigma_m)
# #         nf = sigma_a / Se + sigma_m / (Sut - 1)
# #     elif sigma_a == 0 and sigma_m != 0:
# #         ny = Sy / abs(sigma_m)
# #         nf = None
# #     else:
# #         ny, nf = None, None
# #     return ny, nf
# #
# # def calculate_finite_life(nf, sigma_a, sigma_m, Sut, Se):
# #     if nf and nf < 1:
# #         if sigma_m == 0:
# #             a = (1.06 - 2.8e-3 * Sut + 6.9e-6 * Sut**2) * Sut**2 / Se
# #             b = -1 / 3 * math.log10((1.06 - 2.8e-3 * Sut + 6.9e-6 * Sut**2) * Sut / Se)
# #             N = (sigma_a / a) ** (1 / b)
# #         else:
# #             sigma_ar = sigma_a / (1 - sigma_m / Sut)
# #             a = (1.06 - 2.8e-3 * Sut + 6.9e-6 * Sut**2) * Sut**2 / Se
# #             b = -1 / 3 * math.log10((1.06 - 2.8e-3 * Sut + 6.9e-6 * Sut**2) * Sut / Se)
# #             N = (sigma_ar / a) ** (1 / b)
# #         return N
# #     return None
#
# # TODO: Next script
#
# """
# Solid Mechanics Factor of Safety and Life Calculation Script
# Author: Theresia Heimer
# Usage: Define input values in the `main()` function and run the script.
# """
#
# # TODO: it's different based on the units you use... like in (for kpsi) or mm (for MPA)
# def calculate_stress_state(sigma_max, sigma_min, Kt, p, a, r):
#     """Step 1: Calculate alternating and mean stress"""
#     Kf = 1 + (Kt - 1) / (1 + p * (a / r))
#     sigma_a = Kf * (sigma_max - sigma_min) / 2
#     sigma_m = Kf * (sigma_max + sigma_min) / 2
#     return sigma_a, sigma_m
#
#
# def calculate_surface_factor(Sut, surface_type):
#     """Step 2: Determine surface factor ka"""
#     surface_factors = {
#         "Ground": (1.21, 1.38, -0.067),
#         "Machined/CD": (2.00, 3.04, -0.217),
#         "Hot-Rolled": (11.0, 38.6, -0.650),
#         "As-forged": (12.7, 54.9, -0.758)
#     }
#     a, _, b = surface_factors[surface_type]
#     ka = a * (Sut ** b)
#     return ka
#
#
# def calculate_size_factor(d, rotating=True):
#     """Step 2: Determine size factor kb"""
#     if rotating:
#         return 0.879 * (d ** -0.107) if 0.3 <= d <= 2 else 0.91 * (d ** -0.157)
#     else:
#         return 1  # kb = 1 for axial loading
#
#
# def calculate_load_factor(load_type):
#     """Step 2: Determine load factor kc"""
#     return {"bending": 1, "axial": 0.85, "torsion": 0.59}.get(load_type, 1)
#
#
# def calculate_temperature_factor(T, unit="F"):
#     """Step 2: Determine temperature factor kd"""
#     if unit == "F":
#         return 0.98 + 3.5e-4 * T - 6.3e-7 * (T ** 2)
#     return 0.99 + 5.9e-4 * T - 2.1e-6 * (T ** 2)
#
#
# def calculate_reliability_factor(za):
#     """Step 2: Determine reliability factor ke"""
#     return 1 - 0.08 * za
#
#
# def calculate_endurance_limit(Sut, ka, kb, kc, kd, ke):
#     """Step 2: Calculate Endurance Limit Se"""
#     Se_prime = 0.5 * Sut if Sut <= 200 else 100  # kpsi condition
#     return ka * kb * kc * kd * ke * Se_prime
#
#
# def calculate_factor_of_safety(sigma_a, sigma_m, Se, Sut, Sy):
#     """Step 3: Calculate yield and fatigue factor of safety"""
#     if sigma_a > 0 and sigma_m <= 0:
#         ny = Sy / (sigma_a - sigma_m)
#         nf = Se / sigma_a
#     elif sigma_a > 0 and sigma_m > 0:
#         ny = Sy / (sigma_a + sigma_m)
#         nf = 1 / ((sigma_a / Se) + (sigma_m / (Sut - 1)))
#     else:
#         ny = Sy / abs(sigma_m)
#         nf = None
#     return ny, nf
#
#
# def calculate_finite_life(nf, sigma_a, sigma_m, Sut, Se):
#     """Step 4: Calculate finite life cycles N if nf < 1"""
#     f = 1.06 - 2.8e-3 * Sut + 6.9e-6 * (Sut ** 2) if 70 < Sut < 200 else 1.06 - 4.1e-4 * Sut + 1.5e-7 * (Sut ** 2)
#     a = (f * Sut) ** 2 / Se
#     b = -1 / 3 * (f * Sut / Se)
#     sigma_ar = sigma_a / (1 - sigma_m / Sut)
#     return (sigma_ar / a) ** (1 / b)
#
#
# def main():
#     """Main function to define inputs and run calculations."""
#     # Define input values here (Lines 79-89)
#     sigma_max = 300
#     sigma_min = 100
#     Kt = 1.5
#     p = 0.2
#     a = 0.5
#     r = 1.0
#     Sut = 150  # Ultimate tensile strength in kpsi
#     Sy = 250  # Yield strength in kpsi
#     d = 1.5  # Diameter in inches
#     load_type = "bending"
#     T = 100  # Temperature in Fahrenheit
#     surface_type = "Machined/CD"
#     za = 1.288  # 90% reliability
#
#     # Step 1: Calculate stress state
#     sigma_a, sigma_m = calculate_stress_state(sigma_max, sigma_min, Kt, p, a, r)
#
#     # Step 2: Calculate endurance limit
#     ka = calculate_surface_factor(Sut, surface_type)
#     kb = calculate_size_factor(d, rotating=True)
#     kc = calculate_load_factor(load_type)
#     kd = calculate_temperature_factor(T, unit="F")
#     ke = calculate_reliability_factor(za)
#     Se = calculate_endurance_limit(Sut, ka, kb, kc, kd, ke)
#
#     # Step 3: Calculate factor of safety
#     ny, nf = calculate_factor_of_safety(sigma_a, sigma_m, Se, Sut, Sy)
#
#     # Step 4: Calculate finite life if needed
#     N = calculate_finite_life(nf, sigma_a, sigma_m, Sut, Se) if nf and nf < 1 else "Infinite Life"
#
#     # Output results
#     print(f"Alternating Stress: {sigma_a:.2f} kpsi")
#     print(f"Mean Stress: {sigma_m:.2f} kpsi")
#     print(f"Endurance Limit: {Se:.2f} kpsi")
#     print(f"Yield Factor of Safety (ny): {ny:.2f}")
#     print(f"Fatigue Factor of Safety (nf): {nf:.2f}")
#     print(f"Estimated Life Cycles (N): {N}")
#
#
# if __name__ == "__main__":
#     main()
