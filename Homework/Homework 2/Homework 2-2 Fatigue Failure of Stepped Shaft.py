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
