"""
This file is the main interface for performing calculations, mostly for ease of use (easier than 400 lines right)

CTRL+F for the example you want from the example files and paste it here in main!
Easiest way I promise

Current Modules:
    1. Bending and Deflection Functions:
         Contains functions for bending (normal stress, shear stress, deflections,
         axial and torsional deflections, discontinuity functions, and more)
    2. Stress and Strain Functions:
         Contains functions for stress/strain transformations, principal stress/strain
         calculations, constitutive relations, and plane stress/strain conditions

I'll probably end up making this into some kind of webapp? Doing only python is so clunky :/
I spent a while trying to think of the best way to organize this all lol, and here's what I got
"""

import math
import numpy as np
from sympy import symbols

from Bending_and_Deflection_Functions import *
from Stress_and_Strain_Functions import *


# =============================================================================
# Main function demonstrating usage of the modules.
# =============================================================================
def main():
    """
    Put all your functions in here, ones outside of main won't be run!
    Or I mean, they can, if you modify it

    Some examples, delete these and paste the ones you want
    You can use the same function over and over again btw, go crazy
    """

    # ---------------------------
    # Example Bending Calculation
    # ---------------------------
    # Compute normal bending stress using:  σ = -M * y / I
    M = 500  # Bending moment (Nm)
    y = 0.05  # Distance from the neutral axis (m)
    I = 8e-6  # Moment of inertia (m^4)
    sigma = normal_stress(M, y, I)
    print(f"Inputs: M = {M} Nm, y = {y} m, I = {I} m^4")
    print(f"Normal stress, σ = {sigma:.3f} Pa")

    # ---------------------------
    # Example Moment of a Force
    # ---------------------------
    # Compute the moment (torque) from a force using the cross product
    position = [1, 0, 0]  # Position vector (m)
    force = [0, 10, 0]  # Force vector (N)
    moment_vector = moment_of_force(position, force)
    print(f"Position: {position}, Force: {force}")
    print(f"Resulting moment vector: {moment_vector}")

    # ---------------------------
    # Example Strain Rosette
    # ---------------------------
    # Solves for εx, εy, and γxy from strain rosette measurements.
    # The strain rosette relation is given by:
    #   [1 + cos(2θ1)   1 - cos(2θ1)   2 sin(2θ1)]
    #   [1 + cos(2θ2)   1 - cos(2θ2)   2 sin(2θ2)]  *  [εx, εy, γxy/2]^T  =  [2εa, 2εb, 2εc]^T
    #   [1 + cos(2θ3)   1 - cos(2θ3)   2 sin(2θ3)]
    # Angles in radians
    angles = [0, math.pi / 4, math.pi / 2]
    strains = [0.001, 0.0012, 0.0008]  # measured strains εa, εb, εc
    eps_x_rosette, eps_y_rosette, gamma_xy_rosette = strain_rosette(angles, strains)
    print("Strain Rosette Conversion:")
    print(f"εx = {eps_x_rosette}")
    print(f"εy = {eps_y_rosette}")
    print(f"γxy = {gamma_xy_rosette}")

# =============================================================================
# You can hit the play button here or however you like
# =============================================================================
if __name__ == "__main__":
    main()
