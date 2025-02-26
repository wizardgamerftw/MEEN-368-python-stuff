"""
This is the Table A-20
Deterministic ASTM Minimum Tensile and Yield Strengths for Some Hot-Rolled (HR) and Cold-Drawn (CD) Steels

It's on our google drive too

Instructions:
    1. Run the code
    2. Answer the questions
    3. Bam material info

Inputs:
    - The UNS or SAE AISI number, the 4-6 number thing
    - Hot-Rolled (HR) or Cold-Drawn = (CD)
    - MPa or KPsi
    - Quit (optional if you want out without hitting the button)

Outputs:
    - UNS (to check)
    - SAE/AISI (to check)
    - Processing (to check)
    - Tensile Strength
    - Yield Strength
    - Elongation in 2 in %
    - Reduction in Area
    - Brinell Hardness


I'm trying out the method of asking for inputs instead of changing variables, what do you prefer?
"""
import sys

# ------------------------
# Inputs
# ------------------------
# """
# Manual inputs for the UNS/SAE number and processing
# """
# # MATERIAL_NUMBER can be either the UNS (e.g., "G10180") or the SAE/AISI number (e.g., "1018").
# MATERIAL_NUMBER = "1018"
# # PROCESSING_TYPE should be either "HR" (Hot-Rolled) or "CD" (Cold-Drawn).
# PROCESSING_TYPE = "CD"
print("This file is an input based one, testing it out. Type 'quit' if you ever want out immediately\n------ \n")

def input_or_quit(prompt):
    """
    Gets user input and exits if the user types 'quit'
    """
    user_input = input(prompt).strip()
    if user_input.lower() == "quit":
        sys.exit("User requested exit. Quitting")
    return user_input

def get_user_input():
    """
    Prompts the user for the material number and processing type.
    Checks if the material number exists in the data.
    Raises a ValueError if the material number is not found or the processing type is not 'HR' or 'CD'.
    """
    print("Enter either the UNS number (e.g., G10180) or the SAE/AISI number (e.g., 1018).")
    user_material = input_or_quit("Material number: ")

    # Create a set of all valid material identifiers (both UNS and SAE/AISI)
    valid_materials = {row["UNS"].lower() for row in STEEL_DATA} | {row["SAE_AISI"].lower() for row in STEEL_DATA}
    if user_material.lower() not in valid_materials:
        raise ValueError("Material number not found in the table, is it a steel? Check your input")

    user_processing = input_or_quit("Processing [HR/CD]: ").upper()
    if user_processing not in ("HR", "CD"):
        raise ValueError("Invalid processing type. Enter either 'HR' or 'CD'")
    return user_material, user_processing

def get_unit_input():
    """
    Prompts the user for the desired units for tensile and yield strengths
    Loops until valid input ("MPA" or "KPSI") is provided
    """
    while True:
        unit_input = input_or_quit("Which units do you want for Tensile and Yield Strength? [MPA/KPSI]: ").strip().upper()
        if unit_input in ("MPA", "KPSI"):
            return unit_input
        else:
            print("Invalid units. Enter 'MPA' or 'KPSI'")

# STEEL_DATA dictionary corresponding to each row in the table
# UNS and SAE/AISI numbers, processing, and associated material properties

STEEL_DATA = [
    {
        "UNS": "G10060",
        "SAE_AISI": "1006",
        "Processing": "HR",
        "TensileStrength_MPa": 300,
        "TensileStrength_kpsi": 43,
        "YieldStrength_MPa": 170,
        "YieldStrength_kpsi": 24,
        "Elongation_2in_percent": 30,
        "ReductionArea_percent": 55,
        "BrinellHardness": 86,
    },
    {
        "UNS": "G10060",
        "SAE_AISI": "1006",
        "Processing": "CD",
        "TensileStrength_MPa": 330,
        "TensileStrength_kpsi": 48,
        "YieldStrength_MPa": 280,
        "YieldStrength_kpsi": 41,
        "Elongation_2in_percent": 20,
        "ReductionArea_percent": 45,
        "BrinellHardness": 95,
    },
    {
        "UNS": "G10100",
        "SAE_AISI": "1010",
        "Processing": "HR",
        "TensileStrength_MPa": 320,
        "TensileStrength_kpsi": 47,
        "YieldStrength_MPa": 180,
        "YieldStrength_kpsi": 26,
        "Elongation_2in_percent": 28,
        "ReductionArea_percent": 50,
        "BrinellHardness": 95,
    },
    {
        "UNS": "G10100",
        "SAE_AISI": "1010",
        "Processing": "CD",
        "TensileStrength_MPa": 370,
        "TensileStrength_kpsi": 53,
        "YieldStrength_MPa": 300,
        "YieldStrength_kpsi": 44,
        "Elongation_2in_percent": 20,
        "ReductionArea_percent": 40,
        "BrinellHardness": 105,
    },
    {
        "UNS": "G10150",
        "SAE_AISI": "1015",
        "Processing": "HR",
        "TensileStrength_MPa": 340,
        "TensileStrength_kpsi": 50,
        "YieldStrength_MPa": 190,
        "YieldStrength_kpsi": 27.5,
        "Elongation_2in_percent": 28,
        "ReductionArea_percent": 50,
        "BrinellHardness": 101,
    },
    {
        "UNS": "G10150",
        "SAE_AISI": "1015",
        "Processing": "CD",
        "TensileStrength_MPa": 390,
        "TensileStrength_kpsi": 56,
        "YieldStrength_MPa": 320,
        "YieldStrength_kpsi": 47,
        "Elongation_2in_percent": 18,
        "ReductionArea_percent": 40,
        "BrinellHardness": 111,
    },
    {
        "UNS": "G10180",
        "SAE_AISI": "1018",
        "Processing": "HR",
        "TensileStrength_MPa": 400,
        "TensileStrength_kpsi": 58,
        "YieldStrength_MPa": 220,
        "YieldStrength_kpsi": 32,
        "Elongation_2in_percent": 25,
        "ReductionArea_percent": 50,
        "BrinellHardness": 116,
    },
    {
        "UNS": "G10180",
        "SAE_AISI": "1018",
        "Processing": "CD",
        "TensileStrength_MPa": 440,
        "TensileStrength_kpsi": 64,
        "YieldStrength_MPa": 370,
        "YieldStrength_kpsi": 54,
        "Elongation_2in_percent": 15,
        "ReductionArea_percent": 40,
        "BrinellHardness": 126,
    },
    {
        "UNS": "G10200",
        "SAE_AISI": "1020",
        "Processing": "HR",
        "TensileStrength_MPa": 380,
        "TensileStrength_kpsi": 55,
        "YieldStrength_MPa": 210,
        "YieldStrength_kpsi": 30,
        "Elongation_2in_percent": 25,
        "ReductionArea_percent": 50,
        "BrinellHardness": 111,
    },
    {
        "UNS": "G10200",
        "SAE_AISI": "1020",
        "Processing": "CD",
        "TensileStrength_MPa": 470,
        "TensileStrength_kpsi": 68,
        "YieldStrength_MPa": 390,
        "YieldStrength_kpsi": 57,
        "Elongation_2in_percent": 15,
        "ReductionArea_percent": 40,
        "BrinellHardness": 131,
    },
    {
        "UNS": "G10300",
        "SAE_AISI": "1030",
        "Processing": "HR",
        "TensileStrength_MPa": 470,
        "TensileStrength_kpsi": 68,
        "YieldStrength_MPa": 260,
        "YieldStrength_kpsi": 37.5,
        "Elongation_2in_percent": 20,
        "ReductionArea_percent": 42,
        "BrinellHardness": 137,
    },
    {
        "UNS": "G10300",
        "SAE_AISI": "1030",
        "Processing": "CD",
        "TensileStrength_MPa": 520,
        "TensileStrength_kpsi": 76,
        "YieldStrength_MPa": 440,
        "YieldStrength_kpsi": 64,
        "Elongation_2in_percent": 12,
        "ReductionArea_percent": 35,
        "BrinellHardness": 149,
    },
    {
        "UNS": "G10350",
        "SAE_AISI": "1035",
        "Processing": "HR",
        "TensileStrength_MPa": 500,
        "TensileStrength_kpsi": 72,
        "YieldStrength_MPa": 270,
        "YieldStrength_kpsi": 39.5,
        "Elongation_2in_percent": 18,
        "ReductionArea_percent": 40,
        "BrinellHardness": 143,
    },
    {
        "UNS": "G10350",
        "SAE_AISI": "1035",
        "Processing": "CD",
        "TensileStrength_MPa": 550,
        "TensileStrength_kpsi": 80,
        "YieldStrength_MPa": 460,
        "YieldStrength_kpsi": 67,
        "Elongation_2in_percent": 12,
        "ReductionArea_percent": 35,
        "BrinellHardness": 163,
    },
    {
        "UNS": "G10400",
        "SAE_AISI": "1040",
        "Processing": "HR",
        "TensileStrength_MPa": 520,
        "TensileStrength_kpsi": 76,
        "YieldStrength_MPa": 290,
        "YieldStrength_kpsi": 42,
        "Elongation_2in_percent": 18,
        "ReductionArea_percent": 40,
        "BrinellHardness": 149,
    },
    {
        "UNS": "G10400",
        "SAE_AISI": "1040",
        "Processing": "CD",
        "TensileStrength_MPa": 590,
        "TensileStrength_kpsi": 85,
        "YieldStrength_MPa": 490,
        "YieldStrength_kpsi": 71,
        "Elongation_2in_percent": 12,
        "ReductionArea_percent": 35,
        "BrinellHardness": 170,
    },
    {
        "UNS": "G10450",
        "SAE_AISI": "1045",
        "Processing": "HR",
        "TensileStrength_MPa": 570,
        "TensileStrength_kpsi": 82,
        "YieldStrength_MPa": 310,
        "YieldStrength_kpsi": 45,
        "Elongation_2in_percent": 16,
        "ReductionArea_percent": 40,
        "BrinellHardness": 163,
    },
    {
        "UNS": "G10450",
        "SAE_AISI": "1045",
        "Processing": "CD",
        "TensileStrength_MPa": 630,
        "TensileStrength_kpsi": 91,
        "YieldStrength_MPa": 530,
        "YieldStrength_kpsi": 77,
        "Elongation_2in_percent": 12,
        "ReductionArea_percent": 35,
        "BrinellHardness": 179,
    },
    {
        "UNS": "G10500",
        "SAE_AISI": "1050",
        "Processing": "HR",
        "TensileStrength_MPa": 620,
        "TensileStrength_kpsi": 90,
        "YieldStrength_MPa": 340,
        "YieldStrength_kpsi": 49.5,
        "Elongation_2in_percent": 15,
        "ReductionArea_percent": 35,
        "BrinellHardness": 179,
    },
    {
        "UNS": "G10500",
        "SAE_AISI": "1050",
        "Processing": "CD",
        "TensileStrength_MPa": 690,
        "TensileStrength_kpsi": 100,
        "YieldStrength_MPa": 580,
        "YieldStrength_kpsi": 84,
        "Elongation_2in_percent": 10,
        "ReductionArea_percent": 30,
        "BrinellHardness": 197,
    },
    {
        "UNS": "G10600",
        "SAE_AISI": "1060",
        "Processing": "HR",
        "TensileStrength_MPa": 680,
        "TensileStrength_kpsi": 98,
        "YieldStrength_MPa": 370,
        "YieldStrength_kpsi": 54,
        "Elongation_2in_percent": 12,
        "ReductionArea_percent": 30,
        "BrinellHardness": 201,
    },
    {
        "UNS": "G10800",
        "SAE_AISI": "1080",
        "Processing": "HR",
        "TensileStrength_MPa": 770,
        "TensileStrength_kpsi": 112,
        "YieldStrength_MPa": 420,
        "YieldStrength_kpsi": 61.5,
        "Elongation_2in_percent": 10,
        "ReductionArea_percent": 25,
        "BrinellHardness": 229,
    },
    {
        "UNS": "G10950",
        "SAE_AISI": "1095",
        "Processing": "HR",
        "TensileStrength_MPa": 830,
        "TensileStrength_kpsi": 120,
        "YieldStrength_MPa": 460,
        "YieldStrength_kpsi": 66,
        "Elongation_2in_percent": 10,
        "ReductionArea_percent": 25,
        "BrinellHardness": 248,
    },
]


def get_steel_properties():
    """
    Continuously prompts the user for input until a valid material number and processing type are provided.
    If no matching data is found, raises an error and asks for input again.
    Then prompts for desired units and prints only the strength values in those units.
    """
    while True:
        try:
            user_material, user_processing = get_user_input()

            # Look for matching records in the STEEL_DATA list.
            matches = []
            for row in STEEL_DATA:
                if ((row["UNS"].lower() == user_material.lower() or
                     row["SAE_AISI"].lower() == user_material.lower()) and
                    row["Processing"].upper() == user_processing):
                    matches.append(row)

            if not matches:
                raise ValueError("Material number not found in our data. Please check your input.")

            # If we have valid matches, break out of the loop.
            break

        except ValueError as err:
            print(err)
            print("Something was wrong with the input, try again\n")

    # Prompt for desired units.
    units = get_unit_input()

    # Display the results.
    for m in matches:
        print("\nFound!")
        print(f"UNS: {m['UNS']}, SAE/AISI: {m['SAE_AISI']}, Processing: {m['Processing']}")
        if units == "MPA":
            print(f"Tensile Strength: {m['TensileStrength_MPa']} MPa")
            print(f"Yield Strength:   {m['YieldStrength_MPa']} MPa")
        else:  # units == "KPSI"
            print(f"Tensile Strength: {m['TensileStrength_kpsi']} kpsi")
            print(f"Yield Strength:   {m['YieldStrength_kpsi']} kpsi")
        print(f"Elongation in 2 in: {m['Elongation_2in_percent']}%")
        print(f"Reduction in Area:  {m['ReductionArea_percent']}%")
        print(f"Brinell Hardness:   {m['BrinellHardness']}")

if __name__ == "__main__":
    get_steel_properties()