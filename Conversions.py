"""
I'm traumatized from my 363 exam where the units were really weird

This is for unit conversions, and uses pint

Should work for all of the units we use?

Instructions:
    1. Ensure you have pint installed, if not, in your command line run:
           pip install pint
    2. Modify the 'initial_value' variable to the value and unit you want to convert
    3. Change the conversion target in the 'to()' function as needed
    4. Adjust the 'decimals' variable to set the number of decimal places for the output
    5. Should print when run?

Example:
    Converting 10 meters to feet should output something like:
           10.000 meter = 32.808 foot
"""
import pint
import sys

def main():
    # Create a unit registry
    ureg = pint.UnitRegistry()

    # Define input and output units as strings
    input_unit_str = "meter"
    output_unit_str = "feet"

    # Check if both the input and output units exist in the registry, I'm paranoid lol
    try:
        input_unit = ureg.Unit(input_unit_str)
        output_unit = ureg.Unit(output_unit_str)
        print("Both input and output units exist!")
    except Exception as e:
        print("Error: One of the units does not exist:", e)
        sys.exit(1)

    # Set the number of decimal places you want
    decimals = 3

    # Define the initial quantity to convert
    initial_value = 10 * input_unit

    # Define the conversion unit
    converted_value = initial_value.to(output_unit)

    # Format and print the result
    print(f"{initial_value.magnitude:.{decimals}f} {initial_value.units} = {converted_value.magnitude:.{decimals}f} {converted_value.units}")


# You can press the button to run it
if __name__ == "__main__":
    main()
