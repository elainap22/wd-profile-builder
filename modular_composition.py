import sys
import numpy as np
import mesa_reader as mr
import subprocess
import matplotlib.pyplot as plt
from composition_blend import blend_comps, make_composition_file

def texify(iso):
    "Create a TeX-friendly version of an isotope string."
    element = ''.join([i for i in iso if not i.isdigit()])
    mass_number = ''.join([i for i in iso if i.isdigit()])
    return r"$^{" + str(mass_number) + r"}\mathrm{" + f"{element.title()}" + r"}$"

# Check if the required first argument (initial mass) is provided
if len(sys.argv) < 3:
    print("Error: Initial mass (float) and model file name (string) are required")
    sys.exit(1)

initial_mass = sys.argv[1]
model = mr.MesaData(sys.argv[2]) 

# Function to parse a string of floats in bracket format [1.0,2.0,3.0]
def parse_float_list(arg_string):
    try:
        # Remove brackets and split by comma
        cleaned = arg_string.strip('[]')
        return [float(x.strip()) for x in cleaned.split(',')]
    except ValueError:
        print(f"Error: Invalid float in argument: {arg_string}")
        sys.exit(1)

# Check if optional arguments are provided
if len(sys.argv) >= 5:  # All four arguments provided
    sample_xqs = parse_float_list(sys.argv[3])
    boundary_xqs = parse_float_list(sys.argv[4])
else:  # Only initial mass provided- prompt for both lists
    user_samples = input("Enter sample locations as a list of comma-separated floats within brackets: ")
    sample_xqs = parse_float_list(user_samples)
    
    user_bounds = input("Enter boundary locations as a list of comma-separated floats within brackets: ")
    boundary_xqs = parse_float_list(user_bounds)

if len(boundary_xqs) != len(sample_xqs) - 1:
    print(f"Error: boundary_xqs must be exactly 1 element shorter than sample_xqs")
    sys.exit(1)


xqs = np.cumsum(model.dq)
num_layers = len(sample_xqs)

# Ask user if they want to see a plot of compositions and the lines they have drawn on their given model
plot_bool = input("Would you like to create a plot of your model's composition and sample locations? (y/n): ").lower().strip()
if plot_bool == 'y':
    for iso in ['h1', 'he4', 'c12', 'n14', 'o16', 'ne20']:
        plt.loglog(xqs, model.data(iso), label=texify(iso))
    
    plt.legend(loc='best')
    plt.xlim(1, 1e-8)
    plt.ylim(1.5e-4, 1.5)
    plt.xlabel(r"Exterior Fractional Mass")
    plt.ylabel(r"Mass Fraction")

    # Show vertical dashed lines at sample locations
    for xq in sample_xqs:
        plt.axvline(xq, ls='--', color='k')
    # Show vertical dotted lines at boundary locations
    for xq in boundary_xqs:
        plt.axvline(xq, ls=':', color='lightgray')
    
    plt.legend(loc='best')
    plt.savefig(f'{sys.argv[2]}_plot.pdf', bbox_inches='tight')
    print(f"Plot saved as {sys.argv[2]}_plot.pdf")

# Extract composition at each sample location
isos = model.bulk_names[model.bulk_names.index('h1'):]

dt = np.dtype([(iso, float) for iso in isos])
comps = []


# Sort xqs in descending order
sorted_sample_xqs = sorted(sample_xqs, reverse=True) 
sorted_boundary_xqs = sorted(boundary_xqs, reverse=True)

for i, xq in enumerate(sorted_sample_xqs):
    this_comp = np.zeros(1, dtype=dt)
    idx = np.argmin(np.abs(xqs - xq))
    for iso in isos:
        this_comp[iso] = model.data(iso)[idx]
    comps.append(this_comp)


configs = []
# Sort boundary xqs to be back in ascending order
for i in range(len(sorted_boundary_xqs)):
    boundary_idx = len(sorted_boundary_xqs) - 1 - i  
    configs.append((sorted_boundary_xqs[boundary_idx], 0, comps[i+1]))


# Blend the compositions
blend = blend_comps(comps[0], configs)

# Reverse the order of all columns except xq 
blend_fixed = np.zeros_like(blend)
blend_fixed['xq'] = blend['xq']
for iso in blend.dtype.names[1:]: 
    blend_fixed[iso] = blend[iso][::-1]
blend = blend_fixed

# Generate output filename
output_filename = f'$PATHNAME/compositions/M{initial_mass[0]}P{initial_mass[2:]}_CO_WD.data'

make_composition_file(blend, model.header('net_name'), output_filename)

# Update MESA inlist files
subprocess.run(['shmesa', 'change', 'inlist_wd_builder', 'initial_mass', f'{initial_mass}'])
subprocess.run(['shmesa', 'change', 'inlist_wd_builder', 'relax_composition_filename', f"'{output_filename}'"])
subprocess.run(['shmesa', 'change', 'inlist_wd_builder','save_model_filename', f"'$PATHNAME/outputs/M{initial_mass[0]}P{initial_mass[2:]}_CO_WD.mod'"])

print(f"\nCreated composition file {output_filename}.")