import sys
import numpy as np
import pandas as pd
import subprocess
from composition_blend import blend_comps, make_composition_file
from list_isos import isos_from_net

if len(sys.argv) != 4:
    print("Error: Must provide 4 args-- python manual_composition.py <initial_mass> <csv_file> <network.net>")
    sys.exit(1)

initial_mass = sys.argv[1]
csv_file = sys.argv[2]
net_name = sys.argv[3]

# Read the CSV file & make sure it has an xq column
df = pd.read_csv(csv_file)
if 'xq' not in df.columns:
    print("Error: CSV must have an 'xq' column")
    sys.exit(1)


# Get list of isos 
network_isos = isos_from_net(net_name)

# Get isotope columns from CSV 
csv_isos = [col for col in df.columns if col != 'xq']

# Create the structured array data type with all network isotopes
dt = np.dtype([(iso, float) for iso in network_isos])

# Create comps list
comps = []
sample_xqs = df['xq'].values


for idx, row in df.iterrows():
    this_comp = np.zeros(1, dtype=dt)
    
    # Fill in isotopes from CSV
    for iso in csv_isos:
        if iso in network_isos:
            this_comp[iso] = row[iso]
        else:
            print(f"  Warning: Isotope {iso} from CSV not in network {net_name}")
    
    # Normalize so total mass fraction = 1.0
    total = sum(this_comp[iso] for iso in csv_isos if iso in network_isos)
    if total > 0:
        for iso in network_isos:
            this_comp[iso] = this_comp[iso] / total    
    comps.append(this_comp)


boundary_xqs = []
for i in range(len(sample_xqs)-1):
    boundary_xqs.append(sample_xqs[i])

configs = []

for i in range(len(boundary_xqs)):
    configs.append((boundary_xqs[i], 0, comps[i+1]))


blend = blend_comps(comps[0], configs)

# Generate output filename
output_filename = f'$PATHNAME/compositions/M{initial_mass[0]}P{initial_mass[2:]}_manual.data'
make_composition_file(blend, net_name, output_filename)

# Update MESA inlist files
subprocess.run(['shmesa', 'change', 'inlist_wd_builder', 'initial_mass', f'{initial_mass}'])
subprocess.run(['shmesa', 'change', 'inlist_wd_builder', 'relax_composition_filename', f"'{output_filename}'"])
subprocess.run(['shmesa', 'change', 'inlist_wd_builder', 'save_model_filename', 
               f"'$PATHNAME/outputs/M{initial_mass[0]}P{initial_mass[2:]}_manual.mod'"])


print(f"\nCreated composition file {output_filename}.")


