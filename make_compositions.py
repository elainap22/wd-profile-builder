## Be sure to change PATHNAME to your directory in lines 82, 86 & 87!

from os.path import join, abspath, dirname

import numpy as np

import subprocess

import mesa_reader as mr
from composition_blend import blend_comps, make_composition_file

WOLF_MODELS_DIR = abspath(join('.', 'wolf_2013_models'))

model_0p51 = mr.MesaData(join(WOLF_MODELS_DIR, '0.51Msun_Tc_3e7_init.mod'))
model_0p60 = mr.MesaData(join(WOLF_MODELS_DIR, '0.6Msun_Tc_3e7_init.mod'))
model_0p65 = mr.MesaData(join(WOLF_MODELS_DIR, '0.65Msun_Tc_3e7_init.mod'))
model_0p70 = mr.MesaData(join(WOLF_MODELS_DIR, '0.7Msun_Tc_3e7_init.mod'))
model_0p80 = mr.MesaData(join(WOLF_MODELS_DIR, '0.8Msun_Tc_3e7_init.mod'))
model_0p90 = mr.MesaData(join(WOLF_MODELS_DIR, '0.9Msun_Tc_3e7_init.mod'))
model_1p00 = mr.MesaData(join(WOLF_MODELS_DIR, '1.0Msun_Tc_3e7_init.mod'))
model_1p05 = mr.MesaData(join(WOLF_MODELS_DIR, '1.05Msun_Tc_3e7_init.mod'))
model_1p10 = mr.MesaData(join(WOLF_MODELS_DIR, '1.1Msun_Tc_3e7_init.mod'))
model_1p15 = mr.MesaData(join(WOLF_MODELS_DIR, '1.15Msun_Tc_3e7_init.mod'))
model_1p20 = mr.MesaData(join(WOLF_MODELS_DIR, '1.2Msun_Tc_3e7_init.mod'))
model_1p30 = mr.MesaData(join(WOLF_MODELS_DIR, '1.3Msun_Tc_3e7_init.mod'))
model_1p32 = mr.MesaData(join(WOLF_MODELS_DIR, '1.32Msun_Tc_3e7_init.mod'))
model_1p34 = mr.MesaData(join(WOLF_MODELS_DIR, '1.34Msun_Tc_3e7_init.mod'))
model_1p36 = mr.MesaData(join(WOLF_MODELS_DIR, '1.36Msun_Tc_6e7_init.mod'))

def texify(iso):
    "Create a TeX-friendly version of an isotope string."
    element = ''.join([i for i in iso if not i.isdigit()])
    mass_number = ''.join([i for i in iso if i.isdigit()])
    return r"$^{" + str(mass_number) + r"}\mathrm{" + f"{element.title()}" + r"}$"



import sys

## find corresponding "old" model to use as reference 
initial_mass = float(sys.argv[1])
org_initial_mass = str(initial_mass)
model_masses = [0.51, 0.60, 0.65, 0.70, 0.80, 0.90, 1.00, 1.05, 1.10, 1.15, 1.20, 1.30, 1.32, 1.34, 1.36]
delta_masses = []
for i in range(len(model_masses)):
    delta_masses.append(model_masses[i]-initial_mass)
initial_mass = str(model_masses[delta_masses.index(min(delta_masses, key = abs))])
model = mr.MesaData(join(WOLF_MODELS_DIR, f'{initial_mass[0]}.{initial_mass[2:]}Msun_Tc_3e7_init.mod'))


xqs = np.cumsum(model.dq)


## Specify your layer thicknesses & boundary locations
sample_xqs = []
boundary_xqs = []
samples = {'inner_core': sample_xqs[0], 'outer_core': sample_xqs[1],
           'helium_shell':sample_xqs[2], 'hydrogen_shell': sample_xqs[3]}
boundaries = {'inner_core': boundary_xqs[0], 'outer_core': boundary_xqs[1],
              'helium_core': boundary_xqs[2]}



# extract composition at each sample location
isos = model.bulk_names[model.bulk_names.index('h1'):]

dt = np.dtype([(iso, float) for iso in isos])
comps = []
for xq in sorted(sample_xqs):
    this_comp = np.zeros(1, dtype=dt)
    idx = np.argmin(np.abs(xqs - xq))
    for iso in isos:
        this_comp[iso] = model.data(iso)[idx]
    comps.append(this_comp)

# blend the compositions
blend = blend_comps(comps[0], [(boundaries['helium_core'], 0, comps[1]),
                                (boundaries['outer_core'], 0, comps[2]),
                                (boundaries['inner_core'], 0, comps[3])])


make_composition_file(blend, model.header('net_name'), f'/#PATHNAME/compositions/M{initial_mass[0]}P{initial_mass[2:]}_CO_WD.data')


subprocess.run(['shmesa', 'change', 'inlist_wd_builder', 'initial_mass', f'{initial_mass}'])
subprocess.run(['shmesa', 'change', 'inlist_wd_builder', 'relax_composition_filename', f"'/#PATHNAME/compositions/M{initial_mass[0]}P{initial_mass[2:]}_CO_WD.data'"])
subprocess.run(['shmesa', 'change', 'inlist_wd_builder','save_model_filename', f"'/#PATHNAME/outputs/M{initial_mass[0]}P{initial_mass[2:]}_CO_WD.mod'"])