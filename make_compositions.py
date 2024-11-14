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
 
initial_mass = float(sys.argv[1])


initial_mass = str(initial_mass)
model = mr.MesaData(join(WOLF_MODELS_DIR, f'{initial_mass[0]}.{initial_mass[2:]}Msun_Tc_3e7_init.mod'))

xqs = np.cumsum(model_1p00.dq)


# eyeballed from plot
sample_xqs = [1e-8, 1e-4, 2e-2, 1]
core_shell_xq = 0.45
core_shell_dq = 0.1
core_shell_xm = core_shell_xq * model_1p00.header('M/Msun')
core_shell_dm = core_shell_dq * model_1p00.header('M/Msun')


isos = model_1p00.bulk_names[model_1p00.bulk_names.index('h1'):]

# create composition arrays for each sample location
dt = np.dtype([(iso, float) for iso in isos])
comps = []
for xq in sample_xqs:
    this_comp = np.zeros(1, dtype=dt)
    # find location in model closest to xq
    idx = np.argmin(np.abs(xqs - xq))
    for iso in isos:
        this_comp[iso] = model_1p00.data(iso)[idx]
    comps.append(this_comp)


xqs = np.cumsum(model.dq)

# boundary between core from core burning and core from shell burning
core_shell_xq = core_shell_xm / model.header('M/Msun')
core_shell_dq = core_shell_dm / model.header('M/Msun')
# boundary between he layer and metal "core" (really with the core from shell
# burning)
metal_core_xq = -1
metal_core_dq = 0

# boundary between he layer and hydrogren envelope
he_core_xq = -1
he_core_dq = 0

# determine the h-he boundary and he-core boundary for the 1.3 Msun model
xqs = np.cumsum(model_1p30.dq)
for k, (xq, h1, he4) in enumerate(zip(xqs, model.h1, model.he4)):
    # account for transition zone between h-rich and he-rich layers
    if 1e-3 * model.h1[0] < h1 < 0.98 * model.h1[0] and he_core_xq < 0:
        he_core_dq += model.dq[k]
    if h1 < 0.01 and he_core_xq < 0:
        he_core_xq = xq
    # account for transition zone between he-rich and metal-rich layers
    if he_core_xq > 0 and 0.05 < he4 / max(model.he4) < 0.95:
        metal_core_dq += model.dq[k]
    if he4 < 0.01 and metal_core_xq < 0:
        metal_core_xq = xq
    if metal_core_xq > 0 and he4 < 1e-4:
        break
print(f"Core/shell boundary at xq = {core_shell_xq:.2e} with dq = {core_shell_dq:.2e}")
print(f"Metal core at xq = {metal_core_xq:.2e} with dq = {metal_core_dq:.2e}")
print(f"He core at xq = {he_core_xq:.2e} with dq = {he_core_dq:.2e}")

surf_comp = comps[0]
configs = [
    (he_core_xq, he_core_dq, comps[1]),
    (metal_core_xq, metal_core_dq, comps[2]),
    (core_shell_xq, core_shell_dq, comps[3])
]
blend = blend_comps(surf_comp, configs)


make_composition_file(blend, model.header('net_name'), f'/home/elainaplonis/mesa-r24.03.1/wd_suite/shmesa_suite/compositions/M{initial_mass[0]}P{initial_mass[2:]}_CO_WD.data')


subprocess.run(['shmesa', 'change', 'inlist_wd_builder', 'initial_mass', f'{initial_mass}'])
subprocess.run(['shmesa', 'change', 'inlist_wd_builder', 'relax_composition_filename', f"'compositions/M{initial_mass[0]}P{initial_mass[2:]}_CO_WD.data'"])
subprocess.run(['shmesa', 'change', 'inlist_wd_builder','save_model_filename', f"'outputs/M{initial_mass[0]}P{initial_mass[2:]}_CO_WD.mod'"])