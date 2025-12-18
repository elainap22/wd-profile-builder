# wd-profile-builder
Tool that simplifies creation of '.data' composition files for use in 'wd_builder' module within MESA.

## Installation & Setup

It is assumed that you already have 'wd_builder' installed. See https://github.com/jschwab/wd_builder.git . 
To use this repository, you will need to install 'shMESA' as well. See https://github.com/earlbellinger/shmesa.git .

If you wish to create abundance plots of your WD's composition, you will also need to install 'MesaReader.' See 
https://github/com/wmwolf/py_mesa_reader.git .

To setup this repository, dowload and copy the contents into your 'wd_builder' directory. You will need to update the pathname in both 
'modular_composition.py' and 'manual_composition.py' (see lines 111 & 118 and lines 68 & 74, respectively). 

## Modular Method for Composition Building
You will need an existing MESA model of a white dwarf to use this method. Essentially, you will use the composition of this model as a reference for building the composition of this model.
Four parameters are required:
    initial_mass (float)
    model_name (string)
    sample_xqs (bracketed list of comma-separated floats)
    boundary_xqs (bracketed list of comma-separated floats)

boundary_xqs are the locations of the boundaries between the layers of the WD, and sample_xqs are some "midpoint" between each boundary_xq and the surface (xq=0) and the center (x=1). As a result, boundary_xqs must be exactly 1 element shorter than sample_xqs. The length of sample_xqs corresponds to the number of layers your model will have.

If you do not provide the lists, you will be prompted to enter them. You will also be prompted to create an abundance plot of your model with the boundary/sample locations drawn on to confirm that they match the structure you were intending.

To build the composition profile with this method, run:
      python modular_composition.py initial_mass model_name sample_xqs boundary_xqs

      
## Manual Method for Composition Building
This method requires a .csv file that outines the basic structure of the model you wish to create. See comps.csv for an example of this structure, but it requires an xq column (from 0 to 1), and an column for each isotope in the nuclear network you provide and their respective mass fractions at each xq.
Three parameters are required:
    initial_mass (float)
    csv_file (string)
    network_name.net (string)

Note that the .csv file does not need to be entirely complete; if you do not specify the mass fractions for an isotope that is in the chosen nuclear network, it will be filled in with zeros. If you specify the mass fractions for an isotope that is not in the chosen nucelar network, this mass will be "dumped" into the closest isotope (by my mass). 

To build the composition profile with this method, run:
      python manual_composition.py csv_file network_name.net

## Other Notes
Both of the methods utilize shMESA to update inlist_wd_builder for your convenience. After building the composition profile, you should be able to compile and run your model right away. 

Note that providing the initial mass does not have any physical meaning; it simply changes the initial_mass parameter in inlist_wd_builder, and it is included in the file name of the resulting composition profile & WD model. 
    

