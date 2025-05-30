# wd-profile-builder


Initial Setup:

You will need to have shMESA installed to use this repository. See https://github.com/earlbellinger/shmesa.git to install. 

Open the file 'make_compositions.py' and change the update the #PATHNAME in lines 82, 86 & 87. Once
this is updates to the directory that 'wd-profile-builder' lives in on your machine, all outputs will
automatically be sent to the 'outputs' directory, and the composition files that you build wil be sent to 
the 'compositions' directory.


Building the Compostion:

Before you can create your model, you will need to specify the thicknesses of each layer. To do this, we extract 
the surface, the middle of the He layer, the middle of the outer core, the center of the inner core, and the three transition
points between those four layers from older white dwarf models (found in the 'wolf_2013_models' directory). To see how this is done, look 
look at the 'examples.ipynb' notebook; this will help visualize the structure. 

Once you have decided on what the structure of your star might look like, populate the lists sample_xqs and boundary_xqs in 
lines 55 & 56, respectively. The composition array will be built automatically around these points. 



Running a Model:

To create the most basic of basic white dwarf models, execute the following:

    python make_compositions.py *initial_mass*

All changes to the structure of the composition file can be made in 'make_compositions.py'; all other changes (temperature, luminosity, etc.)
can be make in 'inlist_wd_builder', using either a text editor or shmesa from the command line. 


Once you're ready, execute:

    ./mk && ./rn
