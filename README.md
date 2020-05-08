# ExoSPEEDY
## The ICTP atmospheric general circulation model SPEEDY, adapted for exoplanets

The software in this repository is adapted from the ICTP AGCM software, nicknamed SPEEDY, see:
https://www.ictp.it/research/esp/models/speedy.aspx

The SPEEDY software was developed at ICTP by Franco Molteni and Fred Kucharski and is described in the publications referenced here:
http://users.ictp.it/~kucharsk/speedy-net.html

Please cite these publications when publishing your own work based on this software.

ExoSPEEDY adapts the SPEEDY software for use in climate modeling of exoplanets.
This work was done at Delft University of Technology under the guidance of Daphne Stam.

## Installing and running ExoSPEEDY

Please read doc_ver41.5.pdf.

## Ocean Model
The standard version of Speedy has a vertical slab ocean model, whereby the water columns cannot interact laterally with each other. The changes to the model in the ocean_modle folder make it possible for heat to diffuse horizontally through the ocean from the equator to the pole. This means that surplus incoming radiation can be transported to the poles by both the atmosphere and the ocean.


These changes were made by Frank Selten of the KNMI (www.knmi.nl) for a lecture series ‘Water in the Atmosphere’, in a collaboration with Pier Sibesma and Stephan de Roode from the Delft University of Technology (www.tudelft.nl).


To use this model, copy the content of the ‘ocean_modle’ folder to the ‘update’ folder before running a calculation. The run_exp.s script will then use the contents of the ‘update’ folder for the run.


These ocean model changes may be used freely and are provided ‘as is’, without warranty of any kind.

## Comments in the code
The first step in adapting the code for exoplanets was to identify the flow of earth parameters.
To this purpose, each subroutine is annotated with comment lines starting with the string 'C--IO x' where x is one of:
- r - reads from a file
- w - writes to a file
- h - gets input from a .h file
- s - sets a value in a .f or .h file

## Notes for use with exoplanets
Planet specific parameters that were embedded in the original software, can be found in cpl\_inplanet.h. Other relevant parameters have not been moved from their original cls\_in\*.h file.

Daily mean upper air levels are fixed in ppo\_tminc.f.

The ozone absorption model, polar night cooling in the stratosphere model, stratiform cloud model and energy fractions in LW bands model (assumes 100 to 400K) are defined in phy\_radiat.f.

The water specific model for sea temperature is defined in cpl\_sea.f and cpl\_sea\_model.f

For the ocean model, solar radiation is fixed for the equinox in phy\_radiat.f.

## How the main script works

The main script to run an experiment is speedy\_ver41.5/run/run\_exp.s.
It works by copying input, .f, .h, .s and makefiles from various locations to a speedy\_ver41.5/tmp directory. The copy actions can overwrite previously copied files, so it is important to make any modifications in the correct location. Below is an overview of how the script works. Paths are relative to the speedy\_ver41.5 directory.

**action** | **file** | **from** | **to**
--- | --- | --- | ---
user can edit in situ: | ver41.5.input/cls\_instep.h |  | 
user can edit in situ: | ver41.5.input/cls\_indyns.h |  | 
user can edit in situ: | ver41.5.input/cls\_inphys.h |  | 
user can edit in situ: | ver41.5.input/cls\_inland.h |  | 
user can edit in situ: | ver41.5.input/cls\_insea.h |  | 
user can edit in situ: | ver41.5.input/inpfiles.s |  | 
model parameters written to: | input/exp\_$2/run\_setup |  | 
force remove | * | tmp | 
copy files | makefile, \*.f, \*.h, \*.s | source | tmp
copy file |  | tmp/par\_horres\_$1.h | tmp/atparam.h
copy file |  | tmp/par\_verres.h | tmp/atparam1.h
copy files | \*.h, inpfiles.s | ver41.5.input | tmp
copy files | \*.h, inpfiles.s | ver41.5.input | input/exp\_$2
copy files | \*.f, \*.h, make\* | update | tmp
copy files | \*.f, \*.h, make\* | update | input/exp\_$2
experiment no. and restart file no. written to: | tmp/fort.2 |  | 
output file fort.3 linked to: | output/exp\_$3/atgcm$3.rst |  | 
run inpfiles.s to link fortran units |  |  | 
run make to create: | tmp/imp.exe |  | 
run tmp/imp.exe | |  | 
move | out.lis | tmp | output/exp\_$2/atgcm$2.lis
move | fort.10 | tmp | output/exp\_$2/atgcm$2.rst
move | at\*$2.ctl, at\*$2\_\*.grd | tmp | output/exp\_$2
move | day\*$2.ctl, day\*$2\_\*.grd | tmp | output/exp\_$2


$1 = resolution (eg t21, t30)

$2 = experiment no. (eg 111)

$3 = experiment no. for restart file ( 0 = no restart )
