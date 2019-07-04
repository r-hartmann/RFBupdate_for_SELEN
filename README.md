# RFBupdate_for_SELEN
Repository including a guide and all required files to update SELEN2.9.12/2.9.13
for consideration of rotational feedback (RFB) in the sea level equation (SLE).

################################################################################

SELEN is free software: you can redistribute it and/or modify it under the terms 
of the GNU General Public License as published by the Free Software Foundation, 
either version 3 of the License, or at your option) any later version. 

SELEN is distributed in the hope that it will be useful, but WITHOUT ANY 
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
PARTICULAR PURPOSE. See the GNU General Public License for more details. 

You should have received a copy of the GNU General Public License along with 
SELEN. If not, see <http://www.gnu.org/licenses/>

################################################################################

BUILD INSTRUCTIONS
################################################################################

1. Download and install SELEN 2.9.12 or 2.9.13 from the CIG website, a general 
   installation guide is given in the official SELEN user manual:
   https://geodynamics.org/cig/software/selen/
   
2. Download the repository: 
   https://github.com/r-hartmann/RFBupdate_for_SELEN/src_RFBupdate/
   
3. Copy all .F90/.f90-files from the repository folder ./src_RFBupdate/ into 
   your SELEN source folder ./selen-2.9.12/src/ (replace original files).

4. Replace the Makefile.in in your SELEN directory ./selen-2.9.12/ with the new
   Makefile.in from the repository.
   
5. Replace the config.dat in your SELEN directory ./selen-2.9.12/ with the new 
   config.dat from the repository.

6. DONE!

################################################################################

EXECUTION 
################################################################################

SELEN is still executable in the SELEN directory ./selen-2.9.12 via 
> sh configure
> make run

The consideration of rotational feedback can be switched on/off in option 131 of
config.dat via the parameter 'y'/'n':
>131    Including rotational feedback   'y'

################################################################################
Version 1.0.0
Robert Hartmann (stu200105 AT mail DOT uni-kiel DOT de)
04.07.2019
