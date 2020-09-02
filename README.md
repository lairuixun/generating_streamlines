# generating_streamlines
The code and example data were upload here for the interpolation of the bathymetry of the braided rivers. 

step 1:
Calculating the velocity field of the domain using shallow-water equations. 
01_velocity_field.plt is an example of the velocity field. 
wateredge.shp is an example of water edge file.

step 2:
Generating the streamlines using the Fortran code '02_generation_sl.f90'.
03_one_of_streamline.txt contains the coordinates of one of the generated streamlines.

step 3:
After generation of all the streamlines, the nodes should be transform into shapefile data type for interpolation.
The shapefile 'streamlines.shp' is an example data contains all the streamlines.

step 4:
The measured cross-sectional data is imposed on these streamlines to interpolate bathymetry. The interpolation method was described in Fig. 4.
