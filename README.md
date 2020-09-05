The code and example data were upload here for the interpolation of the bathymetry of the braided rivers. 

step 1:
Generating an initial bathymetry based on represented cross-sections.
TIN_Initial_Bathymetry is the result that can be opened by GIS software such as ArcGIS or QGIS.

step 2:
Calculating the velocity field of the domain using shallow-water equations. 
01_velocity_field.plt is an example of the velocity field. 
shp_wateredge is an example of water edge file.
shp_Mesh is the triangular mesh of the domain.

step 3:
Generating the streamlines using the Fortran code '02_generation_sl.f90'.
03_one_of_streamline.txt contains the coordinates of one of the generated streamlines.

step 4:
After generation of all the streamlines, the nodes should be transform into shapefile data type for interpolation.
The shapefile 'shp_streamlines' is an example data contains all the streamlines.

step 5:
The measured cross-sectional data is imposed on these streamlines to interpolate bathymetry. 
The final bathymetry of the example river can be interpolated using general GIS software such as ArcGIS or QGIS.
TIN_Reconstructed_Bathymetry is the result.
