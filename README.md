# vis3Dvertex
This repository show how to "convert" vertex, face, and cell data from vertex-like models (vertex or Voronoi) into VTK polyhedral unstructured meshes. I
modify three examples from voro++ library to show how this "conversion" is performed. In addition, it shows how to output binary or ASCII files in `.vtu` file format that can be read in [ParaView](https://www.paraview.org/). [`random_points_vtk.cc`](random_points_vtk.cc) shows how to output a series of files as a timeseries such that the user can create movies/animations.

The output files of the examples below are provided in the folder `output`. To create figures 4 and 6 of the manuscript, I provided two ParaView state files (`.pvsm`):
* `threshold_2Dcross_section_filters.pvsm`: Figure 4
* `calculator_filter.pvsm`: Figure 6
For more details on how to implement, refer to manuscript under [Citation](#citation).

## Citation
Please, if you use this in your publication, cite the following publication

## Dependencies
The examples require the libraries
* [voro++](http://math.lbl.gov/voro++/) (see note below about voro++'s dependencies)
* [Eigen3](https://eigen.tuxfamily.org/index.php?title=Main_Page)
* [VTK](https://vtk.org/)

## Singularity container
Examples were written and tested in a Linux machine. Codes were compiled and run using a Singularity container where the host machine had Singularity 2.6.0-dist installed. You can directly download the Singularity container (create link here) or use the Singularity deifinition file `voroVTKPy.def` from this repository to create on your container. 

Singularity containers allows reproducibility of computational results in different machines. It also allows one to avoid going through the pain of installing all of the libraries' dependencies. voro++ has a few dependencies, so to simplify things, it is probably easier to the Singularity container.

## Compilation
1. Clone this repository
2. Then `cd vis3Dvertex`
3. Compile with the command `make all`

## Examples

### `cell_statistics_vtk.cc`
This code was based on `voro++-0.4.6/examples/custom/cell_statistics.cc`. 

This example shows how to create convert cell, face, and vertex data of a single cell into VTK polyhedral unstructured grid and outputs a \*.vtu file to be opened in ParaView.

### `random_points_vtk.cc`
This code was based on `voro++-0.4.6/examples/basic/random_points.cc`. 

In addition to converting cell, face, and vextex data into a VTK polyhedral unstructured mesh, it creates the following scalar VTK `CellData`:
* `cellID`: the index of each cell
* `cellFaces`: the number of faces (i.e. neighbors) of each cell
* `cellVolume`: the volume of each cell
* `cellSurfaceArea`: the surface area of the cell

It also creates a timeseries file `timeseries.pvd` such that the user can see the time evolution of a simulation in ParaView and then export movies.

### `import_vtk.cc`
This code was based on `voro++-0.4.6/examples/basic/import.cc`.

In addition to some of the items included in the examples above, this examples shows how to add vectorial `CellData` in the polyhedral unstructured mesh. The following `CellData` are included:
* `cellID` (scalar): the index of each cell
* `cellVolume` (scalar): the volume of each cell
* `cellPosition` (vectorial): the position of the cell center in global coordinates
* `cellMajorRadius` (scalar): the radius of the major axis of the ellipsoid that was fitted for the a polyhedral cell
* `cellMajorAxis` (vectorial): the vector that describes the major axis of the ellipsoid that was fitted for the polyhedral cell; vector is in local (the cell's) coordinate system.

### `fitEllipsoid.cc`
This is not really an example of using the VTK library. It shows how to calculate the point-mass moment of inertia of a polyhedron [Dobrovolskis, A. R. (1996). "Inertia of Any Polyhedron." Icarus 124(2): 698-704.](https://doi.org/10.1006/icar.1996.0243). The user can access the major, intermediate, and minor radii and axes with the getter functions:
* `majorRadius()`, `majorAxis()`
* `interRadius()`, `interAxis()`
* `minorRadius()`, `minorAxis()`



