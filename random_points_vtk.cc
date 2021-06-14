// Voronoi calculation example code
//
// Original author (voro++ part): Chris H. Rycroft (LBL / UC Berkeley)
// Mofication author (vtk part) : Paula C. Sanematsu (Syracuse University)

#include "voro++.hh"
#include <vtk-7.1/vtkCellArray.h>
#include <vtk-7.1/vtkCellData.h>
#include <vtk-7.1/vtkDataArray.h>
#include <vtk-7.1/vtkDataSetMapper.h>
#include <vtk-7.1/vtkDoubleArray.h>
#include <vtk-7.1/vtkIntArray.h>
#include <vtk-7.1/vtkIdList.h>
#include <vtk-7.1/vtkNamedColors.h>
#include <vtk-7.1/vtkPointData.h>
#include <vtk-7.1/vtkPoints.h>
#include <vtk-7.1/vtkPolyhedron.h>
#include <vtk-7.1/vtkProperty.h>
#include <vtk-7.1/vtkSmartPointer.h>
#include <vtk-7.1/vtkUnstructuredGrid.h>
#include <vtk-7.1/vtkVersion.h>
#include <vtk-7.1/vtkXMLUnstructuredGridWriter.h>
using namespace voro;

// Set up constants for the container geometry
const double x_min=-1,x_max=1;
const double y_min=-1,y_max=1;
const double z_min=-1,z_max=1;
const double cvol=(x_max-x_min)*(y_max-y_min)*(x_max-x_min);

// Set up the number of blocks that the container is divided into
const int n_x=6,n_y=6,n_z=6;

// Set the number of particles that are going to be randomly introduced
const int particles=20;

// This function returns a random double between 0 and 1
double rnd() {return double(rand())/RAND_MAX;}

int main() {
	int i;
	double x,y,z;

	// Create a container with the geometry given above, and make it
	// non-periodic in each of the three coordinates. Allocate space for
	// eight particles within each computational block
	container con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
			false,false,false,8);

	// Randomly add particles into the container
	for(i=0;i<particles;i++) {
		x=x_min+rnd()*(x_max-x_min);
		y=y_min+rnd()*(y_max-y_min);
		z=z_min+rnd()*(z_max-z_min);
		con.put(i,x,y,z);
	}

   // create unstructured grid and points (i.e. vertices)
   vtkSmartPointer<vtkUnstructuredGrid> uGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
   vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

	// Sum up the volumes, and check that this matches the container volume
	double vvol=con.sum_cell_volumes();
	printf("Container volume : %g\n"
	       "Voronoi volume   : %g\n"
	       "Difference       : %g\n",cvol,vvol,vvol-cvol);

	// Output the particle positions in gnuplot format
	con.draw_particles("random_points_p.gnu");

	// Output the Voronoi cells in gnuplot format
	con.draw_cells_gnuplot("random_points_v.gnu");
}
