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

   // total number of vertices in the container with duplicates
   int containerNumberOfVerticesWithDups = 0;

   // initialize variables used in cell loop
   voronoicell_neighbor c;        // create a Voronoi cell with neighbor information
   std::vector<int> neigh;        // neighbors' information
   std::vector<int> f_vert;       // vertices indexes
   std::vector<int> f_order;      // number of vertices per face
   std::vector<double> v;

   // create container loop object
   c_loop_all cellLoop(con);

   // loop through cells
   if( cellLoop.start() ) {
      int counter = 0;
      do {
         if( con.compute_cell(c,cellLoop) ) {  // get Voronoi cell information
            std::cout << counter << std::endl;

            // gather information about the computed Voronoi cell
            cellLoop.pos(x,y,z);
            c.neighbors(neigh);
            c.face_vertices(f_vert);
            c.vertices(x,y,z,v);

            // print neighbors
            std::cout << "neighbors: ";
            for( std::vector<int>::const_iterator i = neigh.begin(); i != neigh.end(); ++i) {
               std::cout << *i << ' ';
            }
            std::cout << endl;

            // print vertices' indexes
            std::cout << "vertices: ";
            for( std::vector<int>::const_iterator i = f_vert.begin(); i != f_vert.end(); ++i) {
               std::cout << *i << ' ';
            }
            std::cout << endl;
            c.output_face_vertices();
            std::cout << endl;

            // print vertices' positions
            std::cout << "vertex position: ";
            for( std::vector<double>::const_iterator i = v.begin(); i != v.end(); ++i) {
               std::cout << *i << ' ';
            }
            std::cout << endl;
            std::cout << "number of vertices: " << (v.size()/3.0) << std::endl;
            c.output_vertices(x,y,z);
            std::cout << endl;

            // loop vertices and store their position; update container counter
            for( unsigned int i=0 ; i<(v.size()/3) ; i+=3 ) {
               points->InsertNextPoint(v[i], v[i+1], v[i+2]);
               ++containerNumberOfVerticesWithDups;
            }

            // vtk faces
            // [numberOfCellFaces, (numberOfPointsOfFace0, pointId0, pointId1, … ),
            //                     (numberOfPointsOfFace1, pointId0, pointId1, …), … ].
            // create faces ID list
            vtkSmartPointer<vtkIdList> vtkFaces = vtkSmartPointer<vtkIdList>::New();
            vtkFaces->InsertNextId(neigh.size());

            // loop over all faces of the Voronoi cell
            int face_counter = 0;
            for( unsigned int i=0 ; i<neigh.size() ; ++i ) {
               for( unsigned int j=0 ; j<f_vert.size() ; ++j ) {
                  vtkFaces->InsertNextId(f_vert[j]);
               }
               std::cout << "    " << face_counter << std::endl;
               face_counter++;
            } // end face loop

            counter++;
         } // end neighbor compute if
      } while(cellLoop.inc()); // end do loop
   } // end cell loop if


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
