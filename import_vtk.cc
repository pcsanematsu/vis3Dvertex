// File import example code
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
#include <Eigen/Dense>
#include "fitEllipsoid.h"
using namespace voro;

// -------------------------------------------------------------------------------------------------------------
// functions to create cell data objects using vtkSmartPointer
vtkSmartPointer<vtkIntArray> createCellAttributeInt(int nComponents, int nTuples, const char* attName) {
   vtkSmartPointer<vtkIntArray> cellAttribute = vtkSmartPointer<vtkIntArray>::New();
   cellAttribute->SetNumberOfComponents(nComponents);
   cellAttribute->SetNumberOfTuples(nTuples);
   cellAttribute->SetName(attName);
   return cellAttribute;
}

vtkSmartPointer<vtkDoubleArray> createCellAttributeDouble(int nComponents, int nTuples, const char* attName) {
   vtkSmartPointer<vtkDoubleArray> cellAttribute = vtkSmartPointer<vtkDoubleArray>::New();
   cellAttribute->SetNumberOfComponents(nComponents);
   cellAttribute->SetNumberOfTuples(nTuples);
   cellAttribute->SetName(attName);
   return cellAttribute;
}
// -------------------------------------------------------------------------------------------------------------


// Set up constants for the container geometry
const double x_min=-5,x_max=5;
const double y_min=-5,y_max=5;
const double z_min=0,z_max=10;

// Set up the number of blocks that the container is divided into
const int n_x=6,n_y=6,n_z=6;

// Number of particles in the pack_ten_cube container
const int particles = 1000;

int main() {
	double x,y,z;   // particle position
   double tempArray[3];

	// Create a container with the geometry given above, and make it
	// non-periodic in each of the three coordinates. Allocate space for
	// eight particles within each computational block
	container con(x_min,x_max,y_min,y_max,z_min,z_max,n_x,n_y,n_z,
			false,false,false,8);

	//Randomly add particles into the container
	con.import("pack_ten_cube");

   // ------------------------------------- VTK declarations ---------------------------------------
   // create unstructured grid and points (i.e. vertices)
   vtkSmartPointer<vtkUnstructuredGrid> uGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
   vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

   // create cell attributes
   vtkSmartPointer<vtkIntArray> cellID = createCellAttributeInt(1, particles, "cellID");
   vtkSmartPointer<vtkDoubleArray> cellVolume = createCellAttributeDouble(1, particles, "cellVolume");
   vtkSmartPointer<vtkDoubleArray> cellPosition = createCellAttributeDouble(3, particles, "cellPosition");
   vtkSmartPointer<vtkDoubleArray> cellMajorRadius = createCellAttributeDouble(1, particles, "cellMajorRadius");
   vtkSmartPointer<vtkDoubleArray> cellMajorAxis = createCellAttributeDouble(3, particles, "cellMajorAxis");
   // ---------------------------------- end VTK declarations --------------------------------------

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

            // gather information about the computed Voronoi cell
            cellLoop.pos(x,y,z);
            c.neighbors(neigh);
            c.face_vertices(f_vert);
            c.vertices(x,y,z,v);

            // loop vertices and store their position
            for( unsigned int i=0 ; i<v.size() ; i+=3 ) {
               points->InsertNextPoint(v[i], v[i+1], v[i+2]);
            }

            // vtk faces
            // [numberOfCellFaces, (numberOfPointsOfFace0, pointId0, pointId1, … ),
            //                     (numberOfPointsOfFace1, pointId0, pointId1, …), … ].
            // create faces ID list
            vtkSmartPointer<vtkIdList> vtkFaces = vtkSmartPointer<vtkIdList>::New();
            vtkFaces->InsertNextId(neigh.size());

            // number of vertices in current cell
            int numberOfVertices = v.size()/3;

            // update total number of vertices in container
            containerNumberOfVerticesWithDups += numberOfVertices;

            // update starting index for the current cell in the container
            int containerVertexStartIndex = containerNumberOfVerticesWithDups - numberOfVertices;

            // loop over all faces of the Voronoi cell and populate vtkFaces with
            // numberOfVerticesPerFace and their vertex indices
            int j,k=0;
            int numberOfVerticesPerFace;
            while( (unsigned int)k<f_vert.size() ) {
               numberOfVerticesPerFace = f_vert[k++];
               vtkFaces->InsertNextId(numberOfVerticesPerFace);  // number of vertices in 1 face

               j = k+numberOfVerticesPerFace;
               while( k<j ) {
                  int containerIndex = f_vert[k++] + containerVertexStartIndex;
                  vtkFaces->InsertNextId(containerIndex);
               } // end single face loop
            } // end vertices loop

            // add cell to unstructure grid
            uGrid->InsertNextCell(VTK_POLYHEDRON,vtkFaces);

            // add attributes to cell
            cellID->InsertValue(counter, counter);
            cellVolume->InsertValue(counter, c.volume());

            // store cell positions into double array so they can be added to vtk array
            tempArray[0] = x; tempArray[1] = y; tempArray[2] = z;
            cellPosition->InsertTuple(counter, tempArray);

            // fit ellipsoid
            fitEllipsoid ellipsoid = fitEllipsoid(c);

            // get major radius and insert into vtk object
            cellMajorRadius->InsertValue(counter, ellipsoid.majorRadius());

            // get major axis, store in double array, and insert into vtk object
            std::vector<double> ma = ellipsoid.majorAxis();
            tempArray[0] = ma[0]; tempArray[1] = ma[1]; tempArray[2] = ma[2];
            cellMajorAxis->InsertTuple(counter, tempArray);

            counter++;
         } // end neighbor compute if
      } while(cellLoop.inc()); // end do loop

   } // end cell loop if

   // add attributes to cellData
   uGrid->GetCellData()->AddArray(cellID);
   uGrid->GetCellData()->AddArray(cellVolume);
   uGrid->GetCellData()->AddArray(cellPosition);
   uGrid->GetCellData()->AddArray(cellMajorRadius);
   uGrid->GetCellData()->AddArray(cellMajorAxis);

   // add points to unstructured grid
   uGrid->SetPoints(points);

   // output unstructured grid
   std::cout << "Writing .vtu file..." << std::endl;
   vtkSmartPointer<vtkXMLUnstructuredGridWriter> vtkWriter = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
   vtkWriter->SetInputData(uGrid);
   vtkWriter->SetFileName("pack_ten_cube.vtu");
   vtkWriter->SetDataModeToBinary();  // much smaller files and faster
   vtkWriter->Update();

	// Save the Voronoi network of all the particles to text files
	// in gnuplot and POV-Ray formats
	con.draw_cells_gnuplot("pack_ten_cube.gnu");
	con.draw_cells_pov("pack_ten_cube_v.pov");

	// Output the particles in POV-Ray format
	con.draw_particles_pov("pack_ten_cube_p.pov");
}
