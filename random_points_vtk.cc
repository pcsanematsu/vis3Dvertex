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

bool is_big_endian(void) {
    union {
        uint32_t i;
        char c[4];
    } bint = {0x01020304};

    return bint.c[0] == 1;
}

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
   ofstream vtkTimeseries ("timeseries.pvd");

   // determine system's endianness that is necessary for timeseries file
   std::string endianness;
   if( is_big_endian() ) {
      endianness = "BigEndian";
   }
   else {
      endianness = "LittleEndian";
   }
   std::cout << "Endianness: " << endianness << std::endl;

   // header of time series file
   vtkTimeseries << "<?xml version=\"1.0\"?>" << endl;
   vtkTimeseries << "<VTKFile type=\"Collection\" version=\"0.1\"" << endl;
   vtkTimeseries << "         byte_order=\"" << endianness << "\"" << endl;
   vtkTimeseries << "         compressor=\"vtkZLibDataCompressor\">" << endl;
   vtkTimeseries << "  <Collection>" << endl;

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

   // ------------------------------------- VTK declarations ---------------------------------------
   // create unstructured grid and points (i.e. vertices)
   vtkSmartPointer<vtkUnstructuredGrid> uGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
   vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

   // create cell attributes
   vtkSmartPointer<vtkIntArray> cellID = createCellAttributeInt(1, particles, "cellID");
   vtkSmartPointer<vtkIntArray> cellFaces = createCellAttributeInt(1, particles, "cellFaces");
   vtkSmartPointer<vtkDoubleArray> cellVolume = createCellAttributeDouble(1, particles, "cellVolume");
   vtkSmartPointer<vtkDoubleArray> cellSurfaceArea = createCellAttributeDouble(1, particles, "cellSurfaceArea");
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
            std::cout << "cell: " << counter << std::endl;

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
            cellFaces->InsertValue(counter, neigh.size());
            cellSurfaceArea->InsertValue(counter, c.surface_area());


            counter++;
         } // end neighbor compute if
      } while(cellLoop.inc()); // end do loop

   } // end cell loop if

   // add attributes to cellData
   uGrid->GetCellData()->AddArray(cellID);
   uGrid->GetCellData()->AddArray(cellVolume);
   uGrid->GetCellData()->AddArray(cellFaces);
   uGrid->GetCellData()->AddArray(cellSurfaceArea);

   // add points to unstructured grid
   uGrid->SetPoints(points);

   // output unstructured grid
   std::cout << "Writing .vtu file..." << std::endl;
   vtkSmartPointer<vtkXMLUnstructuredGridWriter> vtkWriter = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
   vtkWriter->SetInputData(uGrid);
   vtkWriter->SetFileName("random_points.vtu");
   vtkWriter->SetDataModeToAscii();
   //vtkWriter->SetDataModeToBinary();  // much smaller files and faster
   vtkWriter->Update();

   // end of timeseries file
   vtkTimeseries << "  </Collection>" << endl;
   vtkTimeseries << "</VTKFile>" << endl;
   vtkTimeseries.close();

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
