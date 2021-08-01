#include "fitEllipsoid.h"
#include <iostream>
#include <Eigen/Dense>

/*
 * Here, the moment of inertia considers a point mass at each vertex. Other formulations of
 * moment of inertia exist, but here we use the simplest which seems to work very well for our
 * needs of defining an ellipsoid for each cell. For more details, see Dobrovolskis (1996) - Inertia
 * of any polyhedron.
*/
fitEllipsoid::fitEllipsoid(voro::voronoicell_neighbor& cell) {
   std::vector<double> myaxes;
   std::vector<double> v;

   // six terms of symmetric moment of intertia tensor
   double Ixx = 0;
   double Ixy = 0;
   double Ixz = 0;
   double Iyy = 0;
   double Iyz = 0;
   double Izz = 0;

   // gather information about the computed Voronoi cell
   cell.vertices(v);
   int nVertices = v.size();

   // loop vertices and store their position
   for( int i=0 ; i<nVertices ; i+=3 ) {
      double x = v[i];
      double y = v[i+1];
      double z = v[i+2];

      Ixx += y*y + z*z;
      Iyy += x*x + z*z;
      Izz += x*x + y*y;
      Ixy += x*y;
      Ixz += x*z;
      Iyz += y*z;  
   }

   Ixx /= double(nVertices);
   Iyy /= double(nVertices);
   Izz /= double(nVertices);
   Ixy /= double(nVertices);
   Ixz /= double(nVertices);
   Iyz /= double(nVertices);

   // create moment of inertia matrix to obtain eigenvalues
   Eigen::Matrix3d matrix;
   matrix <<  Ixx, -Ixy, -Ixz,
             -Ixy,  Iyy, -Iyz,
             -Ixz, -Iyz,  Izz;

   Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(matrix);
   Eigen::Vector3d evalues = eigensolver.eigenvalues();
   Eigen::Matrix3d evectors = eigensolver.eigenvectors();
   
   // store ellipsoid radii
   _majorRadius = sqrt(2.5*(evalues(1) + evalues(2) - evalues(0)));
   _interRadius = sqrt(2.5*(evalues(0) + evalues(2) - evalues(1)));
   _minorRadius = sqrt(2.5*(evalues(0) + evalues(1) - evalues(2)));

   // store ellipsoid axes
   for( int i=0 ; i < 3 ; ++i) {
      _majorAxis.push_back(evectors(i,0));
      _interAxis.push_back(evectors(i,0));
      _minorAxis.push_back(evectors(0,2));
   }

}
