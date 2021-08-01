#ifndef FITELLIPSOID_H    // To make sure you don't declare the function more than once by including the header multiple times.
#define FITELLIPSOID_H

#include "voro++.hh"
#include <vector>

class fitEllipsoid {
public:
   fitEllipsoid(voro::voronoicell_neighbor& cell);
   double majorRadius() const { return _majorRadius; }
   double interRadius() const { return _interRadius; }
   double minorRadius() const { return _minorRadius; }
   std::vector<double> majorAxis() const { return _majorAxis; }
   std::vector<double> interAxis() const { return _interAxis; }
   std::vector<double> minorAxis() const { return _minorAxis; }

private:
   double _minorRadius, _interRadius, _majorRadius;
   std::vector<double> _minorAxis, _interAxis, _majorAxis;
};

#endif
