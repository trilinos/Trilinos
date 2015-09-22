#ifndef INC_Projections_h
#define INC_Projections_h


#include <cmath>

namespace GAASP {

///////////////////////////////////////////////////////////////////////////////
//	Projection onto piecewise linear interpolant.
///////////////////////////////////////////////////////////////////////////////
double projection(double t, double x1, double x2, double tstart, double tend);

} // namespace GAASP

#endif // INC_Projections_h

