
#include "Projections.h"
#include <cmath>

namespace GAASP {

using std::pow;

///////////////////////////////////////////////////////////////////////////////
//	Projection onto piecewise linear interpolant.
///////////////////////////////////////////////////////////////////////////////
double projection(double t, double x1, double x2, double tstart, double tend)
{
	double slope, nx;
	slope = (x2-x1)/(tend-tstart);
	nx = x1 + slope*(t-tstart);
	return nx;
}

} // namespace GAASP

