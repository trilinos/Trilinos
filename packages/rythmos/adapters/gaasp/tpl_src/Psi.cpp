#include "GAdjointSolve.h"

namespace GAASP {

double* GAdjointSolve::psi()
{
    double * p = new double[ sysDim ];
	switch ( ErrorFlag )
	{
	    case Error_End:
		for ( int i = 0; i < sysDim; ++i )
			p[i] = 0.0;
		break;

	    case Error_Avg:
		for ( int i = 0; i < sysDim; ++i )
			p[i] = 1.0/endTime;
		break;
	}
	return p;
}

} // namespace GAASP

