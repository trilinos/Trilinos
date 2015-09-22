#include <cmath>
namespace GAASP {
    using std::pow;
//*****************************************************************************
//	function: double* l2Norm(double * const yin, int const n)	
//	Written by:	Jeff Sandelin
//	Computes the L2 norm of the given vector.
//
//	Input:	yin		-	input vector
//		n		-	size of the vector
//	Output:	l2		-	norm
//*****************************************************************************
double l2Norm ( double * const yin, int const n )
{
	double l2 = 0.0;

	for( int i = 0; i < n; ++i )
	{
		l2 += pow ( yin[i], 2 );
	}

	l2 = pow ( l2, 0.5 );
	return l2;
}

} // namespace GAASP

