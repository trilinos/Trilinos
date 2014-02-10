//---------------------------------------------------------------------------//
// Ray.hh
//---------------------------------------------------------------------------//

#ifndef mixed_language_src_Ray_hh
#define mixed_language_src_Ray_hh

#include "Vector_Lite.hh"

namespace tribits_mixed
{

//---------------------------------------------------------------------------//

struct Ray 
{
    // 3-Int Vector (used for 3-dimension indices in XYZ space).
    typedef Vector_Lite<int, 3> Indices;
    
    // 3-Vector defining a point.
    typedef Vector_Lite<double, 3> Space_Vector;

    // Mean-free-path along a ray. 
    double tau;

    // Current position of ray.
    Space_Vector p;

    // Current direction of ray.
    Space_Vector dir;

    // (I,J,K) cell-indices of ray.
    Indices index;
};

} // namespace tribits_mixed

#endif

//---------------------------------------------------------------------------//
// end of Ray.hh
//---------------------------------------------------------------------------//

