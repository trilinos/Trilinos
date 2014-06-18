//---------------------------------------------------------------------------//
// Ray_Tracer.hh
//---------------------------------------------------------------------------//

#ifndef mixed_language_src_Ray_Tracer_hh
#define mixed_language_src_Ray_Tracer_hh

#include <vector>
#include "Ray.hh"

namespace tribits_mixed
{

//---------------------------------------------------------------------------//

class Ray_Tracer
{
  public:
    // Typedefs.
    typedef std::vector<double> Vec_Dbl;
    typedef Ray::Space_Vector   Space_Vector;

  private:

    // Mesh dimensions.
    Vec_Dbl d_x, d_y, d_z;

    // Material opacities.
    Vec_Dbl d_sigma;

  public:
    // Constructor.
    Ray_Tracer(const Vec_Dbl &x, const Vec_Dbl &y, const Vec_Dbl &z, 
               const Vec_Dbl &sigma);

    // Destructor.
    ~Ray_Tracer();

    // Ray-trace.
    double ray_trace(const Space_Vector &begin, const Space_Vector &end);
};

} // namespace tribits_mixed

#endif

//---------------------------------------------------------------------------//
// end of Ray_Tracer.hh
//---------------------------------------------------------------------------//
