//---------------------------------------------------------------------------//
// Ray_Tracer.cc
//---------------------------------------------------------------------------//

#include <cmath>
#include <algorithm>
#include <numeric>

#include "MixedLang_config.h"
#include "Ray_Tracer.hh"

//---------------------------------------------------------------------------//
// Prototypes for FORTRAN function kernels.
//---------------------------------------------------------------------------//

#define BLOCK_TRACER FC_FUNC(block_tracer, BLOCK_TRACER)
#define SETUP_RAY_TRACING FC_FUNC(setup_ray_tracing, SETUP_RAY_TRACING)
#define CLEANUP_RAY_TRACING FC_FUNC(cleanup_ray_tracing, CLEANUP_RAY_TRACING)
#define SET_DIMENSIONS FC_FUNC(set_dimensions, SET_DIMENSIONS)

extern "C"
{
    
    void BLOCK_TRACER(const double *x, const double *y, const double *z,
                      const double *sigma, double *pos, const double *end, 
                      int *ijk, double *tau, int *terminator);

    void SETUP_RAY_TRACING();
    
    void CLEANUP_RAY_TRACING();

    void SET_DIMENSIONS(int *im, int *jm, int *km);
}

namespace tribits_mixed
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR AND DESTRUCTOR
//---------------------------------------------------------------------------//

Ray_Tracer::Ray_Tracer(const Vec_Dbl &x,
                       const Vec_Dbl &y,
                       const Vec_Dbl &z,
                       const Vec_Dbl &sigma)
    : d_x(x)
    , d_y(y)
    , d_z(z)
    , d_sigma(sigma)
{
    // set the dimensions
    int i = d_x.size() - 1;
    int j = d_y.size() - 1;
    int k = d_z.size() - 1;
    SET_DIMENSIONS(&i, &j, &k);

    // setup the ray-tracer
    SETUP_RAY_TRACING();
}

//---------------------------------------------------------------------------//

Ray_Tracer::~Ray_Tracer()
{
    CLEANUP_RAY_TRACING();
}

//---------------------------------------------------------------------------//
// PUBLIC FUNCTIONS
//---------------------------------------------------------------------------//

double Ray_Tracer::ray_trace(const Space_Vector &begin,
                             const Space_Vector &end)
{
    // make a ray
    Ray ray;
    ray.p = begin;

    // calculate index of ray
    int i,j,k;
    Vec_Dbl::const_iterator ijk;

    // x
    ijk = std::lower_bound(d_x.begin(), d_x.end(), ray.p[0]);
    i   = ijk - d_x.begin() - 1;

    // y
    ijk = std::lower_bound(d_y.begin(), d_y.end(), ray.p[1]);
    j   = ijk - d_y.begin() - 1;

    // z
    ijk = std::lower_bound(d_z.begin(), d_z.end(), ray.p[2]);
    k   = ijk - d_z.begin() - 1;
    
    ray.index[0] = i;
    ray.index[1] = j;
    ray.index[2] = k;

    // reset the tau
    ray.tau = 0.0;

    // termination condition
    int terminator = 0;

    // trace through the block using the FORTRAN kernel
    BLOCK_TRACER(&d_x[0], &d_y[0], &d_z[0], &d_sigma[0], &ray.p[0], &end[0],
                 &ray.index[0], &ray.tau, &terminator);

    // add up tau
    return ray.tau;
}

} // namespace tribits_mixed

//---------------------------------------------------------------------------//
// end of Ray_Tracer.cc
//---------------------------------------------------------------------------//

