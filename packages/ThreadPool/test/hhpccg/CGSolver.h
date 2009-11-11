
#ifndef CGSolver_h
#define CGSolver_h

#include <tpi_vector.h>
#include <dcrs_matrix.h>

/*--------------------------------------------------------------------*/

void cgsolve_set_lhs( const struct distributed_crs_matrix * matrix ,
                      const VECTOR_SCALAR * const x ,
                            VECTOR_SCALAR * const b );

/* Solve with fused loops */
void cgsolve( const struct distributed_crs_matrix * matrix ,
              const VECTOR_SCALAR * const b ,
                    VECTOR_SCALAR * const x ,
              const VECTOR_SCALAR tolerance ,
              const int max_iter ,
              const int print_iter ,
                    int    * const iter_count ,
                    VECTOR_SCALAR * const norm_resid ,
                    double * const solve_dt );

/* Solve with blas-like calls */
void cgsolve_blas( const struct distributed_crs_matrix * matrix ,
                   const VECTOR_SCALAR * const b ,
                         VECTOR_SCALAR * const x ,
                   const VECTOR_SCALAR tolerance ,
                   const int max_iter ,
                   const int print_iter ,
                         int    * const iter_count ,
                         VECTOR_SCALAR * const norm_resid ,
                         double * const solve_dt );

/*--------------------------------------------------------------------*/

#endif

