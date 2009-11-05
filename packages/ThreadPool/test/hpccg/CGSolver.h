
#include <tpi_vector.h>

struct cgsolve_data {
  int             nRow ; 
  int           * A_pc ; 
  int           * A_ia ; 
  MATRIX_SCALAR * A_a ; 
  int             max_iter ; 
  int             print_iter ; 
  VECTOR_SCALAR   tolerance ; 

  int     np ; 
  int     ip ; 
  int   * recv_pc ; 
  int   * send_pc ; 
  int   * send_id ; 
}; 

void cgsolve_set_lhs( const struct cgsolve_data * data ,
                      const VECTOR_SCALAR * const x ,
                            VECTOR_SCALAR * const b );

void cgsolve( const struct cgsolve_data * data ,
              const VECTOR_SCALAR * const b ,
                    VECTOR_SCALAR * const x ,
                    int    * const iter_count ,
                    VECTOR_SCALAR * const norm_resid ,
                    double * const dt_mxv ,
                    double * const dt_axpby ,
                    double * const dt_dot );

