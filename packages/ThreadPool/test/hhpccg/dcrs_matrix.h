
#ifndef dcrs_matrix_h
#define dcrs_matrix_h

#include <tpi_vector.h>

struct distributed_crs_matrix {
  /* Global parallel */
  int   p_size ;
  int   p_rank ;
  int * p_recv_pc ; /* [np+1], span of received off-processor elements */
  int * p_send_pc ; /* [np+1], span of sent off-processor elements */
  int * p_send_id ; /* [send_pc[np]], indices of sent elements */

  /* Local and local parallel */
  int   n_local_column ; /* Number of local columns */
  int   n_local_row ;    /* Number of local rows */
  int   n_internal_row ; /* Number of local rows with internal columns */
  int * A_pc ;           /* Offsets into A_ia array for column indices */
  int * A_ia ;
  MATRIX_SCALAR * A_a ;
};

/*  1) communicate off-processor portions of input.
 *  2) apply: output = matrix * input ;
 *  3) return: dot( output , input );
 */
double dcrs_apply_and_dot( const struct distributed_crs_matrix * matrix ,
                           VECTOR_SCALAR * input ,
                           VECTOR_SCALAR * output ,
                           const int overlap_communication );

/*  1) communicate off-processor portions of input.
 *  2) apply: output = matrix * input ;
 */
void dcrs_apply( const struct distributed_crs_matrix * matrix ,
                 VECTOR_SCALAR * input ,
                 VECTOR_SCALAR * output );

#endif

