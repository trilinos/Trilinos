/*------------------------------------------------------------------------*/
/*                    TPI: Thread Pool Interface                          */
/*                Copyright (2008) Sandia Corporation                     */
/*                                                                        */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*                                                                        */
/* Redistribution and use in source and binary forms, with or without     */
/* modification, are permitted provided that the following conditions are */
/* met:                                                                   */
/*                                                                        */
/* 1. Redistributions of source code must retain the above copyright      */
/* notice, this list of conditions and the following disclaimer.          */
/*                                                                        */
/* 2. Redistributions in binary form must reproduce the above copyright   */
/* notice, this list of conditions and the following disclaimer in the    */
/* documentation and/or other materials provided with the distribution.   */
/*                                                                        */
/* 3. Neither the name of the Corporation nor the names of the            */
/* contributors may be used to endorse or promote products derived from   */
/* this software without specific prior written permission.               */
/*                                                                        */
/* THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY        */
/* EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE      */
/* IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR     */
/* PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE    */
/* CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,  */
/* EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,    */
/* PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR     */
/* PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF */
/* LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING   */
/* NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS     */
/* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.           */
/*------------------------------------------------------------------------*/


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

