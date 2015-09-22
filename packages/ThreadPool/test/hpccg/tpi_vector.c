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

#include <stdio.h>

#include <stddef.h>

#include <TPI.h>
#include <tpi_vector.h>

/*--------------------------------------------------------------------*/

struct tpi_work_vector {
        VECTOR_SCALAR alpha ;
        VECTOR_SCALAR beta ;
  const VECTOR_SCALAR * x ;
  const VECTOR_SCALAR * y ;
        VECTOR_SCALAR * w ; 
        int  n ;
};

static void tpi_work_span( TPI_Work * const work , const int n ,
                           int * const iBeg , int * const iEnd )
{
  const int chunk = ( n + work->count - 1 ) / work->count ;
  const int i_end = chunk + ( *iBeg = chunk * work->rank );

  *iEnd = n < i_end ? n : i_end ;
}

/*--------------------------------------------------------------------*/

static void tpi_work_fill( TPI_Work * work )
{
  const struct tpi_work_vector * const h =
    (struct tpi_work_vector *) work->info ;

  const VECTOR_SCALAR alpha = h->alpha ;
  VECTOR_SCALAR * const w = h->w ;

  int i , iEnd ;

  tpi_work_span( work , h->n , & i , & iEnd );

  for ( ; i < iEnd ; ++i ) { w[i] = alpha ; }
}

void tpi_fill( int n , VECTOR_SCALAR alpha , VECTOR_SCALAR * x )
{
  struct tpi_work_vector tmp = { 0.0 , 0.0 , NULL , NULL , NULL , 0 };
  tmp.alpha = alpha ;
  tmp.w = x ;
  tmp.n = n ;
  TPI_Run_threads( tpi_work_fill , & tmp , 0 );
}

/*--------------------------------------------------------------------*/

static void tpi_work_scale( TPI_Work * work )
{
  const struct tpi_work_vector * const h =
    (struct tpi_work_vector *) work->info ;

  const VECTOR_SCALAR beta = h->beta ;
  VECTOR_SCALAR * const w = h->w ;

  int i , iEnd ;

  tpi_work_span( work , h->n , & i , & iEnd );

  for ( ; i < iEnd ; ++i ) { w[i] *= beta ; }
}

void tpi_scale( int n , const VECTOR_SCALAR alpha , VECTOR_SCALAR * x )
{
  struct tpi_work_vector tmp = { 0.0 , 0.0 , NULL , NULL , NULL , 0 };
  tmp.alpha = alpha ;
  tmp.w = x ;
  tmp.n = n ;
  TPI_Run_threads( tpi_work_scale , & tmp , 0 );
}

/*--------------------------------------------------------------------*/

static void tpi_work_copy( TPI_Work * work )
{
  const struct tpi_work_vector * const h =
    (struct tpi_work_vector *) work->info ;

  const VECTOR_SCALAR * const x = h->x ;
  VECTOR_SCALAR * const w = h->w ;

  int i , iEnd ;

  tpi_work_span( work , h->n , & i , & iEnd );

  for ( ; i < iEnd ; ++i ) { w[i] = x[i] ; }
}

void tpi_copy( int n , const VECTOR_SCALAR * x , VECTOR_SCALAR * y )
{
  struct tpi_work_vector tmp = { 0.0 , 0.0 , NULL , NULL , NULL , 0 };
  tmp.x = x ;
  tmp.w = y ;
  tmp.n = n ;
  TPI_Run_threads( tpi_work_copy , & tmp , 0 );
}

/*--------------------------------------------------------------------*/

static void tpi_work_axpby( TPI_Work * work )
{
  const struct tpi_work_vector * const h =
    (struct tpi_work_vector *) work->info ;

  const VECTOR_SCALAR alpha = h->alpha ;
  const VECTOR_SCALAR beta  = h->beta ;
  const VECTOR_SCALAR * const x = h->x ;
  VECTOR_SCALAR * const w = h->w ;

  int i , iEnd ;

  tpi_work_span( work , h->n , & i , & iEnd );

  for ( ; i < iEnd ; ++i ) { w[i] = alpha * x[i] + beta * w[i] ; }
}

void tpi_axpby( int n , VECTOR_SCALAR alpha , const VECTOR_SCALAR * x ,
                        VECTOR_SCALAR beta  ,       VECTOR_SCALAR * y )
{
  struct tpi_work_vector tmp = { 0.0 , 0.0 , NULL , NULL , NULL , 0 };
  tmp.alpha = alpha ;
  tmp.beta  = beta ;
  tmp.x = x ;
  tmp.w = y ;
  tmp.n = n ;

  TPI_Run_threads( tpi_work_axpby , & tmp , 0 );
}

/*--------------------------------------------------------------------*/

static void tpi_work_dot_partial( TPI_Work * work )
{
  const struct tpi_work_vector * const h =
    (struct tpi_work_vector *) work->info ;

  VECTOR_SCALAR * const s = (VECTOR_SCALAR *) work->reduce ;
  const VECTOR_SCALAR * const x = h->x ;
  const VECTOR_SCALAR * const y = h->y ;
  VECTOR_SCALAR tmp = *s ;
  int i , iEnd ;

  tpi_work_span( work , h->n , & i , & iEnd );

  for ( ; i < iEnd ; ++i ) { tmp += x[i] * y[i] ; }

  *s = tmp ;
}

static void tpi_work_dot_partial_self( TPI_Work * work )
{
  const struct tpi_work_vector * const h =
    (struct tpi_work_vector *) work->info ;

  VECTOR_SCALAR * const s = (VECTOR_SCALAR *) work->reduce ;
  const VECTOR_SCALAR * const x = h->x ;
  VECTOR_SCALAR tmp = *s ;

  int i , iEnd ;

  tpi_work_span( work , h->n , & i , & iEnd );

  for ( ; i < iEnd ; ++i ) { const VECTOR_SCALAR d = x[i] ; tmp += d * d ; }

  *s = tmp ;
}

static void tpi_work_dot_join( TPI_Work * work , const void * src  )
{
  *((VECTOR_SCALAR *) ( work->reduce) ) += *((const VECTOR_SCALAR *) src);
}

static void tpi_work_dot_init( TPI_Work * work )
{
  *((VECTOR_SCALAR *) ( work->reduce) ) = 0 ;
}

VECTOR_SCALAR tpi_dot( int n , const VECTOR_SCALAR * x , const VECTOR_SCALAR * y )
{
  struct tpi_work_vector tmp = { 0.0 , 0.0 , NULL , NULL , NULL , 0 };
  VECTOR_SCALAR result = 0.0 ;
  tmp.x = x ;
  tmp.y = y ;
  tmp.n = n ;
  if ( x != y ) {
    TPI_Run_threads_reduce( tpi_work_dot_partial , & tmp ,
                            tpi_work_dot_join , tpi_work_dot_init ,
                            sizeof(result) , & result );
  }
  else {
    TPI_Run_threads_reduce( tpi_work_dot_partial_self , & tmp ,
                            tpi_work_dot_join , tpi_work_dot_init ,
                            sizeof(result) , & result );
  }
  return result ;
}

/*--------------------------------------------------------------------*/

struct tpi_crs_matrix {
        int      nRow ;
  const int    * A_pc ;
  const int    * A_ia ;
  const MATRIX_SCALAR * A_a ;
  const VECTOR_SCALAR * x ;
        VECTOR_SCALAR * y ;
};

static void tpi_work_crs_matrix_apply( TPI_Work * work )
{
  const struct tpi_crs_matrix * const h =
    (struct tpi_crs_matrix *) work->info ;

  const int   * const A_pc = h->A_pc ;
  const int   * const A_ia = h->A_ia ;
  const MATRIX_SCALAR * const A_a  = h->A_a ;
  const VECTOR_SCALAR * const x = h->x ;

  const int nRow  = h->nRow ;
  const int chunk = ( nRow + work->count - 1 ) / work->count ;

  int row    = chunk * work->rank ;
  int rowEnd = chunk + row ;

  if ( nRow < rowEnd ) { rowEnd = nRow ; }

  {
    const int * const pc_end = A_pc + rowEnd ;
    const int *       pc     = A_pc + row ;
    VECTOR_SCALAR *   y      = h->y + row ;

    for ( ; pc != pc_end ; ++pc , ++y ) {
      const int   *       ia    = A_ia + *pc ;
      const MATRIX_SCALAR *       a     = A_a  + *pc ;
      const MATRIX_SCALAR * const a_end = A_a  + pc[1] ;
      VECTOR_SCALAR tmp = 0 ;
      for ( ; a != a_end ; ++a , ++ia ) {
        tmp += *a * x[ *ia ];
      }
      *y = tmp ;
    }
  }
}

/*--------------------------------------------------------------------*/

void tpi_crs_matrix_apply(
  const int      nRow ,
  const int    * A_pc ,
  const int    * A_ia ,
  const MATRIX_SCALAR * A_a ,
  const VECTOR_SCALAR * x ,
        VECTOR_SCALAR * y )
{
  struct tpi_crs_matrix h = { 0 , NULL , NULL , NULL , NULL , NULL };
  h.nRow = nRow ;
  h.A_pc = A_pc ;
  h.A_ia = A_ia ;
  h.A_a  = A_a ;
  h.x    = x ;
  h.y    = y ;
  TPI_Run_threads( tpi_work_crs_matrix_apply , & h , 0 );
}


