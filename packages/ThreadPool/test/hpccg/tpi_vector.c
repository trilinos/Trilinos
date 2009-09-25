#include <stdio.h>

#include <stddef.h>

#include <TPI.h>
#include <tpi_vector.h>

/*--------------------------------------------------------------------*/

struct tpi_work_vector {
        double alpha ;
        double beta ;
  const double * x ;
  const double * y ;
        double * w ; 
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

  const double alpha = h->alpha ;
  double * const w = h->w ;

  int i , iEnd ;

  tpi_work_span( work , h->n , & i , & iEnd );

  for ( ; i < iEnd ; ++i ) { w[i] = alpha ; }
}

void tpi_fill( int n , double alpha , double * x )
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

  const double beta = h->beta ;
  double * const w = h->w ;

  int i , iEnd ;

  tpi_work_span( work , h->n , & i , & iEnd );

  for ( ; i < iEnd ; ++i ) { w[i] *= beta ; }
}

void tpi_scale( int n , const double alpha , double * x )
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

  const double * const x = h->x ;
  double * const w = h->w ;

  int i , iEnd ;

  tpi_work_span( work , h->n , & i , & iEnd );

  for ( ; i < iEnd ; ++i ) { w[i] = x[i] ; }
}

void tpi_copy( int n , const double * x , double * y )
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

  const double alpha = h->alpha ;
  const double beta  = h->beta ;
  const double * const x = h->x ;
  double * const w = h->w ;

  int i , iEnd ;

  tpi_work_span( work , h->n , & i , & iEnd );

  for ( ; i < iEnd ; ++i ) { w[i] = alpha * x[i] + beta * w[i] ; }
}

void tpi_axpby( int n , double alpha , const double * x ,
                        double beta  ,       double * y )
{
  struct tpi_work_vector tmp = { 0.0 , 0.0 , NULL , NULL , NULL , 0 };
  tmp.alpha = alpha ;
  tmp.beta  = beta ;
  tmp.x = x ;
  tmp.w = y ;
  tmp.n = n ;

  if ( 0.0 == alpha ) {
    TPI_Run_threads( tpi_work_scale , & tmp , 0 );
  }
  else if ( 0.0 == beta && 1.0 == alpha ) {
    TPI_Run_threads( tpi_work_copy , & tmp , 0 );
  }
  else {
    TPI_Run_threads( tpi_work_axpby , & tmp , 0 );
  }
}

/*--------------------------------------------------------------------*/

static void tpi_work_dot_partial( TPI_Work * work )
{
  const struct tpi_work_vector * const h =
    (struct tpi_work_vector *) work->info ;

  double * const s = (double *) work->reduce ;
  const double * const x = h->x ;
  const double * const y = h->y ;
  double tmp = *s ;
  int i , iEnd ;

  tpi_work_span( work , h->n , & i , & iEnd );

  for ( ; i < iEnd ; ++i ) { tmp += x[i] * y[i] ; }

  *s = tmp ;
}

static void tpi_work_dot_partial_self( TPI_Work * work )
{
  const struct tpi_work_vector * const h =
    (struct tpi_work_vector *) work->info ;

  double * const s = (double *) work->reduce ;
  const double * const x = h->x ;
  double tmp = *s ;

  int i , iEnd ;

  tpi_work_span( work , h->n , & i , & iEnd );

  for ( ; i < iEnd ; ++i ) { const double d = x[i] ; tmp += d * d ; }

  *s = tmp ;
}

static void tpi_work_dot_join( TPI_Work * work , void * src  )
{
  *((double *) ( work->reduce) ) += *((const double *) src);
}

static void tpi_work_dot_init( TPI_Work * work )
{
  *((double *) ( work->reduce) ) = 0 ;
}

double tpi_dot( int n , const double * x , const double * y )
{
  struct tpi_work_vector tmp = { 0.0 , 0.0 , NULL , NULL , NULL , 0 };
  double result = 0.0 ;
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
  const float  * A_a ;
  const double * x ;
        double * y ;
};

static void tpi_work_crs_matrix_apply( TPI_Work * work )
{
  const struct tpi_crs_matrix * const h =
    (struct tpi_crs_matrix *) work->info ;

  const int   * const A_pc = h->A_pc ;
  const int   * const A_ia = h->A_ia ;
  const float * const A_a  = h->A_a ;
  const double * const x = h->x ;

  const int nRow  = h->nRow ;
  const int chunk = ( nRow + work->count - 1 ) / work->count ;

  int row    = chunk * work->rank ;
  int rowEnd = chunk + row ;

  if ( nRow < rowEnd ) { rowEnd = nRow ; }

  {
    const int * const pc_end = A_pc + rowEnd ;
    const int *       pc     = A_pc + row ;
    double    *       y      = h->y + row ;

    for ( ; pc != pc_end ; ++pc , ++y ) {
      const int   *       ia    = A_ia + *pc ;
      const float *       a     = A_a  + *pc ;
      const float * const a_end = A_a  + pc[1] ;
      double tmp = 0 ;
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
  const float  * A_a ,
  const double * x ,
        double * y )
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


