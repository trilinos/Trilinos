
#include <ThreadPool_config.h>

#ifndef tpi_vector_h
#define tpi_vector_h

#define VECTOR_SCALAR float
#define MATRIX_SCALAR float

void tpi_fill( int n , VECTOR_SCALAR alpha , VECTOR_SCALAR * x );

void tpi_scale( int n , const VECTOR_SCALAR alpha , VECTOR_SCALAR * x );

void tpi_copy( int n , const VECTOR_SCALAR * x , VECTOR_SCALAR * y );

void tpi_axpby( int n , VECTOR_SCALAR alpha , const VECTOR_SCALAR * x ,
                        VECTOR_SCALAR beta  ,       VECTOR_SCALAR * y );

VECTOR_SCALAR tpi_dot( int n , const VECTOR_SCALAR * x ,
                               const VECTOR_SCALAR * y );

void tpi_crs_matrix_apply(
  const int             nRow ,
  const int           * A_pc ,
  const int           * A_ia ,
  const MATRIX_SCALAR * A_a ,
  const VECTOR_SCALAR * x ,
        VECTOR_SCALAR * y );

#endif

