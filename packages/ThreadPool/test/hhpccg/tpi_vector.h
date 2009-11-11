
#ifndef tpi_vector_h
#define tpi_vector_h

#define VECTOR_SCALAR float
#define MATRIX_SCALAR float

void tpi_fill( int n , VECTOR_SCALAR alpha , VECTOR_SCALAR * x );

void tpi_scale( int n , const VECTOR_SCALAR alpha , VECTOR_SCALAR * x );

void tpi_copy( int n , const VECTOR_SCALAR * x , VECTOR_SCALAR * y );

void tpi_xpby( int n , const VECTOR_SCALAR * x ,
                             VECTOR_SCALAR beta  , VECTOR_SCALAR * y );

void tpi_axpy( int n , VECTOR_SCALAR alpha , const VECTOR_SCALAR * x ,
                                                   VECTOR_SCALAR * y );

void tpi_axpby( int n , VECTOR_SCALAR alpha , const VECTOR_SCALAR * x ,
                        VECTOR_SCALAR beta  ,       VECTOR_SCALAR * y );

double tpi_dot( int n , const VECTOR_SCALAR * x ,
                        const VECTOR_SCALAR * y );

void tpi_work_span( TPI_Work * const work , const int n ,
                    int * const iBeg , int * const iEnd );

#endif

