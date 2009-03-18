
void tpi_fill( int n , double alpha , double * x );

void tpi_scale( int n , const double alpha , double * x );

void tpi_copy( int n , const double * x , double * y );

void tpi_axpby( int n , double alpha , const double * x ,
                        double beta  ,       double * y );

double tpi_dot( int n , const double * x , const double * y );

void tpi_crs_matrix_apply(
  const int      nRow ,
  const int    * A_pc ,
  const int    * A_ia ,
  const float  * A_a ,
  const double * x ,
        double * y );

