
#ifndef UNITTEST_GENERATETENSORCOEFFICIENTS
#define UNITTEST_GENERATETENSORCOEFFICIENTS

namespace unit_test {

inline
double generate_matrix_coefficient( const unsigned nFEM ,
                                    const unsigned nStoch ,
                                    const unsigned iRowFEM ,
                                    const unsigned iColFEM ,
                                    const unsigned iStoch )
{
  const double A_fem = ( 10.0 + double(iRowFEM) / double(nFEM) ) +
                       (  5.0 + double(iColFEM) / double(nFEM) );

  const double A_stoch = ( 1.0 + double(iStoch) / double(nStoch) );

  return A_fem + A_stoch ;
}

inline
double generate_vector_coefficient( const unsigned nFEM ,
                                    const unsigned nStoch ,
                                    const unsigned iColFEM ,
                                    const unsigned iStoch )
{
  const double X_fem = 100.0 + double(iColFEM) / double(nFEM);
  const double X_stoch =  1.0 + double(iStoch) / double(nStoch);
  return X_fem + X_stoch ;
}

}

#endif /* #ifndef UNITTEST_GENERATETENSORCOEFFICIENTS */


