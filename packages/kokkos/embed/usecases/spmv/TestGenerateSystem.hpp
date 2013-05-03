
#ifndef TESTGENERATE_GRAPH_HPP
#define TESTGENERATE_GRAPH_HPP

#include <cstddef>
#include <vector>

namespace Test {

template< typename IntType >
inline
IntType map_fem_graph_coord( const IntType & N ,
                             const IntType & i ,
                             const IntType & j ,
                             const IntType & k )
{ 
  return k + N * ( j + N * i );
} 

inline
size_t generate_fem_graph( size_t N , std::vector< std::vector<size_t> > & graph )
{
  graph.resize( N * N * N , std::vector<size_t>() );

  size_t total = 0 ;

  for ( int i = 0 ; i < (int) N ; ++i ) {
  for ( int j = 0 ; j < (int) N ; ++j ) {
  for ( int k = 0 ; k < (int) N ; ++k ) {

    const size_t row = map_fem_graph_coord((int)N,i,j,k);

    graph[row].reserve(27);

    for ( int ii = -1 ; ii < 2 ; ++ii ) {
    for ( int jj = -1 ; jj < 2 ; ++jj ) {
    for ( int kk = -1 ; kk < 2 ; ++kk ) {
      if ( 0 <= i + ii && i + ii < (int) N &&
           0 <= j + jj && j + jj < (int) N &&
           0 <= k + kk && k + kk < (int) N ) {
        size_t col = map_fem_graph_coord((int)N,i+ii,j+jj,k+kk);

        graph[row].push_back(col);
      }
    }}}
    total += graph[row].size();
  }}}

  return total ;
}

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

#endif /* #ifndef TESTGENERATE_GRAPH_HPP */

