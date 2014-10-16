
#ifndef TESTGENERATE_GRAPH_HPP
#define TESTGENERATE_GRAPH_HPP

#include <cstddef>
#include <vector>
#include <algorithm>
#include <Kokkos_Core.hpp>
#include <Kokkos_Array.hpp>
#include <TestCrsMatrix.hpp>

namespace Test {


template< class GraphType , class MeshType >
GraphType create_graph_from_mesh( const MeshType & mesh )
{
  typedef GraphType                         graph_type ;
  typedef MeshType                          mesh_type ;
  typedef typename graph_type::device_type  device_type ;
  typedef typename device_type::size_type   size_type  ;

  static const unsigned ElemNodeCount = mesh_type::element_node_count ;

  typename mesh_type::node_elem_ids_type::HostMirror
    node_elem_ids = create_mirror( mesh.node_elem_ids );

  typename mesh_type::elem_node_ids_type::HostMirror
    elem_node_ids = create_mirror( mesh.elem_node_ids );

  deep_copy( elem_node_ids , mesh.elem_node_ids );
  deep_copy( node_elem_ids.entries , mesh.node_elem_ids.entries );

  const size_t owned_node = mesh.parallel_data_map.count_owned ;

  //------------------------------------
  //  Node->node mapping for the CrsMatrix graph

  std::vector< std::vector< unsigned > > node_node_ids( owned_node );
  std::vector< unsigned > node_node_begin( owned_node );

  size_t offset = 0 ;
  for ( size_t i = 0 ; i < owned_node ; ++i ) {
    const size_t j_end = node_elem_ids.row_map[i+1];
          size_t j     = node_elem_ids.row_map[i];

    node_node_begin[i] = offset ;

    std::vector< unsigned > & work = node_node_ids[i] ;

    for ( ; j < j_end ; ++j ) {
      const size_t elem_id = node_elem_ids.entries(j,0);
      for ( size_t k = 0 ; k < ElemNodeCount ; ++k ) {
        work.push_back( elem_node_ids( elem_id , k ) );
      }
    }

    std::sort( work.begin() , work.end() );

    work.erase( std::unique( work.begin() , work.end() ) , work.end() );

    offset += work.size();
  }

  return Kokkos::create_crsarray< graph_type >( "node_node_ids" , node_node_ids );
}

//----------------------------------------------------------------------------

template< class FEMeshType , class ValueType , class Device >
void fill_linear_system( const FEMeshType & mesh ,
                         const Kokkos::CrsMatrix< ValueType , Device > & matrix ,
                         const Kokkos::View< ValueType* , Kokkos::LayoutRight , Device > & rhs ,
                         const typename Kokkos::View< ValueType* , Kokkos::LayoutRight , Device >::HostMirror & solution )
{
  typedef Kokkos::CrsMatrix< ValueType , Device >                        matrix_type ;
  typedef Kokkos::View< ValueType* , Kokkos::LayoutRight , Device > vector_type ;

  typename matrix_type::graph_type ::HostMirror host_graph  = Kokkos::create_mirror( matrix.graph );
  typename matrix_type::values_type::HostMirror host_values = Kokkos::create_mirror( matrix.values );
  typename FEMeshType::node_coords_type::HostMirror host_node_coords = Kokkos::create_mirror_view( mesh.node_coords );

  typename vector_type::HostMirror host_rhs = Kokkos::create_mirror( rhs );

  Kokkos::deep_copy( host_node_coords , mesh.node_coords );

  for ( unsigned iRow = 0 ; iRow < solution.dimension_0() ; ++iRow ) {
    for ( unsigned j = 0 ; j < solution.dimension_1() ; ++j ) {
      host_rhs(iRow,j) = 0 ;
      solution(iRow,j) = host_node_coords(iRow,0) * 5 +
                         host_node_coords(iRow,1) * 2 +
                         host_node_coords(iRow,2) * 1 +
                         0.01 * ( 1 + double(j) / double(solution.dimension_1()) );
    }
  }

  for ( unsigned iRow = 0 ; iRow < host_graph.row_map.dimension_0() - 1 ; ++iRow ) {
    const unsigned endEntry = host_graph.row_map(iRow+1);

    for ( unsigned iEntry = host_graph.row_map(iRow) ; iEntry < endEntry ; ++iEntry ) {

      const unsigned iCol = host_graph.entries(iEntry);

      for ( unsigned j = 0 ; j < solution.dimension_1() ; ++j ) {

        if ( iRow == iCol ) {
          host_values(iEntry,j) = 27 + 0.01 * double(j) / double(solution.dimension_1());
        }
        else {
          host_values(iEntry,j) = -1 - 0.01 * double(j) / double(solution.dimension_1());
        }

        const double sCol = host_node_coords(iCol,0) * 5 +
                            host_node_coords(iCol,1) * 2 +
                            host_node_coords(iCol,2) * 1 +
                            0.01 * ( 1 + double(j) / double(solution.dimension_1()) );

        host_rhs(iRow,j) += host_values(iEntry,j) * sCol ;
      }
    }
  }

  Kokkos::deep_copy( rhs , host_rhs );
  Kokkos::deep_copy( matrix.values , host_values );
}

//----------------------------------------------------------------------------

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

