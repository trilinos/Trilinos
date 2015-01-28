/*
//@HEADER
// ************************************************************************
//
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_EXAMPLE_FENLFUNCTORS_HPP
#define KOKKOS_EXAMPLE_FENLFUNCTORS_HPP

#include <stdio.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <limits>

#include <Kokkos_Core.hpp>
#include <Kokkos_Pair.hpp>
#include <Kokkos_UnorderedMap.hpp>
#include <Kokkos_StaticCrsGraph.hpp>

#include <impl/Kokkos_Timer.hpp>

#include <BoxElemFixture.hpp>
#include <HexElement.hpp>

#include "Sacado_Kokkos.hpp"

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Example {
namespace FENL {

template< typename ValueType , class Space >
struct CrsMatrix {
  typedef Kokkos::StaticCrsGraph< unsigned , Space , void , unsigned >  StaticCrsGraphType ;
  typedef View< ValueType * , Space > coeff_type ;

  StaticCrsGraphType  graph ;
  coeff_type          coeff ;

  CrsMatrix() : graph(), coeff() {}

  CrsMatrix( const StaticCrsGraphType & arg_graph )
    : graph( arg_graph )
    , coeff( "crs_matrix_coeff" , arg_graph.entries.dimension_0() )
    {}
};

template< class ElemNodeIdView , class CrsGraphType , unsigned ElemNode >
class NodeNodeGraph {
public:

  typedef typename ElemNodeIdView::device_type device_type ;
  typedef pair<unsigned,unsigned> key_type ;

  typedef Kokkos::UnorderedMap< key_type, void , device_type > SetType ;
  typedef typename CrsGraphType::row_map_type::non_const_type  RowMapType ;
  typedef Kokkos::View< unsigned ,  device_type >              UnsignedValue ;

  // Static dimensions of 0 generate compiler warnings or errors.
  typedef Kokkos::View< unsigned*[ElemNode][ElemNode] , device_type >
    ElemGraphType ;

private:

  enum PhaseType { FILL_NODE_SET ,
                   SCAN_NODE_COUNT ,
                   FILL_GRAPH_ENTRIES ,
                   SORT_GRAPH_ENTRIES ,
                   FILL_ELEMENT_GRAPH };

  const unsigned        node_count ;
  const ElemNodeIdView  elem_node_id ;
  UnsignedValue         row_total ;
  RowMapType            row_count ;
  RowMapType            row_map ;
  SetType               node_node_set ;
  PhaseType             phase ;

public:

  CrsGraphType          graph ;
  ElemGraphType         elem_graph ;

  struct Times
  {
    double ratio;
    double fill_node_set;
    double scan_node_count;
    double fill_graph_entries;
    double sort_graph_entries;
    double fill_element_graph;
  };

  NodeNodeGraph( const ElemNodeIdView & arg_elem_node_id ,
                 const unsigned         arg_node_count,
                 Times & results
               )
    : node_count(arg_node_count)
    , elem_node_id( arg_elem_node_id )
    , row_total( "row_total" )
    , row_count(Kokkos::ViewAllocateWithoutInitializing("row_count") , node_count ) // will deep_copy to 0 inside loop
    , row_map( "graph_row_map" , node_count + 1 )
    , node_node_set()
    , phase( FILL_NODE_SET )
    , graph()
    , elem_graph()
   {
      //--------------------------------
      // Guess at capacity required for the map:

      Kokkos::Impl::Timer wall_clock ;

      wall_clock.reset();
      phase = FILL_NODE_SET ;

      // upper bound on the capacity
      size_t set_capacity = (((28ull * node_count) / 2ull)*4ull)/3ull;


      // Increase capacity until the (node,node) map is successfully filled.
      {
        // Zero the row count to restart the fill
        Kokkos::deep_copy( row_count , 0u );

        node_node_set = SetType( set_capacity );

        // May be larger that requested:
        set_capacity = node_node_set.capacity();

        Kokkos::parallel_for( elem_node_id.dimension_0() , *this );
      }

      device_type::fence();
      results.ratio = (double)node_node_set.size() / (double)node_node_set.capacity();
      results.fill_node_set = wall_clock.seconds();
      //--------------------------------

      wall_clock.reset();
      phase = SCAN_NODE_COUNT ;

      // Exclusive scan of row_count into row_map
      // including the final total in the 'node_count + 1' position.
      // Zero the 'row_count' values.
      Kokkos::parallel_scan( node_count , *this );

      // Zero the row count for the fill:
      Kokkos::deep_copy( row_count , 0u );

      unsigned graph_entry_count = 0 ;

      Kokkos::deep_copy( graph_entry_count , row_total );

      // Assign graph's row_map and allocate graph's entries
      graph.row_map = row_map ;
      graph.entries = typename CrsGraphType::entries_type( "graph_entries" , graph_entry_count );

      //--------------------------------
      // Fill graph's entries from the (node,node) set.

      device_type::fence();
      results.scan_node_count = wall_clock.seconds();

      wall_clock.reset();
      phase = FILL_GRAPH_ENTRIES ;
      Kokkos::parallel_for( node_node_set.capacity() , *this );

      device_type::fence();
      results.fill_graph_entries = wall_clock.seconds();

      //--------------------------------
      // Done with the temporary sets and arrays
      wall_clock.reset();
      phase = SORT_GRAPH_ENTRIES ;

      row_total = UnsignedValue();
      row_count = RowMapType();
      row_map   = RowMapType();
      node_node_set.clear();

      //--------------------------------

      Kokkos::parallel_for( node_count , *this );

      device_type::fence();
      results.sort_graph_entries = wall_clock.seconds();

      //--------------------------------
      // Element-to-graph mapping:
      wall_clock.reset();
      phase = FILL_ELEMENT_GRAPH ;
      elem_graph = ElemGraphType("elem_graph", elem_node_id.dimension_0() );
      Kokkos::parallel_for( elem_node_id.dimension_0() , *this );

      device_type::fence();
      results.fill_element_graph = wall_clock.seconds();
    }

  //------------------------------------
  // parallel_for: create map and count row length

  KOKKOS_INLINE_FUNCTION
  void fill_set( const unsigned ielem ) const
  {
    // Loop over element's (row_local_node,col_local_node) pairs:
    for ( unsigned row_local_node = 0 ; row_local_node < elem_node_id.dimension_1() ; ++row_local_node ) {

      const unsigned row_node = elem_node_id( ielem , row_local_node );

      for ( unsigned col_local_node = row_local_node ; col_local_node < elem_node_id.dimension_1() ; ++col_local_node ) {

        const unsigned col_node = elem_node_id( ielem , col_local_node );

        // If either node is locally owned then insert the pair into the unordered map:

        if ( row_node < row_count.dimension_0() || col_node < row_count.dimension_0() ) {

          const key_type key = (row_node < col_node) ? make_pair( row_node, col_node ) : make_pair( col_node, row_node ) ;

          const typename SetType::insert_result result = node_node_set.insert( key );

          if ( result.success() ) {
            if ( row_node < row_count.dimension_0() ) { atomic_fetch_add( & row_count( row_node ) , 1 ); }
            if ( col_node < row_count.dimension_0() && col_node != row_node ) { atomic_fetch_add( & row_count( col_node ) , 1 ); }
          }
        }
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void fill_graph_entries( const unsigned iset ) const
  {
    if ( node_node_set.valid_at(iset) ) {
      const key_type key = node_node_set.key_at(iset) ;
      const unsigned row_node = key.first ;
      const unsigned col_node = key.second ;

      if ( row_node < row_count.dimension_0() ) {
        const unsigned offset = graph.row_map( row_node ) + atomic_fetch_add( & row_count( row_node ) , 1 );
        graph.entries( offset ) = col_node ;
      }

      if ( col_node < row_count.dimension_0() && col_node != row_node ) {
        const unsigned offset = graph.row_map( col_node ) + atomic_fetch_add( & row_count( col_node ) , 1 );
        graph.entries( offset ) = row_node ;
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void sort_graph_entries( const unsigned irow ) const
  {
    typedef typename CrsGraphType::size_type size_type;
    const size_type row_beg = graph.row_map( irow );
    const size_type row_end = graph.row_map( irow + 1 );
    for ( size_type i = row_beg + 1 ; i < row_end ; ++i ) {
      const typename CrsGraphType::data_type col = graph.entries(i);
      size_type j = i ;
      for ( ; row_beg < j && col < graph.entries(j-1) ; --j ) {
        graph.entries(j) = graph.entries(j-1);
      }
      graph.entries(j) = col ;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void fill_elem_graph_map( const unsigned ielem ) const
  {
    typedef typename CrsGraphType::data_type entry_type;
    for ( unsigned row_local_node = 0 ; row_local_node < elem_node_id.dimension_1() ; ++row_local_node ) {

      const unsigned row_node = elem_node_id( ielem , row_local_node );

      for ( unsigned col_local_node = 0 ; col_local_node < elem_node_id.dimension_1() ; ++col_local_node ) {

        const unsigned col_node = elem_node_id( ielem , col_local_node );

        entry_type entry = 0 ;

        if ( row_node + 1 < graph.row_map.dimension_0() ) {

          const entry_type entry_end = static_cast<entry_type> (graph.row_map( row_node + 1 ));

          entry = graph.row_map( row_node );

          for ( ; entry < entry_end && graph.entries(entry) != static_cast<entry_type> (col_node) ; ++entry );

          if ( entry == entry_end ) entry = ~0u ;
        }

        elem_graph( ielem , row_local_node , col_local_node ) = entry ;
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( const unsigned iwork ) const
  {
    if ( phase == FILL_NODE_SET ) {
      fill_set( iwork );
    }
    else if ( phase == FILL_GRAPH_ENTRIES ) {
      fill_graph_entries( iwork );
    }
    else if ( phase == SORT_GRAPH_ENTRIES ) {
      sort_graph_entries( iwork );
    }
    else if ( phase == FILL_ELEMENT_GRAPH ) {
      fill_elem_graph_map( iwork );
    }
  }

  //------------------------------------
  // parallel_scan: row offsets

  typedef unsigned value_type ;

  KOKKOS_INLINE_FUNCTION
  void operator()( const unsigned irow , unsigned & update , const bool final ) const
  {
    // exclusive scan
    if ( final ) { row_map( irow ) = update ; }

    update += row_count( irow );

    if ( final ) {
      if ( irow + 1 == row_count.dimension_0() ) {
        row_map( irow + 1 ) = update ;
        row_total()         = update ;
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void init( unsigned & update ) const { update = 0 ; }

  KOKKOS_INLINE_FUNCTION
  void join( volatile unsigned & update , const volatile unsigned & input ) const { update += input ; }

  //------------------------------------
};

} /* namespace FENL */
} /* namespace Example */
} /* namespace Kokkos  */

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Example {
namespace FENL {

struct ElementComputationConstantCoefficient {
  enum { is_constant = true };

  const double coeff_k ;

  KOKKOS_INLINE_FUNCTION
  double operator()( double /* x */
                  , double /* y */
                  , double /* z */
                  ) const
    { return coeff_k ; }

  ElementComputationConstantCoefficient( const double val )
    : coeff_k( val ) {}

  ElementComputationConstantCoefficient( const ElementComputationConstantCoefficient & rhs )
    : coeff_k( rhs.coeff_k ) {}
};

// Traits class to get the value type in various kinds of arrays
template <typename array_type>
struct ArrayTraits;

template <typename T, typename L, typename D, typename M, typename S>
struct ArrayTraits< Kokkos::View<T,L,D,M,S> > {
  typedef Kokkos::View<T,L,D,M,S> array_type;
  typedef typename array_type::non_const_value_type value_type;
};

template <typename T>
struct ArrayTraits< T* > {
  typedef T* array_type;
  typedef T value_type;
};

template <typename T>
struct ArrayTraits< const T* > {
  typedef const T* array_type;
  typedef T value_type;
};

template< class DeviceType , BoxElemPart::ElemOrder Order ,
          class CoordinateMap , typename ScalarType >
class ElementComputationBase
{
public:

  typedef Kokkos::Example::BoxElemFixture< DeviceType, Order, CoordinateMap >  mesh_type ;
  typedef Kokkos::Example::HexElement_Data< mesh_type::ElemNode >              element_data_type ;

  //------------------------------------

  typedef DeviceType   device_type ;
  typedef ScalarType   scalar_type ;

  typedef CrsMatrix< ScalarType , DeviceType >  sparse_matrix_type ;
  typedef typename sparse_matrix_type::StaticCrsGraphType                                       sparse_graph_type ;
  typedef Kokkos::View< scalar_type* , Kokkos::LayoutLeft, device_type > vector_type ;

  //------------------------------------

  static const unsigned SpatialDim       = element_data_type::spatial_dimension ;
  static const unsigned TensorDim        = SpatialDim * SpatialDim ;
  static const unsigned ElemNodeCount    = element_data_type::element_node_count ;
  static const unsigned FunctionCount    = element_data_type::function_count ;
  static const unsigned IntegrationCount = element_data_type::integration_count ;

  //------------------------------------

  typedef typename mesh_type::node_coord_type                                      node_coord_type ;
  typedef typename mesh_type::elem_node_type                                       elem_node_type ;
  typedef Kokkos::View< scalar_type*[FunctionCount][FunctionCount] , device_type > elem_matrices_type ;
  typedef Kokkos::View< scalar_type*[FunctionCount] ,                device_type > elem_vectors_type ;

  typedef typename NodeNodeGraph< elem_node_type , sparse_graph_type , ElemNodeCount >::ElemGraphType elem_graph_type ;
  typedef ElementComputationConstantCoefficient coeff_function_type;

  //------------------------------------


  //------------------------------------
  // Computational data:

  const element_data_type   elem_data ;
  const elem_node_type      elem_node_ids ;
  const node_coord_type     node_coords ;
  const elem_graph_type     elem_graph ;
  const elem_matrices_type  elem_jacobians ;
  const elem_vectors_type   elem_residuals ;
  const vector_type         solution ;
  const vector_type         residual ;
  const sparse_matrix_type  jacobian ;
  const coeff_function_type coeff_function ;

  ElementComputationBase( const ElementComputationBase & rhs )
    : elem_data()
    , elem_node_ids( rhs.elem_node_ids )
    , node_coords(   rhs.node_coords )
    , elem_graph(    rhs.elem_graph )
    , elem_jacobians( rhs.elem_jacobians )
    , elem_residuals( rhs.elem_residuals )
    , solution( rhs.solution )
    , residual( rhs.residual )
    , jacobian( rhs.jacobian )
    , coeff_function( rhs.coeff_function )
    {}

  ElementComputationBase( const mesh_type          & arg_mesh ,
                          const coeff_function_type& arg_coeff_function ,
                          const vector_type        & arg_solution ,
                          const elem_graph_type    & arg_elem_graph ,
                          const sparse_matrix_type & arg_jacobian ,
                          const vector_type        & arg_residual )
    : elem_data()
    , elem_node_ids( arg_mesh.elem_node() )
    , node_coords(   arg_mesh.node_coord() )
    , elem_graph(    arg_elem_graph )
    , elem_jacobians()
    , elem_residuals()
    , solution( arg_solution )
    , residual( arg_residual )
    , jacobian( arg_jacobian )
    , coeff_function( arg_coeff_function )
    {}

  //------------------------------------

  KOKKOS_INLINE_FUNCTION
  double transform_gradients(
    const double grad[][ FunctionCount ] , // Gradient of bases master element
    const double x[] ,
    const double y[] ,
    const double z[] ,
    double dpsidx[] ,
    double dpsidy[] ,
    double dpsidz[] ) const
  {
    enum { j11 = 0 , j12 = 1 , j13 = 2 ,
           j21 = 3 , j22 = 4 , j23 = 5 ,
           j31 = 6 , j32 = 7 , j33 = 8 };

    // Jacobian accumulation:

    double J[ TensorDim ] = { 0, 0, 0,  0, 0, 0,  0, 0, 0 };

    for( unsigned i = 0; i < FunctionCount ; ++i ) {
      const double x1 = x[i] ;
      const double x2 = y[i] ;
      const double x3 = z[i] ;

      const double g1 = grad[0][i] ;
      const double g2 = grad[1][i] ;
      const double g3 = grad[2][i] ;

      J[j11] += g1 * x1 ;
      J[j12] += g1 * x2 ;
      J[j13] += g1 * x3 ;

      J[j21] += g2 * x1 ;
      J[j22] += g2 * x2 ;
      J[j23] += g2 * x3 ;

      J[j31] += g3 * x1 ;
      J[j32] += g3 * x2 ;
      J[j33] += g3 * x3 ;
    }

    // Inverse jacobian:

    double invJ[ TensorDim ] = {
      static_cast<double>( J[j22] * J[j33] - J[j23] * J[j32] ) ,
      static_cast<double>( J[j13] * J[j32] - J[j12] * J[j33] ) ,
      static_cast<double>( J[j12] * J[j23] - J[j13] * J[j22] ) ,

      static_cast<double>( J[j23] * J[j31] - J[j21] * J[j33] ) ,
      static_cast<double>( J[j11] * J[j33] - J[j13] * J[j31] ) ,
      static_cast<double>( J[j13] * J[j21] - J[j11] * J[j23] ) ,

      static_cast<double>( J[j21] * J[j32] - J[j22] * J[j31] ) ,
      static_cast<double>( J[j12] * J[j31] - J[j11] * J[j32] ) ,
      static_cast<double>( J[j11] * J[j22] - J[j12] * J[j21] ) };

    const double detJ = J[j11] * invJ[j11] +
                       J[j21] * invJ[j12] +
                       J[j31] * invJ[j13] ;

    const double detJinv = 1.0 / detJ ;

    for ( unsigned i = 0 ; i < TensorDim ; ++i ) { invJ[i] *= detJinv ; }

    // Transform gradients:

    for( unsigned i = 0; i < FunctionCount ; ++i ) {
      const double g0 = grad[0][i];
      const double g1 = grad[1][i];
      const double g2 = grad[2][i];

      dpsidx[i] = g0 * invJ[j11] + g1 * invJ[j12] + g2 * invJ[j13];
      dpsidy[i] = g0 * invJ[j21] + g1 * invJ[j22] + g2 * invJ[j23];
      dpsidz[i] = g0 * invJ[j31] + g1 * invJ[j32] + g2 * invJ[j33];
    }

    return detJ ;
  }

  //------------------------------------

  template <typename values_type, typename res_type, typename mat_type>
  KOKKOS_INLINE_FUNCTION
  void gatherSolution(const unsigned ielem,
                      values_type val,
                      unsigned node_index[],
                      double x[], double y[], double z[],
                      res_type res,
                      mat_type mat) const
  {
    for ( unsigned i = 0 ; i < ElemNodeCount ; ++i ) {
      const unsigned ni = elem_node_ids( ielem , i );

      node_index[i] = ni ;

      x[i] = this->node_coords( ni , 0 );
      y[i] = this->node_coords( ni , 1 );
      z[i] = this->node_coords( ni , 2 );

      val[i] = solution( ni ) ;
      res[i] = 0 ;

      for( unsigned j = 0; j < FunctionCount ; j++){
        mat(i,j) = 0 ;
      }
    }
  }

  template <typename values_type, typename res_type, typename mat_type>
  KOKKOS_INLINE_FUNCTION
  void gatherSolutionGlobal(const unsigned ielem,
                            const values_type& val,
                            unsigned node_index[],
                            double x[], double y[], double z[],
                            const res_type& res,
                            const mat_type& mat) const
  {
    for ( unsigned i = 0 ; i < ElemNodeCount ; ++i ) {
      const unsigned ni = elem_node_ids( ielem , i );

      node_index[i] = ni ;

      x[i] = this->node_coords( ni , 0 );
      y[i] = this->node_coords( ni , 1 );
      z[i] = this->node_coords( ni , 2 );

      val(ielem,i) = solution( ni ) ;
      res(ielem,i) = 0 ;

      for( unsigned j = 0; j < FunctionCount ; j++){
        mat(ielem,i,j) = 0 ;
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void gatherSolution(const unsigned ielem,
                      scalar_type val[],
                      unsigned node_index[],
                      double x[], double y[], double z[],
                      scalar_type res[],
                      scalar_type mat[][FunctionCount]) const
  {
    for ( unsigned i = 0 ; i < ElemNodeCount ; ++i ) {
      const unsigned ni = elem_node_ids( ielem , i );

      node_index[i] = ni ;

      x[i] = this->node_coords( ni , 0 );
      y[i] = this->node_coords( ni , 1 );
      z[i] = this->node_coords( ni , 2 );

      val[i] = solution( ni ) ;
      res[i] = 0 ;

      for( unsigned j = 0; j < FunctionCount ; j++){
        mat[i][j] = 0 ;
      }
    }
  }

  template <typename values_type, typename res_type>
  KOKKOS_INLINE_FUNCTION
  void gatherSolutionFad(const unsigned ielem,
                         values_type val,
                         unsigned node_index[],
                         double x[], double y[], double z[],
                         res_type res) const
  {
    typedef typename ArrayTraits<values_type>::value_type fad_scalar_type;
    for ( unsigned i = 0 ; i < ElemNodeCount ; ++i ) {
      const unsigned ni = elem_node_ids( ielem , i );

      node_index[i] = ni ;

      x[i] = this->node_coords( ni , 0 );
      y[i] = this->node_coords( ni , 1 );
      z[i] = this->node_coords( ni , 2 );

      val[i] = fad_scalar_type( FunctionCount, i, solution( ni ) );
      res[i] = 0 ;
    }
  }

  template <typename values_type, typename res_type>
  KOKKOS_INLINE_FUNCTION
  void gatherSolutionFadGlobal(const unsigned ielem,
                               const values_type& val,
                               unsigned node_index[],
                               double x[], double y[], double z[],
                               const res_type& res) const
  {
    typedef typename values_type::value_type fad_scalar_type;
    for ( unsigned i = 0 ; i < ElemNodeCount ; ++i ) {
      const unsigned ni = elem_node_ids( ielem , i );

      node_index[i] = ni ;

      x[i] = this->node_coords( ni , 0 );
      y[i] = this->node_coords( ni , 1 );
      z[i] = this->node_coords( ni , 2 );

      //val(ielem,i) = fad_scalar_type( FunctionCount, i, solution( ni ) );
      val(ielem,i).val() = solution( ni );
      val(ielem,i).fastAccessDx(i) = 1.0;
      res(ielem,i) = 0 ;
    }
  }

  //------------------------------------

  template <typename res_type, typename mat_type>
  KOKKOS_INLINE_FUNCTION
  void scatterResidualView(const unsigned ielem,
                           const unsigned node_index[],
                           res_type res,
                           mat_type mat) const
  {
    for( unsigned i = 0 ; i < FunctionCount ; i++ ) {
      const unsigned row = node_index[i] ;
      if ( row < residual.dimension_0() ) {
        atomic_add( & residual( row ) , res[i] );

        for( unsigned j = 0 ; j < FunctionCount ; j++ ) {
          const unsigned entry = elem_graph( ielem , i , j );
          if ( entry != ~0u ) {
            atomic_add( & jacobian.coeff( entry ) , mat(i,j) );
          }
        }
      }
    }
  }

   template <typename res_type, typename mat_type>
  KOKKOS_INLINE_FUNCTION
  void scatterResidualViewGlobal(const unsigned ielem,
                                 const unsigned node_index[],
                                 const res_type& res,
                                 const mat_type& mat) const
  {
    for( unsigned i = 0 ; i < FunctionCount ; i++ ) {
      const unsigned row = node_index[i] ;
      if ( row < residual.dimension_0() ) {
        atomic_add( & residual( row ) , res(ielem,i) );

        for( unsigned j = 0 ; j < FunctionCount ; j++ ) {
          const unsigned entry = elem_graph( ielem , i , j );
          if ( entry != ~0u ) {
            atomic_add( & jacobian.coeff( entry ) , mat(ielem,i,j) );
          }
        }
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void scatterResidual(const unsigned ielem,
                       const unsigned node_index[],
                       const scalar_type res[],
                       const scalar_type mat[][FunctionCount]) const
  {
    for( unsigned i = 0 ; i < FunctionCount ; i++ ) {
      const unsigned row = node_index[i] ;
      if ( row < residual.dimension_0() ) {
        atomic_add( & residual( row ) , res[i] );

        for( unsigned j = 0 ; j < FunctionCount ; j++ ) {
          const unsigned entry = elem_graph( ielem , i , j );
          if ( entry != ~0u ) {
            atomic_add( & jacobian.coeff( entry ) , mat[i][j] );
          }
        }
      }
    }
  }

  template <typename res_type>
  KOKKOS_INLINE_FUNCTION
  void scatterResidualFad(const unsigned ielem,
                          const unsigned node_index[],
                          res_type res) const
  {
    for( unsigned i = 0 ; i < FunctionCount ; i++ ) {
      const unsigned row = node_index[i] ;
      if ( row < residual.dimension_0() ) {
        atomic_add( & residual( row ) , res[i].val() );

        for( unsigned j = 0 ; j < FunctionCount ; j++ ) {
          const unsigned entry = elem_graph( ielem , i , j );
          if ( entry != ~0u ) {
            atomic_add( & jacobian.coeff( entry ) , res[i].fastAccessDx(j) );
          }
        }
      }
    }
  }

  template <typename res_type>
  KOKKOS_INLINE_FUNCTION
  void scatterResidualFadGlobal(const unsigned ielem,
                                const unsigned node_index[],
                                const res_type& res) const
  {
    for( unsigned i = 0 ; i < FunctionCount ; i++ ) {
      const unsigned row = node_index[i] ;
      if ( row < residual.dimension_0() ) {
        atomic_add( & residual( row ) , res(ielem,i).val() );

        for( unsigned j = 0 ; j < FunctionCount ; j++ ) {
          const unsigned entry = elem_graph( ielem , i , j );
          if ( entry != ~0u ) {
            atomic_add( & jacobian.coeff( entry ) , res(ielem,i).fastAccessDx(j) );
          }
        }
      }
    }
  }

  //------------------------------------

  template <typename values_type, typename res_type>
  KOKKOS_INLINE_FUNCTION
  void computeElementResidual(const values_type dof_values ,
                              const double x[],
                              const double y[],
                              const double z[],
                              res_type  elem_res ) const
  {
    typedef typename ArrayTraits<values_type>::value_type array_scalar_type;

    scalar_type coeff_k = 3.456;
    double dpsidx[ FunctionCount ] ;
    double dpsidy[ FunctionCount ] ;
    double dpsidz[ FunctionCount ] ;
    for ( unsigned i = 0 ; i < IntegrationCount ; ++i ) {

      const double integ_weight = elem_data.weights[i];
      const double* bases_vals = elem_data.values[i];

      const double detJ =
        transform_gradients( elem_data.gradients[i] ,
                             x , y , z ,
                             dpsidx , dpsidy , dpsidz );

      array_scalar_type value_at_pt = 0 ;
      array_scalar_type gradx_at_pt = 0 ;
      array_scalar_type grady_at_pt = 0 ;
      array_scalar_type gradz_at_pt = 0 ;

      for ( unsigned m = 0 ; m < FunctionCount ; m++ ) {
        value_at_pt += dof_values[m] * bases_vals[m] ;
        gradx_at_pt += dof_values[m] * dpsidx[m] ;
        grady_at_pt += dof_values[m] * dpsidy[m] ;
        gradz_at_pt += dof_values[m] * dpsidz[m] ;
      }

      const scalar_type k_detJ_weight =
        coeff_k * detJ * integ_weight ;
      const array_scalar_type res_val =
        std::exp( value_at_pt * value_at_pt ) * detJ * integ_weight ;

    // $$ R_i = \int_{\Omega} \nabla \phi_i \cdot (k \nabla T) + \phi_i T^2 d \Omega $$

      for ( unsigned m = 0; m < FunctionCount; ++m) {
        const double bases_val_m = bases_vals[m];
        const double dpsidx_m    = dpsidx[m] ;
        const double dpsidy_m    = dpsidy[m] ;
        const double dpsidz_m    = dpsidz[m] ;

        elem_res[m] += k_detJ_weight * ( dpsidx_m * gradx_at_pt +
                                         dpsidy_m * grady_at_pt +
                                         dpsidz_m * gradz_at_pt ) +
                       res_val * bases_val_m ;
      }
    }
  }

  template <typename values_type, typename res_type>
  KOKKOS_INLINE_FUNCTION
  void computeElementResidualGlobal(const unsigned ielem,
                                    const values_type& dof_values ,
                                    const double x[],
                                    const double y[],
                                    const double z[],
                                    const res_type&  elem_res ) const
  {
    typedef typename values_type::value_type array_scalar_type;

    scalar_type coeff_k = 3.456;
    double dpsidx[ FunctionCount ] ;
    double dpsidy[ FunctionCount ] ;
    double dpsidz[ FunctionCount ] ;
    for ( unsigned i = 0 ; i < IntegrationCount ; ++i ) {

      const double integ_weight = elem_data.weights[i];
      const double* bases_vals = elem_data.values[i];

      const double detJ =
        transform_gradients( elem_data.gradients[i] ,
                             x , y , z ,
                             dpsidx , dpsidy , dpsidz );

      array_scalar_type value_at_pt = 0 ;
      array_scalar_type gradx_at_pt = 0 ;
      array_scalar_type grady_at_pt = 0 ;
      array_scalar_type gradz_at_pt = 0 ;

      for ( unsigned m = 0 ; m < FunctionCount ; m++ ) {
        value_at_pt += dof_values(ielem,m) * bases_vals[m] ;
        gradx_at_pt += dof_values(ielem,m) * dpsidx[m] ;
        grady_at_pt += dof_values(ielem,m) * dpsidy[m] ;
        gradz_at_pt += dof_values(ielem,m) * dpsidz[m] ;
      }

      const scalar_type k_detJ_weight =
        coeff_k * detJ * integ_weight ;
      const array_scalar_type res_val =
        std::exp( value_at_pt * value_at_pt ) * detJ * integ_weight ;

    // $$ R_i = \int_{\Omega} \nabla \phi_i \cdot (k \nabla T) + \phi_i T^2 d \Omega $$

      for ( unsigned m = 0; m < FunctionCount; ++m) {
        const double bases_val_m = bases_vals[m];
        const double dpsidx_m    = dpsidx[m] ;
        const double dpsidy_m    = dpsidy[m] ;
        const double dpsidz_m    = dpsidz[m] ;

        elem_res(ielem,m) += k_detJ_weight * ( dpsidx_m * gradx_at_pt +
                                               dpsidy_m * grady_at_pt +
                                               dpsidz_m * gradz_at_pt ) +
                             res_val * bases_val_m ;
      }
    }
  }

  template <typename values_type, typename res_type, typename mat_type>
  KOKKOS_INLINE_FUNCTION
  void computeElementResidualJacobianView(const values_type dof_values ,
                                          const double x[],
                                          const double y[],
                                          const double z[],
                                          res_type  elem_res ,
                                          mat_type  elem_mat ) const
  {
    scalar_type coeff_k = 3.456;
    double dpsidx[ FunctionCount ] ;
    double dpsidy[ FunctionCount ] ;
    double dpsidz[ FunctionCount ] ;
    for ( unsigned i = 0 ; i < IntegrationCount ; ++i ) {

      const double integ_weight = elem_data.weights[i];
      const double* bases_vals = elem_data.values[i];

      const double detJ =
        this->transform_gradients( elem_data.gradients[i] ,
                                   x , y , z ,
                                   dpsidx , dpsidy , dpsidz );

      scalar_type value_at_pt = 0 ;
      scalar_type gradx_at_pt = 0 ;
      scalar_type grady_at_pt = 0 ;
      scalar_type gradz_at_pt = 0 ;

      for ( unsigned m = 0 ; m < FunctionCount ; m++ ) {
        value_at_pt += dof_values[m] * bases_vals[m] ;
        gradx_at_pt += dof_values[m] * dpsidx[m] ;
        grady_at_pt += dof_values[m] * dpsidy[m] ;
        gradz_at_pt += dof_values[m] * dpsidz[m] ;
      }

      const scalar_type k_detJ_weight =
        coeff_k * detJ * integ_weight ;
      const scalar_type res_val =
        std::exp( value_at_pt * value_at_pt ) * detJ * integ_weight ;
      const scalar_type mat_val =
        2.0 * value_at_pt * res_val;

      // $$ R_i = \int_{\Omega} \nabla \phi_i \cdot (k \nabla T) + \phi_i T^2 d \Omega $$
      // $$ J_{i,j} = \frac{\partial R_i}{\partial T_j} = \int_{\Omega} k \nabla \phi_i \cdot \nabla \phi_j + 2 \phi_i \phi_j T d \Omega $$

      for ( unsigned m = 0; m < FunctionCount; ++m) {
        const double bases_val_m = bases_vals[m];
        const double dpsidx_m    = dpsidx[m] ;
        const double dpsidy_m    = dpsidy[m] ;
        const double dpsidz_m    = dpsidz[m] ;

        elem_res[m] += k_detJ_weight * ( dpsidx_m * gradx_at_pt +
                                         dpsidy_m * grady_at_pt +
                                         dpsidz_m * gradz_at_pt ) +
                       res_val * bases_val_m ;

        for( unsigned n = 0; n < FunctionCount; n++) {

          elem_mat(m,n) += k_detJ_weight * ( dpsidx_m * dpsidx[n] +
                                             dpsidy_m * dpsidy[n] +
                                             dpsidz_m * dpsidz[n] ) +
                           mat_val * bases_val_m * bases_vals[n];
        }
      }
    }
  }

  template <typename values_type, typename res_type, typename mat_type>
  KOKKOS_INLINE_FUNCTION
  void computeElementResidualJacobianViewGlobal(const unsigned ielem,
                                                const values_type& dof_values ,
                                                const double x[],
                                                const double y[],
                                                const double z[],
                                                const res_type& elem_res ,
                                                const mat_type& elem_mat ) const
  {
    scalar_type coeff_k = 3.456;
    double dpsidx[ FunctionCount ] ;
    double dpsidy[ FunctionCount ] ;
    double dpsidz[ FunctionCount ] ;
    for ( unsigned i = 0 ; i < IntegrationCount ; ++i ) {

      const double integ_weight = elem_data.weights[i];
      const double* bases_vals = elem_data.values[i];

      const double detJ =
        this->transform_gradients( elem_data.gradients[i] ,
                                   x , y , z ,
                                   dpsidx , dpsidy , dpsidz );

      scalar_type value_at_pt = 0 ;
      scalar_type gradx_at_pt = 0 ;
      scalar_type grady_at_pt = 0 ;
      scalar_type gradz_at_pt = 0 ;

      for ( unsigned m = 0 ; m < FunctionCount ; m++ ) {
        value_at_pt += dof_values(ielem,m) * bases_vals[m] ;
        gradx_at_pt += dof_values(ielem,m) * dpsidx[m] ;
        grady_at_pt += dof_values(ielem,m) * dpsidy[m] ;
        gradz_at_pt += dof_values(ielem,m) * dpsidz[m] ;
      }

      const scalar_type k_detJ_weight =
        coeff_k * detJ * integ_weight ;
      const scalar_type res_val =
        std::exp( value_at_pt * value_at_pt ) * detJ * integ_weight ;
      const scalar_type mat_val =
        2.0 * value_at_pt * res_val;

      // $$ R_i = \int_{\Omega} \nabla \phi_i \cdot (k \nabla T) + \phi_i T^2 d \Omega $$
      // $$ J_{i,j} = \frac{\partial R_i}{\partial T_j} = \int_{\Omega} k \nabla \phi_i \cdot \nabla \phi_j + 2 \phi_i \phi_j T d \Omega $$

      for ( unsigned m = 0; m < FunctionCount; ++m) {
        const double bases_val_m = bases_vals[m];
        const double dpsidx_m    = dpsidx[m] ;
        const double dpsidy_m    = dpsidy[m] ;
        const double dpsidz_m    = dpsidz[m] ;

        elem_res(ielem,m) += k_detJ_weight * ( dpsidx_m * gradx_at_pt +
                                               dpsidy_m * grady_at_pt +
                                               dpsidz_m * gradz_at_pt ) +
                             res_val * bases_val_m ;

        for( unsigned n = 0; n < FunctionCount; n++) {

          elem_mat(ielem,m,n) += k_detJ_weight * ( dpsidx_m * dpsidx[n] +
                                                   dpsidy_m * dpsidy[n] +
                                                   dpsidz_m * dpsidz[n] ) +
                                 mat_val * bases_val_m * bases_vals[n];
        }
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void computeElementResidualJacobian(const scalar_type dof_values[] ,
                                      const double x[],
                                      const double y[],
                                      const double z[],
                                      scalar_type  elem_res[] ,
                                      scalar_type  elem_mat[][FunctionCount] ) const
  {
    scalar_type coeff_k = 3.456;
    double dpsidx[ FunctionCount ] ;
    double dpsidy[ FunctionCount ] ;
    double dpsidz[ FunctionCount ] ;
    for ( unsigned i = 0 ; i < IntegrationCount ; ++i ) {

      const double integ_weight = elem_data.weights[i];
      const double* bases_vals = elem_data.values[i];

      const double detJ =
        this->transform_gradients( elem_data.gradients[i] ,
                                   x , y , z ,
                                   dpsidx , dpsidy , dpsidz );

      scalar_type value_at_pt = 0 ;
      scalar_type gradx_at_pt = 0 ;
      scalar_type grady_at_pt = 0 ;
      scalar_type gradz_at_pt = 0 ;

      for ( unsigned m = 0 ; m < FunctionCount ; m++ ) {
        value_at_pt += dof_values[m] * bases_vals[m] ;
        gradx_at_pt += dof_values[m] * dpsidx[m] ;
        grady_at_pt += dof_values[m] * dpsidy[m] ;
        gradz_at_pt += dof_values[m] * dpsidz[m] ;
      }

      const scalar_type k_detJ_weight =
        coeff_k * detJ * integ_weight ;
      const scalar_type res_val =
        std::exp( value_at_pt * value_at_pt ) * detJ * integ_weight ;
      const scalar_type mat_val =
        2.0 * value_at_pt * res_val;

      // $$ R_i = \int_{\Omega} \nabla \phi_i \cdot (k \nabla T) + \phi_i T^2 d \Omega $$
      // $$ J_{i,j} = \frac{\partial R_i}{\partial T_j} = \int_{\Omega} k \nabla \phi_i \cdot \nabla \phi_j + 2 \phi_i \phi_j T d \Omega $$

      for ( unsigned m = 0; m < FunctionCount; ++m) {
        scalar_type * const mat = elem_mat[m] ;
        const double bases_val_m = bases_vals[m];
        const double dpsidx_m    = dpsidx[m] ;
        const double dpsidy_m    = dpsidy[m] ;
        const double dpsidz_m    = dpsidz[m] ;

        elem_res[m] += k_detJ_weight * ( dpsidx_m * gradx_at_pt +
                                         dpsidy_m * grady_at_pt +
                                         dpsidz_m * gradz_at_pt ) +
                       res_val * bases_val_m ;

        for( unsigned n = 0; n < FunctionCount; n++) {

          mat[n] += k_detJ_weight * ( dpsidx_m * dpsidx[n] +
                                      dpsidy_m * dpsidy[n] +
                                      dpsidz_m * dpsidz[n] ) +
                    mat_val * bases_val_m * bases_vals[n];
        }
      }
    }
  }

};

enum AssemblyMethod {
  AnalyticLocal,
  AnalyticLocalView,
  AnalyticGlobalView,
  FadLocal,
  FadLocalView,
  FadGlobalView
};

template< class FiniteElementMeshType ,
          class SparseMatrixType ,
          AssemblyMethod Method
        >
class ElementComputation ;

template< class DeviceType , BoxElemPart::ElemOrder Order ,
          class CoordinateMap , typename ScalarType >
class ElementComputation
  < Kokkos::Example::BoxElemFixture< DeviceType , Order , CoordinateMap > ,
    CrsMatrix< ScalarType , DeviceType > ,
    AnalyticLocal > :
    public ElementComputationBase<DeviceType, Order, CoordinateMap,
                                  ScalarType> {
public:

  typedef ElementComputationBase<DeviceType, Order, CoordinateMap,
                                 ScalarType> base_type;

  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::device_type device_type;

  ElementComputation(const ElementComputation& rhs) : base_type(rhs) {}

  ElementComputation(
    const typename base_type::mesh_type          & arg_mesh ,
    const typename base_type::coeff_function_type& arg_coeff_function ,
    const typename base_type::vector_type        & arg_solution ,
    const typename base_type::elem_graph_type    & arg_elem_graph ,
    const typename base_type::sparse_matrix_type & arg_jacobian ,
    const typename base_type::vector_type        & arg_residual ) :
    base_type(arg_mesh, arg_coeff_function, arg_solution, arg_elem_graph,
              arg_jacobian, arg_residual) {}

  //------------------------------------

  void apply() const
  {
    const size_t nelem = this->elem_node_ids.dimension_0();
    parallel_for( nelem , *this );
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( const unsigned ielem ) const
  {
    double x[ base_type::FunctionCount ] ;
    double y[ base_type::FunctionCount ] ;
    double z[ base_type::FunctionCount ] ;
    unsigned node_index[ base_type::ElemNodeCount ];

    scalar_type val[ base_type::FunctionCount ] ;
    scalar_type elem_res[ base_type::FunctionCount ] ;
    scalar_type elem_mat[ base_type::FunctionCount ][ base_type::FunctionCount ] ;

    // Gather nodal coordinates and solution vector:
    this->gatherSolution(ielem, val, node_index, x, y, z, elem_res, elem_mat);

    // Compute nodal element residual vector and Jacobian matrix
    this->computeElementResidualJacobian( val, x, y, z, elem_res , elem_mat );

    // Scatter nodal element residual and Jacobian in global vector and matrix:
    this->scatterResidual( ielem, node_index, elem_res, elem_mat );
  }
}; /* ElementComputation */

template< class DeviceType , BoxElemPart::ElemOrder Order ,
          class CoordinateMap , typename ScalarType >
class ElementComputation
  < Kokkos::Example::BoxElemFixture< DeviceType , Order , CoordinateMap > ,
    CrsMatrix< ScalarType , DeviceType > ,
    AnalyticLocalView > :
    public ElementComputationBase<DeviceType, Order, CoordinateMap,
                                  ScalarType> {
public:

  typedef ElementComputationBase<DeviceType, Order, CoordinateMap,
                                 ScalarType> base_type;

  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::device_type device_type;

  ElementComputation(const ElementComputation& rhs) : base_type(rhs) {}

  ElementComputation(
    const typename base_type::mesh_type          & arg_mesh ,
    const typename base_type::coeff_function_type& arg_coeff_function ,
    const typename base_type::vector_type        & arg_solution ,
    const typename base_type::elem_graph_type    & arg_elem_graph ,
    const typename base_type::sparse_matrix_type & arg_jacobian ,
    const typename base_type::vector_type        & arg_residual ) :
    base_type(arg_mesh, arg_coeff_function, arg_solution, arg_elem_graph,
              arg_jacobian, arg_residual) {}

  //------------------------------------

  void apply() const
  {
    const size_t nelem = this->elem_node_ids.dimension_0();
    parallel_for( nelem , *this );
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( const unsigned ielem ) const
  {
    typedef Kokkos::View<scalar_type[base_type::FunctionCount],device_type,Kokkos::MemoryUnmanaged> local_element_array_type;
    typedef Kokkos::View<scalar_type[base_type::FunctionCount][base_type::FunctionCount],Kokkos::LayoutRight,device_type,Kokkos::MemoryUnmanaged> local_element_matrix_type;

    double x[ base_type::FunctionCount ] ;
    double y[ base_type::FunctionCount ] ;
    double z[ base_type::FunctionCount ] ;
    unsigned node_index[ base_type::ElemNodeCount ];

    scalar_type val_data[ base_type::FunctionCount ] ;
    scalar_type elem_res_data[base_type::FunctionCount] ;
    scalar_type elem_mat_data[base_type::FunctionCount*base_type::FunctionCount];

    local_element_array_type val(val_data, base_type::FunctionCount);
    local_element_array_type elem_res(elem_res_data, base_type::FunctionCount);
    local_element_matrix_type elem_mat(elem_mat_data, base_type::FunctionCount, base_type::FunctionCount);

    // Gather nodal coordinates and solution vector:
    this->gatherSolution(ielem, val, node_index, x, y, z, elem_res, elem_mat);

    // Compute nodal element residual vector and Jacobian matrix
    this->computeElementResidualJacobianView( val, x, y, z, elem_res , elem_mat );

    // Scatter nodal element residual and Jacobian in global vector and matrix:
    this->scatterResidualView( ielem, node_index, elem_res, elem_mat );
  }
}; /* ElementComputation */

template< class DeviceType , BoxElemPart::ElemOrder Order ,
          class CoordinateMap, typename ScalarType >
class ElementComputation
  < Kokkos::Example::BoxElemFixture< DeviceType , Order , CoordinateMap > ,
    CrsMatrix< ScalarType , DeviceType > ,
    AnalyticGlobalView > :
    public ElementComputationBase<DeviceType, Order, CoordinateMap,
                                  ScalarType> {
public:

  typedef ElementComputationBase<DeviceType, Order, CoordinateMap,
                                 ScalarType> base_type;

  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::device_type device_type;

  typedef Kokkos::View<scalar_type*[base_type::FunctionCount],device_type> local_element_array_type;
    typedef Kokkos::View<scalar_type*[base_type::FunctionCount][base_type::FunctionCount],device_type> local_element_matrix_type;

  local_element_array_type val;
  local_element_array_type elem_res;
  local_element_matrix_type elem_mat;

  ElementComputation(const ElementComputation& rhs) :
    base_type(rhs),
    val(rhs.val),
    elem_res(rhs.elem_res),
    elem_mat(rhs.elem_mat) {}

  ElementComputation(
    const typename base_type::mesh_type          & arg_mesh ,
    const typename base_type::coeff_function_type& arg_coeff_function ,
    const typename base_type::vector_type        & arg_solution ,
    const typename base_type::elem_graph_type    & arg_elem_graph ,
    const typename base_type::sparse_matrix_type & arg_jacobian ,
    const typename base_type::vector_type        & arg_residual ) :
    base_type(arg_mesh, arg_coeff_function, arg_solution, arg_elem_graph,
              arg_jacobian, arg_residual),
    val("val", this->elem_node_ids.dimension_0(),base_type::FunctionCount),
    elem_res("elem_res", this->elem_node_ids.dimension_0(),base_type::FunctionCount),
    elem_mat("elem_mat", this->elem_node_ids.dimension_0(),base_type::FunctionCount,base_type::FunctionCount) {}

  //------------------------------------

  void apply() const
  {
    const size_t nelem = this->elem_node_ids.dimension_0();
    parallel_for( nelem , *this );
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( const unsigned ielem ) const
  {
    double x[ base_type::FunctionCount ] ;
    double y[ base_type::FunctionCount ] ;
    double z[ base_type::FunctionCount ] ;
    unsigned node_index[ base_type::ElemNodeCount ];

    // Gather nodal coordinates and solution vector:
    this->gatherSolutionGlobal(ielem, val, node_index, x, y, z, elem_res, elem_mat);

    // Compute nodal element residual vector and Jacobian matrix
    this->computeElementResidualJacobianViewGlobal( ielem, val, x, y, z, elem_res , elem_mat );

    // Scatter nodal element residual and Jacobian in global vector and matrix:
    this->scatterResidualViewGlobal( ielem, node_index, elem_res, elem_mat );
  }
}; /* ElementComputation */

template< class DeviceType , BoxElemPart::ElemOrder Order ,
          class CoordinateMap , typename ScalarType >
class ElementComputation
  < Kokkos::Example::BoxElemFixture< DeviceType , Order , CoordinateMap > ,
    CrsMatrix< ScalarType , DeviceType > ,
    FadLocal > : public ElementComputationBase<DeviceType, Order, CoordinateMap,
                                               ScalarType> {
public:

  typedef ElementComputationBase<DeviceType, Order, CoordinateMap,
                                 ScalarType> base_type;

  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::device_type device_type;
  typedef Sacado::Fad::SFad<scalar_type,base_type::FunctionCount> fad_scalar_type;

  ElementComputation(const ElementComputation& rhs) : base_type(rhs) {}

  ElementComputation(
    const typename base_type::mesh_type          & arg_mesh ,
    const typename base_type::coeff_function_type& arg_coeff_function ,
    const typename base_type::vector_type        & arg_solution ,
    const typename base_type::elem_graph_type    & arg_elem_graph ,
    const typename base_type::sparse_matrix_type & arg_jacobian ,
    const typename base_type::vector_type        & arg_residual ) :
    base_type(arg_mesh, arg_coeff_function, arg_solution, arg_elem_graph,
              arg_jacobian, arg_residual) {}

  //------------------------------------

  void apply() const
  {
    const size_t nelem = this->elem_node_ids.dimension_0();
    parallel_for( nelem , *this );
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( const unsigned ielem ) const
  {
    double x[ base_type::FunctionCount ] ;
    double y[ base_type::FunctionCount ] ;
    double z[ base_type::FunctionCount ] ;
    unsigned node_index[ base_type::ElemNodeCount ];

    fad_scalar_type val[ base_type::FunctionCount ] ;
    fad_scalar_type elem_res[ base_type::FunctionCount ] ;

    // Gather nodal coordinates and solution vector:
    this->gatherSolutionFad( ielem, val, node_index, x, y, z, elem_res );

    // Compute nodal element residual vector:
    this->computeElementResidual( val, x, y, z, elem_res );

    // Scatter nodal element residual in global vector:
    this->scatterResidualFad( ielem, node_index, elem_res );
  }
}; /* ElementComputation */

template< class DeviceType , BoxElemPart::ElemOrder Order ,
          class CoordinateMap , typename ScalarType >
class ElementComputation
  < Kokkos::Example::BoxElemFixture< DeviceType , Order , CoordinateMap > ,
    CrsMatrix< ScalarType , DeviceType > ,
    FadLocalView > :
    public ElementComputationBase<DeviceType, Order, CoordinateMap,
                                  ScalarType> {
public:

  typedef ElementComputationBase<DeviceType, Order, CoordinateMap,
                                 ScalarType> base_type;

  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::device_type device_type;
  typedef Sacado::Fad::SFad<scalar_type,base_type::FunctionCount> fad_scalar_type;

  static const unsigned data_size =
      sizeof(fad_scalar_type[base_type::FunctionCount])/sizeof(scalar_type);

  ElementComputation(const ElementComputation& rhs) : base_type(rhs) {}

  ElementComputation(
    const typename base_type::mesh_type          & arg_mesh ,
    const typename base_type::coeff_function_type& arg_coeff_function ,
    const typename base_type::vector_type        & arg_solution ,
    const typename base_type::elem_graph_type    & arg_elem_graph ,
    const typename base_type::sparse_matrix_type & arg_jacobian ,
    const typename base_type::vector_type        & arg_residual ) :
    base_type(arg_mesh, arg_coeff_function, arg_solution, arg_elem_graph,
              arg_jacobian, arg_residual) {}

  //------------------------------------

  void apply() const
  {
    const size_t nelem = this->elem_node_ids.dimension_0();
    parallel_for( nelem , *this );
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( const unsigned ielem ) const
  {
    typedef Kokkos::View<fad_scalar_type[base_type::FunctionCount],device_type,Kokkos::MemoryUnmanaged> local_element_array_type;

    double x[ base_type::FunctionCount ] ;
    double y[ base_type::FunctionCount ] ;
    double z[ base_type::FunctionCount ] ;
    unsigned node_index[ base_type::ElemNodeCount ];

    scalar_type val_data[data_size];
    scalar_type elem_res_data[data_size];

    local_element_array_type val(val_data, base_type::FunctionCount, base_type::FunctionCount+1);
    local_element_array_type elem_res(elem_res_data, base_type::FunctionCount, base_type::FunctionCount+1);

    // Gather nodal coordinates and solution vector:
    this->gatherSolutionFad( ielem, val, node_index, x, y, z, elem_res );

    // Compute nodal element residual vector:
    this->computeElementResidual( val, x, y, z, elem_res );

    // Scatter nodal element residual in global vector:
    this->scatterResidualFad( ielem, node_index, elem_res );
  }
}; /* ElementComputation */

template< class DeviceType , BoxElemPart::ElemOrder Order ,
          class CoordinateMap , typename ScalarType >
class ElementComputation
  < Kokkos::Example::BoxElemFixture< DeviceType , Order , CoordinateMap > ,
    CrsMatrix< ScalarType , DeviceType >,
    FadGlobalView > :
    public ElementComputationBase<DeviceType, Order, CoordinateMap,
                                  ScalarType> {
public:

  typedef ElementComputationBase<DeviceType, Order, CoordinateMap,
                                 ScalarType> base_type;

  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::device_type device_type;
  typedef Sacado::Fad::SFad<scalar_type,base_type::FunctionCount> fad_scalar_type;
  typedef Kokkos::View<fad_scalar_type*[base_type::FunctionCount],device_type> local_element_array_type;

  local_element_array_type val;
  local_element_array_type elem_res;

  ElementComputation(const ElementComputation& rhs) :
    base_type(rhs),
    val(rhs.val),
    elem_res(rhs.elem_res) {}

  ElementComputation(
    const typename base_type::mesh_type          & arg_mesh ,
    const typename base_type::coeff_function_type& arg_coeff_function ,
    const typename base_type::vector_type        & arg_solution ,
    const typename base_type::elem_graph_type    & arg_elem_graph ,
    const typename base_type::sparse_matrix_type & arg_jacobian ,
    const typename base_type::vector_type        & arg_residual ) :
    base_type(arg_mesh, arg_coeff_function, arg_solution, arg_elem_graph,
              arg_jacobian, arg_residual),
    val("val", this->elem_node_ids.dimension_0(),base_type::FunctionCount,base_type::FunctionCount+1),
    elem_res("elem_res", this->elem_node_ids.dimension_0(),base_type::FunctionCount,base_type::FunctionCount+1) {}

  //------------------------------------

  void apply() const
  {
    const size_t nelem = this->elem_node_ids.dimension_0();
    parallel_for( nelem , *this );
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( const unsigned ielem ) const
  {
    double x[ base_type::FunctionCount ] ;
    double y[ base_type::FunctionCount ] ;
    double z[ base_type::FunctionCount ] ;
    unsigned node_index[ base_type::ElemNodeCount ];

    // Gather nodal coordinates and solution vector:
    this->gatherSolutionFadGlobal( ielem, val, node_index, x, y, z, elem_res );

    // Compute nodal element residual vector:
    this->computeElementResidualGlobal( ielem, val, x, y, z, elem_res );

    // Scatter nodal element residual in global vector:
    this->scatterResidualFadGlobal( ielem, node_index, elem_res );
  }
}; /* ElementComputation */

} /* namespace FENL */
} /* namespace Example */
} /* namespace Kokkos  */

//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_EXAMPLE_FENLFUNCTORS_HPP */
