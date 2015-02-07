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

#include <Kokkos_Pair.hpp>
#include <Kokkos_UnorderedMap.hpp>

#include <impl/Kokkos_Timer.hpp>

#include <fenl.hpp>
#include <BoxElemFixture.hpp>
#include <HexElement.hpp>

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Example {
namespace FENL {

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
    , row_count(Kokkos::ViewAllocateWithoutInitializing("row_count"), node_count ) // will deep_copy to 0 inside loop
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

      row_total = UnsignedValue();
      row_count = RowMapType();
      row_map   = RowMapType();
      node_node_set.clear();

      //--------------------------------

      phase = SORT_GRAPH_ENTRIES ;
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

template< class FiniteElementMeshType , class SparseMatrixType
        , class CoeffFunctionType = ElementComputationConstantCoefficient
        >
class ElementComputation ;


template< class DeviceType , BoxElemPart::ElemOrder Order , class CoordinateMap ,
          typename ScalarType , typename OrdinalType , class MemoryTraits , typename SizeType ,
          class CoeffFunctionType >
class ElementComputation
  < Kokkos::Example::BoxElemFixture< DeviceType , Order , CoordinateMap >
  , Kokkos::CrsMatrix< ScalarType , OrdinalType , DeviceType , MemoryTraits , SizeType >
  , CoeffFunctionType >
{
public:

  typedef Kokkos::Example::BoxElemFixture< DeviceType, Order, CoordinateMap >  mesh_type ;
  typedef Kokkos::Example::HexElement_Data< mesh_type::ElemNode >              element_data_type ;

  //------------------------------------

  typedef DeviceType   device_type ;
  typedef ScalarType   scalar_type ;

  typedef Kokkos::CrsMatrix< ScalarType , OrdinalType , DeviceType , MemoryTraits , SizeType >  sparse_matrix_type ;
  typedef typename sparse_matrix_type::StaticCrsGraphType                                       sparse_graph_type ;
  typedef typename sparse_matrix_type::values_type matrix_values_type ;
  typedef Kokkos::View< scalar_type* , Kokkos::LayoutLeft, device_type > vector_type ;

  //------------------------------------

  typedef LocalViewTraits< vector_type > local_vector_view_traits;
  typedef LocalViewTraits< matrix_values_type> local_matrix_view_traits;
  typedef typename local_vector_view_traits::local_view_type local_vector_type;
  typedef typename local_matrix_view_traits::local_view_type local_matrix_type;
  typedef typename local_vector_view_traits::local_value_type local_scalar_type;
  static const bool use_team = local_vector_view_traits::use_team;

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

  typedef LocalViewTraits< elem_matrices_type > local_elem_matrices_traits;
  typedef LocalViewTraits< elem_vectors_type > local_elem_vectors_traits;
  typedef typename local_elem_matrices_traits::local_view_type local_elem_matrices_type;
  typedef typename local_elem_vectors_traits::local_view_type local_elem_vectors_type;

  typedef typename NodeNodeGraph< elem_node_type , sparse_graph_type , ElemNodeCount >::ElemGraphType elem_graph_type ;

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
  const CoeffFunctionType   coeff_function ;
  const double              coeff_source ;
  const double              coeff_advection ;
  const Kokkos::DeviceConfig dev_config ;

  ElementComputation( const ElementComputation & rhs )
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
    , coeff_source( rhs.coeff_source )
    , coeff_advection( rhs.coeff_advection )
    , dev_config( rhs.dev_config )
    {}

  // If the element->sparse_matrix graph is provided then perform atomic updates
  // Otherwise fill per-element contributions for subequent gather-add into a residual and jacobian.
  ElementComputation( const mesh_type          & arg_mesh ,
                      const CoeffFunctionType  & arg_coeff_function ,
                      const double             & arg_coeff_source ,
                      const double             & arg_coeff_advection ,
                      const vector_type        & arg_solution ,
                      const elem_graph_type    & arg_elem_graph ,
                      const sparse_matrix_type & arg_jacobian ,
                      const vector_type        & arg_residual ,
                      const Kokkos::DeviceConfig arg_dev_config ,
                      const QuadratureData<DeviceType>& qd =
                        QuadratureData<DeviceType>() )
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
    , coeff_source( arg_coeff_source )
    , coeff_advection( arg_coeff_advection )
    , dev_config( arg_dev_config )
    {}

  //------------------------------------

  void apply() const
  {
    const size_t nelem = elem_node_ids.dimension_0();
    if ( use_team ) {
      const size_t team_size = dev_config.block_dim.x * dev_config.block_dim.y;
      const size_t league_size =
        (nelem + dev_config.block_dim.y-1) / dev_config.block_dim.y;
      Kokkos::TeamPolicy< device_type > config( league_size, team_size );
      parallel_for( config , *this );
    }
    else {
      parallel_for( nelem , *this );
    }
  }

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

  KOKKOS_INLINE_FUNCTION
  void contributeResidualJacobian(
    const local_scalar_type dof_values[] ,
    const double  dpsidx[] ,
    const double  dpsidy[] ,
    const double  dpsidz[] ,
    const double  detJ ,
    const double  integ_weight ,
    const double  bases_vals[] ,
    const local_scalar_type  coeff_k ,
    const double  coeff_src ,
    const double  advection[] ,
    local_scalar_type  elem_res[] ,
    local_scalar_type  elem_mat[][ FunctionCount ] ) const
  {
    local_scalar_type value_at_pt = 0 ;
    local_scalar_type gradx_at_pt = 0 ;
    local_scalar_type grady_at_pt = 0 ;
    local_scalar_type gradz_at_pt = 0 ;

    for ( unsigned m = 0 ; m < FunctionCount ; m++ ) {
      value_at_pt += dof_values[m] * bases_vals[m] ;
      gradx_at_pt += dof_values[m] * dpsidx[m] ;
      grady_at_pt += dof_values[m] * dpsidy[m] ;
      gradz_at_pt += dof_values[m] * dpsidz[m] ;
    }

    const double detJ_weight = detJ * integ_weight;

    const local_scalar_type source_term =
      coeff_src * value_at_pt * value_at_pt;
    const local_scalar_type source_deriv =
      2.0 * coeff_src * value_at_pt;

    const local_scalar_type advection_term =
      advection[0]*gradx_at_pt +
      advection[1]*grady_at_pt +
      advection[2]*gradz_at_pt;

    for ( unsigned m = 0; m < FunctionCount; ++m) {
      local_scalar_type * const mat = elem_mat[m] ;
      const double bases_val_m = bases_vals[m];
      const double dpsidx_m    = dpsidx[m] ;
      const double dpsidy_m    = dpsidy[m] ;
      const double dpsidz_m    = dpsidz[m] ;

      elem_res[m] +=
        detJ_weight * ( coeff_k * ( dpsidx_m * gradx_at_pt +
                                    dpsidy_m * grady_at_pt +
                                    dpsidz_m * gradz_at_pt ) +
                        ( advection_term  + source_term ) * bases_val_m ) ;

      for( unsigned n = 0; n < FunctionCount; n++) {

        mat[n] +=
          detJ_weight * ( coeff_k * ( dpsidx_m * dpsidx[n] +
                                      dpsidy_m * dpsidy[n] +
                                      dpsidz_m * dpsidz[n] ) +
                          ( advection[0] * dpsidx[n] +
                            advection[1] * dpsidy[n] +
                            advection[2] * dpsidz[n] +
                            source_deriv * bases_vals[n] ) * bases_val_m );
      }
    }
  }

  typedef typename Kokkos::TeamPolicy< device_type >::member_type team_member ;
  KOKKOS_INLINE_FUNCTION
  void operator()( const team_member & dev ) const
  {

    const unsigned num_ensemble_threads = dev_config.block_dim.x ;
    const unsigned num_element_threads  = dev_config.block_dim.y ;
    const unsigned element_rank  = dev.team_rank() / num_ensemble_threads ;
    const unsigned ensemble_rank = dev.team_rank() % num_ensemble_threads ;

    const unsigned ielem =
      dev.league_rank() * num_element_threads + element_rank;

    if (ielem >= elem_node_ids.dimension_0())
      return;

    (*this)( ielem, ensemble_rank );

  }

  KOKKOS_INLINE_FUNCTION
  void operator()( const unsigned ielem ,
                   const unsigned ensemble_rank = 0 ) const
  {
    local_vector_type local_solution =
      local_vector_view_traits::create_local_view(solution,
                                                  ensemble_rank);
    local_vector_type local_residual =
      local_vector_view_traits::create_local_view(residual,
                                                  ensemble_rank);
    local_matrix_type local_jacobian_values =
      local_matrix_view_traits::create_local_view(jacobian.values,
                                                  ensemble_rank);

    // Gather nodal coordinates and solution vector:

    double x[ FunctionCount ] ;
    double y[ FunctionCount ] ;
    double z[ FunctionCount ] ;
    local_scalar_type val[ FunctionCount ] ;
    unsigned node_index[ ElemNodeCount ];

    for ( unsigned i = 0 ; i < ElemNodeCount ; ++i ) {
      const unsigned ni = elem_node_ids( ielem , i );

      node_index[i] = ni ;

      x[i] = node_coords( ni , 0 );
      y[i] = node_coords( ni , 1 );
      z[i] = node_coords( ni , 2 );

      val[i] = local_solution( ni );
    }


    local_scalar_type elem_vec[ FunctionCount ] ;
    local_scalar_type elem_mat[ FunctionCount ][ FunctionCount ] ;

    for( unsigned i = 0; i < FunctionCount ; i++ ) {
      elem_vec[i] = 0 ;
      for( unsigned j = 0; j < FunctionCount ; j++){
        elem_mat[i][j] = 0 ;
      }
    }

    // advection = [ 1 , 1 , 1 ] * coeff_advection
    double advection[] = { coeff_advection ,
                           coeff_advection ,
                           coeff_advection };

    for ( unsigned i = 0 ; i < IntegrationCount ; ++i ) {
      double dpsidx[ FunctionCount ] ;
      double dpsidy[ FunctionCount ] ;
      double dpsidz[ FunctionCount ] ;

      const double detJ =
        transform_gradients( elem_data.gradients[i] , x , y , z ,
                             dpsidx , dpsidy , dpsidz );

      // Compute physical coordinates of integration point
      double pt[] = { 0 , 0 , 0 } ;
      for ( unsigned j = 0 ; j < FunctionCount ; ++j ) {
        pt[0] += x[j] * elem_data.values[i][j] ;
        pt[1] += y[j] * elem_data.values[i][j] ;
        pt[2] += z[j] * elem_data.values[i][j] ;
      }

      // Evaluate diffusion coefficient
      local_scalar_type coeff_k = coeff_function(pt, ensemble_rank);

      contributeResidualJacobian( val , dpsidx , dpsidy , dpsidz , detJ ,
                                  elem_data.weights[i] , elem_data.values[i] ,
                                  coeff_k , coeff_source , advection ,
                                  elem_vec , elem_mat );
    }

    for( unsigned i = 0 ; i < FunctionCount ; i++ ) {
      const unsigned row = node_index[i] ;
      if ( row < residual.dimension_0() ) {
        atomic_add( & local_residual( row ) , elem_vec[i] );

        for( unsigned j = 0 ; j < FunctionCount ; j++ ) {
          const unsigned entry = elem_graph( ielem , i , j );
          if ( entry != ~0u ) {
            atomic_add( & local_jacobian_values( entry ) , elem_mat[i][j] );
          }
        }
      }
    }
  }
}; /* ElementComputation */

//----------------------------------------------------------------------------

template< class FixtureType , class SparseMatrixType >
class DirichletComputation ;

template< class DeviceType , BoxElemPart::ElemOrder Order , class CoordinateMap ,
          typename ScalarType , typename OrdinalType , class MemoryTraits , typename SizeType >
class DirichletComputation<
  Kokkos::Example::BoxElemFixture< DeviceType , Order , CoordinateMap > ,
  Kokkos::CrsMatrix< ScalarType , OrdinalType , DeviceType , MemoryTraits , SizeType > >
{
public:

  typedef Kokkos::Example::BoxElemFixture< DeviceType, Order, CoordinateMap >  mesh_type ;
  typedef typename mesh_type::node_coord_type                                  node_coord_type ;
  typedef typename node_coord_type::value_type                                 scalar_coord_type ;

  typedef DeviceType   device_type ;
  typedef ScalarType   scalar_type ;

  typedef Kokkos::CrsMatrix< ScalarType , OrdinalType , DeviceType , MemoryTraits , SizeType >  sparse_matrix_type ;
  typedef typename sparse_matrix_type::StaticCrsGraphType                                       sparse_graph_type ;
  typedef typename sparse_matrix_type::values_type matrix_values_type ;
  typedef Kokkos::View< scalar_type* , device_type > vector_type ;

  //------------------------------------

  typedef LocalViewTraits< vector_type > local_vector_view_traits;
  typedef LocalViewTraits< matrix_values_type> local_matrix_view_traits;
  typedef typename local_vector_view_traits::local_view_type local_vector_type;
  typedef typename local_matrix_view_traits::local_view_type local_matrix_type;
  static const bool use_team = local_vector_view_traits::use_team;

  typedef double       bc_scalar_type ;

  //------------------------------------
  // Computational data:

  const node_coord_type     node_coords ;
  const vector_type         solution ;
  const sparse_matrix_type  jacobian ;
  const vector_type         residual ;
  const bc_scalar_type      bc_lower_value ;
  const bc_scalar_type      bc_upper_value ;
  const scalar_coord_type   bc_lower_limit ;
  const scalar_coord_type   bc_upper_limit ;
  const unsigned            bc_plane ;
  const unsigned            node_count ;
        bool                init ;
  const Kokkos::DeviceConfig dev_config ;


  DirichletComputation( const mesh_type          & arg_mesh ,
                        const vector_type        & arg_solution ,
                        const sparse_matrix_type & arg_jacobian ,
                        const vector_type        & arg_residual ,
                        const unsigned             arg_bc_plane ,
                        const bc_scalar_type       arg_bc_lower_value ,
                        const bc_scalar_type       arg_bc_upper_value ,
                        const Kokkos::DeviceConfig arg_dev_config )
    : node_coords( arg_mesh.node_coord() )
    , solution(    arg_solution )
    , jacobian(    arg_jacobian )
    , residual(    arg_residual )
    , bc_lower_value( arg_bc_lower_value )
    , bc_upper_value( arg_bc_upper_value )
    , bc_lower_limit( std::numeric_limits<scalar_coord_type>::epsilon() )
    , bc_upper_limit( scalar_coord_type(1) - std::numeric_limits<scalar_coord_type>::epsilon() )
    , bc_plane(       arg_bc_plane )
    , node_count( arg_mesh.node_count_owned() )
    , init( false )
    , dev_config( arg_dev_config )
    {
      parallel_for( node_count , *this );
      init = true ;
    }

  void apply() const
  {
    if ( use_team ) {
      const size_t team_size = dev_config.block_dim.x * dev_config.block_dim.y;
      const size_t league_size =
        (node_count + dev_config.block_dim.y-1) / dev_config.block_dim.y;
      Kokkos::TeamPolicy< device_type > config( league_size, team_size );
      parallel_for( config , *this );
    }
    else
      parallel_for( node_count , *this );
  }

  //------------------------------------

  typedef typename Kokkos::TeamPolicy< device_type >::member_type team_member ;
  KOKKOS_INLINE_FUNCTION
  void operator()( const team_member & dev ) const
  {

    const unsigned num_ensemble_threads = dev_config.block_dim.x ;
    const unsigned num_node_threads     = dev_config.block_dim.y ;
    const unsigned node_rank     = dev.team_rank() / num_ensemble_threads ;
    const unsigned ensemble_rank = dev.team_rank() % num_ensemble_threads ;

    const unsigned inode = dev.league_rank() * num_node_threads + node_rank;

    if (inode >= node_count)
      return;

    (*this)( inode, ensemble_rank );

  }

  KOKKOS_INLINE_FUNCTION
  void operator()( const unsigned inode ,
                   const unsigned ensemble_rank = 0) const
  {
    local_vector_type local_residual =
      local_vector_view_traits::create_local_view(residual,
                                                  ensemble_rank);
    local_matrix_type local_jacobian_values =
      local_matrix_view_traits::create_local_view(jacobian.values,
                                                  ensemble_rank);

    //  Apply dirichlet boundary condition on the Solution and Residual vectors.
    //  To maintain the symmetry of the original global stiffness matrix,
    //  zero out the columns that correspond to boundary conditions, and
    //  update the residual vector accordingly

    const unsigned iBeg = jacobian.graph.row_map[inode];
    const unsigned iEnd = jacobian.graph.row_map[inode+1];

    const scalar_coord_type c = node_coords(inode,bc_plane);
    const bool bc_lower = c <= bc_lower_limit ;
    const bool bc_upper = bc_upper_limit <= c ;

    if ( ! init ) {
      solution(inode) = bc_lower ? bc_lower_value : (
                        bc_upper ? bc_upper_value : 0 );
    }
    else {
      if ( bc_lower || bc_upper ) {

        local_residual(inode) = 0 ;

        //  zero each value on the row, and leave a one
        //  on the diagonal

        for( unsigned i = iBeg ; i < iEnd ; ++i ) {
          local_jacobian_values(i) = int(inode) == int(jacobian.graph.entries(i)) ? 1 : 0 ;
        }
      }
      else {

        //  Find any columns that are boundary conditions.
        //  Clear them and adjust the residual vector

        for( unsigned i = iBeg ; i < iEnd ; ++i ) {
          const unsigned       cnode = jacobian.graph.entries(i) ;
          const scalar_coord_type cc = node_coords(cnode,bc_plane);

          if ( ( cc <= bc_lower_limit ) || ( bc_upper_limit <= cc ) ) {
            local_jacobian_values(i) = 0 ;
          }
        }
      }
    }
  }
};

template< typename FixtureType , typename VectorType >
class ResponseComputation
{
public:

  typedef FixtureType fixture_type ;
  typedef VectorType vector_type ;
  typedef typename vector_type::device_type device_type ;
  typedef typename vector_type::value_type value_type ;

  typedef Kokkos::Example::HexElement_Data< fixture_type::ElemNode > element_data_type ;
  static const unsigned SpatialDim       = element_data_type::spatial_dimension ;
  static const unsigned TensorDim        = SpatialDim * SpatialDim ;
  static const unsigned ElemNodeCount    = element_data_type::element_node_count ;
  static const unsigned IntegrationCount = element_data_type::integration_count ;

  //------------------------------------
  // Computational data:

  const element_data_type    elem_data ;
  const fixture_type         fixture ;
  const vector_type          solution ;

  ResponseComputation( const ResponseComputation & rhs )
    : elem_data()
    , fixture( rhs.fixture )
    , solution( rhs.solution )
    {}

  ResponseComputation( const fixture_type& arg_fixture ,
                       const vector_type & arg_solution )
    : elem_data()
    , fixture( arg_fixture )
    , solution( arg_solution )
    {}

  //------------------------------------

  value_type apply() const
  {
    value_type response = 0;
    //Kokkos::parallel_reduce( fixture.elem_count() , *this , response );
    Kokkos::parallel_reduce( solution.dimension_0() , *this , response );
    return response;
  }

  //------------------------------------

   KOKKOS_INLINE_FUNCTION
  double compute_detJ(
    const double grad[][ ElemNodeCount ] , // Gradient of bases master element
    const double x[] ,
    const double y[] ,
    const double z[] ) const
  {
    enum { j11 = 0 , j12 = 1 , j13 = 2 ,
           j21 = 3 , j22 = 4 , j23 = 5 ,
           j31 = 6 , j32 = 7 , j33 = 8 };

    // Jacobian accumulation:

    double J[ TensorDim ] = { 0, 0, 0,  0, 0, 0,  0, 0, 0 };

    for( unsigned i = 0; i < ElemNodeCount ; ++i ) {
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

    return detJ ;
  }

  KOKKOS_INLINE_FUNCTION
  value_type contributeResponse(
    const value_type dof_values[] ,
    const double  detJ ,
    const double  integ_weight ,
    const double  bases_vals[] ) const
  {
    // $$ g_i = \int_{\Omega} T^2 d \Omega $$

    value_type value_at_pt = 0 ;
    for ( unsigned m = 0 ; m < ElemNodeCount ; m++ ) {
      value_at_pt += dof_values[m] * bases_vals[m] ;
    }

    value_type elem_response =
      value_at_pt * value_at_pt * detJ * integ_weight ;

    return elem_response;
  }

  /*
  KOKKOS_INLINE_FUNCTION
  void operator()( const unsigned ielem , value_type& response ) const
  {
    // Gather nodal coordinates and solution vector:

    double x[ ElemNodeCount ] ;
    double y[ ElemNodeCount ] ;
    double z[ ElemNodeCount ] ;
    value_type val[ ElemNodeCount ] ;

    for ( unsigned i = 0 ; i < ElemNodeCount ; ++i ) {
      const unsigned ni = fixture.elem_node( ielem , i );

      x[i] = fixture.node_coord( ni , 0 );
      y[i] = fixture.node_coord( ni , 1 );
      z[i] = fixture.node_coord( ni , 2 );

      val[i] = solution( ni );
    }

    for ( unsigned i = 0 ; i < IntegrationCount ; ++i ) {

      const double detJ = compute_detJ( elem_data.gradients[i] , x , y , z );

      response += contributeResponse( val , detJ , elem_data.weights[i] ,
                                      elem_data.values[i] );
    }
  }
  */

  KOKKOS_INLINE_FUNCTION
  void operator()( const unsigned i , value_type& response ) const
  {
    const value_type& u = solution(i);
    response += (u * u) / fixture.node_count_global();
  }

  KOKKOS_INLINE_FUNCTION
  void init( value_type & response ) const
  { response = 0 ; }

  KOKKOS_INLINE_FUNCTION
  void join( volatile value_type & response ,
             volatile const value_type & input ) const
  { response += input ; }

}; /* ResponseComputation */

} /* namespace FENL */
} /* namespace Example */
} /* namespace Kokkos  */

//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_EXAMPLE_FENLFUNCTORS_HPP */
