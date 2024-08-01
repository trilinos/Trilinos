//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

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

#include <Kokkos_Timer.hpp>

#include <BoxElemFixture.hpp>
#include <HexElement.hpp>

#include "Sacado.hpp"

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

namespace Kokkos {
namespace Example {
namespace FENL {

template< typename ValueType , class Space >
struct CrsMatrix {
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE // Don't remove this until Kokkos has removed the deprecated code path probably around September 2018
  typedef Kokkos::StaticCrsGraph< unsigned , Space , void , unsigned >  StaticCrsGraphType ;
#else
  typedef Kokkos::StaticCrsGraph< unsigned , Space , void , void , unsigned >  StaticCrsGraphType ;
#endif
  typedef View< ValueType * , Space > coeff_type ;

  StaticCrsGraphType  graph ;
  coeff_type          coeff ;

  CrsMatrix() : graph(), coeff() {}

  CrsMatrix( const StaticCrsGraphType & arg_graph )
    : graph( arg_graph )
    , coeff( "crs_matrix_coeff" , arg_graph.entries.extent(0) )
    {}
};

template< class ElemNodeIdView , class CrsGraphType , unsigned ElemNode >
class NodeNodeGraph {
public:

  typedef typename ElemNodeIdView::execution_space  execution_space ;
  typedef pair<unsigned,unsigned> key_type ;

  typedef Kokkos::UnorderedMap< key_type, void , execution_space >  SetType ;
  typedef typename CrsGraphType::row_map_type::non_const_type       RowMapType ;
  typedef Kokkos::View< unsigned ,  execution_space >               UnsignedValue ;

  // Static dimensions of 0 generate compiler warnings or errors.
  typedef Kokkos::View< unsigned*[ElemNode][ElemNode] , execution_space >
    ElemGraphType ;

  struct TagFillNodeSet {};
  struct TagScanNodeCount {};
  struct TagFillGraphEntries {};
  struct TagSortGraphEntries {};
  struct TagFillElementGraph {};

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

      Kokkos::Timer wall_clock ;

      wall_clock.reset();
      phase = FILL_NODE_SET ;

      // upper bound on the capacity
      size_t set_capacity = (28ull * node_count) / 2;
      unsigned failed_insert_count = 0 ;

      do {
        // Zero the row count to restart the fill
        Kokkos::deep_copy( row_count , 0u );

        node_node_set = SetType( ( set_capacity += failed_insert_count ) );

        // May be larger that requested:
        set_capacity = node_node_set.capacity();

        Kokkos::parallel_reduce( Kokkos::RangePolicy<execution_space,TagFillNodeSet>(0,elem_node_id.extent(0))
                               , *this
                               , failed_insert_count );

      } while ( failed_insert_count );

      execution_space().fence();
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

      execution_space().fence();
      results.scan_node_count = wall_clock.seconds();

      wall_clock.reset();
      phase = FILL_GRAPH_ENTRIES ;
      Kokkos::parallel_for( node_node_set.capacity() , *this );

      execution_space().fence();
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

      execution_space().fence();
      results.sort_graph_entries = wall_clock.seconds();

      //--------------------------------
      // Element-to-graph mapping:
      wall_clock.reset();
      phase = FILL_ELEMENT_GRAPH ;
      elem_graph = ElemGraphType("elem_graph", elem_node_id.extent(0) );
      Kokkos::parallel_for( elem_node_id.extent(0) , *this );

      execution_space().fence();
      results.fill_element_graph = wall_clock.seconds();
    }

  //------------------------------------
  // parallel_for: create map and count row length

  KOKKOS_INLINE_FUNCTION
  void operator()( const TagFillNodeSet & , unsigned ielem , unsigned & count ) const
  {
    // Loop over element's (row_local_node,col_local_node) pairs:
    for ( unsigned row_local_node = 0 ; row_local_node < elem_node_id.extent(1) ; ++row_local_node ) {

      const unsigned row_node = elem_node_id( ielem , row_local_node );

      for ( unsigned col_local_node = row_local_node ; col_local_node < elem_node_id.extent(1) ; ++col_local_node ) {

        const unsigned col_node = elem_node_id( ielem , col_local_node );

        // If either node is locally owned then insert the pair into the unordered map:

        if ( row_node < row_count.extent(0) || col_node < row_count.extent(0) ) {

          const key_type key = (row_node < col_node) ? make_pair( row_node, col_node ) : make_pair( col_node, row_node ) ;

          const typename SetType::insert_result result = node_node_set.insert( key );

          // A successfull insert: the first time this pair was added
          if ( result.success() ) {

            // If row node is owned then increment count
            if ( row_node < row_count.extent(0) ) { atomic_increment( & row_count( row_node ) ); }

            // If column node is owned and not equal to row node then increment count
            if ( col_node < row_count.extent(0) && col_node != row_node ) { atomic_increment( & row_count( col_node ) ); }
          }
          else if ( result.failed() ) {
            ++count ;
          }
        }
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void fill_graph_entries( const unsigned iset ) const
  {
    typedef typename std::remove_reference< decltype( row_count(0) ) >::type atomic_incr_type;

    if ( node_node_set.valid_at(iset) ) {
      // Add each entry to the graph entries.

      const key_type key = node_node_set.key_at(iset) ;
      const unsigned row_node = key.first ;
      const unsigned col_node = key.second ;

      if ( row_node < row_count.extent(0) ) {
        const unsigned offset = graph.row_map( row_node ) + atomic_fetch_add( & row_count( row_node ) , atomic_incr_type(1) );
        graph.entries( offset ) = col_node ;
      }

      if ( col_node < row_count.extent(0) && col_node != row_node ) {
        const unsigned offset = graph.row_map( col_node ) + atomic_fetch_add( & row_count( col_node ) , atomic_incr_type(1) );
        graph.entries( offset ) = row_node ;
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void sort_graph_entries( const unsigned irow ) const
  {
    const unsigned row_beg = graph.row_map( irow );
    const unsigned row_end = graph.row_map( irow + 1 );
    for ( unsigned i = row_beg + 1 ; i < row_end ; ++i ) {
      const unsigned col = graph.entries(i);
      unsigned j = i ;
      for ( ; row_beg < j && col < graph.entries(j-1) ; --j ) {
        graph.entries(j) = graph.entries(j-1);
      }
      graph.entries(j) = col ;
    }
  }

  KOKKOS_INLINE_FUNCTION
  void fill_elem_graph_map( const unsigned ielem ) const
  {
    for ( unsigned row_local_node = 0 ; row_local_node < elem_node_id.extent(1) ; ++row_local_node ) {

      const unsigned row_node = elem_node_id( ielem , row_local_node );

      for ( unsigned col_local_node = 0 ; col_local_node < elem_node_id.extent(1) ; ++col_local_node ) {

        const unsigned col_node = elem_node_id( ielem , col_local_node );

        unsigned entry = ~0u ;

        if ( row_node + 1 < graph.row_map.extent(0) ) {

          const unsigned entry_end = graph.row_map( row_node + 1 );

          entry = graph.row_map( row_node );

          for ( ; entry < entry_end && graph.entries(entry) != col_node ; ++entry );

          if ( entry == entry_end ) entry = ~0u ;
        }

        elem_graph( ielem , row_local_node , col_local_node ) = entry ;
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( const unsigned iwork ) const
  {
/*
    if ( phase == FILL_NODE_SET ) {
      operator()( TagFillNodeSet() , iwork );
    }
    else */  
    if ( phase == FILL_GRAPH_ENTRIES ) {
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
      if ( irow + 1 == row_count.extent(0) ) {
        row_map( irow + 1 ) = update ;
        row_total()         = update ;
      }
    }
  }

  // For the reduce phase:
  KOKKOS_INLINE_FUNCTION
  void init( const TagFillNodeSet & , unsigned & update ) const { update = 0 ; }

  KOKKOS_INLINE_FUNCTION
  void join( const TagFillNodeSet & 
           , unsigned & update
           , const unsigned & input ) const { update += input ; }

  // For the scan phase::
  KOKKOS_INLINE_FUNCTION
  void init( unsigned & update ) const { update = 0 ; }

  KOKKOS_INLINE_FUNCTION
  void join( unsigned & update
           , const unsigned & input ) const { update += input ; }

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

template< class ExecutionSpace , BoxElemPart::ElemOrder Order ,
          class CoordinateMap , typename ScalarType >
class ElementComputationBase
{
public:

  typedef Kokkos::Example::BoxElemFixture< ExecutionSpace, Order, CoordinateMap >  mesh_type ;
  typedef Kokkos::Example::HexElement_Data< mesh_type::ElemNode >              element_data_type ;

  //------------------------------------

  typedef ExecutionSpace   execution_space ;
  typedef ScalarType   scalar_type ;

  typedef CrsMatrix< ScalarType , ExecutionSpace >  sparse_matrix_type ;
  typedef typename sparse_matrix_type::StaticCrsGraphType                                       sparse_graph_type ;
  typedef Kokkos::View< scalar_type* , Kokkos::LayoutLeft, execution_space > vector_type ;

  //------------------------------------

  static const unsigned SpatialDim       = element_data_type::spatial_dimension ;
  static const unsigned TensorDim        = SpatialDim * SpatialDim ;
  static const unsigned ElemNodeCount    = element_data_type::element_node_count ;
  static const unsigned FunctionCount    = element_data_type::function_count ;
  static const unsigned IntegrationCount = element_data_type::integration_count ;

  //------------------------------------

  typedef typename mesh_type::node_coord_type                                      node_coord_type ;
  typedef typename mesh_type::elem_node_type                                       elem_node_type ;
  typedef Kokkos::View< scalar_type*[FunctionCount][FunctionCount] , execution_space > elem_matrices_type ;
  typedef Kokkos::View< scalar_type*[FunctionCount] ,                execution_space > elem_vectors_type ;

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
    {}

  ElementComputationBase( const mesh_type          & arg_mesh ,
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

};

enum AssemblyMethod {
  Analytic,
  FadElement,
  FadElementOptimized,
  FadQuadPoint
};

template< class FiniteElementMeshType ,
          class SparseMatrixType ,
          AssemblyMethod Method
        >
class ElementComputation ;

template< class ExecutionSpace , BoxElemPart::ElemOrder Order ,
          class CoordinateMap , typename ScalarType >
class ElementComputation
  < Kokkos::Example::BoxElemFixture< ExecutionSpace , Order , CoordinateMap > ,
    CrsMatrix< ScalarType , ExecutionSpace > ,
    Analytic > :
    public ElementComputationBase<ExecutionSpace, Order, CoordinateMap,
                                  ScalarType> {
public:

  typedef ElementComputationBase<ExecutionSpace, Order, CoordinateMap,
                                 ScalarType> base_type;

  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::execution_space execution_space;

  static const unsigned FunctionCount = base_type::FunctionCount;
  static const unsigned IntegrationCount = base_type::IntegrationCount;
  static const unsigned ElemNodeCount = base_type::ElemNodeCount;

  ElementComputation(const ElementComputation& rhs) : base_type(rhs) {}

  ElementComputation(
    const typename base_type::mesh_type          & arg_mesh ,
    const typename base_type::vector_type        & arg_solution ,
    const typename base_type::elem_graph_type    & arg_elem_graph ,
    const typename base_type::sparse_matrix_type & arg_jacobian ,
    const typename base_type::vector_type        & arg_residual ) :
    base_type(arg_mesh, arg_solution, arg_elem_graph,
              arg_jacobian, arg_residual) {}

  //------------------------------------

  void apply() const
  {
    const size_t nelem = this->elem_node_ids.extent(0);
    parallel_for( nelem , *this );
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
      const unsigned ni = this->elem_node_ids( ielem , i );

      node_index[i] = ni ;

      x[i] = this->node_coords( ni , 0 );
      y[i] = this->node_coords( ni , 1 );
      z[i] = this->node_coords( ni , 2 );

      val[i] = this->solution( ni ) ;
      res[i] = 0 ;

      for( unsigned j = 0; j < FunctionCount ; j++){
        mat[i][j] = 0 ;
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
      if ( row < this->residual.extent(0) ) {
        atomic_add( & this->residual( row ) , res[i] );

        for( unsigned j = 0 ; j < FunctionCount ; j++ ) {
          const unsigned entry = this->elem_graph( ielem , i , j );
          if ( entry != ~0u ) {
            atomic_add( & this->jacobian.coeff( entry ) , mat[i][j] );
          }
        }
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void computeElementResidualJacobian(
    const scalar_type dof_values[] ,
    const double x[],
    const double y[],
    const double z[],
    scalar_type  elem_res[] ,
    scalar_type  elem_mat[][FunctionCount] ) const
  {
    double coeff_k = 3.456;
    double coeff_src = 1.234;
    double advection[] = { 1.1, 1.2, 1.3 };
    double dpsidx[ FunctionCount ] ;
    double dpsidy[ FunctionCount ] ;
    double dpsidz[ FunctionCount ] ;
    for ( unsigned i = 0 ; i < IntegrationCount ; ++i ) {

      const double integ_weight = this->elem_data.weights[i];
      const double* bases_vals = this->elem_data.values[i];
      const double detJ =
        this->transform_gradients( this->elem_data.gradients[i] ,
                                   x , y , z ,
                                   dpsidx , dpsidy , dpsidz );
      const double detJ_weight = detJ * integ_weight;
      const double detJ_weight_coeff_k = detJ_weight * coeff_k;

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

      const scalar_type source_term =
        coeff_src * value_at_pt * value_at_pt ;
      const scalar_type source_deriv =
        2.0 * coeff_src * value_at_pt ;

      const scalar_type advection_x = advection[0];
      const scalar_type advection_y = advection[1];
      const scalar_type advection_z = advection[2];

      const scalar_type advection_term =
        advection_x*gradx_at_pt +
        advection_y*grady_at_pt +
        advection_z*gradz_at_pt ;

      for ( unsigned m = 0; m < FunctionCount; ++m) {
        scalar_type * const mat = elem_mat[m] ;
        const double bases_val_m = bases_vals[m] * detJ_weight ;
        const double dpsidx_m    = dpsidx[m] ;
        const double dpsidy_m    = dpsidy[m] ;
        const double dpsidz_m    = dpsidz[m] ;

        elem_res[m] +=
          detJ_weight_coeff_k * ( dpsidx_m * gradx_at_pt +
                                  dpsidy_m * grady_at_pt +
                                  dpsidz_m * gradz_at_pt ) +
          bases_val_m * ( advection_term  + source_term ) ;

        for( unsigned n = 0; n < FunctionCount; n++) {
          const double dpsidx_n    = dpsidx[n] ;
          const double dpsidy_n    = dpsidy[n] ;
          const double dpsidz_n    = dpsidz[n] ;
          mat[n] +=
            detJ_weight_coeff_k * ( dpsidx_m * dpsidx_n +
                                    dpsidy_m * dpsidy_n +
                                    dpsidz_m * dpsidz_n ) +
            bases_val_m * ( advection_x * dpsidx_n +
                            advection_y * dpsidy_n +
                            advection_z * dpsidz_n +
                            source_deriv * bases_vals[n] )  ;
        }
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( const unsigned ielem ) const
  {
    double x[ FunctionCount ] ;
    double y[ FunctionCount ] ;
    double z[ FunctionCount ] ;
    unsigned node_index[ ElemNodeCount ];

    scalar_type val[ FunctionCount ] ;
    scalar_type elem_res[ FunctionCount ] ;
    scalar_type elem_mat[ FunctionCount ][ FunctionCount ] ;

    // Gather nodal coordinates and solution vector:
    gatherSolution(ielem, val, node_index, x, y, z, elem_res, elem_mat);

    // Compute nodal element residual vector and Jacobian matrix
    computeElementResidualJacobian( val, x, y, z, elem_res , elem_mat );

    // Scatter nodal element residual and Jacobian in global vector and matrix:
    scatterResidual( ielem, node_index, elem_res, elem_mat );
  }
}; /* ElementComputation */

template< class ExecutionSpace , BoxElemPart::ElemOrder Order ,
          class CoordinateMap , typename ScalarType >
class ElementComputation
  < Kokkos::Example::BoxElemFixture< ExecutionSpace , Order , CoordinateMap > ,
    CrsMatrix< ScalarType , ExecutionSpace > ,
    FadElement > : public ElementComputationBase<ExecutionSpace, Order, CoordinateMap,
                                                 ScalarType> {
public:

  typedef ElementComputationBase<ExecutionSpace, Order, CoordinateMap,
                                 ScalarType> base_type;

  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::execution_space execution_space;

  static const unsigned FunctionCount = base_type::FunctionCount;
  static const unsigned IntegrationCount = base_type::IntegrationCount;
  static const unsigned ElemNodeCount = base_type::ElemNodeCount;

  typedef Sacado::Fad::SFad<scalar_type,FunctionCount> fad_scalar_type;

  ElementComputation(const ElementComputation& rhs) : base_type(rhs) {}

  ElementComputation(
    const typename base_type::mesh_type          & arg_mesh ,
    const typename base_type::vector_type        & arg_solution ,
    const typename base_type::elem_graph_type    & arg_elem_graph ,
    const typename base_type::sparse_matrix_type & arg_jacobian ,
    const typename base_type::vector_type        & arg_residual ) :
    base_type(arg_mesh, arg_solution, arg_elem_graph,
              arg_jacobian, arg_residual) {}

  //------------------------------------

  void apply() const
  {
    const size_t nelem = this->elem_node_ids.extent(0);
    parallel_for( nelem , *this );
  }

  KOKKOS_INLINE_FUNCTION
  void gatherSolution(const unsigned ielem,
                      fad_scalar_type val[],
                      unsigned node_index[],
                      double x[], double y[], double z[],
                      fad_scalar_type res[]) const
  {
    for ( unsigned i = 0 ; i < ElemNodeCount ; ++i ) {
      const unsigned ni = this->elem_node_ids( ielem , i );

      node_index[i] = ni ;

      x[i] = this->node_coords( ni , 0 );
      y[i] = this->node_coords( ni , 1 );
      z[i] = this->node_coords( ni , 2 );

      val[i].val() = this->solution( ni );
      val[i].diff( i, FunctionCount );
    }
  }

  KOKKOS_INLINE_FUNCTION
  void scatterResidual(const unsigned ielem,
                       const unsigned node_index[],
                       fad_scalar_type res[]) const
  {
    for( unsigned i = 0 ; i < FunctionCount ; i++ ) {
      const unsigned row = node_index[i] ;
      if ( row < this->residual.extent(0) ) {
        atomic_add( & this->residual( row ) , res[i].val() );

        for( unsigned j = 0 ; j < FunctionCount ; j++ ) {
          const unsigned entry = this->elem_graph( ielem , i , j );
          if ( entry != ~0u ) {
            atomic_add( & this->jacobian.coeff( entry ) ,
                        res[i].fastAccessDx(j) );
          }
        }
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void computeElementResidual(const fad_scalar_type dof_values[] ,
                              const double x[],
                              const double y[],
                              const double z[],
                              fad_scalar_type elem_res[] ) const
  {
    double coeff_k = 3.456;
    double coeff_src = 1.234;
    double advection[] = { 1.1, 1.2, 1.3 };
    double dpsidx[ FunctionCount ] ;
    double dpsidy[ FunctionCount ] ;
    double dpsidz[ FunctionCount ] ;
    for ( unsigned i = 0 ; i < IntegrationCount ; ++i ) {

      const double integ_weight = this->elem_data.weights[i];
      const double* bases_vals = this->elem_data.values[i];
      const double detJ =
        this->transform_gradients( this->elem_data.gradients[i] ,
                                   x , y , z ,
                                   dpsidx , dpsidy , dpsidz );
      const double detJ_weight = detJ * integ_weight;
      const double detJ_weight_coeff_k = detJ_weight * coeff_k;

      fad_scalar_type value_at_pt = 0 ;
      fad_scalar_type gradx_at_pt = 0 ;
      fad_scalar_type grady_at_pt = 0 ;
      fad_scalar_type gradz_at_pt = 0 ;
      for ( unsigned m = 0 ; m < FunctionCount ; m++ ) {
        value_at_pt += dof_values[m] * bases_vals[m] ;
        gradx_at_pt += dof_values[m] * dpsidx[m] ;
        grady_at_pt += dof_values[m] * dpsidy[m] ;
        gradz_at_pt += dof_values[m] * dpsidz[m] ;
      }

      const fad_scalar_type source_term =
        coeff_src * value_at_pt * value_at_pt ;

      const fad_scalar_type advection_term =
        advection[0]*gradx_at_pt +
        advection[1]*grady_at_pt +
        advection[2]*gradz_at_pt;

      for ( unsigned m = 0; m < FunctionCount; ++m) {
        const double bases_val_m = bases_vals[m] * detJ_weight ;
        const double dpsidx_m    = dpsidx[m] ;
        const double dpsidy_m    = dpsidy[m] ;
        const double dpsidz_m    = dpsidz[m] ;

        elem_res[m] +=
          detJ_weight_coeff_k * ( dpsidx_m * gradx_at_pt +
                                  dpsidy_m * grady_at_pt +
                                  dpsidz_m * gradz_at_pt ) +
          bases_val_m * ( advection_term  + source_term ) ;
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( const unsigned ielem ) const
  {
    double x[ FunctionCount ] ;
    double y[ FunctionCount ] ;
    double z[ FunctionCount ] ;
    unsigned node_index[ ElemNodeCount ];

    fad_scalar_type val[ FunctionCount ] ;
    fad_scalar_type elem_res[ FunctionCount ] ; // this zeros elem_res

    // Gather nodal coordinates and solution vector:
    gatherSolution( ielem, val, node_index, x, y, z, elem_res );

    // Compute nodal element residual vector:
    computeElementResidual( val, x, y, z, elem_res );

    // Scatter nodal element residual in global vector:
    scatterResidual( ielem, node_index, elem_res );
  }
}; /* ElementComputation */

template< class ExecutionSpace , BoxElemPart::ElemOrder Order ,
          class CoordinateMap , typename ScalarType >
class ElementComputation
  < Kokkos::Example::BoxElemFixture< ExecutionSpace , Order , CoordinateMap > ,
    CrsMatrix< ScalarType , ExecutionSpace > ,
    FadElementOptimized > :
    public ElementComputation< Kokkos::Example::BoxElemFixture< ExecutionSpace , Order , CoordinateMap > ,
                               CrsMatrix< ScalarType , ExecutionSpace > ,
                               FadElement > {
public:

  typedef ElementComputation< Kokkos::Example::BoxElemFixture< ExecutionSpace , Order , CoordinateMap > ,
                              CrsMatrix< ScalarType , ExecutionSpace > ,
                              FadElement > base_type;

  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::execution_space execution_space;

  static const unsigned FunctionCount = base_type::FunctionCount;
  static const unsigned IntegrationCount = base_type::IntegrationCount;
  static const unsigned ElemNodeCount = base_type::ElemNodeCount;

  typedef Sacado::Fad::SFad<scalar_type,FunctionCount> fad_scalar_type;

  ElementComputation(const ElementComputation& rhs) : base_type(rhs) {}

  ElementComputation(
    const typename base_type::mesh_type          & arg_mesh ,
    const typename base_type::vector_type        & arg_solution ,
    const typename base_type::elem_graph_type    & arg_elem_graph ,
    const typename base_type::sparse_matrix_type & arg_jacobian ,
    const typename base_type::vector_type        & arg_residual ) :
    base_type(arg_mesh, arg_solution, arg_elem_graph,
              arg_jacobian, arg_residual) {}

  //------------------------------------

  void apply() const
  {
    const size_t nelem = this->elem_node_ids.extent(0);
    parallel_for( nelem , *this );
  }

  KOKKOS_INLINE_FUNCTION
  void gatherSolution(const unsigned ielem,
                      scalar_type val[],
                      unsigned node_index[],
                      double x[], double y[], double z[],
                      fad_scalar_type res[]) const
  {
    for ( unsigned i = 0 ; i < ElemNodeCount ; ++i ) {
      const unsigned ni = this->elem_node_ids( ielem , i );

      node_index[i] = ni ;

      x[i] = this->node_coords( ni , 0 );
      y[i] = this->node_coords( ni , 1 );
      z[i] = this->node_coords( ni , 2 );

      val[i] = this->solution( ni );
    }
  }

  KOKKOS_INLINE_FUNCTION
  void computeElementResidual(const scalar_type dof_values[] ,
                              const double x[],
                              const double y[],
                              const double z[],
                              fad_scalar_type elem_res[] ) const
  {
    double coeff_k = 3.456;
    double coeff_src = 1.234;
    double advection[] = { 1.1, 1.2, 1.3 };
    double dpsidx[ FunctionCount ] ;
    double dpsidy[ FunctionCount ] ;
    double dpsidz[ FunctionCount ] ;
    for ( unsigned i = 0 ; i < IntegrationCount ; ++i ) {

      const double integ_weight = this->elem_data.weights[i];
      const double* bases_vals = this->elem_data.values[i];
      const double detJ =
        this->transform_gradients( this->elem_data.gradients[i] ,
                                   x , y , z ,
                                   dpsidx , dpsidy , dpsidz );
      const double detJ_weight = detJ * integ_weight;
      const double detJ_weight_coeff_k = detJ_weight * coeff_k;

      fad_scalar_type value_at_pt(FunctionCount, 0.0, Sacado::NoInitDerivArray) ;
      fad_scalar_type gradx_at_pt(FunctionCount, 0.0, Sacado::NoInitDerivArray) ;
      fad_scalar_type grady_at_pt(FunctionCount, 0.0, Sacado::NoInitDerivArray) ;
      fad_scalar_type gradz_at_pt(FunctionCount, 0.0, Sacado::NoInitDerivArray) ;
      for ( unsigned m = 0 ; m < FunctionCount ; m++ ) {
        value_at_pt.val() += dof_values[m] * bases_vals[m] ;
        value_at_pt.fastAccessDx(m) = bases_vals[m] ;

        gradx_at_pt.val() += dof_values[m] * dpsidx[m] ;
        gradx_at_pt.fastAccessDx(m) = dpsidx[m] ;

        grady_at_pt.val() += dof_values[m] * dpsidy[m] ;
        grady_at_pt.fastAccessDx(m) = dpsidy[m] ;

        gradz_at_pt.val() += dof_values[m] * dpsidz[m] ;
        gradz_at_pt.fastAccessDx(m) = dpsidz[m] ;
      }

      const fad_scalar_type source_term =
        coeff_src * value_at_pt * value_at_pt ;

      const fad_scalar_type advection_term =
        advection[0]*gradx_at_pt +
        advection[1]*grady_at_pt +
        advection[2]*gradz_at_pt;

      for ( unsigned m = 0; m < FunctionCount; ++m) {
        const double bases_val_m = bases_vals[m] * detJ_weight ;
        const double dpsidx_m    = dpsidx[m] ;
        const double dpsidy_m    = dpsidy[m] ;
        const double dpsidz_m    = dpsidz[m] ;

        elem_res[m] +=
          detJ_weight_coeff_k * ( dpsidx_m * gradx_at_pt +
                                  dpsidy_m * grady_at_pt +
                                  dpsidz_m * gradz_at_pt ) +
          bases_val_m * ( advection_term  + source_term )  ;
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( const unsigned ielem ) const
  {
    double x[ FunctionCount ] ;
    double y[ FunctionCount ] ;
    double z[ FunctionCount ] ;
    unsigned node_index[ ElemNodeCount ];

    scalar_type val[ FunctionCount ] ;
    fad_scalar_type elem_res[ FunctionCount ] ;

    // Gather nodal coordinates and solution vector:
    gatherSolution( ielem, val, node_index, x, y, z, elem_res );

    // Compute nodal element residual vector:
    computeElementResidual( val, x, y, z, elem_res );

    // Scatter nodal element residual in global vector:
    this->scatterResidual( ielem, node_index, elem_res );
  }
}; /* ElementComputation */

template< class ExecutionSpace , BoxElemPart::ElemOrder Order ,
          class CoordinateMap , typename ScalarType >
class ElementComputation
  < Kokkos::Example::BoxElemFixture< ExecutionSpace , Order , CoordinateMap > ,
    CrsMatrix< ScalarType , ExecutionSpace > ,
    FadQuadPoint > :
    public ElementComputation< Kokkos::Example::BoxElemFixture< ExecutionSpace , Order , CoordinateMap > ,
                               CrsMatrix< ScalarType , ExecutionSpace > ,
                               Analytic > {
public:

  typedef ElementComputation< Kokkos::Example::BoxElemFixture< ExecutionSpace , Order , CoordinateMap > ,
                              CrsMatrix< ScalarType , ExecutionSpace > ,
                              Analytic > base_type;

  typedef typename base_type::scalar_type scalar_type;
  typedef typename base_type::execution_space execution_space;

  static const unsigned FunctionCount = base_type::FunctionCount;
  static const unsigned IntegrationCount = base_type::IntegrationCount;
  static const unsigned ElemNodeCount = base_type::ElemNodeCount;

  typedef Sacado::Fad::SFad<scalar_type,4> fad_scalar_type;

  ElementComputation(const ElementComputation& rhs) : base_type(rhs) {}

  ElementComputation(
    const typename base_type::mesh_type          & arg_mesh ,
    const typename base_type::vector_type        & arg_solution ,
    const typename base_type::elem_graph_type    & arg_elem_graph ,
    const typename base_type::sparse_matrix_type & arg_jacobian ,
    const typename base_type::vector_type        & arg_residual ) :
    base_type(arg_mesh, arg_solution, arg_elem_graph,
              arg_jacobian, arg_residual) {}

  //------------------------------------

  void apply() const
  {
    const size_t nelem = this->elem_node_ids.extent(0);
    parallel_for( nelem , *this );
  }

  KOKKOS_INLINE_FUNCTION
  void computeElementResidualJacobian(
    const scalar_type dof_values[] ,
    const double x[],
    const double y[],
    const double z[],
    scalar_type  elem_res[] ,
    scalar_type  elem_mat[][FunctionCount] ) const
  {
    double coeff_k = 3.456;
    double coeff_src = 1.234;
    double advection[] = { 1.1, 1.2, 1.3 };
    double dpsidx[ FunctionCount ] ;
    double dpsidy[ FunctionCount ] ;
    double dpsidz[ FunctionCount ] ;

    fad_scalar_type value_at_pt(4, 0, 0.0) ;
    fad_scalar_type gradx_at_pt(4, 1, 0.0) ;
    fad_scalar_type grady_at_pt(4, 2, 0.0) ;
    fad_scalar_type gradz_at_pt(4, 3, 0.0) ;
    for ( unsigned i = 0 ; i < IntegrationCount ; ++i ) {

      const double integ_weight = this->elem_data.weights[i];
      const double* bases_vals = this->elem_data.values[i];
      const double detJ =
        this->transform_gradients( this->elem_data.gradients[i] ,
                                   x , y , z ,
                                   dpsidx , dpsidy , dpsidz );
      const double detJ_weight = detJ * integ_weight;
      const double detJ_weight_coeff_k = detJ_weight * coeff_k;

      value_at_pt.val() = 0.0 ;
      gradx_at_pt.val() = 0.0 ;
      grady_at_pt.val() = 0.0 ;
      gradz_at_pt.val() = 0.0 ;
      for ( unsigned m = 0 ; m < FunctionCount ; m++ ) {
        value_at_pt.val() += dof_values[m] * bases_vals[m] ;
        gradx_at_pt.val() += dof_values[m] * dpsidx[m] ;
        grady_at_pt.val() += dof_values[m] * dpsidy[m] ;
        gradz_at_pt.val() += dof_values[m] * dpsidz[m] ;
      }

      const fad_scalar_type source_term =
        coeff_src * value_at_pt * value_at_pt ;

      const fad_scalar_type advection_term =
        advection[0]*gradx_at_pt +
        advection[1]*grady_at_pt +
        advection[2]*gradz_at_pt;

      for ( unsigned m = 0; m < FunctionCount; ++m) {
        const double bases_val_m = bases_vals[m] * detJ_weight ;
        fad_scalar_type res =
          detJ_weight_coeff_k * ( dpsidx[m] * gradx_at_pt +
                                  dpsidy[m] * grady_at_pt +
                                  dpsidz[m] * gradz_at_pt ) +
          bases_val_m * ( advection_term  + source_term )  ;

        elem_res[m] += res.val();

        scalar_type * const mat = elem_mat[m] ;
        for( unsigned n = 0; n < FunctionCount; n++) {
          mat[n] += res.fastAccessDx(0) * bases_vals[n] +
                    res.fastAccessDx(1) * dpsidx[n] +
                    res.fastAccessDx(2) * dpsidy[n] +
                    res.fastAccessDx(3) * dpsidz[n];
        }
      }
    }
  }

  KOKKOS_INLINE_FUNCTION
  void operator()( const unsigned ielem ) const
  {
    double x[ FunctionCount ] ;
    double y[ FunctionCount ] ;
    double z[ FunctionCount ] ;
    unsigned node_index[ ElemNodeCount ];

    scalar_type val[ FunctionCount ] ;
    scalar_type elem_res[ FunctionCount ] ;
    scalar_type elem_mat[ FunctionCount ][ FunctionCount ] ;

    // Gather nodal coordinates and solution vector:
    this->gatherSolution( ielem, val, node_index, x, y, z, elem_res, elem_mat );

    // Compute nodal element residual vector and Jacobian matrix:
    computeElementResidualJacobian( val, x, y, z, elem_res, elem_mat );

    // Scatter nodal element residual and Jacobian in global vector and matrix:
    this->scatterResidual( ielem, node_index, elem_res, elem_mat );
  }
}; /* ElementComputation */

} /* namespace FENL */
} /* namespace Example */
} /* namespace Kokkos  */

//----------------------------------------------------------------------------

#endif /* #ifndef KOKKOS_EXAMPLE_FENLFUNCTORS_HPP */
