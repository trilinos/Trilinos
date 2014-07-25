/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

/**
 * @author H. Carter Edwards
 */

#include <stdexcept>
#include <iostream>

#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>

#include <app/performance_algorithms.hpp>

#define USE_HARDWIRED_NUMBER_NODES 1

namespace stk_classic {
namespace app {

//----------------------------------------------------------------------

#if USE_HARDWIRED_NUMBER_NODES

void element_mean_value_8x3( const unsigned num_elems ,
                             const double * const * node_values ,
                             double * elem_value )
{
  for ( unsigned i = 0 ; i < num_elems ; ++i ) {

    elem_value[0] = ( node_values[0][0] +
                      node_values[1][0] +
                      node_values[2][0] +
                      node_values[3][0] +
                      node_values[4][0] +
                      node_values[5][0] +
                      node_values[6][0] +
                      node_values[7][0] ) / 8 ;

    elem_value[1] = ( node_values[0][1] +
                      node_values[1][1] +
                      node_values[2][1] +
                      node_values[3][1] +
                      node_values[4][1] +
                      node_values[5][1] +
                      node_values[6][1] +
                      node_values[7][1] ) / 8 ;

    elem_value[2] = ( node_values[0][2] +
                      node_values[1][2] +
                      node_values[2][2] +
                      node_values[3][2] +
                      node_values[4][2] +
                      node_values[5][2] +
                      node_values[6][2] +
                      node_values[7][2] ) / 8 ;

    node_values += 8 ;
    elem_value += 3 ;
  }
}

void element_mean_value_4x3( const unsigned num_elems ,
                             const double * const * node_values ,
                             double * elem_value )
{
  for ( unsigned i = 0 ; i < num_elems ; ++i ) {

    elem_value[0] = ( node_values[0][0] +
                      node_values[1][0] +
                      node_values[2][0] +
                      node_values[3][0] ) / 4 ;

    elem_value[1] = ( node_values[0][1] +
                      node_values[1][1] +
                      node_values[2][1] +
                      node_values[3][1] ) / 4 ;

    elem_value[2] = ( node_values[0][2] +
                      node_values[1][2] +
                      node_values[2][2] +
                      node_values[3][2] ) / 4 ;

    node_values += 4 ;
    elem_value += 3 ;
  }
}

enum { GATHER_COUNT_SIZE = 1000 };

void element_gather_mean_value_8x3( const unsigned num_elems ,
                                    const double * const * node_values ,
                                    double * elem_value )
{
  typedef double value_type[3] ;
  value_type gathered_value[ GATHER_COUNT_SIZE * 8 ];

  // Gather:

  {
    value_type * gv = gathered_value ;

    for ( unsigned i = 0 ; i < num_elems ; ++i ) {

      gv[0][0] = node_values[0][0];
      gv[0][1] = node_values[0][1];
      gv[0][2] = node_values[0][2];

      gv[1][0] = node_values[1][0];
      gv[1][1] = node_values[1][1];
      gv[1][2] = node_values[1][2];

      gv[2][0] = node_values[2][0];
      gv[2][1] = node_values[2][1];
      gv[2][2] = node_values[2][2];

      gv[3][0] = node_values[3][0];
      gv[3][1] = node_values[3][1];
      gv[3][2] = node_values[3][2];

      gv[4][0] = node_values[4][0];
      gv[4][1] = node_values[4][1];
      gv[4][2] = node_values[4][2];

      gv[5][0] = node_values[5][0];
      gv[5][1] = node_values[5][1];
      gv[5][2] = node_values[5][2];

      gv[6][0] = node_values[6][0];
      gv[6][1] = node_values[6][1];
      gv[6][2] = node_values[6][2];

      gv[7][0] = node_values[7][0];
      gv[7][1] = node_values[7][1];
      gv[7][2] = node_values[7][2];

      gv += 8 ;
      node_values += 8 ;
    }
  }

  // Compute:

  {
    value_type * gv = gathered_value ;

    for ( unsigned i = 0 ; i < num_elems ; ++i ) {

      elem_value[0] = ( gv[0][0] + gv[1][0] + gv[2][0] + gv[3][0] +
                        gv[4][0] + gv[5][0] + gv[6][0] + gv[7][0] ) / 8 ;

      elem_value[1] = ( gv[0][1] + gv[1][1] + gv[2][1] + gv[3][1] +
                        gv[4][1] + gv[5][1] + gv[6][1] + gv[7][1] ) / 8 ;

      elem_value[2] = ( gv[0][2] + gv[1][2] + gv[2][2] + gv[3][2] +
                        gv[4][2] + gv[5][2] + gv[6][2] + gv[7][2] ) / 8 ;

      gv += 8 ;
      elem_value += 3 ;
    }
  }
}

void element_gather_mean_value_4x3( const unsigned num_elems ,
                                    const double * const * node_values ,
                                    double * elem_value )
{
  typedef double value_type[3] ;
  value_type gathered_value[ GATHER_COUNT_SIZE * 4 ];

  // Gather:

  {
    value_type * gv = gathered_value ;

    for ( unsigned i = 0 ; i < num_elems ; ++i ) {

      gv[0][0] = node_values[0][0];
      gv[0][1] = node_values[0][1];
      gv[0][2] = node_values[0][2];

      gv[1][0] = node_values[1][0];
      gv[1][1] = node_values[1][1];
      gv[1][2] = node_values[1][2];

      gv[2][0] = node_values[2][0];
      gv[2][1] = node_values[2][1];
      gv[2][2] = node_values[2][2];

      gv[3][0] = node_values[3][0];
      gv[3][1] = node_values[3][1];
      gv[3][2] = node_values[3][2];

      gv += 4 ;
      node_values += 4 ;
    }
  }

  // Compute:

  {
    value_type * gv = gathered_value ;

    for ( unsigned i = 0 ; i < num_elems ; ++i ) {

      elem_value[0] = ( gv[0][0] + gv[1][0] + gv[2][0] + gv[3][0] ) / 4 ;
      elem_value[1] = ( gv[0][1] + gv[1][1] + gv[2][1] + gv[3][1] ) / 4 ;
      elem_value[2] = ( gv[0][2] + gv[1][2] + gv[2][2] + gv[3][2] ) / 4 ;

      gv += 4 ;
      elem_value += 3 ;
    }
  }
}

//----------------------------------------------------------------------
#else
//----------------------------------------------------------------------

void element_mean_value_3( const unsigned num_elems ,
                           const unsigned nodes_per_elem ,
                           const double * const * node_values ,
                           double * elem_value )
{
  for ( unsigned i = 0 ; i < num_elems ; ++i ) {

    double tmp[3] = { 0 , 0 , 0 };

    for ( unsigned j = 0 ; j < nodes_per_elem ; ++j ) {
      const double * const value = * node_values ; ++node_values ;
      tmp[0] += value[0] ;
      tmp[1] += value[1] ;
      tmp[2] += value[2] ;
    }

    elem_value[0] = tmp[0] / nodes_per_elem ;
    elem_value[1] = tmp[1] / nodes_per_elem ;
    elem_value[2] = tmp[2] / nodes_per_elem ;

    elem_value += 3 ;
  }
}

//----------------------------------------------------------------------

enum { GATHER_BUFFER_SIZE = 1000 * 8 * 3 };

void element_gather_mean_value_3( const unsigned num_elems ,
                                  const unsigned nodes_per_elem ,
                                  const double * const * node_values ,
                                  double * elem_value )
{
  double gathered_value[ GATHER_BUFFER_SIZE ];

  // Gather:

  for ( unsigned i = 0 ; i < num_elems ; ++i ) {
    const unsigned elem_offset = i * nodes_per_elem ;

    for ( unsigned j = 0 ; j < nodes_per_elem ; ++j ) {

      const unsigned node_offset       = elem_offset + j ;
      const unsigned node_value_offset = node_offset * 3 ;

      gathered_value[node_value_offset+0] = node_values[node_offset][0];
      gathered_value[node_value_offset+1] = node_values[node_offset][1];
      gathered_value[node_value_offset+2] = node_values[node_offset][2];
    }
  }

  // Compute:

  for ( unsigned i = 0 ; i < num_elems ; ++i ) {
    const unsigned elem_offset = i * nodes_per_elem ;

    double tmp[3] = { 0 , 0 , 0 };

    for ( unsigned j = 0 ; j < nodes_per_elem ; ++j ) {
      const unsigned node_offset = elem_offset + j ;
      const double * const value = & gathered_value[ node_offset * 3 ];

      tmp[0] += value[0] ;
      tmp[1] += value[1] ;
      tmp[2] += value[2] ;
    }

    elem_value[0] = tmp[0] / nodes_per_elem ;
    elem_value[1] = tmp[1] / nodes_per_elem ;
    elem_value[2] = tmp[2] / nodes_per_elem ;

    elem_value += 3 ;
  }
}

//----------------------------------------------------------------------
#endif
//----------------------------------------------------------------------

void scale_sum( const unsigned number ,
                const double a ,
                double * const X ,
                const double b ,
                const double * const Y ,
                const double c ,
                const double * const Z )
{
#if 0
  if ( Z ) {
    for ( unsigned i = 0 ; i < number ; ++i ) {
      X[i] = a * X[i] + b * Y[i] + c * Z[i] ;
    }
  }
  else {
    for ( unsigned i = 0 ; i < number ; ++i ) {
      X[i] = a * X[i] + b * Y[i] ;
    }
  }
#else
  double * const x_end = X + number ;
  double * x = X ;
  const double * y = Y ;
  if ( Z ) {
    const double * z = Z ;
    for ( ; x < x_end ; ++x , ++y , ++z ) {
      *x = a * *x + b * *y + c * *z ;
    }
  }
  else {
    for ( ; x < x_end ; ++x , ++y ) {
      *x = a * *x + b * *y ;
    }
  }
#endif
}

void ElementMeanValue::apply( stk_classic::mesh::Bucket::iterator ibegin ,
                              stk_classic::mesh::Bucket::iterator iend ) const
{
  const stk_classic::mesh::BucketArray< stk_classic::mesh::Field<double*,stk_classic::mesh::ElementNode> >
    elem_node( elem_node_field , ibegin , iend );

  const double * const * const node_ptr = elem_node.contiguous_data();

  if ( node_ptr ) {

    const unsigned nodes_per_elem = elem_node.dimension<0>();
    const unsigned number_elem    = elem_node.dimension<1>();

    double * const elem_data = field_data( elem_field , ibegin );

#if USE_HARDWIRED_NUMBER_NODES
    switch( nodes_per_elem ) {
    case 8 :
      element_mean_value_8x3( number_elem, node_ptr, elem_data );
      break ;
    case 4 :
      element_mean_value_4x3( number_elem, node_ptr, elem_data );
      break ;
    default:
      throw std::logic_error( std::string("BAD nodes-per-element" ) );
    }
#else
    element_mean_value_3( number_elem, nodes_per_elem, node_ptr, elem_data );
#endif
  }
}

void ElementMeanValue_Gather::apply( stk_classic::mesh::Bucket::iterator ibegin ,
                                     stk_classic::mesh::Bucket::iterator iend ) const

{
  const stk_classic::mesh::BucketArray< stk_classic::mesh::Field<double*,stk_classic::mesh::ElementNode> >
    elem_node( elem_node_field , ibegin , iend );

  const double * const * const node_ptr = elem_node.contiguous_data();

  if ( node_ptr ) {

    const unsigned nodes_per_elem = elem_node.dimension<0>();
    const unsigned number_elem    = elem_node.dimension<1>();

    double * const elem_data = field_data( elem_field , ibegin );

#if USE_HARDWIRED_NUMBER_NODES
    switch( nodes_per_elem ) {
    case 8 :
      for ( unsigned i = 0 ; i < number_elem ; ) {
        unsigned n = GATHER_COUNT_SIZE < ( number_elem - i ) ?
                     GATHER_COUNT_SIZE : ( number_elem - i );

        element_gather_mean_value_8x3(
          n, node_ptr + i * 8 , elem_data + i * 3 );
        i += n ;
      }
      break ;
    case 4 :
      for ( unsigned i = 0 ; i < number_elem ; ) {
       unsigned n = GATHER_COUNT_SIZE < ( number_elem - i ) ?
                    GATHER_COUNT_SIZE : ( number_elem - i );

        element_gather_mean_value_4x3(
          n, node_ptr + i * 4 , elem_data + i * 3 );
        i += n ;
      }
      break ;
    default:
       throw std::logic_error( std::string("BAD nodes-per-element" ) );
    }
#else
    const unsigned num_elem_workset =
      GATHER_BUFFER_SIZE / ( nodes_per_elem * 3 );

    for ( unsigned i = 0 ; i < number_elem ; ) {
      unsigned n = num_elem_workset < ( number_elem - i ) ?
                   num_elem_workset : ( number_elem - i );

      element_gather_mean_value_3(
        n, nodes_per_elem, node_ptr + i * nodes_per_elem ,
                           elem_data + i * 3 );
      i += n ;
    }
#endif
  }
}

//----------------------------------------------------------------------

void NodeScaleSum::apply( stk_classic::mesh::Bucket::iterator ibegin ,
                          stk_classic::mesh::Bucket::iterator iend ) const
{
  enum { SpaceDim = 3 };

        double * x_val = field_data( X , ibegin );
  const double * y_val = field_data( Y , ibegin );
  const double * z_val = field_data( Z , ibegin );
  const int size = iend - ibegin ;
  
  scale_sum( size * SpaceDim , a , x_val , b , y_val , c , z_val );
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void verify_elem_node_coord(
  stk_classic::mesh::BulkData & mesh ,
  const stk_classic::mesh::Field<double*,stk_classic::mesh::ElementNode> & elem_node_coord ,
  const stk_classic::mesh::Field<double,stk_classic::mesh::Cartesian>    & node_coord )
{
  typedef stk_classic::mesh::Field<double*,stk_classic::mesh::ElementNode> ElemNodeFieldType ;

  const stk_classic::mesh::EntityRank element_rank = stk_classic::mesh::fem::FEMMetaData::get(mesh).element_rank();

  const std::vector<stk_classic::mesh::Bucket*> & buckets = mesh.buckets( element_rank );

  for ( std::vector<stk_classic::mesh::Bucket*>::const_iterator
        k = buckets.begin() ; k != buckets.end() ; ++k ) {
    stk_classic::mesh::Bucket & bucket = **k ;

    stk_classic::mesh::BucketArray< ElemNodeFieldType > array( elem_node_coord , bucket );

    const unsigned num_node = array.dimension<0>();
    const unsigned size     = array.dimension<1>();

    double * const * elem_data = array.contiguous_data();

    for ( unsigned i = 0 ; i < size ; ++i ) {
      stk_classic::mesh::Entity & elem = bucket[i] ;

      stk_classic::mesh::PairIterRelation rel = elem.relations( stk_classic::mesh::fem::FEMMetaData::NODE_RANK );

      for ( unsigned j = 0 ; j < num_node ; ++j , ++elem_data ) {
        stk_classic::mesh::Entity & node = * rel[j].entity();

        double * const node_data = field_data( node_coord , node );
        if ( *elem_data != node_data ) {
          std::cout << "verify_elem_node_coord ERROR, *elem_data != node_data" << std::endl;
        }
      }
    }
  }
}

}//namespace app
}//namespace stk_classic

