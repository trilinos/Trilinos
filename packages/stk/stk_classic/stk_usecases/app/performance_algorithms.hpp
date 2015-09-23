/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

/**
 * @author H. Carter Edwards
 */

#ifndef stk_mesh_use_cases_performance_algorithms_hpp
#define stk_mesh_use_cases_performance_algorithms_hpp

#include <utility>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>

#include <stk_mesh/fem/CoordinateSystems.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>

namespace stk_classic {
namespace app {

//----------------------------------------------------------------------

struct ElementMeanValue {

  enum { maximum_entity_count = 0 }; /* Don't slice buckets */

  void apply( stk_classic::mesh::Bucket::iterator ibegin ,
              stk_classic::mesh::Bucket::iterator iend ) const;

  ElementMeanValue(
                    const stk_classic::mesh::Field<double,stk_classic::mesh::Cartesian>    & arg_elem_field ,
                    const stk_classic::mesh::Field<double*,stk_classic::mesh::ElementNode> & arg_node_field_ptr )
    : elem_field( arg_elem_field ),
      elem_node_field( arg_node_field_ptr ) {}

  const stk_classic::mesh::Field<double,stk_classic::mesh::Cartesian>    & elem_field ;
  const stk_classic::mesh::Field<double*,stk_classic::mesh::ElementNode> & elem_node_field ;
};

struct ElementMeanValue_Gather {

  enum { maximum_entity_count = 0 }; /* Don't slice buckets */

  void apply( stk_classic::mesh::Bucket::iterator ibegin ,
              stk_classic::mesh::Bucket::iterator iend ) const;

  ElementMeanValue_Gather(
                    const stk_classic::mesh::Field<double,stk_classic::mesh::Cartesian>    & arg_elem_field ,
                    const stk_classic::mesh::Field<double*,stk_classic::mesh::ElementNode> & arg_node_field_ptr )
    : elem_field( arg_elem_field ),
      elem_node_field( arg_node_field_ptr ) {}

  const stk_classic::mesh::Field<double,stk_classic::mesh::Cartesian>    & elem_field ;
  const stk_classic::mesh::Field<double*,stk_classic::mesh::ElementNode> & elem_node_field ;
};

//----------------------------------------------------------------------

struct NodeScaleSum {

  enum { maximum_entity_count = 0 }; /* Don't slice buckets. */
  void apply( stk_classic::mesh::Bucket::iterator ibegin ,
              stk_classic::mesh::Bucket::iterator iend ) const;

  NodeScaleSum(
    double arg_a , const stk_classic::mesh::Field<double,stk_classic::mesh::Cartesian> & arg_X ,
    double arg_b , const stk_classic::mesh::Field<double,stk_classic::mesh::Cartesian> & arg_Y ,
    double arg_c , const stk_classic::mesh::Field<double,stk_classic::mesh::Cartesian> & arg_Z )
  : a( arg_a ), b( arg_b ), c( arg_c ),
    X( arg_X ), Y( arg_Y ), Z( arg_Z ) {}

  const double a ;
  const double b ;
  const double c ;
  const stk_classic::mesh::Field<double,stk_classic::mesh::Cartesian> & X ;
  const stk_classic::mesh::Field<double,stk_classic::mesh::Cartesian> & Y ;
  const stk_classic::mesh::Field<double,stk_classic::mesh::Cartesian> & Z ;
};

//----------------------------------------------------------------------

void verify_elem_node_coord(
  stk_classic::mesh::BulkData & mesh ,
  const stk_classic::mesh::Field<double*,stk_classic::mesh::ElementNode> & elem_node_coord ,
  const stk_classic::mesh::Field<double,stk_classic::mesh::Cartesian>    & node_coord );

}//namespace app
}//namespace stk_classic

#endif

