/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef Stk_Mesh_Use_Cases_UseCase_Common_hpp
#define Stk_Mesh_Use_Cases_UseCase_Common_hpp

#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>

#include <stk_mesh/fem/CoordinateSystems.hpp>

namespace stk {
namespace mesh {

class Entity;
class Part;
class ElementNode;

namespace use_cases {

typedef Field<double,Cartesian>    VectorFieldType ;
typedef Field<double>              ScalarFieldType ;
typedef Field<double*,ElementNode> ElementNodePointerFieldType ;

/**
 * For elem with element node pointer field, verify that the
 * element node pointer field correctly refers to the vector field
 * of its nodes.
 */
bool verify_elem_node_coord(
  Entity & elem ,
  const ElementNodePointerFieldType & elem_node_coord ,
  const VectorFieldType & node_coord ,
  const unsigned node_count );

/**
 * For every element in part and bucket vector verify
 * elem node coord.
 */
bool verify_elem_node_coord_by_part(
    Part & part,
    const std::vector<Bucket *> & bucket_vector,
    const ElementNodePointerFieldType & elem_node_coord,
    const VectorFieldType & node_coord,
    const unsigned node_count );


} //namespace use_cases
} //namespace mesh
} //namespace stk

#endif // Stk_Mesh_Use_Cases_UseCase_Common_hpp
