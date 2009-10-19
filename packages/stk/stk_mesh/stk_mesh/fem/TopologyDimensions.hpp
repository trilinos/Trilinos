/**
 * @author H. Carter Edwards  <hcedwar@sandia.gov>
 * @date   June 2008
 */

#ifndef stk_mesh_TopologyDimensions_hpp
#define stk_mesh_TopologyDimensions_hpp

// #include <stk_util/mdarray/Array.hpp>
#include <Shards_Array.hpp>

#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/fem/Stencils.hpp>

namespace stk {
namespace mesh {


//----------------------------------------------------------------------
/** \ingroup stk_mesh_field_dimension_tags
 *  \brief  Define an array dimension of the number of nodes per element.
 */
class ElementNode : public shards::ArrayDimTag {
public:
  const char * name() const ;
  static const ElementNode & tag(); ///< \brief Singleton
private:
  ElementNode() {}
  ElementNode( const ElementNode & );
  ElementNode & operator = ( const ElementNode & );
};

/** \ingroup stk_mesh_relation_stencil
 *  \brief  An element Field defining an array of values, one value per
 *          node of the element.
 */
typedef Field<double,ElementNode> ElementNodeField ;

/** \ingroup stk_mesh_relation_stencil
 *  \brief  A Field defining an array of pointers
 *          to an element's nodal field data.
 */
typedef Field<double*,ElementNode> ElementNodePointerField ;
typedef Field<int *,ElementNode> ElementNodeLockField ;

/** \ingroup stk_mesh_relation_stencil
 *  \brief  Declare an element-node field.
 */
inline
ElementNodeField &
declare_element_node_field( MetaData & md , const std::string & s )
{

  ElementNodeField & f =
    md.declare_field< ElementNodeField >( s, 1 /* 1 state */ );

  return f ;
}

/** \ingroup stk_mesh_relation_stencil
 *  \brief  Declare an element-to-node-data pointer field.
 */
template< class NodeField >
inline
ElementNodePointerField &
declare_element_node_pointer_field(
  MetaData & md , const std::string & s ,
  NodeField & node_field )
{
  const unsigned num_states = node_field.number_of_states();

  ElementNodePointerField & f =
    md.template declare_field< ElementNodePointerField >( s, num_states );

  for ( unsigned i = 0 ; i < num_states ; ++i ) {
    FieldState state = (FieldState) i;
    md.declare_field_relation(
      f.field_of_state( state ) , & element_node_stencil<void> , node_field.field_of_state( state ) );
  }
  
  return f ;
}

/** \ingroup stk_mesh_relation_stencil
 *  \brief  Declare an element-to-node-data pointer field.
 */
template< class NodeField >
inline
ElementNodeLockField &
declare_element_node_lock_field(
  MetaData & md , const std::string & s ,
  NodeField & node_field )
{
  const unsigned num_states = node_field.number_of_states();

  ElementNodeLockField & f =
    md.template declare_field< ElementNodeLockField >( s, num_states );

  for ( unsigned i = 0 ; i < num_states ; ++i ) {
    FieldState state = (FieldState) i;
    md.declare_field_relation(
      f.field_of_state( state ) , & element_node_stencil<void> , node_field.field_of_state( state ) );
  }
  
  return f ;
}
//----------------------------------------------------------------------
/** \addtogroup stk_mesh_field_dimension_tags
 *  \{
 */

struct QuadratureTag : public shards::ArrayDimTag {
  const char * name() const ;
  static const QuadratureTag & tag(); ///< \brief Singleton
private:
  QuadratureTag() {}
  QuadratureTag( const QuadratureTag & );
  QuadratureTag & operator = ( const QuadratureTag & );
};

struct BasisTag : public shards::ArrayDimTag {
  const char * name() const ;
  static const BasisTag & tag(); ///< \brief Singleton
private:
  BasisTag() {}
  BasisTag( const BasisTag & );
  BasisTag & operator = ( const BasisTag & );
};


/** \} */
//----------------------------------------------------------------------

}//namespace mesh
}//namespace stk

#endif

