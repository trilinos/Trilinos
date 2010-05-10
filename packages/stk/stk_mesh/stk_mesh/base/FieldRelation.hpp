/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#ifndef stk_mesh_FieldRelation_hpp
#define stk_mesh_FieldRelation_hpp

namespace stk {
namespace mesh {

class FieldBase;

//----------------------------------------------------------------------
/** \ingroup stk_mesh_relation_stencil
 *  \brief  A defined entity-relationship between a field of a pointer type
 *          and the field that it should point to.
 *          An internal class that should never need to be
 *          <em> directly </em> used within application code.
 *
 *  <b> Let </b>
 *  - <b> rel </b> be the \ref stk::mesh::Relation "relation" from
 *    \ref stk::mesh::Entity "entity" <b>e1</b> to
 *    \ref stk::mesh::Entity "entity" <b>e2</b>
 *
 *  <b> If </b>
 *  - \ref stk::mesh::Field "field" <b> m_root </b>
 *    has a pointer scalar type 'T *' AND
 *  - \ref stk::mesh::Field "field" <b> m_target </b>
 *    has a scalar type 'T' AND
 *  - \ref stk_mesh_field_data "field_data"( *m_root , e1 ) exits AND
 *  - \ref stk_mesh_field_data "field_data"( *m_target , e2 ) exits AND
 *  - \ref stk::mesh::Relation "relation" <b> rel </b> is in the domain of
 *    \ref stk_mesh_relation_stencil "relation stencil" <b> m_function </b>
 *
 *  <b> then </b>
 *  <PRE>
 *    index = (*m_function)( e1.entity_rank() ,
 *                           e2.entity_rank() ,
 *                           rel.identifier() ,
 *                           rel.kind() );
 *
 *    field_data(*m_root,e1)[index] == field_data(*m_target,e2)
 *  </PRE>
 */
struct FieldRelation {
  /** \brief  relation domain part */
  FieldBase          * m_root ;

  /** \brief  relation range part */
  FieldBase          * m_target ;

  /** \brief  \ref stk_mesh_relation_stencil "relation stencil" */
  relation_stencil_ptr m_function ;

#ifndef DOXYGEN_COMPILE

  FieldRelation() : m_root( NULL ), m_target( NULL ), m_function( NULL ) {}

  FieldRelation( const FieldRelation & rhs )
    : m_root( rhs.m_root ),
      m_target( rhs.m_target ),
      m_function( rhs.m_function ) {}

  FieldRelation & operator = ( const FieldRelation & rhs )
    {
      m_root = rhs.m_root ;
      m_target = rhs.m_target ;
      m_function = rhs.m_function ;
      return *this ;
    }

#endif /* DOXYGEN_COMPILE */
};

} // namespace mesh
} // namespace stk


#endif //stk_mesh_FieldRelation_hpp
