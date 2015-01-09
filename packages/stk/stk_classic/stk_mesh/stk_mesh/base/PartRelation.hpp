#ifndef stk_mesh_PartRelation_hpp
#define stk_mesh_PartRelation_hpp

#include <stk_mesh/base/Types.hpp>

namespace stk_classic {
namespace mesh {

class Part;

/** \addtogroup stk_mesh_module
 *  \{
 */

//----------------------------------------------------------------------
/** \ingroup stk_mesh_relations
 *  \brief  A defined entity-relationship between
 *          \ref stk_classic::mesh::Part "parts".
 *          An internal class that should never need to be
 *          <em> directly </em> used within application code.
 *
 *  <b> If </b> an
 *     \ref stk_classic::mesh::Entity "entity" <b> e1 </b> is a member of
 *     \ref stk_classic::mesh::Part "part" <em> m_root </em> and
 *     there exists a
 *     \ref stk_classic::mesh::Relation "relation"
 *     from entity <b> e1 </b> to entity <b> e2 </b> that
 *     is in the domain of the
 *     \ref stk_mesh_relations "relation stencil"
 *     <em> m_function </em>
 *  <b> then </b> entity <b> e2 </b> is a member of part <em> m_target </em>.
 */
struct PartRelation {
  /** \brief  relation domain part */
  Part * m_root ;

  /** \brief  relation range part */
  Part * m_target ;

  /** \brief  \ref stk_mesh_relations "relation stencil" */
  relation_stencil_ptr m_function ;

#ifndef DOXYGEN_COMPILE

  ~PartRelation() {}

  PartRelation() : m_root( NULL ), m_target( NULL ), m_function( NULL ) {}

  PartRelation( const PartRelation & rhs )
    : m_root( rhs.m_root ),
      m_target( rhs.m_target ),
      m_function( rhs.m_function ) {}

  PartRelation & operator = ( const PartRelation & rhs )
  {
    m_root = rhs.m_root ;
    m_target = rhs.m_target ;
    m_function = rhs.m_function ;
    return *this ;
  }

#endif /* DOXYGEN_COMPILE */

};

/** \} */

} // namespace mesh 
} // namespace stk_classic 

#endif // stk_mesh_PartRelation_hpp
