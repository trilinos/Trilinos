// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#ifndef STK_MESH_RELATION_HPP
#define STK_MESH_RELATION_HPP

#include <stdint.h>                     // for uint32_t
#include <iosfwd>                       // for ostream
#include <stk_mesh/base/Entity.hpp>     // for Entity
#include <stk_mesh/base/Types.hpp>      // for RelationType, EntityRank, etc
#include <stk_mesh/baseImpl/MeshImplUtils.hpp>
#include <stk_util/util/ReportHandler.hpp>  // for ThrowAssertMsg
#include <vector>                       // for vector
namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class Part; } }

namespace stk {
namespace mesh {

/** \addtogroup stk_mesh_module
 *  \{
 */
//----------------------------------------------------------------------
/** \brief  A relation between two mesh entities with a
 *          relation <b> identifier </b> and <b> kind </b>.
 *
 *  Each entity owns a collection of relations to other entities.
 *  Each of these relations has a referenced entity, direction,
 *  relation identifier, and kind.
 *  - The direction specifies whether the relation is from the
 *    entity owning the relation to the referenced entity or
 *    conversely from the referenced entity to the owning entity.
 *  - The identifier provides a local numbering convention for relations,
 *    for example the local numbering convention for the nodes of an element.
 *  - Relations can be given a <em> kind </em> raw_relation_id to differentiate
 *    between topological relations (kind = 0) and other kinds of
 *    relationships such as constraints.
 *
 *  Relations are ordered by their
 *  - kind of relation,
 *  - type of the referenced entity,
 *  - if the relation is forward or converse,
 *  - local relation identifier, and then
 *  - referenced entity's identifier.
 */
class Relation {
public:
  typedef uint32_t raw_relation_id_type ;
  typedef uint32_t attribute_type;

  typedef stk::mesh::RelationType RelationType;

  /** \brief  Constructor */
  Relation();

  /** \brief  Construct a relation from a referenced entity and local identifier
   */
  Relation( Entity entity, EntityRank entityRank , RelationIdentifier identifier );
  // Defined in BulkData.hpp to break circular dependency.

  attribute_type   attribute() const { return m_attribute; }
  void set_attribute(attribute_type attr) const  { m_attribute = attr; }

  /** \brief  The encoded relation raw_relation_id */
  inline static raw_relation_id_type raw_relation_id( EntityRank rank , unsigned id );

  /** \brief  The encoded relation raw_relation_id */
  raw_relation_id_type raw_relation_id() const { return m_raw_relation.value ; }

  /** \brief  The rank of the referenced entity */
  inline EntityRank entity_rank() const ;

  /** \brief  The local relation identifier */
  inline RelationIdentifier relation_ordinal() const ;

  /** \brief  The referenced entity */
  inline Entity entity() const;

  /** \brief  Equality operator */
  inline bool operator == ( const Relation & r ) const;

  /** \brief  Inequality operator */
  bool operator != ( const Relation & r ) const
  { return !(*this == r); }

private:

  enum {
    rank_digits = 8  ,
    id_digits   = 24 ,
    id_mask     = ~(0u) >> rank_digits ,
    fwmk_relation_type_digits = 8,
    fmwk_permutation_digits   = 24,
    fmwk_permutation_mask     = ~(0u) >> fwmk_relation_type_digits
  };

  union RawRelationType {
  public:
    raw_relation_id_type value ;

    struct {
      raw_relation_id_type identifier  : id_digits ;
      raw_relation_id_type entity_rank : rank_digits ;
    } normal_view ;

    struct {
      raw_relation_id_type entity_rank : rank_digits ;
      raw_relation_id_type identifier  : id_digits ;
    } reverse_view ;

    RawRelationType( raw_relation_id_type v ) : value(v) {}
    RawRelationType() : value(0) {}
    RawRelationType( const RawRelationType & rhs ) : value( rhs.value ) {}
    RawRelationType & operator = ( const RawRelationType & rhs );
  };

  RawRelationType        m_raw_relation ;
  mutable attribute_type m_attribute ;

  Entity                 m_target_entity;

 public:
  static raw_relation_id_type max_id() { return (1 << id_digits) - 1;}


// Issue: Framework supports relation types (parent, child, etc) that STK_Mesh
// does not support, so these relations will have to be managed by framework
// until:
// A) STK_Mesh is rewritten to support these extra relation types
// B) An extra data-structure is added to manage these relation types and
// this data structure works together with STK_mesh under the generic API
// to manage all relation-types.
// C) We transition to using a mesh that supports everything framework
// supports and this problem goes away.
//
// The problem is that, with framework managing some relations and STK_Mesh
// managing others, the type of the relation descriptor is different depending
// on what type of relation you're dealing with. This can be addressed with
// templates, but this makes the code very ugly. Instead...
//
// Solution: Have framework and STK_Mesh use the same type as its relation_descriptor.
// The code below is designed to make this class compatible with the fmwk
// Relation class.
 public:

  // Moved this to enum in struct RelationType.
  //   /**
  //   * Predefined identifiers for mesh object relationship types.
  //   */
  //  enum RelationType {
  //    USES	= 0 ,
  //    USED_BY	= 1 ,
  //    INVALID     = 10
  //  };

  enum {
    POLARITY_MASK       = 0x80,
    POLARITY_POSITIVE   = 0x80,
    POLARITY_NEGATIVE   = 0x00,
    POLARITY_IDENTITY   = 0x80
  };

  /**
   * Construct filled-out relation, fmwk-style
   */
  // Only needed by Framework and Framework-based apps.
  Relation(EntityRank, Entity obj, const unsigned relation_type, const unsigned ordinal, const unsigned permut = 0);

  // Only needed by Framework and Framework-based apps.
  inline void setMeshObj(Entity object, EntityRank object_rank);

  RelationType  getRelationType() const {
    return static_cast<RelationType::relation_type_t >(attribute() >> fmwk_permutation_digits);
  }

  void setRelationType(RelationType relation_type) {
    set_attribute( (relation_type << fmwk_permutation_digits) | getPermutation() );
  }

  ConnectivityOrdinal getOrdinal() const {
    return static_cast<ConnectivityOrdinal>(relation_ordinal());
  }

  void setOrdinal(RelationIdentifier ordinal) {
    m_raw_relation = Relation::raw_relation_id( entity_rank(), ordinal );
  }

  Permutation getPermutation() const {
    return static_cast<Permutation>(attribute() & fmwk_permutation_mask);
  }

  void setPermutation(unsigned perm) {
    set_attribute( (getRelationType() << fmwk_permutation_digits) | perm );
  }

  /**
   * Query polarity of the related mesh object.
   * A 'true' polarity indicates that the related mesh object is aligned.
   * For element-edge or face-edge a 'true' polarity indicates that the
   * nodes of the edge are compatibly ordered with the nodes of the
   * element or faces.
   * For element-face a 'true' polarity indicates that the ordering
   * of face-nodes is compatible with the ordering of element nodes,
   * i.e. the face's normal defined by a clockwise ordering is outward.
   */
  static bool compute_polarity(const stk::topology & topology, unsigned permutation)
  {
    STK_ThrowAssert(permutation < topology.num_permutations());
    const bool polarity = permutation < topology.num_positive_permutations();
    return polarity;
  }


private:
  bool has_fmwk_state() const { return getRelationType() != RelationType::INVALID; }
};

//----------------------------------------------------------------------

inline
Relation::raw_relation_id_type
Relation::raw_relation_id( EntityRank rank , unsigned id )
{
  STK_ThrowAssertMsg( id <= id_mask,
                  "For args rank " << rank << ", id " << id << ": " <<
                  "id " << " > id_mask=" << id_mask );

  return ( static_cast<raw_relation_id_type>(rank) << id_digits ) | id ;
}

inline
EntityRank Relation::entity_rank() const
{ return static_cast<EntityRank>(m_raw_relation.value >> id_digits); }

inline
RelationIdentifier Relation::relation_ordinal() const
{ return static_cast<unsigned>( m_raw_relation.value & id_mask ); }

/** \brief  Fill a vector with mesh entities of the given rank,
 *          that are related (connected) to all of
 *          the input mesh entities.
 */
template<typename Vec1Type, typename Vec2Type>
void get_entities_through_relations(
               const BulkData& mesh,
               const Vec1Type& entities,
                     EntityRank relatedEntitiesRank,
                     Vec2Type& relatedEntities)
{
  impl::find_entities_these_nodes_have_in_common(mesh, relatedEntitiesRank,
                                                 entities.size(), entities.data(),
                                                 relatedEntities);
}

/** \brief  Induce an entity's part membership based upon relationships
 *          from other entities.  Do not include and parts in the 'omit' list.
 */
void induced_part_membership(const BulkData& mesh,
                             const Entity entity ,
                             OrdinalVector & induced_parts);

//----------------------------------------------------------------------

/** \} */

inline
Relation::Relation() :
  m_raw_relation(),
  m_attribute(),
  m_target_entity()
{
  setRelationType(RelationType::INVALID);
}

inline
Relation::Relation( Entity ent, EntityRank entityRank , RelationIdentifier id )
  : m_raw_relation( Relation::raw_relation_id( entityRank , id ) ),
    m_attribute(),
    m_target_entity(ent)
{
  setRelationType(RelationType::INVALID);
}

inline
bool Relation::operator == ( const Relation & rhs ) const
{
  return m_raw_relation.value == rhs.m_raw_relation.value && m_target_entity == rhs.m_target_entity
    // compared fmwk state too
    && m_attribute == rhs.m_attribute
    ;
}

inline
bool same_specification(const Relation& lhs, const Relation& rhs)
{
  return  lhs.entity_rank()     == rhs.entity_rank() &&
          lhs.getRelationType() == rhs.getRelationType() &&
          lhs.getOrdinal()      == rhs.getOrdinal();
}

inline
RelationType
back_relation_type(const RelationType relType)
{
  switch(relType) {
  case RelationType::USES:
    return RelationType::USED_BY;
  case RelationType::USED_BY:
    return RelationType::USES;
  default:
    return relType;
  }
}

// Made publicly available so that MeshObj.C can use it
template <class Iterator>
bool
verify_relation_ordering(Iterator begin, Iterator end)
{
  for (Iterator itr = begin; itr != end; ) {
    Iterator prev = itr;
    ++itr;
    if (itr != end) {

      if (itr->entity_rank() < prev->entity_rank()) {
        return false ;
      }

      if (itr->entity_rank() == prev->entity_rank()) {

        if (itr->getRelationType() < prev->getRelationType()) {
          return false ;
        }

        if (itr->getRelationType() == prev->getRelationType()) {

          if (itr->getOrdinal() < prev->getOrdinal()) {
            return false ;
          }
        }
      }
    }
  }
  return true ;
}

namespace impl {

inline
bool internal_is_handled_generically(const RelationType relation_type)
{ return relation_type == RelationType::USES || relation_type == RelationType::USED_BY; }

}

inline
Entity Relation::entity() const
{
  return m_target_entity;
}

inline
Relation::Relation(EntityRank rel_rank, Entity obj, const unsigned relation_type, const unsigned ordinal, const unsigned permut)
  :
      m_raw_relation( Relation::raw_relation_id(rel_rank, ordinal )),
      m_attribute( (relation_type << fmwk_permutation_digits) | permut ),
      m_target_entity(obj)
{
  STK_ThrowAssertMsg( permut <= fmwk_permutation_mask,
      "permutation " << permut << " exceeds maximum allowed value");
}

inline
void Relation::setMeshObj(Entity object, EntityRank object_rank )
{
  m_raw_relation = Relation::raw_relation_id( object_rank, relation_ordinal() );
  m_target_entity = object;
}

// Find relation from mesh_obj to an Entity of certain rank by relation ordinal
stk::mesh::Entity find_by_ordinal(stk::mesh::Entity mesh_obj, stk::mesh::EntityRank rank, size_t ordinal, const stk::mesh::BulkData& mesh);

} // namespace mesh
} // namespace stk
#endif
