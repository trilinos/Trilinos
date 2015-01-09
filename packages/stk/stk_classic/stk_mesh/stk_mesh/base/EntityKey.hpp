/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_mesh_EntityKey_hpp
#define stk_mesh_EntityKey_hpp

#include <stdint.h>
#include <limits>
#include <boost/functional/hash.hpp>

#include <stk_mesh/base/Types.hpp>

namespace stk_classic {
namespace mesh {


// Note:  EntityRank and EntityId typedefs are defined in Types.hpp

/** \addtogroup stk_mesh_module
 * \{
 */

//----------------------------------------------------------------------
/** \brief  Integer type for the entity keys, which is an encoding
 *          of the entity type and entity identifier.
 *
 *  This type is used to fully order entities and entity pointers
 *  in numerous containers.  Ordering is first by type and second
 *  by identifier.  Values of this type are frequently compared for
 *  ordering and equality, and are frequently communicated between
 *  parallel processes.  Thus the motivation for this type to be
 *  a "plain old data" type.
 *
 * Notes on construction and validity:
 *
 * EntityKey takes constructor arguments entity_rank and entity_id which
 * are restricted to lie in a sub-range of what can be represented by the
 * EntityRank and EntityId types. The sub-range allowed by EntityKey is
 * dictated by the amount of space used to store the internal encoding of
 * those values (described further in the comments for the EntityKey
 * constructor below).
 *
 * The function entity_key_valid attempts to determine that the key was
 * not created by some erroneous operation such as assignment from a
 * smaller type like a 32-bit int.
 *
 * Note that an instance of stk_classic::mesh may place further restrictions on a
 * 'valid' key, such as requiring that
 *         0 <= entity_rank(key) < meta_data.entity_rank_count().
 *
 * Typically stk_classic::mesh does not take EntityKeys as input. A user of
 * stk_classic::mesh would (for instance) request that an Entity be created by
 * specifying an entity-type and entity-id as input. The resulting
 * Entity would then hold a mesh-created EntityKey that could be queried.
 * Thus stk_classic::mesh can control the validity of EntityKeys associated with
 * the mesh.
 */
union EntityKey {
public:
  typedef uint64_t raw_key_type ;

  enum { rank_digits = 8 };

private:

  enum {
    invalid_key = ~raw_key_type(0) ,
    raw_digits  = std::numeric_limits<raw_key_type>::digits ,
    id_digits   = raw_digits - rank_digits ,
    id_mask     = ~raw_key_type(0) >> rank_digits
  };

  raw_key_type key ;

  struct {
    raw_key_type id   : id_digits ;
    raw_key_type rank : rank_digits ;
  } normal_view ;

  struct {
    raw_key_type rank : rank_digits ;
    raw_key_type id   : id_digits ;
  } reverse_view ;

public:
  /** \brief  Destructor */
  ~EntityKey() {}

  /** Default constructor.
   * Note that entity_key_valid(key) == false if key is default-constructed.
   */
  EntityKey() : key(invalid_key) { }

  EntityKey( const EntityKey & rhs ) : key( rhs.key ) {}

  EntityKey & operator = ( const EntityKey & rhs )
    { key = rhs.key ; return *this ; }

  /** Constructor
   *
   * \param entity_rank is required to lie in the range 0 to 255 (which is
   *    the limit of what can be stored in 8 bits). This limit may be raised
   *    if we decide to use more than 8 bits for encoding an entity-type.
   *
   * \param entity_id is required to lie in the range 1 to 2^id_digits.
   *
   * If entity_rank or entity_id lie outside these ranges an exception will
   * be thrown.
   */
  EntityKey( EntityRank entity_rank, raw_key_type entity_id );

  raw_key_type id() const { return key & id_mask ; }

  EntityRank rank() const { return key >> id_digits ; }

  EntityRank type() const { return rank(); }

  bool operator==(const EntityKey &rhs) const {
    return key == rhs.key;
  }

  bool operator!=(const EntityKey &rhs) const {
    return !(key == rhs.key);
  }

  bool operator<(const EntityKey &rhs) const {
    return key < rhs.key;
  }

  bool operator>(const EntityKey &rhs) const {
    return rhs.key < key;
  }

  bool operator<=(const EntityKey &rhs) const {
    return !(key < rhs.key);
  }

  bool operator>=(const EntityKey &rhs) const {
    return !(rhs.key < key);
  }

  //------------------------------
  // As safe and explict a conversion
  // as possible between the raw_key_type and value.

  explicit EntityKey( const raw_key_type * const value )
   : key( *value ) {}

  raw_key_type raw_key() const { return key ; }
};

// Functions for encoding / decoding entity keys.

/** \brief  Given an entity key, return an entity type (rank). */
inline
EntityRank entity_rank( const EntityKey & key ) {
  return key.rank();
}

/** \brief  Given an entity key, return the identifier for the entity.  */
inline
EntityId  entity_id( const EntityKey & key ) {
  return key.id();
}

/** \brief  Query if an entity key is valid */
inline
bool entity_key_valid( const EntityKey & key ) {
  return key != EntityKey();
}

inline
bool entity_id_valid( EntityKey::raw_key_type id ) {
  return 0 < id && id <= EntityKey().id();
}

inline
size_t hash_value( EntityKey key) {
  return boost::hash_value(key.raw_key());
}



} // namespace mesh
} // namespace stk_classic

#endif /* stk_mesh_EntityKey_hpp */

