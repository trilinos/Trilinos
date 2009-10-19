#ifndef stk_mesh_EntityKey_hpp
#define stk_mesh_EntityKey_hpp

#include <stdint.h>
#include <sstream>
#include <stdexcept>
#include <limits>

#include <boost/detail/endian.hpp>

namespace stk {
namespace mesh {

/** \addtogroup stk_mesh_module
 * \{
 */

//----------------------------------------------------------------------
/** \brief  Integer type for the entity identifiers */

// NOTE: Setting EntityId to 64-bits with a corresponding change of
// EntityKey.x.id to 48 bits fails with the pgi builds (version 7.1
// and 8.0).  A bug report has been submitted to PGI, but until this
// is resolved, I am resetting EntityId back to 32 bits and
// EntityKey.x.id to 32 bits also. GDS.
typedef uint32_t EntityId;
typedef unsigned EntityType;
typedef uint64_t EntityKeyValue;
 
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
 * EntityKey takes constructor arguments entity_type and entity_id which
 * are restricted to lie in a sub-range of what can be represented by the
 * EntityType and EntityId types. The sub-range allowed by EntityKey is
 * dictated by the amount of space used to store the internal encoding of
 * those values (described further in the comments for the EntityKey
 * constructor below).
 *
 * The function entity_key_valid attempts to determine that the key was
 * not created by some erroneous operation such as assignment from a
 * smaller type like a 32-bit int.
 *
 * Note that an instance of stk::mesh may place further restrictions on a
 * 'valid' key, such as requiring that
 *         0 <= entity_type(key) < meta_data.entity_type_count().
 *
 * Typically stk::mesh does not take EntityKeys as input. A user of
 * stk::mesh would (for instance) request that an Entity be created by
 * specifying an entity-type and entity-id as input. The resulting
 * Entity would then hold a mesh-created EntityKey that could be queried.
 * Thus stk::mesh can control the validity of EntityKeys associated with
 * the mesh.
 */

union EntityKey
{
  friend EntityType entity_type(EntityKey key);
  friend EntityId  entity_id(EntityKey key);
  friend bool entity_key_valid(EntityKey key);

private:
  EntityKeyValue        key;
  struct X {
#ifdef BOOST_LITTLE_ENDIAN
    EntityId            id    : 32;
    EntityId            fill  : 16;
    uint8_t             valid :  8;
    uint8_t             type  :  8;
#else
    uint8_t             type  :  8;
    uint8_t             valid :  8;
    EntityId            fill  : 16;
    EntityId            id    : 32;
#endif
  } x;

public:
  /** Default constructor.
   * Note that entity_key_valid(key) == false if key is default-constructed.
   */
  EntityKey()
    : key(0)
  { }

  /** Constructor
   * 
   * \param entity_type is required to lie in the range 0 to 255 (which is
   *    the limit of what can be stored in 8 bits). This limit may be raised
   *    if we decide to use more than 8 bits for encoding an entity-type.
   *
   * \param entity_id is required to lie in the range 1 to 2^31 (which is
   *    the limit of what can be stored in a 32-bit int). This limit may be
   *    raised if/when we use more than 32 bits for encoding an entity-id (and
   *    exodus supports types other than int for ids).
   *
   * If entity_type or entity_id lie outside these ranges an exception will
   * be thrown.
   */
  explicit EntityKey(EntityType entity_type, EntityId entity_id)
    : key(0)
  {
    if (entity_type > 255) {
      throw std::runtime_error("EntityKey ctor: entity_type must be <= 255 (to fit in 8-bits");
    }
    EntityId int_max = std::numeric_limits<int>::max();
    if (entity_id > int_max) {
      std::ostringstream msg;
      msg << "EntityKey ctor: entity_id ("<<entity_id
          <<") must be <= std::numeric_limits<int>::max() (which is "
          << std::numeric_limits<int>::max()<<")";
      std::string str = msg.str();
      throw std::runtime_error(str);
    }

    x.id = entity_id;
    x.valid = 1;
    x.type = entity_type;
  }

  EntityType type() const {
    return x.type;
  }
  
  EntityId id() const {
    return x.id;
  }

  bool operator==(const EntityKey &entity_key) const {
    return key == entity_key.key;
  }

  bool operator!=(const EntityKey &entity_key) const {
    return !(key == entity_key.key);
  }

  bool operator<(const EntityKey &entity_key) const {
    return key < entity_key.key;
  }

  bool operator>(const EntityKey &entity_key) const {
    return entity_key.key < key;
  }

  bool operator<=(const EntityKey &entity_key) const {
    return !(key < entity_key.key);
  }

  bool operator>=(const EntityKey &entity_key) const {
    return !(entity_key.key < key);
  }
};

// Functions for encoding / decoding entity keys.

/** \brief  Given an entity key, return an entity type (rank). */
inline
EntityType entity_type( EntityKey key ) {
  return key.x.type;
}

/** \brief  Given an entity key, return the identifier for the entity.  */
inline
EntityId  entity_id( EntityKey key ) {
  return key.x.id;
}

/** \brief  Query if an entity key is valid */
inline
bool entity_key_valid( EntityKey key ) {
  return key.x.valid && key.x.id != 0 && key.x.id <= (EntityId) std::numeric_limits<int>::max();
}

inline
bool entity_id_valid( EntityId id ) {
  return id != 0 && id <= (EntityId) std::numeric_limits<int>::max();
}

} // namespace mesh
} // namespace stk

#endif /* stk_mesh_EntityKey_hpp */

