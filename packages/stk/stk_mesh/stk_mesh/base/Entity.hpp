#ifndef STK_MESH_ENTIY_HPP
#define STK_MESH_ENTIY_HPP

#include <stdint.h>
#include <boost/functional/hash.hpp>

//----------------------------------------------------------------------

namespace sierra {
  namespace Fmwk {
    class MeshObjRoster;
  }
}

namespace stk{
  namespace mesh{

    struct Entity
    {
      enum Entity_t {
	InvalidEntity = 0ull,
	MinEntity = 1ull,
	MaxEntity = ~0ull,
      };

      uint64_t m_value;

      Entity operator=(Entity_t val) { m_value = val; return *this;}

      /** \brief local_offset is this entity's offset into all local entities of the same rank.
       * An entity's local_offset will generally remain unchanged through mesh-modification cycles,
       * which means that the set of local_offsets may not be compact or contiguous if some
       * entities have been deleted. (local_offsets for deleted entities are no longer valid.)
       * Thus, local_offset is not suitable for use as an equation index for linear-system operations.
       * See local_id() below.
       */
      size_t local_offset() const { return static_cast<size_t>(m_value); }

      bool is_local_offset_valid() const { return local_offset() > 0; }

      /** This setter should only be called by the BulkData class when creating entities.
       * Erroneous calls will lead to undefined (and probably disastrous) behavior.
       */
      void set_local_offset(size_t localOffset) {
	m_value = static_cast<Entity_t>(localOffset);
      }

      bool operator==(Entity entity) const { return m_value == entity.m_value; }

      bool operator==(Entity_t val) const { return m_value == val; }

      bool operator!=(Entity entity) const { return m_value != entity.m_value; }

      bool operator<(Entity entity) const { return m_value < entity.m_value; }

    };

  }
}


namespace stk {
  namespace mesh {
    namespace impl {
      class EntityRepository;
    }
    std::ostream & operator << ( std::ostream & , const Entity & );
    //
    // THIS CLASS STAYS BUT NEEDS TO BE DECLARED AFTER Entity.
    //
    class EntityEqual
    {
    public:
      bool operator()(const Entity lhs, const Entity rhs) const {
	return lhs == rhs;
      }
    };

    //
    // Entity
    //

    inline
    size_t hash_value( Entity entity) {
      return boost::hash_value(entity.local_offset());
    }



    /** \} */

  } // namespace mesh
} // namespace stk



#endif
