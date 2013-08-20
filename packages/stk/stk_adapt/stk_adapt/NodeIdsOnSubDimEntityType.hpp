#ifndef stk_adapt_NodeIdsOnSubDimEntityType_hpp
#define stk_adapt_NodeIdsOnSubDimEntityType_hpp

#include <iostream>
#include <stdexcept>
#include <string>
#include <sstream>
#include <vector>

#include <stk_mesh/base/Entity.hpp>
#include <stk_percept/stk_mesh.hpp>
#include <stk_util/environment/CPUTime.hpp>

#include <stk_percept/NoMallocArray.hpp>
#include <stk_percept/PerceptMesh.hpp>
#include <stk_percept/Util.hpp>

#include <stk_percept/PerceptBoostArray.hpp>
#include <stk_adapt/SubDimCell.hpp>

namespace stk {
  namespace adapt {

    // type defining what is stored on the edge
    typedef stk::mesh::Entity NodeIdsOnSubDimEntityTypeQuantum;

    /// data on a sub-dim entity (global node ids on the entity, the owning element's id)
    // FIXME - don't inherit from vector
    struct NodeIdsOnSubDimEntityType : public std::vector<NodeIdsOnSubDimEntityTypeQuantum>
    {
      typedef std::vector<NodeIdsOnSubDimEntityTypeQuantum> base_type;
      typedef std::vector<stk::mesh::EntityId> entity_id_vector_type;
      entity_id_vector_type m_entity_id_vector;
      unsigned m_mark;

      NodeIdsOnSubDimEntityType(unsigned sz=1, NodeIdsOnSubDimEntityTypeQuantum allValues=stk::mesh::Entity(),
                                unsigned mark=0u) : base_type(sz,allValues),
                                                    m_entity_id_vector(sz,0u),
                                                    m_mark(mark)  {}
      void resize(size_t sz)
      {
        m_entity_id_vector.resize(sz);
        base_type::resize(sz);
      }

      void pack(stk::percept::PerceptMesh& eMesh, stk::CommBuffer& buff)
      {
        buff.pack< unsigned > ( this->size() );
        m_entity_id_vector.resize( this->size() );
        for (unsigned ii = 0; ii < this->size(); ii++)
          {
            if ( ! eMesh.is_valid((*this)[ii]) )
              throw std::logic_error("logic err in NodeIdsOnSubDimEntityType::pack");
            stk::mesh::EntityId id = eMesh.identifier((*this)[ii]);
            VERIFY_OP_ON(id, != , 0, "logic err 2 in NodeIdsOnSubDimEntityType::pack");
            m_entity_id_vector[ii] = id;
            buff.pack<stk::mesh::EntityId>( id );
          }
      }
      void unpack(stk::percept::PerceptMesh& eMesh, stk::CommBuffer& buff)
      {
        unsigned sz;
        buff.unpack< unsigned > ( sz );
        this->resize( sz );
        m_entity_id_vector.resize( sz );
        for (unsigned ii = 0; ii < this->size(); ii++)
          {
            stk::mesh::EntityId id=0;
            buff.unpack<stk::mesh::EntityId>( id );
            m_entity_id_vector[ii] = id;
            VERIFY_OP_ON(id, != , 0, "logic err 3 in NodeIdsOnSubDimEntityType::unpack");
          }
      }
    };

  }
}
#endif
