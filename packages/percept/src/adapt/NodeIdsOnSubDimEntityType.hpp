// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef adapt_NodeIdsOnSubDimEntityType_hpp
#define adapt_NodeIdsOnSubDimEntityType_hpp

#include <iostream>
#include <stdexcept>
#include <string>
#include <sstream>
#include <vector>

#include <stk_mesh/base/Entity.hpp>
#include <percept/stk_mesh.hpp>
#include <stk_util/environment/CPUTime.hpp>

#include <percept/NoMallocArray.hpp>
#include <percept/PerceptMesh.hpp>
#include <percept/Util.hpp>

#include <percept/PerceptBoostArray.hpp>
#include <adapt/SubDimCell.hpp>

  namespace percept {

    /// data on a sub-dim entity (global node ids on the entity, the owning element's id)
    // FIXME - don't inherit from vector
    struct NodeIdsOnSubDimEntityType : public std::vector<stk::mesh::Entity>
    {
      std::vector<stk::mesh::EntityId> m_entity_id_vector;
      unsigned m_mark;
      int m_owner_rank;  // not communicated, local only: used for choosing proc ownership of this data

      NodeIdsOnSubDimEntityType(unsigned sz=1, stk::mesh::Entity allValues = stk::mesh::Entity(),
                                unsigned mark=0u) : std::vector<stk::mesh::Entity>(sz,allValues),
                                                  m_entity_id_vector(sz,0u),
                                                  m_mark(mark),m_owner_rank(-1)  {}
      void resize(size_t sz)
      {
        m_entity_id_vector.resize(sz);
        std::vector<stk::mesh::Entity>::resize(sz);
      }

      void pack(percept::PerceptMesh& eMesh, stk::CommBuffer& buff)
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
      void unpack(percept::PerceptMesh& /*eMesh*/, stk::CommBuffer& buff)
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

    inline std::ostream &operator<<(std::ostream& out, const NodeIdsOnSubDimEntityType nids)
    {
        for(unsigned iEnt=0;iEnt<nids.m_entity_id_vector.size();iEnt++)
            out << nids.m_entity_id_vector[iEnt] << " ";
        out << "m_mark : " << nids.m_mark << " ";
        out << "m_owner_rank : " << nids.m_owner_rank << " ";
        return out;
    }

  }

#endif
