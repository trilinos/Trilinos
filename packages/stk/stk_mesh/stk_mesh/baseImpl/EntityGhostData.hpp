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

#ifndef stk_mesh_baseImpl_EntityGhostData_hpp
#define stk_mesh_baseImpl_EntityGhostData_hpp

#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/Entity.hpp>

namespace stk {
namespace mesh {
class BulkData;

namespace impl {

struct EntityGhostData
{
    enum DIRECTION {
        INVALID,
        NONE,
        SEND,
        RECEIVE
    };
    enum GHOST_LEVEL {
        LOCALLY_OWNED = -1,
        SHARED = 0,
        AURA = 1
    };
    DIRECTION direction;
    int ghostingLevel;
    int processor;
    Entity entity;
    const BulkData * bulkData;

    EntityGhostData()
    : direction(INVALID)
    , ghostingLevel(-2)
    , processor(-1)
    , entity()
    , bulkData(NULL) { }

    // insert to any ostream-like s
    template<class OStream> friend inline OStream& operator << (OStream& s, const DIRECTION& dir)
    {    
        switch (dir) {
            case INVALID: 
                s << "INVALID";
                break; 
            case NONE:
                s << "NONE";
                break; 
            case SEND:
                s << "SEND";
                break; 
            case RECEIVE:
                s << "RECEIVE";
                break; 
            default:
                s << "INVALID";
        }
        return s;
    }
    template<class OStream> inline OStream& printGhostLevel(OStream& s, int gl) const
    {    
        switch (gl) {
            case LOCALLY_OWNED:
                s << "LOCALLY_OWNED";
                break; 
            case SHARED:
                s << "SHARED";
                break; 
            case AURA:
                s << "AURA";
                break; 
            default:
                s << "CUSTOM_" << (gl-1);
        }
        return s;
    }
    template<class OStream> friend inline OStream& operator << (OStream& s, const EntityGhostData& egd)
    {    
        if (egd.bulkData != NULL) {
            s << "(Entity_gid="; 
            s << egd.bulkData->identifier(egd.entity)
              << ", rank=" << static_cast<unsigned int>(egd.bulkData->entity_rank(egd.entity));
        }
        else {
            s << "(Entity_lid=";
            s << egd.entity;
        }
        s << ", direction=" << egd.direction
          << ", processor=" << egd.processor
          << ", ghosting level=";
        egd.printGhostLevel(s,egd.ghostingLevel);
        s << ")";
        return s; 
    }
};

} // namespace impl
} // namespace mesh
} // namespace stk

//----------------------------------------------------------------------
//----------------------------------------------------------------------

#endif // stk_mesh_baseImpl_EntityGhostData_hpp

