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

#ifndef STK_MESH_ENTIY_HPP
#define STK_MESH_ENTIY_HPP

#include <stddef.h>                     // for size_t
#include <stdint.h>                     // for uint64_t, uint32_t
#include <iosfwd>                       // for ostream
#include <stk_util/stk_config.h>

#define STK_32BIT_ENTITY

namespace stk{
namespace mesh{

struct Entity
{
    enum Entity_t {
#ifdef STK_32BIT_ENTITY
        InvalidEntity = 0u,
        MinEntity = 1u,
        MaxEntity = ~0u
#else
        InvalidEntity = 0ull,
        MinEntity = 1ull,
        MaxEntity = ~0ull
#endif
    };

#ifdef STK_32BIT_ENTITY
    typedef uint32_t entity_value_type;
#else
    typedef uint64_t entity_value_type;
#endif

    entity_value_type m_value;

    KOKKOS_FUNCTION
    Entity() : m_value(InvalidEntity) {}

    KOKKOS_FUNCTION
    Entity(Entity_t value) : m_value(value) {}

    KOKKOS_FUNCTION
    explicit Entity(entity_value_type value) : m_value(value) {}

    KOKKOS_FUNCTION
    Entity operator=(entity_value_type val) { m_value = val; return *this;}

    /** \brief local_offset is this entity's offset into all local entities of the same rank.
     * An entity's local_offset will generally remain unchanged through mesh-modification cycles,
     * which means that the set of local_offsets may not be compact or contiguous if some
     * entities have been deleted. (local_offsets for deleted entities are no longer valid.)
     * Thus, local_offset is not suitable for use as an equation index for linear-system operations.
     * See local_id() below.
     */
    KOKKOS_FUNCTION
    entity_value_type local_offset() const { return m_value; }

    KOKKOS_FUNCTION
    bool is_local_offset_valid() const { return local_offset() > 0; }

    /** This setter should only be called by the BulkData class when creating entities.
     * Erroneous calls will lead to undefined (and probably disastrous) behavior.
     */
    KOKKOS_FUNCTION
    void set_local_offset(size_t localOffset) {
        m_value = static_cast<entity_value_type>(localOffset);
    }

    KOKKOS_FUNCTION
    bool operator==(Entity entity) const { return m_value == entity.m_value; }

    KOKKOS_FUNCTION
    bool operator==(entity_value_type val) const { return m_value == val; }

    KOKKOS_FUNCTION
    bool operator!=(Entity entity) const { return m_value != entity.m_value; }

    KOKKOS_FUNCTION
    bool operator!=(entity_value_type val) const { return m_value != val; }

    KOKKOS_FUNCTION
    bool operator<(Entity entity) const { return m_value < entity.m_value; }

    KOKKOS_FUNCTION
    bool operator>(Entity entity) const { return m_value > entity.m_value; }
};

std::ostream & operator << ( std::ostream & , const Entity & );

inline
size_t hash_value( Entity entity) {
//        try to use std::hash
    return entity.local_offset();
}

} // namespace mesh
} // namespace stk


#endif
