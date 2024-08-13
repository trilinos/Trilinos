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

#ifndef STK_PARTICLE_UTILS_H
#define STK_PARTICLE_UTILS_H

// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include <map>                       // for map
#include <memory>                    // for unique_ptr, allocator
#include <stk_mesh/base/Entity.hpp>  // for Entity
#include <stk_mesh/base/Types.hpp>   // for EntityId
#include <vector>                    // for vector
namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class Part; } }
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################


namespace stk
{
namespace unit_test_util
{

class Particle
{
public:
    Particle(stk::mesh::Entity node, stk::mesh::Entity sphericalElement, stk::mesh::Entity owningElement)
      : m_node(node),
        m_sphericalElement(sphericalElement),
        m_owningElement(owningElement)
    {}
    ~Particle() = default;

    stk::mesh::Entity owning_element() const { return m_owningElement; }
    stk::mesh::Entity spherical_element() const { return m_sphericalElement; }
    stk::mesh::Entity node() const { return m_node; }

private:
    stk::mesh::Entity m_node;
    stk::mesh::Entity m_sphericalElement;
    stk::mesh::Entity m_owningElement;
};

typedef std::vector<std::unique_ptr<Particle>> ParticleVector;

class ParticleManager
{
public:
    ParticleManager() = default;
    ~ParticleManager() = default;

    stk::mesh::Entity create_particle_node(stk::mesh::EntityId id, stk::mesh::BulkData & mesh);

    stk::mesh::Entity create_particle_element(stk::mesh::EntityId id, stk::mesh::BulkData & mesh, stk::mesh::Part &particlePart, stk::mesh::Entity connectedNode);

    void add_particle(stk::mesh::Entity owningElement, stk::mesh::Entity newNode, stk::mesh::Entity newElement);

    void create_particle(stk::mesh::EntityId id, stk::mesh::Entity owningElement, stk::mesh::BulkData & mesh, stk::mesh::Part &particlePart);

    ParticleVector & get_particle_vector(stk::mesh::Entity parentElement);

    unsigned count_particles_in_element(stk::mesh::Entity parentElement);

    unsigned count_all_particles();

private:
    std::map<stk::mesh::Entity, ParticleVector> m_particleMap;
};

namespace simple_fields {

class Particle
{
public:
    Particle(stk::mesh::Entity node, stk::mesh::Entity sphericalElement, stk::mesh::Entity owningElement)
      : m_node(node),
        m_sphericalElement(sphericalElement),
        m_owningElement(owningElement)
    {}
    ~Particle() = default;

    stk::mesh::Entity owning_element() const { return m_owningElement; }
    stk::mesh::Entity spherical_element() const { return m_sphericalElement; }
    stk::mesh::Entity node() const { return m_node; }

private:
    stk::mesh::Entity m_node;
    stk::mesh::Entity m_sphericalElement;
    stk::mesh::Entity m_owningElement;
};

typedef std::vector<std::unique_ptr<Particle>> ParticleVector;

class STK_DEPRECATED_MSG("Please use the non-simple_fields-namespaced version of this class instead")
ParticleManager
{
public:
    ParticleManager() = default;
    ~ParticleManager() = default;

    stk::mesh::Entity create_particle_node(stk::mesh::EntityId id, stk::mesh::BulkData & mesh);

    stk::mesh::Entity create_particle_element(stk::mesh::EntityId id, stk::mesh::BulkData & mesh, stk::mesh::Part &particlePart, stk::mesh::Entity connectedNode);

    void add_particle(stk::mesh::Entity owningElement, stk::mesh::Entity newNode, stk::mesh::Entity newElement);

    void create_particle(stk::mesh::EntityId id, stk::mesh::Entity owningElement, stk::mesh::BulkData & mesh, stk::mesh::Part &particlePart);

    ParticleVector & get_particle_vector(stk::mesh::Entity parentElement);

    unsigned count_particles_in_element(stk::mesh::Entity parentElement);

    unsigned count_all_particles();

private:
    std::map<stk::mesh::Entity, ParticleVector> m_particleMap;
};

} // namespace simple_fields

}
}

#endif
