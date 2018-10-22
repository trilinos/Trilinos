// #######################  Start Clang Header Tool Managed Headers ########################
// clang-format off
#include "ParticleUtils.hpp"
#include <memory>                      // for unique_ptr
#include <stk_mesh/base/BulkData.hpp>  // for BulkData
#include <utility>                     // for move, pair
#include "stk_mesh/base/Entity.hpp"    // for Entity
#include "stk_mesh/base/Types.hpp"     // for EntityId, PartVector
namespace stk { namespace mesh { class Part; } }
// clang-format on
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk
{
namespace unit_test_util
{

static int BASE_NODE_OFFSET = 1000000;

stk::mesh::Entity ParticleManager::create_particle_node(stk::mesh::EntityId id, stk::mesh::BulkData & mesh)
{
    const stk::mesh::PartVector noParts;
    return mesh.declare_node(id, noParts);
}

stk::mesh::Entity ParticleManager::create_particle_element(stk::mesh::EntityId id,
                                                           stk::mesh::BulkData & mesh,
                                                           stk::mesh::Part &particlePart,
                                                           stk::mesh::Entity connectedNode)
{
    stk::mesh::Entity newParticle = mesh.declare_element(id, stk::mesh::ConstPartVector{&particlePart});
    mesh.declare_relation(newParticle, connectedNode, 0);
    return newParticle;
}

void ParticleManager::add_particle(stk::mesh::Entity owningElement,
                                   stk::mesh::Entity newNode,
                                   stk::mesh::Entity newElement)
{
    std::unique_ptr<Particle> newParticle(new Particle(newNode, newElement, owningElement));
    m_particleMap[owningElement].push_back(std::move(newParticle));
}

void ParticleManager::create_particle(stk::mesh::EntityId id,
                                      stk::mesh::Entity owningElement,
                                      stk::mesh::BulkData & mesh,
                                      stk::mesh::Part &particlePart)
{
    stk::mesh::Entity newNode = create_particle_node(id+BASE_NODE_OFFSET, mesh);
    stk::mesh::Entity newElement = create_particle_element(id, mesh, particlePart, newNode);
    add_particle(owningElement, newNode, newElement);
}

ParticleVector & ParticleManager::get_particle_vector(stk::mesh::Entity parentElement)
{
    return m_particleMap.at(parentElement);
}

unsigned ParticleManager::count_particles_in_element(stk::mesh::Entity parentElement)
{
    return get_particle_vector(parentElement).size();
}

unsigned ParticleManager::count_all_particles()
{
    unsigned total_count = 0;
    for (auto & map_entry : m_particleMap) {
        total_count += count_particles_in_element(map_entry.first);
    }
    return total_count;
}

}
}
