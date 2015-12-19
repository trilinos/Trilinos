#include <stddef.h>                     // for size_t, nullptr
#include <memory>
#include <map>
#include <vector>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/Types.hpp>

namespace stk { namespace mesh { class BulkData; }}

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
        m_owningElement(owningElement) {};
    ~Particle() {};

    stk::mesh::Entity owning_element() const { return m_owningElement; }
    stk::mesh::Entity spherical_element() const { return m_sphericalElement; }
    stk::mesh::Entity node() const { return m_node; }

private:
    stk::mesh::Entity m_node;
    stk::mesh::Entity m_sphericalElement;
    stk::mesh::Entity m_owningElement;
};

typedef std::vector<std::unique_ptr<Particle> > ParticleVector;

class ParticleManager
{
public:
    ParticleManager() {};
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

}
}
