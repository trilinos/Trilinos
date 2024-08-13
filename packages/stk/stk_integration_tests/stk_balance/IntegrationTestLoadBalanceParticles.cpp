#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <stk_unit_test_utils/ParticleUtils.hpp>
#include "stk_mesh/base/Field.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_unit_test_utils/getOption.h"
#include "stk_unit_test_utils/MeshFixture.hpp"
#include "stk_mesh/base/GetEntities.hpp"
#include <stk_balance/balance.hpp>
#include <stk_balance/balanceUtils.hpp>

namespace {

class StkRebalance : public stk::balance::FieldVertexWeightSettings
{
public:
  StkRebalance(stk::mesh::BulkData &stkMeshBulkData,
               const stk::balance::DoubleFieldType &weightField,
               const double defaultVertexWeight)
    : FieldVertexWeightSettings(stkMeshBulkData, weightField, defaultVertexWeight)
  {}

  virtual ~StkRebalance() = default;

private:
  StkRebalance() = delete;
  StkRebalance(const StkRebalance&) = delete;
  StkRebalance& operator=(const StkRebalance&) = delete;
};

class StkParticleRebalance : public StkRebalance
{
public:
  StkParticleRebalance(stk::unit_test_util::ParticleManager & particleManager,
                       stk::mesh::BulkData &stkMeshBulkData,
                       const stk::balance::DoubleFieldType &weightField,
                       const double defaultVertexWeight)
    : StkRebalance(stkMeshBulkData, weightField, defaultVertexWeight),
      m_particleManager(particleManager)
  {}

  virtual ~StkParticleRebalance() = default;

  virtual void modifyDecomposition(stk::balance::DecompositionChangeList & decomp) const
  {
    delete_particles_from_decomposition(decomp);
    move_particles_with_owning_element(decomp);
  }

private:
  void set_particle_destination(stk::balance::DecompositionChangeList & decomp, stk::mesh::Entity particle, const int destination) const
  {
    EXPECT_TRUE(decomp.get_bulk().is_valid(particle));
    EXPECT_EQ(decomp.get_bulk().bucket(particle).topology(), stk::topology::PARTICLE);
    decomp.set_entity_destination(particle, destination);
  }

  void set_particle_destination_from_owning_element(stk::balance::DecompositionChangeList & decomp, stk::mesh::Entity owner_element, const int destination) const
  {
    stk::unit_test_util::ParticleVector & vec = m_particleManager.get_particle_vector(owner_element);
    for (auto && particlePtr : vec) {
      set_particle_destination(decomp, particlePtr->spherical_element(), destination);
    }
  }

  void move_particles_to_owning_element_processor(stk::balance::DecompositionChangeList & decomp, const stk::mesh::EntityProc & ownerEntityProc) const
  {
    const stk::mesh::Entity ownerElement = ownerEntityProc.first;
    const int destination = ownerEntityProc.second;
    set_particle_destination_from_owning_element(decomp, ownerElement, destination);
  }

  void move_particles_with_owning_element(stk::balance::DecompositionChangeList & decomp) const
  {
    const stk::mesh::EntityProcVec & elementsToMove = decomp.get_all_partition_changes();
    for (auto elementAndDestination : elementsToMove) {
      move_particles_to_owning_element_processor(decomp, elementAndDestination);
    }
  }

  void delete_particles_from_decomposition(stk::balance::DecompositionChangeList & decomp) const
  {
    const stk::mesh::EntityProcVec & elementsToMove = decomp.get_all_partition_changes();
    for (auto elementAndDestination : elementsToMove) {
      const stk::mesh::Entity element = elementAndDestination.first;
      if (decomp.get_bulk().bucket(element).topology() == stk::topology::PARTICLE) {
        decomp.delete_entity(element);
      }
    }
  }

private:
  StkParticleRebalance() = delete;
  StkParticleRebalance(const StkParticleRebalance&) = delete;
  StkParticleRebalance& operator=(const StkParticleRebalance&) = delete;

  stk::unit_test_util::ParticleManager &m_particleManager;
};

class RebalanceParticleMesh : public stk::unit_test_util::MeshFixture
{
protected:
  RebalanceParticleMesh()
    : MeshFixture(),
      m_particleManager(),
      m_particleCountField(nullptr),
      m_particlePart(nullptr) {}
  virtual ~RebalanceParticleMesh() = default;

  void create_parts()
  {
    m_particlePart = & get_meta().declare_part_with_topology("Particles", stk::topology::PARTICLE, true);
  }

  void register_fields()
  {
    double init_value = 0.0;
    m_particleCountField = & get_meta().declare_field<double>(stk::topology::ELEM_RANK, "Particles", 1);
    stk::mesh::put_field_on_mesh(*m_particleCountField, get_meta().universal_part(), &init_value);
  }

  std::vector<int> get_particle_ids(int firstParticleId, int numParticles)
  {
    std::vector<int> particleIds(numParticles);
    for(int i=0; i<numParticles; ++i) {
      particleIds[i] = firstParticleId + i;
    }
    return particleIds;
  }

  std::vector<int> get_particle_ids_for_solid_element_1(int numGlobalElements)
  {
    int firstParticleId = numGlobalElements+1;
    int numParticles = numGlobalElements-1;
    return get_particle_ids(firstParticleId, numParticles);
  }

  std::vector<int> get_particle_info_for_other_solid_elements(int elementId, int numGlobalElements)
  {
    int firstParticleId = (numGlobalElements) + (numGlobalElements-1) + (elementId-1);
    int numParticles = 1;
    return get_particle_ids(firstParticleId, numParticles);
  }

  std::vector<int> get_particle_ids_for_solid_element(int elementId, int numGlobalElements)
  {
    if(elementId == 1)
      return get_particle_ids_for_solid_element_1(numGlobalElements);

    return get_particle_info_for_other_solid_elements(elementId, numGlobalElements);
  }


  void create_particles(const int numLocalElements)
  {
    int numGlobalElements = get_bulk().parallel_size()*numLocalElements;

    for(int elementId=1; elementId<=numGlobalElements; ++elementId)
    {
      std::vector<int> particleIds = get_particle_ids_for_solid_element(elementId, numGlobalElements);
      add_particles_to_element_and_update_weights(particleIds, elementId);
    }
  }

  void create_uneven_particle_distribution(const int numLocalElements)
  {
    get_bulk().modification_begin();
    create_particles(numLocalElements);
    get_bulk().modification_end();
  }

  void create_uneven_weight_distribution(const int numLocalElements)
  {
    int numGlobalElements = get_bulk().parallel_size()*numLocalElements;

    for (int elementId=1; elementId<=numGlobalElements; ++elementId) {
      stk::mesh::Entity element = get_bulk().get_entity(stk::mesh::EntityKey(stk::topology::ELEM_RANK,elementId));
      if (get_bulk().is_valid(element) && get_bulk().bucket(element).owned()) {
        double* weight = stk::mesh::field_data(*m_particleCountField, element);
        EXPECT_TRUE(nullptr != weight);
        if (elementId == 1) {
          *weight = static_cast<double>(numGlobalElements - 1);
        }
        else {
          *weight = 1.0;
        }
      }
    }
  }

  void update_vertex_weight(stk::mesh::Entity elem)
  {
    double* weight = stk::mesh::field_data(*m_particleCountField, elem);
    EXPECT_TRUE(nullptr != weight);
    *weight = static_cast<double> (m_particleManager.count_particles_in_element(elem));
  }

  void add_particles_to_element(const std::vector<int> & particleIds, stk::mesh::Entity elem)
  {
    for (int particleId : particleIds) {
      m_particleManager.create_particle(particleId, elem, get_bulk(), *m_particlePart);
    }
  }

  void add_particles_to_element_and_update_weights(const std::vector<int> & particleIds, int elementId)
  {
    stk::mesh::Entity elem = get_bulk().get_entity(stk::mesh::EntityKey(stk::topology::ELEM_RANK, elementId));
    if (get_bulk().is_valid(elem) && get_bulk().bucket(elem).owned()) {
      add_particles_to_element(particleIds, elem);
      update_vertex_weight(elem);
    }
  }

  void do_stk_particle_rebalance()
  {
    stk::mesh::Selector selector = !(*m_particlePart);
    const double defaultVertexWeight = 0.0;
    StkParticleRebalance graphSettings(m_particleManager, get_bulk(), *m_particleCountField, defaultVertexWeight);
    stk::balance::balanceStkMesh(graphSettings, get_bulk(), {selector});
  }

  void do_stk_rebalance()
  {
    stk::mesh::Selector selector = get_meta().universal_part();
    const double defaultVertexWeight = 0.0;
    StkRebalance graphSettings(get_bulk(), *m_particleCountField, defaultVertexWeight);
    stk::balance::balanceStkMesh(graphSettings, get_bulk(), {selector});
  }

  double get_total_weight_for_these_elements(const stk::mesh::EntityVector & solidElements)
  {
    double totalWeightTheseElements = 0.0;
    for (const stk::mesh::Entity element : solidElements) {
      double* weight = stk::mesh::field_data(*m_particleCountField, element);
      totalWeightTheseElements += (*weight);
    }
    return totalWeightTheseElements;
  }

  double get_total_element_weight_for_this_proc(const int numLocalElements)
  {
    stk::mesh::EntityVector solidElements;
    stk::mesh::Selector solidSelector = (!(*m_particlePart)) & get_meta().locally_owned_part();
    stk::mesh::get_entities(get_bulk(), stk::topology::ELEM_RANK, solidSelector, solidElements);
    return get_total_weight_for_these_elements(solidElements);
  }

  void check_expected_weight_distribution(const int numLocalElements)
  {
    const double expectedWeight = (2 * numLocalElements) - 1;
    const double actualWeightThisProc = get_total_element_weight_for_this_proc(numLocalElements);
    EXPECT_DOUBLE_EQ(expectedWeight, actualWeightThisProc);
  }

  void check_expected_particle_distribution(const int numLocalElements)
  {
    stk::mesh::EntityVector particleElements;
    stk::mesh::Selector particleSelector = (*m_particlePart) & get_meta().locally_owned_part();
    stk::mesh::get_entities(get_bulk(), stk::topology::ELEM_RANK, particleSelector, particleElements);

    const double expectedParticles = (2 * numLocalElements) - 1;
    const double actualParticlesThisProc = get_total_element_weight_for_this_proc(numLocalElements);
    EXPECT_DOUBLE_EQ(expectedParticles, actualParticlesThisProc);
  }

  std::string create_mesh_specification(const int numLocalElements)
  {
    std::ostringstream os;
    os << "generated:1x1x" << numLocalElements*stk::parallel_machine_size(get_comm());
    return os.str();
  }

  void run_stk_particle_rebalance_test(const int numLocalElements, stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    setup_empty_mesh(auraOption);
    create_parts();
    register_fields();
    stk::io::fill_mesh(create_mesh_specification(numLocalElements), get_bulk());
    create_uneven_particle_distribution(numLocalElements);
    do_stk_particle_rebalance();
    check_expected_weight_distribution(numLocalElements);
    check_expected_particle_distribution(numLocalElements);
  }

  void run_stk_rebalance_test(const int numLocalElements, stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    setup_empty_mesh(auraOption);
    create_parts();
    register_fields();
    stk::io::fill_mesh(create_mesh_specification(numLocalElements), get_bulk());
    create_uneven_weight_distribution(numLocalElements);
    do_stk_rebalance();
    check_expected_weight_distribution(numLocalElements);
  }

protected:
  stk::unit_test_util::ParticleManager m_particleManager;
  stk::balance::DoubleFieldType * m_particleCountField;
  stk::mesh::Part * m_particlePart;
};

// Characterization of Fuego rebalance misbehavior under Framework and STK
//
// Initial Mesh (ID.Proc):
// P == particle "owned" by containing element
//
//  o------------o------------o------------o------------o------------o------------o------------o------------o
//  |   P(9.0)   |            |            |            |            |            |            |            |
//  |     .      |   P(16.0)  |  P(17.0)   |   P(18.0)  |   P(19.1)  |   P(20.1)  |   P(21.1)  |   P(22.1)  |
//  |     .      |            |            |            |            |            |            |            |
//  |     .      |            |            |            |            |            |            |            |
//  |   P(15.0)  |            |            |            |            |            |            |            |
//  o------------o------------o------------o------------o------------o------------o------------o------------o
//       1.0          2.0          3.0          4.0          5.1          6.1          7.1          8.1
//               ^                                      ^
//               |                                      |
//               |                                      \------Initial processor boundary
//               \-------------------------------------------- Expected rebalanced processor boundary
//
//
// # Zoltan1 under Balance
// #---------------------------------------------------------------------------------
// With particles:
//  - Does not run due to inability to use Selectors
//
// Without particles:
//  EPP=2,3,4,5,6,7,8,9,10: No change
//  EPP=100: Moved 3 elements to p0
//  EPP=1000: Moved 6 elements to p0
//  EPP=10000: No change
//
//
// # Zoltan2 under Balance
// #---------------------------------------------------------------------------------
// With particles:
//  EPP= 2: Moves everything to p0
//  EPP= 3: Moves everything to p1
//  EPP= 4: Passes
//  EPP= 5: Passes
//  EPP= 6: Passes
//  EPP= 7: Passes
//  EPP= 8: Passes
//  EPP= 9: Passes
//  EPP=10: Passes
//  EPP=100,1000,10000: Passes
//
// Without particles:
//  EPP= 2: Moves everything to p1
//  EPP= 3: Moves everything to p1
//  EPP= 4: Passes
//  EPP= 5: Passes
//  EPP= 6: Passes
//  EPP= 7: Passes
//  EPP= 8: Passes
//  EPP= 9: Passes
//  EPP=10: Passes
//  EPP=100,1000,10000: Passes
//
// Legend
// -----------------------------
// EPP=Element Per Proc

int get_num_local_elements_from_cmdline()
{
  return stk::unit_test_util::get_command_line_option<int>("-nLocal", 100);
}

TEST_F(RebalanceParticleMesh, UnevenParticles2ProcWithAura)
{
  if (2 == stk::parallel_machine_size(get_comm()))
    run_stk_particle_rebalance_test(get_num_local_elements_from_cmdline(), stk::mesh::BulkData::AutomaticAuraOption::AUTO_AURA);
}

TEST_F(RebalanceParticleMesh, UnevenElementWeights2ProcWithAura)
{
  if (2 == stk::parallel_machine_size(get_comm()))
    run_stk_rebalance_test(get_num_local_elements_from_cmdline(), stk::mesh::BulkData::AutomaticAuraOption::AUTO_AURA);
}

TEST_F(RebalanceParticleMesh, UnevenParticles2ProcWithoutAura)
{
  if (2 == stk::parallel_machine_size(get_comm()))
    run_stk_particle_rebalance_test(get_num_local_elements_from_cmdline(), stk::mesh::BulkData::AutomaticAuraOption::NO_AUTO_AURA);
}

TEST_F(RebalanceParticleMesh, UnevenElementWeights2ProcWithoutAura)
{
  if (2 == stk::parallel_machine_size(get_comm()))
    run_stk_rebalance_test(get_num_local_elements_from_cmdline(), stk::mesh::BulkData::AutomaticAuraOption::NO_AUTO_AURA);
}

}
