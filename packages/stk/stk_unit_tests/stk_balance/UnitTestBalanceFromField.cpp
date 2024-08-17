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
#include "MeshFixtureRebalance.hpp"
#include <stk_unit_test_utils/unittestMeshUtils.hpp>
#include <stk_balance/balanceUtils.hpp>
#include <stk_balance/internal/Zoltan2ParallelGraph.hpp>
#include <stk_balance/internal/privateDeclarations.hpp>
#include <stk_balance/internal/Balancer.hpp>
#include <stk_balance/mesh/BalanceMesh.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/environment/EnvData.hpp>
#include "stk_balance/io/BalanceIO.hpp"
#include <vector>
#include <string>

namespace {

class IdAndTimeFieldValueSetter : public stk::unit_test_util::FieldValueSetter
{
public:
    virtual void populate_field(stk::mesh::BulkData &bulk, stk::mesh::FieldBase* field, const unsigned step,
                                const double time) const override
{
    stk::mesh::EntityRank fieldRank = field->entity_rank();

    std::vector<stk::mesh::Entity> entities;
    stk::mesh::get_entities(bulk, fieldRank, entities);

    stk::mesh::FieldVector allTransientFields = stk::io::get_transient_fields(bulk.mesh_meta_data());

    for(stk::mesh::FieldBase * transientField : allTransientFields)
    {
        for(size_t i = 0; i < entities.size(); i++)
        {
            unsigned numEntriesPerEntity = stk::mesh::field_scalars_per_entity(*transientField, entities[i]);
            double value = 100.0 * static_cast<double>(bulk.identifier(entities[i])) + time;
            double *data = static_cast<double*> (stk::mesh::field_data(*transientField, entities[i]));
            for(unsigned j=0; j<numEntriesPerEntity; j++)
                data[j] = value + j;
        }
    }
}
};

class BalanceFromField : public MeshFixtureRebalance
{
public:
  BalanceFromField() {
    stk::EnvData::instance().m_outputP0 = &stk::EnvData::instance().m_outputNull;
    testing::internal::CaptureStderr();
  }

  ~BalanceFromField() override {
    stk::EnvData::instance().m_outputP0 = &std::cout;
    testing::internal::GetCapturedStderr();
  }

  virtual void setup_initial_mesh_with_transient_field_data(const std::string & inputMeshSpec) override
  {
    m_transientTimeSteps = {0.0, 1.0, 2.0};
    m_transientFieldName = "weight_field";
    m_globalVariableName = "global_variable";
    stk::unit_test_util::generated_mesh_with_transient_data_to_file_in_serial(inputMeshSpec,
                                                                                             get_input_file_name(),
                                                                                             m_transientFieldName,
                                                                                             stk::topology::ELEM_RANK,
                                                                                             m_globalVariableName,
                                                                                             m_transientTimeSteps,
                                                                                             IdAndTimeFieldValueSetter());

    allocate_bulk(stk::mesh::BulkData::AUTO_AURA);
    get_meta().set_coordinate_field_name(m_balanceSettings.getCoordinateFieldName());
    m_ioBroker.property_add(Ioss::Property("DECOMPOSITION_METHOD", "RCB"));
    stk::io::fill_mesh_preexisting(m_ioBroker, get_input_file_name(), get_bulk());
  }

  void configure_balance_settings(const std::string & decompMethod = "rcb")
  {
    m_balanceSettings.setVertexWeightMethod(stk::balance::VertexWeightMethod::FIELD);
    m_balanceSettings.setVertexWeightFieldName(m_transientFieldName + "_scalar");
    m_balanceSettings.set_input_filename(get_input_file_name());
    m_balanceSettings.set_output_filename(get_output_file_name());
    m_balanceSettings.set_num_input_processors(1);
    m_balanceSettings.set_num_output_processors(get_parallel_size());
    m_balanceSettings.setDecompMethod(decompMethod);
  }
};

TEST_F(BalanceFromField, 6elems2procs_readLastTimeStepFromFile)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh_with_transient_field_data("1x1x6");
  configure_balance_settings();
  stk::balance::BalanceIO io(get_comm(), m_balanceSettings);
  const stk::balance::Balancer balancer(m_balanceSettings);
  stk::balance::BalanceMesh& initialMesh = io.initial_decomp();
//  balancer.balance(mesh);

  stk::mesh::BulkData & bulk = initialMesh.get_bulk();
  stk::mesh::Field<double> &weightField = *bulk.mesh_meta_data().get_field<double>(stk::topology::ELEM_RANK,
                                                                                   m_transientFieldName + "_scalar");
  const stk::mesh::BucketVector & elemBuckets = bulk.get_buckets(stk::topology::ELEM_RANK,
                                                                 bulk.mesh_meta_data().locally_owned_part());
  for (const stk::mesh::Bucket * bucket : elemBuckets) {
    for (const stk::mesh::Entity elem : *bucket) {
      const stk::mesh::EntityId elemId = bulk.identifier(elem);
      const double * fieldWeight = stk::mesh::field_data(weightField, elem);
      const double expectedFieldWeight = 100 * elemId + 2.0;  // Fixture scales ID by 100 and adds time
      EXPECT_DOUBLE_EQ(*fieldWeight, expectedFieldWeight);
    }
  }
}

TEST_F(BalanceFromField, 6elems2procs_checkGeometricDecomp)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh_with_transient_field_data("1x1x6");
  configure_balance_settings("rcb");
  stk::balance::BalanceIO io(get_comm(), m_balanceSettings);
  const stk::balance::Balancer balancer(m_balanceSettings);
  stk::balance::BalanceMesh& mesh = io.initial_decomp();
  balancer.balance(mesh);

  stk::mesh::BulkData & bulk = mesh.get_bulk();

  std::vector<size_t> counts;
  stk::mesh::count_entities(bulk.mesh_meta_data().locally_owned_part(), bulk, counts);

  if (get_parallel_rank() == 0) {
    EXPECT_EQ(counts[stk::topology::ELEM_RANK], 4u);
  }
  else {
    EXPECT_EQ(counts[stk::topology::ELEM_RANK], 2u);
  }
}

TEST_F(BalanceFromField, 6elems2procs_checkGraphDecomp)
{
  if (stk::parallel_machine_size(get_comm()) != 2) return;

  setup_initial_mesh_with_transient_field_data("1x1x6");
  configure_balance_settings("parmetis");
  stk::balance::BalanceIO io(get_comm(), m_balanceSettings);
  const stk::balance::Balancer balancer(m_balanceSettings);
  stk::balance::BalanceMesh& mesh = io.initial_decomp();
  balancer.balance(mesh);

  stk::mesh::BulkData & bulk = mesh.get_bulk();

  std::vector<size_t> counts;
  stk::mesh::count_entities(bulk.mesh_meta_data().locally_owned_part(), bulk, counts);

  if (get_parallel_rank() == 0) {
    EXPECT_EQ(counts[stk::topology::ELEM_RANK], 2u);
  }
  else {
    EXPECT_EQ(counts[stk::topology::ELEM_RANK], 4u);
  }
}

}
