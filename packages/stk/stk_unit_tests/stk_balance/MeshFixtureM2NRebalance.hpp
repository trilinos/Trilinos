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
#ifndef MESHFIXTUREM2NREBALANCE_HPP
#define MESHFIXTUREM2NREBALANCE_HPP

#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_unit_test_utils/TextMeshToFile.hpp>
#include <stk_balance/m2n/m2nRebalance.hpp>
#include <stk_balance/balanceUtils.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stk_unit_test_utils/BuildMesh.hpp>
#include <Ioss_IOFactory.h>
#include "Ioss_CommSet.h"
#include <vector>
#include <unistd.h>

using NodeIdSharingProc = std::pair<stk::mesh::EntityId, int>;
using stk::unit_test_util::build_mesh;

struct OriginalTopology {
  std::string blockName;
  std::string topologyName;
};

struct AssemblyGrouping {
  std::string assemblyPartName;
  int assemblyId;
  std::vector<std::string> subsetPartNames;
};


class MeshFixtureM2NRebalance : public stk::unit_test_util::MeshFixture
{
protected:
  MeshFixtureM2NRebalance()
    : MeshFixture()
  { }

  ~MeshFixtureM2NRebalance() override = default;

  std::string get_input_file_name() { return "TemporaryInputMesh.g"; }

  std::string get_output_file_name() { return "TemporaryOutputMesh.g"; }

  void read_serial_mesh_with_auto_decomp() {
    allocate_bulk(stk::mesh::BulkData::AUTO_AURA);
    m_ioBroker.property_add(Ioss::Property("DECOMPOSITION_METHOD", "RCB"));
    stk::io::fill_mesh_preexisting(m_ioBroker, get_input_file_name(), get_bulk());
  }

  void setup_initial_mesh(const std::string & inputMeshSpec)
  {
    stk::unit_test_util::generated_mesh_to_file_in_serial(inputMeshSpec, get_input_file_name());
    read_serial_mesh_with_auto_decomp();
  }

  void setup_initial_mesh_textmesh_override_topology(const std::string& inputMeshDesc,
                                                     const std::vector<OriginalTopology>& originalTopologies)
  {
    if (get_parallel_rank() == 0) {
      stk::unit_test_util::TextMeshToFile tMesh(MPI_COMM_SELF, stk::mesh::BulkData::AUTO_AURA);
      tMesh.setup_mesh(inputMeshDesc, get_input_file_name());

      for (const OriginalTopology & ot : originalTopologies) {
        stk::io::set_original_topology_type(*tMesh.get_meta().get_part(ot.blockName), ot.topologyName);
      }

      tMesh.write_mesh();
    }

    read_serial_mesh_with_auto_decomp();
  }

  void setup_initial_mesh_textmesh_add_assemblies(const std::string& inputMeshDesc,
                                                  const std::vector<AssemblyGrouping>& assemblies)
  {
    if (get_parallel_rank() == 0) {
      stk::unit_test_util::TextMeshToFile tMesh(MPI_COMM_SELF, stk::mesh::BulkData::AUTO_AURA);
      tMesh.setup_mesh(inputMeshDesc, get_input_file_name());

      for (const AssemblyGrouping & ag : assemblies) {
        stk::mesh::Part & assemblyPart = tMesh.get_meta().declare_part(ag.assemblyPartName, stk::topology::ELEM_RANK);
        stk::io::put_assembly_io_part_attribute(assemblyPart);
        tMesh.get_meta().set_part_id(assemblyPart, ag.assemblyId);

        for (const std::string & subsetPartName : ag.subsetPartNames) {
          stk::mesh::Part * subsetPart = tMesh.get_meta().get_part(subsetPartName);
          STK_ThrowRequire(subsetPart != nullptr);
          tMesh.get_meta().declare_part_subset(assemblyPart, *subsetPart);
        }
      }

      tMesh.write_mesh();
    }

    read_serial_mesh_with_auto_decomp();
  }

  void setup_initial_mesh_textmesh(const std::string & inputMeshDesc)
  {
    stk::unit_test_util::text_mesh_to_file_in_serial(inputMeshDesc, get_input_file_name());

    allocate_bulk(stk::mesh::BulkData::AUTO_AURA);
    m_ioBroker.property_add(Ioss::Property("DECOMPOSITION_METHOD", "RCB"));
    stk::io::fill_mesh_preexisting(m_ioBroker, get_input_file_name(), get_bulk());
  }

  void setup_initial_mesh_with_transient_field_data(const std::string & inputMeshSpec)
  {
    m_transientTimeSteps = {0.0, 1.0, 2.0};
    m_transientFieldName = "transient_field";
    m_globalVariableName = "global_variable";
    stk::unit_test_util::generated_mesh_with_transient_data_to_file_in_serial(inputMeshSpec,
                                                                                             get_input_file_name(),
                                                                                             m_transientFieldName,
                                                                                             stk::topology::NODE_RANK,
                                                                                             m_globalVariableName,
                                                                                             m_transientTimeSteps,
                                                                                             stk::unit_test_util::IdAndTimeFieldValueSetter());

    allocate_bulk(stk::mesh::BulkData::AUTO_AURA);
    m_ioBroker.property_add(Ioss::Property("DECOMPOSITION_METHOD", "RCB"));
    stk::io::fill_mesh_preexisting(m_ioBroker, get_input_file_name(), get_bulk());
  }

  virtual void rebalance_mesh_m2n(int numFinalProcs) = 0;

  void clean_up_temporary_files()
  {
    if (get_parallel_rank() == 0) {
      unlink(get_input_file_name().c_str());
      for (unsigned i = 0; i < m_numFinalProcs; ++i) {
        unlink(get_subdomain_filename(m_numFinalProcs, i).c_str());
      }
    }
  }

  void test_decomposed_mesh_element_distribution(const std::vector<unsigned> & elemsPerProc)
  {
    // Make sure all procs have written their files before p0 tries to read them
    MPI_Barrier(get_comm());

    ASSERT_EQ(m_numFinalProcs, elemsPerProc.size());

    if (get_parallel_rank() == 0) {
      for (size_t i = 0; i < elemsPerProc.size(); ++i) {
        test_num_elements_this_subdomain(elemsPerProc.size(), i, elemsPerProc[i]);
      }
    }
  }

  void test_decomposed_mesh_element_distribution_and_fields(const std::vector<unsigned> & elemsPerProc)
  {
    // Make sure all procs have written their files before p0 tries to read them
    MPI_Barrier(get_comm());

    ASSERT_EQ(m_numFinalProcs, elemsPerProc.size());

    if (get_parallel_rank() == 0) {
      for (size_t i = 0; i < elemsPerProc.size(); ++i) {
        stk::unit_test_util::MeshFromFile finalMesh(MPI_COMM_SELF);
        finalMesh.fill_from_serial(get_subdomain_filename(elemsPerProc.size(), i));

        stk::unit_test_util::TransientVerifier verifier(MPI_COMM_SELF);
        verifier.verify_time_steps(finalMesh, m_transientTimeSteps);
        verifier.verify_global_variables_at_each_time_step(finalMesh, m_globalVariableName, m_transientTimeSteps);
        verifier.verify_num_transient_fields(finalMesh, 2);
        verifier.verify_transient_field_names(finalMesh, m_transientFieldName);
        verifier.verify_transient_fields(finalMesh);
      }
    }

    test_decomposed_mesh_element_distribution(elemsPerProc);
  }

  void test_decomposed_mesh_element_topology(const std::vector<OriginalTopology>& expectedOriginalTopology)
  {
    if (get_parallel_rank() == 0) {
      for (size_t subdomain = 0; subdomain < m_numFinalProcs; ++subdomain) {
        std::string dbtype("exodusII");
        const std::string subdomainFilename = get_subdomain_filename(m_numFinalProcs, subdomain);
        Ioss::DatabaseIO *dbo = Ioss::IOFactory::create(dbtype, subdomainFilename, Ioss::READ_MODEL, MPI_COMM_SELF);
        Ioss::Region region(dbo, get_subdomain_filename(m_numFinalProcs, subdomain));

        for (const OriginalTopology & ot : expectedOriginalTopology) {
          Ioss::GroupingEntity* ge = region.get_entity(ot.blockName);
          ASSERT_NE(ge, nullptr);

          Ioss::Property topologyProp = ge->get_property("original_topology_type");
          EXPECT_EQ(topologyProp.get_string(), ot.topologyName);
        }
      }
    }
  }

  void test_decomposed_mesh_assemblies(const std::vector<AssemblyGrouping>& expectedAssemblies)
  {
    if (get_parallel_rank() == 0) {
      for (size_t subdomain = 0; subdomain < m_numFinalProcs; ++subdomain) {
        stk::unit_test_util::MeshFromFile finalMesh(MPI_COMM_SELF);
        finalMesh.fill_from_serial(get_subdomain_filename(m_numFinalProcs, subdomain));

        for (const AssemblyGrouping & ag : expectedAssemblies) {
          stk::mesh::Part * assemblyPart = finalMesh.meta.get_part(ag.assemblyPartName);
          ASSERT_NE(assemblyPart, nullptr) << "Missing assembly part " << ag.assemblyPartName;
          EXPECT_TRUE(stk::io::is_part_assembly_io_part(*assemblyPart));
          EXPECT_TRUE(stk::io::is_part_io_part(assemblyPart));
          EXPECT_EQ(assemblyPart->id(), ag.assemblyId);

          for (const std::string & subsetPartName : ag.subsetPartNames) {
            stk::mesh::Part * subsetPart = finalMesh.meta.get_part(subsetPartName);
            ASSERT_NE(subsetPart, nullptr);
            EXPECT_TRUE(assemblyPart->contains(*subsetPart));
          }
        }
      }
    }
  }

  void test_decomposed_mesh_node_sharing(const std::vector<std::vector<NodeIdSharingProc>> & expectedNodeSharingForEachProc)
  {
    ASSERT_EQ(m_numFinalProcs, expectedNodeSharingForEachProc.size());

    if (get_parallel_rank() == 0) {
      for (size_t subdomain = 0; subdomain < m_numFinalProcs; ++subdomain) {
        std::string dbtype("exodusII");
        const std::string subdomainFilename = get_subdomain_filename(m_numFinalProcs, subdomain);
        Ioss::DatabaseIO *dbo = Ioss::IOFactory::create(dbtype, subdomainFilename, Ioss::READ_MODEL, MPI_COMM_SELF);
        Ioss::Region region(dbo, get_subdomain_filename(m_numFinalProcs, subdomain));
        Ioss::CommSet* iocs = region.get_commset("commset_node");

        if (iocs != nullptr) {
          size_t numSharing = iocs->get_field("entity_processor").raw_count();
          std::vector<uint64_t> entityProc;
          iocs->get_field_data("entity_processor", entityProc);

          std::set<NodeIdSharingProc> expectedSharing(expectedNodeSharingForEachProc[subdomain].begin(),
                                                      expectedNodeSharingForEachProc[subdomain].end());

          std::set<NodeIdSharingProc> actualSharing;
          for (unsigned i = 0; i < numSharing; ++i) {
            actualSharing.emplace(entityProc[i*2], entityProc[i*2+1]);
          }

          EXPECT_EQ(expectedSharing, actualSharing) << "Incorrect node sharing in " << subdomainFilename;
        }
      }
    }
  }

  void test_num_elements_this_subdomain(unsigned numProcs, unsigned subdomainId, unsigned expectedNumElements)
  {
    std::shared_ptr<stk::mesh::BulkData> tempBulk = build_mesh(MPI_COMM_SELF);
    std::string subdomainFileName = get_subdomain_filename(numProcs, subdomainId);
    stk::io::fill_mesh(subdomainFileName, *tempBulk);

    stk::mesh::EntityVector entities;
    stk::mesh::get_entities(*tempBulk, stk::topology::ELEM_RANK, entities);

    const unsigned numElements = stk::mesh::count_entities(*tempBulk, stk::topology::ELEM_RANK, tempBulk->mesh_meta_data().universal_part());
    EXPECT_EQ(numElements, expectedNumElements) << "On subdomain " << subdomainId;
  }

  unsigned get_number_width(unsigned number)
  {
    return static_cast<unsigned>(std::log10(static_cast<double>(number))+1);
  }

  std::string get_subdomain_filename(unsigned numProcs, unsigned subdomainId)
  {
    std::ostringstream os;
    os << get_output_file_name();
    if (numProcs > 1) {
      int numberWidth = static_cast<int>(get_number_width(numProcs));
      os << "." << std::setfill('0') << std::setw(numberWidth) << numProcs << "." << subdomainId;
    }
    return os.str();
  }

  unsigned m_numFinalProcs;
  stk::mesh::Field<unsigned> *m_targetDecompField;
  stk::io::StkMeshIoBroker m_ioBroker;
  std::vector<double> m_transientTimeSteps;
  std::string m_transientFieldName;
  std::string m_globalVariableName;
};

#endif // MESHFIXTUREM2NREBALANCE_HPP
