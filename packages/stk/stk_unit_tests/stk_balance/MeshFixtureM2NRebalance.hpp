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
#include <stk_balance/internal/balanceMtoN.hpp>
#include <stk_balance/balanceUtils.hpp>
#include <stk_balance/setup/M2NParser.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <vector>
#include <unistd.h>

class MeshFixtureM2NRebalance : public stk::unit_test_util::MeshFixture
{
protected:
  MeshFixtureM2NRebalance()
    : MeshFixture()
  { }

  ~MeshFixtureM2NRebalance() override = default;

  std::string get_output_file_name() { return "junk.g"; }

  void setup_initial_mesh(const std::string & inputMeshSpec)
  {
    const std::string tempInputFilename = "TemporaryInputMesh.g";
    stk::unit_test_util::generated_mesh_to_file_in_serial(inputMeshSpec, tempInputFilename);

    allocate_bulk(stk::mesh::BulkData::NO_AUTO_AURA);
    create_target_decomp_field_on_entire_mesh();
    m_ioBroker.property_add(Ioss::Property("DECOMPOSITION_METHOD", "RCB"));
    stk::io::fill_mesh_preexisting(m_ioBroker, tempInputFilename, get_bulk());
  }

  void setup_initial_mesh_with_transient_field_data(const std::string & inputMeshSpec)
  {
    const std::string tempInputFilename = "TemporaryTransientInputMesh.g";
    stk::unit_test_util::generated_mesh_with_transient_data_to_file_in_serial(inputMeshSpec,
                                                                              tempInputFilename,
                                                                              "transient_field",
                                                                              "global_variable",
                                                                              {0.0, 1.0, 2.0},
                                                                              stk::unit_test_util::IdAndTimeFieldValueSetter());

    allocate_bulk(stk::mesh::BulkData::NO_AUTO_AURA);
    create_target_decomp_field_on_entire_mesh();
    m_ioBroker.property_add(Ioss::Property("DECOMPOSITION_METHOD", "RCB"));
    stk::io::fill_mesh_preexisting(m_ioBroker, tempInputFilename, get_bulk());
  }

  virtual void rebalance_mesh_m2n(int numFinalProcs) = 0;

  void create_target_decomp_field_on_entire_mesh()
  {
    m_targetDecompField = &get_meta().declare_field<stk::mesh::Field<double>>(stk::topology::ELEMENT_RANK, "TargetDecomp", 1);
    stk::mesh::put_field_on_mesh(*m_targetDecompField, get_meta().universal_part(), static_cast<double*>(nullptr));
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

      for (unsigned i = 0; i < m_numFinalProcs; ++i) {
        unlink(get_subdomain_filename(m_numFinalProcs, i).c_str());
      }
    }
  }

  void test_num_elements_this_subdomain(unsigned numProcs, unsigned subdomainId, unsigned expectedNumElements)
  {
    stk::mesh::MetaData tempMeta;
    stk::mesh::BulkData tempBulk(tempMeta, MPI_COMM_SELF);
    std::string subdomainFileName = get_subdomain_filename(numProcs, subdomainId);
    stk::io::fill_mesh(subdomainFileName, tempBulk);

    stk::mesh::EntityVector entities;
    stk::mesh::get_entities(tempBulk, stk::topology::ELEM_RANK, entities);

    const unsigned numElements = stk::mesh::count_entities(tempBulk, stk::topology::ELEM_RANK, tempMeta.universal_part());
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
  stk::mesh::Field<double> *m_targetDecompField;
  stk::io::StkMeshIoBroker m_ioBroker;
};

#endif // MESHFIXTUREM2NREBALANCE_HPP
