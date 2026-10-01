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

#include <gtest/gtest.h>
#include <stk_util/diag/StringUtil.hpp> // for make_lower
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_mesh/base/Field.hpp>      // for Field
#include <stk_mesh/base/MetaData.hpp>   // for MetaData, put_field
#include "Ioss_DBUsage.h"               // for DatabaseUsage::READ_MODEL, etc
#include "Ioss_ElementTopology.h"       // for NameList
#include "Ioss_Field.h"                 // for Field, etc
#include "Ioss_IOFactory.h"             // for IOFactory
#include "Ioss_NodeBlock.h"             // for NodeBlock
#include "Ioss_Region.h"                // for Region, NodeBlockContainer
#include "Ioss_Utils.h"                 // for Utils
#include "Ioss_SerializeIO.h"
#include "stk_io/DatabasePurpose.hpp"   // for DatabasePurpose::READ_MESH, etc
#include "stk_io/WriteMesh.hpp"
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_unit_test_utils/BuildMesh.hpp"
#include <stk_unit_test_utils/getOption.h>
#include <stk_unit_test_utils/TextMesh.hpp>
#include <stddef.h>
#include <unistd.h>
#include <string>
#include <algorithm>
#include <cctype>
#include <cstdio>

namespace {
  stk::io::EntitySharingInfo get_sharing_info(stk::mesh::BulkData& bulkData)
  {
    stk::mesh::EntityVector sharedNodes;
    const bool sortById = true;
    stk::mesh::get_entities(bulkData, stk::topology::NODE_RANK, bulkData.mesh_meta_data().globally_shared_part(), sharedNodes, sortById);
    stk::io::EntitySharingInfo nodeSharingInfo;
    nodeSharingInfo.reserve(8*sharedNodes.size());

    std::vector<int> sharingProcs;
    for(stk::mesh::Entity sharedNode : sharedNodes)
    {
      bulkData.comm_shared_procs(bulkData.entity_key(sharedNode), sharingProcs);
      for(unsigned j=0;j<sharingProcs.size();++j)
        nodeSharingInfo.push_back(std::make_pair(bulkData.identifier(sharedNode), sharingProcs[j]));
    }

    return nodeSharingInfo;
  }

  template <typename T>
  void set_field_values(const stk::mesh::BulkData &bulkData, stk::mesh::Field<T> & field, T value)
  {
    const stk::mesh::BucketVector & buckets = bulkData.get_buckets(field.entity_rank(), field);
    auto fieldData = field.template data<stk::mesh::ReadWrite>();
    for (stk::mesh::Bucket * bucket : buckets) {
      auto bucketFieldData = fieldData.bucket_values(*bucket);
      for (stk::mesh::EntityIdx nodeIdx : bucket->entities()) {
        bucketFieldData(nodeIdx) = value;
      }
    }
  }

  void write_corrupt_restart_file_for_subdomain(const std::string &baseFilename,
                                                stk::mesh::Field<double> &field,
                                                int indexSubdomain,
                                                int numSubdomains,
                                                int globalNumNodes,
                                                int globalNumElems,
                                                stk::io::OutputParams& params,
                                                const stk::io::EntitySharingInfo &nodeSharingInfo,
                                                const int numSteps, const int skipStep)
  {
    Ioss::DatabaseIO *dbo = stk::io::create_database_for_subdomain(baseFilename, indexSubdomain, numSubdomains, false, Ioss::WRITE_RESULTS);
    Ioss::Region outRegion(dbo, "name");

    STK_ThrowRequireMsg(params.io_region_ptr() == nullptr, "OutputParams argument must have a NULL IORegion");
    params.set_io_region(&outRegion);
    stk::io::add_properties_for_subdomain(params, indexSubdomain, numSubdomains, globalNumNodes, globalNumElems);

    stk::io::write_mesh_data_for_subdomain(params, nodeSharingInfo);
    const stk::mesh::BulkData &bulkData = params.bulk_data();
    for(int i=0; i<=numSteps; i++) {
      if(stk::parallel_machine_rank(bulkData.parallel()) == 1 && (i == skipStep)) {
        continue;
      }

      double time = i;
      set_field_values(bulkData, field, time);
      stk::io::write_transient_data_for_subdomain(params, time);
    }

    stk::io::delete_selector_property(outRegion);
    params.set_io_region(nullptr);
  }
}

namespace stk
{
namespace unit_test_util
{
  std::string get_corrupt_restart_mesh_spec(stk::ParallelMachine communicator)
  {
    std::ostringstream oss;
    oss << "1x1x";
    oss << stk::parallel_machine_size(communicator);
    return oss.str();
  }

  std::string get_corrupt_restart_mesh_filename(stk::ParallelMachine communicator)
  {
    std::ostringstream oss;
    oss << "generated:1x1x";
    oss << stk::parallel_machine_size(communicator);
    return oss.str();
  }


  void create_corrupt_restart(stk::ParallelMachine communicator,
                              const std::string& restartFilename,
                              const std::string& internalClientFieldName,
                              const int nSteps, const int skipStep)
{
  std::string parallelFilename;

  stk::io::StkMeshIoBroker stkIo(communicator);
  const std::string exodusFileName = get_corrupt_restart_mesh_filename(communicator);
  size_t index = stkIo.add_mesh_database(exodusFileName, stk::io::READ_MESH);
  stkIo.set_active_mesh(index);
  stkIo.create_input_mesh();

  stk::mesh::MetaData &stkMeshMetaData = stkIo.meta_data();

  const int numberOfStates = 2;
  stk::mesh::Field<double> &field0 = stkMeshMetaData.declare_field<double>(stk::topology::NODE_RANK,
                                                                           internalClientFieldName,
                                                                           numberOfStates);

  stk::mesh::put_field_on_mesh(field0, stkMeshMetaData.universal_part(), nullptr);

  int numStatesToWrite = std::max(numberOfStates-1, 1);
  for(int state=0; state < numStatesToWrite; state++) {
    stk::mesh::FieldState state_identifier = static_cast<stk::mesh::FieldState>(state);
    stk::mesh::FieldBase *statedField = field0.field_state(state_identifier);
    stk::io::set_field_role(*statedField, Ioss::Field::TRANSIENT);
  }

  stkIo.populate_bulk_data();

  stk::mesh::BulkData &stkMeshBulkData = stkIo.bulk_data();
  stk::io::EntitySharingInfo nodeSharingInfo = get_sharing_info(stkMeshBulkData);

  std::vector<size_t> counts;
  stk::mesh::comm_mesh_counts(stkMeshBulkData, counts);
  int global_num_nodes = counts[stk::topology::NODE_RANK];
  int global_num_elems = counts[stk::topology::ELEM_RANK];

  stk::io::OutputParams params(stkMeshBulkData);
  write_corrupt_restart_file_for_subdomain(restartFilename,
                                           field0,
                                           stkMeshBulkData.parallel_rank(),
                                           stkMeshBulkData.parallel_size(),
                                           global_num_nodes,
                                           global_num_elems,
                                           params,
                                           nodeSharingInfo,
                                           nSteps,
                                           skipStep);
}

} // namespace unit_test_util
} // namespace stk


