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

#ifndef stk_mesh_MeshBuilder_hpp
#define stk_mesh_MeshBuilder_hpp

#include <stk_mesh/base/BulkData.hpp>
#include <memory>
#include <vector>
#include <string>

namespace stk {
namespace mesh {

class MeshBuilder
{
public:
  MeshBuilder();
  explicit MeshBuilder(ParallelMachine comm);
  virtual ~MeshBuilder() = default;

  MeshBuilder& set_spatial_dimension(unsigned spatialDimension);
  MeshBuilder& set_entity_rank_names(const std::vector<std::string>& entityRankNames);

  MeshBuilder& set_communicator(ParallelMachine comm);
  MeshBuilder& set_aura_option(BulkData::AutomaticAuraOption auraOption);
  virtual MeshBuilder& set_add_fmwk_data(bool addFmwkData);
#ifndef STK_HIDE_DEPRECATED_CODE  // Delete after 2025-08-19
  STK_DEPRECATED MeshBuilder& set_field_data_manager(std::unique_ptr<FieldDataManager> fieldDataManager);
#endif

  MeshBuilder& set_bucket_capacity(unsigned bucketCapacity);
  MeshBuilder& set_initial_bucket_capacity(unsigned initialCapacity);
  MeshBuilder& set_maximum_bucket_capacity(unsigned maximumCapacity);

  MeshBuilder& set_upward_connectivity(bool onOrOff);
  MeshBuilder& set_symmetric_ghost_info(bool onOrOff);
  MeshBuilder& set_maintain_local_ids(bool onOrOff);

  virtual std::unique_ptr<BulkData> create();
  virtual std::unique_ptr<BulkData> create(std::shared_ptr<MetaData> metaData);

  virtual std::shared_ptr<MetaData> create_meta_data();

protected:
  std::shared_ptr<impl::AuraGhosting> create_aura_ghosting();

  ParallelMachine m_comm;
  bool m_haveComm;
  BulkData::AutomaticAuraOption m_auraOption;
  bool m_addFmwkData;
  std::unique_ptr<FieldDataManager> m_fieldDataManager;
  unsigned m_initialBucketCapacity;
  unsigned m_maximumBucketCapacity;
  unsigned m_spatialDimension;
  std::vector<std::string> m_entityRankNames;
  bool m_upwardConnectivity;
  bool m_symmetricGhostInfo;
  bool m_maintainLocalIds;
};

}
}

#endif // stk_mesh_MeshBuilder_hpp
