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

#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/baseImpl/AuraGhostingDownwardConnectivity.hpp>

namespace stk {
namespace mesh {

MeshBuilder::MeshBuilder()
 : m_comm(MPI_COMM_NULL),
   m_haveComm(false),
   m_auraOption(BulkData::AUTO_AURA),
   m_addFmwkData(false),
   m_initialBucketCapacity(stk::mesh::get_default_initial_bucket_capacity()),
   m_maximumBucketCapacity(stk::mesh::get_default_maximum_bucket_capacity()),
   m_spatialDimension(0),
   m_entityRankNames(),
   m_upwardConnectivity(true),
   m_symmetricGhostInfo(true),
   m_maintainLocalIds(false)
{
}

MeshBuilder::MeshBuilder(ParallelMachine comm)
 : m_comm(comm),
   m_haveComm(true),
   m_auraOption(BulkData::AUTO_AURA),
   m_addFmwkData(false),
   m_initialBucketCapacity(stk::mesh::get_default_initial_bucket_capacity()),
   m_maximumBucketCapacity(stk::mesh::get_default_maximum_bucket_capacity()),
   m_spatialDimension(0),
   m_entityRankNames(),
   m_upwardConnectivity(true),
   m_symmetricGhostInfo(true),
   m_maintainLocalIds(false)
{
}

MeshBuilder& MeshBuilder::set_spatial_dimension(unsigned spatialDimension)
{
  m_spatialDimension = spatialDimension;
  return *this;
}

MeshBuilder& MeshBuilder::set_entity_rank_names(const std::vector<std::string>& entityRankNames)
{
  m_entityRankNames = entityRankNames;
  return *this;
}

MeshBuilder& MeshBuilder::set_communicator(ParallelMachine comm)
{
  m_comm = comm;
  m_haveComm = true;
  return *this;
}

MeshBuilder& MeshBuilder::set_aura_option(BulkData::AutomaticAuraOption auraOption)
{
  m_auraOption = auraOption;
  return *this;
}

MeshBuilder& MeshBuilder::set_add_fmwk_data(bool addFmwkData)
{
  m_addFmwkData = addFmwkData;
  return *this;
}

#ifndef STK_HIDE_DEPRECATED_CODE  // Delete after 2025-08-19
STK_DEPRECATED MeshBuilder& MeshBuilder::set_field_data_manager(std::unique_ptr<FieldDataManager> fieldDataManager)
{
  m_fieldDataManager = std::move(fieldDataManager);
  return *this;
}
#endif

MeshBuilder& MeshBuilder::set_bucket_capacity(unsigned bucketCapacity)
{
  m_initialBucketCapacity = bucketCapacity;
  m_maximumBucketCapacity = bucketCapacity;
  return *this;
}

MeshBuilder& MeshBuilder::set_initial_bucket_capacity(unsigned initialCapacity)
{
  m_initialBucketCapacity = initialCapacity;
  return *this;
}

MeshBuilder& MeshBuilder::set_maximum_bucket_capacity(unsigned maximumCapacity)
{
  m_maximumBucketCapacity = maximumCapacity;
  return *this;
}

MeshBuilder& MeshBuilder::set_upward_connectivity(bool onOrOff)
{
  m_upwardConnectivity = onOrOff;
  return *this;
}

MeshBuilder& MeshBuilder::set_symmetric_ghost_info(bool onOrOff)
{
  m_symmetricGhostInfo = onOrOff;
  return *this;
}

MeshBuilder& MeshBuilder::set_maintain_local_ids(bool onOrOff)
{
  m_maintainLocalIds = onOrOff;
  return *this;
}

std::shared_ptr<MetaData> MeshBuilder::create_meta_data()
{
  if (m_spatialDimension > 0 || !m_entityRankNames.empty()) {
    return std::make_shared<MetaData>(m_spatialDimension, m_entityRankNames);
  }

  return std::make_shared<MetaData>();
}

std::shared_ptr<impl::AuraGhosting> MeshBuilder::create_aura_ghosting()
{
  if (m_upwardConnectivity) {
    return std::make_shared<impl::AuraGhosting>();
  }
  return std::make_shared<impl::AuraGhostingDownwardConnectivity>();
}

std::unique_ptr<BulkData> MeshBuilder::create(std::shared_ptr<MetaData> metaData)
{
  STK_ThrowRequireMsg(m_haveComm, "MeshBuilder must be given an MPI communicator before creating BulkData");

  std::unique_ptr<BulkData> bulkPtr(new BulkData(metaData,
                                                m_comm,
                                                m_auraOption,
#ifdef SIERRA_MIGRATION
                                                m_addFmwkData,
#endif
                                                std::move(m_fieldDataManager),
                                                m_initialBucketCapacity,
                                                m_maximumBucketCapacity,
                                                create_aura_ghosting(),
                                                m_upwardConnectivity));
  bulkPtr->set_symmetric_ghost_info(m_symmetricGhostInfo);
  bulkPtr->set_maintain_local_ids(m_maintainLocalIds);
  return bulkPtr;
}

std::unique_ptr<BulkData> MeshBuilder::create()
{
  return create(create_meta_data());
}

}
}
