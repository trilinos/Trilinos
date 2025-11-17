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

#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/DataTraits.hpp"  // for DataTraits
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/FieldBase.hpp"  // for FieldBase
#include "stk_mesh/base/Types.hpp"      // for ConnectivityOrdinal, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/parallel/ParallelComm.hpp"  // for CommBuffer
#include "stk_util/util/ReportHandler.hpp"  // for ThrowAssertMsg
#include <algorithm>                       // for operator<<
#include <sstream>                      // for operator<<, basic_ostream
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/base/EntityCommDatabase.hpp>
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/Relation.hpp>   // for Relation
#include <string>                       // for operator<<


namespace stk {
namespace mesh {
namespace impl {

//----------------------------------------------------------------------------

namespace {

unsigned count_parallel_consistent_parts(const MetaData & meta,
                                         const unsigned* first, const unsigned* last)
{
  unsigned count = 0;
  for (unsigned part_index=0; part_index < last-first; ++part_index) {
    const unsigned part_ordinal = first[part_index];
    if ( (part_ordinal != meta.locally_owned_part().mesh_meta_data_ordinal()) &&
         (part_ordinal != meta.globally_shared_part().mesh_meta_data_ordinal()) &&
         (meta.get_parts()[part_ordinal]->entity_membership_is_parallel_consistent()) ) {
      ++count;
    }
  }
  return count;
}

void pack_bucket_part_list(const Bucket & bucket, CommBuffer & buf )
{
    const MetaData & meta = bucket.mesh().mesh_meta_data();
    const std::pair<const unsigned *, const unsigned *>
      part_ordinals = bucket.superset_part_ordinals();
    buf.pack<unsigned>( count_parallel_consistent_parts(meta, part_ordinals.first, part_ordinals.second) );
    unsigned nparts = part_ordinals.second - part_ordinals.first;
    for (unsigned part_index=0; part_index < nparts; ++part_index) {
        const unsigned part_ordinal = part_ordinals.first[part_index];
        if ( (part_ordinal != meta.locally_owned_part().mesh_meta_data_ordinal()) &&
             (part_ordinal != meta.globally_shared_part().mesh_meta_data_ordinal()) &&
             (meta.get_parts()[part_ordinal]->entity_membership_is_parallel_consistent() )) {
            buf.pack<unsigned>(part_ordinal);
        }
    }
}

}

void pack_entity_info(const BulkData& mesh,
                      CommBuffer& buf,
                      const Entity entity,
                      bool onlyPackDownwardRelations)
{
  const EntityKey & key   = mesh.entity_key(entity);
  const unsigned    owner = mesh.parallel_owner_rank(entity);

  buf.pack<EntityKey>( key );
  buf.pack<unsigned>( owner );
  pack_bucket_part_list(mesh.bucket(entity), buf);

  const bool onlyCountDownwardRelations = onlyPackDownwardRelations;
  const unsigned tot_rel = mesh.count_relations(entity, onlyCountDownwardRelations);
  buf.pack<unsigned>( tot_rel );

  Bucket& bucket = mesh.bucket(entity);
  unsigned ebo   = mesh.bucket_ordinal(entity);

  STK_ThrowAssertMsg(mesh.is_valid(entity), "BulkData at " << &mesh << " does not know Entity " << entity.local_offset());
  const EntityRank end_rank = onlyPackDownwardRelations ? mesh.entity_rank(entity) : static_cast<EntityRank>(mesh.mesh_meta_data().entity_rank_count());

  for (EntityRank irank = stk::topology::BEGIN_RANK; irank < end_rank; ++irank)
  {
    const unsigned nrel = bucket.num_connectivity(ebo, irank);
    if (nrel > 0) {
      Entity const *rel_entities = bucket.begin(ebo, irank);
      ConnectivityOrdinal const *rel_ordinals = bucket.begin_ordinals(ebo, irank);
      Permutation const *rel_permutations = bucket.begin_permutations(ebo, irank);
      for ( unsigned i = 0 ; i < nrel ; ++i ) {
        if (mesh.is_valid(rel_entities[i])) {
          STK_ThrowAssert(rel_ordinals);
          buf.pack<EntityKey>( mesh.entity_key(rel_entities[i]) );
          buf.pack<unsigned>( rel_ordinals[i] );
          if (should_store_permutations(bucket.entity_rank(),irank)) {
            STK_ThrowAssert(rel_permutations);
            buf.pack<unsigned>( rel_permutations[i] );
          }
        }
      }
    }
  }
}

void unpack_entity_info(
  CommBuffer       & buf,
  const BulkData   & mesh ,
  EntityKey        & key ,
  int              & owner ,
  PartVector       & parts ,
  RelationVector& relations )
{
  unsigned nparts = 0 ;
  unsigned nrel = 0 ;

  buf.unpack<EntityKey>( key );
  buf.unpack<int>( owner );
  buf.unpack<unsigned>( nparts );

  parts.resize( nparts );

  for ( unsigned i = 0 ; i < nparts ; ++i ) {
    unsigned part_ordinal = ~0u ;
    buf.unpack<unsigned>( part_ordinal );
    parts[i] = & mesh.mesh_meta_data().get_part( part_ordinal );
  }

  buf.unpack( nrel );

  relations.clear();
  relations.reserve( nrel );

  for ( unsigned i = 0 ; i < nrel ; ++i ) {
    EntityKey rel_key ;
    unsigned rel_id = 0 ;
    unsigned rel_attr = 0 ;
    buf.unpack<EntityKey>( rel_key );
    buf.unpack<unsigned>( rel_id );
    if (should_store_permutations(key.rank(), rel_key.rank())) {
      buf.unpack<unsigned>( rel_attr );
    }
    Entity const entity =
      mesh.get_entity( rel_key.rank(), rel_key.id() );
    if ( mesh.is_valid(entity) ) {
      Relation rel(entity, mesh.entity_rank(entity), rel_id );
      rel.set_attribute(rel_attr);
      relations.push_back( rel );
    }
  }
}

//----------------------------------------------------------------------

struct SideSetInfo {
  unsigned partOrdinal;
  bool fromInput;
  std::vector<ConnectivityOrdinal> sideOrdinals;
};

void fill_sideset_info_for_entity(const MetaData& meta, const Entity entity, SideSet* sideset, std::vector<SideSetInfo>& sideSetInfo)
{
  std::vector<SideSetEntry>::iterator lowerBound = std::lower_bound(sideset->begin(), sideset->end(), SideSetEntry(entity, 0));
  std::vector<SideSetEntry>::iterator upperBound = std::upper_bound( sideset->begin(), sideset->end(), SideSetEntry(entity, INVALID_CONNECTIVITY_ORDINAL));
  unsigned distance = std::distance(lowerBound, upperBound);
  if (distance > 0) {
    SideSetInfo sideInfo;
    Part* part = meta.get_part(sideset->get_name());
    sideInfo.partOrdinal = part->mesh_meta_data_ordinal();
    sideInfo.fromInput = sideset->is_from_input();
    for (std::vector<SideSetEntry>::iterator iter = lowerBound; iter != upperBound; ++iter) {
      sideInfo.sideOrdinals.push_back(iter->side);
    }
    sideSetInfo.push_back(sideInfo);
  }
}

std::vector<SideSetInfo> get_sideset_info_for_entity(BulkData& mesh, const Entity entity) {
  const MetaData& meta = mesh.mesh_meta_data();
  std::vector<SideSet*> sidesets = mesh.get_sidesets();
  std::vector<SideSetInfo> sideSetInfo;
  for (SideSet* sideset : sidesets) {
    fill_sideset_info_for_entity(meta, entity, sideset, sideSetInfo);
  }
  return sideSetInfo;
}

void fill_comm_buffer_with_sideset_info(const std::vector<SideSetInfo>& sideSetInfo, CommBuffer& buf) {
  buf.pack<unsigned>(sideSetInfo.size());
  if (sideSetInfo.size() > 0) {
    for (const SideSetInfo& sideInfo : sideSetInfo) {
      buf.pack<unsigned>(sideInfo.partOrdinal);
      buf.pack<bool>(sideInfo.fromInput);
      buf.pack<unsigned>(sideInfo.sideOrdinals.size());
      for (const ConnectivityOrdinal sideOrdinal : sideInfo.sideOrdinals) {
        buf.pack<ConnectivityOrdinal>(sideOrdinal);
      }
    }
  }
}

void pack_sideset_info(BulkData& mesh, CommBuffer & buf, const Entity entity)
{
  if (mesh.entity_rank(entity) == stk::topology::ELEMENT_RANK)
  {
    std::vector<SideSetInfo> sideSetInfo = get_sideset_info_for_entity(mesh, entity);
    fill_comm_buffer_with_sideset_info(sideSetInfo, buf);
  }
}

void unpack_sideset_info(
  CommBuffer & buf,
  BulkData & mesh,
  const Entity entity)
{
  if (mesh.entity_rank(entity) == stk::topology::ELEMENT_RANK)
  {
    unsigned numSideSetInfo;
    buf.unpack<unsigned>( numSideSetInfo );

    if (numSideSetInfo > 0)
    {
      const stk::mesh::MetaData& meta = mesh.mesh_meta_data();
      for (unsigned sideInfo = 0; sideInfo < numSideSetInfo; ++sideInfo)
      {
        unsigned partOrdinal;
        bool fromInput;
        buf.unpack<unsigned>(partOrdinal);
        buf.unpack<bool>(fromInput);

        stk::mesh::Part& sidePart = meta.get_part(partOrdinal);
        if (!mesh.does_sideset_exist(sidePart))
        {
          mesh.create_sideset(sidePart, fromInput);
        }
        SideSet& sideset = mesh.get_sideset(sidePart);

        unsigned numSideOrdinals;
        buf.unpack<unsigned>(numSideOrdinals);
        for (unsigned side = 0; side < numSideOrdinals; ++side)
        {
          ConnectivityOrdinal sideOrdinal;
          buf.unpack<ConnectivityOrdinal>(sideOrdinal);
          sideset.add(entity, sideOrdinal);
        }
      }
    }
  }
}

//----------------------------------------------------------------------

void pack_field_values(const BulkData& mesh, CommBuffer& buf, Entity entity)
{
  if (!mesh.is_field_updating_active()) {
    return;
  }
  const Bucket& bucket = mesh.bucket(entity);
  const MetaData& mesh_meta_data = mesh.mesh_meta_data();
  const std::vector<FieldBase*>& fields = mesh_meta_data.get_fields(bucket.entity_rank());
  for (FieldBase* field : fields) {
    if (field->data_traits().is_pod) {
      auto& fieldDataBytes = field->data_bytes<const std::byte>();
      auto entityBytes = fieldDataBytes.entity_bytes(entity);
      const int numBytes = entityBytes.num_bytes();

#ifndef NDEBUG
      buf.pack<int>(numBytes);
#endif
      if (numBytes) {
        buf.pack_bytes(entityBytes.pointer(), numBytes, entityBytes.bytes_per_scalar(),
                       entityBytes.scalar_byte_stride());
      }
    }
  }
}

bool unpack_field_values(const BulkData& mesh, CommBuffer& buf, Entity entity, [[maybe_unused]] std::ostream& error_msg)
{
  if (!mesh.is_field_updating_active()) {
    return true;
  }
  const Bucket& bucket = mesh.bucket(entity);
  const MetaData& mesh_meta_data = mesh.mesh_meta_data();
  const std::vector<FieldBase*>& fields = mesh_meta_data.get_fields(bucket.entity_rank());
  bool ok = true;
  for (const FieldBase* field : fields) {
    if (field->data_traits().is_pod) {
      auto& fieldDataBytes = field->data_bytes<std::byte>();
      auto entityBytes = fieldDataBytes.entity_bytes(entity);
      const int numBytes = entityBytes.num_bytes();

#ifndef NDEBUG
      int receivedNumBytes = 0;
      buf.unpack<int>(receivedNumBytes);
      if (numBytes != receivedNumBytes) {
        if (ok) {
          ok = false;
          error_msg << mesh.identifier(entity);
        }
        error_msg << " " << field->name() << " " << numBytes << " != " << receivedNumBytes;
        buf.skip<unsigned char>(receivedNumBytes);
      }
#endif

      if (numBytes) {
        buf.unpack_bytes(entityBytes.pointer(), numBytes, entityBytes.bytes_per_scalar(),
                         entityBytes.scalar_byte_stride());
      }
    }
  }
  return ok;
}

} // namespace impl
} // namespace mesh
}

