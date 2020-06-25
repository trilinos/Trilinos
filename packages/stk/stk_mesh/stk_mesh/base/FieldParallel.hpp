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

#ifndef stk_mesh_FieldParallel_hpp
#define stk_mesh_FieldParallel_hpp

#include <stk_util/stk_config.h>
#include <stk_mesh/base/Types.hpp>      // for EntityProc
#include <stk_mesh/base/FieldTraits.hpp>  // for FieldTraits
#include <stk_mesh/base/FieldBase.hpp>  // for FieldBase
#include <stk_mesh/base/BulkData.hpp>
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/CommNeighbors.hpp>
#include <stk_util/util/ReportHandler.hpp>  // for ThrowRequireMsg

#include <stddef.h>                     // for size_t
#include <vector>                       // for vector

namespace stk { namespace mesh { class Ghosting; } }

namespace stk {
namespace mesh {

/**
 * This file contains some helper functions that are part of the Field API.
 */

/** Send field-data from entities to their ghosts, for a specified 'ghosting'.
 * For entities that are ghosted, this function updates field-data from the
 * original entity to the ghosts.
 */
void communicate_field_data(
  const Ghosting                        & ghosts ,
  const std::vector< const FieldBase *> & fields );

/** Copy data for the given fields, from owned entities to shared-but-not-owned entities.
 * I.e., shared-but-not-owned entities get an update of the field-data from the owned entity.
*/
inline
void copy_owned_to_shared( const BulkData& mesh,
                           const std::vector< const FieldBase *> & fields )
{
  communicate_field_data(*mesh.ghostings()[BulkData::SHARED], fields);
}

//----------------------------------------------------------------------

/** Sum/Max/Min (assemble) field-data for the specified fields on shared entities such that each shared entity
 * will have the same field values on each sharing proc.
 */
void parallel_sum(const BulkData& mesh, const std::vector<const FieldBase*>& fields, bool deterministic = true);
void parallel_max(const BulkData& mesh, const std::vector<const FieldBase*>& fields);
void parallel_min(const BulkData& mesh, const std::vector<const FieldBase*>& fields);


inline bool find_proc_before_index(const EntityCommInfoVector& infovec, int proc, int index)
{
    for(int i=0; i<index; ++i) {
        if (proc == infovec[i].proc) {
            return true;
        }
    }
    return false;
}

struct InfoRankLess {
  bool operator()(const EntityCommListInfo& info, EntityRank rank) const
  { return info.key.rank() < rank; }
  bool operator()(EntityRank rank, const EntityCommListInfo& info) const
  { return rank < info.key.rank(); }
};

inline void communicate_field_data(
  const BulkData                        & mesh ,
  const std::vector< const FieldBase *> & fields )
{
  const int parallel_size = mesh.parallel_size();
  if ( fields.empty() || parallel_size == 1) { return; }

  const int numFields = fields.size();

  std::vector<std::vector<unsigned char> > send_data(parallel_size);
  std::vector<std::vector<unsigned char> > recv_data(parallel_size);

  const EntityCommListInfoVector &comm_info_vec = mesh.internal_comm_list();

  std::vector<unsigned> send_sizes(parallel_size, 0);
  std::vector<unsigned> recv_sizes(parallel_size, 0);

  std::vector<std::pair<int,int> > fieldRange(fields.size(), std::make_pair(-1,-1));
  for(int fi=0; fi<numFields; ++fi) {
    fields[fi]->sync_to_host();
    fields[fi]->modify_on_host();
    EntityRank fieldRank = fields[fi]->entity_rank();
    EntityCommListInfoVector::const_iterator startIter = std::lower_bound(comm_info_vec.begin(), comm_info_vec.end(), fieldRank, InfoRankLess());
    EntityCommListInfoVector::const_iterator endIter = std::upper_bound(startIter, comm_info_vec.end(), fieldRank, InfoRankLess());
    fieldRange[fi].first = std::distance(comm_info_vec.begin(), startIter);
    fieldRange[fi].second = std::distance(comm_info_vec.begin(), endIter);
  }

  //this first loop calculates send_sizes and recv_sizes.
  for(int fi=0; fi<numFields; ++fi)
  {
      const FieldBase & f = *fields[fi];
      for(int i = fieldRange[fi].first; i<fieldRange[fi].second; ++i)
      {
          const MeshIndex& meshIndex = mesh.mesh_index(comm_info_vec[i].entity);
          const Bucket& bucket = *meshIndex.bucket;

          unsigned e_size = 0;

          const unsigned bucketId = bucket.bucket_id();
          const unsigned size = field_bytes_per_entity(f, bucketId);
          e_size += size;

          if(e_size == 0)
          {
              continue;
          }

          const bool owned = bucket.owned();
          if(owned)
          {
              const EntityCommInfoVector& infovec = comm_info_vec[i].entity_comm->comm_map;
              const int infovec_size = infovec.size();
              for(int j=0; j<infovec_size; ++j)
              {
                  const int proc = infovec[j].proc;

                  const bool proc_already_found = j>0 && find_proc_before_index(infovec, proc, j); 
                  if (!proc_already_found) {
                      send_sizes[proc] += e_size;
                  }
              }
          }
          else
          {
              const int owner = mesh.parallel_owner_rank(comm_info_vec[i].entity);

              recv_sizes[owner] += e_size;
          }
      }
  }

  //now size the send_data buffers
  size_t max_len = 0;
  for(int p=0; p<parallel_size; ++p)
  {
      if (send_sizes[p] > 0)
      {
          if (send_sizes[p] > max_len)
          {
              max_len = send_sizes[p];
          }
          send_data[p].resize(send_sizes[p]);
          send_sizes[p] = 0;
      }
  }

  //now pack the send buffers
  std::vector<unsigned char> field_data(max_len);
  unsigned char* field_data_ptr = field_data.data();

  for(int fi=0; fi<numFields; ++fi)
  {
      const FieldBase & f = *fields[fi];
      for(int i = fieldRange[fi].first; i<fieldRange[fi].second; ++i)
      {
          const MeshIndex& meshIndex = mesh.mesh_index(comm_info_vec[i].entity);
          const Bucket& bucket = *meshIndex.bucket;

          const bool owned = bucket.owned();

          unsigned e_size = 0;

          {
              const unsigned bucketId = bucket.bucket_id();
              const unsigned size = field_bytes_per_entity(f, bucketId);
              if (owned && size > 0)
              {
                  unsigned char * ptr = reinterpret_cast<unsigned char*>(stk::mesh::field_data(f, bucketId, meshIndex.bucket_ordinal, size));
                  std::memcpy(field_data_ptr+e_size, ptr, size);
              }
              e_size += size;
          }

          if(e_size == 0)
          {
              continue;
          }

          if(owned)
          {
              const EntityCommInfoVector& infovec = comm_info_vec[i].entity_comm->comm_map;
              const int infovec_size = infovec.size();
              for(int j=0; j<infovec_size; ++j)
              {
                  const int proc = infovec[j].proc;
    
                  const bool proc_already_found = j>0 && find_proc_before_index(infovec, proc, j);
                  if (!proc_already_found) {
                      unsigned char* dest_ptr = send_data[proc].data()+send_sizes[proc];
                      const unsigned char* src_ptr = field_data_ptr;
                      std::memcpy(dest_ptr, src_ptr, e_size);
                      send_sizes[proc] += e_size;
                  }
              }
          }
      }
  }

  for(int p=0; p<parallel_size; ++p)
  {
      if (recv_sizes[p] > 0)
      {
          recv_data[p].resize(recv_sizes[p]);
          recv_sizes[p] = 0;
      }
  }

  parallel_data_exchange_nonsym_known_sizes_t(send_data, recv_data, mesh.parallel());

  //now unpack and store the recvd data
  for(int fi=0; fi<numFields; ++fi)
  {
      const FieldBase & f = *fields[fi];

      for(int i = fieldRange[fi].first; i<fieldRange[fi].second; ++i)
      {
          const MeshIndex& meshIndex = mesh.mesh_index(comm_info_vec[i].entity);
          const Bucket& bucket = *meshIndex.bucket;
          if (bucket.owned()) {
              continue;
          }

          const int owner = mesh.parallel_owner_rank(comm_info_vec[i].entity);
          if(recv_data[owner].size() == 0)
          {
              continue;
          }

          {
              const unsigned bucketId = bucket.bucket_id();
              const unsigned size = field_bytes_per_entity(f, bucketId);
              if (size > 0)
              {
                  unsigned char * ptr = reinterpret_cast<unsigned char*>(stk::mesh::field_data(f, bucketId, meshIndex.bucket_ordinal, size));

                  std::memcpy(ptr, &(recv_data[owner][recv_sizes[owner]]), size);
                  recv_sizes[owner] += size;
              }
          }
      }
  }
}

void parallel_sum_including_ghosts(const BulkData & mesh, const std::vector<const FieldBase *> & fields);
void parallel_max_including_ghosts(const BulkData & mesh, const std::vector<const FieldBase *> & fields);
void parallel_min_including_ghosts(const BulkData & mesh, const std::vector<const FieldBase *> & fields);

} // namespace mesh
} // namespace stk

#endif

