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

#include <stk_util/stk_config.h>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_util/parallel/ParallelComm.hpp>
#include <stk_util/parallel/CommNeighbors.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>  // for Reduce, ReduceSum, etc
#include <stk_util/parallel/Parallel.hpp>  // for parallel_machine_rank, etc
#include <stk_util/util/PairIter.hpp>   // for PairIter
#include <stk_util/util/SortAndUnique.hpp>   // for PairIter

#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Entity.hpp>     // for Entity
#include <stk_mesh/base/Ghosting.hpp>   // for Ghosting
#include <stk_mesh/base/Types.hpp>      // for PairIterEntityComm, etc
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/baseImpl/MeshCommImplUtils.hpp>

#include <utility>                      // for pair
#include <sstream>                      // for basic_ostream::operator<<, etc
#include <set>

namespace stk {
namespace mesh {

bool find_proc_before_index(const EntityCommInfoVector& infovec, int proc, int index)
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

void compute_field_entity_ranges(const EntityCommListInfoVector& commListInfo,
                                 const std::vector<const FieldBase*>& fields,
                                 std::vector<std::pair<int,int>>& fieldRange)
{
  int numFields = fields.size();
  for(int fi=0; fi<numFields; ++fi) {
    EntityRank fieldRank = fields[fi]->entity_rank();
    EntityCommListInfoVector::const_iterator startIter = std::lower_bound(commListInfo.begin(), commListInfo.end(), fieldRank, InfoRankLess());
    EntityCommListInfoVector::const_iterator endIter = std::upper_bound(startIter, commListInfo.end(), fieldRank, InfoRankLess());
    fieldRange[fi].first = std::distance(commListInfo.begin(), startIter);
    fieldRange[fi].second = std::distance(commListInfo.begin(), endIter);
  }
}

void do_initial_sync_to_host(const std::vector<const FieldBase*>& fields, Selector& selector)
{
  for(const FieldBase* fptr : fields) {
    fptr->sync_to_host();
    fptr->modify_on_host(selector);
  }
}

//----------------------------------------------------------------------
void communicate_field_data(const Ghosting& ghosts, const std::vector<const FieldBase*>& fields)
{
  if ( fields.empty() ) { return; }

  const BulkData & mesh = ghosts.mesh();
  const int parallel_size = mesh.parallel_size();
  const int parallel_rank = mesh.parallel_rank();
  const unsigned ghost_id = ghosts.ordinal();

  Selector ghostSelector = !mesh.mesh_meta_data().locally_owned_part();
  do_initial_sync_to_host(fields, ghostSelector);

  // Sizing for send and receive

  const unsigned zero = 0 ;
  std::vector<unsigned> send_size( parallel_size , zero );
  std::vector<unsigned> recv_size( parallel_size , zero );

  const EntityCommDatabase& commDB = mesh.internal_comm_db();

  for ( const EntityCommListInfo& ecli : mesh.internal_comm_list()) {
    Entity e = ecli.entity;
    const MeshIndex meshIdx = mesh.mesh_index(e);
    const unsigned bucketId = meshIdx.bucket->bucket_id();
    EntityRank erank = meshIdx.bucket->entity_rank();

    const bool owned = meshIdx.bucket->owned();

    unsigned e_size = 0 ;
    for ( const FieldBase* fptr : fields) {
      const FieldBase & f = *fptr ;

      if(is_matching_rank(f, erank)) {
        e_size += field_bytes_per_entity( f , bucketId );
      }
    }

    if (e_size == 0) {
      continue;
    }

    const int owner = mesh.parallel_owner_rank(ecli.entity);
    PairIterEntityComm info = commDB.comm(ecli.entity_comm);
    if ( owned ) {
      for (; !info.empty(); ++info) {
        if (info->ghost_id == ghost_id) {
          send_size[ info->proc ] += e_size ;
        }
      }
    }
    else {
      for (; !info.empty(); ++info) {
        if (info->ghost_id == ghost_id) {
          recv_size[ owner ] += e_size ;
          break;//jump out since we know we're only recving 1 msg for this entity from the 1-and-only owner
        }
      }
    }
  }

  std::vector<int> send_procs, recv_procs;
  for(int p=0; p<mesh.parallel_size(); ++p) {
      if (send_size[p] > 0) {
          send_procs.push_back(p);
      }
      if (recv_size[p] > 0) {
          recv_procs.push_back(p);
      }
  }

  CommNeighbors sparse(mesh.parallel(), send_procs, recv_procs);

  for(int p=0; p<mesh.parallel_size(); ++p) {
      if (send_size[p] > 0) {
          sparse.send_buffer(p).reserve(send_size[p]);
      }
  }

  // Send packing:
  for (int phase = 0; phase < 2; ++phase) {

    for ( const EntityCommListInfo& ecli : mesh.internal_comm_list()) {
      const int owner = mesh.parallel_owner_rank(ecli.entity);
      if ( (phase == 0 && owner == parallel_rank) || (phase == 1 && owner != parallel_rank) ) {
        Entity e = ecli.entity;
        const MeshIndex meshIdx = mesh.mesh_index(e);
        const unsigned bucketId = meshIdx.bucket->bucket_id();
        EntityRank erank = meshIdx.bucket->entity_rank();

        for (const FieldBase* fptr : fields) {
          const FieldBase & f = *fptr ;

          if(!is_matching_rank(f, erank)) continue;

          const unsigned size = field_bytes_per_entity( f , bucketId );

          if ( size ) {
            unsigned char * ptr =
              reinterpret_cast<unsigned char *>(stk::mesh::field_data( f , bucketId, meshIdx.bucket_ordinal, size ));

            PairIterEntityComm info = commDB.comm(ecli.entity_comm);
            if (phase == 0) { // send
              for (; !info.empty(); ++info) {
                if (info->ghost_id == ghost_id) {
                  CommBufferV & b = sparse.send_buffer( info->proc );
                  b.pack<unsigned char>( ptr , size );
                }
              }
            }
            else { //recv
              for (; !info.empty(); ++info) {
                if (info->ghost_id == ghost_id) {
                  CommBufferV & b = sparse.recv_buffer( owner );
                  b.unpack<unsigned char>( ptr , size );
                  break;
                }
              }
            }
          }
        }
      }
    }
    if (phase == 0) { sparse.communicate(); }
  }
}

void communicate_field_data(const BulkData& mesh ,
                            const std::vector< const FieldBase *> & fields)
{
  const int parallel_size = mesh.parallel_size();
  if ( fields.empty() || parallel_size == 1) {
    return;
  }

  const int myRank = mesh.parallel_rank();
  const int numFields = fields.size();

  const EntityCommListInfoVector &comm_info_vec = mesh.internal_comm_list();

  Selector ghostSelector = !mesh.mesh_meta_data().locally_owned_part();
  do_initial_sync_to_host(fields, ghostSelector);

  static std::vector<int> int_buffer;
  int_buffer.resize(4*parallel_size+2);
  std::fill(int_buffer.begin(), int_buffer.end(), 0);

  int* send_sizes = int_buffer.data();
  int* sendOffsets = send_sizes + parallel_size;
  int* recv_sizes = sendOffsets + parallel_size + 1;
  int* recvOffsets = recv_sizes + parallel_size;
  STK_ThrowAssertMsg(static_cast<size_t>(std::distance(send_sizes, recvOffsets+parallel_size+1)) == int_buffer.size(),
             "communicate_field_data: sizing error in poor-man's memory-pool.");

  std::vector<std::pair<int,int> > fieldRange(fields.size(), std::make_pair(-1,-1));
  compute_field_entity_ranges(comm_info_vec, fields, fieldRange);

  const EntityCommDatabase& commDB = mesh.internal_comm_db();
  std::vector<int> commProcs;
  //this first loop calculates send_sizes and recv_sizes.
  for(int fi=0; fi<numFields; ++fi)
  {
      const FieldBase & f = *fields[fi];
      for(int i = fieldRange[fi].first; i<fieldRange[fi].second; ++i)
      {
          const MeshIndex& meshIndex = mesh.mesh_index(comm_info_vec[i].entity);
          const Bucket& bucket = *meshIndex.bucket;

          const unsigned bucketId = bucket.bucket_id();
          const unsigned e_size = field_bytes_per_entity(f, bucketId);

          if (e_size == 0) {
              continue;
          }

          const bool owned = bucket.owned();
          if(owned)
          {
              impl::fill_sorted_procs(commDB.comm(comm_info_vec[i].entity_comm), commProcs);
              for(int proc : commProcs) {
                send_sizes[proc] += e_size;
              }
          }
          else {
              const int owner = mesh.parallel_owner_rank(comm_info_vec[i].entity);
              recv_sizes[owner] += e_size;
          }
      }
  }

  size_t totalSendSize = 0;
  size_t totalRecvSize = 0;
  for (int p=0; p<parallel_size; ++p) {
    sendOffsets[p] = totalSendSize;
    totalSendSize += send_sizes[p];
    send_sizes[p] = 0;

    recvOffsets[p] = totalRecvSize;
    totalRecvSize += recv_sizes[p];
    recv_sizes[p] = 0;
  }

  sendOffsets[parallel_size] = totalSendSize;
  recvOffsets[parallel_size] = totalRecvSize;

  std::unique_ptr<unsigned char[]> uninitializedData(new unsigned char[totalSendSize + totalRecvSize]);

  unsigned char* send_data = uninitializedData.get();
  unsigned char* recv_data = send_data+totalSendSize;

  //now pack the send buffers
  for(int fi=0; fi<numFields; ++fi)
  {
      const FieldBase & f = *fields[fi];
      for(int i = fieldRange[fi].first; i<fieldRange[fi].second; ++i)
      {
          const MeshIndex& meshIndex = mesh.mesh_index(comm_info_vec[i].entity);
          const Bucket& bucket = *meshIndex.bucket;

          if (bucket.owned()) {
              const unsigned bucketId = bucket.bucket_id();
              const unsigned e_size = field_bytes_per_entity(f, bucketId);
              if (e_size > 0) {
                const unsigned char* field_data_ptr = reinterpret_cast<unsigned char*>(stk::mesh::field_data(f, bucketId, meshIndex.bucket_ordinal, e_size));

                impl::fill_sorted_procs(commDB.comm(comm_info_vec[i].entity_comm), commProcs);
                for(int proc : commProcs) {
                  unsigned char* dest_ptr = &(send_data[sendOffsets[proc]+send_sizes[proc]]);
                  std::memcpy(dest_ptr, field_data_ptr, e_size);
                  send_sizes[proc] += e_size;
                }
              }
          }
      }
  }

  parallel_data_exchange_nonsym_known_sizes_t(sendOffsets, send_data,
                                              recvOffsets, recv_data, mesh.parallel(),
                                              mesh.is_mesh_consistency_check_on());

  //now unpack and store the recvd data
  for(int fi=0; fi<numFields; ++fi)
  {
      const FieldBase & f = *fields[fi];

      for(int i = fieldRange[fi].first; i<fieldRange[fi].second; ++i)
      {
          const int owner = mesh.parallel_owner_rank(comm_info_vec[i].entity);
          if (owner == myRank) {
            continue;
          }

          const int recvSize = recvOffsets[owner+1]-recvOffsets[owner];
          if(recvSize == 0) {
              continue;
          }

          const MeshIndex& meshIndex = mesh.mesh_index(comm_info_vec[i].entity);
          const Bucket& bucket = *meshIndex.bucket;

          const unsigned bucketId = bucket.bucket_id();
          const unsigned size = field_bytes_per_entity(f, bucketId);
          if (size > 0)
          {
              unsigned char * ptr = reinterpret_cast<unsigned char*>(stk::mesh::field_data(f, bucketId, meshIndex.bucket_ordinal, size));

              std::memcpy(ptr, &(recv_data[recvOffsets[owner]+recv_sizes[owner]]), size);
              recv_sizes[owner] += size;
          }
      }
  }
}

//----------------------------------------------------------------------

/** Sum (assemble) field-data for the specified fields on shared entities such that each shared entity
 * will have the same field values on each sharing proc.
 */

namespace {

template <typename T, Operation OP>
struct DoOp;

template <typename T>
struct DoOp<T, Operation::SUM>
{
  T operator()(T lhs, T rhs) const
  { return lhs + rhs; }
};

template <typename T>
struct DoOp<T, Operation::MIN>
{
  T operator()(T lhs, T rhs) const
  { return lhs < rhs ? lhs : rhs; }
};

template <typename T>
struct DoOp<T, Operation::MAX>
{
  T operator()(T lhs, T rhs) const
  { return lhs > rhs ? lhs : rhs; }
};

template <typename T, Operation OP>
void parallel_op_impl(const BulkData& mesh, std::vector<const FieldBase*> fields, bool deterministic = false)
{
  if (fields.empty()) {
    return;
  }

  std::vector<int> comm_procs = mesh.all_sharing_procs(fields[0]->entity_rank());
  stk::mesh::EntityRank first_field_rank = fields[0]->entity_rank();
  for(size_t i=1; i<fields.size(); ++i) {
      const FieldBase* f = fields[i];
      if (f->entity_rank() != first_field_rank) {
          const std::vector<int>& sharing_procs = mesh.all_sharing_procs(f->entity_rank());
          for(int p : sharing_procs) {
              stk::util::insert_keep_sorted_and_unique(p, comm_procs);
          }
      }
  }

  auto msgPacker = [&fields, &mesh](int proc, std::vector<T>& send_data)
  {
    send_data.clear();
    size_t reserve_len = 0;
    for (size_t j = 0 ; j < fields.size() ; ++j ) {
        const FieldBase& f = *fields[j];
        STK_ThrowRequireMsg(f.type_is<T>(),
                      "Please don't mix fields with different primitive types in the same parallel assemble operation");

        f.sync_to_host();
        HostCommMapIndices commMapIndices = mesh.volatile_fast_shared_comm_map(f.entity_rank(), proc);
        for(size_t i=0; i<commMapIndices.extent(0); ++i) {
            const unsigned bucket = commMapIndices(i).bucket_id;
            const int num_bytes_per_entity = field_bytes_per_entity( f , bucket );
            if (num_bytes_per_entity > 0) {
                const int num_Ts_per_field = num_bytes_per_entity / sizeof(T);
                reserve_len += num_Ts_per_field;
            }
        }
    }
    send_data.reserve(reserve_len);

    for (size_t j = 0 ; j < fields.size() ; ++j ) {
        const FieldBase& f = *fields[j];
        HostCommMapIndices commMapIndices = mesh.volatile_fast_shared_comm_map(f.entity_rank(),proc);

        for(size_t i=0; i<commMapIndices.extent(0); ++i) {
            const unsigned bucket = commMapIndices(i).bucket_id;
            const unsigned offset = commMapIndices(i).bucket_ord;
            const int num_bytes_per_entity = field_bytes_per_entity( f , bucket );
            if (num_bytes_per_entity > 0) {
                const int num_Ts_per_field = num_bytes_per_entity / sizeof(T);
                const T* data = reinterpret_cast<T*>(stk::mesh::field_data( f , bucket)) + offset*num_Ts_per_field;

                for (int d = 0; d < num_Ts_per_field; ++d) {
                   send_data.push_back(data[d]);
                }
            }
        }
    }
  };

  auto msgUnpacker = [&fields, &mesh](int iproc, std::vector<T>& recv_data)
  {
    DoOp<T, OP> do_op;

    unsigned recv_offset = 0;
    for (size_t j = 0 ; j < fields.size() ; ++j ) {
        const FieldBase& f = *fields[j] ;
        HostCommMapIndices commMapIndices = mesh.volatile_fast_shared_comm_map(f.entity_rank(), iproc);

        f.sync_to_host();
        f.modify_on_host();

        for(size_t i=0; i<commMapIndices.extent(0); ++i) {
            const unsigned bucket = commMapIndices(i).bucket_id;
            const unsigned offset = commMapIndices(i).bucket_ord;
            const int num_bytes_per_entity = field_bytes_per_entity( f , bucket );
            if (num_bytes_per_entity > 0) {
                const int num_Ts_per_field = num_bytes_per_entity / sizeof(T);
                T* data = reinterpret_cast<T*>(stk::mesh::field_data( f , bucket))+offset*num_Ts_per_field;

                for (int d = 0; d < num_Ts_per_field; ++d) {
                    data[d] = do_op(data[d], recv_data[recv_offset++]);
                }
            }
        }
    }
  };

  MPI_Comm comm = mesh.parallel();
  stk::parallel_data_exchange_sym_pack_unpack<T>(comm, comm_procs, msgPacker, msgUnpacker, deterministic);
}

template <Operation OP>
inline
void parallel_op(const BulkData& mesh, const std::vector<const FieldBase*>& fields, bool deterministic)
{
  if (mesh.parallel_size() == 1 || fields.empty()) return;

  if (fields[0]->type_is<long double>()) {
    parallel_op_impl<long double, OP>(mesh, fields, deterministic);
  }
  else if (fields[0]->type_is<double>()) {
    parallel_op_impl<double, OP>(mesh, fields, deterministic);
  }
  else if (fields[0]->type_is<float>()) {
    parallel_op_impl<float, OP>(mesh, fields, deterministic);
  }
  else if (fields[0]->type_is<int>()) {
    parallel_op_impl<int, OP>(mesh, fields, deterministic);
  }
  else if (fields[0]->type_is<unsigned long>()) {
    parallel_op_impl<unsigned long, OP>(mesh, fields, deterministic);
  }
  else {
    STK_ThrowRequireMsg(false, "Error, parallel_op only operates on fields of type long double, double, float, int, or unsigned long.");
  }
}

}

void parallel_sum(const BulkData& mesh, const std::vector<const FieldBase*>& fields, bool deterministic)
{
  parallel_op<Operation::SUM>(mesh, fields, deterministic);
}

//----------------------------------------------------------------------

/** Communicate and take the maximum value of field-data for the specified fields
 * on shared entities such that each shared entity
 * will have the same (maximum) field values on each sharing proc.
 */
void parallel_max(const BulkData& mesh, const std::vector<const FieldBase*>& fields)
{
  parallel_op<Operation::MAX>(mesh, fields, false);
}

/** Communicate and take the minimum value of field-data for the specified fields
 * on shared entities such that each shared entity
 * will have the same (minimum) field values on each sharing proc.
 */
void parallel_min(const BulkData& mesh, const std::vector<const FieldBase*>& fields)
{
  parallel_op<Operation::MIN>(mesh, fields, false);
}

template<Operation OP, typename FIELD_DATA_TYPE>
inline void send_or_recv_field_data_for_assembly(stk::CommNeighbors& sparse, int phase, const stk::mesh::FieldBase& f, int owner, const std::vector<int>& procs, unsigned scalars_per_entity, unsigned bucketId, unsigned bucket_ordinal)
{
    FIELD_DATA_TYPE * ptr =
      reinterpret_cast<FIELD_DATA_TYPE *>(stk::mesh::field_data( f , bucketId, bucket_ordinal, scalars_per_entity*sizeof(FIELD_DATA_TYPE) ));

    if (phase == 0) { // send
        CommBufferV & b = sparse.send_buffer( owner );
        for(unsigned i=0; i<scalars_per_entity; ++i)
        {
            b.pack<FIELD_DATA_TYPE>( ptr[i] );
        }
    }
    else { //recv
        DoOp<FIELD_DATA_TYPE, OP> do_op;
        for (int proc : procs) {
          // Only unpack one time from each proc
          CommBufferV & b = sparse.recv_buffer( proc );
          for(unsigned i=0; i<scalars_per_entity; ++i) {
              FIELD_DATA_TYPE recvd_value;
              b.unpack<FIELD_DATA_TYPE>( recvd_value );
              ptr[i] = do_op(ptr[i], recvd_value);
          }
        }
    }
}

template<typename FIELD_DATA_TYPE>
inline void send_or_recv_assembled_data(stk::CommNeighbors& sparse, int phase, const stk::mesh::FieldBase& f, int owner, const std::vector<int>& procs, unsigned scalars_per_entity, unsigned bucketId, unsigned bucket_ordinal)
{
    FIELD_DATA_TYPE * ptr =
      reinterpret_cast<FIELD_DATA_TYPE *>(stk::mesh::field_data( f , bucketId, bucket_ordinal, scalars_per_entity*sizeof(FIELD_DATA_TYPE) ));

    if (phase == 0) { // send
        for (int proc : procs) {
          CommBufferV & b = sparse.send_buffer( proc );
          b.pack<FIELD_DATA_TYPE>( ptr, scalars_per_entity );
        }
    }
    else { //recv
        CommBufferV & b = sparse.recv_buffer( owner );
        b.unpack<FIELD_DATA_TYPE>( ptr, scalars_per_entity );
    }
}

template <Operation OP>
void assemble_to_owner(const stk::mesh::BulkData& mesh,
                       const std::vector<const FieldBase*>& fields,
                       const std::vector<std::pair<int,int>>& fieldRange,
                       const EntityCommListInfoVector& comm_info_vec,
                       CommNeighbors& sparse)
{
  int parallel_rank = mesh.parallel_rank();
  int numFields = fields.size();
  std::vector<int> commProcs;
  for (int phase = 0; phase < 2; ++phase) {

    for ( int fi = 0 ; fi < numFields ; ++fi ) {
      const FieldBase & f = *fields[fi] ;

      unsigned prevBucketId = INVALID_BUCKET_ID;
      unsigned scalarsPerEntity = 0;
      for (int i=fieldRange[fi].first; i<fieldRange[fi].second; ++i) {
        const Entity entity = comm_info_vec[i].entity;
        const bool owned = mesh.parallel_owner_rank(entity) == parallel_rank;

        if ( (!owned && phase == 0) || (owned && phase == 1) )
        {
            const MeshIndex& meshIndex = mesh.mesh_index(entity);
            const unsigned bucketId = meshIndex.bucket->bucket_id();
            if (bucketId != prevBucketId) {
              prevBucketId = bucketId;
              scalarsPerEntity = field_scalars_per_entity(f, bucketId);
            }
            if ( scalarsPerEntity > 0 ) {
              const size_t bucket_ordinal = meshIndex.bucket_ordinal;
              const int owner = mesh.parallel_owner_rank(comm_info_vec[i].entity);
              mesh.comm_procs(comm_info_vec[i].entity, commProcs);

              if (f.data_traits().is_floating_point && f.data_traits().size_of == 8)
              {
                  send_or_recv_field_data_for_assembly<OP, double>(sparse, phase, f, owner, commProcs, scalarsPerEntity, bucketId, bucket_ordinal);
              }
              else if (f.data_traits().is_floating_point && f.data_traits().size_of == 4)
              {
                  send_or_recv_field_data_for_assembly<OP, float>(sparse, phase, f, owner, commProcs, scalarsPerEntity, bucketId, bucket_ordinal);
              }
              else if (f.data_traits().is_integral && f.data_traits().size_of == 4 && f.data_traits().is_unsigned)
              {
                  send_or_recv_field_data_for_assembly<OP, unsigned>(sparse, phase, f, owner, commProcs, scalarsPerEntity, bucketId, bucket_ordinal);
              }
              else if (f.data_traits().is_integral && f.data_traits().size_of == 8 && f.data_traits().is_unsigned)
              {
                  send_or_recv_field_data_for_assembly<OP, unsigned long>(sparse, phase, f, owner, commProcs, scalarsPerEntity, bucketId, bucket_ordinal);
              }
              else if (f.data_traits().is_integral && f.data_traits().size_of == 4 && f.data_traits().is_signed)
              {
                  send_or_recv_field_data_for_assembly<OP, int>(sparse, phase, f, owner, commProcs, scalarsPerEntity, bucketId, bucket_ordinal);
              }
              else
              {
                  STK_ThrowRequireMsg(false,"Unsupported field type in parallel_sum_including_ghosts");
              }
            }
          }
        }
      }

      if (phase == 0) { sparse.communicate(); }
  }
}

void send_back_to_non_owners(const stk::mesh::BulkData& mesh,
                             const std::vector<const FieldBase*>& fields,
                             const std::vector<std::pair<int,int>>& fieldRange,
                             const EntityCommListInfoVector& comm_info_vec,
                             CommNeighbors& sparse)
{
  int parallel_rank = mesh.parallel_rank();
  int numFields = fields.size();
  std::vector<int> commProcs;
  for (int phase = 0; phase < 2; ++phase) {

    for ( int fi = 0 ; fi < numFields ; ++fi ) {
      const FieldBase & f = *fields[fi] ;

      unsigned prevBucketId = INVALID_BUCKET_ID;
      unsigned scalarsPerEntity = 0;
      for (int i=fieldRange[fi].first; i<fieldRange[fi].second; ++i) {
        const Entity entity = comm_info_vec[i].entity;
        const bool owned = mesh.parallel_owner_rank(entity) == parallel_rank;

        if ( (owned && phase == 0) || (!owned && phase == 1) )
        {
            const MeshIndex& meshIndex = mesh.mesh_index(entity);
            const unsigned bucketId = meshIndex.bucket->bucket_id();
            if (bucketId != prevBucketId) {
              prevBucketId = bucketId;
              scalarsPerEntity = field_scalars_per_entity(f, bucketId);
            }
            if ( scalarsPerEntity > 0 ) {
              const size_t bucket_ordinal = meshIndex.bucket_ordinal;
              const int owner = mesh.parallel_owner_rank(comm_info_vec[i].entity);
              mesh.comm_procs(comm_info_vec[i].entity, commProcs);

              if (f.data_traits().is_floating_point && f.data_traits().size_of == 8)
              {
                  send_or_recv_assembled_data<double>(sparse, phase, f, owner, commProcs, scalarsPerEntity, bucketId, bucket_ordinal);
              }
              else if (f.data_traits().is_floating_point && f.data_traits().size_of == 4)
              {
                  send_or_recv_assembled_data<float>(sparse, phase, f, owner, commProcs, scalarsPerEntity, bucketId, bucket_ordinal);
              }
              else if (f.data_traits().is_integral && f.data_traits().size_of == 4 && f.data_traits().is_unsigned)
              {
                  send_or_recv_assembled_data<unsigned>(sparse, phase, f, owner, commProcs, scalarsPerEntity, bucketId, bucket_ordinal);
              }
              else if (f.data_traits().is_integral && f.data_traits().size_of == 8 && f.data_traits().is_unsigned)
              {
                  send_or_recv_assembled_data<unsigned long>(sparse, phase, f, owner, commProcs, scalarsPerEntity, bucketId, bucket_ordinal);
              }
              else if (f.data_traits().is_integral && f.data_traits().size_of == 4 && f.data_traits().is_signed)
              {
                  send_or_recv_assembled_data<int>(sparse, phase, f, owner, commProcs, scalarsPerEntity, bucketId, bucket_ordinal);
              }
              else
              {
                  STK_ThrowRequireMsg(false,"Unsupported field type in parallel_sum_including_ghosts");
              }
            }
          }
        }
      }

      if (phase == 0) { sparse.communicate(); }
  }
}

template <Operation OP>
void parallel_op_including_ghosts_impl(const BulkData & mesh, const std::vector<const FieldBase *> & fields)
{
  if ( fields.empty() ) { return; }

  const int parallel_size = mesh.parallel_size();

  const EntityCommListInfoVector& comm_info_vec = mesh.internal_comm_list();
  const EntityCommDatabase& commDB = mesh.internal_comm_db();
  std::vector<std::pair<int,int> > fieldRange(fields.size(), std::make_pair(-1,-1));
  compute_field_entity_ranges(comm_info_vec, fields, fieldRange);

  // Sizing for send and receive

  const unsigned zero = 0;
  std::vector<unsigned> send_size( parallel_size , zero );
  std::vector<unsigned> recv_size( parallel_size , zero );

  int numFields = fields.size();
  for (int fi = 0 ; fi < numFields ; ++fi ) {
    const FieldBase & f = *fields[fi];
    f.sync_to_host();
    f.modify_on_host();

    for (int i=fieldRange[fi].first; i<fieldRange[fi].second; ++i) {
      STK_ThrowAssertMsg(mesh.is_valid(comm_info_vec[i].entity),"parallel_sum_including_ghosts found invalid entity");
      const MeshIndex& meshIndex = mesh.mesh_index(comm_info_vec[i].entity);
      const Bucket& bucket = *meshIndex.bucket;

      const unsigned bucketId = bucket.bucket_id();
      const unsigned fieldDataSize = field_bytes_per_entity( f , bucketId );

      if (fieldDataSize == 0) {
        continue;
      }

      const bool owned = bucket.owned();

      if ( !owned ) {
         send_size[ mesh.parallel_owner_rank(comm_info_vec[i].entity) ] += fieldDataSize ;
      }
      else {
          PairIterEntityComm info = commDB.comm(comm_info_vec[i].entity_comm);
          for (; !info.empty(); ++info) {
              recv_size[ info->proc ] += fieldDataSize ;
          }
      }
    }
  }

  std::vector<int> send_procs, recv_procs;
  for(int p=0; p<mesh.parallel_size(); ++p) {
      if (send_size[p] > 0) {
          send_procs.push_back(p);
      }
      if (recv_size[p] > 0) {
          recv_procs.push_back(p);
      }
  }

  CommNeighbors sparse(mesh.parallel(), send_procs, recv_procs);

  for(int p=0; p<mesh.parallel_size(); ++p) {
    if (send_size[p] > 0) {
      sparse.send_buffer(p).reserve(send_size[p]);
    }
  }

  assemble_to_owner<OP>(mesh, fields, fieldRange, comm_info_vec, sparse);

  sparse.reset_buffers();
  sparse.swap_send_recv();

  send_back_to_non_owners(mesh, fields, fieldRange, comm_info_vec, sparse);
}

void parallel_sum_including_ghosts(const BulkData & mesh, const std::vector<const FieldBase *> & fields)
{
  parallel_op_including_ghosts_impl<Operation::SUM>(mesh, fields);
}

void parallel_max_including_ghosts(const BulkData & mesh, const std::vector<const FieldBase *> & fields)
{
  parallel_op_including_ghosts_impl<Operation::MAX>(mesh, fields);
}

void parallel_min_including_ghosts(const BulkData & mesh, const std::vector<const FieldBase *> & fields)
{
  parallel_op_including_ghosts_impl<Operation::MIN>(mesh, fields);
}


} // namespace mesh
} // namespace stk
