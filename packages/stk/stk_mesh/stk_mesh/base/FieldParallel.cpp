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
#include <stk_mesh/base/Entity.hpp>     // for Entity
#include <stk_mesh/base/Ghosting.hpp>   // for Ghosting
#include <stk_mesh/base/Types.hpp>      // for PairIterEntityComm, etc

#include <utility>                      // for pair
#include <sstream>                      // for basic_ostream::operator<<, etc
#include <set>

namespace stk {
namespace mesh {


//----------------------------------------------------------------------
void communicate_field_data(
  const Ghosting                        & ghosts ,
  const std::vector< const FieldBase *> & fields )
{
  if ( fields.empty() ) { return; }

  const BulkData & mesh = ghosts.mesh();
  const int parallel_size = mesh.parallel_size();
  const int parallel_rank = mesh.parallel_rank();
  const unsigned ghost_id = ghosts.ordinal();

  for ( const FieldBase* fptr : fields) {
    fptr->sync_to_host();
    fptr->modify_on_host();
  }

  // Sizing for send and receive

  const unsigned zero = 0 ;
  std::vector<unsigned> send_size( parallel_size , zero );
  std::vector<unsigned> recv_size( parallel_size , zero );

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
    const EntityCommInfoVector& infovec = ecli.entity_comm->comm_map;
    if ( owned ) {
      for (const EntityCommInfo& ec : infovec) {
        if (ec.ghost_id == ghost_id) {
          send_size[ ec.proc ] += e_size ;
        }
      }
    }
    else {
      for (const EntityCommInfo& ec : infovec) {
        if (ec.ghost_id == ghost_id) {
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

            const EntityCommInfoVector& infovec = ecli.entity_comm->comm_map;
            if (phase == 0) { // send
              for (const EntityCommInfo& ec : infovec) {
                if (ec.ghost_id == ghost_id) {
                  CommBufferV & b = sparse.send_buffer( ec.proc );
                  b.pack<unsigned char>( ptr , size );
                }
              }
            }
            else { //recv
              for (const EntityCommInfo& ec : infovec) {
                if (ec.ghost_id == ghost_id) {
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
        ThrowRequireMsg(f.type_is<T>(),
                      "Please don't mix fields with different primitive types in the same parallel assemble operation");

        f.sync_to_host();
        const BucketIndices& bktIndices = mesh.volatile_fast_shared_comm_map(f.entity_rank())[proc];
        for(size_t i=0; i<bktIndices.bucket_info.size(); ++i) {
            unsigned bucket = bktIndices.bucket_info[i].bucket_id;
            const int num_bytes_per_entity = field_bytes_per_entity( f , bucket );
            if (num_bytes_per_entity > 0) {
                const int num_Ts_per_field = num_bytes_per_entity / sizeof(T);
                reserve_len += bktIndices.bucket_info[i].num_entities_this_bucket*num_Ts_per_field;
            }
        }
    }
    send_data.reserve(reserve_len);

    for (size_t j = 0 ; j < fields.size() ; ++j ) {
        const FieldBase& f = *fields[j];
        const BucketIndices& bktIndices = mesh.volatile_fast_shared_comm_map(f.entity_rank())[proc];

        const std::vector<unsigned>& ords = bktIndices.ords;
        unsigned offset = 0;

        for(size_t i=0; i<bktIndices.bucket_info.size(); ++i) {
            unsigned bucket = bktIndices.bucket_info[i].bucket_id;
            unsigned num_entities = bktIndices.bucket_info[i].num_entities_this_bucket;
            const int num_bytes_per_entity = field_bytes_per_entity( f , bucket );
            if (num_bytes_per_entity > 0) {
                const int num_Ts_per_field = num_bytes_per_entity / sizeof(T);
                const T* data = reinterpret_cast<T*>(stk::mesh::field_data( f , bucket));

                for (unsigned iord=0; iord<num_entities; ++iord) {
                    const unsigned idx = ords[offset++]*num_Ts_per_field;

                    for (int d = 0; d < num_Ts_per_field; ++d) {
                      send_data.push_back(data[idx+d]);
                    }
                }
            }
            else {
                offset += num_entities;
            }
        }
    }
  };

  auto msgUnpacker = [&fields, &mesh](int iproc, std::vector<T>& recv_data)
  {
    DoOp<T, OP> do_op;

    unsigned offset = 0;
    for (size_t j = 0 ; j < fields.size() ; ++j ) {
        const FieldBase& f = *fields[j] ;
        const BucketIndices& bktIndices = mesh.volatile_fast_shared_comm_map(f.entity_rank())[iproc];
        const std::vector<unsigned>& ords = bktIndices.ords;

        f.sync_to_host();
        f.modify_on_host();

        size_t ords_offset = 0;
        for(size_t i=0; i<bktIndices.bucket_info.size(); ++i) {
            unsigned bucket = bktIndices.bucket_info[i].bucket_id;
            unsigned num_entities = bktIndices.bucket_info[i].num_entities_this_bucket;
            const int num_bytes_per_entity = field_bytes_per_entity( f , bucket );
            if (num_bytes_per_entity > 0) {
                const int num_Ts_per_field = num_bytes_per_entity / sizeof(T);
                T* data = reinterpret_cast<T*>(stk::mesh::field_data( f , bucket));

                for (unsigned iord=0; iord<num_entities; ++iord) {
                    const unsigned idx = ords[ords_offset++]*num_Ts_per_field;

                    for (int d = 0; d < num_Ts_per_field; ++d) {
                      data[idx+d] = do_op(data[idx+d], recv_data[offset + d]);
                    }
                    offset += num_Ts_per_field;
                }
            }
            else {
                ords_offset += num_entities;
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
    ThrowRequireMsg(false, "Error, parallel_op only operates on fields of type long double, double, float, int, or unsigned long.");
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
inline void send_or_recv_field_data_for_assembly(stk::CommNeighbors& sparse, int phase, const stk::mesh::FieldBase& f, int owner, const EntityCommInfoVector& infovec, unsigned scalars_per_entity, unsigned bucketId, unsigned bucket_ordinal)
{
    FIELD_DATA_TYPE * ptr =
      reinterpret_cast<FIELD_DATA_TYPE *>(stk::mesh::field_data( f , bucketId, bucket_ordinal, scalars_per_entity*sizeof(FIELD_DATA_TYPE) ));

    if (phase == 0)
    { // send
        CommBufferV & b = sparse.send_buffer( owner );
        for(unsigned i=0; i<scalars_per_entity; ++i)
        {
            b.pack<FIELD_DATA_TYPE>( ptr[i] );
        }
    }
    else
    { //recv
        DoOp<FIELD_DATA_TYPE, OP> do_op;
        std::set<int> receivedProcs;
        PairIterEntityComm ec(infovec.begin(), infovec.end());
        for ( ; !ec.empty() ; ++ec )
        {
            auto haveNotReceivedFromProc = receivedProcs.insert(ec->proc);
            if (haveNotReceivedFromProc.second) {
                // Only unpack one time from each proc
                CommBufferV & b = sparse.recv_buffer( ec->proc );
                for(unsigned i=0; i<scalars_per_entity; ++i)
                {
                    FIELD_DATA_TYPE recvd_value;
                    b.unpack<FIELD_DATA_TYPE>( recvd_value );
                    ptr[i] = do_op(ptr[i], recvd_value);
                }
            }
        }
    }
}


template <Operation OP>
void parallel_op_including_ghosts_impl(const BulkData & mesh, const std::vector<const FieldBase *> & fields)
{
  if ( fields.empty() ) { return; }

  const int parallel_size = mesh.parallel_size();
  const int parallel_rank = mesh.parallel_rank();

  const std::vector<const FieldBase *>::const_iterator fe = fields.end();
  const std::vector<const FieldBase *>::const_iterator fb = fields.begin();
        std::vector<const FieldBase *>::const_iterator fi ;

  // Sizing for send and receive

  const unsigned zero = 0 ;
  std::vector<int> send_size( parallel_size , zero );
  std::vector<int> recv_size( parallel_size , zero );

  const EntityCommListInfoVector& comm_info_vec = mesh.internal_comm_list();
  size_t comm_info_vec_size = comm_info_vec.size();
  for ( fi = fb ; fi != fe ; ++fi ) {
    const FieldBase & f = **fi ;
    f.sync_to_host();
    f.modify_on_host();

    for (size_t i=0; i<comm_info_vec_size; ++i) {
        if (!mesh.is_valid(comm_info_vec[i].entity))
        {
            ThrowAssertMsg(mesh.is_valid(comm_info_vec[i].entity),"parallel_sum_including_ghosts found invalid entity");
        }
      const MeshIndex& meshIndex = mesh.mesh_index(comm_info_vec[i].entity);
      const Bucket& bucket = *meshIndex.bucket;

      unsigned e_size = 0 ;
      if(is_matching_rank(f, bucket)) {
        const unsigned bucketId = bucket.bucket_id();
        e_size += field_bytes_per_entity( f , bucketId );
      }

      if (e_size == 0) {
        continue;
      }

      const bool owned = bucket.owned();

      if ( !owned ) {
         send_size[ mesh.parallel_owner_rank(comm_info_vec[i].entity) ] += e_size ;
      }
      else {
          const EntityCommInfoVector& infovec = comm_info_vec[i].entity_comm->comm_map;
          size_t info_vec_size = infovec.size();
          for (size_t j=0; j<info_vec_size; ++j ) {
              recv_size[ infovec[j].proc ] += e_size ;
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

    for ( fi = fb ; fi != fe ; ++fi ) {
      const FieldBase & f = **fi ;

      for (size_t i=0; i<comm_info_vec_size; ++i) {
        const bool owned = mesh.parallel_owner_rank(comm_info_vec[i].entity) == parallel_rank;
        if ( (!owned && phase == 0) || (owned && phase == 1) )
        {
            const MeshIndex& meshIndex = mesh.mesh_index(comm_info_vec[i].entity);
            const Bucket& bucket = *meshIndex.bucket;

            if(!is_matching_rank(f, bucket)) continue;

            const unsigned bucketId = bucket.bucket_id();
            const size_t bucket_ordinal = meshIndex.bucket_ordinal;
            const unsigned scalars_per_entity = field_scalars_per_entity(f, bucketId);

            if ( scalars_per_entity > 0 ) {
              const int owner = mesh.parallel_owner_rank(comm_info_vec[i].entity);

              if (f.data_traits().is_floating_point && f.data_traits().size_of == 8)
              {
                  send_or_recv_field_data_for_assembly<OP, double>(sparse, phase, f, owner, comm_info_vec[i].entity_comm->comm_map, scalars_per_entity, bucketId, bucket_ordinal);
              }
              else if (f.data_traits().is_floating_point && f.data_traits().size_of == 4)
              {
                  send_or_recv_field_data_for_assembly<OP, float>(sparse, phase, f, owner, comm_info_vec[i].entity_comm->comm_map, scalars_per_entity, bucketId, bucket_ordinal);
              }
              else if (f.data_traits().is_integral && f.data_traits().size_of == 4 && f.data_traits().is_unsigned)
              {
                  send_or_recv_field_data_for_assembly<OP, unsigned>(sparse, phase, f, owner, comm_info_vec[i].entity_comm->comm_map, scalars_per_entity, bucketId, bucket_ordinal);
              }
              else if (f.data_traits().is_integral && f.data_traits().size_of == 8 && f.data_traits().is_unsigned)
              {
                  send_or_recv_field_data_for_assembly<OP, unsigned long>(sparse, phase, f, owner, comm_info_vec[i].entity_comm->comm_map, scalars_per_entity, bucketId, bucket_ordinal);
              }
              else if (f.data_traits().is_integral && f.data_traits().size_of == 4 && f.data_traits().is_signed)
              {
                  send_or_recv_field_data_for_assembly<OP, int>(sparse, phase, f, owner, comm_info_vec[i].entity_comm->comm_map, scalars_per_entity, bucketId, bucket_ordinal);
              }
              else
              {
                  ThrowRequireMsg(false,"Unsupported field type in parallel_sum_including_ghosts");
              }
            }
          }
        }
      }

      if (phase == 0) { sparse.communicate(); }
  }

  communicate_field_data(mesh, fields);
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
