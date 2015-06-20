// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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
#include <stk_util/parallel/ParallelComm.hpp>  // for CommAll, CommBuffer
#include <stk_util/parallel/ParallelReduce.hpp>  // for Reduce, ReduceSum, etc
#include <stk_util/parallel/Parallel.hpp>  // for parallel_machine_rank, etc
#include <stk_util/util/PairIter.hpp>   // for PairIter

#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/base/Entity.hpp>     // for Entity
#include <stk_mesh/base/Ghosting.hpp>   // for Ghosting
#include <stk_mesh/base/Types.hpp>      // for PairIterEntityComm, etc

#include <utility>                      // for pair
#include <sstream>                      // for basic_ostream::operator<<, etc

namespace stk {
namespace mesh {

void communicate_field_data(
  const BulkData                        & mesh ,
  const std::vector< const FieldBase *> & fields )
{
  if ( fields.empty() ) { return; }

  const int parallel_size = mesh.parallel_size();
  const int parallel_rank = mesh.parallel_rank();

  const std::vector<const FieldBase *>::const_iterator fe = fields.end();
  const std::vector<const FieldBase *>::const_iterator fb = fields.begin();
        std::vector<const FieldBase *>::const_iterator fi ;

  // Sizing for send and receive

  const unsigned zero = 0 ;
  std::vector<unsigned> send_size( parallel_size , zero );
  std::vector<unsigned> recv_size( parallel_size , zero );

  std::vector<int> visited_procs;
  visited_procs.reserve(128);
  for ( EntityCommListInfoVector::const_iterator
        i =  mesh.internal_comm_list().begin() , iend = mesh.internal_comm_list().end(); i != iend ; ++i ) {
    Entity e = i->entity;
    const MeshIndex meshIdx = mesh.mesh_index(e);
    const unsigned bucketId = meshIdx.bucket->bucket_id();

    const bool owned = i->owner == parallel_rank ;

    unsigned e_size = 0 ;
    for ( fi = fb ; fi != fe ; ++fi ) {
      const FieldBase & f = **fi ;

      if(is_matching_rank(f, *meshIdx.bucket)) {
        e_size += field_bytes_per_entity( f , bucketId );
      }
    }

    if (e_size == 0) {
      continue;
    }

    if ( owned ) {
      const EntityCommInfoVector& infovec = i->entity_comm->comm_map;
      PairIterEntityComm ec(infovec.begin(), infovec.end());
      visited_procs.clear();
      for ( ; ! ec.empty() ; ++ec ) {
          if (std::find(visited_procs.begin(),visited_procs.end(),ec->proc)==visited_procs.end()) {
              send_size[ ec->proc ] += e_size ;
              visited_procs.push_back(ec->proc);
          }
      }
    }
    else {
      recv_size[ i->owner ] += e_size ;
    }
  }

  // Allocate send and receive buffers:

  CommAll sparse ;

  {
    const unsigned * const snd_size = & send_size[0] ;
    const unsigned * const rcv_size = & recv_size[0] ;
    sparse.allocate_buffers( mesh.parallel(), snd_size, rcv_size);
  }

  // Send packing:

  for (int phase = 0; phase < 2; ++phase) {

    for ( EntityCommListInfoVector::const_iterator i =  mesh.internal_comm_list().begin(), iend = mesh.internal_comm_list().end() ; i != iend ; ++i ) {
      if ( (i->owner == parallel_rank && phase == 0) ||
           (i->owner != parallel_rank && phase == 1) ) {
        Entity e = i->entity;
        const MeshIndex meshIdx = mesh.mesh_index(e);
        const unsigned bucketId = meshIdx.bucket->bucket_id();

        for ( fi = fb ; fi != fe ; ++fi ) {
          const FieldBase & f = **fi ;

          if(!is_matching_rank(f, *meshIdx.bucket)) continue;

          const unsigned size = field_bytes_per_entity( f , bucketId );

          if ( size ) {
            unsigned char * ptr =
              reinterpret_cast<unsigned char *>(stk::mesh::field_data( f , bucketId, meshIdx.bucket_ordinal, size ));

            if (phase == 0) { // send
              const EntityCommInfoVector& infovec = i->entity_comm->comm_map;
              PairIterEntityComm ec(infovec.begin(), infovec.end());
              visited_procs.clear();
              for ( ; !ec.empty() ; ++ec ) {
                  if (std::find(visited_procs.begin(),visited_procs.end(),ec->proc)==visited_procs.end()) {
                      CommBuffer & b = sparse.send_buffer( ec->proc );
                      b.pack<unsigned char>( ptr , size );
                      visited_procs.push_back(ec->proc);
                  }
              }
            }
            else { //recv
              CommBuffer & b = sparse.recv_buffer( i->owner );
              b.unpack<unsigned char>( ptr , size );
            }
          }
        }
      }
    }
    if (phase == 0) { sparse.communicate(); }
  }
}

//----------------------------------------------------------------------
void communicate_field_data(
  const Ghosting                        & ghosts ,
  const std::vector< const FieldBase *> & fields )
{
  if ( fields.empty() ) { return; }

  const BulkData & mesh = BulkData::get(ghosts);
  const int parallel_size = mesh.parallel_size();
  const int parallel_rank = mesh.parallel_rank();
  const unsigned ghost_id = ghosts.ordinal();

  const std::vector<const FieldBase *>::const_iterator fe = fields.end();
  const std::vector<const FieldBase *>::const_iterator fb = fields.begin();
        std::vector<const FieldBase *>::const_iterator fi ;

  // Sizing for send and receive

  const unsigned zero = 0 ;
  std::vector<unsigned> send_size( parallel_size , zero );
  std::vector<unsigned> recv_size( parallel_size , zero );

  for ( EntityCommListInfoVector::const_iterator
        i =  mesh.internal_comm_list().begin() , iend = mesh.internal_comm_list().end(); i != iend ; ++i ) {
    Entity e = i->entity;
    const MeshIndex meshIdx = mesh.mesh_index(e);
    const unsigned bucketId = meshIdx.bucket->bucket_id();

    const bool owned = i->owner == parallel_rank ;

    unsigned e_size = 0 ;
    for ( fi = fb ; fi != fe ; ++fi ) {
      const FieldBase & f = **fi ;

      if(is_matching_rank(f, *meshIdx.bucket)) {
        e_size += field_bytes_per_entity( f , bucketId );
      }
    }

    if (e_size == 0) {
      continue;
    }

    const EntityCommInfoVector& infovec = i->entity_comm->comm_map;
    PairIterEntityComm ec(infovec.begin(), infovec.end());
    if ( owned ) {
      for ( ; ! ec.empty() ; ++ec ) {
        if (ec->ghost_id == ghost_id) {
          send_size[ ec->proc ] += e_size ;
        }
      }
    }
    else {
      for ( ; ! ec.empty() ; ++ec ) {
        if (ec->ghost_id == ghost_id) {
          recv_size[ i->owner ] += e_size ;
          break;//jump out since we know we're only recving 1 msg from the 1-and-only owner
        }
      }
    }
  }

  // Allocate send and receive buffers:

  CommAll sparse ;

  {
    const unsigned * const snd_size = & send_size[0] ;
    const unsigned * const rcv_size = & recv_size[0] ;
    sparse.allocate_buffers( mesh.parallel(), snd_size, rcv_size);
  }

  // Send packing:

  for (int phase = 0; phase < 2; ++phase) {

    for ( EntityCommListInfoVector::const_iterator i =  mesh.internal_comm_list().begin(), iend = mesh.internal_comm_list().end() ; i != iend ; ++i ) {
      if ( (i->owner == parallel_rank && phase == 0) || (i->owner != parallel_rank && phase == 1) ) {
        Entity e = i->entity;
        const MeshIndex meshIdx = mesh.mesh_index(e);
        const unsigned bucketId = meshIdx.bucket->bucket_id();

        for ( fi = fb ; fi != fe ; ++fi ) {
          const FieldBase & f = **fi ;

          if(!is_matching_rank(f, e)) continue;

          const unsigned size = field_bytes_per_entity( f , e );

          if ( size ) {
            unsigned char * ptr =
              reinterpret_cast<unsigned char *>(stk::mesh::field_data( f , bucketId, meshIdx.bucket_ordinal, size ));

            const EntityCommInfoVector& infovec = i->entity_comm->comm_map;
            PairIterEntityComm ec(infovec.begin(), infovec.end());
            if (phase == 0) { // send
              for ( ; !ec.empty() ; ++ec ) {
                if (ec->ghost_id == ghost_id) {
                  CommBuffer & b = sparse.send_buffer( ec->proc );
                  b.pack<unsigned char>( ptr , size );
                }
              }
            }
            else { //recv
              for ( ; !ec.empty(); ++ec ) {
                if (ec->ghost_id == ghost_id) {
                  CommBuffer & b = sparse.recv_buffer( i->owner );
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

enum Operation
{
  SUM,
  MIN,
  MAX
};

template <typename T, Operation OP>
struct DoOp;

template <typename T>
struct DoOp<T, SUM>
{
  T operator()(T lhs, T rhs) const
  { return lhs + rhs; }
};

template <typename T>
struct DoOp<T, MIN>
{
  T operator()(T lhs, T rhs) const
  { return lhs < rhs ? lhs : rhs; }
};

template <typename T>
struct DoOp<T, MAX>
{
  T operator()(T lhs, T rhs) const
  { return lhs > rhs ? lhs : rhs; }
};


template <typename T, Operation OP>
void parallel_op_impl(const BulkData& mesh, std::vector<FieldBase*> fields)
{
  const int parallel_size = mesh.parallel_size();

  std::vector<std::vector<T> > send_data(parallel_size);
  std::vector<std::vector<T> > recv_data(parallel_size);

  for (size_t j = 0 ; j < fields.size() ; ++j ) {
    const FieldBase& f = *fields[j];
    ThrowRequireMsg(f.type_is<T>(),
                    "Please don't mix fields with different primitive types in the same parallel assemble operation");

    VolatileFastSharedCommMapOneRank const& fast_comm_map = mesh.volatile_fast_shared_comm_map(f.entity_rank());

    for (int iproc=0; iproc<parallel_size; ++iproc) {
      // Not enough for multidimensional fields, but better than nothing
      send_data[iproc].reserve(fast_comm_map[iproc].size());

      for (size_t idata=0, idata_end = fast_comm_map[iproc].size(); idata < idata_end; ++idata) {
        unsigned const bucket = fast_comm_map[iproc][idata].bucket_id;
        unsigned const ord    = fast_comm_map[iproc][idata].bucket_ord;

        const int num_bytes_per_field = field_bytes_per_entity( f , bucket );
        const int num_Ts_per_field = num_bytes_per_field / sizeof(T);
        if (num_Ts_per_field > 0) {
          T* data = reinterpret_cast<T*>(stk::mesh::field_data( f , bucket, ord, num_bytes_per_field ));
          for (int d = 0; d < num_Ts_per_field; ++d) {
            send_data[iproc].push_back(data[d]);
          }
        }
      }
    }
  }

  MPI_Comm comm = mesh.parallel();
  parallel_data_exchange_sym_t(send_data, recv_data, comm);

  DoOp<T, OP> do_op;

  std::vector<unsigned> offset(parallel_size, 0);
  for (size_t j = 0 ; j < fields.size() ; ++j ) {
    const FieldBase& f = *fields[j] ;
    stk::mesh::VolatileFastSharedCommMapOneRank const& fast_comm_map = mesh.volatile_fast_shared_comm_map(f.entity_rank());

    for (int iproc=0; iproc<parallel_size; ++iproc) {

      for (size_t idata=0, idata_end = fast_comm_map[iproc].size(); idata < idata_end; ++idata) {
        unsigned const bucket = fast_comm_map[iproc][idata].bucket_id;
        unsigned const ord    = fast_comm_map[iproc][idata].bucket_ord;

        const int num_bytes_per_field = field_bytes_per_entity( f , bucket );
        const int num_Ts_per_field = num_bytes_per_field / sizeof(T);
        if (num_Ts_per_field > 0) {
          T* data = reinterpret_cast<T*>(stk::mesh::field_data( f , bucket, ord, num_bytes_per_field ));
          for (int d = 0; d < num_Ts_per_field; ++d) {
            data[d] = do_op(data[d], recv_data[iproc][offset[iproc] + d]);
          }
          offset[iproc] += num_Ts_per_field;
        }
      }
    }
  }
}

template <Operation OP>
inline
void parallel_op(const BulkData& mesh, const std::vector<FieldBase*>& fields)
{
  if (mesh.parallel_size() == 1 || fields.empty()) return;

  if (fields[0]->type_is<double>()) {
    parallel_op_impl<double, OP>(mesh, fields);
  }
  else if (fields[0]->type_is<float>()) {
    parallel_op_impl<float, OP>(mesh, fields);
  }
  else if (fields[0]->type_is<int>()) {
    parallel_op_impl<int, OP>(mesh, fields);
  }
  else {
    ThrowRequireMsg(false, "Error, parallel_max only operates on fields of type double, float or int.");
  }
}

}

void parallel_sum(const BulkData& mesh, const std::vector<FieldBase*>& fields)
{
  parallel_op<SUM>(mesh, fields);
}

//----------------------------------------------------------------------

/** Communicate and take the maximum value of field-data for the specified fields
 * on shared entities such that each shared entity
 * will have the same (maximum) field values on each sharing proc.
 */
void parallel_max(const BulkData& mesh, const std::vector<FieldBase*>& fields)
{
  parallel_op<MAX>(mesh, fields);
}

/** Communicate and take the minimum value of field-data for the specified fields
 * on shared entities such that each shared entity
 * will have the same (minimum) field values on each sharing proc.
 */
void parallel_min(const BulkData& mesh, const std::vector<FieldBase*>& fields)
{
  parallel_op<MIN>(mesh, fields);
}

//
//  Determine the number of items each other process will send to the current processor
//
std::vector<int> ComputeReceiveList(std::vector<int>& sendSizeArray, MPI_Comm &mpi_communicator) {
  const int msg_tag = 10240;
  int num_procs = sendSizeArray.size();
  int my_proc;
  MPI_Comm_rank(mpi_communicator, &my_proc);
  std::vector<int> receiveSizeArray(num_procs, 0);
  //
  //  Determine the total number of messages every processor will receive
  //
  std::vector<int> local_number_to_receive(num_procs, 0);
  std::vector<int> global_number_to_receive(num_procs, 0);
  for(int iproc = 0; iproc < num_procs; ++iproc) {
    if(sendSizeArray[iproc] > 0) local_number_to_receive[iproc] = 1;
  }
  MPI_Allreduce(&local_number_to_receive[0], &global_number_to_receive[0], num_procs, MPI_INT, MPI_SUM, mpi_communicator);
  MPI_Barrier(mpi_communicator);
  //
  //  Now each processor knows how many messages it will recive, but does not know the message lengths or where
  //  the messages will be recived from.  Need to extract this information.
  //  Post a recieve for each expected message.
  //
  std::vector<MPI_Request> recv_handles(num_procs);
  int num_to_recv = global_number_to_receive[my_proc];
  std::vector<int> recv_size_buffers(num_to_recv);
  for(int imsg = 0; imsg < num_to_recv; ++imsg) {
    int *recv_buffer = &(recv_size_buffers[imsg]);
    MPI_Irecv(recv_buffer, 1, MPI_INT, MPI_ANY_SOURCE,
              msg_tag, mpi_communicator, &recv_handles[imsg]);
  }
  MPI_Barrier(mpi_communicator);
  //
  //  Send message lengths
  //
  for(int iproc = 0; iproc < num_procs; ++iproc) {
    if(sendSizeArray[iproc] > 0) {
      int send_length = sendSizeArray[iproc];
      MPI_Send(&send_length, 1, MPI_INT, iproc, msg_tag, mpi_communicator);
    }
  }
  //
  //  Get each message and place the length in the proper place in the length array
  //
  for(int imsg = 0; imsg < num_to_recv; ++imsg) {
    MPI_Status status;
    MPI_Wait(&recv_handles[imsg], &status);
    receiveSizeArray[status.MPI_SOURCE] = recv_size_buffers[imsg];
  }

  return receiveSizeArray;
}
} // namespace mesh
} // namespace stk
