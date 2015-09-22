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

#ifndef stk_mesh_FieldParallel_hpp
#define stk_mesh_FieldParallel_hpp

#include <stk_util/stk_config.h>
#include <stk_mesh/base/Types.hpp>      // for EntityProc
#include <stk_mesh/base/FieldTraits.hpp>  // for FieldTraits
#include <stk_mesh/base/FieldBase.hpp>  // for FieldBase
#include <stk_mesh/base/BulkData.hpp>
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine
#include <stk_util/parallel/ParallelComm.hpp>  // for CommAll
#include <stk_util/environment/ReportHandler.hpp>  // for ThrowRequireMsg

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
  communicate_field_data(*mesh.ghostings()[0], fields);
}

//----------------------------------------------------------------------

/** Sum/Max/Min (assemble) field-data for the specified fields on shared entities such that each shared entity
 * will have the same field values on each sharing proc.
 */
void parallel_sum(const BulkData& mesh, const std::vector<FieldBase*>& fields);
void parallel_max(const BulkData& mesh, const std::vector<FieldBase*>& fields);
void parallel_min(const BulkData& mesh, const std::vector<FieldBase*>& fields);


//
//  Generalized comm plans
//
//  This plan assumes the send and recv lists have identical sizes so no extra sizing communications are needed
//
template<typename T>
void parallel_data_exchange_sym_t(std::vector< std::vector<T> > &send_lists,
                                  std::vector< std::vector<T> > &recv_lists,
                                  MPI_Comm &mpi_communicator )
{
  //
  //  Determine the number of processors involved in this communication
  //
#if defined( STK_HAS_MPI)
  const int msg_tag = 10242;
  int num_procs = stk::parallel_machine_size(mpi_communicator);
  int class_size = sizeof(T);

  //
  //  Send the actual messages as raw byte streams.
  //
  std::vector<MPI_Request> recv_handles(num_procs);
  for(int iproc = 0; iproc < num_procs; ++iproc) {
    recv_lists[iproc].resize(send_lists[iproc].size());
    if(recv_lists[iproc].size() > 0) {
      char* recv_buffer = (char*)&recv_lists[iproc][0];
      int recv_size = recv_lists[iproc].size()*class_size;
      MPI_Irecv(recv_buffer, recv_size, MPI_CHAR,
                iproc, msg_tag, mpi_communicator, &recv_handles[iproc]);
    }
  }
  MPI_Barrier(mpi_communicator);
  for(int iproc = 0; iproc < num_procs; ++iproc) {
    if(send_lists[iproc].size() > 0) {
      char* send_buffer = (char*)&send_lists[iproc][0];
      int send_size = send_lists[iproc].size()*class_size;
      MPI_Send(send_buffer, send_size, MPI_CHAR,
               iproc, msg_tag, mpi_communicator);
    }
  }
  for(int iproc = 0; iproc < num_procs; ++iproc) {
    if(recv_lists[iproc].size() > 0) {
      MPI_Status status;
      MPI_Wait( &recv_handles[iproc], &status );
    }
  }
#endif
}

template<typename T>
inline void parallel_data_exchange_nonsym_known_sizes_t(std::vector< std::vector<T> > &send_lists,
                                                std::vector< std::vector<T> > &recv_lists,
                                                MPI_Comm mpi_communicator )
{
#if defined( STK_HAS_MPI)
  const int msg_tag = 10243; //arbitrary tag value, anything less than 32768 is legal
  int num_procs = stk::parallel_machine_size(mpi_communicator);
  int class_size = sizeof(T);

  //
  //  Send the actual messages as raw byte streams.
  //
  std::vector<MPI_Request> recv_handles(num_procs);
  for(int iproc = 0; iproc < num_procs; ++iproc) {
    if(recv_lists[iproc].size() > 0) {
      char* recv_buffer = (char*)&recv_lists[iproc][0];
      int recv_size = recv_lists[iproc].size()*class_size;
      MPI_Irecv(recv_buffer, recv_size, MPI_CHAR, iproc, msg_tag, mpi_communicator, &recv_handles[iproc]);
    }
  }

  MPI_Barrier(mpi_communicator);

  for(int iproc = 0; iproc < num_procs; ++iproc) {
    if(send_lists[iproc].size() > 0) {
      char* send_buffer = (char*)&send_lists[iproc][0];
      int send_size = send_lists[iproc].size()*class_size;
      MPI_Send(send_buffer, send_size, MPI_CHAR, iproc, msg_tag, mpi_communicator);
    }
  }

  for(int iproc = 0; iproc < num_procs; ++iproc) {
    if(recv_lists[iproc].size() > 0) {
      MPI_Status status;
      MPI_Wait( &recv_handles[iproc], &status );
    }
  }
#endif
}

inline bool find_proc_before_index(const EntityCommInfoVector& infovec, int proc, int index)
{
    for(int i=0; i<index; ++i) {
        if (proc == infovec[i].proc) {
            return true;
        }
    }
    return false;
}

inline void communicate_field_data(
  const BulkData                        & mesh ,
  const std::vector< const FieldBase *> & fields )
{
  if ( fields.empty() ) { return; }

  const int parallel_size = mesh.parallel_size();
  const int parallel_rank = mesh.parallel_rank();

  const std::vector<const FieldBase *>::const_iterator fe = fields.end();
  const std::vector<const FieldBase *>::const_iterator fb = fields.begin();
        std::vector<const FieldBase *>::const_iterator fi ;

  std::vector<std::vector<unsigned char> > send_data(parallel_size);
  std::vector<std::vector<unsigned char> > recv_data(parallel_size);

  const EntityCommListInfoVector &comm_info_vec = mesh.internal_comm_list();
  size_t comm_info_vec_size = comm_info_vec.size();

  std::vector<unsigned> send_sizes(parallel_size, 0);
  std::vector<unsigned> recv_sizes(parallel_size, 0);

  //this first loop calculates send_sizes and recv_sizes.
  for(fi = fb; fi != fe; ++fi)
  {
      const FieldBase & f = **fi;
      for(size_t i = 0; i<comm_info_vec_size; ++i)
      {
          const Bucket* bucket = comm_info_vec[i].bucket;

          int owner = comm_info_vec[i].owner;
          const bool owned = (owner == parallel_rank);

          unsigned e_size = 0;

          if(is_matching_rank(f, *bucket))
          {
              const unsigned bucketId = bucket->bucket_id();
              unsigned size = field_bytes_per_entity(f, bucketId);
              e_size += size;
          }

          if(e_size == 0)
          {
              continue;
          }

          if(owned)
          {
              const EntityCommInfoVector& infovec = comm_info_vec[i].entity_comm->comm_map;
              size_t infovec_size = infovec.size();
              for(size_t j=0; j<infovec_size; ++j)
              {
                  int proc = infovec[j].proc;

                  bool proc_already_found = find_proc_before_index(infovec, proc, j); 
                  if (!proc_already_found) {
                      send_sizes[proc] += e_size;
                  }
              }
          }
          else
          {
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

  for(fi = fb; fi != fe; ++fi)
  {
      const FieldBase & f = **fi;
      for(size_t i = 0; i<comm_info_vec_size; ++i)
      {
          const Bucket* bucket = comm_info_vec[i].bucket;

          int owner = comm_info_vec[i].owner;
          const bool owned = (owner == parallel_rank);

          unsigned e_size = 0;

          if(is_matching_rank(f, *bucket))
          {
              const unsigned bucketId = bucket->bucket_id();
              unsigned size = field_bytes_per_entity(f, bucketId);
              if (owned && size > 0)
              {
                  unsigned char * ptr = reinterpret_cast<unsigned char*>(stk::mesh::field_data(f, bucketId, comm_info_vec[i].bucket_ordinal, size));
                  std::memcpy(field_data_ptr+e_size, ptr, size);
 //                 field_data.insert(field_data.end(), ptr, ptr+size);
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
              size_t infovec_size = infovec.size();
              for(size_t j=0; j<infovec_size; ++j)
              {
                  int proc = infovec[j].proc;
    
                  bool proc_already_found = find_proc_before_index(infovec, proc, j);
                  if (!proc_already_found) {
                      unsigned char* dest_ptr = send_data[proc].data()+send_sizes[proc];
                      unsigned char* src_ptr = field_data_ptr;
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
  for(fi = fb; fi != fe; ++fi)
  {
      const FieldBase & f = **fi;

      for(size_t i=0; i<comm_info_vec_size; ++i)
      {
          int owner = comm_info_vec[i].owner;
          const bool owned = (owner == parallel_rank);

          if(owned || recv_data[owner].size() == 0)
          {
              continue;
          }

          const Bucket* bucket = comm_info_vec[i].bucket;

          if(is_matching_rank(f, *bucket))
          {
              const unsigned bucketId = bucket->bucket_id();
              unsigned size = field_bytes_per_entity(f, bucketId);
              if (size > 0)
              {
                  unsigned char * ptr = reinterpret_cast<unsigned char*>(stk::mesh::field_data(f, bucketId, comm_info_vec[i].bucket_ordinal, size));

                  std::memcpy(ptr, &(recv_data[owner][recv_sizes[owner]]), size);
                  recv_sizes[owner] += size;
              }
          }
      }
  }
}

template<typename FIELD_DATA_TYPE>
inline void send_or_recv_field_data_for_assembly(stk::CommAll& sparse, int phase, const stk::mesh::FieldBase& f, int owner, const EntityCommInfoVector& infovec, unsigned scalars_per_entity, unsigned bucketId, unsigned bucket_ordinal)
{
    FIELD_DATA_TYPE * ptr =
      reinterpret_cast<FIELD_DATA_TYPE *>(stk::mesh::field_data( f , bucketId, bucket_ordinal, scalars_per_entity*sizeof(FIELD_DATA_TYPE) ));

    if (phase == 0)
    { // send
        CommBuffer & b = sparse.send_buffer( owner );
        for(unsigned i=0; i<scalars_per_entity; ++i)
        {
          b.pack<FIELD_DATA_TYPE>( ptr[i] );
        }
    }
    else
    { //recv
        PairIterEntityComm ec(infovec.begin(), infovec.end());
        for ( ; !ec.empty() ; ++ec )
        {
            CommBuffer & b = sparse.recv_buffer( ec->proc );
            for(unsigned i=0; i<scalars_per_entity; ++i)
            {
                FIELD_DATA_TYPE recvd_value;
                b.unpack<FIELD_DATA_TYPE>( recvd_value );
                ptr[i] += recvd_value;
            }
        }
    }
}

inline void parallel_sum_including_ghosts(
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

  const EntityCommListInfoVector& comm_info_vec = mesh.internal_comm_list();
  size_t comm_info_vec_size = comm_info_vec.size();
  for ( fi = fb ; fi != fe ; ++fi ) {
    const FieldBase & f = **fi ;

    for (size_t i=0; i<comm_info_vec_size; ++i) {
        if (!mesh.is_valid(comm_info_vec[i].entity))
        {
            ThrowAssertMsg(mesh.is_valid(comm_info_vec[i].entity),"parallel_sum_including_ghosts found invalid entity");
        }
      const Bucket* bucket = comm_info_vec[i].bucket;

      unsigned e_size = 0 ;
      if(is_matching_rank(f, *bucket)) {
        const unsigned bucketId = bucket->bucket_id();
        e_size += field_bytes_per_entity( f , bucketId );
      }

      if (e_size == 0) {
        continue;
      }

      const bool owned = comm_info_vec[i].owner == parallel_rank ;

      if ( !owned ) {
         send_size[ comm_info_vec[i].owner ] += e_size ;
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

  // Allocate send and receive buffers:

  CommAll sparse ;

  {
    const unsigned * const snd_size = & send_size[0] ;
    const unsigned * const rcv_size = & recv_size[0] ;
    sparse.allocate_buffers( mesh.parallel(), snd_size, rcv_size);
  }

  // Send packing:

  for (int phase = 0; phase < 2; ++phase) {

    for ( fi = fb ; fi != fe ; ++fi ) {
      const FieldBase & f = **fi ;

      for (size_t i=0; i<comm_info_vec_size; ++i) {
        const bool owned = comm_info_vec[i].owner == parallel_rank;
        if ( (!owned && phase == 0) || (owned && phase == 1) )
        {
            const Bucket* bucket = comm_info_vec[i].bucket;

            if(!is_matching_rank(f, *bucket)) continue;

            const unsigned bucketId = bucket->bucket_id();
            const size_t bucket_ordinal = comm_info_vec[i].bucket_ordinal;
            const unsigned scalars_per_entity = field_scalars_per_entity(f, bucketId);

            if ( scalars_per_entity > 0 ) {
              int owner = comm_info_vec[i].owner;

              if (f.data_traits().is_floating_point && f.data_traits().size_of == 8)
              {
                  send_or_recv_field_data_for_assembly<double>(sparse, phase, f, owner, comm_info_vec[i].entity_comm->comm_map, scalars_per_entity, bucketId, bucket_ordinal);
              }
              else if (f.data_traits().is_floating_point && f.data_traits().size_of == 4)
              {
                  send_or_recv_field_data_for_assembly<float>(sparse, phase, f, owner, comm_info_vec[i].entity_comm->comm_map, scalars_per_entity, bucketId, bucket_ordinal);
              }
              else if (f.data_traits().is_integral && f.data_traits().size_of == 4 && f.data_traits().is_unsigned)
              {
                  send_or_recv_field_data_for_assembly<unsigned>(sparse, phase, f, owner, comm_info_vec[i].entity_comm->comm_map, scalars_per_entity, bucketId, bucket_ordinal);
              }
              else if (f.data_traits().is_integral && f.data_traits().size_of == 4 && f.data_traits().is_signed)
              {
                  send_or_recv_field_data_for_assembly<int>(sparse, phase, f, owner, comm_info_vec[i].entity_comm->comm_map, scalars_per_entity, bucketId, bucket_ordinal);
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

//
//  This plan assumes the send and recv lists are matched, but that the actual ammount of data to send is unknown.
//  A processor knows which other processors it will be receiving data from, but does not know who much data.
//  Thus the comm plan is known from the inputs, but an additional message sizing call must be done.
//
template<typename T>
void parallel_data_exchange_sym_unknown_size_t(std::vector< std::vector<T> > &send_lists,
                                               std::vector< std::vector<T> > &recv_lists,
                                               MPI_Comm &mpi_communicator )
{
#if defined( STK_HAS_MPI)
  const int msg_tag = 10242;
  int num_procs = stk::parallel_machine_size(mpi_communicator);
  int class_size = sizeof(T);

  //
  //  Send the message sizes
  //
  std::vector<int> send_msg_sizes(num_procs);
  std::vector<int> recv_msg_sizes(num_procs);
  std::vector<MPI_Request> recv_handles(num_procs);

  for(int iproc = 0; iproc < num_procs; ++iproc) {
    send_msg_sizes[iproc] = send_lists[iproc].size();
  }    
  for(int iproc = 0; iproc < num_procs; ++iproc) {
    if(recv_lists[iproc].size()>0) {
      MPI_Irecv(&recv_msg_sizes[iproc], 1, MPI_INT, iproc, msg_tag, mpi_communicator, &recv_handles[iproc]);
    }
  }
  MPI_Barrier(mpi_communicator);
  for(int iproc = 0; iproc < num_procs; ++iproc) {
    if(send_lists[iproc].size()>0) {
      MPI_Send(&send_msg_sizes[iproc], 1, MPI_INT, iproc, msg_tag, mpi_communicator);
    }
  }
  for(int iproc = 0; iproc < num_procs; ++iproc) {
    if(recv_lists[iproc].size() > 0) {
      MPI_Status status;
      MPI_Wait( &recv_handles[iproc], &status );
      recv_lists[iproc].resize(recv_msg_sizes[iproc]);
    }
  }
  //
  //  Send the actual messages as raw byte streams.
  //
  for(int iproc = 0; iproc < num_procs; ++iproc) {
    if(recv_lists[iproc].size() > 0) {
      char* recv_buffer = (char*)&recv_lists[iproc][0];
      int recv_size = recv_lists[iproc].size()*class_size;
      MPI_Irecv(recv_buffer, recv_size, MPI_CHAR,
                iproc, msg_tag, mpi_communicator, &recv_handles[iproc]);
    }
  }
  MPI_Barrier(mpi_communicator);
  for(int iproc = 0; iproc < num_procs; ++iproc) {
    if(send_lists[iproc].size() > 0) {
      char* send_buffer = (char*)&send_lists[iproc][0];
      int send_size = send_lists[iproc].size()*class_size;
      MPI_Send(send_buffer, send_size, MPI_CHAR,
               iproc, msg_tag, mpi_communicator);
    }
  }
  for(int iproc = 0; iproc < num_procs; ++iproc) {
    if(recv_lists[iproc].size() > 0) {
      MPI_Status status;
      MPI_Wait( &recv_handles[iproc], &status );
    }
  }
#endif
}

std::vector<int> ComputeReceiveList(std::vector<int>& sendSizeArray, MPI_Comm &mpi_communicator);

//
//  Parallel_Data_Exchange: General object exchange template with unknown comm plan
//
template<typename T>
void parallel_data_exchange_t(std::vector< std::vector<T> > &send_lists,
                              std::vector< std::vector<T> > &recv_lists,
                              MPI_Comm &mpi_communicator ) {
  //
  //  Determine the number of processors involved in this communication
  //
  const int msg_tag = 10242;
  int num_procs;
  MPI_Comm_size(mpi_communicator, &num_procs);
  int my_proc;
  MPI_Comm_rank(mpi_communicator, &my_proc);
  ThrowRequire((unsigned int) num_procs == send_lists.size() && (unsigned int) num_procs == recv_lists.size());
  int class_size = sizeof(T);
  //
  //  Determine number of items each other processor will send to the current processor
  //
  std::vector<int> global_number_to_send(num_procs);
  for(int iproc=0; iproc<num_procs; ++iproc) {
    global_number_to_send[iproc] = send_lists[iproc].size();
  }
  std::vector<int> numToRecvFrom = ComputeReceiveList(global_number_to_send, mpi_communicator);
  //
  //  Send the actual messages as raw byte streams.
  //
  std::vector<MPI_Request> recv_handles(num_procs);
  for(int iproc = 0; iproc < num_procs; ++iproc) {
    recv_lists[iproc].resize(numToRecvFrom[iproc]);
    if(recv_lists[iproc].size() > 0) {
      char* recv_buffer = (char*)&recv_lists[iproc][0];
      int recv_size = recv_lists[iproc].size()*class_size;
      MPI_Irecv(recv_buffer, recv_size, MPI_CHAR,
                iproc, msg_tag, mpi_communicator, &recv_handles[iproc]);
    }
  }
  MPI_Barrier(mpi_communicator);
  for(int iproc = 0; iproc < num_procs; ++iproc) {
    if(send_lists[iproc].size() > 0) {
      char* send_buffer = (char*)&send_lists[iproc][0];
      int send_size = send_lists[iproc].size()*class_size;
      MPI_Send(send_buffer, send_size, MPI_CHAR,
               iproc, msg_tag, mpi_communicator);
    }
  }
  for(int iproc = 0; iproc < num_procs; ++iproc) {
    if(recv_lists[iproc].size() > 0) {
      MPI_Status status;
      MPI_Wait( &recv_handles[iproc], &status );
    }
  }
}

} // namespace mesh
} // namespace stk

#endif

