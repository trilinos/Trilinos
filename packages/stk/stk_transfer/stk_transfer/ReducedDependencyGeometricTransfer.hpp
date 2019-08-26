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

#ifndef  ReducedDependecy_STK_GEOMETRICTRANSFER_HPP
#define  ReducedDependecy_STK_GEOMETRICTRANSFER_HPP

#include <algorithm>
#include <set>
#include <string>
#include <vector>
#include <memory>

#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <stk_util/util/StaticAssert.hpp>

#include <stk_search/CoarseSearch.hpp>
#include <stk_search/SearchMethod.hpp>
#include <stk_transfer/GeometricTransferImpl.hpp>
#include <stk_transfer/TransferBase.hpp>

#include<impl/Kokkos_Timer.hpp>


namespace stk {
namespace transfer {


struct ReducedDependencyCommData
{
  std::vector<std::pair <int, int>> offset_and_num_keys_to_mesh;
  std::vector<std::pair <int, int>> offset_and_num_keys_from_mesh;
  int numToMeshCommunications;
  int numFromMeshCommunications;
  std::vector<int> uniqueFromProcVec;
  std::vector<int> uniqueToProcVec;
  stk::ParallelMachine m_shared_comm = MPI_COMM_WORLD;
};

template <class INTERPOLATE> class ReducedDependencyGeometricTransfer : public TransferBase {

public :

  typedef INTERPOLATE                                     InterpolateClass;
  typedef typename InterpolateClass::MeshA                MeshA;
  typedef typename InterpolateClass::MeshB                MeshB;
  typedef typename InterpolateClass::EntityKeyA           EntityKeyA;
  typedef typename InterpolateClass::EntityKeyB           EntityKeyB;
  typedef typename InterpolateClass::EntityKeyMap         EntityKeyMap;

  typedef typename InterpolateClass::EntityProcA          EntityProcA;
  typedef typename InterpolateClass::EntityProcB          EntityProcB;

  typedef typename InterpolateClass::EntityProcRelation       EntityProcRelation;
  typedef typename InterpolateClass::EntityProcRelationVec    EntityProcRelationVec;

  typedef typename std::set     <EntityKeyA>              EntityKeySetA;
  typedef typename std::set     <EntityKeyB>              EntityKeySetB;
  typedef typename MeshA::BoundingBox                     BoundingBoxA;
  typedef typename MeshB::BoundingBox                     BoundingBoxB;

  typedef typename MeshB::Point             ToPoint;
  typedef typename MeshB::ToPointsContainer ToPointsContainer;
  typedef typename MeshB::ToPointsDistanceContainer ToPointsDistanceContainer;
  typedef stk::search::Point<double>         Point;

  enum {Dimension = 3};

  ReducedDependencyGeometricTransfer(std::shared_ptr<MeshA> &mesha,
                    std::shared_ptr<MeshB> &meshb,
                    const std::string &name,
                    stk::ParallelMachine pm,
                    const double expansion_factor = 1.5,
                    const stk::search::SearchMethod search_method = stk::search::BOOST_RTREE);
  virtual ~ReducedDependencyGeometricTransfer(){};
  void coarse_search() override;
  void communication() override;
  void local_search() override;
  void apply() override;

  void determine_entities_to_copy(typename MeshB::EntityProcVec   &entities_to_copy_to,
                                  typename MeshA::EntityProcVec   &entities_to_copy_from ) const;

  const std::shared_ptr<MeshA> mesha() const {return m_mesha;}
  const std::shared_ptr<MeshB> meshb() const {return m_meshb;}

private :

  void buildExchangeLists(typename MeshB::EntityProcVec entity_key_proc_to, typename MeshA::EntityProcVec entity_key_proc_from);
  void communicate_distances ();
  void filter_to_nearest (typename MeshB::EntityProcVec to_entity_keys, typename MeshA::EntityProcVec from_entity_keys );

  std::shared_ptr<MeshA>               m_mesha;
  std::shared_ptr<MeshB>               m_meshb;


  const std::string     m_name;
  const double          m_expansion_factor;
  const stk::search::SearchMethod m_search_method;

  EntityProcRelationVec m_global_range_to_domain;
  InterpolateClass      m_interpolate;
  ToPointsContainer to_points_on_to_mesh;
  ToPointsContainer to_points_on_from_mesh;
  ToPointsDistanceContainer to_points_distance_on_from_mesh;
  ToPointsDistanceContainer to_points_distance_on_to_mesh;
  ReducedDependencyCommData m_comm_data;
  typename MeshB::EntityProcVec to_entity_keys_masked;
  typename MeshA::EntityProcVec from_entity_keys_masked;
};

namespace impl {

template <class Mesh>
void get_unique_procs_from_entity_keys(const typename Mesh::EntityProcVec & vec, std::vector<int> & uniqueProcVec)
{
  //recalculate uniqueProcs using masked values
  uniqueProcVec.reserve(vec.size());
  std::transform(vec.begin(), vec.end(), std::back_inserter(uniqueProcVec),
      [](typename Mesh::EntityProc a) {return a.proc();});
  std::sort(uniqueProcVec.begin(), uniqueProcVec.end());
  auto last = std::unique(uniqueProcVec.begin(), uniqueProcVec.end());
  uniqueProcVec.erase(last, uniqueProcVec.end());
}

template <class VecType>
void copy_according_to_mask(const VecType & src, VecType & dest, const std::vector<int> & mask)
{
  int offset = 0;
  for (unsigned int ii = 0; ii < src.size(); ++ii )
  {
    if (mask[ii])
    {
      dest[offset]= src[ii];
      offset++;
    }
  }
}


template <class VecType>
void create_offset_and_num_key(const std::vector<int> & uniqueProcVec,
    const VecType & entity_key_proc,
    std::vector<std::pair <int, int>> & offset_and_num_keys)
{
  for(unsigned iproc = 0; iproc < uniqueProcVec.size(); ++iproc)
  {
    const unsigned procid = uniqueProcVec[iproc];
    int offset = -1;
    int numKeys = 0;
    for(int jj = 0; jj < (int) entity_key_proc.size(); ++jj)
    {
      if(entity_key_proc[jj].proc() == procid)
      {
        if(offset < 0)
          offset = jj;
        numKeys += 1;
      }
    }
    offset_and_num_keys[iproc] = std::make_pair(offset, numKeys);
  }
}
}

template <class MeshA, class MeshB>
void construct_comm_data(const typename MeshA::EntityProcVec & entity_key_proc_from,
    const typename MeshB::EntityProcVec entity_key_proc_to,
    ReducedDependencyCommData & comm)
{
  impl::get_unique_procs_from_entity_keys<MeshA>(entity_key_proc_from, comm.uniqueFromProcVec);
  comm.numFromMeshCommunications = comm.uniqueFromProcVec.size();

  impl::get_unique_procs_from_entity_keys<MeshB>(entity_key_proc_to, comm.uniqueToProcVec);
  comm.numToMeshCommunications = comm.uniqueToProcVec.size();

  comm.offset_and_num_keys_to_mesh.resize(comm.numToMeshCommunications);
  impl::create_offset_and_num_key(comm.uniqueToProcVec, entity_key_proc_to, comm.offset_and_num_keys_to_mesh );

  comm.offset_and_num_keys_from_mesh.resize(comm.numFromMeshCommunications);
  impl::create_offset_and_num_key(comm.uniqueFromProcVec, entity_key_proc_from, comm.offset_and_num_keys_from_mesh );
}

//Send data from mesh b to mesh a
template <class MeshAVec, class MeshBVec>
void do_reverse_communication(const ReducedDependencyCommData & comm_data, MeshAVec & a_vec, const MeshBVec & b_vec)
{
  static_assert(std::is_same<typename MeshAVec::value_type, typename MeshBVec::value_type>::value, "Incompatible types for src and dest vector in do_communication");

  std::vector<MPI_Request> receiveRequests(comm_data.numFromMeshCommunications);
  std::vector<MPI_Request> sendRequests(comm_data.numToMeshCommunications);

  // Communicate Mask
  for (int ii = 0; ii < comm_data.numFromMeshCommunications; ++ii)
  {
    int source = comm_data.uniqueFromProcVec[ii];
    const int recv_size = comm_data.offset_and_num_keys_from_mesh[ii].second;
    const int recv_offset = comm_data.offset_and_num_keys_from_mesh[ii].first;
    int recvMessageSize = recv_size*sizeof(typename MeshAVec::value_type);

    MPI_Irecv(&a_vec[recv_offset], recvMessageSize, MPI_BYTE, source,
              MPI_ANY_TAG, comm_data.m_shared_comm, &receiveRequests[ii]);
  }

  for (int ii = 0; ii < comm_data.numToMeshCommunications; ++ii)
  {
    int destination = comm_data.uniqueToProcVec[ii];
    const int send_size = comm_data.offset_and_num_keys_to_mesh[ii].second;
    const int send_offset = comm_data.offset_and_num_keys_to_mesh[ii].first;
    int sendMessageSize = send_size*sizeof(typename MeshBVec::value_type);

    MPI_Isend(&b_vec[send_offset], sendMessageSize, MPI_BYTE, destination,
              0, comm_data.m_shared_comm, &sendRequests[ii]);
  }

  std::vector<MPI_Status> receiveStati(receiveRequests.size());
  MPI_Waitall(receiveRequests.size(), &receiveRequests[0], &receiveStati[0] );

  std::vector<MPI_Status> sendStati(sendRequests.size());
  MPI_Waitall(sendRequests.size(), &sendRequests[0], &sendStati[0] );
}


//Send data from mesh a to mesh b
template <class MeshAVec, class MeshBVec>
void do_communication(const ReducedDependencyCommData & comm_data, const MeshAVec & a_vec, MeshBVec & b_vec, int stride = 1)
{
  static_assert(std::is_same<typename MeshAVec::value_type, typename MeshBVec::value_type>::value, "Incompatible types for src and dest vector in do_communication");
  std::vector<MPI_Request> receiveRequests(comm_data.numToMeshCommunications);
  std::vector<MPI_Request> sendRequests(comm_data.numFromMeshCommunications);


  for (int ii = 0; ii < comm_data.numToMeshCommunications; ++ii)
  {
    int source = comm_data.uniqueToProcVec[ii];
    const int recv_size = comm_data.offset_and_num_keys_to_mesh[ii].second * stride;
    const int recv_offset = comm_data.offset_and_num_keys_to_mesh[ii].first * stride;
    int recvMessageSize = recv_size*sizeof(typename MeshAVec::value_type);

    MPI_Irecv(&b_vec[recv_offset], recvMessageSize, MPI_BYTE, source,
              MPI_ANY_TAG, comm_data.m_shared_comm, &receiveRequests[ii]);
  }

  for (int ii = 0; ii < comm_data.numFromMeshCommunications; ++ii)
  {
    int destination = comm_data.uniqueFromProcVec[ii];
    const int send_size = comm_data.offset_and_num_keys_from_mesh[ii].second * stride;
    const int send_offset = comm_data.offset_and_num_keys_from_mesh[ii].first * stride;

    int sendMessageSize = send_size*sizeof(typename MeshBVec::value_type);

    MPI_Isend(&a_vec[send_offset], sendMessageSize, MPI_BYTE, destination,
              0, comm_data.m_shared_comm, &sendRequests[ii]);
  }

  std::vector<MPI_Status> receiveStati(receiveRequests.size());
  MPI_Waitall(receiveRequests.size(), &receiveRequests[0], &receiveStati[0] );

  std::vector<MPI_Status> sendStati(sendRequests.size());
  MPI_Waitall(sendRequests.size(), &sendRequests[0], &sendStati[0] );
}



template <class INTERPOLATE> ReducedDependencyGeometricTransfer<INTERPOLATE>::ReducedDependencyGeometricTransfer
(std::shared_ptr<MeshA> &mesha,
 std::shared_ptr<MeshB> &meshb,
 const std::string        &name,
 stk::ParallelMachine pm,
 const double              expansion_factor,
 const stk::search::SearchMethod search_method) :
      m_mesha(mesha), m_meshb(meshb), m_name(name), m_expansion_factor(expansion_factor), m_search_method(search_method)
  {
    //In an mpmd program, there's no guarantee that the types specified for the entity keys are honored by each program,
    //so for now, enforce that the types are 64bit for consistency during mpi comms
    static_assert(8 == sizeof(typename InterpolateClass::EntityKeyA), "Size of EntityKeyA needs to be 64 bit");
    static_assert(8 == sizeof(typename InterpolateClass::EntityKeyB), "Size of EntityKeyB needs to be 64 bit");
    m_comm_data.m_shared_comm = pm;
    ThrowRequire(mesha || meshb);
  }

template <class INTERPOLATE> void ReducedDependencyGeometricTransfer<INTERPOLATE>::coarse_search() {

  m_global_range_to_domain.clear();
   impl::coarse_search_impl<INTERPOLATE>(m_global_range_to_domain,
                m_comm_data.m_shared_comm,
                m_mesha.get(),
                m_meshb.get(),
                m_search_method,
                m_expansion_factor);
}

template <class INTERPOLATE> void ReducedDependencyGeometricTransfer<INTERPOLATE>::communication() {

  typename MeshB::EntityProcVec to_entity_keys;
  typename MeshA::EntityProcVec from_entity_keys;

  determine_entities_to_copy(to_entity_keys, from_entity_keys);

  if (m_meshb)
    m_meshb->get_to_points_coordinates(to_entity_keys, to_points_on_to_mesh);

  buildExchangeLists(to_entity_keys, from_entity_keys);
  to_points_distance_on_from_mesh.resize(from_entity_keys.size());
  if(m_mesha)
    m_interpolate.obtain_parametric_coords(from_entity_keys, *m_mesha, to_points_on_from_mesh, to_points_distance_on_from_mesh);
  communicate_distances();
  filter_to_nearest(to_entity_keys, from_entity_keys);
}

template <class INTERPOLATE> void ReducedDependencyGeometricTransfer<INTERPOLATE>::local_search()
{
  //no op (filtering done at communication)
}

template <class INTERPOlATE>
void ReducedDependencyGeometricTransfer<INTERPOlATE>::filter_to_nearest(typename MeshB::EntityProcVec to_entity_keys,
    typename MeshA::EntityProcVec from_entity_keys )
{

  // Find the winner
  std::map<EntityKeyB, std::pair<double, int> > filterMap;

  for (unsigned int ii = 0; ii < m_comm_data.offset_and_num_keys_to_mesh.size(); ++ii )
  {
    int offset = m_comm_data.offset_and_num_keys_to_mesh[ii].first;
    for(int jj =0; jj < m_comm_data.offset_and_num_keys_to_mesh[ii].second; ++jj){
      ThrowRequire(offset+jj < (int)to_points_distance_on_to_mesh.size());
      std::pair<double,int> dist_and_to_entity_index = std::make_pair(to_points_distance_on_to_mesh[offset+jj], offset+jj);
      auto key = to_entity_keys[offset+jj].id();
      if ( filterMap.find(key) == filterMap.end() )
      {
        filterMap[key] = dist_and_to_entity_index;
      }
      else
      {
        auto current_dist_and_index = filterMap.at(key);
        if (to_points_distance_on_to_mesh[offset+jj] < current_dist_and_index.first) {
          filterMap[key] = dist_and_to_entity_index;
        }
      }
    }
  }


  // Create Mask of to_entity_keys winners, if the key is in the map give the mask 1 otherwise 0
  std::vector<int> FilterMaskTo(to_entity_keys.size(),0);

  for (unsigned int ii = 0; ii < to_entity_keys.size(); ++ii )
  {
    auto key = to_entity_keys[ii].id();
    auto iter = filterMap.find(key);
    //key must match and must be the same entity index
    if( iter != filterMap.end() && iter->second.second == static_cast<int>(ii))
    {
      FilterMaskTo[ii] = 1;
    }
    else
    {
      FilterMaskTo[ii] = 0;
    }
  }

  std::vector<int> FilterMaskFrom(from_entity_keys.size(),0);
  do_reverse_communication(m_comm_data, FilterMaskFrom, FilterMaskTo);

  const int to_count = std::count(FilterMaskTo.begin(), FilterMaskTo.end(), 1);
  const int from_count = std::count(FilterMaskFrom.begin(), FilterMaskFrom.end(), 1);
  m_interpolate.mask_parametric_coords(FilterMaskFrom, from_count);


  from_entity_keys_masked.resize(from_count);
  impl::copy_according_to_mask(from_entity_keys, from_entity_keys_masked, FilterMaskFrom);

  to_entity_keys_masked.resize(to_count);
  impl::copy_according_to_mask(to_entity_keys, to_entity_keys_masked, FilterMaskTo);

  construct_comm_data<MeshA, MeshB>(from_entity_keys_masked, to_entity_keys_masked, m_comm_data);
}



template <class INTERPOLATE> void ReducedDependencyGeometricTransfer<INTERPOLATE>::apply()
{
  if (m_mesha)
    m_mesha->update_values();

  m_interpolate.apply(m_meshb.get(),
      m_mesha.get(),
      to_entity_keys_masked,
      from_entity_keys_masked,
      m_comm_data);

  if(m_meshb)
    m_meshb->update_values();
}

template <class INTERPOLATE> void ReducedDependencyGeometricTransfer<INTERPOLATE>::determine_entities_to_copy(
                         typename MeshB::EntityProcVec   &entities_to_copy_to,
                         typename MeshA::EntityProcVec   &entities_to_copy_from ) const {

  entities_to_copy_to.clear();
  entities_to_copy_from.clear();

  ParallelMachine comm = m_comm_data.m_shared_comm;
  const unsigned my_rank = parallel_machine_rank(comm);

  for (auto && elem : m_global_range_to_domain)
  {
    const unsigned domain_owning_rank = elem.second.proc();
    const unsigned range_owning_rank = elem.first.proc();

    if (range_owning_rank == my_rank) {
      const EntityKeyB entity = elem.first.id();
      const typename MeshB::EntityProc ep(entity, domain_owning_rank);
      entities_to_copy_to.push_back(ep);
    }
    if (domain_owning_rank == my_rank) {
      const EntityKeyA entity = elem.second.id();
      const typename MeshA::EntityProc ep (entity, range_owning_rank);
      entities_to_copy_from.push_back(ep);
    }
  }
  std::sort(entities_to_copy_to.begin(), entities_to_copy_to.end());

}

template<class INTERPOLATE>
void
ReducedDependencyGeometricTransfer<INTERPOLATE>::buildExchangeLists(typename MeshB::EntityProcVec entity_key_proc_to,
    typename MeshA::EntityProcVec entity_key_proc_from)
{

  construct_comm_data<MeshA, MeshB>(entity_key_proc_from, entity_key_proc_to, m_comm_data);

  //allocate buffers
  std::vector<int> receiveSizesBuffers(m_comm_data.numFromMeshCommunications);
  std::vector<int> sendSizesBuffers(m_comm_data.numToMeshCommunications);
  std::vector<MPI_Request> receiveRequests(m_comm_data.numFromMeshCommunications);
  std::vector<MPI_Request> sendRequests(m_comm_data.numToMeshCommunications);

  //communicate buffer sizes
  for(int ii = 0; ii < m_comm_data.numToMeshCommunications; ++ii)
  {
    sendSizesBuffers[ii] = m_comm_data.offset_and_num_keys_to_mesh[ii].second;
  }

  for(int ii = 0; ii < m_comm_data.numFromMeshCommunications; ++ii)
  {
    int source = m_comm_data.uniqueFromProcVec[ii];
    MPI_Irecv(&receiveSizesBuffers[ii], 1, MPI_INT, source, MPI_ANY_TAG, m_comm_data.m_shared_comm, &receiveRequests[ii]);
  }

  for(int ii = 0; ii < m_comm_data.numToMeshCommunications; ++ii)
  {
    int destination = m_comm_data.uniqueToProcVec[ii];
    MPI_Isend(&sendSizesBuffers[ii], 1, MPI_INT, destination, 0, m_comm_data.m_shared_comm, &sendRequests[ii]);
  }

  std::vector<MPI_Status> receiveStati(receiveRequests.size());
  MPI_Waitall(receiveRequests.size(), &receiveRequests[0], &receiveStati[0]);

  std::vector<MPI_Status> sendStati(sendRequests.size());
  MPI_Waitall(sendRequests.size(), &sendRequests[0], &sendStati[0]);

  //communicate coordinates
  int numRecvPoints = 0;
  for(unsigned ii = 0; ii < receiveSizesBuffers.size(); ++ii)
  {
    numRecvPoints += receiveSizesBuffers[ii];
  }
  to_points_on_from_mesh.resize(numRecvPoints);

  do_reverse_communication(m_comm_data, to_points_on_from_mesh, to_points_on_to_mesh);
}

template <class INTERPOLATE>  void ReducedDependencyGeometricTransfer<INTERPOLATE>::communicate_distances ()
{
  to_points_distance_on_to_mesh.resize(to_points_on_to_mesh.size());
  do_communication(m_comm_data, to_points_distance_on_from_mesh, to_points_distance_on_to_mesh);
}


}
}


#endif

