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
#include <limits>
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
#include <stk_transfer/TransferUtil.hpp>
#include <stk_transfer/ReducedDependencyCommData.hpp>
#include "stk_util/parallel/CouplingVersions.hpp"
#include "stk_util/parallel/DataExchangeKnownPatternUserDataNonBlocking.hpp"
#include "stk_util/parallel/ParallelReduceBool.hpp"


namespace stk {
namespace transfer {


template <class INTERPOLATE> class ReducedDependencyGeometricTransfer : public TransferBase {

public :

  typedef INTERPOLATE                                     InterpolateClass;
  typedef typename InterpolateClass::MeshA                MeshA;
  typedef typename InterpolateClass::MeshB                MeshB;
  typedef typename InterpolateClass::EntityKeyA           EntityKeyA;
  typedef typename InterpolateClass::EntityKeyB           EntityKeyB;

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

  ReducedDependencyGeometricTransfer(std::shared_ptr<MeshA> mesha,
      std::shared_ptr<MeshB> meshb,
      const std::string &name,
      stk::ParallelMachine pm,
      const double expansion_factor = 1.5,
      const stk::search::SearchMethod search_method = stk::search::KDTREE);

  ReducedDependencyGeometricTransfer(std::shared_ptr<MeshA> mesha,
      std::shared_ptr<MeshB> meshb,
      const std::string &name,
      stk::ParallelMachine pm,
      InterpolateClass interpolate,
      const double expansion_factor = 1.5,
      const stk::search::SearchMethod search_method = stk::search::KDTREE);

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
  void communicate_distances();
  void exchange_transfer_ids();
  void filter_to_nearest(typename MeshB::EntityProcVec to_entity_keys, typename MeshA::EntityProcVec from_entity_keys );

  std::shared_ptr<MeshA>               m_mesha;
  std::shared_ptr<MeshB>               m_meshb;


  const std::string     m_name;
  const double          m_expansion_factor;
  const stk::search::SearchMethod m_search_method;

  EntityProcRelationVec m_domain_to_range;
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
  uniqueProcVec.clear();
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
    const int procid = uniqueProcVec[iproc];
    int offset = -1;
    int numKeys = 0;
    for(int jj = 0; jj < (int) entity_key_proc.size(); ++jj)
    {
      if(static_cast<int>(entity_key_proc[jj].proc()) == procid)
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
  comm.m_otherTransferId.assign(comm.numToMeshCommunications, 0);//zeros for now, will be replaced during exchange_transfer_ids().

  comm.offset_and_num_keys_to_mesh.resize(comm.numToMeshCommunications);
  impl::create_offset_and_num_key(comm.uniqueToProcVec, entity_key_proc_to, comm.offset_and_num_keys_to_mesh );

  comm.offset_and_num_keys_from_mesh.resize(comm.numFromMeshCommunications);
  impl::create_offset_and_num_key(comm.uniqueFromProcVec, entity_key_proc_from, comm.offset_and_num_keys_from_mesh );
}

// Send data from mesh b (domain) to mesh a (range)
template <class MeshAVec, class MeshBVec>
void do_reverse_communication(const ReducedDependencyCommData &comm_data, MeshAVec &a_vec, MeshBVec &b_vec)
{
  using value_t = typename MeshAVec::value_type;

  static_assert(std::is_same<value_t, typename MeshBVec::value_type>::value,
      "Incompatible types for src and dest vector in do_reverse_communication");

  if (stk::util::get_common_coupling_version() > 14) {
    auto exchange = DataExchangeKnownPatternUserDataNonBlocking(comm_data.m_shared_comm);

    std::vector<PointerAndSize> recvs;
    std::vector<int> source_ranks;
    for (int i = 0; i < comm_data.numFromMeshCommunications; ++i) {
      const auto &[recv_offset, recv_size] = comm_data.offset_and_num_keys_from_mesh[i];
      recvs.emplace_back(reinterpret_cast<unsigned char *>(a_vec.data() + recv_offset), sizeof(value_t) * recv_size);
      source_ranks.emplace_back(comm_data.uniqueFromProcVec[i]);
    }

    std::vector<PointerAndSize> sends;
    std::vector<int> dest_ranks;
    for (int i = 0; i < comm_data.numToMeshCommunications; ++i) {
      const auto &[send_offset, send_size] = comm_data.offset_and_num_keys_to_mesh[i];
      sends.emplace_back(reinterpret_cast<unsigned char *>(b_vec.data() + send_offset), sizeof(value_t) * send_size);
      dest_ranks.emplace_back(comm_data.uniqueToProcVec[i]);
    }

    exchange.start_nonblocking(sends, dest_ranks, recvs, source_ranks);
    exchange.complete_receives(recvs, source_ranks, [](int, PointerAndSize){});
    exchange.complete_sends();

  } else {
    // route through integer overflow path
    do_reverse_communication_max_int(comm_data, a_vec, b_vec);
  }
}

// Send data from mesh a (range) to mesh b (domain)
template <class MeshAVec, class MeshBVec>
void do_communication(const ReducedDependencyCommData &comm_data, MeshAVec &a_vec, MeshBVec &b_vec, int stride = 1)
{
  using value_t = typename MeshAVec::value_type;

  static_assert(std::is_same<value_t, typename MeshBVec::value_type>::value,
      "Incompatible types for src and dest vector in do_reverse_communication");

  if (stk::util::get_common_coupling_version() > 14) {
    auto exchange = DataExchangeKnownPatternUserDataNonBlocking(comm_data.m_shared_comm);
    std::vector<PointerAndSize> recvs;
    std::vector<int> source_ranks;

    for (int i = 0; i < comm_data.numToMeshCommunications; ++i){
      const auto source_rank = comm_data.uniqueToProcVec[i];
      const auto & [recv_offset, recv_size] = comm_data.offset_and_num_keys_to_mesh[i];
      recvs.emplace_back(reinterpret_cast<unsigned char *>(b_vec.data() + recv_offset * stride),
          sizeof(value_t) * (recv_size * stride));
      source_ranks.emplace_back(source_rank);
    }

    std::vector<PointerAndSize> sends;
    std::vector<int> dest_ranks;
    for (int i = 0; i < comm_data.numFromMeshCommunications; ++i) {
      const auto dest_rank = comm_data.uniqueFromProcVec[i];
      const auto &[send_offset, send_size] = comm_data.offset_and_num_keys_from_mesh[i];
      sends.emplace_back(reinterpret_cast<unsigned char *>(a_vec.data() + send_offset * stride),
          sizeof(value_t) * (send_size * stride));
      dest_ranks.emplace_back(dest_rank);
    }

    exchange.start_nonblocking(sends, dest_ranks, recvs, source_ranks);
    exchange.complete_receives(recvs, source_ranks, [](int, PointerAndSize){});
    exchange.complete_sends();
  } else {
    // route through integer overflow path
    do_communication_max_int(comm_data, a_vec, b_vec, stride);
  }
}

// These paths potentially hit integer overflow in gpu runs where we partition problem such that we have
// orders of magnitude more nodes / elems per rank (i.e., several million vs tens of thousands). reason for the
// overflow is not only the increased # domain points but also 1) and 2)

// 1) the coarse search can return a large # of avg. candidates. this is especially problematic for vol -> vol xfers
// on high aspect ratio meshes i.e., the bounding boxes are currently not aligned to a local coordinate system

// 2) various buffers (coords, dist, filter mask) scale w/ the avg. number of candidates which is problematic.
// moving forward we should refactor this class such that scaling w/ candidates is avoided where possible, punt 
// on this for now.
template <class MeshAVec, class MeshBVec>
void do_reverse_communication_max_int(
    const ReducedDependencyCommData &comm_data, MeshAVec &a_vec, const MeshBVec &b_vec)
{
  constexpr size_t max_int = std::numeric_limits<int>::max();

  stk::util::print_unsupported_version_warning(14, __LINE__, __FILE__);
  static_assert(std::is_same<typename MeshAVec::value_type, typename MeshBVec::value_type>::value, "Incompatible types for src and dest vector in do_communication");

  std::vector<MPI_Request> receiveRequests(comm_data.numFromMeshCommunications);
  std::vector<MPI_Request> sendRequests(comm_data.numToMeshCommunications);

  // Communicate Mask
  for (int ii = 0; ii < comm_data.numFromMeshCommunications; ++ii)
  {
    int source = comm_data.uniqueFromProcVec[ii];
    const int recv_size = comm_data.offset_and_num_keys_from_mesh[ii].second;
    const int recv_offset = comm_data.offset_and_num_keys_from_mesh[ii].first;
    auto recvMessageSize = recv_size*sizeof(typename MeshAVec::value_type);
    STK_ThrowRequireMsg(recvMessageSize <= max_int,
        "Integer overflow detected during recv in do_reverse_communication_max_int(), use a more recent STK coupling version");

    MPI_Irecv(&a_vec[recv_offset], recvMessageSize, MPI_BYTE, source, MPI_ANY_TAG, comm_data.m_shared_comm,
        &receiveRequests[ii]);
  }

  for (int ii = 0; ii < comm_data.numToMeshCommunications; ++ii)
  {
    int destination = comm_data.uniqueToProcVec[ii];
    const int send_size = comm_data.offset_and_num_keys_to_mesh[ii].second;
    const int send_offset = comm_data.offset_and_num_keys_to_mesh[ii].first;
    auto sendMessageSize = send_size*sizeof(typename MeshBVec::value_type);

    STK_ThrowRequireMsg(sendMessageSize <= max_int,
        "Integer overflow detected during send in do_reverse_communication_max_int(), use a more recent STK coupling version");

    MPI_Isend(&b_vec[send_offset], sendMessageSize, MPI_BYTE, destination,
              0, comm_data.m_shared_comm, &sendRequests[ii]);
  }

  std::vector<MPI_Status> receiveStati(receiveRequests.size());
  MPI_Waitall(receiveRequests.size(), receiveRequests.data(), receiveStati.data());

  std::vector<MPI_Status> sendStati(sendRequests.size());
  MPI_Waitall(sendRequests.size(), sendRequests.data(), sendStati.data());
}


//Send data from mesh a to mesh b
template <class MeshAVec, class MeshBVec>
void do_communication_max_int(const ReducedDependencyCommData & comm_data, const MeshAVec & a_vec, MeshBVec & b_vec, int stride = 1)
{
  constexpr size_t max_int = std::numeric_limits<int>::max();

  stk::util::print_unsupported_version_warning(14, __LINE__, __FILE__);
  static_assert(std::is_same<typename MeshAVec::value_type, typename MeshBVec::value_type>::value, "Incompatible types for src and dest vector in do_communication");
  std::vector<MPI_Request> receiveRequests(comm_data.numToMeshCommunications);
  std::vector<MPI_Request> sendRequests(comm_data.numFromMeshCommunications);

  int sendTag = comm_data.m_transferId;

  for (int ii = 0; ii < comm_data.numToMeshCommunications; ++ii)
  {
    int source = comm_data.uniqueToProcVec[ii];
    int recvTag = comm_data.m_otherTransferId[ii];
    const int recv_size = comm_data.offset_and_num_keys_to_mesh[ii].second * stride;
    const int recv_offset = comm_data.offset_and_num_keys_to_mesh[ii].first * stride;
    auto recvMessageSize = recv_size*sizeof(typename MeshAVec::value_type);

    STK_ThrowRequireMsg(recvMessageSize <= max_int,
        "Integer overflow detected during recv in do_communication_max_int(), use a more recent STK coupling version");

    MPI_Irecv(&b_vec[recv_offset], recvMessageSize, MPI_BYTE, source,
              recvTag, comm_data.m_shared_comm, &receiveRequests[ii]);
  }

  for (int ii = 0; ii < comm_data.numFromMeshCommunications; ++ii)
  {
    int destination = comm_data.uniqueFromProcVec[ii];
    const int send_size = comm_data.offset_and_num_keys_from_mesh[ii].second * stride;
    const int send_offset = comm_data.offset_and_num_keys_from_mesh[ii].first * stride;

    auto sendMessageSize = send_size*sizeof(typename MeshBVec::value_type);

    STK_ThrowRequireMsg(sendMessageSize <= max_int,
        "Integer overflow detected during send in do_communication_max_int(), use a more recent STK coupling version");

    MPI_Isend(&a_vec[send_offset], sendMessageSize, MPI_BYTE, destination,
              sendTag, comm_data.m_shared_comm, &sendRequests[ii]);
  }

  std::vector<MPI_Status> receiveStati(receiveRequests.size());
  MPI_Waitall(receiveRequests.size(), receiveRequests.data(), MPI_STATUSES_IGNORE);

  std::vector<MPI_Status> sendStati(sendRequests.size());
  MPI_Waitall(sendRequests.size(), sendRequests.data(), sendStati.data());
}

template <class INTERPOLATE> ReducedDependencyGeometricTransfer<INTERPOLATE>::ReducedDependencyGeometricTransfer
(std::shared_ptr<MeshA> mesha,
 std::shared_ptr<MeshB> meshb,
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
    STK_ThrowRequire(mesha || meshb);

  }

  template <class INTERPOLATE>
  ReducedDependencyGeometricTransfer<INTERPOLATE>::ReducedDependencyGeometricTransfer(std::shared_ptr<MeshA> mesha,
      std::shared_ptr<MeshB> meshb,
      const std::string &name,
      stk::ParallelMachine pm,
      InterpolateClass interpolate,
      const double expansion_factor,
      const stk::search::SearchMethod search_method)
      : m_mesha(mesha),
        m_meshb(meshb),
        m_name(name),
        m_expansion_factor(expansion_factor),
        m_search_method(search_method),
        m_interpolate(std::move(interpolate))

  {
    //In an mpmd program, there's no guarantee that the types specified for the entity keys are honored by each program,
    //so for now, enforce that the types are 64bit for consistency during mpi comms
    static_assert(8 == sizeof(typename InterpolateClass::EntityKeyA), "Size of EntityKeyA needs to be 64 bit");
    static_assert(8 == sizeof(typename InterpolateClass::EntityKeyB), "Size of EntityKeyB needs to be 64 bit");
    m_comm_data.m_shared_comm = pm;
    STK_ThrowRequire(mesha || meshb);
  }

  template <class INTERPOLATE>
  void ReducedDependencyGeometricTransfer<INTERPOLATE>::coarse_search()
  {
    // for a given input (domain), we want the output (range) candidates
    // e.g., nodes of mesh B -> candidate elements of mesh A

    // in the context of MPMD
    // for the domain app domain_to_range is pairs of (local domain ids, off rank candidate range ids)
    // for the range app domain_to_range is pairs of (off rank domain ids, local candidate range ids)

    impl::coarse_search_impl<INTERPOLATE>(m_domain_to_range, m_comm_data.m_shared_comm, m_mesha.get(),
        m_meshb.get(), m_search_method, m_expansion_factor);
  }

template <class INTERPOLATE> void ReducedDependencyGeometricTransfer<INTERPOLATE>::communication() {  
  typename MeshB::EntityProcVec to_entity_keys;
  typename MeshA::EntityProcVec from_entity_keys;

  determine_entities_to_copy(to_entity_keys, from_entity_keys);

  {
    // domain to range is potentially very large e.g. # domain points * avg candidates * 32 byte
    // do not hold on to it as we don't currently use it after this point
    EntityProcRelationVec().swap(m_domain_to_range);
  }

  if (m_meshb)
    m_meshb->get_to_points_coordinates(to_entity_keys, to_points_on_to_mesh);

  buildExchangeLists(to_entity_keys, from_entity_keys);
  to_points_distance_on_from_mesh.resize(from_entity_keys.size());
  if(m_mesha)
    m_interpolate.obtain_parametric_coords(from_entity_keys, *m_mesha, to_points_on_from_mesh, to_points_distance_on_from_mesh);
  communicate_distances();
  filter_to_nearest(to_entity_keys, from_entity_keys);
  
  const auto coupling_version = stk::util::get_common_coupling_version();
  if (coupling_version >= 11 and coupling_version <= 14) {
    exchange_transfer_ids();
  }
}

template <class INTERPOLATE> void ReducedDependencyGeometricTransfer<INTERPOLATE>::local_search()
{
  //no op (filtering done at communication)
}

template <class INTERPOLATE>
void ReducedDependencyGeometricTransfer<INTERPOLATE>::exchange_transfer_ids()
{
  impl::exchange_transfer_ids(m_comm_data);
}

template <class INTERPOlATE>
void ReducedDependencyGeometricTransfer<INTERPOlATE>::filter_to_nearest(typename MeshB::EntityProcVec to_entity_keys,
    typename MeshA::EntityProcVec from_entity_keys )
{
  // Find the winner
  std::map<EntityKeyB, std::pair<double, int> > filterMap;

  for (const auto &[offset, num_keys] : m_comm_data.offset_and_num_keys_to_mesh) {
    for(int jj =0; jj < num_keys; ++jj){
      STK_ThrowRequireMsg(offset + jj < (int) to_points_distance_on_to_mesh.size(),
          "'offset+jj' (" << offset << "+" << jj << ") required to be less than to_points_distance_on_to_mesh.size() ("
                          << to_points_distance_on_to_mesh.size() << ")");
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

  if (stk::util::get_common_coupling_version() >= 15)
  {
    bool exception_on_local_proc = false;
    std::string error_msg;
    try {
      m_interpolate.mask_parametric_coords(FilterMaskFrom, from_count);
      error_msg = "error on another proc";
    } catch (std::exception& e)
    {
      exception_on_local_proc = true;
      error_msg = e.what();
    }
    
    bool exception_on_any_proc = stk::is_true_on_any_proc(m_comm_data.m_shared_comm, exception_on_local_proc);
    STK_ThrowRequireMsg(!exception_on_any_proc, error_msg);
  } else
  {
    m_interpolate.mask_parametric_coords(FilterMaskFrom, from_count);
  }


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

template <class INTERPOLATE>
void ReducedDependencyGeometricTransfer<INTERPOLATE>::determine_entities_to_copy(
    typename MeshB::EntityProcVec &entities_to_copy_to, typename MeshA::EntityProcVec &entities_to_copy_from) const
{
  entities_to_copy_to.clear();
  entities_to_copy_from.clear();

  ParallelMachine comm = m_comm_data.m_shared_comm;
  const unsigned my_rank = parallel_machine_rank(comm);
  for (const auto &[domain, range] : m_domain_to_range) {
    const unsigned range_rank = range.proc();
    const unsigned domain_rank = domain.proc();

    if (domain_rank == my_rank) {
      entities_to_copy_to.emplace_back(domain.id(), range_rank);
    }

    if (range_rank == my_rank) {
      entities_to_copy_from.emplace_back(range.id(), domain_rank);
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

  int sendTag = m_comm_data.m_transferId;

  for(int ii = 0; ii < m_comm_data.numToMeshCommunications; ++ii)
  {
    int destination = m_comm_data.uniqueToProcVec[ii];
    MPI_Isend(&sendSizesBuffers[ii], 1, MPI_INT, destination, sendTag, m_comm_data.m_shared_comm, &sendRequests[ii]);
  }

  std::vector<MPI_Status> receiveStati;
  receiveStati.resize(receiveRequests.size());
  MPI_Waitall(receiveRequests.size(), receiveRequests.data(), receiveStati.data());

  std::vector<MPI_Status> sendStati;
  sendStati.resize(sendRequests.size());
  MPI_Waitall(sendRequests.size(), sendRequests.data(), sendStati.data());

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

