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

#ifndef COMMON_SEARCH_UTIL_H_
#define COMMON_SEARCH_UTIL_H_

#ifdef _OPENMP
#include <omp.h>
#endif

#include "stk_util/environment/WallTime.hpp"
#include "stk_util/parallel/Parallel.hpp"
#include "stk_util/parallel/CommSparse.hpp"
#include "stk_util/parallel/ParallelComm.hpp"
#include "stk_util/util/SortAndUnique.hpp"
#include "stk_search/BoxIdent.hpp"
#include "stk_search/kdtree/KDTree_BoundingBox.hpp"
#include "stk_search/kdtree/KDTree.hpp"
#include "DeviceMPIUtils.hpp"
#include <any>
#include <memory>

namespace stk::search {

template <typename DomainBox>
inline void global_box_combine(DomainBox &box_array, MPI_Comm &communicator)
{
  typedef typename DomainBox::value_type::coordinate_t coordinate_t;

  int num_boxes = static_cast<int>(box_array.size());
  //
  //  Allocate a common set of arrays to perform the reductions on
  //
  int array_length = num_boxes * 3;
  std::vector<coordinate_t> all_box_min_local (array_length);
  std::vector<coordinate_t> all_box_max_local (array_length);
  std::vector<coordinate_t> all_box_min_global(array_length);
  std::vector<coordinate_t> all_box_max_global(array_length);
  //
  //  Fill the local arrays
  //
  for(int ibox = 0; ibox < num_boxes; ++ibox) {
    all_box_min_local[ibox * 3 + 0] = box_array[ibox].GetBox().get_x_min();
    all_box_min_local[ibox * 3 + 1] = box_array[ibox].GetBox().get_y_min();
    all_box_min_local[ibox * 3 + 2] = box_array[ibox].GetBox().get_z_min();
    all_box_max_local[ibox * 3 + 0] = box_array[ibox].GetBox().get_x_max();
    all_box_max_local[ibox * 3 + 1] = box_array[ibox].GetBox().get_y_max();
    all_box_max_local[ibox * 3 + 2] = box_array[ibox].GetBox().get_z_max();
  }
  //
  //  Perform the global MPI reductions
  //
  MPI_Datatype floatType;
  if(sizeof(coordinate_t) == sizeof(float)) {
    floatType = MPI_FLOAT;
  }
  else if (sizeof(coordinate_t) == sizeof(double)) {
    floatType = MPI_DOUBLE;
  }
  else {
    floatType = MPI_DOUBLE;
  }
  MPI_Allreduce(all_box_min_local.data(), all_box_min_global.data(), array_length, floatType, MPI_MIN, communicator);
  MPI_Allreduce(all_box_max_local.data(), all_box_max_global.data(), array_length, floatType, MPI_MAX, communicator);
  //
  //  Scatter the local arrays back to the boxes
  //
  for(int ibox = 0; ibox < num_boxes; ++ibox) {
    box_array[ibox].GetBox().set_box(all_box_min_global[ibox * 3 + 0],
                                     all_box_min_global[ibox * 3 + 1],
                                     all_box_min_global[ibox * 3 + 2],
                                     all_box_max_global[ibox * 3 + 0],
                                     all_box_max_global[ibox * 3 + 1],
                                     all_box_max_global[ibox * 3 + 2]);
  }
}

//
//  Exchange boxes so that current proc local box is sent to the global box on all processors.
//
template <typename DomainBox>
inline void all_gather_helper(const DomainBox& localBox, std::vector<DomainBox> &global_box_array, MPI_Comm &comm)
{
  int numProc;
  MPI_Comm_size(comm, &numProc);
  global_box_array.clear();
  global_box_array.resize(numProc);
  const char* localDataConst = reinterpret_cast<const char*>(&localBox);
  //  NKC, hack to support old MPI version used by Goodyear
  char* localData = const_cast<char*>(localDataConst);

  MPI_Allgather(localData, sizeof(DomainBox), MPI_CHAR, global_box_array.data(), sizeof(DomainBox), MPI_CHAR, comm);
}

template <typename DomainBox, typename ExecutionSpace>
void all_gather_helper(const DomainBox& localBox, Kokkos::View<DomainBox*, ExecutionSpace> global_box_array,
                       MPI_Comm comm)
{
  int commSize = stk::parallel_machine_size(comm);
  if (static_cast<int>(global_box_array.extent(0)) < commSize)
  {
    Kokkos::resize(global_box_array, commSize);
  }
  auto global_box_array_host = Kokkos::create_mirror_view(global_box_array);

  impl::check_view_c_layout(global_box_array_host);

  MPI_Allgather(&localBox, sizeof(DomainBox), MPI_CHAR,
                global_box_array_host.data(), sizeof(DomainBox), MPI_CHAR, comm);
  Kokkos::deep_copy(global_box_array, global_box_array_host);
}

template <typename DataType>
inline void concatenate_thread_lists(const std::vector<std::vector<DataType>> &vectorIn,
                                     std::vector<DataType> &vectorOut)
{
  const unsigned numThreadLists = vectorIn.size();
  std::vector<unsigned> offsets(numThreadLists);
  unsigned totSize = 0;
  for (unsigned i = 0; i < numThreadLists; ++i) {
    offsets[i] = totSize;
    totSize += vectorIn[i].size();
  }

  vectorOut.resize(totSize);

#ifdef _OPENMP
#pragma omp parallel default(shared)
  {
    const unsigned ithread = omp_get_thread_num();
    const std::vector<DataType> &data = vectorIn[ithread];
    std::copy(data.begin(), data.end(), &vectorOut[offsets[ithread]]);
  }
#else
  for (unsigned ithread = 0; ithread < numThreadLists; ++ithread) {
    const std::vector<DataType> &data = vectorIn[ithread];
    std::copy(data.begin(), data.end(), &vectorOut[offsets[ithread]]);
  }
#endif

}

template <typename DomainKey, typename RangeKey>
void communicate_vector(stk::ParallelMachine arg_comm,
                        std::vector<std::pair<DomainKey, RangeKey>> & search_relations,
                        bool enforceSearchResultSymmetry = true)
{
  typedef std::pair<DomainKey, RangeKey> ValueType;

  CommSparse commSparse(arg_comm);

  const int p_rank = commSparse.parallel_rank();
  const int p_size = commSparse.parallel_size();

  if (1 == p_size) {
    return;
  }

  typename std::vector< ValueType >::const_iterator i;

  size_t numLocal = 0;
  for (i = search_relations.begin(); i != search_relations.end(); ++i) {
    const ValueType & val = *i;
    if (static_cast<int>(val.first.proc()) == p_rank ||
        (enforceSearchResultSymmetry && static_cast<int>(val.second.proc()) == p_rank)) {
      ++numLocal;
    }
    if (static_cast<int>(val.first.proc()) != p_rank) {
      CommBuffer & buf = commSparse.send_buffer(val.first.proc());
      buf.skip<ValueType>(1);
    }
    if (enforceSearchResultSymmetry) {
      if (static_cast<int>(val.second.proc()) != p_rank && val.second.proc() != val.first.proc()) {
        CommBuffer & buf = commSparse.send_buffer(val.second.proc());
        buf.skip<ValueType>(1);
      }
    }
  }

  commSparse.allocate_buffers();

  for (i = search_relations.begin(); i != search_relations.end(); ++i) {
    const ValueType & val = *i;
    if (static_cast<int>(val.first.proc()) != p_rank) {
      CommBuffer & buf = commSparse.send_buffer(val.first.proc());
      buf.pack<ValueType>(val);
    }
    if (enforceSearchResultSymmetry) {
      if (static_cast<int>(val.second.proc()) != p_rank && val.second.proc() != val.first.proc()) {
        CommBuffer & buf = commSparse.send_buffer(val.second.proc());
        buf.pack<ValueType>(val);
      }
    }
  }

  commSparse.communicate();

  size_t numRecvd = 0;
  for (int p = 0; p < p_size; ++p) {
    CommBuffer & buf = commSparse.recv_buffer(p);
    numRecvd += (buf.remaining()/sizeof(ValueType));
  }
  search_relations.reserve(numLocal+numRecvd);

  size_t keep = 0;
  for (size_t j=0; j<search_relations.size(); ++j) {
    const ValueType & val = search_relations[j];
    if (static_cast<int>(val.first.proc()) == p_rank ||
        (enforceSearchResultSymmetry && static_cast<int>(val.second.proc()) == p_rank)) {
      if (j > keep) {
        search_relations[keep] = val;
      }
      ++keep;
    }
  }

  for (int p = 0; p < p_size; ++p) {
    CommBuffer & buf = commSparse.recv_buffer(p);
    while (buf.remaining()) {
      ValueType val;
      buf.unpack<ValueType>(val);
      search_relations.push_back(val);
    }
  }
}

template <typename SearchRelationType>
void communicate_views(stk::ParallelMachine arg_comm,
                       SearchRelationType& search_relations,
                       bool enforceSearchResultSymmetry = true)
{
  using ValueType = typename SearchRelationType::value_type;
  CommSparse commSparse(arg_comm);

  const int p_rank = commSparse.parallel_rank();
  const int p_size = commSparse.parallel_size();

  if (1 == p_size) {
    return;
  }

  size_t numLocal = 0;
  for (unsigned i = 0; i != search_relations.extent(0); ++i) {
    const ValueType & val = search_relations(i);
    if (static_cast<int>(val.domainIdentProc.proc()) == p_rank ||
        (enforceSearchResultSymmetry && static_cast<int>(val.rangeIdentProc.proc()) == p_rank)) {
      ++numLocal;
    }
    if (static_cast<int>(val.domainIdentProc.proc()) != p_rank) {
      CommBuffer & buf = commSparse.send_buffer(val.domainIdentProc.proc());
      buf.skip<ValueType>(1);
    }
    if (enforceSearchResultSymmetry) {
      if (static_cast<int>(val.rangeIdentProc.proc()) != p_rank && val.rangeIdentProc.proc() != val.domainIdentProc.proc()) {
        CommBuffer & buf = commSparse.send_buffer(val.rangeIdentProc.proc());
        buf.skip<ValueType>(1);
      }
    }
  }

  commSparse.allocate_buffers();

  for (unsigned i = 0; i != search_relations.extent(0); ++i) {
    const ValueType & val = search_relations(i);
    if (static_cast<int>(val.domainIdentProc.proc()) != p_rank) {
      CommBuffer & buf = commSparse.send_buffer(val.domainIdentProc.proc());
      buf.pack<ValueType>(val);
    }
    if (enforceSearchResultSymmetry) {
      if (static_cast<int>(val.rangeIdentProc.proc()) != p_rank && val.rangeIdentProc.proc() != val.domainIdentProc.proc()) {
        CommBuffer & buf = commSparse.send_buffer(val.rangeIdentProc.proc());
        buf.pack<ValueType>(val);
      }
    }
  }

  commSparse.communicate();

  size_t numRecvd = 0;
  for (int p = 0; p < p_size; ++p) {
    CommBuffer & buf = commSparse.recv_buffer(p);
    numRecvd += (buf.remaining()/sizeof(ValueType));
  }

  auto keep = search_relations.extent(0);
  Kokkos::resize(Kokkos::WithoutInitializing, search_relations, numLocal+numRecvd);

  for (int p = 0; p < p_size; ++p) {
    CommBuffer & buf = commSparse.recv_buffer(p);
    while (buf.remaining()) {
      ValueType val;
      buf.unpack<ValueType>(val);
      search_relations(keep++) = val;
    }
  }
}

namespace impl {
template <typename T>
bool constexpr is_stk_box =
    std::is_same_v<T, Box<typename T::value_type>> || std::is_base_of_v<Box<typename T::value_type>, T>;

template <typename T>
bool constexpr is_stk_sphere =
    std::is_same_v<T, Sphere<typename T::value_type>> || std::is_base_of_v<Sphere<typename T::value_type>, T>;

template <typename T>
bool constexpr is_stk_point =
    std::is_same_v<T, Point<typename T::value_type>> || std::is_base_of_v<Point<typename T::value_type>, T>;
}

template <typename SearchResultsType, typename ExecutionSpace>
class SearchResultCommunication
{
    using ValueType = typename SearchResultsType::value_type;
    using DeviceBufferAppender = impl::DeviceMPIBufferAppender<ValueType, ExecutionSpace>;
    using DeviceBuffers = impl::DeviceMPIBuffers<ValueType, ExecutionSpace>;

  public:
    struct FillBuffer {};
    struct UnpackBuffer {};


    //Note: unlike communicate_views(), this assumes search_relations is
    //      exactly the size it needs to be (no unused entries at the end)
    SearchResultCommunication(stk::ParallelMachine comm,
                              SearchResultsType searchResults,
                              ExecutionSpace execSpace,
                              bool enforceSearchResultSymmetry = true) :
      m_comm(comm),
      m_searchResults(searchResults),
      m_execSpace(execSpace),
      m_enforceSearchResultSymmetry(enforceSearchResultSymmetry),
      m_bufferAppender(stk::parallel_machine_size(comm), execSpace),
      m_commRank(stk::parallel_machine_rank(comm)),
      m_numLocalResults(searchResults.size())
    {}

    SearchResultsType run()
    {
      Kokkos::Profiling::pushRegion("Filling MPI Buffers");
      Kokkos::RangePolicy<FillBuffer, ExecutionSpace> packPolicy(m_execSpace, 0, m_searchResults.extent(0));
      Kokkos::parallel_for("mpi_buffer_sizing", packPolicy, *this);
      m_execSpace.fence();

      m_bufferAppender.allocate_buffers();

      Kokkos::parallel_for("mpi_buffer_fill", packPolicy, *this);
      m_execSpace.fence();
      Kokkos::Profiling::popRegion();

      Kokkos::Profiling::pushRegion("Parallel Communication");
      impl::DeviceDataExchangeUnknownPattern<ValueType, ExecutionSpace> exchanger(m_bufferAppender.getBuffers(), m_execSpace, m_comm);
      m_recvBuffers = exchanger.communicate();
      Kokkos::Profiling::popRegion();

      Kokkos::Profiling::pushRegion("Unpacking buffers");
      Kokkos::resize(Kokkos::WithoutInitializing, m_searchResults, m_numLocalResults + m_recvBuffers.buffers.extent(0));
      Kokkos::RangePolicy<UnpackBuffer, ExecutionSpace> unpackPolicy(m_execSpace, 0, m_recvBuffers.buffers.extent(0));
      Kokkos::parallel_for("mpi_buffer_unpack", unpackPolicy, *this);
      m_execSpace.fence();
      Kokkos::Profiling::popRegion();

      return m_searchResults;
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(FillBuffer /*tag*/, int idx) const
    {
      ValueType searchResult = m_searchResults(idx);
      int domainProc = static_cast<int>(searchResult.domainIdentProc.proc());

      if (domainProc != m_commRank) {
        m_bufferAppender.push_back(domainProc, searchResult);
      }

      if (m_enforceSearchResultSymmetry)
      {
        int rangeProc = static_cast<int>(searchResult.rangeIdentProc.proc());
        if (rangeProc != m_commRank &&
            rangeProc != domainProc) {
          m_bufferAppender.push_back(rangeProc, searchResult);
        }
      }
    }

    KOKKOS_INLINE_FUNCTION
    void operator()(UnpackBuffer /*tag*/, int idx) const
    {
      m_searchResults(idx + m_numLocalResults) = m_recvBuffers.buffers(idx);
    }


  private:
    stk::ParallelMachine m_comm;
    SearchResultsType m_searchResults;
    ExecutionSpace m_execSpace;
    bool m_enforceSearchResultSymmetry;

    DeviceBufferAppender m_bufferAppender;
    int m_commRank;
    int m_numLocalResults;
    DeviceBuffers m_recvBuffers;
};

struct SearchData
{
  std::any data;
};

} // end namespace stk::search

#endif
