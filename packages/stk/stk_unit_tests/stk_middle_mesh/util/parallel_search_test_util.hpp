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

#ifndef  STK_MIDDLE_MESH_PARALLEL_SEARCH_TEST_UTIL_HPP
#define  STK_MIDDLE_MESH_PARALLEL_SEARCH_TEST_UTIL_HPP

#include "stk_middle_mesh/bounding_box_search.hpp"

namespace {

using namespace stk::middle_mesh::search;

class SplitCommTestUtil {

public:
  SplitCommTestUtil(int nProcsSendMesh, int nProcsRecvMesh)
  : m_numProcsSendMesh(nProcsSendMesh), m_numProcsRecvMesh(nProcsRecvMesh)
  {
    m_color = get_color();

    MPI_Comm_split(MPI_COMM_WORLD, static_cast<int>(m_color), 0, &m_splitComm);
  }

  ~SplitCommTestUtil() {
    MPI_Comm_free(&m_splitComm);
  }

  SplitCommTestUtil() = delete;

  SplitCommColor get_color() const
  {
    int myRankWorld;
    MPI_Comm_rank(MPI_COMM_WORLD, &myRankWorld);

    SplitCommColor color = myRankWorld < m_numProcsSendMesh ? SplitCommColor::SEND : SplitCommColor::RECV;

    return color;
  }

  bool get_status() const
  {
    int commSizeWorld;
    MPI_Comm_size(MPI_COMM_WORLD, &commSizeWorld);

    if (commSizeWorld != (m_numProcsSendMesh + m_numProcsRecvMesh))
      return false;

    return true;
  }

  MPI_Comm get_comm() const
  {
    return m_splitComm;
  }

  int get_local_comm_rank(int globalRank) const
  {
    int commSizeWorld, commSizeSplit;
    MPI_Comm_size(MPI_COMM_WORLD, &commSizeWorld);
    MPI_Comm_size(m_splitComm, &commSizeSplit);

    int commSizeSend = (m_color == SplitCommColor::SEND) ? commSizeSplit : (commSizeWorld - commSizeSplit);

    int localRank = (globalRank < commSizeSend) ? globalRank : (globalRank - commSizeSend);

    return localRank;
  }

  int get_global_comm_rank(int localRankInOtherComm) const
  {
    int commSizeWorld, commSizeSplit;
    MPI_Comm_size(MPI_COMM_WORLD, &commSizeWorld);
    MPI_Comm_size(m_splitComm, &commSizeSplit);

    int commSizeSend = (m_color == SplitCommColor::SEND) ? commSizeSplit : (commSizeWorld - commSizeSplit);

    int globalRank = (m_color == SplitCommColor::SEND) ? (commSizeSend + localRankInOtherComm) : localRankInOtherComm;

    return globalRank;
  }

private:
  int m_numProcsSendMesh;
  int m_numProcsRecvMesh;
  SplitCommColor m_color;
  MPI_Comm m_splitComm;
};

}

#endif
