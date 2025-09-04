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

#ifndef stk_util_CommNeighbors_hpp
#define stk_util_CommNeighbors_hpp

#include "stk_util/stk_config.h"              // for STK_HAS_MPI
#include "stk_util/parallel/CommBufferV.hpp"  // for CommBufferV
#include "stk_util/parallel/Parallel.hpp"     // for ParallelMachine, OMPI_MAJOR_VERSION, ompi_c...
#include "stk_util/util/ReportHandler.hpp"
#include <vector>                             // for vector

//------------------------------------------------------------------------
//
#if defined( STK_HAS_MPI )

#ifdef STK_MPI_SUPPORTS_NEIGHBOR_COMM
#undef STK_MPI_SUPPORTS_NEIGHBOR_COMM
#endif

#if MPI_VERSION >= 3
#define STK_MPI_SUPPORTS_NEIGHBOR_COMM
#endif

#ifdef OMPI_MAJOR_VERSION
//OpenMPI 3.1.x seems to have a bug in the MPI_Neighbor* functions.
#if OMPI_MAJOR_VERSION == 3 && OMPI_MINOR_VERSION == 1
#undef STK_MPI_SUPPORTS_NEIGHBOR_COMM
#endif
//OpenMPI 2.x.y doesn't seem to support MPI_Neighbor* functions either...
#if OMPI_MAJOR_VERSION == 2
#undef STK_MPI_SUPPORTS_NEIGHBOR_COMM
#endif
#if OMPI_MAJOR_VERSION == 10
//This is the IBM Spectrum MPI, which also has slow MPI_Neighbor functions.
#undef STK_MPI_SUPPORTS_NEIGHBOR_COMM
#endif

#endif

//the MPI_Neighbor functions seem to be unacceptably slow with intel mpi
#ifdef I_MPI_VERSION
#undef STK_MPI_SUPPORTS_NEIGHBOR_COMM
#endif

// Needs to be reassessed with OneAPI?
#if (defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER))
#undef STK_MPI_SUPPORTS_NEIGHBOR_COMM
#endif

//Finally: if the user explicitly enables or disables mpi-neighbor-comm
//by defining one of the following macros (e.g., with cmake option or with
//-D on compile line etc), then they take precedence over anything that
//happened in the ifdef logic above.
//
#ifdef STK_DISABLE_MPI_NEIGHBOR_COMM
#undef STK_MPI_SUPPORTS_NEIGHBOR_COMM
#endif

#ifdef STK_ENABLE_MPI_NEIGHBOR_COMM
#define STK_MPI_SUPPORTS_NEIGHBOR_COMM
#endif

#endif

//------------------------------------------------------------------------

namespace stk {

class CommNeighbors {
public:

  stk::ParallelMachine parallel()      const { return m_comm ; }
  int             parallel_size() const { return m_size ; }
  int             parallel_rank() const { return m_rank ; }

  /** Obtain the message buffer for a given processor */
  CommBufferV & send_buffer( int p )
  {
    STK_ThrowAssertMsg(m_size > p, "CommNeighbors::send_buffer p="<<p<<" out of range.");
    return m_send[p] ;
  }

  /** Obtain the message buffer for a given processor */
  CommBufferV & recv_buffer( int p )
  {
    STK_ThrowAssertMsg(m_size > p, "CommNeighbors::recv_buffer p="<<p<<" out of range.");
    return m_recv[p] ;
  }

  //############################################################################
  //  Do not use CommNeighbors in new code!  It will be deprecated and removed
  //  soon.  Please use CommSparse as an equivalent capability.
  //############################################################################
  CommNeighbors( stk::ParallelMachine comm, const std::vector<int>& neighbor_procs );
  CommNeighbors( stk::ParallelMachine comm, const std::vector<int>& send_procs, const std::vector<int>& recv_procs );

  //----------------------------------------
  /** Communicate send buffers to receive buffers.  */
  void communicate();

  /** Reset, but do not reallocate, message buffers for reprocessing.
   *  Sets 'size() == 0' and 'remaining() == capacity()'.
   */
  void reset_buffers();

  /** send procs become recv procs, send buffers become recv buffers...
  */
  void swap_send_recv();

  virtual ~CommNeighbors();

  const std::vector<int>& send_procs() const { return m_send_procs; }
  const std::vector<int>& recv_procs() const { return m_recv_procs; }

protected:

  virtual stk::ParallelMachine setup_neighbor_comm(stk::ParallelMachine fullComm,
                                                  const std::vector<int>& sendProcs,
                                                  const std::vector<int>& recvProcs);

  virtual void perform_neighbor_communication(MPI_Comm neighborComm,
                                              const std::vector<unsigned char>& sendBuf,
                                              const std::vector<int>& sendCounts,
                                              const std::vector<int>& sendDispls,
                                                    std::vector<unsigned char>& recvBuf,
                                                    std::vector<int>& recvCounts,
                                                    std::vector<int>& recvDispls);
  //----------------------------------------
  /** default Constructor not allowed
   */
  CommNeighbors() = delete;

  CommNeighbors( const CommNeighbors & );
  CommNeighbors & operator = ( const CommNeighbors & );

  void sort_procs_and_resize_buffers();

  stk::ParallelMachine m_comm ;
  bool m_created_dist_graph;
  int             m_size ;
  int             m_rank ;
  std::vector<CommBufferV> m_send;
  std::vector<CommBufferV> m_recv;
  std::vector<unsigned char> m_send_data;
  std::vector<unsigned char> m_recv_data;
  std::vector<int> m_send_procs;
  std::vector<int> m_recv_procs;
};

}

//----------------------------------------------------------------------

#endif

