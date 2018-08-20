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

#ifndef stk_util_CommNeighbors_hpp
#define stk_util_CommNeighbors_hpp

#include <cstddef>                      // for size_t, ptrdiff_t
#include <vector>
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine
#include <stk_util/parallel/CommBufferV.hpp>

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
#ifndef NDEBUG
    if ( m_size <= p ) { rank_error("send_buffer",p); }
#endif
    return m_send[p] ;
  }

  /** Obtain the message buffer for a given processor */
  CommBufferV & recv_buffer( int p )
  {
#ifndef NDEBUG
    if ( m_size <= p ) { rank_error("recv_buffer",p); }
#endif
    return m_recv[p] ;
  }

  //----------------------------------------
  /** Construct for a to-be-sized communication.
   *  Comm scenario:
   *  1) constructor argument 'neighbor_procs' specifies the processors
   *     which may be communicated with (neighbors==send-procs==recv-procs).
   *  2) send-buffers are packed with data to be sent
   *     All processors sent to, must be members of neighbor_procs.
   *  3) communicate() performs the communication and stores recvd data
   *     in recv_buffers.
   *     All processors recvd from, are members of neighbor_procs.
   */
  CommNeighbors( stk::ParallelMachine comm, const std::vector<int>& neighbor_procs );
  CommNeighbors( stk::ParallelMachine comm, const std::vector<int>& send_procs, const std::vector<int>& recv_procs );

  //----------------------------------------
  /** Communicate send buffers to receive buffers.  */
  void communicate();

  /** Reset, but do not reallocate, message buffers for reprocessing.
   *  Sets 'size() == 0' and 'remaining() == capacity()'.
   */
  void reset_buffers();

  ~CommNeighbors();

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

  void rank_error( const char * , int ) const ;
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

