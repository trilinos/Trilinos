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

#ifndef stk_util_parallel_CommSparse_hpp
#define stk_util_parallel_CommSparse_hpp

#include <cstddef>                      // for size_t, ptrdiff_t
#include <vector>
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine
#include <stk_util/parallel/ParallelComm.hpp>  // for CommBuffer etc
#include <boost/static_assert.hpp>

//------------------------------------------------------------------------

namespace stk {

/** Given the send sizes determine the receive sizes and also
 * set output vectors of send-procs and recv-procs.
 * Send and receive size arrays are dimensioned to
 * the size of the parallel machine. (i.e., the number
 * of MPI processor ranks.)
 * Output vectors for send-procs and recv-procs will have
 * length num-send-procs and num-recv-procs respectively.
 */
void comm_recv_procs_and_msg_sizes(ParallelMachine comm,
                     const unsigned * const send_size,
                     unsigned * const recv_size,
                     std::vector<int>& output_send_procs,
                     std::vector<int>& output_recv_procs);

void comm_recv_procs_and_msg_sizes(ParallelMachine comm ,
                                   const std::vector<CommBuffer>& send_bufs ,
                                         std::vector<CommBuffer>& recv_bufs,
                                   std::vector<int>& send_procs,
                                   std::vector<int>& recv_procs);

/** Given send sizes (of length number-of-MPI-processor-ranks) and
 * send-procs and recv-procs (of length number-of-procs-to-send/recv-with),
 * set recv sizes (recv_size array has length number-of-MPI-processor-ranks).
 */
void comm_recv_msg_sizes(ParallelMachine comm ,
                     const unsigned * const send_size ,
                     const std::vector<int>& send_procs,
                     const std::vector<int>& recv_procs,
                     unsigned * const recv_size);

void comm_recv_msg_sizes(ParallelMachine comm ,
                     const std::vector<int>& send_procs,
                     const std::vector<int>& recv_procs,
                     const std::vector<CommBuffer>& send_bufs,
                     std::vector<CommBuffer>& recv_bufs);

class CommSparse {
public:

  ParallelMachine parallel()      const { return m_comm ; }
  int             parallel_size() const { return m_size ; }
  int             parallel_rank() const { return m_rank ; }

  /** Obtain the message buffer for a given processor */
  CommBuffer & send_buffer( int p )
  {
#ifndef NDEBUG
    if ( m_size <= p ) { rank_error("send_buffer",p); }
#endif
    return m_send[p] ;
  }

  /** Obtain the message buffer for a given processor */
  CommBuffer & recv_buffer( int p )
  {
#ifndef NDEBUG
    if ( m_size <= p ) { rank_error("recv_buffer",p); }
#endif
    return m_recv[p] ;
  }

  //----------------------------------------
  /** Construct for undefined communication.
   *  No buffers are allocated.
   */
  CommSparse();
 
  //----------------------------------------
  /** Construct for a to-be-sized communication.
   *  Allocate surrogate send buffers to enable
   *  no-op packing for the purpose of send sizing.
   *  Surrogate send scenario:
   *  1) Surrogate send buffers are "packed" for sizing where
   *     packing sizes are recorded but no data is copied.
   *  2) 'allocate_buffers()' is called to allocate
   *     buffers.
   *  3) Send buffers are identically packed; however, this
   *     packing copies data into the send buffers.
   */
  explicit CommSparse( ParallelMachine );

  /** Allocate communication buffers based upon
   *  sizing from the surrogate send buffer packing.
   */
  void allocate_buffers();
  /** Allocate communication buffers based upon
   *  sizing from the surrogate send buffer packing, with user-specified
   *  vectors of procs to send and receive with. Allowing the send/recv
   *  procs to be specified allows the CommSparse class to avoid a
   *  call to MPI_Reduce_scatter.
   */
  void allocate_buffers(const std::vector<int>& send_procs, const std::vector<int>& recv_procs);

  //----------------------------------------
  /** Communicate send buffers to receive buffers.  */
  void communicate();

  //----------------------------------------
  /** Swap send and receive buffers leading to reversed communication. */
  void swap_send_recv();

  /** Reset, but do not reallocate, message buffers for reprocessing.
   *  Sets 'size() == 0' and 'remaining() == capacity()'.
   */
  void reset_buffers();

  ~CommSparse();

private:

  CommSparse( const CommAll & );
  CommSparse & operator = ( const CommAll & );

  void rank_error( const char * , int ) const ;

  void allocate_data(std::vector<CommBuffer>& bufs, std::vector<unsigned char>& data);

  ParallelMachine m_comm ;
  int             m_size ;
  int             m_rank ;
  std::vector<CommBuffer> m_send;
  std::vector<CommBuffer> m_recv;
  std::vector<unsigned char> m_send_data;
  std::vector<unsigned char> m_recv_data;
  std::vector<int> m_send_procs;
  std::vector<int> m_recv_procs;
};

}

//----------------------------------------------------------------------

#endif

