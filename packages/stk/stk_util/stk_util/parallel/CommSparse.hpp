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

#ifndef stk_util_parallel_CommSparse_hpp
#define stk_util_parallel_CommSparse_hpp

#include <cstddef>                      // for size_t, ptrdiff_t
#include <vector>
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine
#include <stk_util/parallel/ParallelComm.hpp>  // for CommBuffer etc

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

  /** Obtain the message buffer for a given processor */
  const CommBuffer & recv_buffer( int p ) const
  {
#ifndef NDEBUG
    if ( m_size <= p ) { rank_error("recv_buffer",p); }
#endif
    return m_recv[p] ;
  }

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
  explicit CommSparse( ParallelMachine comm)
    : m_comm( comm ),
      m_size( parallel_machine_size( comm ) ),
      m_rank( parallel_machine_rank( comm ) ),
      m_send(m_size),
      m_recv(m_size),
      m_send_data(),
      m_recv_data(),
      m_send_procs(),
      m_recv_procs()
  {
  }

  /** Allocate communication buffers based upon
   *  sizing from the surrogate send buffer packing.
   *  Returns true if the local processor is actually
   *  sending or receiving.
   */
  bool allocate_buffers();

  /** Allocate communication buffers based upon
   *  sizing from the surrogate send buffer packing, with user-specified
   *  vectors of procs to send and receive with. Allowing the send/recv
   *  procs to be specified allows the CommSparse class to avoid a
   *  round of communication.
   */
  void allocate_buffers(const std::vector<int>& send_procs, const std::vector<int>& recv_procs);

  /** Communicate send buffers to receive buffers.  */
  void communicate();

  /** Swap send and receive buffers leading to reversed communication. */
  void swap_send_recv();

  /** Reset, but do not reallocate, message buffers for reprocessing.
   *  Sets 'size() == 0' and 'remaining() == capacity()'.
   */
  void reset_buffers();

  ~CommSparse()
  {
    m_comm = parallel_machine_null();
    m_size = 0 ;
    m_rank = 0 ;
    m_send.clear();
    m_recv.clear();
  }
private:

  /** Construct for undefined communication.
   *  No buffers are allocated.
   */
  CommSparse()
    : m_comm( parallel_machine_null() ),
      m_size( 0 ), 
      m_rank( 0 ),
      m_send(),
      m_recv(),
      m_send_data(),
      m_recv_data(),
      m_send_procs(),
      m_recv_procs()
  {}

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

template<typename COMM, typename PACK_ALGORITHM>
void pack_and_communicate(COMM & comm, const PACK_ALGORITHM & algorithm)
{
    algorithm();
    const bool actuallySendingOrReceiving = comm.allocate_buffers();
    if (actuallySendingOrReceiving) {
        algorithm();
        comm.communicate();
    }
}

template<typename COMM, typename UNPACK_ALGORITHM>
void unpack_communications(COMM & comm, const UNPACK_ALGORITHM & algorithm)
{
    for(int proc_id=0; proc_id<comm.parallel_size(); ++proc_id)
    {
        if (proc_id != comm.parallel_rank())
        {
            while(comm.recv_buffer(proc_id).remaining())
            {
                algorithm(proc_id);
            }
        }
    }
}

template <typename T>
void pack_vector_to_proc(stk::CommSparse& comm, const T& data, int otherProc)
{
    comm.send_buffer(otherProc).pack<unsigned>(data.size());
    for(size_t i=0; i<data.size(); ++i)
        comm.send_buffer(otherProc).pack<typename T::value_type>(data[i]);
}

template <typename T>
void unpack_vector_from_proc(stk::CommSparse& comm, T& data, int fromProc)
{
    unsigned num_items = 0;
    comm.recv_buffer(fromProc).unpack<unsigned>(num_items);
    data.resize(num_items);
    for(unsigned i=0;i<num_items;++i)
        comm.recv_buffer(fromProc).unpack<typename T::value_type>(data[i]);
}

}

//----------------------------------------------------------------------

#endif

