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

#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/parallel/ParallelComm.hpp>
#include <sstream>                      // for basic_ostream::operator<<, etc
#include <stdexcept>                    // for runtime_error
#include <vector>                       // for vector
#include "stk_util/parallel/Parallel.hpp"  // for ParallelMachine, etc
#include "stk_util/stk_config.h"        // for STK_HAS_MPI

namespace stk {

#if defined( STK_HAS_MPI )

/// @todo REFACTOR The per-processor separator needs to be addressed.  Already created a duplicate
/// for deferred error messages which could be replacement.

void all_write_string( ParallelMachine arg_comm ,
                       std::ostream &  arg_root_os ,
                       const std::string & arg_msg )
{
  const int i_zero = 0 ;
  const int p_root = 0 ;
  const unsigned p_size = parallel_machine_size( arg_comm );
  const unsigned p_rank = parallel_machine_rank( arg_comm );

  int result = MPI_SUCCESS ;

  // Gather the send counts on root processor

  int send_count = arg_msg.size();

  std::vector<int> recv_count( p_size , i_zero );

  int * const recv_count_ptr = & recv_count[0] ;

  result = MPI_Gather( & send_count , 1 , MPI_INT ,
                       recv_count_ptr , 1 , MPI_INT ,
                       p_root , arg_comm );

  if ( MPI_SUCCESS != result ) {
    std::ostringstream msg ;
    msg << "stk::all_write FAILED: MPI_Gather = " << result ;
    throw std::runtime_error( msg.str() );
  }

  // Receive counts are only non-zero on the root processor:

  std::vector<int> recv_displ( p_size + 1 , i_zero );

  for ( unsigned i = 0 ; i < p_size ; ++i ) {
    recv_displ[i+1] = recv_displ[i] + recv_count[i] ;
  }

  const unsigned recv_size = static_cast<unsigned>(recv_displ[ p_size ]);

  std::vector<char> buffer( recv_size );

  {
    const char * const send_ptr = arg_msg.c_str();
    char * const recv_ptr = recv_size ? & buffer[0] : NULL ;
    int * const recv_displ_ptr = & recv_displ[0] ;

    result = MPI_Gatherv( const_cast<char*>(send_ptr), send_count, MPI_CHAR ,
                          recv_ptr, recv_count_ptr, recv_displ_ptr, MPI_CHAR,
                          p_root, arg_comm );
  }

  if ( MPI_SUCCESS != result ) {
    std::ostringstream msg ;
    msg << "stk::all_write FAILED: MPI_Gatherv = " << result ;
    throw std::runtime_error( msg.str() );
  }

  if ( p_root == static_cast<int>(p_rank) ) {
//    arg_root_os << std::endl ;
    for ( unsigned i = 0 ; i < p_size ; ++i ) {
      if ( recv_count[i] ) {
        char * const ptr = & buffer[ recv_displ[i] ];
        arg_root_os.write( ptr , recv_count[i] );
        arg_root_os << std::endl ;
      }
    }
    arg_root_os.flush();
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------

void all_reduce( ParallelMachine  arg_comm ,
                 ParallelReduceOp arg_op ,
                 void           * arg_in ,
                 void           * arg_out ,
                 unsigned         arg_len )
{
  MPI_Op mpi_op = MPI_OP_NULL ;

  MPI_Op_create( arg_op , 0 , & mpi_op );

  // The SUN was buggy when combining an
  // MPI_Allreduce with a user defined operator,
  // use reduce/broadcast instead.
/*
  const int result =
    MPI_Allreduce(arg_in,arg_out,arg_len,MPI_BYTE,mpi_op,arg_comm);
*/

  const int result_reduce =
    MPI_Reduce(arg_in,arg_out,arg_len,MPI_BYTE,mpi_op,0,arg_comm);

  const int result_bcast =
    MPI_Bcast(arg_out,arg_len,MPI_BYTE,0,arg_comm);

  MPI_Op_free( & mpi_op );

  if ( MPI_SUCCESS != result_reduce || MPI_SUCCESS != result_bcast ) {
    std::ostringstream msg ;
    msg << "stk::all_reduce FAILED: MPI_Reduce = " << result_reduce
        << " MPI_Bcast = " << result_bcast ;
    throw std::runtime_error( msg.str() );
  }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
void all_reduce_impl( ParallelMachine comm ,
                     const size_t * local , size_t * global , unsigned count, MPI_Op op)
{
  size_t * tmp = const_cast<size_t*>( local );

  if ( sizeof(size_t) == sizeof(unsigned) ) {
    MPI_Allreduce( tmp , global , count , MPI_UNSIGNED , op , comm );

  }
  else if ( sizeof(size_t) == sizeof(unsigned long) ) {
    MPI_Allreduce( tmp , global , count , MPI_UNSIGNED_LONG , op , comm );
  }
  else {
    unsigned long * const in  = new unsigned long[ count ];
    unsigned long * const out = new unsigned long[ count ];

    for ( unsigned i = 0 ; i < count ; ++i ) { in[i] = local[i] ; }
    MPI_Allreduce( in , out , count , MPI_UNSIGNED_LONG , op , comm );
    for ( unsigned i = 0 ; i < count ; ++i ) { global[i] = out[i] ; }

    delete[] in ;
    delete[] out ;
  }
}

void all_reduce_bor( ParallelMachine comm ,
                     const unsigned * local ,
                     unsigned * global , unsigned count )
{
  unsigned * tmp = const_cast<unsigned*>( local );
  MPI_Allreduce( tmp , global , count , MPI_UNSIGNED , MPI_BOR , comm );
}
//----------------------------------------------------------------------
//----------------------------------------------------------------------

#else

void all_write_string( ParallelMachine ,
                       std::ostream &  arg_root_os ,
                       const std::string & arg_msg )
{
  arg_root_os << arg_msg ;
}

void all_reduce_bor( ParallelMachine ,
                     const unsigned * local ,
                     unsigned * global , unsigned count )
{
  for ( unsigned i = 0 ; i < count ; ++i ) { global[i] = local[i] ; }
}

//----------------------------------------------------------------------

void all_reduce( ParallelMachine ,
                 ParallelReduceOp ,
                 void   * arg_in ,
                 void   * arg_out ,
                 unsigned arg_len )
{
  unsigned char * i = reinterpret_cast<unsigned char *>( arg_in );
  unsigned char * o = reinterpret_cast<unsigned char *>( arg_out );
  for ( unsigned char * const e = i + arg_len ; e != i ; ++i , ++o ) {
    *o = *i ;
  }
}

#endif

}

