// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
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

#include <stk_util/parallel/Parallel.hpp>
#include "stk_util/stk_config.h"        // for STK_HAS_MPI
#if defined( STK_HAS_MPI )
#  include "mpi.h"                      // for MPI_COMM_NULL, etc
#endif

/*--------------------------------------------------------------------*/
/* Parallel operations */

#if defined( STK_HAS_MPI )

namespace stk {

int parallel_machine_size( ParallelMachine parallel_machine )
{
  int value = 1 ;
  if (parallel_machine != MPI_COMM_NULL) {
    if ( MPI_SUCCESS != MPI_Comm_size( parallel_machine , &value ) ) {
      value = 1 ;
    }
  }
  return value ;
}

int parallel_machine_rank( ParallelMachine parallel_machine )
{
  int value = 0 ;
  if (parallel_machine != MPI_COMM_NULL) {
    if ( MPI_SUCCESS != MPI_Comm_rank( parallel_machine , &value ) ) {
      value = 0 ;
    }
  }
  return value ;
}

void parallel_machine_barrier( ParallelMachine parallel_machine )
{
  if (parallel_machine != MPI_COMM_NULL) {
    MPI_Barrier( parallel_machine );
  }
}

}

#else

namespace stk {

int parallel_machine_size( ParallelMachine parallel_machine) { return 1 ; }

int parallel_machine_rank( ParallelMachine parallel_machine) { return 0 ; }

void parallel_machine_barrier( ParallelMachine parallel_machine) {}

}

#endif
