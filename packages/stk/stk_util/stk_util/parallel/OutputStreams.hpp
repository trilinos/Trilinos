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

#ifndef STK_UTIL_PARALLEL_OutputStreams_hpp
#define STK_UTIL_PARALLEL_OutputStreams_hpp

#include <stk_util/parallel/Parallel.hpp>
#include <iosfwd>

namespace stk {

// outputP0() defaults to std::cout on MPI rank 0 (of comm-world)
//            and a 'null' stream (which discards everything) on other ranks.
//
// If default stream (or comm-world) is not suitable, it can be altered
// with set_outputP0(..) below.
//
std::ostream& outputP0();

std::ostream& outputNull();

// output() is an ostream that holds output until output_flush() is
// called.
//
std::ostream& output();

// output_flush() sends all output in the 'output()' stream on each
// mpi rank to 'outputP0()'.
// output_flush() must be called on each mpi rank in the communicator
//
void output_flush();

void set_outputP0(std::ostream* ostreamPtr, ParallelMachine comm = MPI_COMM_WORLD); // CHECK: ALLOW MPI_COMM_WORLD

// if you fear that the output streams are in a bad state,
// this resets them to their defaults
void reset_default_output_streams(ParallelMachine comm = MPI_COMM_WORLD); // CHECK: ALLOW MPI_COMM_WORLD

}

#endif // STK_UTIL_PARALLEL_OutputStreams_hpp

