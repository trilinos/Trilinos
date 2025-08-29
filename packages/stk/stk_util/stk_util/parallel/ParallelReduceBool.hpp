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

#ifndef stk_util_parallel_ParallelReduceBool_hpp
#define stk_util_parallel_ParallelReduceBool_hpp

#include "stk_util/parallel/ParallelReduce.hpp"  // for ParallelMachine, etc
#include "stk_util/parallel/CouplingVersions.hpp"

//------------------------------------------------------------------------

namespace stk {

inline bool is_true_on_all_procs(ParallelMachine comm , const bool truthValue)
{
#ifdef STK_HAS_MPI
  int truthValueInt = truthValue, globalResult;
  MPI_Allreduce(&truthValueInt, &globalResult, 1, MPI_INT, MPI_LAND, comm);
  return globalResult;
#else
  return truthValue;
#endif
}

inline bool is_true_on_any_proc(ParallelMachine comm , const bool truthValue)
{
#ifdef STK_HAS_MPI
  int truthValueInt = truthValue, globalResult;
  MPI_Allreduce(&truthValueInt, &globalResult, 1, MPI_INT, MPI_LOR, comm);
  return globalResult;
#else
  return truthValue;
#endif
}

}

#endif
