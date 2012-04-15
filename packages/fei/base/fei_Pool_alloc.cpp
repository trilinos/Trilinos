/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/


#include "fei_macros.hpp"

#include "fei_Pool_alloc.hpp"
#include "fei_Pool.hpp"

#include <cstdlib>
#include <map>

struct fei_map_holder {
  std::map<size_t,fei_Pool*> fei_mem_pools;

  ~fei_map_holder()
  {
    std::map<size_t,fei_Pool*>::iterator
      iter = fei_mem_pools.begin(),
      iter_end = fei_mem_pools.end();
    for(; iter != iter_end; ++iter) {
      delete iter->second;
    }
  }
};

static fei_map_holder fmh;

fei_Pool* get_fei_mem_pool(size_t n)
{
  std::map<size_t,fei_Pool*>::iterator
    iter = fmh.fei_mem_pools.find(n);

  if (iter == fmh.fei_mem_pools.end()) {
    fei_Pool* new_pool = new fei_Pool(n);
    fmh.fei_mem_pools.insert(std::make_pair(n,new_pool));
    return(new_pool);
  }

  return(iter->second);
}

