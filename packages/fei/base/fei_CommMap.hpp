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

#ifndef _fei_CommMap_hpp_
#define _fei_CommMap_hpp_

#include <fei_macros.hpp>
#include <fei_ArrayUtils.hpp>
#include <set>
#include <map>

namespace fei {

/** Map that maps processors to associated vectors of items to be sent or recv'd. */
template<typename T>
struct CommMap {
  typedef std::map<int,std::vector<T> > Type;
};

/** Given a proc and an array of items, add the mapping
       proc -> items
   to the given comm_map.
   Optionally ensure that the comm_map's vector of items for proc remains
   sorted and unique.
*/
template<typename T>
void addItemsToCommMap(int proc, size_t numItems, const T* items,
                       typename CommMap<T>::Type& comm_map,
                       bool keep_sorted_and_unique = true)
{
  typename CommMap<T>::Type::iterator iter = comm_map.find(proc);
  if (iter == comm_map.end()) {
    iter = comm_map.insert(std::make_pair(proc,std::vector<T>())).first;
  }

  std::vector<T>& comm_items = iter->second;

  if (keep_sorted_and_unique) {
    for(size_t i=0; i<numItems; ++i) {
      fei::sortedListInsert(items[i], comm_items);
    }
  }
  else {
    for(size_t i=0; i<numItems; ++i) {
      comm_items.push_back(items[i]);
    }
  }
}

} //namespace fei

#endif // _fei_CommMap_hpp_

