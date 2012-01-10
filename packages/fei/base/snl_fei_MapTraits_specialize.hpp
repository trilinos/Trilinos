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

#ifndef _snl_fei_MapTraits_specialize_hpp_
#define _snl_fei_MapTraits_specialize_hpp_


#include <fei_macros.hpp>

#include <snl_fei_MapTraits.hpp>
#include <snl_fei_MapContig.hpp>

namespace snl_fei {
#if defined(FEI_HASH_MAP) && defined(FEI_HASH_SET)
template<>
struct MapTraits<FEI_HASH_MAP<int,FEI_HASH_SET<int>*> > {
  static FEI_HASH_MAP<int,FEI_HASH_SET<int>*>::iterator
    lower_bound(FEI_HASH_MAP<int,FEI_HASH_SET<int>*>& map_obj,
                int item)
  { return( map_obj.find(item) ); }

  static void insert(FEI_HASH_MAP<int,FEI_HASH_SET<int>*>& map_obj,
                     FEI_HASH_MAP<int,FEI_HASH_SET<int>*>::iterator& pos,
                     FEI_HASH_MAP<int,FEI_HASH_SET<int>*>::value_type& val)
  { map_obj.insert(val); }
};
#endif

#if defined(FEI_HASH_MAP)

template<>
struct MapTraits<FEI_HASH_MAP<int,fei::ctg_set<int>*> >
{
  static FEI_HASH_MAP<int,fei::ctg_set<int>*>::iterator
    lower_bound(FEI_HASH_MAP<int,fei::ctg_set<int>*>& map_obj, int item)
  { return( map_obj.find(item) ); }

  static void insert(FEI_HASH_MAP<int,fei::ctg_set<int>*>& map_obj,
                     FEI_HASH_MAP<int,fei::ctg_set<int>*>::iterator& pos,
                     FEI_HASH_MAP<int,fei::ctg_set<int>*>::value_type& val)
 { map_obj.insert(val); }};

template<>
struct MapTraits<FEI_HASH_MAP<int,std::set<int>*> > {
  static FEI_HASH_MAP<int,std::set<int>*>::iterator
    lower_bound(FEI_HASH_MAP<int,std::set<int>*>& map_obj,
                int item)
  { return( map_obj.find(item) ); }

  static void insert(FEI_HASH_MAP<int,std::set<int>*>& map_obj,
                     FEI_HASH_MAP<int,std::set<int>*>::iterator& pos,
                     FEI_HASH_MAP<int,std::set<int>*>::value_type& val)
  { map_obj.insert(val); }
};
#endif

}//namespace snl_fei
#endif

