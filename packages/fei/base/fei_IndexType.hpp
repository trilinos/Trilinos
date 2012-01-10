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


#ifndef _fei_IndexType_hpp_
#define _fei_IndexType_hpp_

#include <fei_macros.hpp>

#include <map>
#include <algorithm>
#include <functional>

#ifdef HAVE_FEI_BOOST
#include <boost/unordered_map.hpp>
#endif


namespace fei {

template< typename T_Key, typename T_Value>
class IndexType : public
#ifdef HAVE_FEI_BOOST
   boost::unordered_map<T_Key, T_Value>
#else
  std::map<T_Key, T_Value>
#endif
{
  typedef std::map<T_Key, T_Value> mapRep;
public:
#ifdef HAVE_FEI_BOOST
public:
  typedef boost::unordered_map<T_Key, T_Value> base;
  bool isStdMap() const {return false;}
  const mapRep& asMap(mapRep &map_rep) const {
    resyncToMap(map_rep);
    return map_rep;
  }
  mapRep& asMap(mapRep &map_rep) {
    resyncToMap(map_rep);
    return map_rep;
  }

  T_Key getMaxKey() const {
    if (this->size() == 0 ) 
      return T_Value(0);
    T_Key k = this->begin()->first;
    typename base::const_iterator iter=++(this->begin()),end=this->end();
    for ( ; iter != end; ++iter) 
      k = std::max(k,iter->first);
    return k;
  }

  T_Key getMinKey() const {
    if (this->size() == 0 ) 
      return T_Value(0);
    T_Key k = this->begin()->first;
    typename base::const_iterator iter=++(this->begin()),end=this->end();
    for ( ; iter != end; ++iter) 
      k = std::min(k,iter->first);
    return k;
  }

#else
  typedef std::map<T_Key, T_Value> base;
  bool isStdMap() const {return true;}
  mapRep &asMap(mapRep &rep) const{ return *static_cast<base *>(&rep);}
  const mapRep &asMap(const mapRep&rep) const { return *static_cast<const base *>( &rep);}

  T_Key getMaxKey() const {
    IndexType<int,int>::const_iterator highest = this->end();
    if (this->size() > 0) 
      return (--highest)->first;
    return T_Key(0);
  }

  T_Key getMinKey() const {
    if (this->size() > 0) 
      return this->begin()->first;
    return T_Key(0);
  }

#endif
  
  bool lightWeightCompare(const mapRep &map_rep) const {
    return (this->size() == map_rep.size());
  }

  bool heavyWeightCompare(const mapRep &map_rep) const {

    if (this->size() != map_rep.size())
      return false;
    typename base::const_iterator iter=this->begin(),end=this->end();
    for ( ; iter != end; ++iter) {
      typename mapRep::const_iterator lookup = map_rep.find(iter->first);
      if ( lookup == map_rep.end())
	return false;
      if (lookup->second  != iter->second)
	return false;
    }
    return true;
  }

  void resyncToMap( mapRep &map_rep) const {
    if ( !heavyWeightCompare(map_rep) ) {
      map_rep.clear();
      typename base::const_iterator iter=this->begin(), end=this->end();
      for ( ; iter != end; ++iter) 
	map_rep.insert(*iter);
    }
  }

  void resyncFromMap(const mapRep &map_rep) {
    if ( !heavyWeightCompare(map_rep) ) {
      this->clear();
      typename mapRep::const_iterator iter=map_rep.begin(), end=map_rep.end();
      for ( ; iter != end; ++iter) 
	this->insert(*iter);
    }
  }  
}; 
}
#endif
