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


#ifndef _fei_FillableVec_hpp_
#define _fei_FillableVec_hpp_

#include "fei_macros.hpp"
#include "fei_Pool_alloc.hpp"
#include <map>

namespace fei {

class FillableVec {
 public:
  FillableVec();
  virtual ~FillableVec();

  void addEntry(int index, double coef);

  void addEntries(unsigned numEntries,
                  const double* coefs,
                  const int* indices);

  void putEntry(int index, double coef);

  void putEntries(unsigned numEntries,
                  const double* coefs,
                  const int* indices);

  void setValues(double value);

  unsigned size() const;

  bool hasEntry(int index) const;

  /** Return coef for index if index is present, otherwise throw.
  */
  double getEntry(int index) const;

  void removeEntry(int index);

  void clear();

  typedef
    std::map<int,double,std::less<int>,
     fei_Pool_alloc<std::pair<const int,double> > > feipoolmap;

  typedef feipoolmap::iterator iterator;
  typedef feipoolmap::const_iterator const_iterator;

  iterator begin() {return vecdata_.begin();}
  iterator end() {return vecdata_.end();}

  const_iterator begin() const {return vecdata_.begin();}
  const_iterator end() const {return vecdata_.end();}

  bool operator==(const FillableVec& rhs) const;

  bool operator!=(const FillableVec& rhs) const;

 private:
  feipoolmap vecdata_;
};//class FillableVec

inline void
FillableVec::addEntry(int index, double coef)
{
  feipoolmap::iterator iter = vecdata_.lower_bound(index);
  if (iter == vecdata_.end() || iter->first != index) {
    vecdata_.insert(iter, std::make_pair(index, coef));
  }
  else {
    iter->second += coef;
  }
}

inline void
FillableVec::putEntry(int index, double coef)
{
  feipoolmap::iterator iter = vecdata_.lower_bound(index);
  if (iter == vecdata_.end() || iter->first != index) {
    vecdata_.insert(iter, std::make_pair(index, coef));
  }
  else {
    iter->second = coef;
  }
}

inline
bool
FillableVec::operator==(const FillableVec& rhs) const
{
  return vecdata_ == rhs.vecdata_;
}

inline
bool
FillableVec::operator!=(const FillableVec& rhs) const
{
  return vecdata_ != rhs.vecdata_;
}

}//namespace fei

#endif

