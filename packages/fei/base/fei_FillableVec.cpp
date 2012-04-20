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


#include "fei_FillableVec.hpp"

namespace fei {

FillableVec::FillableVec()
 : vecdata_()
{
}

FillableVec::~FillableVec()
{
}

void
FillableVec::addEntries(unsigned numEntries,
                  const double* coefs,
                  const int* indices)
{
  for(unsigned i=0; i<numEntries; ++i) {
    addEntry(indices[i], coefs[i]);
  }
}

void
FillableVec::putEntries(unsigned numEntries,
                  const double* coefs,
                  const int* indices)
{
  for(unsigned i=0; i<numEntries; ++i) {
    putEntry(indices[i], coefs[i]);
  }
}

void
FillableVec::setValues(double value)
{
  if (size() > 0) {
    feipoolmap::iterator iter = vecdata_.begin(), iter_end = vecdata_.end();
    for(; iter != iter_end; ++iter) {
      iter->second = value;
    }
  }
}

void
FillableVec::clear()
{
  vecdata_.clear();
}

unsigned
FillableVec::size() const
{
  return vecdata_.size();
}

bool
FillableVec::hasEntry(int index) const
{
  feipoolmap::const_iterator iter = find(index);
  return iter != vecdata_.end();
}

double
FillableVec::getEntry(int index) const
{
  feipoolmap::const_iterator iter = find(index);
  if (iter == vecdata_.end()) {
    throw std::runtime_error("FillableVec::getEntry: index not found.");
  }

  return iter->second;
}

void
FillableVec::removeEntry(int index)
{
  feipoolmap::iterator iter = find(index);
  if (iter != vecdata_.end()) {
    vecdata_.erase(iter);
  }
}

}//namespace fei

