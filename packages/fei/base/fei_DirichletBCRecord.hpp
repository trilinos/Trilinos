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

#ifndef _fei_DirichletBCRecord_hpp_
#define _fei_DirichletBCRecord_hpp_


#include <fei_macros.hpp>

namespace fei {

struct DirichletBCRecord {
  int IDType;
  int ID;
  int fieldID;
  int whichComponent;
  double prescribedValue;

  bool operator!=(const DirichletBCRecord& rhs) const
  {
    return IDType != rhs.IDType || ID != rhs.ID || 
           fieldID != rhs.fieldID || whichComponent != rhs.whichComponent;
  }
};

class less_DirichletBCRecord {
 public:
  less_DirichletBCRecord(){}
  ~less_DirichletBCRecord(){}

  bool operator()(const DirichletBCRecord& lhs,
                  const DirichletBCRecord& rhs)
  {
    if (lhs.IDType < rhs.IDType) return true;
    if (lhs.IDType > rhs.IDType) return false;

    if (lhs.ID < rhs.ID) return true;
    if (lhs.ID > rhs.ID) return false;

    if (lhs.fieldID < rhs.fieldID) return true;
    if (lhs.fieldID > rhs.fieldID) return false;

    if (lhs.whichComponent < rhs.whichComponent) return true;
    if (lhs.whichComponent > rhs.whichComponent) return false;

    return false;
  }
};

}//namespace fei

#endif

