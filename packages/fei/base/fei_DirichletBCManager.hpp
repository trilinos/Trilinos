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

#ifndef _fei_DirichletBCManager_hpp_
#define _fei_DirichletBCManager_hpp_


#include <fei_DirichletBCRecord.hpp>
#include <SNL_FEI_Structure.hpp>
#include <fei_VectorSpace.hpp>

#include <fei_Pool_alloc.hpp>
#include <map>

class NodeDatabase;
class EqnBuffer;

namespace fei {
class Matrix;

class DirichletBCManager {
 public:
  DirichletBCManager(SNL_FEI_Structure* structure)
   : structure_(structure), vecSpace_() {}

  DirichletBCManager(fei::SharedPtr<fei::VectorSpace> vecspace)
   : structure_(NULL), vecSpace_(vecspace) {}

  ~DirichletBCManager(){}

  void addBCRecords(int numBCs,
                    int IDType,
                    int fieldID,
                    int offsetIntoField,
                    const int* IDs,
                    const double* prescribedValues);

  void addBCRecords(int numBCs,
                    int IDType,
                    int fieldID,
                    const int* IDs,
                    const int* offsetsIntoField,
                    const double* prescribedValues);

  int finalizeBCEqns(fei::Matrix& matrix,
                     bool throw_if_bc_slave_conflict=false);

  int finalizeBCEqns(EqnBuffer& bcEqns);

  size_t getNumBCRecords() const;

  void clearAllBCs();

 private:
  int getEqnNumber(int IDType, int ID, int fieldID, int offsetIntoField);

  SNL_FEI_Structure* structure_;
  fei::SharedPtr<fei::VectorSpace> vecSpace_;

  typedef std::map<int,double,std::less<int>,
            fei_Pool_alloc<std::pair<const int, double> > > bc_map;
  bc_map bcs_;
};//class DirichletBCManager
}//namespace fei
#endif

