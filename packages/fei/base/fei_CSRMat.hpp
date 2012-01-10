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

#ifndef _fei_CSRMat_hpp_
#define _fei_CSRMat_hpp_


#include "fei_macros.hpp"
#include "fei_FillableMat.hpp"
#include "fei_SparseRowGraph.hpp"
#include "fei_CSVec.hpp"

namespace fei {

/** Compressed Sparse Row Matrix object.
*/
class CSRMat {
 public:
  CSRMat();
  CSRMat(const FillableMat& fmat);
  virtual ~CSRMat();

  SparseRowGraph& getGraph() {return srg_;}
  const SparseRowGraph& getGraph() const {return srg_;}
 
  std::vector<double>& getPackedCoefs() {return packedcoefs_;}
  const std::vector<double>& getPackedCoefs() const {return packedcoefs_;}

  unsigned getNumRows() const {return srg_.rowNumbers.size();}

  CSRMat& operator=(const FillableMat& src);

  CSRMat& operator+=(const CSRMat& src);

  bool operator==(const CSRMat& rhs) const;

  bool operator!=(const CSRMat& rhs) const;

 private:
  SparseRowGraph srg_;
  std::vector<double> packedcoefs_;
};//class CSRMat

/** form y = A*x */
void multiply_CSRMat_CSVec(const CSRMat& A, const CSVec& x, CSVec& y);

/** form y = A^T*x */
void multiply_trans_CSRMat_CSVec(const CSRMat& A, const CSVec& x, CSVec& y);

/** form C = A*B */
void multiply_CSRMat_CSRMat(const CSRMat& A, const CSRMat& B, CSRMat& C,
                            bool storeResultZeros=false);

/** form C = A^T*B */
void multiply_trans_CSRMat_CSRMat(const CSRMat& A, const CSRMat& B, CSRMat& C,
                                  bool storeResultZeros=false);

void add_CSRMat_to_FillableMat(const CSRMat& csrm, FillableMat& fm);

}//namespace fei

#endif

