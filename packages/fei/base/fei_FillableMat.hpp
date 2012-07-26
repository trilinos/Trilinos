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


#ifndef _fei_FillableMat_hpp_
#define _fei_FillableMat_hpp_

#include <vector>
#include <fei_CSVec.hpp>
#include <fei_EqnBuffer.hpp>

namespace fei {

class FillableMat {
 public:
  FillableMat();
  FillableMat(EqnBuffer& eqnbuf);
  virtual ~FillableMat();

  FillableMat& operator=(const FillableMat& src);

  void setValues(double value);

  void createPosition(int row, int col);

  void sumInCoef(int row, int col, double coef);
  void putCoef(int row, int col, double coef);

  void sumInRow(int row, const int* cols, const double* coefs, unsigned len);
  void putRow(int row, const int* cols, const double* coefs, unsigned len);

  unsigned getNumRows() const;

  bool hasRow(int row) const;

  const CSVec* getRow(int row) const;
  CSVec* create_or_getRow(int row);

  typedef std::map<int, CSVec*, std::less<int>,
                fei_Pool_alloc<std::pair<const int,CSVec*> > > feipoolmat;

  typedef feipoolmat::iterator iterator;
  typedef feipoolmat::const_iterator const_iterator;

  iterator begin() {return matdata_.begin();}
  iterator end() {return matdata_.end();}

  const_iterator begin() const {return matdata_.begin();}
  const_iterator end() const {return matdata_.end();}

  void clear();

  bool operator==(const FillableMat& rhs) const;

  bool operator!=(const FillableMat& rhs) const;

 private:
  feipoolmat matdata_;
  fei_Pool_alloc<CSVec> vecpool_;
}; //class FillableMat

/** Print the contents of a FillableMat to the given std::ostream. */
void print(std::ostream& os, const FillableMat& mat);

/** Return the number of nonzeros in mat. */
int count_nnz(const FillableMat& mat);

/** Fill a std::vector with the row-numbers from the given matrix. */
void get_row_numbers(const FillableMat& mat, std::vector<int>& rows);

}//namespace fei

#endif

