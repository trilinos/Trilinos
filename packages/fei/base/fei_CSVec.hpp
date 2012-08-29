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

#ifndef _fei_CSVec_hpp_
#define _fei_CSVec_hpp_


#include <fei_macros.hpp>
#include <vector>
#include <algorithm>

namespace fei {

/** 'Compressed Sparse Vector' stored as two std::vectors: a vector of ints for the
    indices and a vector of doubles for the coefficients.

   Non-member functions add_entry and put_entry maintain sortedness of the
   vector when inserting new entries.
*/
class CSVec {
 public:
  CSVec(unsigned sz=0);
  virtual ~CSVec();

  CSVec& operator=(const CSVec& invec);

  std::vector<int>& indices() {return indices_;}
  const std::vector<int>& indices() const {return indices_;}
  std::vector<double>& coefs() {return coefs_;}
  const std::vector<double>& coefs() const {return coefs_;}

  size_t size() const {return indices_.size();}

  void clear() { indices_.clear(); coefs_.clear(); }

  bool operator==(const CSVec& rhs) const {
    return indices_==rhs.indices_ && coefs_==rhs.coefs_;
  }

  bool operator!=(const CSVec& rhs) const {
    return indices_!=rhs.indices_ || coefs_!=rhs.coefs_;
  }

 private:
  std::vector<int> indices_;
  std::vector<double> coefs_;
};//class CSVec

inline
void add_entry(CSVec& vec, int eqn, double coef)
{
  std::vector<int>& v_ind = vec.indices();
  std::vector<double>& v_coef = vec.coefs();

  std::vector<int>::iterator
    iter = std::lower_bound(v_ind.begin(), v_ind.end(), eqn);

  size_t offset = iter - v_ind.begin();

  if (iter == v_ind.end() || *iter != eqn) {
    v_ind.insert(iter, eqn);
    v_coef.insert(v_coef.begin()+offset, coef);
  }
  else {
    v_coef[offset] += coef;
  }
}


void add_entries(CSVec& vec, int num, const int* eqns, const double* coefs);

void put_entry(CSVec& vec, int eqn, double coef);

double get_entry(const CSVec& vec, int eqn);

void remove_entry(CSVec& vec, int eqn);

void set_values(CSVec& vec, double scalar);

/** form v = v + u
*/
void add_CSVec_CSVec(const CSVec& u, CSVec& v);

}//namespace fei

#endif

