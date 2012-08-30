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



#include "fei_CSRMat.hpp"
#include <fei_impl_utils.hpp>
#include "fei_ArrayUtils.hpp"
#include <limits>
#include <cmath>

namespace fei {

CSRMat::CSRMat()
 : srg_(),
   packedcoefs_()
{
}

CSRMat::CSRMat(const FillableMat& fmat)
 : srg_(),
   packedcoefs_()
{
  *this = fmat;
}

CSRMat::~CSRMat()
{
}

CSRMat&
CSRMat::operator=(const fei::FillableMat& src)
{
  FillableMat::const_iterator iter = src.begin(), iter_end = src.end();

  unsigned nrows = src.getNumRows();

  srg_.rowNumbers.resize(nrows);
  srg_.rowOffsets.resize(nrows+1);

  unsigned nnz = 0;
  unsigned i = 0;
  for(; iter != iter_end; ++iter, ++i) {
    srg_.rowNumbers[i] = iter->first;
    srg_.rowOffsets[i] = nnz;
    nnz += iter->second->size();
  }

  srg_.rowOffsets[nrows] = nnz;

  srg_.packedColumnIndices.resize(nnz);
  packedcoefs_.resize(nnz);

  int* colind_ptr = (srg_.packedColumnIndices.size()
    ? &(srg_.packedColumnIndices[0]) : 0);
  double* coef_ptr = (packedcoefs_.size()
    ? &(packedcoefs_[0]) : 0);

  iter = src.begin();

  unsigned offset = 0;
  for(; iter != iter_end; ++iter) {
    const CSVec* v = iter->second;
    const std::vector<int>& v_ind = v->indices();
    const std::vector<double>& v_coef = v->coefs();
    for(size_t i=0; i<v_ind.size(); ++i) {
      colind_ptr[offset] = v_ind[i];
      coef_ptr[offset++] = v_coef[i];
    }
  }

  return *this;
}

CSRMat&
CSRMat::operator+=(const CSRMat& src)
{
  FillableMat tmp;
  add_CSRMat_to_FillableMat(*this, tmp);
  add_CSRMat_to_FillableMat(src, tmp);
  *this = tmp;
  return *this;
}

bool
CSRMat::operator==(const CSRMat& rhs) const
{
  if (getGraph() != rhs.getGraph()) return false;
  return getPackedCoefs() == rhs.getPackedCoefs();
}

bool
CSRMat::operator!=(const CSRMat& rhs) const
{
  return !(*this == rhs);
}

void multiply_CSRMat_CSVec(const CSRMat& A, const CSVec& x, CSVec& y)
{
  //This function is unit-tested in fei/utest_cases/fei_unit_CSRMat_CSVec.cpp

  const std::vector<int>& rows = A.getGraph().rowNumbers;
  const int* rowoffs = &(A.getGraph().rowOffsets[0]);
  const std::vector<int>& colinds = A.getGraph().packedColumnIndices;
  const double* Acoef = &(A.getPackedCoefs()[0]);

  const std::vector<int>& xind = x.indices();
  const std::vector<double>& xcoef = x.coefs();

  const double* xcoef_ptr = &xcoef[0];
  const int* xind_ptr = &xind[0];
  int xlen = xcoef.size();

  std::vector<int>& yind = y.indices();
  std::vector<double>& ycoef = y.coefs();

  unsigned nrows = A.getNumRows();

  yind.resize(nrows);
  ycoef.resize(nrows);

  int* yind_ptr = &yind[0];
  double* ycoef_ptr = &ycoef[0];

  int jbeg = *rowoffs++;
  for(unsigned i=0; i<nrows; ++i) {
    int jend = *rowoffs++;

    double sum = 0.0;
    while(jbeg<jend) {
      int xoff = fei::binarySearch(colinds[jbeg], xind_ptr, xlen);

      if (xoff > -1) {
        sum += Acoef[jbeg]*xcoef_ptr[xoff];
      }
      ++jbeg;
    }

    yind_ptr[i] = rows[i];
    ycoef_ptr[i] = sum;
  }
}

void multiply_trans_CSRMat_CSVec(const CSRMat& A, const CSVec& x, CSVec& y)
{
  const std::vector<int>& rows = A.getGraph().rowNumbers;
  const int* rowoffs = &(A.getGraph().rowOffsets[0]);
  const int* colinds = &(A.getGraph().packedColumnIndices[0]);
  const double* Acoef = &(A.getPackedCoefs()[0]);

  const std::vector<int>& xind = x.indices();
  const std::vector<double>& xcoef = x.coefs();

  const double* xcoef_ptr = &xcoef[0];

  unsigned nrows = A.getNumRows();

  std::vector<int> offsets;
  fei::impl_utils::find_offsets(rows, xind, offsets);
  const int* offsetsptr = &offsets[0];

  fei::CSVec fy;

  int jbeg = *rowoffs++;
  for(unsigned i=0; i<nrows; ++i) {
    int jend = *rowoffs++;

    int xoff = offsetsptr[i];
    if (xoff < 0) {
      jbeg = jend;
      continue;
    }

    double xcoeff = xcoef_ptr[xoff];

    while(jbeg<jend) {
      add_entry(fy, colinds[jbeg],Acoef[jbeg]*xcoeff);
      ++jbeg;
    }
  }

  y = fy;
}

void multiply_CSRMat_CSRMat(const CSRMat& A, const CSRMat& B, CSRMat& C,
                            bool storeResultZeros)
{
  //This function is unit-tested in fei/utest_cases/fei_unit_CSRMat_CSVec.cpp

  fei::FillableMat fc;

  const std::vector<int>& Arows = A.getGraph().rowNumbers;
  const std::vector<int>& Brows = B.getGraph().rowNumbers;
  if (Arows.size() < 1 || Brows.size() < 1) {
    C = fc;
    return;
  }
  const int* Arowoffs = &(A.getGraph().rowOffsets[0]);
  const int* Acols = &(A.getGraph().packedColumnIndices[0]);
  const double* Acoefs = &(A.getPackedCoefs()[0]);

  const int* Browoffs = &(B.getGraph().rowOffsets[0]);
  const std::vector<int>& Bcols = B.getGraph().packedColumnIndices;
  const double* Bcoefs = &(B.getPackedCoefs()[0]);

  static double fei_min = std::numeric_limits<double>::min();

  int jbeg = *Arowoffs++;
  for(size_t i=0; i<Arows.size(); ++i) {
    int row = Arows[i];
    int jend = *Arowoffs++;

    fei::CSVec* fc_row = NULL;
    if (storeResultZeros) {
      fc_row = fc.create_or_getRow(row);
    }
    else {
      fc_row = fc.hasRow(row) ? fc.create_or_getRow(row) : NULL;
    }

    while(jbeg<jend) {
      ++jbeg;
      int Acol = *Acols++;
      double Acoef = *Acoefs++;

      int Brow_offset = fei::binarySearch(Acol, &Brows[0], Brows.size());

      if (Brow_offset < 0) {
        continue;
      }

      if (!storeResultZeros) {
        if (std::abs(Acoef) < fei_min) {
          continue;
        }
      }

      const int* Brow_cols = &(Bcols[Browoffs[Brow_offset]]);
      const double* Brow_coefs = &(Bcoefs[Browoffs[Brow_offset]]);
      int Brow_len = Browoffs[Brow_offset+1]-Browoffs[Brow_offset];

      for(int k=0; k<Brow_len; ++k) {
        double resultCoef = Acoef*Brow_coefs[k];
        int resultCol = Brow_cols[k];

        if (!storeResultZeros) {
          if (std::abs(resultCoef) < fei_min) {
            continue;
          }
        }

        if (fc_row == NULL) {
          fc_row = fc.create_or_getRow(row);
        }

        add_entry(*fc_row, resultCol, resultCoef);
      }
    }
  }

  C = fc;
}

void multiply_trans_CSRMat_CSRMat(const CSRMat& A, const CSRMat& B, CSRMat& C,
                                  bool storeResultZeros)
{
  //This function is unit-tested in fei/utest_cases/fei_unit_CSRMat_CSVec.cpp

  fei::FillableMat fc;

  const std::vector<int>& Arows = A.getGraph().rowNumbers;
  const std::vector<int>& Brows = B.getGraph().rowNumbers;
  if (Arows.size() < 1 || Brows.size() < 1) {
    C = fc;
    return;
  }

  const size_t numArows = Arows.size();
  const int* Arowoffs = &(A.getGraph().rowOffsets[0]);
  const int* Acols = &(A.getGraph().packedColumnIndices[0]);
  const double* Acoefs = &(A.getPackedCoefs()[0]);

  const int* Browoffs = &(B.getGraph().rowOffsets[0]);
  const std::vector<int>& Bcols = B.getGraph().packedColumnIndices;
  const double* Bcoefs = &(B.getPackedCoefs()[0]);

  std::vector<double> row_coefs;

  static double fei_min = std::numeric_limits<double>::min();

  std::vector<int> offsets;
  fei::impl_utils::find_offsets(Arows, Brows, offsets);

  int jbeg = *Arowoffs++;
  for(size_t i=0; i<numArows; ++i) {
    int jend = *Arowoffs++;

    int Brow_offset = offsets[i];
    if (Brow_offset < 0) {
      jbeg = jend;
      continue;
    }

    const int* Brow_cols = &(Bcols[Browoffs[Brow_offset]]);
    const double* Brow_coefs = &(Bcoefs[Browoffs[Brow_offset]]);
    int Brow_len = Browoffs[Brow_offset+1]-Browoffs[Brow_offset];

    if ((int)row_coefs.size() < Brow_len) row_coefs.resize(Brow_len*2);
    double* row_coefs_ptr = &row_coefs[0];

    while(jbeg<jend) {
      int Acol = Acols[jbeg];
      double Acoef = Acoefs[jbeg++];

      if (std::abs(Acoef) < fei_min && !storeResultZeros) {
        continue;
      }

      for(int k=0; k<Brow_len; ++k) {
        row_coefs_ptr[k] = Acoef*Brow_coefs[k];
      }

      fc.sumInRow(Acol, Brow_cols, row_coefs_ptr, Brow_len);
    }
  }

  C = fc;
}

void add_CSRMat_to_FillableMat(const CSRMat& csrm, FillableMat& fm)
{
  const std::vector<int>& rows = csrm.getGraph().rowNumbers;
  const int* rowoffs = &(csrm.getGraph().rowOffsets[0]);
  const std::vector<int>& cols = csrm.getGraph().packedColumnIndices;
  const double* coefs = &(csrm.getPackedCoefs()[0]);

  for(size_t i=0; i<rows.size(); ++i) {
    int row = rows[i];

    for(int j=rowoffs[i]; j<rowoffs[i+1]; ++j) {
      fm.sumInCoef(row, cols[j], coefs[j]);
    }
  }
}

}//namespace fei

