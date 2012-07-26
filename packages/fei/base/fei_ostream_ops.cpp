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


#include <fei_macros.hpp>
#include <fei_ostream_ops.hpp>

#include <fei_Vector.hpp>
#include <fei_Matrix.hpp>
#include <fei_FillableMat.hpp>
#include <fei_FillableVec.hpp>
#include <fei_CSRMat.hpp>
#include <fei_CSVec.hpp>


FEI_OSTREAM& operator<<(FEI_OSTREAM& os, fei::Vector& vec)
{
  vec.writeToStream(os);
  return(os);
}

FEI_OSTREAM& operator<<(FEI_OSTREAM& os, fei::Matrix& mat)
{
  mat.writeToStream(os);
  return(os);
}

FEI_OSTREAM& operator<<(FEI_OSTREAM& os, fei::FillableVec& vec)
{
  fei::FillableVec::iterator
    iter = vec.begin(), iter_end = vec.end();

  os << "   numEntries: " << vec.size() << FEI_ENDL;

  for(; iter!=iter_end; ++iter) {
    os << "     " << iter->first<< ": "<<iter->second << FEI_ENDL;
  }

  return(os);
}

FEI_OSTREAM& operator<<(FEI_OSTREAM& os, fei::FillableMat& mat)
{
  os << "num rows: " << mat.getNumRows() << FEI_ENDL;
  fei::FillableMat::iterator
    iter = mat.begin(), iter_end = mat.end();

  for(; iter!=iter_end; ++iter) {
    int row = iter->first;
    const fei::CSVec* v = iter->second;
    const std::vector<int>& v_ind = v->indices();
    const std::vector<double>& v_coef = v->coefs();
    os << row << ": ";
    for(size_t i=0; i<v_ind.size(); ++i) {
      os << "("<<v_ind[i]<<","<<v_coef[i]<<") ";
    }
    os << FEI_ENDL;
  }

  return(os);
}

FEI_OSTREAM& operator<<(FEI_OSTREAM& os, fei::CSVec& vec)
{
  size_t len = vec.size();

  os << "   numEntries: " << len << FEI_ENDL;

  for(size_t i=0; i<len; ++i) {
    os << "     " << vec.indices()[i]<< ": "<<vec.coefs()[i] << FEI_ENDL;
  }

  return(os);
}

FEI_OSTREAM& operator<<(FEI_OSTREAM& os, fei::CSRMat& mat)
{
  os << "num rows: " << mat.getNumRows() << FEI_ENDL;

  const std::vector<int>& rows = mat.getGraph().rowNumbers;
  const int* rowoffs = &(mat.getGraph().rowOffsets[0]);
  const std::vector<int>& cols = mat.getGraph().packedColumnIndices;
  const double* coefs = &(mat.getPackedCoefs()[0]);

  for(size_t i=0; i<rows.size(); ++i) {
    int row = rows[i];

    os << row << ": ";
    for(int j=rowoffs[i]; j<rowoffs[i+1]; ++j) {
      os << "("<<cols[j]<<","<<coefs[j]<<") ";
    }
    os << FEI_ENDL;
  }

  return(os);
}

