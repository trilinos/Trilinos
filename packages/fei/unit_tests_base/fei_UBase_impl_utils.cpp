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


#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include <fei_impl_utils.hpp>
#include <fei_FillableMat.hpp>
#include <fei_CSVec.hpp>
#include <fei_CommUtils.hpp>

#include <fei_iostream.hpp>

#include <cmath>
#include <limits>

namespace {

bool verify_offsets(const std::vector<int>& src, const std::vector<int>& tgt,
                    const std::vector<int>& offsets)
{
  if (offsets.size() != src.size()) return false;

  for(size_t i=0; i<src.size(); ++i) {
    if (offsets[i] < 0) continue;

    if (src[i] != tgt[offsets[i]]) return false;
  }

  return true;
}

TEUCHOS_UNIT_TEST(impl_utils, find_offsets)
{
  std::vector<int> s1(5), t1(5), offs1;

  for(size_t i=0; i<5; ++i) {
    s1[i] = i;
    t1[i] = i;
  }

  fei::impl_utils::find_offsets(s1, t1, offs1);

  if (!verify_offsets(s1, t1, offs1)) {
    throw std::runtime_error("failed test 1.");
  }

  for(size_t i=0; i<5; ++i) {
    s1[i] = i;
    t1[i] = i + 2;
  }

  fei::impl_utils::find_offsets(s1, t1, offs1);

  if (!verify_offsets(s1, t1, offs1)) {
    throw std::runtime_error("failed test 2.");
  }

  for(size_t i=0; i<5; ++i) {
    s1[i] = i;
    t1[i] = i - 2;
  }

  fei::impl_utils::find_offsets(s1, t1, offs1);

  if (!verify_offsets(s1, t1, offs1)) {
    throw std::runtime_error("failed test 3.");
  }

  t1 = s1;
  t1.resize(t1.size()*2);
  int val = s1[s1.size()-1];
  for(size_t i=s1.size(); i<t1.size(); ++i) t1[i] = val+i;

  fei::impl_utils::find_offsets(s1, t1, offs1);

  if (!verify_offsets(s1, t1, offs1)) {
    throw std::runtime_error("failed test 4.");
  }
}

TEUCHOS_UNIT_TEST(impl_utils, pack_unpack_FillableMat)
{
  fei::FillableMat fm, fm2;

  fm.putCoef(0, 0, 0.0);
  fm.putCoef(0, 1, 0.1);
  fm.putCoef(1, 1, 1.1);
  fm.putCoef(1, 2, 1.2);
  fm.putCoef(2, 1, 2.1);
  fm.putCoef(2, 2, 2.2);

  std::vector<char> data(fei::impl_utils::num_bytes_FillableMat(fm));

  fei::impl_utils::pack_FillableMat(fm, &data[0]); 

  fei::impl_utils::unpack_FillableMat(&data[0], &data[0]+data.size(), fm2);

  if (fm.getNumRows() != fm2.getNumRows()) {
    throw std::runtime_error("pack/unpack FillableMat, wrong number of rows");
  }

  if (fei::count_nnz(fm) != fei::count_nnz(fm2)) {
    throw std::runtime_error("pack/unpack FillableMat, wrong number of nonzeros");
  }

  bool result = fm == fm2;
  TEUCHOS_TEST_EQUALITY(result, true, out, success);
}

TEUCHOS_UNIT_TEST(impl_utils, separateBCEqns)
{
  fei::FillableMat fm;

  fm.putCoef(0, 0, 0.0);
  fm.putCoef(2, 2, 2.0);
  fm.putCoef(3, 3, 3.0);
  fm.putCoef(4, 4, 4.0);
  fm.putCoef(5, 5, 5.0);

  std::vector<int> bcEqns;
  std::vector<double> bcVals;

  fei::impl_utils::separate_BC_eqns(fm, bcEqns, bcVals);

  if (bcEqns.size() != 5 || bcEqns.size() != bcVals.size()) {
    throw std::runtime_error("separate_BC_eqns test 1 failed.");
  }

  const double eps = std::numeric_limits<double>::epsilon();

  if (bcEqns[1] != 2 || std::abs(bcVals[1] - 2.0) > eps) {
    throw std::runtime_error("separate_BC_eqns test 2 failed.");
  }

  fei::FillableMat fm2;

  fm2.putCoef(1, 1, 1.0);

  fei::impl_utils::separate_BC_eqns(fm2, bcEqns, bcVals);

  if (bcEqns.size() != 6 || bcEqns.size() != bcVals.size()) {
    throw std::runtime_error("separate_BC_eqns test 3 failed.");
  }

  if (bcEqns[2] != 2 || std::abs(bcVals[2] - 2.0) > eps) {
    throw std::runtime_error("separate_BC_eqns test 4 failed.");
  }
}

TEUCHOS_UNIT_TEST(impl_utils, create_col_to_row_map)
{
  fei::FillableMat fm;

  fm.putCoef(0, 0, 0.0);
  fm.putCoef(0, 1, 0.1);
  fm.putCoef(1, 1, 1.1);
  fm.putCoef(2, 2, 2.2);

  std::multimap<int,int> crmap;

  fei::impl_utils::create_col_to_row_map(fm, crmap);

  if (crmap.size() != 4) {
    FEI_COUT << "ERROR, crmap.size()=="<<crmap.size()<<", expected 4."<<FEI_ENDL;
    throw std::runtime_error("create_col_to_row_map failed 1");
  }

  //next make sure that the col-to-row-map indicates that col 1 appears in 2 rows.
  typedef std::multimap<int,int>::iterator MM_Iter;
  std::pair<MM_Iter,MM_Iter> mm = crmap.equal_range(1);
  int num = 0;
  for(; mm.first!=mm.second; ++mm.first) ++num;

  if (num != 2) {
    FEI_COUT << "ERROR, size of equal_range(1)=="<<num<<", expected 2."<<FEI_ENDL;
    throw std::runtime_error("create_col_to_row_map failed 2");
  }
}

TEUCHOS_UNIT_TEST(impl_utils, remove_couplings)
{
  fei::FillableMat fm;

  fm.putCoef(2,  0, 0.5);
  fm.putCoef(2, 10, 0.5);
  fm.putCoef(8,  2, 0.5);
  fm.putCoef(8, 10, 0.5);

  int levels = fei::impl_utils::remove_couplings(fm);

  if (levels != 1) {
    throw std::runtime_error("remove_couplings test failed 1.");
  }

  //after remove_couplings, the matrix-row for 8 should have
  //2 column-indices, and they should be 0 and 10. Also, the
  //coefficients should be 0.25 and 0.75.
  const fei::CSVec* matrow = fm.getRow(8);

  if (matrow->size() != 2) {
    throw std::runtime_error("matrow 8 has wrong length");
  }

  const std::vector<int>& indices = matrow->indices();
  const std::vector<double>& coefs = matrow->coefs();
  if (indices[0] != 0 || indices[1] != 10 ||
      std::abs(coefs[0] -0.25) > 1.e-49 || std::abs(coefs[1] -0.75) > 1.e-49) {
    throw std::runtime_error("matrow 8 has wrong contents after remove_couplings");
  }

  levels = fei::impl_utils::remove_couplings(fm);
  if (levels > 0) {
    throw std::runtime_error("remove_couplings test 2 failed");
  }
}

TEUCHOS_UNIT_TEST(impl_utils, global_union_mat)
{
  MPI_Comm comm = MPI_COMM_WORLD;

  int numProcs = fei::numProcs(comm);
  int localProc = fei::localProc(comm);

  int numlocalrows = 5;
  int rowlen = 5;

  fei::FillableMat globalmat0;
  fei::FillableMat localmat;
  int row=0;
  for(int p=0; p<numProcs; ++p) {
    for(int i=0; i<numlocalrows; ++i) {
      for(int j=0; j<rowlen; ++j) {
        globalmat0.putCoef(row, j, 1.0);
        if (p == localProc) {
          localmat.putCoef(row, j, 1.0);
        }
      }
      ++row;
    }
  }

  fei::FillableMat globalmat;

  fei::impl_utils::global_union(comm, localmat, globalmat);

  if (globalmat != globalmat0) {
    throw std::runtime_error("globalUnion test (mat) failed");
  }
}

}//namespace <anonymous>

