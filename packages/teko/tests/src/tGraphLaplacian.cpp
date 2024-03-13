/*
// @HEADER
//
// ***********************************************************************
//
//      Teko: A package for block and physics based preconditioning
//                  Copyright 2010 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Eric C. Cyr (eccyr@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

*/

#include "tGraphLaplacian.hpp"

#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"

#include "Teko_Utilities.hpp"

using Teuchos::RCP;
using Teuchos::rcp;

namespace Teko {
namespace Test {

#define hard_point(i, m, n, o) \
  {                            \
    x[i] = m;                  \
    y[i] = n;                  \
    z[i] = o;                  \
  }

static double exact[5][5] = {{1.0, -1.0, 0.0, 0.0, 0.0},
                             {-1.0, 1.12038585308577, 0.0, -0.120385853085769, 0.0},
                             {0.0, 0.0, 0.141421356237310, 0.0, -0.141421356237310},
                             {0.0, -0.120385853085769, 0.0, 0.120385853085769, 0.0},
                             {0.0, -0.110431526074847, -0.14142135623731, 0.0, 0.251852882312156}};

inline void simple_point(int i, std::vector<double>& vec, double x, double y, double z) {
  vec[3 * i + 0] = x;
  vec[3 * i + 1] = y;
  vec[3 * i + 2] = z;
}

void tGraphLaplacian::initializeTest() { tolerance_ = 1e-12; }

Teuchos::RCP<Epetra_CrsMatrix> stencil(const Epetra_Comm& comm) {
  Teuchos::RCP<Epetra_Map> map       = rcp(new Epetra_Map(5, 0, comm));
  Teuchos::RCP<Epetra_CrsMatrix> mat = rcp(new Epetra_CrsMatrix(Copy, *map, 0));

  // arrays for indicies and values
  int indicies[3];
  double values[3] = {1, 1, 1};

  // build linear system
  indicies[0] = 0;
  indicies[1] = 1;
  indicies[2] = 3;
  mat->InsertGlobalValues(1, 3, values, indicies);

  indicies[0] = 4;
  mat->InsertGlobalValues(2, 1, values, indicies);

  indicies[0] = 1;
  indicies[1] = 3;
  mat->InsertGlobalValues(3, 2, values, indicies);

  indicies[0] = 1;
  indicies[1] = 2;
  indicies[2] = 4;
  mat->InsertGlobalValues(4, 3, values, indicies);

  indicies[0] = 0;
  indicies[1] = 1;
  indicies[2] = 3;
  values[2]   = 0.0;
  mat->InsertGlobalValues(0, 3, values, indicies);

  mat->FillComplete();

  return mat;
}

static void coords(std::vector<double>& vec) {
  int count = 5;
  int dim   = 3;

  vec.clear();
  vec.resize(dim * count);

  simple_point(0, vec, 0.0, 0.0, 0.0);
  simple_point(1, vec, 0.0, -1.0, 0.0);
  simple_point(2, vec, 1.0, 0.0, 2.0);
  simple_point(3, vec, -2.0, 3.0, 7.0);
  simple_point(4, vec, 0.0, 0.0, 9.0);
}

static void coords(std::vector<double>& x, std::vector<double>& y, std::vector<double>& z) {
  int count = 5;
  // int dim = 3;

  x.clear();
  x.resize(count);
  y.clear();
  y.resize(count);
  z.clear();
  z.resize(count);

  hard_point(0, 0.0, 0.0, 0.0);
  hard_point(1, 0.0, -1.0, 0.0);
  hard_point(2, 1.0, 0.0, 2.0);
  hard_point(3, -2.0, 3.0, 7.0);
  hard_point(4, 0.0, 0.0, 9.0);
}

int tGraphLaplacian::runTest(int verbosity, std::ostream& stdstrm, std::ostream& failstrm,
                             int& totalrun) {
  bool allTests = true;
  bool status;
  int failcount = 0;

  failstrm << "tGraphLaplacian";

  status = test_single_array(verbosity, failstrm);
  Teko_TEST_MSG(stdstrm, 1, "   \"single_array\" ... PASSED", "   \"single_array\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_multi_array(verbosity, failstrm);
  Teko_TEST_MSG(stdstrm, 1, "   \"multi_array\" ... PASSED", "   \"multi_array\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = allTests;
  if (verbosity >= 10) {
    Teko_TEST_MSG(failstrm, 0, "tGraphLaplacian...PASSED", "tGraphLaplacian...FAILED");
  } else {  // Normal Operatoring Procedures (NOP)
    Teko_TEST_MSG(failstrm, 0, "...PASSED", "tGraphLaplacian...FAILED");
  }

  return failcount;
}

bool tGraphLaplacian::compareMatrix(const Epetra_CrsMatrix& gl, const std::string& name,
                                    int verbosity, std::ostream& os) const {
  bool status    = false;
  bool allPassed = true;

  int count;
  int indicies[5];
  double values[5];
  for (int i = 0; i < 5; i++) {
    gl.ExtractGlobalRowCopy(i, 5, count, values, indicies);

    for (int j = 0; j < count; j++) {
      int col     = indicies[j];
      double diff = 0.0;
      if (exact[i][col] == 0.0)
        diff = std::fabs((exact[i][col] - values[j]));
      else
        diff = std::fabs((exact[i][col] - values[j]) / exact[i][col]);
      TEST_ASSERT(diff <= tolerance_,
                  "\n   tGraphLaplacian::" << name << ": " << toString(status) << "\n"
                                           << "      (row,col) = ( " << i << ", " << col << " )\n"
                                           << "      exact = " << exact[i][col] << "\n"
                                           << "      gl = " << values[j] << "\n"
                                           << "      rel err = " << diff << "<= " << tolerance_
                                           << "\n");
    }
  }

  return allPassed;
}

bool tGraphLaplacian::test_single_array(int verbosity, std::ostream& os) {
  bool status    = false;
  bool allPassed = true;
  // double diff;
  std::vector<double> points;

  // build coordinates and the stencil
  RCP<Epetra_CrsMatrix> sten = stencil(*GetComm());
  coords(points);

  // build the graph laplacian
  RCP<Epetra_CrsMatrix> gl = buildGraphLaplacian(3, &points[0], *sten);

  TEST_ASSERT(
      compareMatrix(*gl, "test_single_array", verbosity, os),
      "\n   tGraphLaplacian::test_single_array: " << toString(status) << "\n"
                                                  << "      checked single array laplacian matrix");

  return allPassed;
}

bool tGraphLaplacian::test_multi_array(int verbosity, std::ostream& os) {
  bool status    = false;
  bool allPassed = true;
  // double diff;
  std::vector<double> x, y, z;
  std::vector<double> points;

  // build coordinates and the stencil
  RCP<Epetra_CrsMatrix> sten = stencil(*GetComm());
  coords(x, y, z);
  coords(points);

  // build the graph laplacian
  RCP<Epetra_CrsMatrix> gl = buildGraphLaplacian(&x[0], &y[0], &z[0], 1, *sten);

  TEST_ASSERT(compareMatrix(*gl, "test_multi_array", verbosity, os),
              "\n   tGraphLaplacian::test_multi_array: "
                  << toString(status) << "\n"
                  << "      checked multi array laplacian matrix, unit stride");

  // build the graph laplacian
  RCP<Epetra_CrsMatrix> gl2 = buildGraphLaplacian(&points[0], &points[1], &points[2], 3, *sten);

  TEST_ASSERT(compareMatrix(*gl2, "test_multi_array", verbosity, os),
              "\n   tGraphLaplacian::test_multi_array: "
                  << toString(status) << "\n"
                  << "      checked multi array laplacian matrix, 3 stride");

  return allPassed;
}

}  // namespace Test
}  // end namespace Teko
