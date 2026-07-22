// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "tGraphLaplacian_tpetra.hpp"

#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Map.hpp"

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

static ST exact[5][5] = {{1.0, -1.0, 0.0, 0.0, 0.0},
                         {-1.0, 1.12038585308577, 0.0, -0.120385853085769, 0.0},
                         {0.0, 0.0, 0.141421356237310, 0.0, -0.141421356237310},
                         {0.0, -0.120385853085769, 0.0, 0.120385853085769, 0.0},
                         {0.0, -0.110431526074847, -0.14142135623731, 0.0, 0.251852882312156}};

inline void simple_point(int i, std::vector<ST>& vec, ST x, ST y, ST z) {
  vec[3 * i + 0] = x;
  vec[3 * i + 1] = y;
  vec[3 * i + 2] = z;
}

void tGraphLaplacian_tpetra::initializeTest() { tolerance_ = 1e-12; }

Teuchos::RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > stencil(const Teuchos::Comm<int>& comm) {
  Teuchos::RCP<Tpetra::Map<LO, GO, NT> > map =
      rcp(new Tpetra::Map<LO, GO, NT>(5, 0, rcpFromRef(comm)));
  Teuchos::RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > mat =
      rcp(new Tpetra::CrsMatrix<ST, LO, GO, NT>(map, 6));
  // KDD 8/2019:  The value 6 above is needed because this test inserts
  // three indices twice into row 4 of the matrix.
  // Tpetra is robust enough to correctly handle the double insertion.
  // But if this double insertion is, in fact, an error, the 6 can probably
  // be changed to 3.  Note that this test doesn't insert any nonzeros into
  // row 0 of the matrix.

  // arrays for indicies and values
  GO indicies[3];
  ST values[3] = {1, 1, 1};

  // build linear system
  indicies[0] = 0;
  indicies[1] = 1;
  indicies[2] = 3;
  mat->insertGlobalValues(1, Teuchos::ArrayView<GO>(indicies, 3),
                          Teuchos::ArrayView<ST>(values, 3));

  indicies[0] = 4;
  mat->insertGlobalValues(2, Teuchos::ArrayView<GO>(indicies, 1),
                          Teuchos::ArrayView<ST>(values, 1));

  indicies[0] = 1;
  indicies[1] = 3;
  mat->insertGlobalValues(3, Teuchos::ArrayView<GO>(indicies, 2),
                          Teuchos::ArrayView<ST>(values, 2));

  indicies[0] = 1;
  indicies[1] = 2;
  indicies[2] = 4;
  mat->insertGlobalValues(4, Teuchos::ArrayView<GO>(indicies, 3),
                          Teuchos::ArrayView<ST>(values, 3));

  indicies[0] = 0;
  indicies[1] = 1;
  indicies[2] = 3;
  values[2]   = 0.0;
  mat->insertGlobalValues(4, Teuchos::ArrayView<GO>(indicies, 3),
                          Teuchos::ArrayView<ST>(values, 3));

  mat->fillComplete();

  return mat;
}

static void coords(std::vector<ST>& vec) {
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

static void coords(std::vector<ST>& x, std::vector<ST>& y, std::vector<ST>& z) {
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

int tGraphLaplacian_tpetra::runTest(int verbosity, std::ostream& stdstrm, std::ostream& failstrm,
                                    int& totalrun) {
  bool allTests = true;
  bool status;
  int failcount = 0;

  failstrm << "tGraphLaplacian_tpetra";

  status = test_single_array(verbosity, failstrm);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"single_array\" ... PASSED",
                       "   \"single_array\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = test_multi_array(verbosity, failstrm);
  Teko_TEST_MSG_tpetra(stdstrm, 1, "   \"multi_array\" ... PASSED",
                       "   \"multi_array\" ... FAILED");
  allTests &= status;
  failcount += status ? 0 : 1;
  totalrun++;

  status = allTests;
  if (verbosity >= 10) {
    Teko_TEST_MSG_tpetra(failstrm, 0, "tGraphLaplacian_tpetra...PASSED",
                         "tGraphLaplacian_tpetra...FAILED");
  } else {  // Normal Operatoring Procedures (NOP)
    Teko_TEST_MSG_tpetra(failstrm, 0, "...PASSED", "tGraphLaplacian_tpetra...FAILED");
  }

  return failcount;
}

bool tGraphLaplacian_tpetra::compareMatrix(const Tpetra::CrsMatrix<ST, LO, GO, NT>& gl,
                                           const std::string& name, int verbosity,
                                           std::ostream& os) const {
  bool status    = false;
  bool allPassed = true;

  size_t count;
  auto indices = typename Tpetra::CrsMatrix<ST, LO, GO, NT>::nonconst_global_inds_host_view_type(
      Kokkos::ViewAllocateWithoutInitializing("rowIndices"), 5);
  auto values = typename Tpetra::CrsMatrix<ST, LO, GO, NT>::nonconst_values_host_view_type(
      Kokkos::ViewAllocateWithoutInitializing("rowIndices"), 5);
  for (GO i = 0; i < 5; i++) {
    gl.getGlobalRowCopy(i, indices, values, count);

    for (size_t j = 0; j < count; j++) {
      GO col  = indices(j);
      ST diff = 0.0;
      if (exact[i][col] == 0.0)
        diff = std::fabs((exact[i][col] - values(j)));
      else
        diff = std::fabs((exact[i][col] - values(j)) / exact[i][col]);
      TEST_ASSERT(diff <= tolerance_, "\n   tGraphLaplacian_tpetra::"
                                          << name << ": " << toString(status) << "\n"
                                          << "      (row,col) = ( " << i << ", " << col << " )\n"
                                          << "      exact = " << exact[i][col] << "\n"
                                          << "      gl = " << values(j) << "\n"
                                          << "      rel err = " << diff << "<= " << tolerance_
                                          << "\n");
    }
  }

  return allPassed;
}

bool tGraphLaplacian_tpetra::test_single_array(int verbosity, std::ostream& os) {
  bool status    = false;
  bool allPassed = true;
  // double diff;
  std::vector<ST> points;

  // build coordinates and the stencil
  RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > sten = stencil(*GetComm_tpetra());
  coords(points);

  // build the graph laplacian
  RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > gl = buildGraphLaplacian(3, &points[0], *sten);

  TEST_ASSERT(compareMatrix(*gl, "test_single_array", verbosity, os),
              "\n   tGraphLaplacian_tpetra::test_single_array: "
                  << toString(status) << "\n"
                  << "      checked single array laplacian matrix");

  return allPassed;
}

bool tGraphLaplacian_tpetra::test_multi_array(int verbosity, std::ostream& os) {
  bool status    = false;
  bool allPassed = true;
  // double diff;
  std::vector<ST> x, y, z;
  std::vector<ST> points;

  // build coordinates and the stencil
  RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > sten = stencil(*GetComm_tpetra());
  coords(x, y, z);
  coords(points);

  // build the graph laplacian
  RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > gl = buildGraphLaplacian(&x[0], &y[0], &z[0], 1, *sten);

  TEST_ASSERT(compareMatrix(*gl, "test_multi_array", verbosity, os),
              "\n   tGraphLaplacian_tpetra::test_multi_array: "
                  << toString(status) << "\n"
                  << "      checked multi array laplacian matrix, unit stride");

  // build the graph laplacian
  RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > gl2 =
      buildGraphLaplacian(&points[0], &points[1], &points[2], 3, *sten);

  TEST_ASSERT(compareMatrix(*gl2, "test_multi_array", verbosity, os),
              "\n   tGraphLaplacian_tpetra::test_multi_array: "
                  << toString(status) << "\n"
                  << "      checked multi array laplacian matrix, 3 stride");

  return allPassed;
}

}  // namespace Test
}  // end namespace Teko
