/*
//@HEADER
// ***********************************************************************
//
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

/// \file Ifpack2_UnitTestEquilibration.cpp
/// \brief Unit test for Ifpack2's equilibration

#include "Teuchos_UnitTestHarness.hpp"
#include "Tpetra_Details_EquilibrationInfo.hpp"
#include "Tpetra_computeRowAndColumnOneNorms.hpp"
#include "Tpetra_computeRowAndColumnOneNorms_def.hpp" // unfortunate, but needed for test below; USERS DON'T DO THIS!
#include "Tpetra_leftAndOrRightScaleCrsMatrix.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_RowMatrix.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Details_gathervPrint.hpp"
#include "Teuchos_CommHelpers.hpp"
#include <array>
#include <functional>
#include <type_traits>

namespace { // (anonymous)

template<class SC, const bool isComplex = Kokkos::ArithTraits<SC>::is_complex>
struct NaughtyValues {
  static SC infinity ();
  static SC quiet_NaN ();
};

template<class SC>
struct NaughtyValues<SC, false> {
  static SC infinity () {
    return std::numeric_limits<SC>::infinity ();
  }
  static SC quiet_NaN () {
    return std::numeric_limits<SC>::quiet_NaN ();
  }
};

// mfh 15 Jun 2018: std::numeric_limits' infinity() and quiet_NaN()
// methods don't seem to worki right for std::complex or
// Kokkos::complex.  The code below works around that.  This test
// passes whether just the real part is (Inf or NaN), or whether both
// real and imaginary parts are (Inf or NaN).
template<class SC>
struct NaughtyValues<SC, true> {
  using mag_type = typename Kokkos::ArithTraits<SC>::mag_type;

  static SC infinity () {
    return SC (std::numeric_limits<mag_type>::infinity (), 0.0);
  }
  static SC quiet_NaN () {
    return SC (std::numeric_limits<mag_type>::quiet_NaN (), 0.0);
  }
};

template<class SC, class LO, class GO, class NT>
Tpetra::Vector<SC, LO, GO, NT>
createVectorFromCopyOf1DView (const Teuchos::RCP<const Tpetra::Map<LO, GO, NT> >& map,
                              const Kokkos::View<SC*, typename NT::device_type>& inputView)
{
  using dual_view_type = typename Tpetra::Vector<SC, LO, GO, NT>::dual_view_type;
  using dev_memory_space = typename NT::device_type::memory_space;

  dual_view_type dv ("MV::DualView", inputView.extent (0), 1);
  dv.template modify<dev_memory_space> ();

  auto outputView_2d = dv.template view<dev_memory_space> ();
  auto outputView_1d = Kokkos::subview (outputView_2d, Kokkos::ALL (), 0);

  Kokkos::deep_copy (outputView_1d, inputView);
  return Tpetra::Vector<SC, LO, GO, NT> (map, dv);
}

template<class ValueType>
bool
near (const ValueType& x,
      const ValueType& y,
      const typename Kokkos::ArithTraits<ValueType>::mag_type& factor)
{
  using KAT = Kokkos::ArithTraits<ValueType>;

  if (KAT::isNan (x) && KAT::isNan (y)) {
    // Any comparison involving a NaN will fail, so we need a special case.
    return true;
  }
  else if (KAT::isInf (x) && KAT::isInf (y)) {
    // Inf - Inf -> NaN, so we need a special case.
    return true;
  }
  else {
    const auto eps = KAT::eps ();
    const auto absDiff = KAT::abs (x - y);
    return absDiff <= factor * eps;
  }
}

template<class SC, class LO, class GO, class NT>
Teuchos::RCP<Tpetra::CrsMatrix<SC, LO, GO, NT> >
deepCopyFillCompleteCrsMatrix (const Tpetra::CrsMatrix<SC, LO, GO, NT>& A)
{
  using Teuchos::RCP;
  using crs_matrix_type = Tpetra::CrsMatrix<SC, LO, GO, NT>;

  TEUCHOS_TEST_FOR_EXCEPTION
    (! A.isFillComplete (), std::invalid_argument,
     "deepCopyFillCompleteCrsMatrix: Input matrix A must be fillComplete.");
  RCP<crs_matrix_type> A_copy (new crs_matrix_type (A.getCrsGraph ()));
  auto A_copy_lcl = A_copy->getLocalMatrix ();
  auto A_lcl = A.getLocalMatrix ();
  Kokkos::deep_copy (A_copy_lcl.values, A_lcl.values);
  A_copy->fillComplete (A.getDomainMap (), A.getRangeMap ());
  return A_copy;
}

template<class SC, class LO, class GO, class NT>
void
testCrsMatrixEquality (bool& success,
                       Teuchos::FancyOStream& out,
                       const Tpetra::CrsMatrix<SC, LO, GO, NT>& A_expected,
                       const Tpetra::CrsMatrix<SC, LO, GO, NT>& A_actual)
{
  using std::endl;
  using mag_type = typename Kokkos::ArithTraits<SC>::mag_type;
  const mag_type toleranceFactor = 10.0; // factor of eps

  auto A_expected_lcl = A_expected.getLocalMatrix ();
  auto ptr_h = Kokkos::create_mirror_view (A_expected_lcl.graph.row_map);
  Kokkos::deep_copy (ptr_h, A_expected_lcl.graph.row_map);

  auto expected_val_h = Kokkos::create_mirror_view (A_expected_lcl.values);
  Kokkos::deep_copy (expected_val_h, A_expected_lcl.values);

  auto A_actual_lcl = A_actual.getLocalMatrix ();
  auto actual_val_h = Kokkos::create_mirror_view (A_actual_lcl.values);
  Kokkos::deep_copy (actual_val_h, A_actual_lcl.values);

  using size_type = typename decltype (A_actual_lcl.graph)::size_type;

  if (A_expected_lcl.numRows () != A_actual_lcl.numRows ()) {
    TEST_EQUALITY( A_expected_lcl.numRows (), A_actual_lcl.numRows () );
    return;
  }
  const LO lclNumRows = A_expected_lcl.numRows ();

  bool ok = true;
  for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
    for (size_type k = ptr_h[lclRow]; k < ptr_h[lclRow+1]; ++k) {
      if (! near (actual_val_h[k],
                  expected_val_h[k],
                  toleranceFactor)) {
        ok = false;
        break;
      }
    }
  }
  if (! ok) {
    for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
      out << "lclRow: " << lclRow << endl;
      Teuchos::OSTab tab2 (out);

      using size_type = typename decltype (A_actual_lcl.graph)::size_type;
      out << "Expected values: [";
      for (size_type k = ptr_h[lclRow]; k < ptr_h[lclRow+1]; ++k) {
        out << expected_val_h[k];
        if (k + size_type (1) < ptr_h[lclRow+1]) {
          out << ", ";
        }
      }
      out << "]" << endl
          << "Actual values: [";
      for (size_type k = ptr_h[lclRow]; k < ptr_h[lclRow+1]; ++k) {
        out << actual_val_h[k];
        if (k + size_type (1) < ptr_h[lclRow+1]) {
          out << ", ";
        }
      }
      out << "]" << endl;
    }
  }
}

template<class SC, class LO, class GO, class NT>
struct EquilibrationTest {
  using map_type = Tpetra::Map<LO, GO, NT>;
  using crs_graph_type = Tpetra::CrsGraph<LO, GO, NT>;
  using crs_matrix_type = Tpetra::CrsMatrix<SC, LO, GO, NT>;
  using val_type = typename Kokkos::ArithTraits<SC>::val_type;
  using mag_type = typename Kokkos::ArithTraits<val_type>::mag_type;

  Teuchos::RCP<crs_matrix_type> A;
  std::vector<mag_type> lclRowNorms;
  std::vector<val_type> lclRowDiagonalEntries;
  std::vector<mag_type> lclColNorms;
  std::vector<val_type> lclColDiagonalEntries;
  std::vector<mag_type> lclRowScaledColNorms;
  std::vector<mag_type> gblRowNorms;
  std::vector<mag_type> gblColNorms;
  std::vector<mag_type> gblRowScaledColNorms;
  Teuchos::RCP<crs_matrix_type> A_leftScaled;
  Teuchos::RCP<crs_matrix_type> A_rightScaled;
  bool lclFoundInf;
  bool lclFoundNan;
  bool lclFoundZeroDiag;
  bool lclFoundZeroRowNorm;
  bool gblFoundInf;
  bool gblFoundNan;
  bool gblFoundZeroDiag;
  bool gblFoundZeroRowNorm;
};

template<class SC, class LO, class GO, class NT>
void
testEquilibration (Teuchos::FancyOStream& out,
                   bool& success,
                   const EquilibrationTest<SC, LO, GO, NT>& test,
                   const bool assumeSymmetric)
{
  using Teuchos::outArg;
  using Teuchos::RCP;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  using std::endl;
  using map_type = Tpetra::Map<LO, GO, NT>;
  using row_matrix_type = Tpetra::RowMatrix<SC, LO, GO, NT>;
  using mag_type = typename Kokkos::ArithTraits<SC>::mag_type;

  const mag_type toleranceFactor = 10.0; // factor of eps
  int lclSuccess = success ? 1 : 0; // for reduceAll (see below)
  int gblSuccess = 0; // output argument for reduceAll (see below)
  auto comm = test.A->getMap ()->getComm ();

  out << "Test equilibration: assumeSymmetric="
      << (assumeSymmetric ? "true" : "false") << endl;
  Teuchos::OSTab tab1 (out);

  const LO lclNumRows =
    static_cast<LO> (test.A->getRowMap ()->getNodeNumElements ());
  RCP<const map_type> colMap = test.A->getColMap ();
  const LO lclNumCols = static_cast<LO> (colMap->getNodeNumElements ());

  // Test computeLocalRowAndColumnOneNorms (CrsMatrix)
  {
    out << "Test computeLocalRowAndColumnOneNorms (CrsMatrix)" << endl;
    Teuchos::OSTab tab2 (out);
    auto result0 = Tpetra::Details::computeLocalRowAndColumnOneNorms (* (test.A), assumeSymmetric);

    {
      out << "Test detection of local error conditions" << endl;
      Teuchos::OSTab tab3 (out);

      TEST_EQUALITY( test.lclFoundInf, result0.foundInf );
      TEST_EQUALITY( test.lclFoundNan, result0.foundNan );
      TEST_EQUALITY( test.lclFoundZeroDiag, result0.foundZeroDiag );
      TEST_EQUALITY( test.lclFoundZeroRowNorm, result0.foundZeroRowNorm );
    }

    {
      out << "Test local row norms" << endl;
      Teuchos::OSTab tab3 (out);

      auto rowNorms_h = Kokkos::create_mirror_view (result0.rowNorms);
      Kokkos::deep_copy (rowNorms_h, result0.rowNorms);

      std::ostringstream os;
      os << "Expected local rowNorms: [";
      for (LO k = 0; k < lclNumRows; ++k) {
        os << test.lclRowNorms[k];
        if (k + LO (1) < lclNumRows) {
          os << ", ";
        }
      }
      os << "]" << endl << "Actual local rowNorms: [";
      for (LO k = 0; k < lclNumRows; ++k) {
        os << rowNorms_h[k];
        if (k + LO (1) < lclNumRows) {
          os << ", ";
        }
      }
      os << "]" << endl;
      out << os.str ();

      bool ok = true;
      for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
        if (! near (rowNorms_h[lclRow],
                    test.lclRowNorms[lclRow],
                    toleranceFactor)) {
          ok = false;
          break;
        }
      }
      TEST_ASSERT( ok );
    } // test local row norms

    {
      out << "Test local row diagonal entries" << endl;
      Teuchos::OSTab tab3 (out);

      auto rowDiagonalEntries_h =
        Kokkos::create_mirror_view (result0.rowDiagonalEntries);
      Kokkos::deep_copy (rowDiagonalEntries_h, result0.rowDiagonalEntries);

      std::ostringstream os;
      os << "Expected local rowDiagonalEntries: [";
      for (LO k = 0; k < lclNumRows; ++k) {
        os << test.lclRowDiagonalEntries[k];
        if (k + LO (1) < lclNumRows) {
          os << ", ";
        }
      }
      os << "]" << endl << "Actual local rowDiagonalEntries: [";
      for (LO k = 0; k < lclNumRows; ++k) {
        os << rowDiagonalEntries_h[k];
        if (k + LO (1) < lclNumRows) {
          os << ", ";
        }
      }
      os << "]" << endl;
      out << os.str ();

      bool ok = true;
      for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
        if (! near (rowDiagonalEntries_h[lclRow],
                    test.lclRowDiagonalEntries[lclRow],
                    toleranceFactor)) {
          ok = false;
          break;
        }
      }
      TEST_ASSERT( ok );
    } // test local row diagonal entries

    if (! assumeSymmetric) {
      out << "Test local column norms" << endl;
      Teuchos::OSTab tab3 (out);

      auto colNorms_h = Kokkos::create_mirror_view (result0.colNorms);
      Kokkos::deep_copy (colNorms_h, result0.colNorms);

      std::ostringstream os;
      os << "Expected local colNorms: [";
      for (LO k = 0; k < lclNumCols; ++k) {
        os << test.lclColNorms[k];
        if (k + LO (1) < lclNumCols) {
          os << ", ";
        }
      }
      os << "]" << endl << "Actual local colNorms: [";
      for (LO k = 0; k < lclNumCols; ++k) {
        os << colNorms_h[k];
        if (k + LO (1) < lclNumCols) {
          os << ", ";
        }
      }
      os << "]" << endl;
      out << os.str ();

      bool ok = true;
      for (LO lclCol = 0; lclCol < lclNumCols; ++lclCol) {
        if (! near (colNorms_h[lclCol],
                    test.lclColNorms[lclCol],
                    toleranceFactor)) {
          ok = false;
          break;
        }
      }
      TEST_ASSERT( ok );
    } // test local column norms

    if (! assumeSymmetric) {
      out << "Test local column diagonal entries" << endl;
      Teuchos::OSTab tab3 (out);

      auto colDiagonalEntries_h =
        Kokkos::create_mirror_view (result0.colDiagonalEntries);
      Kokkos::deep_copy (colDiagonalEntries_h, result0.colDiagonalEntries);

      std::ostringstream os;
      os << "Expected local colDiagonalEntries: [";
      for (LO k = 0; k < lclNumCols; ++k) {
        os << test.lclColDiagonalEntries[k];
        if (k + LO (1) < lclNumCols) {
          os << ", ";
        }
      }
      os << "]" << endl << "Actual local colDiagonalEntries: [";
      for (LO k = 0; k < lclNumCols; ++k) {
        os << colDiagonalEntries_h[k];
        if (k + LO (1) < lclNumCols) {
          os << ", ";
        }
      }
      os << "]" << endl;
      out << os.str ();

      bool ok = true;
      for (LO lclCol = 0; lclCol < lclNumCols; ++lclCol) {
        if (! near (colDiagonalEntries_h[lclCol],
                    test.lclColDiagonalEntries[lclCol],
                    toleranceFactor)) {
          ok = false;
          break;
        }
      }
      TEST_ASSERT( ok );
    } // test local column diagonal entries

    out << "Test globalizing row one norms" << endl;
    Tpetra::Details::globalizeRowOneNorms (result0, * (test.A));

    {
      out << "Test detection of global error conditions" << endl;
      Teuchos::OSTab tab3 (out);

      TEST_EQUALITY( test.gblFoundInf, result0.foundInf );
      TEST_EQUALITY( test.gblFoundNan, result0.foundNan );
      TEST_EQUALITY( test.gblFoundZeroDiag, result0.foundZeroDiag );
      TEST_EQUALITY( test.gblFoundZeroRowNorm, result0.foundZeroRowNorm );
    }

    {
      out << "Test global row norms" << endl;
      Teuchos::OSTab tab3 (out);

      auto rowNorms_h = Kokkos::create_mirror_view (result0.rowNorms);
      Kokkos::deep_copy (rowNorms_h, result0.rowNorms);

      std::ostringstream os;
      os << "Expected global rowNorms: [";
      for (LO k = 0; k < lclNumRows; ++k) {
        os << test.gblRowNorms[k];
        if (k + LO (1) < lclNumRows) {
          os << ", ";
        }
      }
      os << "]" << endl << "Actual global rowNorms: [";
      for (LO k = 0; k < lclNumRows; ++k) {
        os << rowNorms_h[k];
        if (k + LO (1) < lclNumRows) {
          os << ", ";
        }
      }
      os << "]" << endl;
      out << os.str ();

      bool ok = true;
      for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
        if (! near (rowNorms_h[lclRow],
                    test.gblRowNorms[lclRow],
                    toleranceFactor)) {
          ok = false;
          break;
        }
      }
      TEST_ASSERT( ok );
    } // test global row norms

    if (! assumeSymmetric) {
      out << "Test computeLocalRowScaledColumnNorms" << endl;
      Tpetra::Details::computeLocalRowScaledColumnNorms (result0, * (test.A));

      auto rowScaledColNorms_h =
        Kokkos::create_mirror_view (result0.rowScaledColNorms);
      Kokkos::deep_copy (rowScaledColNorms_h, result0.rowScaledColNorms);

      std::ostringstream os;
      os << "Expected local rowScaledColNorms: [";
      for (LO k = 0; k < lclNumCols; ++k) {
        os << test.lclRowScaledColNorms[k];
        if (k + LO (1) < lclNumCols) {
          os << ", ";
        }
      }
      os << "]" << endl << "Actual local rowScaledColNorms: [";
      for (LO k = 0; k < lclNumCols; ++k) {
        os << rowScaledColNorms_h[k];
        if (k + LO (1) < lclNumCols) {
          os << ", ";
        }
      }
      os << "]" << endl;
      out << os.str ();

      bool ok = true;
      for (LO lclCol = 0; lclCol < lclNumCols; ++lclCol) {
        if (! near (rowScaledColNorms_h[lclCol],
                    test.lclRowScaledColNorms[lclCol],
                    toleranceFactor)) {
          ok = false;
          break;
        }
      }
      TEST_ASSERT( ok );
    } // test computeLocalRowScaledColumnNorms

    out << "Test globalizing column one norms" << endl;
    Tpetra::Details::globalizeColumnOneNorms (result0, * (test.A), assumeSymmetric);

    {
      out << "Test global column norms" << endl;
      Teuchos::OSTab tab3 (out);

      auto colNorms_h = Kokkos::create_mirror_view (result0.colNorms);
      Kokkos::deep_copy (colNorms_h, result0.colNorms);

      std::ostringstream os;
      os << "Expected global colNorms: [";
      for (LO k = 0; k < lclNumCols; ++k) {
        os << test.gblColNorms[k];
        if (k + LO (1) < lclNumCols) {
          os << ", ";
        }
      }
      os << "]" << endl << "Actual global colNorms: [";
      for (LO k = 0; k < lclNumCols; ++k) {
        os << colNorms_h[k];
        if (k + LO (1) < lclNumCols) {
          os << ", ";
        }
      }
      os << "]" << endl;
      out << os.str ();

      bool ok = true;
      for (LO lclCol = 0; lclCol < lclNumCols; ++lclCol) {
        if (! near (colNorms_h[lclCol],
                    test.gblColNorms[lclCol],
                    toleranceFactor)) {
          ok = false;
          break;
        }
      }
      TEST_ASSERT( ok );
    } // test global column norms

    if (! assumeSymmetric) {
      out << "Test global row-scaled column norms" << endl;
      Teuchos::OSTab tab3 (out);

      auto rowScaledColNorms_h =
        Kokkos::create_mirror_view (result0.rowScaledColNorms);
      Kokkos::deep_copy (rowScaledColNorms_h, result0.rowScaledColNorms);

      std::ostringstream os;
      os << "Expected global rowScaledColNorms: [";
      for (LO k = 0; k < lclNumCols; ++k) {
        os << test.gblRowScaledColNorms[k];
        if (k + LO (1) < lclNumCols) {
          os << ", ";
        }
      }
      os << "]" << endl << "Actual global rowScaledColNorms: [";
      for (LO k = 0; k < lclNumCols; ++k) {
        os << rowScaledColNorms_h[k];
        if (k + LO (1) < lclNumCols) {
          os << ", ";
        }
      }
      os << "]" << endl;
      out << os.str ();

      bool ok = true;
      for (LO lclCol = 0; lclCol < lclNumCols; ++lclCol) {
        if (! near (rowScaledColNorms_h[lclCol],
                    test.gblRowScaledColNorms[lclCol],
                    toleranceFactor)) {
          ok = false;
          break;
        }
      }
      TEST_ASSERT( ok );
    } // test global row-scaled column norms

    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY( gblSuccess, 1 );
    if (gblSuccess != 1) {
      out << "Test FAILED on some process!" << endl;
      return;
    }
  } // test computeLocalRowAndColumnNorms (CrsMatrix)

  // Test computeLocalRowAndColumnNorms (RowMatrix)
  {
    out << "Test computeLocalRowAndColumnOneNorms (RowMatrix)" << endl;
    Teuchos::OSTab tab2 (out);

    using Tpetra::Details::computeLocalRowAndColumnOneNorms_RowMatrix;
    auto result1 =
      computeLocalRowAndColumnOneNorms_RowMatrix (* (test.A), assumeSymmetric);

    {
      out << "Test detection of local error conditions" << endl;
      Teuchos::OSTab tab3 (out);

      TEST_EQUALITY( test.lclFoundInf, result1.foundInf );
      TEST_EQUALITY( test.lclFoundNan, result1.foundNan );
      TEST_EQUALITY( test.lclFoundZeroDiag, result1.foundZeroDiag );
      TEST_EQUALITY( test.lclFoundZeroRowNorm, result1.foundZeroRowNorm );
    }

    {
      out << "Test local row norms" << endl;
      Teuchos::OSTab tab3 (out);
      auto rowNorms_h = Kokkos::create_mirror_view (result1.rowNorms);
      Kokkos::deep_copy (rowNorms_h, result1.rowNorms);

      std::ostringstream os;
      os << "Expected local rowNorms: [";
      for (LO k = 0; k < lclNumRows; ++k) {
        os << test.lclRowNorms[k];
        if (k + LO (1) < lclNumRows) {
          os << ", ";
        }
      }
      os << "]" << endl << "Actual local rowNorms: [";
      for (LO k = 0; k < lclNumRows; ++k) {
        os << rowNorms_h[k];
        if (k + LO (1) < lclNumRows) {
          os << ", ";
        }
      }
      os << "]" << endl;
      out << os.str ();

      bool ok = true;
      for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
        if (! near (rowNorms_h[lclRow],
                    test.lclRowNorms[lclRow],
                    toleranceFactor)) {
          ok = false;
          break;
        }
      }
      TEST_ASSERT( ok );
    } // test local row norms

    {
      out << "Test local row diagonal entries" << endl;
      Teuchos::OSTab tab3 (out);
      auto rowDiagonalEntries_h =
        Kokkos::create_mirror_view (result1.rowDiagonalEntries);
      Kokkos::deep_copy (rowDiagonalEntries_h, result1.rowDiagonalEntries);

      bool ok = true;
      for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
        if (! near (rowDiagonalEntries_h[lclRow],
                    test.lclRowDiagonalEntries[lclRow],
                    toleranceFactor)) {
          ok = false;
          break;
        }
      }
      TEST_ASSERT( ok );
    } // test local row diagonal entries

    if (! assumeSymmetric) {
      out << "Test local column norms" << endl;
      Teuchos::OSTab tab3 (out);
      auto colNorms_h = Kokkos::create_mirror_view (result1.colNorms);
      Kokkos::deep_copy (colNorms_h, result1.colNorms);

      bool ok = true;
      for (LO lclCol = 0; lclCol < lclNumCols; ++lclCol) {
        if (! near (colNorms_h[lclCol],
                    test.lclColNorms[lclCol],
                    toleranceFactor)) {
          ok = false;
          break;
        }
      }
      TEST_ASSERT( ok );
    } // ! assumeSymmetric

    if (! assumeSymmetric) {
      out << "Test local column diagonal entries" << endl;
      Teuchos::OSTab tab3 (out);
      auto colDiagonalEntries_h =
        Kokkos::create_mirror_view (result1.colDiagonalEntries);
      Kokkos::deep_copy (colDiagonalEntries_h, result1.colDiagonalEntries);

      bool ok = true;
      for (LO lclCol = 0; lclCol < lclNumCols; ++lclCol) {
        if (! near (colDiagonalEntries_h[lclCol],
                    test.lclColDiagonalEntries[lclCol],
                    toleranceFactor)) {
          ok = false;
          break;
        }
      }
      TEST_ASSERT( ok );
    } // ! assumeSymmetric

    // We've already tested globalize{Row,Column}OneNorms above;
    // neither depends on whether the matrix is a CrsMatrix.

    reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
    TEST_EQUALITY( gblSuccess, 1 );
    if (gblSuccess != 1) {
      out << "Test FAILED on some process!" << endl;
      return;
    }
  } // test computeLocalRowAndColumnNorms (RowMatrix)

  // Test computeRowAndColumnOneNorms
  {
    out << "Test computeRowAndColumnOneNorms" << endl;
    Teuchos::OSTab tab2 (out);
    // FIXME (mfh 24 May 2018) Why the cast?
    auto result2 =
      Tpetra::computeRowAndColumnOneNorms (static_cast<row_matrix_type&> (* (test.A)),
                                           assumeSymmetric);
    {
      out << "Test whether assumeSymmetric got communicated" << endl;
      Teuchos::OSTab tab3 (out);
      TEST_EQUALITY( assumeSymmetric, result2.assumeSymmetric );
    }

    {
      out << "Test detection of global error conditions" << endl;
      Teuchos::OSTab tab3 (out);

      TEST_EQUALITY( test.gblFoundInf, result2.foundInf );
      TEST_EQUALITY( test.gblFoundNan, result2.foundNan );
      TEST_EQUALITY( test.gblFoundZeroDiag, result2.foundZeroDiag );
      TEST_EQUALITY( test.gblFoundZeroRowNorm, result2.foundZeroRowNorm );
    }

    {
      out << "Test global row norms" << endl;
      Teuchos::OSTab tab3 (out);
      auto rowNorms_h = Kokkos::create_mirror_view (result2.rowNorms);
      Kokkos::deep_copy (rowNorms_h, result2.rowNorms);

      bool ok = true;
      for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
        if (! near (rowNorms_h[lclRow],
                    test.gblRowNorms[lclRow],
                    toleranceFactor)) {
          ok = false;
          break;
        }
      }
      TEST_ASSERT( ok );
    } // test global row norms

    {
      out << "Test global column norms" << endl;
      Teuchos::OSTab tab3 (out);
      auto colNorms_h = Kokkos::create_mirror_view (result2.colNorms);
      Kokkos::deep_copy (colNorms_h, result2.colNorms);

      bool ok = true;
      for (LO lclCol = 0; lclCol < lclNumCols; ++lclCol) {
        if (! near (colNorms_h[lclCol],
                    test.gblColNorms[lclCol],
                    toleranceFactor)) {
          ok = false;
          break;
        }
      }
      TEST_ASSERT( ok );
    } // test global column norms

    if (! assumeSymmetric) {
      out << "Test global row-scaled column norms" << endl;
      Teuchos::OSTab tab3 (out);
      auto rowScaledColNorms_h =
        Kokkos::create_mirror_view (result2.rowScaledColNorms);
      Kokkos::deep_copy (rowScaledColNorms_h, result2.rowScaledColNorms);

      std::ostringstream os;
      os << "Expected global rowScaledColNorms: [";
      for (LO k = 0; k < lclNumCols; ++k) {
        os << test.gblRowScaledColNorms[k];
        if (k + LO (1) < lclNumCols) {
          os << ", ";
        }
      }
      os << "]" << endl << "Actual global rowScaledColNorms: [";
      for (LO k = 0; k < lclNumCols; ++k) {
        os << rowScaledColNorms_h[k];
        if (k + LO (1) < lclNumCols) {
          os << ", ";
        }
      }
      os << "]" << endl;
      out << os.str ();

      bool ok = true;
      for (LO lclCol = 0; lclCol < lclNumCols; ++lclCol) {
        if (! near (rowScaledColNorms_h[lclCol],
                    test.gblRowScaledColNorms[lclCol],
                    toleranceFactor)) {
          ok = false;
          break;
        }
      }
      TEST_ASSERT( ok );
    } // test global row-scaled column norms

    {
      out << "Test deepCopyFillCompleteCrsMatrix" << endl;
      Teuchos::OSTab tab3 (out);
      auto A_copy = deepCopyFillCompleteCrsMatrix (*(test.A));
      testCrsMatrixEquality (success, out, *(test.A), *A_copy);
    }

    Tpetra::Vector<mag_type, LO, GO, NT> rowNormsVec =
      createVectorFromCopyOf1DView (test.A->getRowMap (), result2.rowNorms);
    Tpetra::Vector<mag_type, LO, GO, NT> colNormsVec =
      createVectorFromCopyOf1DView (test.A->getColMap (),
                                    result2.assumeSymmetric ?
                                      result2.colNorms :
                                      result2.rowScaledColNorms);
    {
      out << "Test left-scaling CrsMatrix with Kokkos::View" << endl;
      Teuchos::OSTab tab3 (out);
      auto A_copy = deepCopyFillCompleteCrsMatrix (*(test.A));
      Tpetra::leftAndOrRightScaleCrsMatrix (*A_copy,
                                            result2.rowNorms,
                                            result2.assumeSymmetric ?
                                              result2.colNorms :
                                              result2.rowScaledColNorms, // ignored
                                            true, false,
                                            result2.assumeSymmetric,
                                            Tpetra::SCALING_DIVIDE);
      testCrsMatrixEquality (success, out, *(test.A_leftScaled), *A_copy);
    }
    {
      out << "Test left-scaling CrsMatrix with Vector" << endl;
      Teuchos::OSTab tab3 (out);
      auto A_copy = deepCopyFillCompleteCrsMatrix (*(test.A));

      Tpetra::leftAndOrRightScaleCrsMatrix (*A_copy,
                                            rowNormsVec,
                                            colNormsVec, // ignored
                                            true, false,
                                            result2.assumeSymmetric,
                                            Tpetra::SCALING_DIVIDE);
      testCrsMatrixEquality (success, out, *(test.A_leftScaled), *A_copy);
    }
    {
      out << "Test right-scaling CrsMatrix with Kokkos::View" << endl;
      Teuchos::OSTab tab3 (out);
      auto A_copy = deepCopyFillCompleteCrsMatrix (*(test.A));
      Tpetra::leftAndOrRightScaleCrsMatrix (*A_copy,
                                            result2.rowNorms, // ignored
                                            result2.assumeSymmetric ?
                                              result2.colNorms :
                                              result2.rowScaledColNorms,
                                            false, true,
                                            result2.assumeSymmetric,
                                            Tpetra::SCALING_DIVIDE);
      testCrsMatrixEquality (success, out, *(test.A_rightScaled), *A_copy);
    }
    {
      out << "Test right-scaling CrsMatrix with Tpetra::Vector" << endl;
      Teuchos::OSTab tab3 (out);
      auto A_copy = deepCopyFillCompleteCrsMatrix (*(test.A));
      Tpetra::leftAndOrRightScaleCrsMatrix (*A_copy,
                                            rowNormsVec, // ignored
                                            colNormsVec,
                                            false, true,
                                            result2.assumeSymmetric,
                                            Tpetra::SCALING_DIVIDE);
      testCrsMatrixEquality (success, out, *(test.A_rightScaled), *A_copy);
    }
  } // test computeRowAndColumnOneNorms

  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  TEST_EQUALITY( gblSuccess, 1 );
  if (gblSuccess != 1) {
    out << "Test FAILED on some process!" << endl;
  }
}

// NOTE (mfh 23 May 2018) This function is only for the test.
template<class SC, class LO, class GO, class NT>
EquilibrationTest<SC, LO, GO, NT>
makeSymmetricPositiveDefiniteTridiagonalMatrixTest (Teuchos::FancyOStream& out,
                                                    const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                                                    const bool assumeSymmetric)
{
  using Teuchos::RCP;
  using std::endl;
  using GST = Tpetra::global_size_t;
  using map_type = Tpetra::Map<LO, GO, NT>;
  using crs_graph_type = Tpetra::CrsGraph<LO, GO, NT>;
  using crs_matrix_type = Tpetra::CrsMatrix<SC, LO, GO, NT>;
  using KAT = Kokkos::ArithTraits<SC>;
  using val_type = typename KAT::val_type;
  using mag_type = typename Kokkos::ArithTraits<val_type>::mag_type;
  using KAM = Kokkos::ArithTraits<mag_type>;

  const LO lclNumRows = 5;
  const GO gblNumRows =
    static_cast<GO> (comm->getSize ()) * static_cast<GO> (lclNumRows);
  const GO indexBase = 0;

  out << "Create symmetric positive definite tridiagonal matrix problem"
      << endl;
  Teuchos::OSTab tab0 (out);

  out << "Create Maps" << endl;
  RCP<const map_type> rowMap =
    rcp (new map_type (static_cast<GST> (gblNumRows),
                       static_cast<size_t> (lclNumRows),
                       indexBase, comm));
  RCP<const map_type> domMap = rowMap;
  RCP<const map_type> ranMap = rowMap;

  out << "Create CrsGraph" << endl;
  const size_t maxNumEntPerRow = 3;
  RCP<crs_graph_type> G =
    rcp (new crs_graph_type (rowMap, maxNumEntPerRow, Tpetra::StaticProfile));
  std::vector<GO> globalIndices (maxNumEntPerRow);

  for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
    const GO gblRow = rowMap->getGlobalElement (lclRow);

    if (gblRow == 0) {
      const LO numEnt = 2;
      globalIndices[0] = gblRow;
      globalIndices[1] = gblRow+1;
      G->insertGlobalIndices (gblRow, numEnt, globalIndices.data ());
    }
    else if (gblRow == gblNumRows - GO (1)) {
      const LO numEnt = 2;
      globalIndices[0] = gblRow-1;
      globalIndices[1] = gblRow;
      G->insertGlobalIndices (gblRow, numEnt, globalIndices.data ());
    }
    else {
      const LO numEnt = 3;
      globalIndices[0] = gblRow-1;
      globalIndices[1] = gblRow;
      globalIndices[2] = gblRow+1;
      G->insertGlobalIndices (gblRow, numEnt, globalIndices.data ());
    }
  }
  G->fillComplete (domMap, ranMap);

  const SC diagVal {2.0};
  const SC offDiagVal {-1.0};
  std::vector<SC> firstRowValues {diagVal, offDiagVal};
  std::vector<SC> middleRowValues {offDiagVal, diagVal, offDiagVal};
  std::vector<SC> lastRowValues {offDiagVal, diagVal};

  out << "Create test CrsMatrix" << endl;
  RCP<crs_matrix_type> A = rcp (new crs_matrix_type (G));
  for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
    const GO gblRow = rowMap->getGlobalElement (lclRow);

    if (gblRow == 0) {
      const LO numEnt = 2;
      globalIndices[0] = gblRow;
      globalIndices[1] = gblRow+1;
      A->replaceGlobalValues (gblRow, numEnt, firstRowValues.data (), globalIndices.data ());
    }
    else if (gblRow == gblNumRows - GO (1)) {
      const LO numEnt = 2;
      globalIndices[0] = gblRow-1;
      globalIndices[1] = gblRow;
      A->replaceGlobalValues (gblRow, numEnt, lastRowValues.data (), globalIndices.data ());
    }
    else {
      const LO numEnt = 3;
      globalIndices[0] = gblRow-1;
      globalIndices[1] = gblRow;
      globalIndices[2] = gblRow+1;
      A->replaceGlobalValues (gblRow, numEnt, middleRowValues.data (), globalIndices.data ());
    }
  }
  A->fillComplete (domMap, ranMap);

  RCP<const map_type> colMap = G->getColMap ();
  const LO lclNumCols = static_cast<LO> (colMap->getNodeNumElements ());
  const GO gblNumCols = static_cast<GO> (G->getDomainMap ()->getGlobalNumElements ());

  const mag_type diagAbsVal = KAT::abs (diagVal);
  const mag_type offDiagAbsVal = KAT::abs (offDiagVal);

  out << "Compute local row norms and diagonal entries" << endl;
  std::vector<val_type> lclRowDiagonalEntries (lclNumRows);
  std::vector<mag_type> lclRowNorms (lclNumRows);
  for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
    const GO gblRow = rowMap->getGlobalElement (lclRow);
    const mag_type expectedRowNorm =
      (gblRow == 0 || gblRow == gblNumRows - GO (1)) ?
      offDiagAbsVal + diagAbsVal :
      offDiagAbsVal + diagAbsVal + offDiagAbsVal;

    lclRowDiagonalEntries[lclRow] = diagVal;
    lclRowNorms[lclRow] = expectedRowNorm;
  }

  // For this matrix, the global row norms are the same as the local
  // row norms, since the matrix's row Map and range Map are the same.
  // This may not be the case for matrices with an overlapping or
  // permuted row Map.
  out << "Compute global row norms" << endl;
  std::vector<mag_type> gblRowNorms (lclRowNorms.begin (), lclRowNorms.end ());

  Teuchos::RCP<crs_matrix_type> A_leftScaled = [&] () {
    Teuchos::RCP<crs_matrix_type> A_copy = deepCopyFillCompleteCrsMatrix (*A);
    A_copy->resumeFill ();

    auto A_lcl = A_copy->getLocalMatrix ();
    auto val_h = Kokkos::create_mirror_view (A_lcl.values);
    Kokkos::deep_copy (val_h, A_lcl.values);
    auto ptr_h = Kokkos::create_mirror_view (A_lcl.graph.row_map);
    Kokkos::deep_copy (ptr_h, A_lcl.graph.row_map);
    auto ind_h = Kokkos::create_mirror_view (A_lcl.graph.entries);
    Kokkos::deep_copy (ind_h, A_lcl.graph.entries);

    for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
      const GO gblRow = rowMap->getGlobalElement (lclRow);

      using size_type = typename decltype (A_lcl.graph)::size_type;
      for (size_type k = ptr_h[lclRow]; k < ptr_h[lclRow+1]; ++k) {
        const LO lclCol = ind_h[k];
        const GO gblCol = colMap->getGlobalElement (lclCol);
        const val_type expectedUnscaledVal =
          (gblRow == gblCol) ? diagVal : offDiagVal;
        const mag_type rowNorm = gblRowNorms[lclRow];
        const mag_type scalingFactor = assumeSymmetric ?
          KAM::sqrt (rowNorm) : rowNorm;
        val_h[k] = expectedUnscaledVal / scalingFactor;
      }
    }
    Kokkos::deep_copy (A_lcl.values, val_h);
    A_copy->fillComplete (A_copy->getDomainMap (), A_copy->getRangeMap ());
    return A_copy;
  } ();

  if (assumeSymmetric) {
    out << "assumeSymmetric=true: Skip local (column norms, "
      "diagonal entries, and row-scaled column norms)" << endl;
  }
  else {
    out << "assumeSymmetric=false: Compute local (column norms, "
      "diagonal entries, and row-scaled column norms)" << endl;
  }

  std::vector<mag_type> lclColNorms
    (assumeSymmetric ? LO (0) : lclNumCols);
  std::vector<val_type> lclColDiagonalEntries
    (assumeSymmetric ? LO (0) : lclNumCols);
  std::vector<mag_type> lclRowScaledColNorms
    (assumeSymmetric ? LO (0) : lclNumCols);

  if (! assumeSymmetric) {
    // Columns are a little more complicated, since in the usual
    // distributed case, the local column norms may not be the same as
    // the global column norms.
    for (LO lclCol = 0; lclCol < lclNumCols; ++lclCol) {
      const GO gblCol = colMap->getGlobalElement (lclCol);
      const GO gblRow = gblCol;

      mag_type expectedColNorm {0.0};
      val_type expectedColDiagonalVal {0.0};

      // These are local column norms, not global column norms.
      if (rowMap->isNodeGlobalElement (gblRow-1)) {
        expectedColNorm += offDiagAbsVal;
      }
      if (rowMap->isNodeGlobalElement (gblRow)) {
        expectedColNorm += diagAbsVal;
        expectedColDiagonalVal += diagVal;
      }
      if (rowMap->isNodeGlobalElement (gblRow+1)) {
        expectedColNorm += offDiagAbsVal;
      }

      lclColNorms[lclCol] = expectedColNorm;
      lclColDiagonalEntries[lclCol] = expectedColDiagonalVal;
    } // for each local column

    for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
      const GO gblRow = rowMap->getGlobalElement (lclRow);
      const GO gblCol = gblRow;
      mag_type lclRowScaledColNorm {0.0};

      if (gblRow == rowMap->getMinAllGlobalIndex ()) {
        const LO diagLclColInd = colMap->getLocalElement (gblCol);
        const LO rightLclColInd = colMap->getLocalElement (gblCol + GO (1));

        if (diagLclColInd != Tpetra::Details::OrdinalTraits<LO>::invalid ()) {
          lclRowScaledColNorms[diagLclColInd] += diagAbsVal / gblRowNorms[lclRow];
        }
        if (rightLclColInd != Tpetra::Details::OrdinalTraits<LO>::invalid ()) {
          lclRowScaledColNorms[rightLclColInd] += offDiagAbsVal / gblRowNorms[lclRow];
        }
      }
      else if (gblRow == rowMap->getMaxAllGlobalIndex ()) {
        const LO leftLclColInd = colMap->getLocalElement (gblCol - GO (1));
        const LO diagLclColInd = colMap->getLocalElement (gblCol);

        if (leftLclColInd != Tpetra::Details::OrdinalTraits<LO>::invalid ()) {
          lclRowScaledColNorms[leftLclColInd] += offDiagAbsVal / gblRowNorms[lclRow];
        }
        if (diagLclColInd != Tpetra::Details::OrdinalTraits<LO>::invalid ()) {
          lclRowScaledColNorms[diagLclColInd] += diagAbsVal / gblRowNorms[lclRow];
        }
      }
      else {
        const LO leftLclColInd = colMap->getLocalElement (gblCol - GO (1));
        const LO diagLclColInd = colMap->getLocalElement (gblCol);
        const LO rightLclColInd = colMap->getLocalElement (gblCol + GO (1));

        if (leftLclColInd != Tpetra::Details::OrdinalTraits<LO>::invalid ()) {
          lclRowScaledColNorms[leftLclColInd] += offDiagAbsVal / gblRowNorms[lclRow];
        }
        if (diagLclColInd != Tpetra::Details::OrdinalTraits<LO>::invalid ()) {
          lclRowScaledColNorms[diagLclColInd] += diagAbsVal / gblRowNorms[lclRow];
        }
        if (rightLclColInd != Tpetra::Details::OrdinalTraits<LO>::invalid ()) {
          lclRowScaledColNorms[rightLclColInd] += offDiagAbsVal / gblRowNorms[lclRow];
        }
      }
    }
  } // ! assumeSymmetric

  out << "Compute global column norms" << endl;
  std::vector<mag_type> gblColNorms (lclNumCols);
  // The matrix is symmetric, so this holds regardless of assumeSymmetric.
  for (LO lclCol = 0; lclCol < lclNumCols; ++lclCol) {
    const GO gblCol = colMap->getGlobalElement (lclCol);

    mag_type expectedColNorm {0.0};
    if (gblCol == 0 || gblCol == gblNumCols - GO (1)) {
      expectedColNorm = diagAbsVal + offDiagAbsVal;
    }
    else {
      expectedColNorm = offDiagAbsVal + diagAbsVal + offDiagAbsVal;
    }
    gblColNorms[lclCol] = expectedColNorm;
  }

  if (assumeSymmetric) {
    out << "Skip global row-scaled column norms" << endl;
  }
  else {
    out << "Compute global row-scaled column norms" << endl;
  }
  std::vector<mag_type> gblRowScaledColNorms
    (assumeSymmetric ? LO (0) : lclNumCols);
  if (! assumeSymmetric && colMap->getGlobalNumElements () != 0) {
    const GO gblMinGblCol = domMap->getMinAllGlobalIndex ();
    const GO gblMaxGblCol = domMap->getMaxAllGlobalIndex ();

    for (LO lclCol = 0; lclCol < lclNumCols; ++lclCol) {
      const GO gblCol = colMap->getGlobalElement (lclCol);
      if (gblCol == gblMinGblCol || gblCol == gblMaxGblCol) {
        gblRowScaledColNorms[lclCol] = 11.0 / 12.0; // 2/3 + 1/4
      }
      else if (gblCol == gblMinGblCol + GO (1) ||
               gblCol == gblMaxGblCol - GO (1)) {
        gblRowScaledColNorms[lclCol] = 13.0 / 12.0; // 1/3 + 1/2 + 1/4
      }
      else {
        gblRowScaledColNorms[lclCol] = 1.0; // 1/4 + 1/2 + 1/4
      }
    }
  }

  out << "Compute right-scaled matrix" << endl;
  Teuchos::RCP<crs_matrix_type> A_rightScaled = [&] () {
    Teuchos::RCP<crs_matrix_type> A_copy = deepCopyFillCompleteCrsMatrix (*A);
    A_copy->resumeFill ();

    auto A_lcl = A_copy->getLocalMatrix ();
    auto val_h = Kokkos::create_mirror_view (A_lcl.values);
    Kokkos::deep_copy (val_h, A_lcl.values);
    auto ptr_h = Kokkos::create_mirror_view (A_lcl.graph.row_map);
    Kokkos::deep_copy (ptr_h, A_lcl.graph.row_map);
    auto ind_h = Kokkos::create_mirror_view (A_lcl.graph.entries);
    Kokkos::deep_copy (ind_h, A_lcl.graph.entries);

    for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
      using size_type = typename decltype (A_lcl.graph)::size_type;
      for (size_type k = ptr_h[lclRow]; k < ptr_h[lclRow+1]; ++k) {
        const LO lclCol = ind_h[k];
        const mag_type colNorm = assumeSymmetric ?
          gblColNorms[lclCol] : gblRowScaledColNorms[lclCol];
        // const mag_type rowScalingFactor = assumeSymmetric ?
        //   KAM::sqrt (gblRowNorms[lclRow]) : gblRowNorms[lclRow];
        const mag_type colScalingFactor = assumeSymmetric ?
          KAM::sqrt (colNorm) : colNorm;
        // const mag_type scalingFactor = rowScalingFactor * colScalingFactor;
        const mag_type scalingFactor = colScalingFactor;
        val_h[k] = val_h[k] / scalingFactor;
      }
    }
    Kokkos::deep_copy (A_lcl.values, val_h);
    A_copy->fillComplete (A_copy->getDomainMap (), A_copy->getRangeMap ());
    return A_copy;
  } ();

  const bool lclFoundInf = false;
  const bool lclFoundNan = false;
  const bool lclFoundZeroDiag = false;
  const bool lclFoundZeroRowNorm = false;
  const bool gblFoundInf = false;
  const bool gblFoundNan = false;
  const bool gblFoundZeroDiag = false;
  const bool gblFoundZeroRowNorm = false;

  return EquilibrationTest<SC, LO, GO, NT> {
    A,
    lclRowNorms,
    lclRowDiagonalEntries,
    lclColNorms,
    lclColDiagonalEntries,
    lclRowScaledColNorms,
    gblRowNorms,
    gblColNorms,
    gblRowScaledColNorms,
    A_leftScaled,
    A_rightScaled,
    lclFoundInf,
    lclFoundNan,
    lclFoundZeroDiag,
    lclFoundZeroRowNorm,
    gblFoundInf,
    gblFoundNan,
    gblFoundZeroDiag,
    gblFoundZeroRowNorm
  };
}

// NOTE (mfh 23 May 2018) This function is only for the test.
template<class SC, class LO, class GO, class NT>
EquilibrationTest<SC, LO, GO, NT>
makeMatrixTestWithExplicitZeroDiag (Teuchos::FancyOStream& out,
                                    const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                                    const bool assumeSymmetric)
{
  using Teuchos::RCP;
  using std::endl;
  using GST = Tpetra::global_size_t;
  using map_type = Tpetra::Map<LO, GO, NT>;
  using crs_graph_type = Tpetra::CrsGraph<LO, GO, NT>;
  using crs_matrix_type = Tpetra::CrsMatrix<SC, LO, GO, NT>;
  using KAT = Kokkos::ArithTraits<SC>;
  using val_type = typename KAT::val_type;
  using mag_type = typename Kokkos::ArithTraits<val_type>::mag_type;
  using KAM = Kokkos::ArithTraits<mag_type>;

  const int myRank = comm->getRank ();
  const int numProcs = comm->getSize ();

  const LO lclNumRows = 1;
  const GO gblNumRows = GO (numProcs) * GO (lclNumRows);
  const GO indexBase = 0;

  out << "Create diagonal matrix problem with explicit zero "
    "diagonal entry on one process, not Process 0" << endl;
  Teuchos::OSTab tab0 (out);

  out << "Create Maps" << endl;
  RCP<const map_type> rowMap =
    rcp (new map_type (GST (gblNumRows), std::size_t (lclNumRows),
                       indexBase, comm));
  RCP<const map_type> domMap = rowMap;
  RCP<const map_type> ranMap = rowMap;

  out << "Create CrsGraph" << endl;
  const size_t maxNumEntPerRow = 1;
  RCP<crs_graph_type> G =
    rcp (new crs_graph_type (rowMap, maxNumEntPerRow, Tpetra::StaticProfile));
  std::vector<GO> globalIndices (maxNumEntPerRow);

  for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
    const GO gblRow = rowMap->getGlobalElement (lclRow);
    const LO numEnt = 1;
    globalIndices[0] = gblRow;
    G->insertGlobalIndices (gblRow, numEnt, globalIndices.data ());
  }
  G->fillComplete (domMap, ranMap);

  out << "Create test CrsMatrix" << endl;
  RCP<crs_matrix_type> A = rcp (new crs_matrix_type (G));
  std::vector<SC> values (maxNumEntPerRow);
  for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
    const GO gblRow = rowMap->getGlobalElement (lclRow);

    const LO numEnt = 1;
    globalIndices[0] = gblRow;
    if (myRank == 1) {
      values[0] = 0.0; // explicitly stored zero diagonal
    }
    else {
      values[0] = 1.0;
    }
    A->replaceGlobalValues (gblRow, numEnt, values.data (), globalIndices.data ());
  }
  A->fillComplete (domMap, ranMap);

  RCP<const map_type> colMap = G->getColMap ();
  const LO lclNumCols = static_cast<LO> (colMap->getNodeNumElements ());

  out << "Compute local row norms and diagonal entries" << endl;
  std::vector<val_type> lclRowDiagonalEntries (lclNumRows);
  std::vector<mag_type> lclRowNorms (lclNumRows);
  for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
    const mag_type expectedRowNorm = myRank == 1 ? 0.0 : 1.0;
    lclRowDiagonalEntries[lclRow] = myRank == 1 ? 0.0 : 1.0;
    lclRowNorms[lclRow] = expectedRowNorm;
  }

  // For this matrix, the global row norms are the same as the local
  // row norms, since the matrix's row Map and range Map are the same.
  // This may not be the case for matrices with an overlapping or
  // permuted row Map.
  out << "Compute global row norms" << endl;
  std::vector<mag_type> gblRowNorms (lclRowNorms.begin (), lclRowNorms.end ());

  Teuchos::RCP<crs_matrix_type> A_leftScaled = [&] () {
    Teuchos::RCP<crs_matrix_type> A_copy = deepCopyFillCompleteCrsMatrix (*A);
    A_copy->resumeFill ();

    auto A_lcl = A_copy->getLocalMatrix ();
    auto val_h = Kokkos::create_mirror_view (A_lcl.values);
    Kokkos::deep_copy (val_h, A_lcl.values);
    auto ptr_h = Kokkos::create_mirror_view (A_lcl.graph.row_map);
    Kokkos::deep_copy (ptr_h, A_lcl.graph.row_map);
    auto ind_h = Kokkos::create_mirror_view (A_lcl.graph.entries);
    Kokkos::deep_copy (ind_h, A_lcl.graph.entries);

    for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
      using size_type = typename decltype (A_lcl.graph)::size_type;
      for (size_type k = ptr_h[lclRow]; k < ptr_h[lclRow+1]; ++k) {
        const val_type expectedUnscaledVal = myRank == 1 ? 0.0 : 1.0;
        const mag_type rowNorm = gblRowNorms[lclRow];
        const mag_type scalingFactor = assumeSymmetric ?
          KAM::sqrt (rowNorm) : rowNorm;
        val_h[k] = expectedUnscaledVal / scalingFactor; // may be NaN in this case
      }
    }
    Kokkos::deep_copy (A_lcl.values, val_h);
    A_copy->fillComplete (A_copy->getDomainMap (), A_copy->getRangeMap ());
    return A_copy;
  } ();

  if (assumeSymmetric) {
    out << "assumeSymmetric=true: Skip local (column norms, "
      "diagonal entries, and row-scaled column norms)" << endl;
  }
  else {
    out << "assumeSymmetric=false: Compute local (column norms, "
      "diagonal entries, and row-scaled column norms)" << endl;
  }

  std::vector<mag_type> lclColNorms
    (assumeSymmetric ? LO (0) : lclNumCols);
  std::vector<val_type> lclColDiagonalEntries
    (assumeSymmetric ? LO (0) : lclNumCols);
  std::vector<mag_type> lclRowScaledColNorms
    (assumeSymmetric ? LO (0) : lclNumCols);

  if (! assumeSymmetric) {
    // Columns are a little more complicated, since in the usual
    // distributed case, the local column norms may not be the same as
    // the global column norms.
    for (LO lclCol = 0; lclCol < lclNumCols; ++lclCol) {
      lclColNorms[lclCol] = myRank == 1 ? 0.0 : 1.0;
      lclColDiagonalEntries[lclCol] = myRank == 1 ? 0.0 : 1.0;
    } // for each local column

    for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
      const GO gblRow = rowMap->getGlobalElement (lclRow);
      const GO gblCol = gblRow;
      const LO lclCol = colMap->getLocalElement (gblCol);
      mag_type lclRowScaledColNorm {0.0};

      const mag_type matrixAbsVal = myRank == 1 ? 0.0 : 1.0;
      lclRowScaledColNorms[lclCol] = myRank == 1 ?
        NaughtyValues<mag_type>::quiet_NaN () : 1.0;
    }
  } // ! assumeSymmetric

  out << "Compute global column norms" << endl;
  std::vector<mag_type> gblColNorms (lclNumCols);
  // The matrix is symmetric, so this holds regardless of assumeSymmetric.
  for (LO lclCol = 0; lclCol < lclNumCols; ++lclCol) {
    gblColNorms[lclCol] = myRank == 1 ? 0.0 : 1.0;
  }

  if (assumeSymmetric) {
    out << "Skip global row-scaled column norms" << endl;
  }
  else {
    out << "Compute global row-scaled column norms" << endl;
  }
  std::vector<mag_type> gblRowScaledColNorms
    (assumeSymmetric ? LO (0) : lclNumCols);
  if (! assumeSymmetric && colMap->getGlobalNumElements () != 0) {
    for (LO lclCol = 0; lclCol < lclNumCols; ++lclCol) {
      gblRowScaledColNorms[lclCol] = myRank == 1 ?
        NaughtyValues<mag_type>::quiet_NaN () : 1.0;
    }
  }

  Teuchos::RCP<crs_matrix_type> A_rightScaled = [&] () {
    Teuchos::RCP<crs_matrix_type> A_copy = deepCopyFillCompleteCrsMatrix (*A);
    A_copy->resumeFill ();

    auto A_lcl = A_copy->getLocalMatrix ();
    auto val_h = Kokkos::create_mirror_view (A_lcl.values);
    Kokkos::deep_copy (val_h, A_lcl.values);
    auto ptr_h = Kokkos::create_mirror_view (A_lcl.graph.row_map);
    Kokkos::deep_copy (ptr_h, A_lcl.graph.row_map);
    auto ind_h = Kokkos::create_mirror_view (A_lcl.graph.entries);
    Kokkos::deep_copy (ind_h, A_lcl.graph.entries);

    for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
      using size_type = typename decltype (A_lcl.graph)::size_type;
      for (size_type k = ptr_h[lclRow]; k < ptr_h[lclRow+1]; ++k) {
        const LO lclCol = ind_h[k];
        const mag_type colNorm = assumeSymmetric ?
          gblColNorms[lclCol] : gblRowScaledColNorms[lclCol];
        // const mag_type rowScalingFactor = assumeSymmetric ?
        //   KAM::sqrt (gblRowNorms[lclRow]) : gblRowNorms[lclRow];
        const mag_type colScalingFactor = assumeSymmetric ?
          KAM::sqrt (colNorm) : colNorm;
        // const mag_type scalingFactor = rowScalingFactor * colScalingFactor;
        const mag_type scalingFactor = colScalingFactor;
        val_h[k] = val_h[k] / scalingFactor;
      }
    }
    Kokkos::deep_copy (A_lcl.values, val_h);
    A_copy->fillComplete (A_copy->getDomainMap (), A_copy->getRangeMap ());
    return A_copy;
  } ();

  const bool lclFoundInf = false;
  const bool lclFoundNan = false;
  const bool lclFoundZeroDiag = myRank == 1 ? true : false;
  const bool lclFoundZeroRowNorm = myRank == 1 ? true : false;
  const bool gblFoundInf = false;
  const bool gblFoundNan = false;
  const bool gblFoundZeroDiag = numProcs == 1 ? false : true;
  const bool gblFoundZeroRowNorm = numProcs == 1 ? false : true;

  return EquilibrationTest<SC, LO, GO, NT> {
    A,
    lclRowNorms,
    lclRowDiagonalEntries,
    lclColNorms,
    lclColDiagonalEntries,
    lclRowScaledColNorms,
    gblRowNorms,
    gblColNorms,
    gblRowScaledColNorms,
    A_leftScaled,
    A_rightScaled,
    lclFoundInf,
    lclFoundNan,
    lclFoundZeroDiag,
    lclFoundZeroRowNorm,
    gblFoundInf,
    gblFoundNan,
    gblFoundZeroDiag,
    gblFoundZeroRowNorm
  };
}

// NOTE (mfh 14 Jun 2018) This function is only for the test.
template<class SC, class LO, class GO, class NT>
EquilibrationTest<SC, LO, GO, NT>
makeMatrixTestWithImplicitZeroDiag (Teuchos::FancyOStream& out,
                                    const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                                    const bool assumeSymmetric)
{
  using Teuchos::RCP;
  using std::endl;
  using GST = Tpetra::global_size_t;
  using map_type = Tpetra::Map<LO, GO, NT>;
  using crs_graph_type = Tpetra::CrsGraph<LO, GO, NT>;
  using crs_matrix_type = Tpetra::CrsMatrix<SC, LO, GO, NT>;
  using KAT = Kokkos::ArithTraits<SC>;
  using val_type = typename KAT::val_type;
  using mag_type = typename Kokkos::ArithTraits<val_type>::mag_type;
  using KAM = Kokkos::ArithTraits<mag_type>;

  const int myRank = comm->getRank ();
  const int numProcs = comm->getSize ();

  const LO lclNumRows = 1;
  const GO gblNumRows = GO (numProcs) * GO (lclNumRows);
  const GO indexBase = 0;

  out << "Create diagonal matrix problem with implicit zero "
    "diagonal entry on one process, not Process 0" << endl;
  Teuchos::OSTab tab0 (out);

  out << "Create Maps" << endl;
  RCP<const map_type> rowMap =
    rcp (new map_type (GST (gblNumRows), std::size_t (lclNumRows),
                       indexBase, comm));
  RCP<const map_type> domMap = rowMap;
  RCP<const map_type> ranMap = rowMap;

  out << "Create CrsGraph" << endl;
  const size_t maxNumEntPerRow = 1;
  RCP<crs_graph_type> G =
    rcp (new crs_graph_type (rowMap, maxNumEntPerRow, Tpetra::StaticProfile));

  // Process 1 gets an implicit zero diagonal entry.
  std::vector<GO> globalIndices (maxNumEntPerRow);
  if (myRank != 1) {
    for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
      const GO gblRow = rowMap->getGlobalElement (lclRow);
      const LO numEnt = 1;
      globalIndices[0] = gblRow;
      G->insertGlobalIndices (gblRow, numEnt, globalIndices.data ());
    }
  }
  G->fillComplete (domMap, ranMap);

  out << "Create test CrsMatrix" << endl;
  RCP<crs_matrix_type> A = rcp (new crs_matrix_type (G));

  std::vector<SC> values (maxNumEntPerRow);
  if (myRank != 1) {
    for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
      const GO gblRow = rowMap->getGlobalElement (lclRow);
      const LO numEnt = 1;
      globalIndices[0] = gblRow;
      values[0] = 1.0;
      A->replaceGlobalValues (gblRow, numEnt, values.data (), globalIndices.data ());
    }
  }
  A->fillComplete (domMap, ranMap);

  RCP<const map_type> colMap = G->getColMap ();
  const LO lclNumCols = static_cast<LO> (colMap->getNodeNumElements ());

  out << "Compute local row norms and diagonal entries" << endl;
  std::vector<val_type> lclRowDiagonalEntries (lclNumRows);
  std::vector<mag_type> lclRowNorms (lclNumRows);
  for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
    const mag_type expectedRowNorm = myRank == 1 ? 0.0 : 1.0;
    lclRowDiagonalEntries[lclRow] = myRank == 1 ? 0.0 : 1.0;
    lclRowNorms[lclRow] = expectedRowNorm;
  }

  // For this matrix, the global row norms are the same as the local
  // row norms, since the matrix's row Map and range Map are the same.
  // This may not be the case for matrices with an overlapping or
  // permuted row Map.
  out << "Compute global row norms" << endl;
  std::vector<mag_type> gblRowNorms (lclRowNorms.begin (), lclRowNorms.end ());

  Teuchos::RCP<crs_matrix_type> A_leftScaled = [&] () {
    Teuchos::RCP<crs_matrix_type> A_copy = deepCopyFillCompleteCrsMatrix (*A);
    A_copy->resumeFill ();

    auto A_lcl = A_copy->getLocalMatrix ();
    auto val_h = Kokkos::create_mirror_view (A_lcl.values);
    Kokkos::deep_copy (val_h, A_lcl.values);
    auto ptr_h = Kokkos::create_mirror_view (A_lcl.graph.row_map);
    Kokkos::deep_copy (ptr_h, A_lcl.graph.row_map);
    auto ind_h = Kokkos::create_mirror_view (A_lcl.graph.entries);
    Kokkos::deep_copy (ind_h, A_lcl.graph.entries);

    for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
      using size_type = typename decltype (A_lcl.graph)::size_type;
      for (size_type k = ptr_h[lclRow]; k < ptr_h[lclRow+1]; ++k) {
        const val_type expectedUnscaledVal = 1.0; // Proc 1 simply has no entry here
        const mag_type rowNorm = gblRowNorms[lclRow];
        const mag_type scalingFactor = assumeSymmetric ?
          KAM::sqrt (rowNorm) : rowNorm;
        val_h[k] = expectedUnscaledVal / scalingFactor;
      }
    }
    Kokkos::deep_copy (A_lcl.values, val_h);
    A_copy->fillComplete (A_copy->getDomainMap (), A_copy->getRangeMap ());
    return A_copy;
  } ();

  if (assumeSymmetric) {
    out << "assumeSymmetric=true: Skip local (column norms, "
      "diagonal entries, and row-scaled column norms)" << endl;
  }
  else {
    out << "assumeSymmetric=false: Compute local (column norms, "
      "diagonal entries, and row-scaled column norms)" << endl;
  }

  std::vector<mag_type> lclColNorms
    (assumeSymmetric ? LO (0) : lclNumCols);
  std::vector<val_type> lclColDiagonalEntries
    (assumeSymmetric ? LO (0) : lclNumCols);
  std::vector<mag_type> lclRowScaledColNorms
    (assumeSymmetric ? LO (0) : lclNumCols);

  if (! assumeSymmetric) {
    // Columns are a little more complicated, since in the usual
    // distributed case, the local column norms may not be the same as
    // the global column norms.
    for (LO lclCol = 0; lclCol < lclNumCols; ++lclCol) {
      lclColNorms[lclCol] = myRank == 1 ? 0.0 : 1.0;
      lclColDiagonalEntries[lclCol] = myRank == 1 ? 0.0 : 1.0;
    } // for each local column

    for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
      const GO gblRow = rowMap->getGlobalElement (lclRow);
      const GO gblCol = gblRow;
      const LO lclCol = colMap->getLocalElement (gblCol);
      if (lclCol != Tpetra::Details::OrdinalTraits<LO>::invalid ()) {
        mag_type lclRowScaledColNorm {0.0};

        const mag_type matrixAbsVal = myRank == 1 ? 0.0 : 1.0;
        lclRowScaledColNorms[lclCol] = myRank == 1 ?
          NaughtyValues<mag_type>::quiet_NaN () : 1.0;
      }
    }
  } // ! assumeSymmetric

  out << "Compute global column norms" << endl;
  std::vector<mag_type> gblColNorms (lclNumCols);
  // The matrix is symmetric, so this holds regardless of assumeSymmetric.
  for (LO lclCol = 0; lclCol < lclNumCols; ++lclCol) {
    gblColNorms[lclCol] = myRank == 1 ? 0.0 : 1.0;
  }

  if (assumeSymmetric) {
    out << "Skip global row-scaled column norms" << endl;
  }
  else {
    out << "Compute global row-scaled column norms" << endl;
  }
  std::vector<mag_type> gblRowScaledColNorms
    (assumeSymmetric ? LO (0) : lclNumCols);
  if (! assumeSymmetric && colMap->getGlobalNumElements () != 0) {
    for (LO lclCol = 0; lclCol < lclNumCols; ++lclCol) {
      gblRowScaledColNorms[lclCol] = myRank == 1 ?
        NaughtyValues<mag_type>::quiet_NaN () : 1.0;
    }
  }

  out << "Compute right-scaled matrix" << endl;
  Teuchos::RCP<crs_matrix_type> A_rightScaled = [&] () {
    Teuchos::RCP<crs_matrix_type> A_copy = deepCopyFillCompleteCrsMatrix (*A);
    A_copy->resumeFill ();

    auto A_lcl = A_copy->getLocalMatrix ();
    auto val_h = Kokkos::create_mirror_view (A_lcl.values);
    Kokkos::deep_copy (val_h, A_lcl.values);
    auto ptr_h = Kokkos::create_mirror_view (A_lcl.graph.row_map);
    Kokkos::deep_copy (ptr_h, A_lcl.graph.row_map);
    auto ind_h = Kokkos::create_mirror_view (A_lcl.graph.entries);
    Kokkos::deep_copy (ind_h, A_lcl.graph.entries);

    for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
      using size_type = typename decltype (A_lcl.graph)::size_type;
      for (size_type k = ptr_h[lclRow]; k < ptr_h[lclRow+1]; ++k) {
        const LO lclCol = ind_h[k];
        const mag_type colNorm = assumeSymmetric ?
          gblColNorms[lclCol] : gblRowScaledColNorms[lclCol];
        // const mag_type rowScalingFactor = assumeSymmetric ?
        //   KAM::sqrt (gblRowNorms[lclRow]) : gblRowNorms[lclRow];
        const mag_type colScalingFactor = assumeSymmetric ?
          KAM::sqrt (colNorm) : colNorm;
        // const mag_type scalingFactor = rowScalingFactor * colScalingFactor;
        const mag_type scalingFactor = colScalingFactor;
        val_h[k] = val_h[k] / scalingFactor;
      }
    }
    Kokkos::deep_copy (A_lcl.values, val_h);
    A_copy->fillComplete (A_copy->getDomainMap (), A_copy->getRangeMap ());
    return A_copy;
  } ();

  const bool lclFoundInf = false;
  const bool lclFoundNan = false;
  const bool lclFoundZeroDiag = myRank == 1 ? true : false;
  const bool lclFoundZeroRowNorm = myRank == 1 ? true : false;
  const bool gblFoundInf = false;
  const bool gblFoundNan = false;
  const bool gblFoundZeroDiag = numProcs == 1 ? false : true;
  const bool gblFoundZeroRowNorm = numProcs == 1 ? false : true;

  return EquilibrationTest<SC, LO, GO, NT> {
    A,
    lclRowNorms,
    lclRowDiagonalEntries,
    lclColNorms,
    lclColDiagonalEntries,
    lclRowScaledColNorms,
    gblRowNorms,
    gblColNorms,
    gblRowScaledColNorms,
    A_leftScaled,
    A_rightScaled,
    lclFoundInf,
    lclFoundNan,
    lclFoundZeroDiag,
    lclFoundZeroRowNorm,
    gblFoundInf,
    gblFoundNan,
    gblFoundZeroDiag,
    gblFoundZeroRowNorm
  };
}

// NOTE (mfh 23 May 2018) This function is only for the test.
template<class SC, class LO, class GO, class NT>
EquilibrationTest<SC, LO, GO, NT>
makeMatrixTestWithExplicitInfAndNan (Teuchos::FancyOStream& out,
                                     const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                                     const bool assumeSymmetric)
{
  using Teuchos::RCP;
  using std::endl;
  using GST = Tpetra::global_size_t;
  using map_type = Tpetra::Map<LO, GO, NT>;
  using crs_graph_type = Tpetra::CrsGraph<LO, GO, NT>;
  using crs_matrix_type = Tpetra::CrsMatrix<SC, LO, GO, NT>;
  using KAT = Kokkos::ArithTraits<SC>;
  using val_type = typename KAT::val_type;
  using mag_type = typename Kokkos::ArithTraits<val_type>::mag_type;
  using KAM = Kokkos::ArithTraits<mag_type>;

  const int myRank = comm->getRank ();
  const int numProcs = comm->getSize ();

  const LO lclNumRows = 1;
  const GO gblNumRows = GO (numProcs) * GO (lclNumRows);
  const GO indexBase = 0;

  out << "Create diagonal matrix problem with an explicit Inf entry on "
    "one process, and an explicit NaN entry on another process" << endl;
  Teuchos::OSTab tab0 (out);

  out << "Create Maps" << endl;
  RCP<const map_type> rowMap =
    rcp (new map_type (GST (gblNumRows), std::size_t (lclNumRows),
                       indexBase, comm));
  RCP<const map_type> domMap = rowMap;
  RCP<const map_type> ranMap = rowMap;

  out << "Create CrsGraph" << endl;
  const size_t maxNumEntPerRow = 1;
  RCP<crs_graph_type> G =
    rcp (new crs_graph_type (rowMap, maxNumEntPerRow, Tpetra::StaticProfile));
  std::vector<GO> globalIndices (maxNumEntPerRow);

  for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
    const GO gblRow = rowMap->getGlobalElement (lclRow);
    const LO numEnt = 1;
    globalIndices[0] = gblRow;
    G->insertGlobalIndices (gblRow, numEnt, globalIndices.data ());
  }
  G->fillComplete (domMap, ranMap);

  out << "Create test CrsMatrix" << endl;
  RCP<crs_matrix_type> A = rcp (new crs_matrix_type (G));
  std::vector<SC> values (maxNumEntPerRow);
  for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
    const GO gblRow = rowMap->getGlobalElement (lclRow);

    const LO numEnt = 1;
    globalIndices[0] = gblRow;
    if (myRank == 1) {
      values[0] = NaughtyValues<SC>::infinity ();
    }
    else if (myRank == 2) {
      values[0] = NaughtyValues<SC>::quiet_NaN ();
    }
    else {
      values[0] = 1.0;
    }
    A->replaceGlobalValues (gblRow, numEnt, values.data (), globalIndices.data ());
  }
  A->fillComplete (domMap, ranMap);

  RCP<const map_type> colMap = G->getColMap ();
  const LO lclNumCols = static_cast<LO> (colMap->getNodeNumElements ());

  out << "Compute local row norms and diagonal entries" << endl;
  std::vector<val_type> lclRowDiagonalEntries (lclNumRows);
  std::vector<mag_type> lclRowNorms (lclNumRows);
  for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
    mag_type expectedRowNorm = 1.0;
    val_type diagEnt = 1.0;
    if (myRank == 1) {
      expectedRowNorm = NaughtyValues<mag_type>::infinity ();
      diagEnt = NaughtyValues<val_type>::infinity ();
    }
    else if (myRank == 2) {
      expectedRowNorm = NaughtyValues<mag_type>::quiet_NaN ();
      diagEnt = NaughtyValues<val_type>::quiet_NaN ();
    }
    lclRowDiagonalEntries[lclRow] = diagEnt;
    lclRowNorms[lclRow] = expectedRowNorm;
  }

  // For this matrix, the global row norms are the same as the local
  // row norms, since the matrix's row Map and range Map are the same.
  // This may not be the case for matrices with an overlapping or
  // permuted row Map.
  out << "Compute global row norms" << endl;
  std::vector<mag_type> gblRowNorms (lclRowNorms.begin (), lclRowNorms.end ());

  Teuchos::RCP<crs_matrix_type> A_leftScaled = [&] () {
    Teuchos::RCP<crs_matrix_type> A_copy = deepCopyFillCompleteCrsMatrix (*A);
    A_copy->resumeFill ();

    auto A_lcl = A_copy->getLocalMatrix ();
    auto val_h = Kokkos::create_mirror_view (A_lcl.values);
    Kokkos::deep_copy (val_h, A_lcl.values);
    auto ptr_h = Kokkos::create_mirror_view (A_lcl.graph.row_map);
    Kokkos::deep_copy (ptr_h, A_lcl.graph.row_map);
    auto ind_h = Kokkos::create_mirror_view (A_lcl.graph.entries);
    Kokkos::deep_copy (ind_h, A_lcl.graph.entries);

    for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
      using size_type = typename decltype (A_lcl.graph)::size_type;
      for (size_type k = ptr_h[lclRow]; k < ptr_h[lclRow+1]; ++k) {
        val_type expectedUnscaledVal = 1.0;
        if (myRank == 1) {
          expectedUnscaledVal = NaughtyValues<val_type>::infinity ();
        }
        else if (myRank == 2) {
          expectedUnscaledVal = NaughtyValues<val_type>::quiet_NaN ();
        }
        const mag_type rowNorm = gblRowNorms[lclRow];
        const mag_type scalingFactor = assumeSymmetric ?
          KAM::sqrt (rowNorm) : rowNorm;
        val_h[k] = expectedUnscaledVal / scalingFactor;
      }
    }
    Kokkos::deep_copy (A_lcl.values, val_h);
    A_copy->fillComplete (A_copy->getDomainMap (), A_copy->getRangeMap ());
    return A_copy;
  } ();

  if (assumeSymmetric) {
    out << "assumeSymmetric=true: Skip local (column norms, "
      "diagonal entries, and row-scaled column norms)" << endl;
  }
  else {
    out << "assumeSymmetric=false: Compute local (column norms, "
      "diagonal entries, and row-scaled column norms)" << endl;
  }

  std::vector<mag_type> lclColNorms
    (assumeSymmetric ? LO (0) : lclNumCols);
  std::vector<val_type> lclColDiagonalEntries
    (assumeSymmetric ? LO (0) : lclNumCols);
  std::vector<mag_type> lclRowScaledColNorms
    (assumeSymmetric ? LO (0) : lclNumCols);

  if (! assumeSymmetric) {
    // Columns are a little more complicated, since in the usual
    // distributed case, the local column norms may not be the same as
    // the global column norms.
    for (LO lclCol = 0; lclCol < lclNumCols; ++lclCol) {
      mag_type expectedNorm = 1.0;
      val_type diagEnt = 1.0;
      if (myRank == 1) {
        expectedNorm = NaughtyValues<mag_type>::infinity ();
        diagEnt = NaughtyValues<val_type>::infinity ();
      }
      else if (myRank == 2) {
        expectedNorm = NaughtyValues<mag_type>::quiet_NaN ();
        diagEnt = NaughtyValues<val_type>::quiet_NaN ();
      }
      lclColNorms[lclCol] = expectedNorm;
      lclColDiagonalEntries[lclCol] = diagEnt;
    } // for each local column

    for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
      const GO gblRow = rowMap->getGlobalElement (lclRow);
      const GO gblCol = gblRow;
      const LO lclCol = colMap->getLocalElement (gblCol);
      mag_type lclRowScaledColNorm {0.0};

      mag_type matrixAbsVal = 1.0;
      if (myRank == 1) {
        matrixAbsVal = NaughtyValues<mag_type>::infinity ();
      }
      else if (myRank == 2) {
        matrixAbsVal = NaughtyValues<mag_type>::quiet_NaN ();
      }
      lclRowScaledColNorms[lclCol] = matrixAbsVal / gblRowNorms[lclRow];
    }
  } // ! assumeSymmetric

  out << "Compute global column norms" << endl;
  std::vector<mag_type> gblColNorms (lclNumCols);
  // The matrix is symmetric, so this holds regardless of assumeSymmetric.
  for (LO lclCol = 0; lclCol < lclNumCols; ++lclCol) {
    mag_type expectedNorm = 1.0;
    if (myRank == 1) {
      expectedNorm = NaughtyValues<mag_type>::infinity ();
    }
    else if (myRank == 2) {
      expectedNorm = NaughtyValues<mag_type>::quiet_NaN ();
    }
    gblColNorms[lclCol] = expectedNorm;
  }

  if (assumeSymmetric) {
    out << "Skip global row-scaled column norms" << endl;
  }
  else {
    out << "Compute global row-scaled column norms" << endl;
  }
  std::vector<mag_type> gblRowScaledColNorms
    (assumeSymmetric ? LO (0) : lclNumCols);
  if (! assumeSymmetric && colMap->getGlobalNumElements () != 0) {
    for (LO lclCol = 0; lclCol < lclNumCols; ++lclCol) {
      mag_type expectedNorm = 1.0;
      if (myRank == 1) {
        expectedNorm = NaughtyValues<mag_type>::quiet_NaN ();
      }
      else if (myRank == 2) {
        expectedNorm = NaughtyValues<mag_type>::quiet_NaN ();
      }
      gblRowScaledColNorms[lclCol] = expectedNorm;
    }
  }

  Teuchos::RCP<crs_matrix_type> A_rightScaled = [&] () {
    Teuchos::RCP<crs_matrix_type> A_copy = deepCopyFillCompleteCrsMatrix (*A);
    A_copy->resumeFill ();

    auto A_lcl = A_copy->getLocalMatrix ();
    auto val_h = Kokkos::create_mirror_view (A_lcl.values);
    Kokkos::deep_copy (val_h, A_lcl.values);
    auto ptr_h = Kokkos::create_mirror_view (A_lcl.graph.row_map);
    Kokkos::deep_copy (ptr_h, A_lcl.graph.row_map);
    auto ind_h = Kokkos::create_mirror_view (A_lcl.graph.entries);
    Kokkos::deep_copy (ind_h, A_lcl.graph.entries);

    for (LO lclRow = 0; lclRow < lclNumRows; ++lclRow) {
      using size_type = typename decltype (A_lcl.graph)::size_type;
      for (size_type k = ptr_h[lclRow]; k < ptr_h[lclRow+1]; ++k) {
        const LO lclCol = ind_h[k];
        const mag_type colNorm = assumeSymmetric ?
          gblColNorms[lclCol] : gblRowScaledColNorms[lclCol];
        // const mag_type rowScalingFactor = assumeSymmetric ?
        //   KAM::sqrt (gblRowNorms[lclRow]) : gblRowNorms[lclRow];
        const mag_type colScalingFactor = assumeSymmetric ?
          KAM::sqrt (colNorm) : colNorm;
        // const mag_type scalingFactor = rowScalingFactor * colScalingFactor;
        const mag_type scalingFactor = colScalingFactor;
        val_h[k] = val_h[k] / scalingFactor;
      }
    }
    Kokkos::deep_copy (A_lcl.values, val_h);
    A_copy->fillComplete (A_copy->getDomainMap (), A_copy->getRangeMap ());
    return A_copy;
  } ();

  const bool lclFoundInf = myRank == 1 ? true : false;
  const bool lclFoundNan = myRank == 2 ? true : false;
  const bool lclFoundZeroDiag = false;
  const bool lclFoundZeroRowNorm = false;
  const bool gblFoundInf = numProcs > 1;
  const bool gblFoundNan = numProcs > 2;
  const bool gblFoundZeroDiag = false;
  const bool gblFoundZeroRowNorm = false;

  return EquilibrationTest<SC, LO, GO, NT> {
    A,
    lclRowNorms,
    lclRowDiagonalEntries,
    lclColNorms,
    lclColDiagonalEntries,
    lclRowScaledColNorms,
    gblRowNorms,
    gblColNorms,
    gblRowScaledColNorms,
    A_leftScaled,
    A_rightScaled,
    lclFoundInf,
    lclFoundNan,
    lclFoundZeroDiag,
    lclFoundZeroRowNorm,
    gblFoundInf,
    gblFoundNan,
    gblFoundZeroDiag,
    gblFoundZeroRowNorm
  };
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Equilibration, Test0, SC, LO, GO, NT)
{
  // We are now in a class method declared by the above macro.
  // The method has these input arguments:
  // (Teuchos::FancyOStream& out, bool& success)
  using Teuchos::Comm;
  using Teuchos::FancyOStream;
  using Teuchos::outArg;
  using Teuchos::RCP;
  using Teuchos::rcpFromRef;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  using std::endl;

  //const bool debugMode = true;
  const bool debugMode = false;
  RCP<FancyOStream> errStrmPtr = debugMode ?
    Teuchos::getFancyOStream (rcpFromRef (std::cerr)) : Teuchos::null;
  FancyOStream& testOut = debugMode ? *errStrmPtr : out;

  testOut << "Test equilibration" << endl;
  Teuchos::OSTab tab0 (testOut);
  auto comm = Tpetra::getDefaultComm ();

  using test_return_type = EquilibrationTest<SC, LO, GO, NT>;
  using test_type =
    std::function<test_return_type (FancyOStream&,
                                    const RCP<const Comm<int> >&,
                                    const bool)>;
  std::array<test_type, 4> tests {
    makeSymmetricPositiveDefiniteTridiagonalMatrixTest<SC, LO, GO, NT>,
    makeMatrixTestWithExplicitZeroDiag<SC, LO, GO, NT>,
    makeMatrixTestWithExplicitInfAndNan<SC, LO, GO, NT>,
    makeMatrixTestWithImplicitZeroDiag<SC, LO, GO, NT>
  };

  for (bool assumeSymmetric : {false, true}) {
    for (auto&& currentTest : tests) {
      // The default FancyOStream 'out' only prints to Process 0 by
      // default.  Thus, we gather up output from each process into a
      // single string, and print it on Process 0 at the end.
      RCP<std::ostringstream> osPtr (new std::ostringstream);
      RCP<FancyOStream> curOutPtr = Teuchos::getFancyOStream (osPtr);
      FancyOStream& curOut = *curOutPtr;

      curOut << endl << endl
             << ">>> Process " << comm->getRank () << ":" << endl << endl
             << "assumeSymmetric=" << (assumeSymmetric ? "true" : "false")
             << endl << endl;
      Teuchos::OSTab tab1 (curOut);
      auto test0 = currentTest (curOut, comm, assumeSymmetric);

      bool curSuccess = true;
      testEquilibration (curOut, curSuccess, test0, assumeSymmetric);
      success = success && curSuccess;

      int lclSuccess = success ? 1 : 0;
      int gblSuccess = 0; // output argument
      reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
      TEST_EQUALITY( gblSuccess, 1 );
      if (gblSuccess != 1) {
        // Don't actually print the test output unless the test failed.
        Tpetra::Details::gathervPrint (testOut, osPtr->str (), *comm);
        out << "Test FAILED on some process!" << endl;
        return;
      }
    }
  }
}

// Define the set of unit tests to instantiate in this file.
#define UNIT_TEST_GROUP( SC, LO, GO, NT ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( Equilibration, Test0, SC, LO, GO, NT )

#include "TpetraCore_ETIHelperMacros.h"

TPETRA_ETI_MANGLING_TYPEDEFS()

TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR( UNIT_TEST_GROUP )

} // namespace (anonymous)

