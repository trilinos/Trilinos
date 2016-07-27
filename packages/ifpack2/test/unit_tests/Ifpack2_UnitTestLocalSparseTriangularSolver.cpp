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

/// \file Ifpack2_UnitTestLocalSparseTriangularSolver.cpp
/// \brief Unit tests for Ifpack2::LocalSparseTriangularSolver

#include "Teuchos_UnitTestHarness.hpp"
#include "Ifpack2_LocalSparseTriangularSolver.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_DefaultPlatform.hpp"

namespace { // (anonymous)

using Teuchos::outArg;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::REDUCE_MIN;
using Teuchos::reduceAll;
using std::endl;
typedef Tpetra::global_size_t GST;

template<class CrsMatrixType, class MultiVectorType>
void
referenceApply (const CrsMatrixType& A,
                const MultiVectorType& X,
                MultiVectorType& Y,
                const Teuchos::ETransp mode,
                const typename MultiVectorType::scalar_type& alpha,
                const typename MultiVectorType::scalar_type& beta)
{
  typedef typename MultiVectorType::scalar_type ST;
  typedef Teuchos::ScalarTraits<ST> STS;
  typedef MultiVectorType MV;

  if (beta == STS::zero ()) {
    if (alpha == STS::zero ()) {
      Y.putScalar (STS::zero ()); // Y := 0 * Y (ignore contents of Y)
    }
    else { // alpha != 0
      A.template localSolve<ST, ST> (X, Y, mode);
      if (alpha != STS::one ()) {
        Y.scale (alpha);
      }
    }
  }
  else { // beta != 0
    if (alpha == STS::zero ()) {
      Y.scale (beta); // Y := beta * Y
    }
    else { // alpha != 0
      MV Y_tmp (Y, Teuchos::Copy);
      A.template localSolve<ST, ST> (X, Y_tmp, mode); // Y_tmp := M * X
      Y.update (alpha, Y_tmp, beta); // Y := beta * Y + alpha * Y_tmp
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(LocalSparseTriangularSolver, CompareToLocalSolve, Scalar, LocalOrdinal, GlobalOrdinal)
{
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;
  typedef Tpetra::Map<LO, GO> map_type;
  typedef typename map_type::device_type device_type;
  typedef Tpetra::CrsMatrix<Scalar, LO, GO> crs_matrix_type;
  typedef Tpetra::RowMatrix<Scalar, LO, GO> row_matrix_type;
  typedef Tpetra::MultiVector<Scalar, LO, GO> mv_type;
  typedef Ifpack2::LocalSparseTriangularSolver<row_matrix_type> solver_type;
  typedef Teuchos::ScalarTraits<Scalar> STS;
  typedef typename crs_matrix_type::mag_type mag_type;
  int lclSuccess = 1;
  int gblSuccess = 1;

  out << "Ifpack2::LocalSparseTriangularSolver CompareToLocalSolve" << endl;
  Teuchos::OSTab tab0 (out);

  auto comm = Tpetra::DefaultPlatform::getDefaultPlatform ().getComm ();

  const LO lclNumRows = 5;
  const GO gblNumRows = comm->getSize () * lclNumRows;
  const GO indexBase = 0;
  RCP<const map_type> rowMap =
    rcp (new map_type (static_cast<GST> (gblNumRows),
                       static_cast<size_t> (lclNumRows),
                       indexBase, comm));

  // If we construct an upper or lower triangular matrix with an
  // implicit unit diagonal, then we need to specify the column Map
  // explicitly.  Otherwise, the matrix will report having the wrong
  // number of columns.
  RCP<const map_type> colMap = rowMap;

  const GO maxLowerBandwidth = 1;
  const GO maxUpperBandwidth = 1;
  const GO maxNumEntPerRow =
    maxLowerBandwidth + static_cast<GO> (1) + maxUpperBandwidth;

  Teuchos::Array<GO> gblColIndsBuf (maxNumEntPerRow);
  Teuchos::Array<Scalar> valsBuf (maxNumEntPerRow);

  const bool isUpperTriangularValues[] = {false, true};
  const bool implicitUnitDiagValues[] = {false, true};

  // Make sure the matrix is diagonally dominant.  Cast to mag_type
  // first, since some Scalar types (like std::complex<T> for T =
  // float, double) don't have a direct conversion from integer types.
  const Scalar diagVal = static_cast<Scalar> (static_cast<mag_type> (2 * maxNumEntPerRow));
  const Scalar offDiagVal = STS::one () / diagVal;

  // For this test, the opposite of "is upper triangular" is "is lower
  // triangular."
  for (bool isUpperTriangular : isUpperTriangularValues) {
    out << (isUpperTriangular ? "Upper" : "Lower") << " triangular" << endl;
    Teuchos::OSTab tab1 (out);

    for (bool implicitUnitDiag : implicitUnitDiagValues) {
      out << (implicitUnitDiag ? "Implicit" : "Explicit") << " unit diagonal" << endl;
      Teuchos::OSTab tab2 (out);

      out << "Construct the test matrix A" << endl;
      RCP<crs_matrix_type> A =
        rcp (new crs_matrix_type (rowMap, colMap,
                                  static_cast<size_t> (maxNumEntPerRow),
                                  Tpetra::StaticProfile));

      out << "Fill the test matrix A" << endl;
      // Pick the entries of the matrix so that triangular solves are
      // guaranteed not to get Inf or NaN results.

      const LO minLclRow = 0;
      const LO maxLclRow = lclNumRows;

      for (LO lclRow = minLclRow; lclRow < maxLclRow; ++lclRow) {
        const GO gblRow = rowMap->getGlobalElement (lclRow);
        TEST_ASSERT( gblRow != Tpetra::Details::OrdinalTraits<GO>::invalid () );
        if (! success) {
          break;
        }

        // Don't subtract at all when computing these bounds (even the
        // upper one), in case LO is unsigned.
        //
        //if (lclRow - maxLowerBandwidth < minLclRow) { ... }
        // if (lclRow + maxUpperBandwidth >= maxLclRow) { ... }

        LO lowerBandwidth = (lclRow < minLclRow + maxLowerBandwidth) ?
          lclRow :
          maxLowerBandwidth;
        LO upperBandwidth = (lclRow + maxUpperBandwidth >= maxLclRow) ?
          (maxLclRow - static_cast<LO> (1)) - lclRow :
          maxUpperBandwidth;

        TEST_ASSERT( lclRow >= lowerBandwidth + minLclRow );
        TEST_ASSERT( lclRow + upperBandwidth < maxLclRow );
        if (! success) {
          // If the test messed up the bounds, don't continue; go fix
          // the test.  We don't want to crash with an illegal array
          // access.
          break;
        }

        LO numEnt = 0;
        if (! implicitUnitDiag) {
          numEnt++;
        }
        if (isUpperTriangular) {
          numEnt += upperBandwidth;
        }
        else { // is lower triangular
          numEnt += lowerBandwidth;
        }

        TEST_ASSERT( numEnt <= maxNumEntPerRow );
        if (! success) {
          break;
        }

        Teuchos::ArrayView<GO> gblColInds = gblColIndsBuf.view (0, numEnt);
        Teuchos::ArrayView<Scalar> vals = valsBuf.view (0, numEnt);
        TEST_EQUALITY( static_cast<LO> (gblColInds.size ()), numEnt );
        TEST_EQUALITY( static_cast<LO> (vals.size ()), numEnt );
        if (! success) {
          break;
        }

        LO curPos = 0;
        if (! implicitUnitDiag) {
          gblColInds[curPos] = gblRow;
          vals[curPos] = diagVal;
          curPos++;
        }
        if (isUpperTriangular) {
          for (LO k = 0; k < upperBandwidth; ++k) {
            gblColInds[curPos] = gblRow + static_cast<GO> (k+1);
            vals[curPos] = offDiagVal;
            curPos++;
          }
        }
        else { // is lower triangular
          for (LO k = 0; k < lowerBandwidth; ++k) {
            gblColInds[curPos] = gblRow - static_cast<GO> (k+1);
            vals[curPos] = offDiagVal;
            curPos++;
          }
        }

        TEST_NOTHROW( A->insertGlobalValues (gblRow, gblColInds, vals) );
      }

      lclSuccess = success ? 1 : 0;
      gblSuccess = 0; // to be revised
      reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
      TEST_EQUALITY( gblSuccess, 1 );
      if (! gblSuccess) {
        out << "Aborting test" << endl;
        return;
      }

      out << "Call fillComplete on the test matrix A" << endl;

      TEST_NOTHROW( A->fillComplete () );
      lclSuccess = success ? 1 : 0;
      gblSuccess = 0; // to be revised
      reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
      TEST_EQUALITY( gblSuccess, 1 );
      if (! gblSuccess) {
        out << "Aborting test" << endl;
        return;
      }

      out << "Make a deep copy of A" << endl;
      // Make a true deep copy of A.  This will help us test whether
      // LocalSparseTriangularSolver modified the entries of the matrix
      // (it shouldn't).
      RCP<crs_matrix_type> A_copy;
      {
        typedef typename crs_matrix_type::local_matrix_type local_matrix_type;
        typedef typename crs_matrix_type::local_graph_type local_graph_type;

        local_matrix_type A_lcl = A->getLocalMatrix ();

        typename local_matrix_type::row_map_type::non_const_type ptr ("A_copy.ptr", A_lcl.graph.row_map.dimension_0 ());
        Kokkos::deep_copy (ptr, A_lcl.graph.row_map);
        typename local_graph_type::entries_type::non_const_type ind ("A_copy.ind", A_lcl.graph.entries.dimension_0 ());
        Kokkos::deep_copy (ind, A_lcl.graph.entries);
        typename local_matrix_type::values_type::non_const_type val ("A_copy.val", A_lcl.values.dimension_0 ());
        Kokkos::deep_copy (val, A_lcl.values);

        TEST_NOTHROW( A_copy = rcp (new crs_matrix_type (rowMap, A->getColMap (), ptr, ind, val)) );
      }
      TEST_ASSERT( ! A_copy.is_null () );
      lclSuccess = success ? 1 : 0;
      gblSuccess = 0; // to be revised
      reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
      TEST_EQUALITY( gblSuccess, 1 );
      if (! gblSuccess) {
        out << "Aborting test" << endl;
        return;
      }

      out << "Call fillComplete on the deep copy of A" << endl;
      TEST_NOTHROW( A_copy->fillComplete (A->getDomainMap (), A->getRangeMap ()) );
      lclSuccess = success ? 1 : 0;
      gblSuccess = 0; // to be revised
      reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
      TEST_EQUALITY( gblSuccess, 1 );
      if (! gblSuccess) {
        out << "Aborting test" << endl;
        return;
      }

      out << "Create the solver" << endl;
      RCP<solver_type> solver;
      TEST_NOTHROW( solver = rcp (new solver_type (A)) );
      TEST_ASSERT( ! solver.is_null () );
      lclSuccess = success ? 1 : 0;
      gblSuccess = 0; // to be revised
      reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
      TEST_EQUALITY( gblSuccess, 1 );
      if (! gblSuccess) {
        out << "Aborting test" << endl;
        return;
      }

      out << "Set up the solver" << endl;
      // Set up the solver.  For now, we only do this once.  Another test
      // should exercise the case where we do repeated solves, between
      // which we change the matrix itself, its graph, or its values.
      TEST_NOTHROW( solver->initialize () );
      if (success) {
        TEST_NOTHROW( solver->compute () );
      }
      lclSuccess = success ? 1 : 0;
      gblSuccess = 0; // to be revised
      reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
      TEST_EQUALITY( gblSuccess, 1 );
      if (! gblSuccess) {
        out << "Aborting test" << endl;
        return;
      }

      out << "Test the solver" << endl;

      // For different combinations of alpha, beta, and mode, apply
      // the solver.  Test against Tpetra::CrsMatrix::localSolve with
      // A_copy, a deep copy of the original matrix A.  This tests
      // whether the solver changes A (it shouldn't).

      const size_t numVecValues[] = {static_cast<size_t> (1), static_cast<size_t> (3)};
      const Scalar ZERO = STS::zero ();
      const Scalar ONE = STS::one ();
      const Scalar THREE = ONE + ONE + ONE;
      const Scalar alphaValues[] = {ZERO, ONE, -ONE, THREE, -THREE};
      const Scalar betaValues[] = {ZERO, ONE, -ONE, THREE, -THREE};
      const Teuchos::ETransp modes[] = {Teuchos::NO_TRANS, Teuchos::TRANS, Teuchos::CONJ_TRANS};

      for (size_t numVecs : numVecValues) {
        out << "numVecs: " << numVecs << endl;
        Teuchos::OSTab tab3 (out);

        typename Kokkos::View<mag_type*, device_type>::HostMirror norms ("norms", numVecs);

        for (Teuchos::ETransp mode : modes) {
          out << "mode: ";
          if (mode == Teuchos::NO_TRANS) {
            out << "NO_TRANS";
          }
          else if (mode == Teuchos::TRANS) {
            out << "TRANS";
          }
          else if (mode == Teuchos::CONJ_TRANS) {
            out << "CONJ_TRANS";
          }
          else {
            TEST_ASSERT( false );
          }
          out << endl;
          Teuchos::OSTab tab4 (out);

          auto domainMap = mode == Teuchos::NO_TRANS ?
            A_copy->getDomainMap () :
            A_copy->getRangeMap ();
          auto rangeMap = mode == Teuchos::NO_TRANS ?
            A_copy->getRangeMap () :
            A_copy->getDomainMap ();
          mv_type X (domainMap, numVecs);
          mv_type Y (rangeMap, numVecs);
          mv_type X_copy (X.getMap (), numVecs);
          mv_type Y_copy (Y.getMap (), numVecs);

          for (Scalar alpha : alphaValues) {
            out << "alpha: " << alpha << endl;
            Teuchos::OSTab tab5 (out);

            for (Scalar beta : betaValues) {
              out << "beta: " << beta << endl;
              Teuchos::OSTab tab6 (out);

              X.randomize ();
              Y.randomize ();
              Tpetra::deep_copy (X_copy, X);
              Tpetra::deep_copy (Y_copy, Y);

              TEST_NOTHROW( solver->apply (X, Y, mode, alpha, beta) );

              // Don't continue unless no processes threw.  Otherwise, the
              // test may deadlock if run with multiple MPI processes,
              // since comparing the vectors requires all-reduces.
              lclSuccess = success ? 1 : 0;
              gblSuccess = 0; // to be revised
              reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
              TEST_EQUALITY( gblSuccess, 1 );
              if (! gblSuccess) {
                out << "Aborting test" << endl;
                return;
              }

              // Test against a reference implementation.
              TEST_NOTHROW( referenceApply (*A_copy, X_copy, Y_copy, mode, alpha, beta) );

              // Don't continue unless no processes threw.  Otherwise, the
              // test may deadlock if run with multiple MPI processes,
              // since comparing the vectors requires all-reduces.
              lclSuccess = success ? 1 : 0;
              gblSuccess = 0; // to be revised
              reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
              TEST_EQUALITY( gblSuccess, 1 );
              if (! gblSuccess) {
                out << "Aborting test" << endl;
                return;
              }

              // Compare Y and Y_copy.  Compute difference in Y_copy.
              Y_copy.update (ONE, Y, -ONE); // Y_copy := Y - Y_copy
              Y_copy.normInf (norms);

              const mag_type tol = static_cast<mag_type> (gblNumRows) * STS::eps ();
              for (size_t j = 0; j < numVecs; ++j) {
                TEST_ASSERT( norms(j) <= tol );
                if (norms(j) > tol) {
                  out << "norms(" << j << ") = " << norms(j) << " > tol = " << tol << endl;
                }
              }
            }
          }
        }
      }
    }
  }

  // Make sure that all processes succeeded.
  lclSuccess = success ? 1 : 0;
  gblSuccess = 0; // to be revised
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  TEST_EQUALITY( gblSuccess, 1 );
  if (gblSuccess != 1) {
    out << "The test FAILED on at least one MPI process." << endl;
  }
}

#define UNIT_TEST_GROUP_SC_LO_GO(SC, LO, GO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(LocalSparseTriangularSolver, CompareToLocalSolve, SC, LO, GO)

#include "Ifpack2_ETIHelperMacros.h"

IFPACK2_ETI_MANGLING_TYPEDEFS()

// Test all enabled combinations of Scalar (SC), LocalOrdinal (LO),
// and GlobalOrdinal (GO) types, where Scalar is real.

IFPACK2_INSTANTIATE_SLG_REAL( UNIT_TEST_GROUP_SC_LO_GO )

} // namespace (anonymous)

