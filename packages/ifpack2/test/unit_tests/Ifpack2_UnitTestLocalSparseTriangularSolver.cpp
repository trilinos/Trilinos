// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// \file Ifpack2_UnitTestLocalSparseTriangularSolver.cpp
/// \brief Unit tests for Ifpack2::LocalSparseTriangularSolver

#include "Teuchos_UnitTestHarness.hpp"
#include "Ifpack2_LocalSparseTriangularSolver.hpp"
#include "Tpetra_Details_gathervPrint.hpp"
#include "Tpetra_BlockView.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_RowMatrixTransposer.hpp"
#include "KokkosSparse_trsv.hpp"

namespace { // (anonymous)

using Teuchos::outArg;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::REDUCE_MIN;
using Teuchos::reduceAll;
using std::endl;
typedef Tpetra::global_size_t GST;

template<class CrsMatrixType>
void
localSolve (Tpetra::MultiVector<
              typename CrsMatrixType::scalar_type,
              typename CrsMatrixType::local_ordinal_type,
              typename CrsMatrixType::global_ordinal_type,
              typename CrsMatrixType::node_type>& X,
            const CrsMatrixType& A,
            const Tpetra::MultiVector<
              typename CrsMatrixType::scalar_type,
              typename CrsMatrixType::local_ordinal_type,
              typename CrsMatrixType::global_ordinal_type,
              typename CrsMatrixType::node_type>& Y,
            const bool isUpperTriangular, // opposite is "lower triangular"
            const bool implicitUnitDiag, // opposite is "explicitly stored, possibly non-unit diagonal"
            Teuchos::ETransp mode)
{
  using Teuchos::CONJ_TRANS;
  using Teuchos::NO_TRANS;
  using Teuchos::TRANS;
  using scalar_type = typename CrsMatrixType::scalar_type;
  using STS = Teuchos::ScalarTraits<scalar_type>;
  const char prefix[] = "localSolve: ";

  TEUCHOS_TEST_FOR_EXCEPTION
    (! A.isFillComplete (), std::runtime_error,
     prefix << "The matrix is not fill complete.");
  TEUCHOS_TEST_FOR_EXCEPTION
    (! X.isConstantStride () || ! Y.isConstantStride (), std::invalid_argument,
     prefix << "X and Y must be constant stride.");
  TEUCHOS_TEST_FOR_EXCEPTION
    (STS::isComplex && mode == TRANS, std::logic_error, prefix << "This "
     "function does not currently support non-conjugated transposed solve "
     "(mode == Teuchos::TRANS) for complex scalar types.");

  // FIXME (mfh 27 Aug 2014) Tpetra has always made the odd decision
  // that if _some_ diagonal entries are missing locally, then it
  // assumes that the matrix has an implicitly stored unit diagonal.
  // Whether the matrix has an implicit unit diagonal or not should
  // be up to the user to decide.  What if the graph has no diagonal
  // entries, and the user wants it that way?  The only reason this
  // matters, though, is for the triangular solve, and in that case,
  // missing diagonal entries will cause trouble anyway.  However,
  // it would make sense to warn the user if they ask for a
  // triangular solve with an incomplete diagonal.  Furthermore,
  // this code should only assume an implicitly stored unit diagonal
  // if the matrix has _no_ explicitly stored diagonal entries.

  const std::string uplo = isUpperTriangular ? "U" : "L";
  const std::string trans = (mode == Teuchos::CONJ_TRANS) ? "C" :
    (mode == Teuchos::TRANS ? "T" : "N");
  const std::string diag = implicitUnitDiag ? "U" : "N";

  auto A_lcl = A.getLocalMatrixHost ();

  if (X.isConstantStride () && Y.isConstantStride ()) {
    auto X_lcl = X.getLocalViewHost (Tpetra::Access::OverwriteAll);
    auto Y_lcl = Y.getLocalViewHost (Tpetra::Access::ReadOnly);
    KokkosSparse::trsv (uplo.c_str (), trans.c_str (), diag.c_str (),
                        A_lcl, Y_lcl, X_lcl);
  }
  else {
    const size_t numVecs =
      std::min (X.getNumVectors (), Y.getNumVectors ());
    for (size_t j = 0; j < numVecs; ++j) {
      auto X_j = X.getVectorNonConst (j);
      auto Y_j = Y.getVector (j);
      auto X_lcl = X_j->getLocalViewHost (Tpetra::Access::OverwriteAll);
      auto Y_lcl = Y_j->getLocalViewHost (Tpetra::Access::ReadOnly);
      KokkosSparse::trsv (uplo.c_str (), trans.c_str (),
                          diag.c_str (), A_lcl, Y_lcl, X_lcl);
    }
  }
}

template<class CrsMatrixType, class MultiVectorType>
void
referenceApply (MultiVectorType& Y,
                const CrsMatrixType& A,
                const MultiVectorType& X,
                const bool isUpperTriangular, // opposite is "lower triangular"
                const bool implicitUnitDiag, // opposite is "explicitly stored, possibly non-unit diagonal"
                const Teuchos::ETransp mode,
                const typename MultiVectorType::scalar_type& alpha,
                const typename MultiVectorType::scalar_type& beta)
{
  typedef typename MultiVectorType::scalar_type ST;
  typedef Teuchos::ScalarTraits<ST> STS;
  typedef MultiVectorType MV;

  if (beta == STS::zero ()) {
    Y.putScalar (STS::zero ()); // Y := 0 * Y (ignore contents of Y)
    if (alpha != STS::zero ()) {
      localSolve (Y, A, X, isUpperTriangular, implicitUnitDiag, mode);
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
      localSolve (Y_tmp, A, X, isUpperTriangular, implicitUnitDiag, mode); // Y_tmp := M * X
      Y.update (alpha, Y_tmp, beta); // Y := beta * Y + alpha * Y_tmp
    }
  }
}

struct TrisolverDetails {
  enum Enum { Internal, HTS, KSPTRSV };
};

static bool isGblSuccess (const bool success, Teuchos::FancyOStream& out)
{
  auto comm = Tpetra::getDefaultComm ();
  const int lclSuccess = success ? 1 : 0;
  int gblSuccess = 0; // to be revised
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  if (! gblSuccess) {
    out << "Aborting test" << endl;
    return false;
  }
  return true;
}

template<typename Scalar, typename LocalOrdinal, typename GlobalOrdinal>
void testCompareToLocalSolve (bool& success, Teuchos::FancyOStream& out,
                              const TrisolverDetails::Enum trisolverType)
{
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;
  typedef Tpetra::Map<LO, GO> map_type;
  typedef typename map_type::device_type device_type;
  typedef Tpetra::CrsGraph<LO, GO> crs_graph_type;
  typedef Tpetra::CrsMatrix<Scalar, LO, GO> crs_matrix_type;
  typedef Tpetra::RowMatrix<Scalar, LO, GO> row_matrix_type;
  typedef Tpetra::MultiVector<Scalar, LO, GO> mv_type;
  typedef Ifpack2::LocalSparseTriangularSolver<row_matrix_type> solver_type;
  typedef Teuchos::ScalarTraits<Scalar> STS;
  typedef typename crs_matrix_type::mag_type mag_type;

  out << "Ifpack2::LocalSparseTriangularSolver CompareToLocalSolve" << endl;
  Teuchos::OSTab tab0 (out);

  auto comm = Tpetra::getDefaultComm ();

  const LO lclNumRows = 5;
  const GO gblNumRows = comm->getSize () * lclNumRows;
  const GO indexBase = 0;
  RCP<const map_type> rowMap =
    rcp (new map_type (static_cast<GST> (gblNumRows),
                       static_cast<std::size_t> (lclNumRows),
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
                                  static_cast<std::size_t> (maxNumEntPerRow)));

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

      if (! isGblSuccess (success, out)) return;

      out << "Call fillComplete on the test matrix A" << endl;

      TEST_NOTHROW( A->fillComplete () );
      if (! isGblSuccess (success, out)) return;

      out << "Make a deep copy of A" << endl;
      // Make a true deep copy of A.  This will help us test whether
      // LocalSparseTriangularSolver modified the entries of the matrix
      // (it shouldn't).
      RCP<crs_matrix_type> A_copy;
      {
        typedef typename crs_matrix_type::local_matrix_device_type local_matrix_type;
        typedef typename crs_graph_type::local_graph_device_type local_graph_type;

        local_matrix_type A_lcl = A->getLocalMatrixDevice ();

        typename local_matrix_type::row_map_type::non_const_type ptr ("A_copy.ptr", A_lcl.graph.row_map.extent (0));
        Kokkos::deep_copy (ptr, A_lcl.graph.row_map);
        typename local_graph_type::entries_type::non_const_type ind ("A_copy.ind", A_lcl.graph.entries.extent (0));
        Kokkos::deep_copy (ind, A_lcl.graph.entries);
        typename local_matrix_type::values_type::non_const_type val ("A_copy.val", A_lcl.values.extent (0));
        Kokkos::deep_copy (val, A_lcl.values);

        TEST_NOTHROW( A_copy = rcp (new crs_matrix_type (rowMap, A->getColMap (), ptr, ind, val)) );
      }
      TEST_ASSERT( ! A_copy.is_null () );
      if (! isGblSuccess (success, out)) return;

      out << "Call fillComplete on the deep copy of A" << endl;
      TEST_NOTHROW( A_copy->fillComplete (A->getDomainMap (), A->getRangeMap ()) );
      if (! isGblSuccess (success, out)) return;

      out << "Create the solver" << endl;
      RCP<solver_type> solver;
      TEST_NOTHROW( solver = rcp (new solver_type (A)) );
      TEST_ASSERT( ! solver.is_null () );
      if (! isGblSuccess (success, out)) return;

      if (trisolverType == TrisolverDetails::HTS) {
        out << "Set solver parameters" << endl;
        // This line can throw. It should throw if HTS is not built in.
        Teuchos::ParameterList pl ("LocalSparseTriangularSolver parameters");
        pl.set ("trisolver: type", "HTS");
        try {
          solver->setParameters (pl);
        } catch (...) {
#ifdef HAVE_IFPACK2_SHYLU_NODEHTS
          // This should not happen.
          isGblSuccess (false, out);
          return;
#else
          // This should happen. Continue with the default solver.
#endif
        }
      }
      else if (trisolverType == TrisolverDetails::KSPTRSV) {
        out << "Set solver parameters" << endl;
        // This line can throw. It should throw if HTS is not built in.
        Teuchos::ParameterList pl ("LocalSparseTriangularSolver parameters");
        pl.set ("trisolver: type", "KSPTRSV");
        try {
          solver->setParameters (pl);
        } catch (...) {
          isGblSuccess (false, out);
          return;
        }
      }

      out << "Set up the solver" << endl;
      TEST_NOTHROW( solver->initialize () );
      if (success) {
        TEST_NOTHROW( solver->compute () );
      }
      if (! isGblSuccess (success, out)) return;

      out << "Make sure that we can call setMatrix with a null input matrix, "
        "and that this resets the solver" << endl;
      TEST_NOTHROW( solver->setMatrix (Teuchos::null) );
      TEST_ASSERT( ! solver->isInitialized () );
      TEST_ASSERT( ! solver->isComputed () );
      if (! isGblSuccess (success, out)) return;

      out << "Set up the solver again with the original input matrix A" << endl;
      TEST_NOTHROW( solver->setMatrix (A) );
      TEST_NOTHROW( solver->initialize () );
      if (success) {
        TEST_NOTHROW( solver->compute () );
      }
      if (! isGblSuccess (success, out)) return;

      out << "Call compute() on the solver again, testing reuse of the nonzero pattern" << endl;
      TEST_NOTHROW( solver->compute () );

      out << "Test the solver" << endl;

      // For different combinations of alpha, beta, and mode, apply
      // the solver.  Test against Tpetra::CrsMatrix::localSolve with
      // A_copy, a deep copy of the original matrix A.  This tests
      // whether the solver changes A (it shouldn't).

      const std::size_t numVecValues[] = {
        static_cast<std::size_t> (1),
        static_cast<std::size_t> (3)
      };
      const Scalar ZERO = STS::zero ();
      const Scalar ONE = STS::one ();
      const Scalar THREE = ONE + ONE + ONE;
      const Scalar alphaValues[] = {ZERO, ONE, -ONE, THREE, -THREE};
      const Scalar betaValues[] = {ZERO, ONE, -ONE, THREE, -THREE};
      const Teuchos::ETransp modes[] = {Teuchos::NO_TRANS, Teuchos::TRANS, Teuchos::CONJ_TRANS};

      for (std::size_t numVecs : numVecValues) {
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
              if (! isGblSuccess (success, out)) return;

              // Test against a reference implementation.

              TEST_NOTHROW( referenceApply (Y_copy, *A_copy, X_copy,
                                            isUpperTriangular,
                                            implicitUnitDiag,
                                            mode, alpha, beta) );
              if (! isGblSuccess (success, out)) return;

              // Compare Y and Y_copy.  Compute difference in Y_copy.
              Y_copy.update (ONE, Y, -ONE); // Y_copy := Y - Y_copy
              Y_copy.normInf (norms);

              const mag_type tol = static_cast<mag_type> (gblNumRows) * STS::eps ();
              for (std::size_t j = 0; j < numVecs; ++j) {
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
  if (! isGblSuccess (success, out)) return;
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(LocalSparseTriangularSolver, CompareInternalToLocalSolve, Scalar, LocalOrdinal, GlobalOrdinal)
{
  testCompareToLocalSolve<Scalar, LocalOrdinal, GlobalOrdinal> (success, out, TrisolverDetails::Internal);
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(LocalSparseTriangularSolver, CompareKSPTRSVToLocalSolve, Scalar, LocalOrdinal, GlobalOrdinal)
{
  testCompareToLocalSolve<Scalar, LocalOrdinal, GlobalOrdinal> (success, out, TrisolverDetails::KSPTRSV);
}

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(LocalSparseTriangularSolver, CompareHTSToLocalSolve, Scalar, LocalOrdinal, GlobalOrdinal)
{
#ifdef HAVE_IFPACK2_SHYLU_NODEHTS
  testCompareToLocalSolve<Scalar, LocalOrdinal, GlobalOrdinal> (success, out, TrisolverDetails::HTS);
#else

#endif
}

// Consider a real arrow matrix (arrow pointing down and right) with
// diagonal entries d and other nonzero entries 1.  Here is a 4 x 4
// example:
//
// [d     1]
// [  d   1]
// [    d 1]
// [1 1 1 d]
//
// Compute its LU factorization without pivoting, assuming that all
// the values exist:
//
// [1            ] [d        1      ]
// [    1        ] [   d     1      ]
// [        1    ] [      d  1      ]
// [1/d 1/d 1/d 1] [         d - 3/d]
//
// Generalize the pattern: the off-diagonal nonzero entries of the L
// factor are all 1/d, and the lower right entry of U is d - (n-1)/d,
// where the original matrix A is n by n.  If d is positive and big
// enough, say d >= 2n, then all diagonal entries of U will be
// sufficiently large for this factorization to make sense.
// Furthermore, if d is a power of 2, then 1/d is exact in binary
// floating-point arithmetic (if it doesn't overflow), as is (1/d +
// 1/d).  This lets us easily check our work.
//
// Suppose that we want to solve Ax=b for b = [1 2 ... n]^T.
// For c = Ux, first solve Lc = b:
//
// c = [1, 2, ..., n-1, n - n(n-1)/(2d)]^T
//
// and then solve Ux = c.  First,
//
// x_n = c_n / (d - (n-1)/d).
//
// Then, for k = 1, ..., n-1, dx_k + x_n = k, so
//
// x_k = (k - x_n) / d, for k = 1, ..., n-1.
//
// Now, multiply b through by d - (n-1)/d.  This completely avoids
// rounding error, as long as no quantities overflow.  To get the
// right answer, multiply both c and x through by the same scaling
// factor.

template<class SC, class LO, class DT>
void
testArrowMatrixWithDense (bool& success, Teuchos::FancyOStream& out, const LO lclNumRows)
{
  using val_type = typename Kokkos::ArithTraits<SC>::val_type;
  using mag_type = typename Kokkos::ArithTraits<val_type>::mag_type;
  using host_execution_space =
    typename Kokkos::View<val_type**, DT>::HostMirror::execution_space;
  using host_memory_space = Kokkos::HostSpace;
  using HDT = Kokkos::Device<host_execution_space, host_memory_space>;

  Teuchos::OSTab tab0 (out);
  out << "Test arrow matrix problem using dense matrices" << endl;
  Teuchos::OSTab tab1 (out);

  const LO lclNumCols = lclNumRows;
  out << "Test with " << lclNumRows << " x " << lclNumCols << " matrices" << endl;

  // Kokkos Views fill with zeros by default.
  Kokkos::View<val_type**, HDT> A ("A", lclNumRows, lclNumCols);
  Kokkos::View<val_type**, HDT> L ("L", lclNumRows, lclNumCols);
  Kokkos::View<val_type**, HDT> U ("U", lclNumRows, lclNumCols);

  const val_type ZERO = Kokkos::ArithTraits<val_type>::zero ();
  const val_type ONE = Kokkos::ArithTraits<val_type>::one ();
  const val_type TWO = ONE + ONE;
  const val_type N = static_cast<val_type> (static_cast<mag_type> (lclNumRows));
  const val_type d = TWO * N;

  out << "Construct test problem (d = " << d << ")" << endl;

  for (LO i = 0; i < lclNumRows; ++i) {
    A(i, i) = d;
    if (i + 1 < lclNumRows) {
      A(i, lclNumCols-1) = ONE;
    }
    else if (i + 1 == lclNumRows) {
      for (LO j = 0; j + 1 < lclNumCols; ++j) {
        A(i, j) = ONE;
      }
    }
  }

  for (LO i = 0; i < lclNumRows; ++i) {
    L(i, i) = ONE;
    if (i + 1 == lclNumRows) {
      for (LO j = 0; j + 1 < lclNumCols; ++j) {
        L(i, j) = ONE / d;
      }
    }
  }

  for (LO i = 0; i < lclNumRows; ++i) {
    if (i + 1 < lclNumRows) {
      U(i, i) = d;
      U(i, lclNumCols-1) = ONE;
    }
    else if (i + 1 == lclNumRows) {
      U(i, i) = d - (N - ONE) / d;
    }
  }

  out << "Use dense matrix-matrix multiply to check that A == L*U" << endl;

  Kokkos::View<val_type**, HDT> A_copy ("A_copy", lclNumRows, lclNumCols);
  Tpetra::GEMM ("N", "N", ONE, L, U, ZERO, A_copy);
  for (LO i = 0; i < lclNumRows; ++i) {
    out << "Row " << i << endl;
    for (LO j = 0; j < lclNumCols; ++j) {
      TEST_EQUALITY( A(i,j), A_copy(i,j) );
    }
  }

  out << "Check that the LU factorization of A is LU" << endl;

  Kokkos::deep_copy (A_copy, A);
  Kokkos::View<LO*, HDT> ipiv ("ipiv", lclNumRows);
  int info = 0;
  Tpetra::GETF2 (A_copy, ipiv, info);
  TEST_EQUALITY( info, 0 );

  for (LO i = 0; i < lclNumRows; ++i) {
    out << "Row " << i << endl;
    for (LO j = 0; j < lclNumCols; ++j) {
      if (j < i) { // lower triangle (L) part of result
        TEST_EQUALITY( L(i,j), A_copy(i,j) );
      }
      else { // upper triangle (U) part of result
        TEST_EQUALITY( U(i,j), A_copy(i,j) );
      }
    }
  }
  for (LO i = 0; i < lclNumRows; ++i) {
    // LAPACK pivots are one-based.
    TEST_EQUALITY( static_cast<LO> (ipiv[i]), i + static_cast<LO> (1) );
  }

  // Test our exact solution to Ax=b.
  Kokkos::View<val_type*, HDT> x ("x", lclNumCols);
  Kokkos::View<val_type*, HDT> c ("c", lclNumRows);
  Kokkos::View<val_type*, HDT> b ("b", lclNumRows);

  const val_type scalingFactor = d - (N - ONE) / d;
  for (LO i = 0; i < lclNumRows; ++i) {
    b(i) = scalingFactor * static_cast<val_type> (static_cast<mag_type> (i+1));
    //b(i) = static_cast<val_type> (static_cast<mag_type> (i+1));
  }
  // GETRS overwrites the input right-hand side with the solution.
  Kokkos::deep_copy (x, b);

  Tpetra::GETRS ("N", A_copy, ipiv, x, info);
  TEST_EQUALITY( info, 0 );

  const val_type c_n_unscaled_expected = N - ((N - ONE)*N) / (TWO * d);
  const val_type c_n_expected = scalingFactor * c_n_unscaled_expected;
  const val_type x_n_unscaled_expected = c_n_unscaled_expected / scalingFactor;
  const val_type x_n_expected = c_n_unscaled_expected;

  //TEST_EQUALITY( x(lclNumRows-1), x_n_unscaled_expected );
  TEST_EQUALITY( x(lclNumRows-1), x_n_expected );
  for (LO i = 0; i + 1 < lclNumRows; ++i) {
    const val_type K = static_cast<val_type> (static_cast<mag_type> (i+1));
    //const val_type x_i_unscaled_expected = (K - x_n_unscaled_expected) / d;
    //TEST_EQUALITY( x(i), x_i_unscaled_expected );

    const val_type x_i_expected = (scalingFactor * (K - x_n_unscaled_expected)) / d;
    TEST_EQUALITY( x(i), x_i_expected );
  }

  // Now repeat the test to see if we can do a triangular solve with
  // L.  We do this to test whether solving Lc = b works correctly.
  Kokkos::deep_copy (x, b);
  Kokkos::deep_copy (c, ZERO);

  // Compute c, the result of solving Lc = b
  c(0) = b(0);
  for (LO i = 1; i < lclNumRows; ++i) {
    val_type c_i = b(i);
    for (LO j = 0; j < i; ++j) {
      c_i = c_i - L(i,j) * c(j);
    }
    c(i) = c_i / L(i,i);
  }

  // Test the resulting vector c
  for (LO i = 0; i + 1 < lclNumRows; ++i) {
    const val_type K = static_cast<val_type> (static_cast<mag_type> (i+1));
    // const val_type c_i_unscaled_expected = K; // unused
    const val_type c_i_expected = scalingFactor * K;
    TEST_EQUALITY( c(i), c_i_expected );
  }
  TEST_EQUALITY( c(lclNumRows-1), c_n_expected );
}


template<class crs_matrix_type, class map_type>
bool
testArrowMatrixAssembly(const int lclNumRows,
                        const bool explicitlyStoreUnitDiagonalOfL,
                        RCP<const map_type> rowMap,
                        RCP<const map_type> colMap,
                        RCP<const map_type> domMap,
                        RCP<const map_type> ranMap,
                        RCP<crs_matrix_type> & L,
                        RCP<crs_matrix_type> & U,
                        Teuchos::FancyOStream& out)
{
  int gblSuccess=1, lclSuccess=1;
  bool success=true;
  using LO = typename crs_matrix_type::local_ordinal_type;
  using SC = typename crs_matrix_type::scalar_type;

  typedef Kokkos::ArithTraits<SC> KAT;
  typedef typename KAT::val_type IST;
  typedef typename KAT::mag_type mag_type;
  typedef typename crs_matrix_type::local_graph_device_type local_graph_type;
  typedef typename crs_matrix_type::local_matrix_device_type local_matrix_type;
  typedef typename local_matrix_type::row_map_type::non_const_type row_offsets_type;
  typedef typename local_graph_type::entries_type::non_const_type col_inds_type;
  typedef typename local_matrix_type::values_type::non_const_type values_type;

  const LO lclNumCols = lclNumRows;

  auto comm = rowMap->getComm();

  //
  // The suffix _d here stands for (GPU) "device," and the suffix _h
  // stands for (CPU) "host."  
  //
  row_offsets_type L_ptr_d ("ptr", lclNumRows + 1);
  auto L_ptr_h = Kokkos::create_mirror_view (L_ptr_d);
  row_offsets_type U_ptr_d ("ptr", lclNumRows + 1);
  auto U_ptr_h = Kokkos::create_mirror_view (U_ptr_d);

  // The local number of _entries_ could in theory require 64 bits
  // even if LO is 32 bits.  This example doesn't require it, but why
  // not be general if there is no serious cost?  We use ptrdiff_t
  // because it is signed.
  const ptrdiff_t L_lclNumEnt = explicitlyStoreUnitDiagonalOfL ?
    (2*lclNumRows - 1) :
    (lclNumRows - 1);
  const ptrdiff_t U_lclNumEnt = 2*lclNumRows - 1;

  col_inds_type L_ind_d ("ind", L_lclNumEnt);
  auto L_ind_h = Kokkos::create_mirror_view (L_ind_d);
  values_type L_val_d ("val", L_lclNumEnt);
  auto L_val_h = Kokkos::create_mirror_view (L_val_d);

  col_inds_type U_ind_d ("ind", U_lclNumEnt);
  auto U_ind_h = Kokkos::create_mirror_view (U_ind_d);
  values_type U_val_d ("val", U_lclNumEnt);
  auto U_val_h = Kokkos::create_mirror_view (U_val_d);

  const IST ONE = KAT::one ();
  const IST TWO = KAT::one () + KAT::one ();
  // Don't cast directly from an integer type to IST,
  // since if IST is complex, that cast may not exist.
  const IST N = static_cast<IST> (static_cast<mag_type> (lclNumRows));
  const IST d = TWO * N;

  ptrdiff_t L_curPos = 0;
  for (LO i = 0; i < lclNumRows; ++i) {
    L_ptr_h[i] = L_curPos;

    if (i + 1 == lclNumRows) {
      // Last row: Add the off-diagonal entries
      for (LO j = 0; j + 1 < lclNumCols; ++j) {
        L_ind_h[L_curPos] = j;
        L_val_h[L_curPos] = ONE / d;
        ++L_curPos;
      }
    }
    if (explicitlyStoreUnitDiagonalOfL) {
      // Add the diagonal entry
      L_ind_h[L_curPos] = i;
      L_val_h[L_curPos] = ONE;
      ++L_curPos;
    }
  }
  L_ptr_h[lclNumRows] = L_curPos;

  ptrdiff_t U_curPos = 0;
  for (LO i = 0; i < lclNumRows; ++i) {
    U_ptr_h[i] = U_curPos;

    if (i + 1 < lclNumRows) {
      // Add the diagonal entry (first in the row)
      U_ind_h[U_curPos] = i;
      U_val_h[U_curPos] = d;
      ++U_curPos;

      // Add the last entry in the row
      U_ind_h[U_curPos] = lclNumCols - 1;
      U_val_h[U_curPos] = ONE;
      ++U_curPos;
    }
    else if (i + 1 == lclNumRows) {
      // Add the last row's diagonal entry (only entry in this row)
      U_ind_h[U_curPos] = lclNumCols - 1;
      U_val_h[U_curPos] = d - (N - ONE) / d;
      ++U_curPos;
    }
  }
  U_ptr_h[lclNumRows] = U_curPos;

  // Make sure that we counted the number of entries correctly.
  TEST_ASSERT( L_curPos == L_lclNumEnt );
  TEST_ASSERT( U_curPos == U_lclNumEnt );
  lclSuccess = success ? 1 : 0;
  gblSuccess = 0; // output argument
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  TEST_EQUALITY( gblSuccess, 1 );
  if (! gblSuccess) {
    out << "Aborting test" << endl;
    return gblSuccess;
  }

  Kokkos::deep_copy (L_ptr_d, L_ptr_h);
  Kokkos::deep_copy (L_ind_d, L_ind_h);
  Kokkos::deep_copy (L_val_d, L_val_h);

  Kokkos::deep_copy (U_ptr_d, U_ptr_h);
  Kokkos::deep_copy (U_ind_d, U_ind_h);
  Kokkos::deep_copy (U_val_d, U_val_h);

  out << "Create the lower triangular Tpetra::CrsMatrix L" << endl;
  TEST_NOTHROW( L = rcp (new crs_matrix_type (rowMap, colMap, L_ptr_d, L_ind_d, L_val_d)) );
  TEST_ASSERT( ! L.is_null () );
  lclSuccess = success ? 1 : 0;
  gblSuccess = 0; // output argument
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  TEST_EQUALITY( gblSuccess, 1 );
  if (! gblSuccess) {
    out << "Aborting test" << endl;
    return gblSuccess;
  }
  out << "Call fillComplete on the lower triangular matrix L" << endl;
  TEST_NOTHROW( L->fillComplete (domMap, ranMap) );
  lclSuccess = success ? 1 : 0;
  gblSuccess = 0; // output argument
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  TEST_EQUALITY( gblSuccess, 1 );
  if (! gblSuccess) {
    out << "Aborting test" << endl;
    return gblSuccess;
  }

  out << "Create the upper triangular Tpetra::CrsMatrix U" << endl;
  TEST_NOTHROW( U = rcp (new crs_matrix_type (rowMap, colMap, U_ptr_d, U_ind_d, U_val_d)) );
  TEST_ASSERT( ! U.is_null () );
  lclSuccess = success ? 1 : 0;
  gblSuccess = 0; // output argument
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  TEST_EQUALITY( gblSuccess, 1 );
  if (! gblSuccess) {
    out << "Aborting test" << endl;
    return gblSuccess;
  }
  out << "Call fillComplete on the upper triangular matrix U" << endl;
  TEST_NOTHROW( U->fillComplete (domMap, ranMap) );
  lclSuccess = success ? 1 : 0;
  gblSuccess = 0; // output argument
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  TEST_EQUALITY( gblSuccess, 1 );
  if (! gblSuccess) {
    out << "Aborting test" << endl;
    return gblSuccess;
  }
  return gblSuccess;
}


template<class SC = Tpetra::Vector<>::scalar_type,
         class LO = Tpetra::Vector<>::local_ordinal_type,
         class GO = Tpetra::Vector<>::global_ordinal_type>
void testArrowMatrix (bool& success, Teuchos::FancyOStream& out)
{
  typedef Tpetra::Map<LO, GO> map_type;
  typedef typename map_type::device_type device_type;
  typedef Tpetra::CrsMatrix<SC, LO, GO> crs_matrix_type;
  typedef Tpetra::RowMatrix<SC, LO, GO> row_matrix_type;
  typedef Tpetra::Vector<SC, LO, GO> vec_type;
  typedef Ifpack2::LocalSparseTriangularSolver<row_matrix_type> solver_type;
  typedef Kokkos::ArithTraits<SC> KAT;
  typedef typename KAT::val_type IST;
  typedef typename KAT::mag_type mag_type;
  int lclSuccess = 1;
  int gblSuccess = 1;

  const bool explicitlyStoreUnitDiagonalOfL = false;

  Teuchos::OSTab tab0 (out);
  out << "Ifpack2::LocalSparseTriangularSolver: Test with arrow matrix" << endl;
  Teuchos::OSTab tab1 (out);

  auto comm = Tpetra::getDefaultComm ();

  const LO lclNumRows = 8; // power of two (see above)
  const LO lclNumCols = lclNumRows;
  const GO gblNumRows = comm->getSize () * lclNumRows;
  const GO indexBase = 0;
  RCP<const map_type> rowMap =
    rcp (new map_type (static_cast<GST> (gblNumRows),
                       static_cast<std::size_t> (lclNumRows),
                       indexBase, comm));

  // At this point, we know Kokkos has been initialized, so test the
  // dense version of the problem.
  testArrowMatrixWithDense<SC, LO, device_type> (success, out, lclNumRows);

  // If we construct an upper or lower triangular matrix with an
  // implicit unit diagonal, then we need to specify the column Map
  // explicitly.  Otherwise, the matrix will report having the wrong
  // number of columns.  In this case, the local matrix is square and
  // every column is populated, so we can set column Map = row Map.
  RCP<const map_type> colMap = rowMap;
  RCP<const map_type> domMap = rowMap;
  RCP<const map_type> ranMap = rowMap;

  // All of the matrix assembly stuff had to get hived off into a different
  // scope to keep the later accessors from violating the "you can't have a 
  // host and a device view at the same time" assumption
  RCP<crs_matrix_type> L, U;

  gblSuccess=testArrowMatrixAssembly(lclNumRows,
                                     explicitlyStoreUnitDiagonalOfL,
                                     rowMap,colMap,domMap,ranMap,
                                     L,U,out);
  if(!gblSuccess) return;

  typedef typename crs_matrix_type::local_inds_host_view_type const_local_inds_type;
  typedef typename crs_matrix_type::values_host_view_type const_values_type;

  const IST ONE = KAT::one ();
  const IST TWO = KAT::one () + KAT::one ();
  // Don't cast directly from an integer type to IST,
  // since if IST is complex, that cast may not exist.
  const IST N = static_cast<IST> (static_cast<mag_type> (lclNumRows));
  const IST d = TWO * N;

  out << "Make sure that the last row of L is correct" << endl;
  {
    Teuchos::OSTab tab2 (out);

    const_local_inds_type lclColInds;
    const_values_type vals;

    L->getLocalRowView (lclNumRows - 1, lclColInds, vals);
    if (explicitlyStoreUnitDiagonalOfL) {
      TEST_EQUALITY( static_cast<LO> (lclColInds.size ()), lclNumCols );
      TEST_EQUALITY( static_cast<LO> (vals.size ()), lclNumCols );
    }
    else {
      TEST_EQUALITY( static_cast<LO> (lclColInds.size ()),
                     lclNumCols - static_cast<LO> (1) );
      TEST_EQUALITY( static_cast<LO> (vals.size ()),
                     lclNumCols - static_cast<LO> (1) );
    }
    if (success) {
      // FIXME (mfh 23 Aug 2016) This depends on the entries being
      // sorted.  They should be, since they were sorted on input,
      // using a KokkosSparse::CrsMatrix input.
      for (LO j = 0; j + 1 < lclNumCols; ++j) {
        TEST_EQUALITY( lclColInds[j], j );
        TEST_EQUALITY( static_cast<IST> (vals[j]), ONE / d );
      }
      if (explicitlyStoreUnitDiagonalOfL) {
        TEST_EQUALITY( lclColInds[lclNumCols-1], lclNumCols - 1 );
        TEST_EQUALITY( static_cast<IST> (vals[lclNumCols-1]), ONE );
      }
    }
  }
  lclSuccess = success ? 1 : 0;
  gblSuccess = 0; // output argument
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  TEST_EQUALITY( gblSuccess, 1 );
  if (! gblSuccess) {
    out << "Aborting test" << endl;
    return;
  }


  out << "Create the solver for L" << endl;
  RCP<solver_type> L_solver;
  TEST_NOTHROW( L_solver = rcp (new solver_type (L, Teuchos::rcpFromRef (out))) );
  TEST_ASSERT( ! L_solver.is_null () );
  lclSuccess = success ? 1 : 0;
  gblSuccess = 0; // to be revised
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  TEST_EQUALITY( gblSuccess, 1 );
  if (! gblSuccess) {
    out << "Aborting test" << endl;
    return;
  }

  out << "Set up the solver for L" << endl;
  TEST_NOTHROW( L_solver->initialize () );
  if (success) {
    TEST_NOTHROW( L_solver->compute () );
  }
  lclSuccess = success ? 1 : 0;
  gblSuccess = 0; // to be revised
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  TEST_EQUALITY( gblSuccess, 1 );
  if (! gblSuccess) {
    out << "Aborting test" << endl;
    return;
  }

  out << "Create the solver for U" << endl;
  RCP<solver_type> U_solver;
  TEST_NOTHROW( U_solver = rcp (new solver_type (U, Teuchos::rcpFromRef (out))) );
  TEST_ASSERT( ! U_solver.is_null () );
  lclSuccess = success ? 1 : 0;
  gblSuccess = 0; // to be revised
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  TEST_EQUALITY( gblSuccess, 1 );
  if (! gblSuccess) {
    out << "Aborting test" << endl;
    return;
  }

  out << "Set up the solver for U" << endl;
  TEST_NOTHROW( U_solver->initialize () );
  if (success) {
    TEST_NOTHROW( U_solver->compute () );
  }
  lclSuccess = success ? 1 : 0;
  gblSuccess = 0; // to be revised
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  TEST_EQUALITY( gblSuccess, 1 );
  if (! gblSuccess) {
    out << "Aborting test" << endl;
    return;
  }

  const IST scalingFactor = d - (N - ONE) / d;
  const bool scaleProblem = true;//false;

  // Set up the right-hand side b.
  vec_type b (ranMap);
  {
    auto b_lcl_2d = b.getLocalViewHost (Tpetra::Access::OverwriteAll);
    auto b_lcl_1d = Kokkos::subview (b_lcl_2d, Kokkos::ALL (), 0);

    for (LO i = 0; i < lclNumRows; ++i) {
      // Don't cast directly from an integer type to IST,
      // since if IST is complex, that cast may not exist.
      const IST K = static_cast<IST> (static_cast<mag_type> (i+1));
      if (scaleProblem) {
        b_lcl_1d(i) = scalingFactor * K;
      }
      else {
        b_lcl_1d(i) = K;
      }
    }
  }

  // We solve Ax=b (with A = LU) by first solving Lc = b, and then
  // solving Ux = c.  Thus, c's Map is the same as U's range Map.
  vec_type c (U->getRangeMap ());
  vec_type x (domMap);

  out << "Solve Lc = b for c" << endl;
  int lclNotThrew = 1; // to be set below
  int gblNotThrew = 0; // output argument
  {
    std::ostringstream errStrm;
    bool threw = false;
    try {
      L_solver->apply (b, c);
    }
    catch (std::exception& e) {
      errStrm << "L_solver->apply(b,c) threw an exception: " << e.what () << std::endl;
      threw = true;
    }
    catch (...) {
      errStrm << "L_solver->apply(b,c) threw an exception not a subclass of "
        "std::exception" << std::endl;
      threw = true;
    }
    lclNotThrew = threw ? 0 : 1;
    reduceAll<int, int> (*comm, REDUCE_MIN, lclNotThrew, outArg (gblNotThrew));
    TEST_EQUALITY( gblNotThrew, 1 );
    if (gblNotThrew != 1) {
      Tpetra::Details::gathervPrint (out, errStrm.str (), *comm);
    }
  }
  if (gblNotThrew != 1) {
    out << "Aborting test" << endl;
    return;
  }

  // Test the entries of c for correctness.  These are EXACT tests,
  // which we may do since the solves should not have committed any
  // rounding error.  See discussion above.
  //
  // NOTE (mfh 21 Aug 2016) This won't work if we accept approximate
  // sparse triangular solves.

  const IST c_n_unscaled_expected = N - ((N - ONE)*N) / (TWO * d);
  const IST c_n_expected = scaleProblem ?
    (scalingFactor * c_n_unscaled_expected) :
    c_n_unscaled_expected;
  const IST x_n_unscaled_expected = c_n_unscaled_expected / scalingFactor;
  const IST x_n_expected = scaleProblem ? c_n_unscaled_expected : x_n_unscaled_expected;

  out << "Test entries of c (solution of Lc=b)" << endl;
  {
    Teuchos::OSTab tab2 (out);

    auto c_lcl_2d = c.getLocalViewHost (Tpetra::Access::ReadOnly);
    auto c_lcl_1d = Kokkos::subview (c_lcl_2d, Kokkos::ALL (), 0);

    for (LO i = 0; i + 1 < lclNumRows; ++i) {
      // Don't cast directly from an integer type to IST,
      // since if IST is complex, that cast may not exist.
      const IST K = static_cast<IST> (static_cast<mag_type> (i+1));
      const IST c_i_expected = scaleProblem ? (scalingFactor * K) : K;
      TEST_EQUALITY( c_lcl_1d(i), c_i_expected );
    }
    TEST_EQUALITY( c_lcl_1d(lclNumRows-1), c_n_expected );
  }
  // lclSuccess = success ? 1 : 0;
  // gblSuccess = 0; // to be revised
  // reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  // TEST_EQUALITY( gblSuccess, 1 );
  // if (! gblSuccess) {
  //   out << "Aborting test" << endl;
  //   return;
  // }

  out << "Solve Ux = c for x" << endl;
  TEST_NOTHROW( U_solver->apply (c, x) );
  // lclSuccess = success ? 1 : 0;
  // gblSuccess = 0; // to be revised
  // reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  // TEST_EQUALITY( gblSuccess, 1 );
  // if (! gblSuccess) {
  //   out << "Aborting test" << endl;
  //   return;
  // }

  // Test the entries of x for correctness.  These are EXACT tests,
  // which we may do since the solves should not have committed any
  // rounding error.  See discussion above.
  //
  // NOTE (mfh 21 Aug 2016) This won't work if we accept approximate
  // sparse triangular solves.

  out << "Test entries of x (solution of Ux=c)" << endl;
  {
    Teuchos::OSTab tab2 (out);

    auto x_lcl_2d = x.getLocalViewHost (Tpetra::Access::ReadOnly);
    auto x_lcl_1d = Kokkos::subview (x_lcl_2d, Kokkos::ALL (), 0);

    for (LO i = 0; i + 1 < lclNumRows; ++i) {
      // Don't cast directly from an integer type to IST,
      // since if IST is complex, that cast may not exist.
      const IST K = static_cast<IST> (static_cast<mag_type> (i+1));
      const IST x_i_expected = scaleProblem ?
        ((scalingFactor * (K - x_n_unscaled_expected)) / d) :
        ((K - x_n_unscaled_expected) / d);
      TEST_EQUALITY( x_lcl_1d(i), x_i_expected );
    }
    TEST_EQUALITY( x_lcl_1d(lclNumRows-1), x_n_expected );
  }

  out << "Test against a reference sparse triangular solver" << endl;

  c.putScalar (Teuchos::ScalarTraits<SC>::zero ());
  x.putScalar (Teuchos::ScalarTraits<SC>::zero ());

  const std::string unitDiagL = explicitlyStoreUnitDiagonalOfL ?
    "No unit diagonal" : "Unit diagonal";
  localSolve (c, *L, b, false, ! explicitlyStoreUnitDiagonalOfL, Teuchos::NO_TRANS);
  out << "Test entries of c (solution of Lc=b)" << endl;
  {
    Teuchos::OSTab tab2 (out);

    auto c_lcl_2d = c.getLocalViewHost (Tpetra::Access::ReadOnly);
    auto c_lcl_1d = Kokkos::subview (c_lcl_2d, Kokkos::ALL (), 0);

    for (LO i = 0; i + 1 < lclNumRows; ++i) {
      // Don't cast directly from an integer type to IST,
      // since if IST is complex, that cast may not exist.
      const IST K = static_cast<IST> (static_cast<mag_type> (i+1));
      const IST c_i_expected = scaleProblem ? (scalingFactor * K) : K;
      TEST_EQUALITY( c_lcl_1d(i), c_i_expected );
    }
    TEST_EQUALITY( c_lcl_1d(lclNumRows-1), c_n_expected );
  }

  localSolve (x, *U, c, true, false, Teuchos::NO_TRANS);
  out << "Test entries of x (solution of Ux=c)" << endl;
  {
    Teuchos::OSTab tab2 (out);

    auto x_lcl_2d = x.getLocalViewHost (Tpetra::Access::ReadOnly);
    auto x_lcl_1d = Kokkos::subview (x_lcl_2d, Kokkos::ALL (), 0);

    for (LO i = 0; i + 1 < lclNumRows; ++i) {
      // Don't cast directly from an integer type to IST,
      // since if IST is complex, that cast may not exist.
      const IST K = static_cast<IST> (static_cast<mag_type> (i+1));
      const IST x_i_expected = scaleProblem ?
        ((scalingFactor * (K - x_n_unscaled_expected)) / d) :
        ((K - x_n_unscaled_expected) / d);
      TEST_EQUALITY( x_lcl_1d(i), x_i_expected );
    }
    TEST_EQUALITY( x_lcl_1d(lclNumRows-1), x_n_expected );
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

// This test is really only useful for Scalar = double.
#ifdef HAVE_TPETRA_INST_DOUBLE
TEUCHOS_UNIT_TEST(LocalSparseTriangularSolver, ArrowMatrix)
{
  testArrowMatrix<double> (success, out);
}
#endif // HAVE_TPETRA_INST_DOUBLE

#define UNIT_TEST_GROUP_SC_LO_GO(SC, LO, GO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(LocalSparseTriangularSolver, CompareInternalToLocalSolve, SC, LO, GO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(LocalSparseTriangularSolver, CompareKSPTRSVToLocalSolve, SC, LO, GO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(LocalSparseTriangularSolver, CompareHTSToLocalSolve, SC, LO, GO)


#include "Ifpack2_ETIHelperMacros.h"

IFPACK2_ETI_MANGLING_TYPEDEFS()

// Test all enabled combinations of Scalar (SC), LocalOrdinal (LO),
// and GlobalOrdinal (GO) types, where Scalar is real.

IFPACK2_INSTANTIATE_SLG_REAL( UNIT_TEST_GROUP_SC_LO_GO )

} // namespace (anonymous)

