// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#include "Teuchos_UnitTestHarness.hpp"
#include "Ifpack2_LocalSparseTriangularSolver.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Details_determineLocalTriangularStructure.hpp"
#include "Teuchos_CommHelpers.hpp"
#include <type_traits>

namespace {

  template<class LO, class GO, class NT>
  Tpetra::Details::LocalTriangularStructureResult<LO>
  getLocalTriangularStructure (const Tpetra::RowGraph<LO, GO, NT>& G)
  {
    using Tpetra::Details::determineLocalTriangularStructure;
    using crs_graph_type = Tpetra::CrsGraph<LO, GO, NT>;

    const crs_graph_type& G_crs = dynamic_cast<const crs_graph_type&> (G);

    auto G_lcl = G_crs.getLocalGraphDevice ();
    auto lclRowMap = G.getRowMap ()->getLocalMap ();
    auto lclColMap = G.getColMap ()->getLocalMap ();
    return determineLocalTriangularStructure (G_lcl, lclRowMap, lclColMap, true);
  }

  template<class SC, class LO, class GO, class NT>
  Tpetra::Details::LocalTriangularStructureResult<LO>
  getLocalTriangularStructure (const Tpetra::CrsMatrix<SC, LO, GO, NT>& A)
  {
    return getLocalTriangularStructure (* (A.getGraph ()));
  }

  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::as;
  using Teuchos::Comm;
  using Teuchos::outArg;
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::REDUCE_SUM;
  using Teuchos::reduceAll;
  using Teuchos::tuple;

  using Teuchos::ETransp;
  using Teuchos::CONJ_TRANS;
  using Teuchos::NO_TRANS;
  using Teuchos::TRANS;

  using Teuchos::EDiag;
  using Teuchos::UNIT_DIAG;
  using Teuchos::NON_UNIT_DIAG;
  using Teuchos::EUplo;
  using Teuchos::LOWER_TRI;
  using Teuchos::UNDEF_TRI;
  using Teuchos::UPPER_TRI;

  using std::endl;

  using GST = Tpetra::global_size_t;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
  }

  //
  // UNIT TESTS
  //

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, EmptyTriSolve, Scalar, LO, GO, Node )
  {
    using crs_matrix_type = Tpetra::CrsMatrix<Scalar, LO, GO, Node>;
    using row_matrix_type = Tpetra::RowMatrix<Scalar, LO, GO, Node>;
    using solver_type = Ifpack2::LocalSparseTriangularSolver<row_matrix_type>;
    using STS = Teuchos::ScalarTraits<Scalar>;
    using MV = Tpetra::MultiVector<Scalar,LO,GO,Node>;
    using mag_type = typename STS::magnitudeType;
    using STM = Teuchos::ScalarTraits<mag_type>;
    using map_type = Tpetra::Map<LO, GO, Node>;

    const size_t numLocal = 13, numVecs = 7;
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
    // get a comm
    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
    // create a Map
    RCP<const map_type> map =
      Tpetra::createContigMapWithNode<LO, GO, Node> (INVALID, numLocal, comm);

    /* Create a triangular matrix with no entries, for testing implicit diagonals.
      We test with Transpose and Non-Transpose application solve (these should be equivalent for the identity matrix)
    */

    MV X(map,numVecs), B(map,numVecs), Xhat(map,numVecs);
    X.setObjectLabel("X");
    B.setObjectLabel("B");
    Xhat.setObjectLabel("Xhat");
    X.randomize();
    for (size_t tnum=0; tnum < 2; ++tnum) {
      ETransp trans     = ((tnum & 1) == 1 ? CONJ_TRANS        : NO_TRANS);
      RCP<solver_type> ZeroIOp;
      {
        RCP<crs_matrix_type> ZeroMat;
        // must explicitly provide the column map for implicit diagonals
        ZeroMat = rcp(new crs_matrix_type(map,map,0));
        RCP<ParameterList> params = parameterList();
        RCP<ParameterList> fillparams = sublist(params,"Local Sparse Ops");
        fillparams->set("Prepare Solve", true);
        fillparams->set("Prepare Transpose Solve", true);
        fillparams->set("Prepare Conjugate Transpose Solve", true);
        ZeroMat->fillComplete(params);

        auto lclTri = getLocalTriangularStructure (*ZeroMat);
        TEST_ASSERT( lclTri.couldBeLowerTriangular );
        TEST_ASSERT( lclTri.couldBeUpperTriangular );
        GO gblDiagCount = 0;
        reduceAll<int, GO> (*comm, REDUCE_SUM,
                            static_cast<GO> (lclTri.diagCount),
                            outArg (gblDiagCount));
        TEST_EQUALITY( gblDiagCount, static_cast<GO> (0) );

        ZeroIOp = rcp<solver_type> (new solver_type (ZeroMat.getConst ()));
        ZeroIOp->initialize ();
        ZeroIOp->compute ();
      }
      X = B;
      Xhat.randomize();
      ZeroIOp->apply(B,Xhat,trans);
      //
      Xhat.update(-STS::one(),X,STS::one());
      Array<mag_type> errnrms(numVecs), normsB(numVecs), zeros(numVecs, STM::zero());
      Xhat.norm1(errnrms());
      B.norm1(normsB());
      mag_type maxBnrm = *std::max_element( normsB.begin(), normsB.end() );
      if (std::is_integral<Scalar>::value) {
        TEST_COMPARE_ARRAYS(errnrms, zeros);
      }
      else {
        TEST_COMPARE_FLOATING_ARRAYS( errnrms, zeros, maxBnrm );
      }
    }
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, TriSolve, Scalar, LO, GO, Node )
  {
    out << "Testing Tpetra::CrsMatrix triangular solve with nonempty matrices"
        << endl;
    Teuchos::OSTab tab0 (out);

    // mfh 26 Feb 2014: Organizing the if-else in this way avoids a
    // build warning for "dynamic initialization in unreachable code."
    if (std::is_integral<Scalar>::value) {
      out << "Skipping testing for the integral type Scalar="
          << Teuchos::TypeNameTraits<Scalar>::name () << "." << endl;
      return;
    }
    else {
      using crs_matrix_type = Tpetra::CrsMatrix<Scalar,LO,GO,Node>;
      using row_matrix_type = Tpetra::RowMatrix<Scalar,LO,GO,Node>;
      using solver_type = Ifpack2::LocalSparseTriangularSolver<row_matrix_type>;
      using map_type = Tpetra::Map<LO, GO, Node>;
      using STS = Teuchos::ScalarTraits<Scalar>;
      using MV = Tpetra::MultiVector<Scalar,LO,GO,Node>;
      using mag_type = typename STS::magnitudeType;
      using STM = Teuchos::ScalarTraits<mag_type>;

      const size_t numLocal = 13, numVecs = 7;
      const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
      RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
      // Create a row Map for the matrix.
      // This will be the same as the domain and range Maps.
      RCP<const map_type> map =
        Tpetra::createContigMapWithNode<LO, GO, Node> (INVALID, numLocal, comm);
      Scalar SONE = STS::one ();

      /* Create one of the following locally triangular matries:

         0  [1 2       ]
         1  [  1 3     ]
         .  [    .  .  ] = U
        n-2 [       1 n]
        n-1 [         1]

         0  [1           ]
         1  [2 1         ]
         .  [   .  .     ] = L
        n-2 [     n-1 1  ]
        n-1 [         n 1]

        The resulting global matrices are diag(U,U,...,U) resp. diag(L,L,...,L).

        For each of these (upper or lower triangular), we test all
        16 combinations of the following options:
        - Explicit or implicit unit diagonal
        - With and without the (conjugate) transpose
        - Optimized or nonoptimized storage
      */

      MV X (map, numVecs), B (map, numVecs), Xhat (map, numVecs);
      X.setObjectLabel("X");
      B.setObjectLabel("B");
      Xhat.setObjectLabel("Xhat");
      X.randomize();

      // Sanity check for X.
      Array<mag_type> normsX (numVecs);
      {
        X.norm1 (normsX ());
        Array<size_t> badColumns;
        for (size_t j = 0; j < numVecs; ++j) {
          if (STS::isnaninf (normsX[j])) {
            badColumns.push_back (j);
          }
        }
        TEUCHOS_TEST_FOR_EXCEPTION(
          badColumns.size () > 0,
          std::runtime_error,
          "Columns " << Teuchos::toString (badColumns) << " of the input X "
          "have a 1-norm either Inf or NaN.  That suggests that randomize() "
          "is broken.  Here are the 1-norms of each column: "
          << Teuchos::toString (normsX));
      }

      RCP<ParameterList> params = parameterList();
      // Test all 16 combinations of options.
      for (size_t tnum = 0; tnum < 16; ++tnum) {
        const EUplo   uplo  = ((tnum & 1) == 1 ? UPPER_TRI  : LOWER_TRI);
        const EDiag   diag  = ((tnum & 2) == 2 ? UNIT_DIAG  : NON_UNIT_DIAG);
        const ETransp trans = ((tnum & 8) == 8 ? CONJ_TRANS : NO_TRANS);
        const bool optimizeStorage = (tnum & 4) == 4;

        std::string diagStr;
        if (diag == UNIT_DIAG) {
          diagStr = "UNIT_DIAG";
        } else if (diag == NON_UNIT_DIAG) {
          diagStr = "NON_UNIT_DIAG";
        } else {
          diagStr = "UNDEFINED";
        }
        std::string uploStr;
        if (uplo == LOWER_TRI) {
          uploStr = "LOWER_TRI";
        } else if (uplo == UPPER_TRI) {
          uploStr = "UPPER_TRI";
        } else if (uplo == UNDEF_TRI) {
          uploStr = "UNDEF_TRI";
        } else {
          uploStr = "UNDEFINED";
        }
        std::string transStr;
        if (trans == CONJ_TRANS) {
          transStr = "CONJ_TRANS";
        } else if (trans == TRANS) {
          transStr = "TRANS";
        } else if (trans == NO_TRANS) {
          transStr = "NO_TRANS";
        } else {
          transStr = "UNDEFINED";
        }

        out << "Test " << (tnum+1) << " of " << 16 << ":" << endl;
        Teuchos::OSTab tab1 (out);
        {
          out << "Parameters:" << endl;
          Teuchos::OSTab tab2 (out);
          out << "uplo: " << uploStr << endl
              << "diag: " << diagStr << endl
              << "trans: " << transStr << endl
              << "optimizeStorage: " << (optimizeStorage ? "true" : "false") << endl
              << endl;
        }

        params->set ("Optimize Storage", optimizeStorage);
        RCP<ParameterList> fillparams = sublist (params, "Local Sparse Ops");
        fillparams->set ("Prepare Solve", true);
        fillparams->set ("Prepare Transpose Solve", true);
        fillparams->set ("Prepare Conjugate Transpose Solve", true);

        RCP<solver_type> AIOp;
        RCP<crs_matrix_type> AMat;
        {
          if (diag == UNIT_DIAG) {
            // must explicitly specify the column map
            AMat = rcp(new crs_matrix_type(map,map,2));
          }
          else {
            // can let the matrix compute a column map
            AMat = rcp(new crs_matrix_type(map,2));
          }
          // fill the matrix
          if (uplo == UPPER_TRI) {
            if (diag == UNIT_DIAG) {
              for (GO gid=map->getMinGlobalIndex(); gid <= map->getMaxGlobalIndex(); ++gid) {
                if (gid == map->getMaxGlobalIndex()) {
                  // do nothing
                }
                else {
                  AMat->insertGlobalValues (gid, tuple<GO> (gid+1), tuple<Scalar> (as<Scalar> (gid+2)));
                }
              }
            }
            else {
              for (GO gid=map->getMinGlobalIndex(); gid <= map->getMaxGlobalIndex(); ++gid) {
                if (gid == map->getMaxGlobalIndex()) {
                  AMat->insertGlobalValues (gid, tuple<GO> (gid), tuple<Scalar> (SONE));
                }
                else {
                  AMat->insertGlobalValues (gid, tuple<GO> (gid,gid+1), tuple<Scalar> (SONE, as<Scalar>(gid+2)));
                }
              }
            }
          }
          else { // uplo == LOWER_TRI
            if (diag == UNIT_DIAG) {
              for (GO gid=map->getMinGlobalIndex(); gid <= map->getMaxGlobalIndex(); ++gid) {
                if (gid == map->getMinGlobalIndex()) {
                  // do nothing
                }
                else {
                  AMat->insertGlobalValues (gid, tuple<GO> (gid-1), tuple<Scalar> (as<Scalar> (gid+1)));
                }
              }
            }
            else {
              for (GO gid=map->getMinGlobalIndex(); gid <= map->getMaxGlobalIndex(); ++gid) {
                if (gid == map->getMinGlobalIndex()) {
                  AMat->insertGlobalValues (gid, tuple<GO> (gid), tuple<Scalar> (SONE));
                }
                else {
                  AMat->insertGlobalValues (gid, tuple<GO> (gid-1, gid), tuple<Scalar> (as<Scalar> (gid+1), SONE));
                }
              }
            }
          }
          AMat->fillComplete(params);

          auto lclTri = getLocalTriangularStructure (*AMat);
          TEST_EQUALITY( lclTri.couldBeLowerTriangular, uplo == LOWER_TRI );
          TEST_EQUALITY( lclTri.couldBeUpperTriangular, uplo == UPPER_TRI );

          GO gblDiagCount = 0;
          reduceAll<int, GO> (*comm, REDUCE_SUM,
                              static_cast<GO> (lclTri.diagCount),
                              outArg (gblDiagCount));
          TEST_EQUALITY( gblDiagCount == static_cast<GO> (0), diag == UNIT_DIAG );

          // AIOp.apply (X,B,trans) solves op(A) X=B for X locally,
          // using a triangular solve.  op(A) is just A if
          // trans==NO_TRANS, else A^H (Hermitian transpose) if
          // trans==CONJ_TRANS.

          AIOp = rcp<solver_type> (new solver_type (AMat.getConst ()));
          AIOp->initialize ();
          AIOp->compute ();
        }
        B.randomize ();
        AIOp->apply (X, B, trans);
        if (diag == UNIT_DIAG) {
          // we want (I+A)*X -> B
          // A*X -> B needs to be augmented with X
          B.update (STS::one (), X, STS::one());
        }
        Array<mag_type> normsB (numVecs);
        B.norm1 (normsB ());
        {
          Array<size_t> badColumns;
          for (size_t j = 0; j < numVecs; ++j) {
            if (STS::isnaninf (normsB[j])) {
              badColumns.push_back (j);
            }
          }
          if (badColumns.size () > 0) {
            out << "Result of triangular solve contains Inf or NaN:" << endl;
            Teuchos::OSTab tab2 (out);
            B.normInf (normsB ());
            out << "Inf-norms of each column of B = A \\ X: "
                << Teuchos::toString (normsB) << endl;
            B.norm2 (normsB ());
            out << "2-norms of each column of B = A \\ X: "
                << Teuchos::toString (normsB) << endl;
            B.norm1 (normsB ());
            out << "1-norms of each column of B = A \\ X: "
                << Teuchos::toString (normsB) << endl;
            out << "Input MV X:" << endl;
            X.describe (out, Teuchos::VERB_EXTREME);
            out << "Output MV B, on output:" << endl;
            B.describe (out, Teuchos::VERB_EXTREME);
          }
        }

        Xhat.randomize ();
        AIOp->apply (B, Xhat, trans);
        Xhat.update (-STS::one (), X, STS::one ());
        Array<mag_type> errnrms (numVecs), zeros (numVecs, STM::zero ());
        Xhat.norm1 (errnrms ());
        B.norm1 (normsB ());
        mag_type maxBnrm = *std::max_element (normsB.begin (), normsB.end ());
        if (std::is_integral<Scalar>::value) {
          TEST_COMPARE_ARRAYS( errnrms, zeros );
        } else {
          TEST_COMPARE_FLOATING_ARRAYS( errnrms, zeros, maxBnrm );
        }
      }
    }
  }

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, EmptyTriSolve, SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, TriSolve, SCALAR, LO, GO, NODE )

#include "Ifpack2_ETIHelperMacros.h"

  IFPACK2_ETI_MANGLING_TYPEDEFS()

  IFPACK2_INSTANTIATE_SLGN( UNIT_TEST_GROUP )

}


