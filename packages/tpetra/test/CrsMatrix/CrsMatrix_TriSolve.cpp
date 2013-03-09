/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER
*/

#include <Tpetra_TestingUtilities.hpp>

#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_CrsMatrixSolveOp.hpp>

namespace {

  using Tpetra::TestingUtilities::getNode;
  using Tpetra::TestingUtilities::getDefaultComm;

  using std::endl;
  using std::string;

  using Teuchos::as;
  using Teuchos::FancyOStream;
  using Teuchos::RCP;
  using Teuchos::ArrayRCP;
  using Teuchos::rcp;
  using Teuchos::arcp;
  using Teuchos::outArg;
  using Teuchos::arcpClone;
  using Teuchos::arrayView;
  using Teuchos::broadcast;
  using Teuchos::OrdinalTraits;
  using Teuchos::ScalarTraits;
  using Teuchos::Comm;
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::tuple;
  using Teuchos::null;
  using Teuchos::VERB_NONE;
  using Teuchos::VERB_LOW;
  using Teuchos::VERB_MEDIUM;
  using Teuchos::VERB_HIGH;
  using Teuchos::VERB_EXTREME;
  using Teuchos::ETransp;
  using Teuchos::NO_TRANS;
  using Teuchos::TRANS;
  using Teuchos::CONJ_TRANS;
  using Teuchos::EDiag;
  using Teuchos::UNIT_DIAG;
  using Teuchos::NON_UNIT_DIAG;
  using Teuchos::EUplo;
  using Teuchos::UPPER_TRI;
  using Teuchos::LOWER_TRI;
  using Teuchos::ParameterList;
  using Teuchos::parameterList;

  using Tpetra::Map;
  using Tpetra::MultiVector;
  using Tpetra::Vector;
  using Tpetra::Operator;
  using Tpetra::CrsMatrix;
  using Tpetra::CrsGraph;
  using Tpetra::RowMatrix;
  using Tpetra::Import;
  using Tpetra::global_size_t;
  using Tpetra::createNonContigMapWithNode;
  using Tpetra::createUniformContigMapWithNode;
  using Tpetra::createContigMapWithNode;
  using Tpetra::createLocalMapWithNode;
  using Tpetra::createCrsMatrixSolveOp;
  using Tpetra::createVector;
  using Tpetra::createCrsMatrix;
  using Tpetra::DefaultPlatform;
  using Tpetra::ProfileType;
  using Tpetra::StaticProfile;
  using Tpetra::DynamicProfile;
  using Tpetra::OptimizeOption;
  using Tpetra::DoOptimizeStorage;
  using Tpetra::DoNotOptimizeStorage;
  using Tpetra::GloballyDistributed;
  using Tpetra::INSERT;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.addOutputSetupOptions(true);
    clp.setOption(
        "test-mpi", "test-serial", &Tpetra::TestingUtilities::testMpi,
        "Test MPI (if available) or force test of serial.  In a serial build,"
        " this option is ignored and a serial comm is always used." );
  }

  //
  // UNIT TESTS
  //

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, EmptyTriSolve, LO, GO, Scalar, Node )
  {
    RCP<Node> node = getNode<Node>();
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef Operator<Scalar,LO,GO,Node>  OP;
    typedef ScalarTraits<Scalar> ST;
    typedef MultiVector<Scalar,LO,GO,Node> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const size_t numLocal = 13, numVecs = 7;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    // get a comm
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a Map
    RCP<const Map<LO,GO,Node> > map = createContigMapWithNode<LO,GO>(INVALID,numLocal,comm,node);

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
      RCP<OP> ZeroIOp;
      {
        RCP<MAT> ZeroMat;
        // must explicitly provide the column map for implicit diagonals
        ZeroMat = rcp(new MAT(map,map,0));
        RCP<ParameterList> params = parameterList();
        RCP<ParameterList> fillparams = sublist(params,"Local Sparse Ops");
        fillparams->set("Prepare Solve", true);
        fillparams->set("Prepare Transpose Solve", true);
        fillparams->set("Prepare Conjugate Transpose Solve", true);
        ZeroMat->fillComplete(params);
        TEST_EQUALITY_CONST(ZeroMat->isUpperTriangular(), true);
        TEST_EQUALITY_CONST(ZeroMat->isLowerTriangular(), true);
        TEST_EQUALITY_CONST(ZeroMat->getGlobalNumDiags(), 0);
        ZeroIOp = createCrsMatrixSolveOp<Scalar>(ZeroMat.getConst());
      }
      X = B;
      Xhat.randomize();
      ZeroIOp->apply(B,Xhat,trans);
      //
      Xhat.update(-ST::one(),X,ST::one());
      Array<Mag> errnrms(numVecs), normsB(numVecs), zeros(numVecs, MT::zero());
      Xhat.norm1(errnrms());
      B.norm1(normsB());
      Mag maxBnrm = *std::max_element( normsB.begin(), normsB.end() );
      if (ST::isOrdinal) {
        TEST_COMPARE_ARRAYS(errnrms, zeros);
      } else {
        TEST_COMPARE_FLOATING_ARRAYS( errnrms, zeros, maxBnrm );
      }
    }
  }

  ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( CrsMatrix, TriSolve, LO, GO, Scalar, Node )
  {
    using std::endl;

    if (Teuchos::ScalarTraits<Scalar>::isOrdinal) {
      out << "Skipping testing for the integral type Scalar=" 
	  << Teuchos::TypeNameTraits<Scalar>::name () << "." << endl;
      return;
    }
    RCP<Node> node = getNode<Node>();
    typedef CrsMatrix<Scalar,LO,GO,Node> MAT;
    typedef Operator<Scalar,LO,GO,Node>  OP;
    typedef ScalarTraits<Scalar> ST;
    typedef MultiVector<Scalar,LO,GO,Node> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const size_t numLocal = 13, numVecs = 7;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid ();
    RCP<const Comm<int> > comm = getDefaultComm (); 
    // Create a row Map for the matrix.
    // This will be the same as the domain and range Maps.
    RCP<const Map<LO,GO,Node> > map = 
      createContigMapWithNode<LO, GO, Node> (INVALID, numLocal, comm, node);
    Scalar SONE = ST::one ();

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
    Array<Mag> normsX (numVecs);
    {
      X.norm1 (normsX ());
      Array<size_t> badColumns;
      for (size_t j = 0; j < numVecs; ++j) {
	if (ST::isnaninf (normsX[j])) {
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

      params->set ("Optimize Storage", optimizeStorage);
      RCP<ParameterList> fillparams = sublist (params, "Local Sparse Ops");
      fillparams->set ("Prepare Solve", true);
      fillparams->set ("Prepare Transpose Solve", true);
      fillparams->set ("Prepare Conjugate Transpose Solve", true);

      RCP<OP> AIOp;
      RCP<MAT> AMat;
      {
        if (diag == UNIT_DIAG) {
          // must explicitly specify the column map
          AMat = rcp(new MAT(map,map,2));
        }
        else {
          // can let the matrix compute a column map
          AMat = rcp(new MAT(map,2));
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
        TEST_EQUALITY(AMat->isUpperTriangular(), uplo == UPPER_TRI);
        TEST_EQUALITY(AMat->isLowerTriangular(), uplo == LOWER_TRI);
        TEST_EQUALITY(AMat->getGlobalNumDiags() == 0, diag == UNIT_DIAG);
	// AIOp.apply (X,B,trans) solves op(A) X=B for X locally,
	// using a triangular solve.  op(A) is just A if
	// trans==NO_TRANS, else A^H (Hermitian transpose) if
	// trans==CONJ_TRANS.
        AIOp = createCrsMatrixSolveOp<Scalar> (AMat.getConst ());
      }
      B.randomize ();
      if (diag != UNIT_DIAG && uplo == LOWER_TRI && 
	  trans == CONJ_TRANS && ! optimizeStorage) {
	const std::string diagStr = (diag == UNIT_DIAG) ? 
	  "UNIT_DIAG" : "NON_UNIT_DIAG";
	std::string uploStr;
	if (uplo == LOWER_TRI) {
	  uploStr = "LOWER_TRI";
	} else if (uplo == UPPER_TRI) {
	  uploStr = "UPPER_TRI";
	} else {
	  uploStr = "NEITHER";
	}
	std::string transStr;
	if (trans == CONJ_TRANS) {
	  transStr = "CONJ_TRANS";
	} else if (trans == TRANS) {
	  transStr = "TRANS";
	} else {
	  transStr = "NO_TRANS";
	}
	out << endl
	    << "================================" << endl
	    << "HERE IS THE INTERESTING USE CASE" << endl
	    << "================================" << endl
	    << "  uplo: " << uploStr << endl
	    << "  diag: " << diagStr << endl
	    << "  trans: " << transStr << endl
	    << "  optimizeStorage: " << optimizeStorage << endl
	    << endl;
	out << "Before the solve:" << endl
	    << "  Input MV X:" << endl;
	X.describe (out, Teuchos::VERB_EXTREME);
	out << "  Output MV B:" << endl;
	B.describe (out, Teuchos::VERB_EXTREME);
	out << "  Sparse matrix A:" << endl;
	AMat->describe (out, Teuchos::VERB_EXTREME);
	out << "Time for the solve!" << endl;
      }
      AIOp->apply(X,B,trans);
      if (diag == UNIT_DIAG) {
        // we want (I+A)*X -> B
        // A*X -> B needs to be augmented with X
        B.update(ST::one(),X,ST::one());
      }
      Array<Mag> normsB (numVecs);
      B.norm1 (normsB ());
      {
	Array<size_t> badColumns;
        for (size_t j = 0; j < numVecs; ++j) {
	  if (ST::isnaninf (normsB[j])) {
	    badColumns.push_back (j);
	  }
	}
	if (badColumns.size () > 0) {
	  out << endl << "*** TRIANGULAR SOLVE FAILED ***" << endl << endl;
	  B.normInf (normsB ());
	  out << "Here are the inf-norms of each column of B = A \\ X: " 
	       << Teuchos::toString (normsB) << endl;
	  B.norm2 (normsB ());
	  out << "Here are the 2-norms of each column of B = A \\ X: " 
	       << Teuchos::toString (normsB) << endl;
	  B.norm1 (normsB ());
	  out << "Here are the 1-norms of each column of B = A \\ X: " 
	       << Teuchos::toString (normsB) << endl;
	  out << "Here is the input MV X:" << endl;
	  X.describe (out, Teuchos::VERB_EXTREME);
	  out << "Here is the output MV B, on output:" << endl;
	  B.describe (out, Teuchos::VERB_EXTREME);
	}
        TEUCHOS_TEST_FOR_EXCEPTION(
          badColumns.size () > 0,
          std::runtime_error,
   	  "Columns " << Teuchos::toString (badColumns) << " of B = A \\ X "
	  "have a 1-norm either Inf or NaN." << std::endl 
	  << "That suggests the triangular solve is broken." << endl
	  << "Here are the 1-norms of each column: " 
	  << Teuchos::toString (normsB));
      }

      Xhat.randomize();
      AIOp->apply(B,Xhat,trans);
      //
      Xhat.update(-ST::one(),X,ST::one());
      Array<Mag> errnrms(numVecs), zeros(numVecs, MT::zero());
      Xhat.norm1(errnrms());
      B.norm1(normsB());
      Mag maxBnrm = *std::max_element( normsB.begin(), normsB.end() );
      if (ST::isOrdinal) {
        TEST_COMPARE_ARRAYS(errnrms, zeros);
      } else {
        TEST_COMPARE_FLOATING_ARRAYS( errnrms, zeros, maxBnrm );
      }
    }
  }

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, EmptyTriSolve, LO, GO, SCALAR, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, TriSolve, LO, GO, SCALAR, NODE )

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_SLGN( UNIT_TEST_GROUP )

}
