// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <string>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_as.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_TestingHelpers.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_ParameterXMLFileReader.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <MatrixMarket_Tpetra.hpp> // for loading matrices from file

#include "Amesos2.hpp"
#include "Amesos2_Meta.hpp"


namespace {

  using std::cout;
  using std::endl;
  using std::string;

  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::ArrayRCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  using Teuchos::Comm;
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::tuple;
  using Teuchos::ScalarTraits;
  using Teuchos::OrdinalTraits;
  using Teuchos::FancyOStream;
  using Teuchos::VerboseObjectBase;
  using Teuchos::ETransp;
  using Teuchos::EUplo;
  using Teuchos::LOWER_TRI;
  using Teuchos::UPPER_TRI;
  using Teuchos::CONJ_TRANS;
  using Teuchos::TRANS;
  using Teuchos::NO_TRANS;


  using Tpetra::global_size_t;
  using Tpetra::CrsMatrix;
  using Tpetra::MultiVector;
  using Tpetra::Map;
  using Tpetra::createContigMap;
  using Tpetra::createUniformContigMap;

  // using Amesos2::MatrixAdapter;
  // using Amesos2::MultiVecAdapter;
  using Amesos2::PardisoMKL;

  typedef Tpetra::Map<>::node_type Node;

  bool testMpi = true;

  // Where to look for input files
  string filedir;

  RCP<const Comm<int> > getDefaultComm()
  {
    RCP<const Comm<int> > ret;
    if( testMpi ){
      ret = Tpetra::getDefaultComm();
    } else {
      ret = rcp(new Teuchos::SerialComm<int>());
    }
    return ret;
  }

  //RCP<FancyOStream> getDefaultOStream()
  //{
  //  return( VerboseObjectBase::getDefaultOStream() );
  //}

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.setOption("filedir",&filedir,"Directory of matrix files.");
    clp.addOutputSetupOptions(true);
    clp.setOption("test-mpi", "test-serial", &testMpi,
                  "Test MPI by default or force serial test.  In a serial build,"
                  " this option is ignored and a serial comm is always used." );
  }


  /*
   * UNIT TESTS
   */

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( PardisoMKL, Partial1, SCALAR, LO, GO )
  {
    typedef CrsMatrix<SCALAR,LO,GO,Node> MAT;
    typedef ScalarTraits<SCALAR> ST;
    typedef MultiVector<SCALAR,LO,GO,Node> MV;
    typedef typename ST::magnitudeType Mag;
    const size_t numVecs = 1;

    for (int test = 0; test < 2; test++) {
      // test-0 with default-comm, matrix will be gather to mpi-0 and only mpi-0 computes Schur complement
      // test-1 with serial-comm, all processes calls pardiso and computes schur complement independently
      RCP<const Comm<int> > comm = (test == 0 ? Tpetra::getDefaultComm() :
                                                Teuchos::rcp_dynamic_cast<const Comm<int>>(rcp(new Teuchos::SerialComm<int>())));
      const size_t myRank = comm->getRank();
      const size_t numRanks = comm->getSize();

      RCP<const Comm<int> > global_comm = Tpetra::getDefaultComm(); // just for printing
      if (global_comm->getRank()==0) {
        std::cout << std::endl
                  << " >> UnitTest for PardisoMKL:Partial1 with Scalar = " << ST::name()
                  << " and " << (test == 0 ? "DefaultComm" : "SerialComm") << " <<" << std::endl << std::endl;
      }

      // Construct matrix
      // NOTE: with serial-comm, every process construct its own 10-by-10 matrix, independently
      //       with global-comm, all the process construct a (10*numRanks)-by-(10*numRanks) matrix, jointly
      // Schur complement will half of the global matrix
      const SCALAR one  = ST::one();
      const SCALAR mone = -one;
      const SCALAR two  = one + one;
      const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
      const size_t numLocal = 10;
      const size_t numGlobal = numLocal*numRanks;
      RCP<Map<LO,GO,Node> > map = rcp( new Map<LO,GO,Node>(INVALID,numLocal,0,comm) );
      RCP<MAT> A = rcp( new MAT(map, 3) );
      GO base = numLocal*myRank;
      for( size_t i = 0; i < numLocal; ++i ){
        if (base+i > 0) {
          A->insertGlobalValues(base+i,tuple<GO>(base+i-1),tuple<SCALAR>(mone));
        }
        A->insertGlobalValues(base+i,tuple<GO>(base+i),tuple<SCALAR>(two));
        if (base+i < numGlobal-1) {
          A->insertGlobalValues(base+i,tuple<GO>(base+i+1),tuple<SCALAR>(mone));
        }
      }
      A->fillComplete();

      RCP<const Map<LO,GO,Node> > dmnmap = A->getDomainMap();
      RCP<const Map<LO,GO,Node> > rngmap = A->getRangeMap();

      RCP<MV> X = rcp(new MV(dmnmap,numVecs));
      RCP<MV> B = rcp(new MV(rngmap,numVecs));
      RCP<MV> Xhat = rcp(new MV(dmnmap,numVecs));
      X->setObjectLabel("X");
      B->setObjectLabel("B");
      Xhat->setObjectLabel("Xhat");
      X->randomize();

      A->apply(*X,*B);            // no transpose

      Xhat->randomize();

      {
        // Create PardisoMKL solver
        RCP<Amesos2::Solver<MAT,MV> > solver
          = Amesos2::create<MAT,MV>("PARDISOMKL", A, Xhat, B );

        // Parameters
        Teuchos::ParameterList amesos2_paramlist;
        amesos2_paramlist.setName("Amesos2");
        Teuchos::ParameterList & pardiso_paramlist = amesos2_paramlist.sublist("PARDISOMKL");
        // partial factorization currently requires at least two threads
        pardiso_paramlist.set("PartialFacto", 1, "Partial Factorization");
        // Schur part has odd row IDs
        Teuchos::Array<LO> schurPart(numGlobal);
        for( size_t i = 0; i < numGlobal; i++) {
          if (i%2 == 0) schurPart[i] = 0;
          if (i%2 == 1) schurPart[i] = 1;
        }
        pardiso_paramlist.set("SchurPart", (const LO*)schurPart.getRawPtr());
        pardiso_paramlist.set("MessageLevel", (global_comm->getRank() == 0 ? 1 : 0));
        solver->setParameters(Teuchos::rcpFromRef(amesos2_paramlist));

        // Solve A*Xhat = B for Xhat using the PardisoMKL solver
        solver->symbolicFactorization();
        solver->numericFactorization();
        solver->solve();

        // Check result of solve
        Array<Mag> xhatnorms(numVecs), xnorms(numVecs);
        Xhat->norm2(xhatnorms());
        X->norm2(xnorms());
        if (global_comm->getRank()==0) {
          for (int i=0; i<xnorms.size(); i++)
            std::cout << "err[" << i << "]  = " << xnorms[i] << " - " << xhatnorms[i]
                      << " = " << xnorms[i]-xhatnorms[i] << std::endl;
        }
        TEST_COMPARE_FLOATING_ARRAYS( xhatnorms, xnorms, 0.005 );
      }
    }
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( PardisoMKL, Partial2, SCALAR, LO, GO )
  {
    typedef ScalarTraits<SCALAR> ST;
    typedef CrsMatrix<SCALAR,LO,GO,Node> MAT;
    typedef MultiVector<SCALAR,LO,GO,Node> MV;

    for (int test = 0; test < 2; test++) {
      // test-0 with default-comm, matrix will be gather to mpi-0 and only mpi-0 computes Schur complement
      // test-1 with serial-comm, all processes calls pardiso and computes schur complement independently
      RCP<const Comm<int> > comm = (test == 0 ? Tpetra::getDefaultComm() :
                                                Teuchos::rcp_dynamic_cast<const Comm<int>>(rcp(new Teuchos::SerialComm<int>())));
      const size_t myRank = comm->getRank();

      RCP<const Comm<int> > global_comm = Tpetra::getDefaultComm(); // just for printing
      if (global_comm->getRank()==0) {
        std::cout << std::endl
                  << " >> UnitTest for PardisoMKL::Partial2 with Scalar = " << ST::name()
                  << " and " << (test == 0 ? "DefaultComm" : "SerialComm") << " <<" << std::endl << std::endl;
      }
      const size_t numRanks = comm->getSize();
      const size_t numVecs = 1;

      // Construct matrix
      // NOTE: with serial-comm, every process construct its own 10-by-10 matrix, independently
      //       with global-comm, all the process construct a (10*numRanks)-by-(10*numRanks) matrix, jointly
      // Schur complement will half of the global matrix
      const SCALAR one  = ST::one();
      const SCALAR mone = -one;
      const SCALAR two  = one + one;
      const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
      const size_t numLocal = 10;
      const size_t numGlobal = numLocal*numRanks;
      RCP<Map<LO,GO,Node> > map = rcp( new Map<LO,GO,Node>(INVALID,numLocal,0,comm) );
      RCP<MAT> A = rcp( new MAT(map, 3) );
      GO base = numLocal*myRank;
      for( size_t i = 0; i < numLocal; ++i ){
        if (base+i > 0) {
          A->insertGlobalValues(base+i,tuple<GO>(base+i-1),tuple<SCALAR>(mone));
        }
        A->insertGlobalValues(base+i,tuple<GO>(base+i),tuple<SCALAR>(two));
        if (base+i < numGlobal-1) {
          A->insertGlobalValues(base+i,tuple<GO>(base+i+1),tuple<SCALAR>(mone));
        }
      }
      A->fillComplete();

      {
        // Create PardisoMKL solver
        RCP<Amesos2::Solver<MAT,MV> > solver
          = Amesos2::create<MAT,MV>("PARDISOMKL", A );

        // Parameters
        Teuchos::ParameterList amesos2_paramlist;
        amesos2_paramlist.setName("Amesos2");
        Teuchos::ParameterList & pardiso_paramlist = amesos2_paramlist.sublist("PARDISOMKL");
        pardiso_paramlist.set("PartialFacto", 2, "Partial Factorization");
        // Schur part has odd row IDs
        Teuchos::Array<LO> schurPart(numGlobal);
        for( size_t i = 0; i < numGlobal; i++) {
          if (i%2 == 0) schurPart[i] = 0;
          if (i%2 == 1) schurPart[i] = 1;
        }
        const size_t numSchur = numGlobal/2;
        Teuchos::Array<SCALAR> schurOut(numSchur*numSchur);
        pardiso_paramlist.set("SchurPart", (const LO*)schurPart.getRawPtr());
        pardiso_paramlist.set("SchurOut", (SCALAR*)schurOut.getRawPtr());
        pardiso_paramlist.set("MessageLevel", (global_comm->getRank() == 0 ? 1 : 0));
        solver->setParameters(Teuchos::rcpFromRef(amesos2_paramlist));

        // Perform Partial Facto to compute Schur complement
        solver->symbolicFactorization();
        solver->numericFactorization();

        // Check Schur complement
        const SCALAR zero = ST::zero();
        const SCALAR half  = one / two;
        Teuchos::Array<SCALAR> schur(numSchur*numSchur, zero);
        // MPI-0 on comm (either serial or default) will compute Schur complement
        // NOTE: with serial comm, every one has myRank == 0, and compute Schur complement independently
        if (myRank == 0) {
          for( size_t i = 0; i < numSchur; ++i ){
            if (i > 0) {
              schur[i-1 + i*numSchur] = -half;
            }
            if (i == numSchur-1) {
              schur[i + i*numSchur] = one+half;
            } else {
              schur[i + i*numSchur] = one;
            }
            if (base+i < numSchur-1) {
              schur[i+1 + i*numSchur] = -half;
            }
          }
        }
        if (global_comm->getRank() == 0) {
          std::cout << "[" << std::endl;
          for( size_t i = 0; i < numSchur; ++i ){
            std::cout << "  ";
            for( size_t j = 0; j < numSchur; ++j ) std::cout << schurOut[i+j*numSchur] << " ";
            std::cout << std::endl;
          }
          std::cout << "];" << std::endl;
          std::cout << "[" << std::endl;
          for( size_t i = 0; i < numSchur; ++i ){
            std::cout << "  ";
            for( size_t j = 0; j < numSchur; ++j ) std::cout << schur[i+j*numSchur] << " ";
            std::cout << std::endl;
          }
          std::cout << "];" << std::endl;
        }
        TEST_COMPARE_FLOATING_ARRAYS( schurOut, schur, 0.005 );
      }
    }
  }

  /*
   * Instantiations
   */

#ifdef HAVE_TPETRA_INST_FLOAT
#  define UNIT_TEST_GROUP_ORDINAL_FLOAT( LO, GO )       \
  UNIT_TEST_GROUP_ORDINAL_SCALAR( LO, GO, float )
#else
#  define UNIT_TEST_GROUP_ORDINAL_FLOAT( LO, GO )
#endif
#ifdef HAVE_TPETRA_INST_DOUBLE
#  define UNIT_TEST_GROUP_ORDINAL_DOUBLE( LO, GO )      \
  UNIT_TEST_GROUP_ORDINAL_SCALAR( LO, GO, double )
#else
#  define UNIT_TEST_GROUP_ORDINAL_DOUBLE( LO, GO )
#endif

#define UNIT_TEST_GROUP_ORDINAL_SCALAR( LO, GO, SCALAR )                \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( PardisoMKL, Partial1, SCALAR, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( PardisoMKL, Partial2, SCALAR, LO, GO )

#define UNIT_TEST_GROUP_ORDINAL_ORDINAL( LO, GO )     \
  UNIT_TEST_GROUP_ORDINAL_FLOAT(LO, GO)               \
  UNIT_TEST_GROUP_ORDINAL_DOUBLE(LO, GO)              \

#define UNIT_TEST_GROUP_ORDINAL( ORDINAL )            \
  UNIT_TEST_GROUP_ORDINAL_ORDINAL( ORDINAL, ORDINAL )

  //Add JDB (10-19-215)
#ifndef HAVE_AMESOS2_EXPLICIT_INSTANTIATION
  UNIT_TEST_GROUP_ORDINAL(int)
  typedef long int LongInt;
  UNIT_TEST_GROUP_ORDINAL_ORDINAL( int, LongInt )
  #ifdef HAVE_TPETRA_INT_LONG_LONG
  typedef long long int LongLongInt;
  UNIT_TEST_GROUP_ORDINAL_ORDINAL( int, LongLongInt )
  #endif
#else  //ETI
  #ifdef HAVE_TPETRA_INST_INT_INT
  UNIT_TEST_GROUP_ORDINAL(int)
  #endif
  #ifdef HAVE_TPETRA_INST_INT_LONG
  typedef long int LongInt;
  UNIT_TEST_GROUP_ORDINAL_ORDINAL(int,LongInt)
  #endif
  #ifdef HAVE_TPETRA_INST_INT_LONG_LONG
  typedef long long int LongLongInt;
  UNIT_TEST_GROUP_ORDINAL_ORDINAL(int,LongLongInt)
  #endif
#endif  // EXPL-INST

} // end anonymous namespace
