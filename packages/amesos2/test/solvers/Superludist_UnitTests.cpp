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

#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <MatrixMarket_Tpetra.hpp> // for loading matrices from file

#include "Amesos2.hpp"
#include "Amesos2_Meta.hpp"

#if defined(HAVE_AMESOS2_EPETRA)
#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif
#include <Epetra_Map.h>
#include <Epetra_MultiVector.h>
#include <Epetra_CrsMatrix.h>
#endif

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

  using Tpetra::global_size_t;
  using Tpetra::CrsMatrix;
  using Tpetra::MultiVector;
  using Tpetra::Map;
  using Tpetra::createContigMap;
  using Tpetra::createUniformContigMap;

  using Amesos2::Superludist;

  typedef Tpetra::Map<>::node_type Node;

  bool testMpi = true;

  // Where to look for input files
  string filedir;

  TEUCHOS_STATIC_SETUP()
  {
    Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
    clp.setOption("filedir",&filedir,"Directory of matrix files.");
    clp.addOutputSetupOptions(true);
    clp.setOption("test-mpi", "test-serial", &testMpi,
                  "Test MPI by default or force serial test.  In a serial build,"
                  " this option is ignored and a serial comm is always used." );
  }

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


  /*
   * UNIT TESTS
   */

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( Superludist, Initialization, SCALAR, LO, GO )
  {
    /* Test correct initialization of the Solver instance
     *
     * - All Constructors
     * - Correct initialization of class members
     * - Correct typedefs
     */
    typedef ScalarTraits<SCALAR> ST;
    typedef CrsMatrix<SCALAR,LO,GO,Node> MAT;
    typedef MultiVector<SCALAR,LO,GO,Node> MV;
    typedef Superludist<MAT,MV> SOLVER;

    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    RCP<const Comm<int> > comm = getDefaultComm();
    //const size_t numprocs = comm->getSize();
    const size_t rank     = comm->getRank();
    // create a Map
    const size_t numLocal = 10;
    RCP<Map<LO,GO,Node> > map = rcp( new Map<LO,GO,Node>(INVALID,numLocal,0,comm) );
    RCP<MAT> eye = rcp( new MAT(map,1) );
    GO base = numLocal*rank;
    for( size_t i = 0; i < numLocal; ++i ){
      eye->insertGlobalValues(base+i,tuple<GO>(base+i),tuple<SCALAR>(ST::one()));
    }
    eye->fillComplete();

    // Create X
    RCP<MV> X = rcp(new MV(map,11));
    X->randomize();

    // Create B
    RCP<MV> B = rcp(new MV(map,11));
    B->randomize();

    // Constructor from Factory.  Every overloaded Factory method
    // eventually calls the same/only Solver constructor
    RCP<Amesos2::Solver<MAT,MV> > solver = Amesos2::create<MAT,MV>("SuperLU_DIST",eye,X,B);

    TEST_ASSERT( solver->getStatus().getNumSymbolicFact() == 0 );
    TEST_ASSERT( solver->getStatus().getNumNumericFact() == 0 );
    TEST_ASSERT( solver->getStatus().getNumSolve() == 0 );

    // The following should all pass at compile time
    TEST_ASSERT( (std::is_same_v<MAT,typename SOLVER::matrix_type>) );
    TEST_ASSERT( (std::is_same_v<MV,typename SOLVER::vector_type>) );
    TEST_ASSERT( (std::is_same_v<SCALAR,typename SOLVER::scalar_type>) );
    TEST_ASSERT( (std::is_same_v<LO,typename SOLVER::local_ordinal_type>) );
    TEST_ASSERT( (std::is_same_v<GO,typename SOLVER::global_ordinal_type>) );
    TEST_ASSERT( (std::is_same_v<global_size_t,typename SOLVER::global_size_type>) );
    // TEST_ASSERT( (std::is_same_v<Node,typename SOLVER::node_type>) );
  }


  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( Superludist, SymbolicFactorization, SCALAR, LO, GO )
  {
    typedef ScalarTraits<SCALAR> ST;
    typedef CrsMatrix<SCALAR,LO,GO,Node> MAT;
    typedef MultiVector<SCALAR,LO,GO,Node> MV;

    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    RCP<const Comm<int> > comm = getDefaultComm();
    //const size_t numprocs = comm->getSize();
    const size_t rank     = comm->getRank();
    // create a Map
    const size_t numLocal = 10;
    RCP<Map<LO,GO,Node> > map = rcp( new Map<LO,GO,Node>(INVALID,numLocal,0,comm) );
    RCP<MAT> eye = rcp( new MAT(map,1) );
    GO base = numLocal*rank;
    for( size_t i = 0; i < numLocal; ++i ){
      eye->insertGlobalValues(base+i,tuple<GO>(base+i),tuple<SCALAR>(ST::one()));
    }
    eye->fillComplete();

    // Create X
    RCP<MV> X = rcp(new MV(map,11));
    X->randomize();

    // Create B
    RCP<MV> B = rcp(new MV(map,11));
    B->randomize();

    // Constructor from Factory
    RCP<Amesos2::Solver<MAT,MV> > solver = Amesos2::create<MAT,MV>("SuperLU_DIST",eye,X,B);

    solver->symbolicFactorization();
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( Superludist, NumericFactorization, SCALAR, LO, GO )
  {
    typedef ScalarTraits<SCALAR> ST;
    typedef CrsMatrix<SCALAR,LO,GO,Node> MAT;
    typedef MultiVector<SCALAR,LO,GO,Node> MV;

    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    RCP<const Comm<int> > comm = getDefaultComm();
    //const size_t numprocs = comm->getSize();
    const size_t rank     = comm->getRank();
    // create a Map
    const size_t numLocal = 10;
    RCP<Map<LO,GO,Node> > map = rcp( new Map<LO,GO,Node>(INVALID,numLocal,0,comm) );
    RCP<MAT> eye = rcp( new MAT(map,1) );
    GO base = numLocal*rank;
    for( size_t i = 0; i < numLocal; ++i ){
      eye->insertGlobalValues(base+i,tuple<GO>(base+i),tuple<SCALAR>(ST::one()));
    }
    eye->fillComplete();

    // Create X
    RCP<MV> X = rcp(new MV(map,11));
    X->randomize();

    // Create B
    RCP<MV> B = rcp(new MV(map,11));
    B->randomize();

    // Constructor from Factory
    RCP<Amesos2::Solver<MAT,MV> > solver = Amesos2::create<MAT,MV>("SuperLU_DIST",eye,X,B);

    solver->symbolicFactorization().numericFactorization();

    // Good way to check the factors L and U?  Probs not, since they are private members
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( Superludist, NonContgGID, SCALAR, LO, GO )
  {
    typedef CrsMatrix<SCALAR,LO,GO,Node> MAT;
    typedef ScalarTraits<SCALAR> ST;
    typedef MultiVector<SCALAR,LO,GO,Node> MV;
    typedef typename ST::magnitudeType Mag;

    using Tpetra::global_size_t;
    using Teuchos::tuple;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Scalar = SCALAR;

    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();

    size_t myRank = comm->getRank();
    const global_size_t numProcs = comm->getSize();

    // Unit test created for 2 processes
    if ( numProcs == 2 ) {

      const global_size_t numVectors = 1;
      const global_size_t nrows = 6;

      const GO numGlobalEntries = nrows;
      const LO numLocalEntries = nrows / numProcs;

      // Create non-contiguous Map
      // This example: np 2 leads to GIDS: proc0 - 0,2,4  proc 1 - 1,3,5
      Teuchos::Array<GO> elementList(numLocalEntries);
      for ( LO k = 0; k < numLocalEntries; ++k ) {
        elementList[k] = myRank + k*numProcs + 4*myRank;
      }

      typedef Tpetra::Map<LO,GO>  map_type;
      RCP< const map_type > map = rcp( new map_type(numGlobalEntries, elementList, 0, comm) );
      TEUCHOS_TEST_FOR_EXCEPTION(
          comm->getSize () > 1 && map->isContiguous (),
          std::logic_error,
          "SuperLU_DIST NonContigGID Test: The non-contiguous map claims to be contiguous.");

      RCP<MAT> A = rcp( new MAT(map,3) ); // max of three entries in a row
      A->setObjectLabel("A");

      /*
       * We will solve a system with a known solution, for which we will be using
       * the following matrix:
       *
       *  GID  0   2   4   5   7   9
       * [ 0 [ 7,  0, -3,  0, -1,  0 ]
       *   2 [ 2,  8,  0,  0,  0,  0 ]
       *   4 [ 0,  0,  1,  0,  0,  0 ]
       *   5 [-3,  0,  0,  5,  0,  0 ]
       *   7 [ 0, -1,  0,  0,  4,  0 ]
       *   9 [ 0,  0,  0, -2,  0,  6 ] ]
       *
       */

      // Construct matrix
      if( myRank == 0 ){
        A->insertGlobalValues(0,tuple<GO>(0,4,7),tuple<Scalar>(7,-3,-1));
        A->insertGlobalValues(2,tuple<GO>(0,2),tuple<Scalar>(2,8));
        A->insertGlobalValues(4,tuple<GO>(4),tuple<Scalar>(1));
        A->insertGlobalValues(5,tuple<GO>(0,5),tuple<Scalar>(-3,5));
        A->insertGlobalValues(7,tuple<GO>(2,7),tuple<Scalar>(-1,4));
        A->insertGlobalValues(9,tuple<GO>(5,9),tuple<Scalar>(-2,6));
      }

      A->fillComplete();

      TEUCHOS_TEST_FOR_EXCEPTION(
          comm->getSize () > 1 && A->getMap()->isContiguous (),
          std::logic_error,
          "SuperLU_DIST NonContigGID Test: The non-contiguous map of A claims to be contiguous.");


      // Create X with random values
      RCP<MV> X = rcp(new MV(map,numVectors));
      X->setObjectLabel("X");
      X->randomize();

      /* Create B, use same GIDs
       *
       * Use RHS:
       *
       *  [[-7]
       *   [18]
       *   [ 3]
       *   [17]
       *   [18]
       *   [28]]
       */
      RCP<MV> B = rcp(new MV(map,numVectors));
      B->setObjectLabel("B");
      GO rowids[nrows] = {0,2,4,5,7,9};
      Scalar data[nrows] = {-7,18,3,17,18,28};
      for( global_size_t i = 0; i < nrows; ++i ){
        if( B->getMap()->isNodeGlobalElement(rowids[i]) ){
          B->replaceGlobalValue(rowids[i],0,data[i]);
        }
      }

      TEUCHOS_TEST_FOR_EXCEPTION(
          comm->getSize () > 1 && X->getMap()->isContiguous () && B->getMap()->isContiguous (),
          std::logic_error,
          "SuperLU_DIST NonContigGID Test: The non-contiguous maps of X or B claims to be contiguous.");


      // Create solver interface with Amesos2 factory method
      RCP<Amesos2::Solver<MAT,MV> > solver = Amesos2::create<MAT,MV>("SuperLU_DIST", A, X, B);

      // Create a Teuchos::ParameterList to hold solver parameters
      Teuchos::ParameterList amesos2_params("Amesos2");
      amesos2_params.sublist("SuperLU_DIST").set("IsContiguous", false, "Are GIDs Contiguous");

      solver->setParameters( Teuchos::rcpFromRef(amesos2_params) );

      solver->symbolicFactorization().numericFactorization().solve();


      /* Check the solution
       *
       * Should be:
       *
       *  [[1]
       *   [2]
       *   [3]
       *   [4]
       *   [5]
       *   [6]]
       */
      // Solution Vector for later comparison
      RCP<MV> Xhat = rcp(new MV(map,numVectors));
      Xhat->setObjectLabel("Xhat");
      GO rowids_soln[nrows] = {0,2,4,5,7,9};
      Scalar data_soln[nrows] = {1,2,3,4,5,6};
      for( global_size_t i = 0; i < nrows; ++i ){
        if( Xhat->getMap()->isNodeGlobalElement(rowids_soln[i]) ){
          Xhat->replaceGlobalValue(rowids_soln[i],0,data_soln[i]);
        }
      }

      // Check result of solve
      Array<Mag> xhatnorms(numVectors), xnorms(numVectors);
      Xhat->norm2(xhatnorms());
      X->norm2(xnorms());
      TEST_COMPARE_FLOATING_ARRAYS( xhatnorms, xnorms, 0.005 );
    } // end if numProcs = 2
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( Superludist, ContgGID, SCALAR, LO, GO )
  {
    typedef CrsMatrix<SCALAR,LO,GO,Node> MAT;
    typedef ScalarTraits<SCALAR> ST;
    typedef MultiVector<SCALAR,LO,GO,Node> MV;
    typedef typename ST::magnitudeType Mag;

    using Tpetra::global_size_t;
    using Teuchos::tuple;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Scalar = SCALAR;

    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();

    size_t myRank = comm->getRank();
    const global_size_t numProcs = comm->getSize();

    printf( " Contig test\n" );
    // Unit test created for 2 processes
    if ( numProcs == 2 ) {

      const global_size_t numVectors = 1;
      const global_size_t nrows = 6;

      const GO numGlobalEntries = nrows;
      const LO numLocalEntries = (myRank == 0 ? nrows : 0);

      // Create non-contiguous Map
      // This example: np 2 leads to GIDS: proc0 - 0,2,4  proc 1 - 1,3,5
      Teuchos::Array<GO> elementList(numLocalEntries);
      for ( LO k = 0; k < numLocalEntries; ++k ) {
        elementList[k] = k;
      }

      typedef Tpetra::Map<LO,GO>  map_type;
      RCP< const map_type > map = rcp( new map_type(numGlobalEntries, elementList, 0, comm) );

      RCP<MAT> A = rcp( new MAT(map,3) ); // max of three entries in a row
      A->setObjectLabel("A");

      /*
       * We will solve a system with a known solution, for which we will be using
       * the following 2D Laplace matrix:
       */

      // Construct matrix
      if( myRank == 0 ){
        A->insertGlobalValues(0,tuple<GO>(  0,1),tuple<Scalar>(2,-1));
        A->insertGlobalValues(1,tuple<GO>(0,1,2),tuple<Scalar>(-1,2,-1));
        A->insertGlobalValues(2,tuple<GO>(1,2,3),tuple<Scalar>(-1,2,-1));
        A->insertGlobalValues(3,tuple<GO>(2,3,4),tuple<Scalar>(-1,2,-1));
        A->insertGlobalValues(4,tuple<GO>(3,4,5),tuple<Scalar>(-1,2,-1));
        A->insertGlobalValues(5,tuple<GO>(4,5),  tuple<Scalar>(-1,2));
      }
      A->fillComplete();
      //{ RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(rcpFromRef(std::cout)); A->describe(*fancy,Teuchos::VERB_EXTREME); }

      // Create X with random values
      RCP<MV> X = rcp(new MV(map,numVectors));
      X->setObjectLabel("X");
      X->randomize();
      //{ RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(rcpFromRef(std::cout)); X->describe(*fancy,Teuchos::VERB_EXTREME); }

      /* Create B, use same GIDs
       *
       * Use RHS:
       *
       *  [[1]
       *   [0]
       *   [0]
       *   [0]
       *   [0]
       *   [1]]
       */
      RCP<MV> B = rcp(new MV(map,numVectors));
      B->setObjectLabel("B");
      GO rowids[nrows] = {0,1,2,3,4,5};
      Scalar data[nrows] = {1,0,0,0,0,1};
      for( global_size_t i = 0; i < nrows; ++i ){
        if( B->getMap()->isNodeGlobalElement(rowids[i]) ){
          B->replaceGlobalValue(rowids[i],0,data[i]);
        }
      }

      // Create solver interface with Amesos2 factory method
      RCP<Amesos2::Solver<MAT,MV> > solver = Amesos2::create<MAT,MV>("SuperLU_DIST", A, X, B);

      // Create a Teuchos::ParameterList to hold solver parameters
      solver->symbolicFactorization().numericFactorization().solve();
      //{ RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(rcpFromRef(std::cout)); B->describe(*fancy,Teuchos::VERB_EXTREME); }
      //{ RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(rcpFromRef(std::cout)); X->describe(*fancy,Teuchos::VERB_EXTREME); }


      /* Check the solution
       *
       * Should be:
       *
       *  [[1]
       *   [1]
       *   [1]
       *   [1]
       *   [1]
       *   [1]]
       */
      // Solution Vector for later comparison
      RCP<MV> Xhat = rcp(new MV(map,numVectors));
      Xhat->setObjectLabel("Xhat");
      GO rowids_soln[nrows] = {0,1,2,3,4,5};
      Scalar data_soln[nrows] = {1,1,1,1,1,1};
      for( global_size_t i = 0; i < nrows; ++i ){
        if( Xhat->getMap()->isNodeGlobalElement(rowids_soln[i]) ){
          Xhat->replaceGlobalValue(rowids_soln[i],0,data_soln[i]);
        }
      }

      // Check result of solve
      Array<Mag> xhatnorms(numVectors), xnorms(numVectors);
      Xhat->norm2(xhatnorms());
      X->norm2(xnorms());
      TEST_COMPARE_FLOATING_ARRAYS( xhatnorms, xnorms, 0.005 );
    } // end if numProcs = 2
  }

#if defined(HAVE_AMESOS2_EPETRA)
  //! @test Test for Epetra non-contiguous GIDs.
  TEUCHOS_UNIT_TEST( Superludist, NonContigGIDEpetra)
  {
    using SCALAR = double;
    using LO = int;
    using GO = int;
    typedef Epetra_CrsMatrix MAT;
    typedef ScalarTraits<SCALAR> ST;
    typedef Epetra_MultiVector MV;
    typedef typename ST::magnitudeType Mag;
    typedef Epetra_Map map_type;

    using Tpetra::global_size_t;
    using Teuchos::tuple;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Scalar = SCALAR;

    RCP<const Comm<int> > t_comm = Tpetra::getDefaultComm();
    size_t myRank = t_comm->getRank();
    const global_size_t numProcs = t_comm->getSize();

    // Unit test created for 2 processes
    if ( numProcs == 2 ) {

      const global_size_t numVectors = 1;
      const global_size_t nrows = 6;

      const int numGlobalEntries = nrows;
      const LO numLocalEntries = nrows / numProcs;

      // Create non-contiguous Map
      // This example: np 2 leads to GIDS: proc0 - 0,2,4  proc 1 - 1,3,5
      Teuchos::Array<int> elementList(numLocalEntries);
      for ( LO k = 0; k < numLocalEntries; ++k ) {
        elementList[k] = myRank + k*numProcs + 4*myRank;
      }

#ifdef HAVE_MPI
      const Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
      const Epetra_SerialComm comm;
#endif

      RCP< const map_type > map = rcp( new Epetra_Map( numGlobalEntries, numLocalEntries, elementList.data(), 0, comm ));
      RCP<MAT> A = rcp( new MAT(Epetra_DataAccess::Copy, *map, numLocalEntries) );

      /*
       * We will solve a system with a known solution, for which we will be using
       * the following matrix:
       *
       *  GID  0   2   4   5   7   9
       * [ 0 [ 7,  0, -3,  0, -1,  0 ]
       *   2 [ 2,  8,  0,  0,  0,  0 ]
       *   4 [ 0,  0,  1,  0,  0,  0 ]
       *   5 [-3,  0,  0,  5,  0,  0 ]
       *   7 [ 0, -1,  0,  0,  4,  0 ]
       *   9 [ 0,  0,  0, -2,  0,  6 ] ]
       *
       */

      // Construct matrix
      if( myRank == 0 ){
        A->InsertGlobalValues(0,3,tuple<Scalar>(7,-3,-1).data(),tuple<GO>(0,4,7).data());
        A->InsertGlobalValues(2,2,tuple<Scalar>(2,8).data(),tuple<GO>(0,2).data());
        A->InsertGlobalValues(4,1,tuple<Scalar>(1).data(),tuple<GO>(4).data());
      } else {
        A->InsertGlobalValues(5,2,tuple<Scalar>(-3,5).data(),tuple<GO>(0,5).data());
        A->InsertGlobalValues(7,2,tuple<Scalar>(-1,4).data(),tuple<GO>(2,7).data());
        A->InsertGlobalValues(9,2,tuple<Scalar>(-2,6).data(),tuple<GO>(5,9).data());
      }
      A->FillComplete();

      // Create X with random values
      RCP<MV> X = rcp(new MV(*map, numVectors));
      X->Random();

      /* Create B, use same GIDs
       *
       * Use RHS:
       *
       *  [[-7]
       *   [18]
       *   [ 3]
       *   [17]
       *   [18]
       *   [28]]
       */
      RCP<MV> B = rcp(new MV(*map, numVectors));
      GO rowids[nrows] = {0,2,4,5,7,9};
      Scalar data[nrows] = {-7,18,3,17,18,28};
      for( global_size_t i = 0; i < nrows; ++i ){
        if( B->Map().MyGID(rowids[i]) ){
          B->ReplaceGlobalValue(rowids[i],0,data[i]);
        }
      }

      // Create solver interface with Amesos2 factory method
      RCP<Amesos2::Solver<MAT,MV> > solver = Amesos2::create<MAT,MV>("SuperLU_DIST", A, X, B);

      // Create a Teuchos::ParameterList to hold solver parameters
      Teuchos::ParameterList amesos2_params("Amesos2");
      amesos2_params.sublist("SuperLU_DIST").set("IsContiguous", false, "Are GIDs Contiguous");

      solver->setParameters( Teuchos::rcpFromRef(amesos2_params) );

      solver->symbolicFactorization().numericFactorization().solve();

      /* Check the solution
       *
       * Should be:
       *
       *  [[1]
       *   [2]
       *   [3]
       *   [4]
       *   [5]
       *   [6]]
       */
      // Solution Vector for later comparison
      RCP<MV> Xhat = rcp(new MV(*map, numVectors));
      GO rowids_soln[nrows] = {0,2,4,5,7,9};
      Scalar data_soln[nrows] = {1,2,3,4,5,6};
      for( global_size_t i = 0; i < nrows; ++i ){
        if( Xhat->Map().MyGID(rowids_soln[i]) ){
          Xhat->ReplaceGlobalValue(rowids_soln[i],0,data_soln[i]);
        }
      }

      // Check result of solve
      Array<Mag> xhatnorms(numVectors), xnorms(numVectors);
      Xhat->Norm2(xhatnorms().data());
      X->Norm2(xnorms().data());
      TEST_COMPARE_FLOATING_ARRAYS( xhatnorms, xnorms, 0.005 );
    } // end if numProcs = 2
  }
#endif


  /*
   * Instantiations
   */
#ifdef HAVE_TEUCHOS_COMPLEX
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_SCALAR(LO, GO, SCALAR)        \
  typedef std::complex<SCALAR>  Complex##SCALAR;                        \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Superludist, Initialization, Complex##SCALAR, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Superludist, SymbolicFactorization, Complex##SCALAR, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Superludist, NumericFactorization, Complex##SCALAR, LO, GO ) \

#  ifdef HAVE_TPETRA_INST_COMPLEX_FLOAT
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(LO, GO) \
  UNIT_TEST_GROUP_ORDINAL_COMPLEX_SCALAR(LO, GO, float)
#  else
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(LO, GO)
#  endif

#  ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_DOUBLE(LO, GO)        \
  UNIT_TEST_GROUP_ORDINAL_COMPLEX_SCALAR(LO, GO, double)
#  else
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_DOUBLE(LO, GO)
#  endif

#else  // !(defined HAVE_TEUCHOS_COMPLEX
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(LO, GO)
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_DOUBLE(LO, GO)
#endif

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

  // Uncomment this for really fast development cycles but make sure to comment
  // it back again before checking in so that we can test all the types.
  // #define FAST_DEVELOPMENT_UNIT_TEST_BUILD

#define UNIT_TEST_GROUP_ORDINAL_SCALAR( LO, GO, SCALAR )                \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Superludist, Initialization, SCALAR, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Superludist, SymbolicFactorization, SCALAR, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Superludist, NumericFactorization, SCALAR, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Superludist, NonContgGID, SCALAR, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Superludist, ContgGID, SCALAR, LO, GO ) \

#define UNIT_TEST_GROUP_ORDINAL( ORDINAL )              \
  UNIT_TEST_GROUP_ORDINAL_ORDINAL( ORDINAL, ORDINAL )

#ifdef FAST_DEVELOPMENT_UNIT_TEST_BUILD
#  define UNIT_TEST_GROUP_ORDINAL_ORDINAL( LO, GO )     \
  UNIT_TEST_GROUP_ORDINAL_SCALAR( LO, GO, double)       \
  UNIT_TEST_GROUP_ORDINAL(int)
  UNIT_TEST_GROUP_ORDINAL_COMPLEX(LO,GO,float)

#else // not FAST_DEVELOPMENT_UNIT_TEST_BUILD

#  define UNIT_TEST_GROUP_ORDINAL_ORDINAL( LO, GO )     \
  UNIT_TEST_GROUP_ORDINAL_FLOAT(LO, GO)                 \
  UNIT_TEST_GROUP_ORDINAL_DOUBLE(LO, GO)                \
  UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(LO,GO)          \
  UNIT_TEST_GROUP_ORDINAL_COMPLEX_DOUBLE(LO,GO)

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

#endif // FAST_DEVELOPMENT_UNIT_TEST_BUILD

} // end anonymous namespace
