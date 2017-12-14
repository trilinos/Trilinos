// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package
//                  Copyright 2011 Sandia Corporation
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
// ***********************************************************************
//
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

#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <MatrixMarket_Tpetra.hpp> // for loading matrices from file

#include "Amesos2.hpp"
#include "Amesos2_Meta.hpp"

//#include "Amesos2_Basker_decl.hpp"
//#include "Amesos2_Basker_def.hpp"


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
  using Amesos2::Basker;
  using Amesos2::Meta::is_same;

  typedef Tpetra::DefaultPlatform::DefaultPlatformType Platform;
  typedef Platform::NodeType Node;

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
      ret = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
    } else {
      ret = rcp(new Teuchos::SerialComm<int>());
    }
    return ret;
  }

  RCP<FancyOStream> getDefaultOStream()
  {
    return( VerboseObjectBase::getDefaultOStream() );
  }


  /*
   * UNIT TESTS
   */

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( Basker, Initialization, SCALAR, LO, GO )
  {
    /* Test correct initialization of the Solver instance
     *
     * - All Constructors
     * - Correct initialization of class members
     * - Correct typedefs ( using Amesos2::is_same<> )
     */
    typedef ScalarTraits<SCALAR> ST;
    typedef CrsMatrix<SCALAR,LO,GO,Node> MAT;
    typedef MultiVector<SCALAR,LO,GO,Node> MV;
    //typedef Basker<MAT,MV> SOLVER;

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
    RCP<Amesos2::Solver<MAT,MV> > solver = Amesos2::create<MAT,MV>("Basker",eye,X,B);

    TEST_ASSERT( solver->getStatus().getNumSymbolicFact() == 0 );
    TEST_ASSERT( solver->getStatus().getNumNumericFact() == 0 );
    TEST_ASSERT( solver->getStatus().getNumSolve() == 0 );

    // The following should all pass at compile time
    //TEST_ASSERT( (is_same<MAT,typename SOLVER::matrix_type>::value) );
    //TEST_ASSERT( (is_same<MV,typename SOLVER::vector_type>::value) );
    //TEST_ASSERT( (is_same<SCALAR,typename SOLVER::scalar_type>::value) );
    //TEST_ASSERT( (is_same<LO,typename SOLVER::local_ordinal_type>::value) );
    //TEST_ASSERT( (is_same<GO,typename SOLVER::global_ordinal_type>::value) );
    //TEST_ASSERT( (is_same<global_size_t,typename SOLVER::global_size_type>::value) );
    // TEST_ASSERT( (is_same<Node,typename SOLVER::node_type>::value) );
  }


  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( Basker, SymbolicFactorization, SCALAR, LO, GO )
  {
    typedef ScalarTraits<SCALAR> ST;
    typedef CrsMatrix<SCALAR,LO,GO,Node> MAT;
    typedef MultiVector<SCALAR,LO,GO,Node> MV;
    //typedef Basker<MAT,MV> SOLVER;

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
    RCP<Amesos2::Solver<MAT,MV> > solver = Amesos2::create<MAT,MV>("Basker",eye,X,B);

    solver->symbolicFactorization();
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( Basker, NumericFactorization, SCALAR, LO, GO )
  {
    typedef ScalarTraits<SCALAR> ST;
    typedef CrsMatrix<SCALAR,LO,GO,Node> MAT;
    typedef MultiVector<SCALAR,LO,GO,Node> MV;
    //typedef Basker<MAT,MV> SOLVER;

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
    RCP<Amesos2::Solver<MAT,MV> > solver = Amesos2::create<MAT,MV>("Basker",eye,X,B);

    solver->symbolicFactorization().numericFactorization();

    if ( rank == 0 ) {
      TEST_ASSERT( solver->getStatus().getNnzLU() != 0 );
      TEST_ASSERT( solver->getStatus().getNnzLU() == numLocal*static_cast<const size_t>(2*comm->getSize()) );
      // Good way to check the factors L and U?  Probs not, since they are private members
    }
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( Basker, Solve, SCALAR, LO, GO )
  {
    typedef CrsMatrix<SCALAR,LO,GO,Node> MAT;
    typedef ScalarTraits<SCALAR> ST;
    typedef MultiVector<SCALAR,LO,GO,Node> MV;
    typedef typename ST::magnitudeType Mag;
    //typedef ScalarTraits<Mag> MT;
    const size_t numVecs = 1;

    Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
    RCP<const Comm<int> > comm = platform.getComm();
    RCP<Node>             node = platform.getNode();

    RCP<MAT> A =
      Tpetra::MatrixMarket::Reader<MAT>::readSparseFile("../matrices/amesos2_test_mat1.mtx",comm,node);


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
    Xhat->describe(*(getDefaultOStream()), Teuchos::VERB_EXTREME);
    X->describe(*(getDefaultOStream()), Teuchos::VERB_EXTREME);
    B->describe(*(getDefaultOStream()), Teuchos::VERB_EXTREME);


    // Solve A*Xhat = B for Xhat using the Bakser solver
    RCP<Amesos2::Solver<MAT,MV> > solver
      = Amesos2::create<MAT,MV>("Basker", A, Xhat, B );

    solver->symbolicFactorization();
    solver->numericFactorization();
    solver->solve();

    Xhat->describe(*(getDefaultOStream()), Teuchos::VERB_EXTREME);
    X->describe(*(getDefaultOStream()), Teuchos::VERB_EXTREME);
    B->describe(*(getDefaultOStream()), Teuchos::VERB_EXTREME);

    // Check result of solve
    Array<Mag> xhatnorms(numVecs), xnorms(numVecs);
    Xhat->norm2(xhatnorms());
    X->norm2(xnorms());
    TEST_COMPARE_FLOATING_ARRAYS( xhatnorms, xnorms, 0.005 );
  }

 /* TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( KLU2, SolveTrans, SCALAR, LO, GO )
  {
    typedef CrsMatrix<SCALAR,LO,GO,Node> MAT;
    typedef ScalarTraits<SCALAR> ST;
    typedef MultiVector<SCALAR,LO,GO,Node> MV;
    typedef typename ST::magnitudeType Mag;
    typedef ScalarTraits<Mag> MT;
    const size_t numVecs = 7;

    Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
    RCP<const Comm<int> > comm = platform.getComm();
    RCP<Node>             node = platform.getNode();

    RCP<MAT> A =
      Tpetra::MatrixMarket::Reader<MAT>::readSparseFile("../matrices/amesos2_test_mat1.mtx",comm,node);

    RCP<const Map<LO,GO,Node> > dmnmap = A->getDomainMap();
    RCP<const Map<LO,GO,Node> > rngmap = A->getRangeMap();

    RCP<MV> X = rcp(new MV(rngmap,numVecs));
    RCP<MV> B = rcp(new MV(dmnmap,numVecs));
    RCP<MV> Xhat = rcp(new MV(rngmap,numVecs));
    X->setObjectLabel("X");
    B->setObjectLabel("B");
    Xhat->setObjectLabel("Xhat");
    X->randomize();

    A->apply(*X,*B,Teuchos::TRANS);

    Xhat->randomize();

    // Solve A*Xhat = B for Xhat using the KLU2 solver
    cout <<"I am in solvetrans create" << endl;
    RCP<Amesos2::Solver<MAT,MV> > solver
      = Amesos2::create<MAT,MV>("KLU2", A, Xhat, B );

    Teuchos::ParameterList amesos2_params("Amesos2");
    amesos2_params.sublist("KLU2").set("Trans","TRANS","Solve with transpose");

    cout <<"Setting parameters" << amesos2_params << endl;
    solver->setParameters( rcpFromRef(amesos2_params) );
    cout <<"Calling everything up to solve" << endl;
    solver->symbolicFactorization().numericFactorization().solve();

    Xhat->describe(out, Teuchos::VERB_EXTREME);
    X->describe(out, Teuchos::VERB_EXTREME);

    // Check result of solve
    Array<Mag> xhatnorms(numVecs), xnorms(numVecs);
    Xhat->norm2(xhatnorms());
    X->norm2(xnorms());
    TEST_COMPARE_FLOATING_ARRAYS( xhatnorms, xnorms, 0.005 );
  }*/

  /*
   * Unit Tests for Complex data types
   */


  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( Basker, NonContgGID, SCALAR, LO, GO )
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

    Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
    RCP<const Comm<int> > comm = platform.getComm();
    RCP<Node>             node = platform.getNode();

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
          "Basker NonContigGID Test: The non-contiguous map claims to be contiguous.");

      //RCP<MAT> A = rcp( new MAT(map,3) ); // max of three entries in a row
      RCP<MAT> A = rcp( new MAT(map,0) );
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
          "Basker NonContigGID Test: The non-contiguous map of A claims to be contiguous.");


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
          "Basker NonContigGID Test: The non-contiguous maps of X or B claims to be contiguous.");


      // Create solver interface with Amesos2 factory method
      RCP<Amesos2::Solver<MAT,MV> > solver = Amesos2::create<MAT,MV>("Basker", A, X, B);

      // Create a Teuchos::ParameterList to hold solver parameters
      Teuchos::ParameterList amesos2_params("Amesos2");
      amesos2_params.sublist("Basker").set("IsContiguous", false, "Are GIDs Contiguous");

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

      A->describe(out, Teuchos::VERB_EXTREME);
      B->describe(out, Teuchos::VERB_EXTREME);
      Xhat->describe(out, Teuchos::VERB_EXTREME);
      X->describe(out, Teuchos::VERB_EXTREME);

      //A->describe(*(getDefaultOStream()), Teuchos::VERB_EXTREME);
      //B->describe(*(getDefaultOStream()), Teuchos::VERB_EXTREME);
      //Xhat->describe(*(getDefaultOStream()), Teuchos::VERB_EXTREME);
      //X->describe(*(getDefaultOStream()), Teuchos::VERB_EXTREME);

      // Check result of solve
      Array<Mag> xhatnorms(numVectors), xnorms(numVectors);
      Xhat->norm2(xhatnorms());
      X->norm2(xnorms());
      TEST_COMPARE_FLOATING_ARRAYS( xhatnorms, xnorms, 0.005 );
    } // end if numProcs = 2
  }


  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( Basker, ComplexSolve, SCALAR, LO, GO )
  {
    typedef std::complex<SCALAR> cmplx;
    typedef CrsMatrix<cmplx,LO,GO,Node> MAT;
    typedef ScalarTraits<cmplx> ST;
    typedef MultiVector<cmplx,LO,GO,Node> MV;
    typedef typename ST::magnitudeType Mag;
    //typedef ScalarTraits<Mag> MT;
    Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
    RCP<const Comm<int> > comm = platform.getComm();
    RCP<Node>             node = platform.getNode();

    RCP<MAT> A =
      Tpetra::MatrixMarket::Reader<MAT>::readSparseFile("../matrices/amesos2_test_mat4.mtx",comm,node);

    RCP<const Map<LO,GO,Node> > dmnmap = A->getDomainMap();
    RCP<const Map<LO,GO,Node> > rngmap = A->getRangeMap();

    // Create the know-solution vector
    std::map<GO,cmplx> xValues;
    xValues[0] = cmplx(0.552267,1.22818);
    xValues[1] = cmplx(-0.290371,-0.600974);
    xValues[2] = cmplx(-0.629824,-0.340952);

    typename std::map<GO,cmplx>::iterator it;
    RCP<MV> X = rcp(new MV(dmnmap, 1));
    X->setObjectLabel("X");

    for( it = xValues.begin(); it != xValues.end(); ++it ){
      if( rngmap->isNodeGlobalElement( (*it).first ) ){
        out << "replacing global row " << (*it).first << " with " << (*it).second << std::endl;
        X->replaceGlobalValue( (*it).first, 0, (*it).second );
      }
    }

    // Create the rhs vector B
    std::map<GO,cmplx> bValues;
    bValues[0] = cmplx(1,1);
    bValues[1] = cmplx(3,-2);
    bValues[2] = cmplx(1,-2);

    RCP<MV> B = rcp(new MV(rngmap, 1));
    B->setObjectLabel("B");

    for( it = bValues.begin(); it != bValues.end(); ++it ){
      if( rngmap->isNodeGlobalElement( (*it).first ) ){
        out << "replacing global row " << (*it).first << " with " << (*it).second << std::endl;
        B->replaceGlobalValue( (*it).first, 0, (*it).second );
      }
    }

    // Create the solution vector
    RCP<MV> Xhat = rcp(new MV(dmnmap,1));
    Xhat->setObjectLabel("Xhat");

    // Create solver interface to Basker through Amesos2 factory method
    RCP<Amesos2::Solver<MAT,MV> > solver = Amesos2::create<MAT,MV>("Basker",A,Xhat,B);

    solver->symbolicFactorization().numericFactorization().solve();

    Xhat->describe(out, Teuchos::VERB_EXTREME);
    X->describe(out, Teuchos::VERB_EXTREME);

    // Check result of solve
    Array<Mag> xhatnorms(1), xnorms(1);
    Xhat->norm2(xhatnorms());
    X->norm2(xnorms());
    TEST_COMPARE_FLOATING_ARRAYS( xhatnorms, xnorms, 0.005 );
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( Basker, ComplexSolve2, SCALAR, LO, GO )
  {
    typedef std::complex<SCALAR> cmplx;
    typedef CrsMatrix<cmplx,LO,GO,Node> MAT;
    typedef ScalarTraits<cmplx> ST;
    typedef MultiVector<cmplx,LO,GO,Node> MV;
    typedef typename ST::magnitudeType Mag;
    //typedef ScalarTraits<Mag> MT;
    const size_t numVecs = 7;

    Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
    RCP<const Comm<int> > comm = platform.getComm();
    RCP<Node>             node = platform.getNode();

    RCP<MAT> A =
      Tpetra::MatrixMarket::Reader<MAT>::readSparseFile("../matrices/amesos2_test_mat2.mtx",comm,node);

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

    // Solve A*Xhat = B for Xhat using the Basker solver
    RCP<Amesos2::Solver<MAT,MV> > solver
      = Amesos2::create<MAT,MV>("Basker", A, Xhat, B);

    solver->symbolicFactorization().numericFactorization().solve();

    Xhat->describe(out, Teuchos::VERB_EXTREME);
    X->describe(out, Teuchos::VERB_EXTREME);

    // Check result of solve
    Array<Mag> xhatnorms(numVecs), xnorms(numVecs);
    Xhat->norm2(xhatnorms());
    X->norm2(xnorms());
    TEST_COMPARE_FLOATING_ARRAYS( xhatnorms, xnorms, 0.005 );
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( Basker, ComplexSolve2Trans, SCALAR, LO, GO )
  {
    typedef std::complex<SCALAR> cmplx;
    typedef CrsMatrix<cmplx,LO,GO,Node> MAT;
    typedef ScalarTraits<cmplx> ST;
    typedef MultiVector<cmplx,LO,GO,Node> MV;
    typedef typename ST::magnitudeType Mag;
    //typedef ScalarTraits<Mag> MT;
    const size_t numVecs = 7;

    Platform &platform = Tpetra::DefaultPlatform::getDefaultPlatform();
    RCP<const Comm<int> > comm = platform.getComm();
    RCP<Node>             node = platform.getNode();

    RCP<MAT> A =
      Tpetra::MatrixMarket::Reader<MAT>::readSparseFile("../matrices/amesos2_test_mat3.mtx",comm,node);

    RCP<const Map<LO,GO,Node> > dmnmap = A->getDomainMap();
    RCP<const Map<LO,GO,Node> > rngmap = A->getRangeMap();

    RCP<MV> X = rcp(new MV(dmnmap,numVecs));
    RCP<MV> B = rcp(new MV(rngmap,numVecs));
    RCP<MV> Xhat = rcp(new MV(dmnmap,numVecs));
    X->setObjectLabel("X");
    B->setObjectLabel("B");
    Xhat->setObjectLabel("Xhat");
    X->randomize();

    A->apply(*X,*B,Teuchos::CONJ_TRANS); // use conjugate-transpose

    Xhat->randomize();

    // Solve A*Xhat = B for Xhat using the Basker solver
    RCP<Amesos2::Solver<MAT,MV> > solver
      = Amesos2::create<MAT,MV>("Basker", A, Xhat, B);

    Teuchos::ParameterList amesos2_params("Amesos2");
    amesos2_params.sublist("Basker").set("Trans","CONJ","Solve with conjugate-transpose");

    solver->setParameters( rcpFromRef(amesos2_params) );
    solver->symbolicFactorization().numericFactorization().solve();

    Xhat->describe(out, Teuchos::VERB_EXTREME);
    X->describe(out, Teuchos::VERB_EXTREME);

    // Check result of solve
    Array<Mag> xhatnorms(numVecs), xnorms(numVecs);
    Xhat->norm2(xhatnorms());
    X->norm2(xnorms());
    TEST_COMPARE_FLOATING_ARRAYS( xhatnorms, xnorms, 0.005 );
  }


  /*
   * Instantiations
   */

#if defined(HAVE_TEUCHOS_COMPLEX)

  // mfh 11 Jan 2016: Clang 3.7 doesn't like the following
  // commented-out pragma.  It says: "error: pragma message requires
  // parenthesized string".  I' m not sure how to fix that (the string
  // looks parenthesized to me!), so I'll comment this out for now,
  // just so that I can get Basker to build with Clang 3.7.
  //
  //#pragma message("T COMPLEX");

#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_SCALAR(LO, GO, SCALAR)        \
  typedef std::complex<SCALAR>  Complex##SCALAR;                        \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Basker, Initialization, Complex##SCALAR, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Basker, SymbolicFactorization, Complex##SCALAR, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Basker, NumericFactorization, Complex##SCALAR, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Basker, ComplexSolve, SCALAR, LO, GO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Basker, ComplexSolve2, SCALAR, LO, GO)

#  ifdef HAVE_TPETRA_INST_COMPLEX_FLOAT
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(LO, GO) \
  UNIT_TEST_GROUP_ORDINAL_COMPLEX_SCALAR(LO, GO, float)
#  else
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(LO, GO)
#  endif//end have complex_flox

#  ifdef HAVE_TPETRA_INST_COMPLEX_DOUBLE
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_DOUBLE(LO, GO)        \
  UNIT_TEST_GROUP_ORDINAL_COMPLEX_SCALAR(LO, GO, double)
#  else
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_DOUBLE(LO, GO)
#  endif//end complex_double
#else  // !(defined HAVE_TEUCHOS_COMPLEX
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_FLOAT(LO, GO)
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_DOUBLE(LO, GO)
#endif
  //#endif

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
  //TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( KLU2, SolveTrans, SCALAR, LO, GO )

#define UNIT_TEST_GROUP_ORDINAL_SCALAR( LO, GO, SCALAR )                \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Basker, Initialization, SCALAR, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Basker, SymbolicFactorization, SCALAR, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Basker, NumericFactorization, SCALAR, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Basker, Solve, SCALAR, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Basker, NonContgGID, SCALAR, LO, GO )

#define UNIT_TEST_GROUP_ORDINAL( ORDINAL )              \
  UNIT_TEST_GROUP_ORDINAL_ORDINAL( ORDINAL, ORDINAL )

#ifdef FAST_DEVELOPMENT_UNIT_TEST_BUILD
#  define UNIT_TEST_GROUP_ORDINAL_ORDINAL( LO, GO )     \
  UNIT_TEST_GROUP_ORDINAL_SCALAR( LO, GO, double)       \
  UNIT_TEST_GROUP_ORDINAL(int)

#else // not FAST_DEVELOPMENT_UNIT_TEST_BUILD

#  define UNIT_TEST_GROUP_ORDINAL_ORDINAL( LO, GO )     \
  UNIT_TEST_GROUP_ORDINAL_FLOAT(LO, GO)                 \
  UNIT_TEST_GROUP_ORDINAL_DOUBLE(LO, GO)                \
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
#endif  // EXPL-INST


#endif // FAST_DEVELOPMENT_UNIT_TEST_BUILD

} // end anonymous namespace
