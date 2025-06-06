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
  using Amesos2::TachoSolver;

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
  RCP<FancyOStream> getDefaultOStream()
  {
    return( VerboseObjectBase::getDefaultOStream() );
  }
  */


  /*
   * UNIT TESTS
   */

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( Tacho, Initialization, SCALAR, LO, GO )
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
    typedef TachoSolver<MAT,MV> SOLVER;

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
    RCP<Amesos2::Solver<MAT,MV> > solver = Amesos2::create<MAT,MV>("Tacho",eye,X,B);

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


  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( Tacho, SymbolicFactorization, SCALAR, LO, GO )
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
    RCP<Amesos2::Solver<MAT,MV> > solver = Amesos2::create<MAT,MV>("Tacho",eye,X,B);

    solver->symbolicFactorization();
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( Tacho, NumericFactorization, SCALAR, LO, GO )
  {
    typedef ScalarTraits<SCALAR> ST;
    typedef CrsMatrix<SCALAR,LO,GO,Node> MAT;
    typedef MultiVector<SCALAR,LO,GO,Node> MV;

    const SCALAR ONE = ST::one();
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
      eye->insertGlobalValues(base+i,tuple<GO>(base+i),tuple<SCALAR>(ONE));
    }
    eye->fillComplete();

    // Create X
    RCP<MV> X = rcp(new MV(map,11));
    X->randomize();

    // Create B
    RCP<MV> B = rcp(new MV(map,11));
    B->randomize();

    // Constructor from Factory
    RCP<Amesos2::Solver<MAT,MV> > solver = Amesos2::create<MAT,MV>("Tacho",eye,X,B);

    solver->symbolicFactorization().numericFactorization();

    {
      // Modify matrix
      // NOTE: Running with one MPI, previously Amesos2/Tacho holds on to device_nzvals_views of eye
      //                             which leads to "annot access data on host while a device view is alive"
      //                             when replacing values of eye (PR 14019)
      eye->resumeFill();
      for( size_t i = 0; i < numLocal; ++i ){
        eye->replaceLocalValues(base+i,tuple<LO>(i),tuple<SCALAR>(ONE+ONE));
      }
      eye->fillComplete();

      // Refactor
      solver->numericFactorization();
    }
    // Good way to check the factors L and U?  Probs not, since they are private members
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( Tacho, Solve, SCALAR, LO, GO )
  {
    typedef CrsMatrix<SCALAR,LO,GO,Node> MAT;
    typedef ScalarTraits<SCALAR> ST;
    typedef MultiVector<SCALAR,LO,GO,Node> MV;
    typedef typename ST::magnitudeType Mag;
    const size_t numVecs = 7;

    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();

    out << "Reading file" << std::endl;

    RCP<MAT> A =
      Tpetra::MatrixMarket::Reader<MAT>::readSparseFile("../matrices/amesos2_test_mat1.mtx",comm);

    out << "Done reading file" << std::endl;

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

    out << "Creating solver" << std::endl;

    // Solve A*Xhat = B for Xhat using the Tacho solver
    RCP<Amesos2::Solver<MAT,MV> > solver
      = Amesos2::create<MAT,MV>("Tacho", A, Xhat, B );

    // Use LU
    Teuchos::ParameterList amesos2_params("Amesos2");
    amesos2_params.sublist("Tacho").set("method", "lu", "Factorization Time");
    solver->setParameters( rcpFromRef(amesos2_params) );
    out << "Done creating solver" << std::endl;

    solver->symbolicFactorization();
    solver->numericFactorization();
    out << "Doing solve" << std::endl;
    solver->solve();
    out << "Solve finished" << std::endl;

    Xhat->describe(out, Teuchos::VERB_EXTREME);
    X->describe(out, Teuchos::VERB_EXTREME);

    // Check result of solve
    Array<Mag> xhatnorms(numVecs), xnorms(numVecs);
    Xhat->norm2(xhatnorms());
    X->norm2(xnorms());
    TEST_COMPARE_FLOATING_ARRAYS( xhatnorms, xnorms, 0.005 );
  }

  /*
   * Unit Tests for Complex data types
   */

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( Tacho, ComplexSolve, SCALAR, LO, GO )
  {
    typedef std::complex<SCALAR> cmplx;
    typedef CrsMatrix<cmplx,LO,GO,Node> MAT;
    typedef ScalarTraits<cmplx> ST;
    typedef MultiVector<cmplx,LO,GO,Node> MV;
    typedef typename ST::magnitudeType Mag;
    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();

    RCP<MAT> A =
      Tpetra::MatrixMarket::Reader<MAT>::readSparseFile("../matrices/amesos2_test_mat4.mtx",comm);

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

    // Create solver interface to Tacho through Amesos2 factory method
    RCP<Amesos2::Solver<MAT,MV> > solver = Amesos2::create<MAT,MV>("Tacho",A,Xhat,B);

    // Use LU
    Teuchos::ParameterList amesos2_params("Amesos2");
    amesos2_params.sublist("Tacho").set("method", "lu", "Factorization Time");
    solver->setParameters( rcpFromRef(amesos2_params) );
    out << "Done creating solver" << std::endl;

    solver->symbolicFactorization().numericFactorization().solve();

    Xhat->describe(out, Teuchos::VERB_EXTREME);
    X->describe(out, Teuchos::VERB_EXTREME);

    // Check result of solve
    Array<Mag> xhatnorms(1), xnorms(1);
    Xhat->norm2(xhatnorms());
    X->norm2(xnorms());
    TEST_COMPARE_FLOATING_ARRAYS( xhatnorms, xnorms, 0.005 );
  }

  TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL( Tacho, ComplexSolve2, SCALAR, LO, GO )
  {
    typedef std::complex<SCALAR> cmplx;
    typedef CrsMatrix<cmplx,LO,GO,Node> MAT;
    typedef ScalarTraits<cmplx> ST;
    typedef MultiVector<cmplx,LO,GO,Node> MV;
    typedef typename ST::magnitudeType Mag;
    const size_t numVecs = 7;

    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();

    RCP<MAT> A =
      Tpetra::MatrixMarket::Reader<MAT>::readSparseFile("../matrices/amesos2_test_mat2.mtx",comm);

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

    // Solve A*Xhat = B for Xhat using the Tacho solver
    RCP<Amesos2::Solver<MAT,MV> > solver
      = Amesos2::create<MAT,MV>("Tacho", A, Xhat, B);

    // Use LU
    Teuchos::ParameterList amesos2_params("Amesos2");
    amesos2_params.sublist("Tacho").set("method", "lu", "Factorization Time");
    solver->setParameters( rcpFromRef(amesos2_params) );
    out << "Done creating solver" << std::endl;

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
#ifdef HAVE_TEUCHOS_COMPLEX
#  define UNIT_TEST_GROUP_ORDINAL_COMPLEX_SCALAR(LO, GO, SCALAR)        \
  typedef std::complex<SCALAR>  Complex##SCALAR;                        \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Tacho, Initialization, Complex##SCALAR, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Tacho, SymbolicFactorization, Complex##SCALAR, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Tacho, NumericFactorization, Complex##SCALAR, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Tacho, ComplexSolve, SCALAR, LO, GO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Tacho, ComplexSolve2, SCALAR, LO, GO)

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
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Tacho, Initialization, SCALAR, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Tacho, SymbolicFactorization, SCALAR, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Tacho, NumericFactorization, SCALAR, LO, GO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT( Tacho, Solve, SCALAR, LO, GO )


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
