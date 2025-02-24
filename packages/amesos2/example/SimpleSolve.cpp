// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/**
   \file   SimpleSolve.cpp
   \author Eric Bavier <etbavie@sandia.gov>
   \date   Sat Jul 17 10:35:39 2010

   \brief  Simple example of Amesos2 usage.

   This example solves a simple sparse system of linear equations using the
   Amesos2 interface to the Superlu solver.
*/

#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_StackedTimer.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_ParameterXMLFileReader.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include "Amesos2.hpp"
#include "Amesos2_Version.hpp"
#if defined(HAVE_AMESOS2_XPETRA) && defined(HAVE_AMESOS2_GALERI)
 #include "Galeri_XpetraMaps.hpp"
 #include "Galeri_XpetraProblemFactory.hpp"
#endif

int main(int argc, char *argv[]) {
  Tpetra::ScopeGuard tpetraScope(&argc,&argv);

  typedef Tpetra::CrsMatrix<>::scalar_type Scalar;
  typedef Tpetra::Map<>::local_ordinal_type LO;
  typedef Tpetra::Map<>::global_ordinal_type GO;

  typedef Tpetra::Map<LO,GO> MAP;
  typedef Tpetra::CrsMatrix<Scalar,LO,GO> MAT;
  typedef Tpetra::MultiVector<Scalar,LO,GO> MV;

  using Tpetra::global_size_t;
  using Teuchos::tuple;
  using Teuchos::RCP;
  using Teuchos::rcp;

  bool verbose = false;
  GO nx = 1;
  std::string GaleriName3D {"Laplace3D"};
  std::string solvername("Superlu");
  std::string xml_filename("");
  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("solvername",&solvername,"Name of solver.");
  cmdp.setOption("xml_filename",&xml_filename,"XML Filename for Solver parameters.");
  cmdp.setOption("nx",&nx,"Dimension of 3D problem.");
  cmdp.setOption ("galeriMatrixName", &GaleriName3D, "Name of 3D Galeri Matrix");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }

  // Before we do anything, check that SuperLU is enabled
  if( !Amesos2::query(solvername) ){
    std::cerr << solvername << " not enabled.  Exiting..." << std::endl;
    return EXIT_SUCCESS;        // Otherwise CTest will pick it up as
                                // failure, which it isn't really
  }

  Teuchos::RCP<const Teuchos::Comm<int> > comm =
    Tpetra::getDefaultComm();

  size_t myRank = comm->getRank();

  std::ostream &out = std::cout;

  out << Amesos2::version() << std::endl << std::endl;

  const size_t numVectors = 1;

  RCP<MAT> A;
  RCP<const MAP> map;
  if (nx > 0) {
    #if defined(HAVE_AMESOS2_XPETRA) && defined(HAVE_AMESOS2_GALERI)
    typedef Galeri::Xpetra::Problem<MAP, MAT, MV> Galeri_t;

    Teuchos::ParameterList galeriList;
    Tpetra::global_size_t nGlobalElements = nx * nx * nx;
    galeriList.set("nx", nx);
    galeriList.set("ny", nx);
    galeriList.set("nz", nx);
    if (GaleriName3D == "Elasticity3D") {
      GO mx = 1;
      galeriList.set("mx", mx);
      galeriList.set("my", mx);
      galeriList.set("mz", mx);
      nGlobalElements *= 3;
    }
    map = rcp(new MAP(nGlobalElements, 0, comm));
    RCP<Galeri_t> galeriProblem =
                  Galeri::Xpetra::BuildProblem<Scalar, LO, GO, MAP, MAT, MV>
                                     (GaleriName3D, map, galeriList);
    A = galeriProblem->BuildMatrix();
    #else
    std::cerr << "Galeri or Xpetra not enabled.  Exiting..." << std::endl;
    return EXIT_SUCCESS;        // Otherwise CTest will pick it up as
    #endif
  } else {
    // create a Map
    global_size_t nrows = 6;
    map = rcp( new MAP(nrows,0,comm) );

    A = rcp( new MAT(map,3) ); // max of three entries in a row

    /*
     * We will solve a system with a known solution, for which we will be using
     * the following matrix:
     *
     * [ [ 7,  0,  -3, 0,  -1, 0 ]
     *   [ 2,  8,  0,  0,  0,  0 ]
     *   [ 0,  0,  1,  0,  0,  0 ]
     *   [ -3, 0,  0,  5,  0,  0 ]
     *   [ 0,  -1, 0,  0,  4,  0 ]
     *   [ 0,  0,  0,  -2, 0,  6 ] ]
     *
     */
    // Construct matrix
    if( myRank == 0 ){
      A->insertGlobalValues(0,tuple<GO>(0,2,4),tuple<Scalar>(7,-3,-1));
      A->insertGlobalValues(1,tuple<GO>(0,1),tuple<Scalar>(2,8));
      A->insertGlobalValues(2,tuple<GO>(2),tuple<Scalar>(1));
      A->insertGlobalValues(3,tuple<GO>(0,3),tuple<Scalar>(-3,5));
      A->insertGlobalValues(4,tuple<GO>(1,4),tuple<Scalar>(-1,4));
      A->insertGlobalValues(5,tuple<GO>(3,5),tuple<Scalar>(-2,6));
    }
    A->fillComplete();
  }

  // Create X
  RCP<MV> X = rcp(new MV(map,numVectors));
  X->putScalar(1);

  /* Create B */
  RCP<MV> B = rcp(new MV(map,numVectors));
  A->apply(*X, *B);
  X->randomize();

  // Create solver interface with Amesos2 factory method
  RCP<Amesos2::Solver<MAT,MV> > solver = Amesos2::create<MAT,MV>(solvername, A, X, B);
  if (xml_filename != "") {
    Teuchos::ParameterList test_params =
      Teuchos::ParameterXMLFileReader(xml_filename).getParameters();
    Teuchos::ParameterList& amesos2_params = test_params.sublist("Amesos2");
    solver->setParameters( Teuchos::rcpFromRef(amesos2_params) );
  }

  RCP<Teuchos::StackedTimer> stackedTimer;
  stackedTimer = rcp(new Teuchos::StackedTimer("Amesos2 SimpleSolve-File"));
  Teuchos::TimeMonitor::setStackedTimer(stackedTimer);
  {
    solver->symbolicFactorization().numericFactorization().solve();
  }
  stackedTimer->stopBaseTimer();
  {
    Teuchos::StackedTimer::OutputOptions options;
    options.num_histogram=3;
    options.print_warnings = false;
    options.output_histogram = true;
    options.output_fraction=true;
    options.output_minmax = true;
    stackedTimer->report(std::cout, comm, options);
  }
  if (verbose) {
    /* Print the solution
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
    RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));
    *fos << "Solution :" << std::endl;
    X->describe(*fos,Teuchos::VERB_EXTREME);
    *fos << std::endl;
  }
  // We are done.
  return 0;
}
