// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_StackedTimer.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_ParameterXMLFileReader.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>

// I/O for Matrix-Market files
#include <MatrixMarket_Tpetra.hpp>
#include <Tpetra_Import.hpp>

#include <Amesos2.hpp>
#include <Amesos2_Version.hpp>

#if defined(HAVE_AMESOS2_XPETRA) && defined(HAVE_AMESOS2_ZOLTAN2)
# include <Zoltan2_OrderingProblem.hpp>
# include <Zoltan2_PartitioningProblem.hpp>
# include <Zoltan2_TpetraRowGraphAdapter.hpp>
# include <Zoltan2_TpetraCrsMatrixAdapter.hpp>
# include <Zoltan2_XpetraMultiVectorAdapter.hpp>
#endif


int main(int argc, char *argv[]) {
  Tpetra::ScopeGuard tpetraScope(&argc,&argv);

  typedef Tpetra::CrsMatrix<>::scalar_type Scalar;
  typedef Tpetra::Map<>::local_ordinal_type LO;
  typedef Tpetra::Map<>::global_ordinal_type GO;
  typedef Tpetra::Map<>::node_type NO;

  typedef Tpetra::RowGraph<LO, GO, NO> Graph;
  typedef Tpetra::CrsMatrix<Scalar,LO,GO> MAT;
  typedef Tpetra::MultiVector<Scalar,LO,GO> MV;

  using Tpetra::global_size_t;
  using Tpetra::Map;
  using Tpetra::Import;
  using Teuchos::RCP;
  using Teuchos::rcp;


  //
  // Get the default communicator
  //
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int myRank = comm->getRank();

  Teuchos::oblackholestream blackhole;

  bool printMatrix     = false;
  bool printSolution   = false;
  bool checkSolution   = false;
  bool printTiming     = false;
  bool useStackedTimer = false;
  bool allprint        = false;
  bool verbose = (myRank==0);
  bool useZoltan2 = false;
  bool useParMETIS = false;
  std::string mat_filename("arc130.mtx");
  std::string rhs_filename("");
  std::string solvername("Superlu");
  std::string xml_filename("");
  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
  cmdp.setOption("filename",&mat_filename,"Filename for Matrix-Market test matrix.");
  cmdp.setOption("rhs_filename",&rhs_filename,"Filename for Matrix-Market right-hand-side.");
  cmdp.setOption("solvername",&solvername,"Name of solver.");
  cmdp.setOption("xml_filename",&xml_filename,"XML Filename for Solver parameters.");
  cmdp.setOption("print-matrix","no-print-matrix",&printMatrix,"Print the full matrix after reading it.");
  cmdp.setOption("print-solution","no-print-solution",&printSolution,"Print solution vector after solve.");
  cmdp.setOption("check-solution","no-check-solution",&checkSolution,"Check solution vector after solve.");
  cmdp.setOption("use-zoltan2","no-zoltan2",&useZoltan2,"Use Zoltan2 (Hypergraph) for repartitioning");
  cmdp.setOption("use-parmetis","no-parmetis",&useParMETIS,"Use ParMETIS for repartitioning");
  cmdp.setOption("print-timing","no-print-timing",&printTiming,"Print solver timing statistics");
  cmdp.setOption("use-stacked-timer","no-stacked-timer",&useStackedTimer,"Use StackedTimer to print solver timing statistics");
  cmdp.setOption("all-print","root-print",&allprint,"All processors print to out");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }

  std::ostream& out = ( (allprint || (myRank == 0)) ? std::cout : blackhole );
  RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));

  // Say hello
  out << myRank << " : " << Amesos2::version() << std::endl << std::endl;

  const size_t numVectors = 1;

  // Read matrix
  RCP<const MAT> A = Tpetra::MatrixMarket::Reader<MAT>::readSparseFile(mat_filename, comm);

  // get the map (Range Map used for both X & B)
  RCP<const Map<LO,GO> > rngmap = A->getRangeMap();
  RCP<const Map<LO,GO> > dmnmap = A->getDomainMap();
  GO nrows = A->getGlobalNumRows();

  // Create random X
  RCP<MV> X = rcp(new MV(dmnmap,numVectors));
  X->randomize();

  // Create B
  RCP<MV> B = rcp(new MV(rngmap,numVectors));
  if (rhs_filename == "") {
    /*
     * Use RHS:
     *
     *  [[10]
     *   [10]
     *   [10]
     *   [10]
     *   [10]
     *   [10]]
     */
    B->putScalar(10);
  } else {
    B = Tpetra::MatrixMarket::Reader<MAT>::readDenseFile (rhs_filename, comm, rngmap);
  }

  if (useZoltan2 || useParMETIS) {
#if defined(HAVE_AMESOS2_XPETRA) && defined(HAVE_AMESOS2_ZOLTAN2)
    // Specify partitioning parameters
    Teuchos::ParameterList zoltan_params;
    zoltan_params.set("partitioning_approach", "partition");
    //
    if (useParMETIS) {
      if (comm->getRank() == 0) {
        std::cout << "Using Zoltan2(ParMETIS)" << std::endl;
      }
      zoltan_params.set("algorithm", "parmetis");
      zoltan_params.set("symmetrize_input", "transpose");
      zoltan_params.set("partitioning_objective", "minimize_cut_edge_weight");
    } else {
      if (comm->getRank() == 0) {
        std::cout << "Using Zoltan2(HyperGraph)" << std::endl;
      }
      zoltan_params.set("algorithm", "phg");
    }

    // Create an input adapter for the Tpetra matrix.
    Zoltan2::TpetraRowGraphAdapter<Graph>
      zoltan_graph(A->getGraph());

    // Create and solve partitioning problem
    Zoltan2::PartitioningProblem<Zoltan2::TpetraRowGraphAdapter<Graph>>
      problem(&zoltan_graph, &zoltan_params);
    problem.solve();

    // Redistribute matrix
    RCP<MAT> zoltan_A;
    Zoltan2::TpetraCrsMatrixAdapter<MAT> zoltan_matrix(A);
    zoltan_matrix.applyPartitioningSolution (*A, zoltan_A, problem.getSolution());
    // Set it as coefficient matrix, and update range map
    A = zoltan_A;
    rngmap = A->getRangeMap();

    // Redistribute RHS
    RCP<MV> zoltan_b;
    Zoltan2::XpetraMultiVectorAdapter<MV> adapterRHS(rcpFromRef (*B));
    adapterRHS.applyPartitioningSolution (*B, zoltan_b, problem.getSolution());
    // Set it as RHS
    B = zoltan_b;

    // Redistribute Sol
    RCP<MV> zoltan_x;
    Zoltan2::XpetraMultiVectorAdapter<MV> adapterSol(rcpFromRef (*X));
    adapterSol.applyPartitioningSolution (*X, zoltan_x, problem.getSolution());
    // Set it as Sol
    X = zoltan_x;
#else
    TEUCHOS_TEST_FOR_EXCEPTION(
      useZoltan2, std::invalid_argument,
      "Both Xpetra and Zoltan2 are needed to use Zoltan2.");
#endif
  }
  if( printMatrix ){
    A->describe(*fos, Teuchos::VERB_EXTREME);
  }
  else if( verbose ){
    std::cout << std::endl << A->description() << std::endl << std::endl;
  }

  // Constructor from Factory
  RCP<Amesos2::Solver<MAT,MV> > solver;
  if( !Amesos2::query(solvername) ){
    *fos << solvername << " solver not enabled.  Exiting..." << std::endl;
    return EXIT_SUCCESS;
  }

  solver = Amesos2::create<MAT,MV>(solvername, A, X, B);
  if (xml_filename != "") {
    Teuchos::ParameterList test_params =
      Teuchos::ParameterXMLFileReader(xml_filename).getParameters();
    Teuchos::ParameterList& amesos2_params = test_params.sublist("Amesos2");
    *fos << amesos2_params.currentParametersString() << std::endl;
    solver->setParameters( Teuchos::rcpFromRef(amesos2_params) );
  }

  RCP<Teuchos::StackedTimer> stackedTimer;
  if(useStackedTimer) {
    stackedTimer = rcp(new Teuchos::StackedTimer("Amesos2 SimpleSolve-File"));
    Teuchos::TimeMonitor::setStackedTimer(stackedTimer);
  }
  solver->symbolicFactorization().numericFactorization().solve();
  if(useStackedTimer) {
    stackedTimer->stopBaseTimer();
  }

  if( printSolution ){
    // Print the solution
    RCP<Map<LO,GO> > root_map
      = rcp( new Map<LO,GO>(nrows,myRank == 0 ? nrows : 0,0,comm) );
    RCP<MV> Xhat = rcp( new MV(root_map,numVectors) );
    RCP<Import<LO,GO> > importer = rcp( new Import<LO,GO>(rngmap,root_map) );
    if( allprint ){
      if( myRank == 0 ) *fos << "Solution :" << std::endl;
      Xhat->describe(*fos,Teuchos::VERB_EXTREME);
      *fos << std::endl;
    } else {
      Xhat->doImport(*X,*importer,Tpetra::REPLACE);
      if( myRank == 0 ){
        *fos << "Solution :" << std::endl;
        Xhat->describe(*fos,Teuchos::VERB_EXTREME);
        *fos << std::endl;
      }
    }
  }

  if( checkSolution ){
    const Scalar one = Teuchos::ScalarTraits<Scalar>::one ();
    RCP<MV> R = rcp(new MV(rngmap,numVectors));
    A->apply(*X, *R);
    R->update(one, *B, -one);
    for (size_t j = 0; j < numVectors; ++j) {
      auto Rj = R->getVector(j);
      auto Bj = B->getVector(j);
      auto r_norm = Rj->norm2();
      auto b_norm = Bj->norm2();
      if (myRank == 0) {
        *fos << "Relative Residual norm = " << r_norm << " / " << b_norm << " = "
             << r_norm / b_norm << std::endl;
      }
    }
    if (myRank == 0) *fos << std::endl;
  }

  if(useStackedTimer) {
    Teuchos::StackedTimer::OutputOptions options;
    options.num_histogram=3;
    options.print_warnings = false;
    options.output_histogram = true;
    options.output_fraction=true;
    options.output_minmax = true;
    stackedTimer->report(std::cout, comm, options);
  } else if( printTiming ){
    // Print some timing statistics
    solver->printTiming(*fos);
  } else {
    Teuchos::TimeMonitor::summarize();
  }

  // We are done.
  return 0;
}
