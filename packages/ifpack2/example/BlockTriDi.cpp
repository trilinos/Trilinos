#include <Ifpack2_Factory.hpp>
#include <Ifpack2_BlockTriDiContainer.hpp>
#include <BelosTpetraAdapter.hpp>
#include <BelosSolverFactory.hpp>
#include <MatrixMarket_Tpetra.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_BlockCrsMatrix.hpp>
#include <Tpetra_BlockCrsMatrix_Helpers.hpp>

#include <Tpetra_Core.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterXMLFileReader.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include "Teuchos_StackedTimer.hpp"
#include <Teuchos_StandardCatchMacros.hpp>

namespace { // (anonymous)

// Values of command-line arguments.
struct CmdLineArgs {
  CmdLineArgs ():blockSize(-1),numIters(10),tol(1e-12),nx(172),lpp(10),useStackedTimer(false){}

  std::string mapFilename;
  std::string matrixFilename;
  std::string rhsFilename;
  std::string lineFilename;
  int blockSize;
  int numIters;
  double tol;
  int nx;
  int lpp;
  bool useStackedTimer;
  std::string problemName;
};

// Read in values of command-line arguments.
bool
getCmdLineArgs (CmdLineArgs& args, int argc, char* argv[])
{
  Teuchos::CommandLineProcessor cmdp (false, true);
  cmdp.setOption ("mapFilename", &args.mapFilename, "Name of the map "
                  "file for the right-hand side vector(s) B");
  cmdp.setOption ("matrixFilename", &args.matrixFilename, "Name of Matrix "
                  "Market file with the sparse matrix A");
  cmdp.setOption ("rhsFilename", &args.rhsFilename, "Name of Matrix Market "
                  "file with the right-hand side vector(s) B");
  cmdp.setOption ("lineFilename", &args.lineFilename, "Name of Matrix Market "
                  "file with the lineid of each node listed");
  cmdp.setOption ("blockSize", &args.blockSize, "Size of block to use");
  cmdp.setOption ("numIters", &args.numIters, "Number of iterations");
  cmdp.setOption ("tol", &args.tol, "Solver tolerance");
  cmdp.setOption ("nx", &args.nx, "If using inline meshing, number of nodes in the x direction per proc");
  cmdp.setOption ("lpp", &args.lpp, "If using inline meshing, number of lines per proc");
  cmdp.setOption ("withStackedTimer", "withoutStackedTimer", &args.useStackedTimer,
      "Whether to run with a StackedTimer and print the timer tree at the end (and try to output Watchr report)");
  cmdp.setOption("problemName", &args.problemName, "Human-readable problem name for Watchr plot");
  auto result = cmdp.parse (argc, argv);
  return result == Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL;
}

} // namespace (anonymous)

// Xpetra / Galeri
#if defined(HAVE_IFPACK2_XPETRA) 
#include "Xpetra_ConfigDefs.hpp"
#include "Xpetra_DefaultPlatform.hpp"
#include "Xpetra_Parameters.hpp"
#include "Xpetra_MapFactory.hpp"
#include "Xpetra_TpetraMap.hpp"
#include "Xpetra_CrsMatrix.hpp"
#include "Xpetra_TpetraCrsMatrix.hpp"
#include "Galeri_XpetraProblemFactory.hpp"
#include "Galeri_XpetraMatrixTypes.hpp"
#include "Galeri_XpetraParameters.hpp"
#include "Galeri_XpetraUtils.hpp"
#include "Galeri_XpetraMaps.hpp"


// Create a matrix as specified by parameter list options
template <class SC, class LO, class GO, class NO>
static Teuchos::RCP<Xpetra::Matrix<SC,LO,GO,NO> > BuildMatrix(Teuchos::ParameterList &matrixList, Teuchos::RCP<const Teuchos::Comm<int> > & comm) {
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  Xpetra::UnderlyingLib lib = Xpetra::UseTpetra;
  using Map = Xpetra::Map<LO,GO,NO>;
  using Matrix = Xpetra::Matrix<SC,LO,GO,NO>;
  using CrsMatrixWrap = Xpetra::CrsMatrixWrap<SC,LO,GO,NO>;
  using MapFactory = Xpetra::MapFactory<LO,GO,NO>;
  using MultiVector = Xpetra::MultiVector<SC,LO,GO,NO>;
  
  GO nx,ny,nz;
  nx = ny = nz = 5;
  nx = matrixList.get("nx", nx);
  ny = matrixList.get("ny", ny);
  nz = matrixList.get("nz", nz);
  
  std::string matrixType = matrixList.get("matrixType","Laplace1D");
  GO numGlobalElements; //global_size_t
  if (matrixType == "Laplace1D")
    numGlobalElements = nx;
  else if (matrixType == "Laplace2D" || matrixType == "Star2D" || matrixType == "Cross2D")
    numGlobalElements = nx*ny;
  else if(matrixType == "Elasticity2D")
    numGlobalElements = 2*nx*ny;
  else if (matrixType == "Laplace3D" || matrixType == "Brick3D")
    numGlobalElements = nx*ny*nz;
  else if  (matrixType == "Elasticity3D")
    numGlobalElements = 3*nx*ny*nz;
  else {
    std::string msg = matrixType + " is unsupported (in unit testing)";
    throw std::runtime_error(msg);
  }
  
  RCP<const Map> map = MapFactory::Build(lib, numGlobalElements, 0, comm);
  RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr =
    Galeri::Xpetra::BuildProblem<SC,LO,GO,Map,CrsMatrixWrap,MultiVector>(matrixType, map, matrixList);
  RCP<Matrix> Op = Pr->BuildMatrix();
  
  return Op;
} // BuildMatrix()



template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Tpetra::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > BuildBlockMatrix(Teuchos::ParameterList &matrixList,Teuchos:: RCP<const Teuchos::Comm<int> > & comm) {
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;

  // Thanks for the code, Travis!
  
  // Make the graph
  RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > FirstMatrix = BuildMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(matrixList,comm);
  RCP<const Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node> > FGraph = FirstMatrix->getCrsGraph();
  
  int blocksize = 3;
  RCP<const Xpetra::TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node> > TGraph = rcp_dynamic_cast<const Xpetra::TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node> >(FGraph);
  RCP<const Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > TTGraph = TGraph->getTpetra_CrsGraph();
  
  RCP<Tpetra::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > bcrsmatrix = rcp(new Tpetra::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> (*TTGraph, blocksize));
  
  const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>& meshRowMap = *bcrsmatrix->getRowMap();
  const Scalar zero   = Teuchos::ScalarTraits<Scalar>::zero();
  const Scalar one   = Teuchos::ScalarTraits<Scalar>::one();
  const Scalar two   = one+one;
  const Scalar three = two+one;
  
  Teuchos::Array<Scalar> basematrix(blocksize*blocksize, zero);
  basematrix[0] = two;
  basematrix[2] = three;
  basematrix[3] = three;
  basematrix[4] = two;
  basematrix[7] = three;
  basematrix[8] = two;
  Teuchos::Array<Scalar> offmatrix(blocksize*blocksize, zero);
  offmatrix[0]=offmatrix[4]=offmatrix[8]=-1;
  
  Teuchos::Array<LocalOrdinal> lclColInds(1);
  for (LocalOrdinal lclRowInd = meshRowMap.getMinLocalIndex (); lclRowInd <= meshRowMap.getMaxLocalIndex(); ++lclRowInd) {
    lclColInds[0] = lclRowInd;
    bcrsmatrix->replaceLocalValues(lclRowInd, lclColInds.getRawPtr(), &basematrix[0], 1);
    
    // Off diagonals
    if(lclRowInd > meshRowMap.getMinLocalIndex ()) {
      lclColInds[0] = lclRowInd - 1;
      bcrsmatrix->replaceLocalValues(lclRowInd, lclColInds.getRawPtr(), &offmatrix[0], 1);
    }
    if(lclRowInd < meshRowMap.getMaxLocalIndex ()) {
      lclColInds[0] = lclRowInd + 1;
      bcrsmatrix->replaceLocalValues(lclRowInd, lclColInds.getRawPtr(), &offmatrix[0], 1);
    }
    
  }

  return bcrsmatrix; 
} // BuildBlockMatrix()

#endif // HAVE_IFPACK2_XPETRA



int
main (int argc, char* argv[])
{
  using Teuchos::Comm;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::Time;
  using Teuchos::StackedTimer;
  using std::cerr;
  using std::endl;
  typedef Tpetra::CrsMatrix<> crs_matrix_type;
  typedef Tpetra::BlockCrsMatrix<> block_crs_matrix_type;
  typedef Tpetra::Map<> map_type;
  typedef Tpetra::MultiVector<> MV;
  typedef Tpetra::RowMatrix<> row_matrix_type;

  using SC = typename MV::scalar_type;
  using LO = typename MV::local_ordinal_type;
  using GO = typename MV::global_ordinal_type;
  using NO = typename MV::node_type;
  using MT = Teuchos::ScalarTraits<SC>::magnitudeType; 

  typedef Tpetra::Vector<LO,LO,GO,NO> IV;
  typedef Tpetra::MatrixMarket::Reader<crs_matrix_type> reader_type;
  typedef Tpetra::MatrixMarket::Reader<Tpetra::CrsMatrix<LO,LO,GO,NO> > LO_reader_type;
  typedef Ifpack2::BlockTriDiContainer<row_matrix_type> BTDC;

  Tpetra::ScopeGuard tpetraScope (&argc, &argv);

  RCP<const Comm<int> > comm = Tpetra::getDefaultComm ();
  bool rank0 = comm->getRank() == 0;

  // Get command-line arguments.
  CmdLineArgs args;
  const bool gotCmdLineArgs = getCmdLineArgs (args, argc, argv);
  if (! gotCmdLineArgs) {
    if(rank0) cerr << "Failed to get command-line arguments!" << endl;
    return EXIT_FAILURE;
  }

  // If using StackedTimer, then do not time Galeri or matrix I/O because that will dominate the total time.

  RCP<StackedTimer> stackedTimer;
  RCP<Time> totalTime;
  RCP<Teuchos::TimeMonitor> totalTimeMon;
  RCP<Time> precSetupTime = Teuchos::TimeMonitor::getNewTimer ("Preconditioner setup");
  RCP<Time> solveTime = Teuchos::TimeMonitor::getNewTimer ("Solve");
  if(!args.useStackedTimer)
  {
    totalTime = Teuchos::TimeMonitor::getNewTimer ("Total");
    totalTimeMon = rcp(new Teuchos::TimeMonitor(*totalTime));
  }

  bool inline_matrix = false;
#if defined(HAVE_IFPACK2_XPETRA)
  // If we have Xpetra/Galeri, we can use inline matrix generation.  If we're doing that, we
  // also reset everything else to the empty string to make the later code easier
  if (args.matrixFilename == "") {
     args.mapFilename   = "";
     args.rhsFilename  = "";
     args.lineFilename = "";
     args.blockSize = 3; 
     inline_matrix = true;
  }
#endif
  if(inline_matrix == false) {
    if (args.mapFilename == "") {
      if (rank0) cerr << "Must specify filename for loading the map of the right-hand side(s)!" << endl;
      return EXIT_FAILURE;
    }
    if (args.matrixFilename == "") {
      if (rank0) cerr << "Must specify sparse matrix filename!" << endl;
      return EXIT_FAILURE;
    }
    if (args.rhsFilename == "") {
      if (rank0) cerr << "Must specify filename for loading right-hand side(s)!" << endl;
      return EXIT_FAILURE;
    }
    if (args.lineFilename == "") {
      if (rank0) cerr << "Must specify filename for loading line information!" << endl;
      return EXIT_FAILURE;    
    }
    if (args.blockSize <= 0) {
      if (rank0) cerr << "Must specify block size!" << endl;
      return EXIT_FAILURE;    
    }
  }

  // Read sparse matrix A from Matrix Market file.
  RCP<crs_matrix_type> A;
  RCP<block_crs_matrix_type> Ablock;
  RCP<MV> B,X;
  RCP<IV> line_info;
#if defined(HAVE_IFPACK2_XPETRA)
  if(args.matrixFilename == "") {
    // matrix
    Teuchos::ParameterList plist;
    plist.set("matrixType","Laplace1D");    
    plist.set("nx", (GO)args.nx*comm->getSize());
    Ablock = BuildBlockMatrix<SC,LO,GO,NO>(plist,comm);

    //rhs 
    B = rcp(new MV(Ablock->getRangeMap(),1));
    B->putScalar(Teuchos::ScalarTraits<SC>::one());

    // line info (lpp lines per proc)
    int line_length = std::max(1, args.nx  / args.lpp);
    line_info = rcp(new IV(Ablock->getRowMap()));
    auto line_ids = line_info->get1dViewNonConst();
    for(LO i=0; i<(LO)line_ids.size(); i++)
      line_ids[i] = i / line_length;     

    if(rank0)
      std::cout<<"Using block_size = "<<args.blockSize<<" # lines per proc= "<<args.lpp<<" and average line length = "<<line_length<<std::endl;
  }
  else
#endif 
    {
      // Read map
      if(rank0) std::cout<<"Reading map file..."<<std::endl;
      RCP<const map_type> point_map = reader_type::readMapFile(args.mapFilename, comm);
      if(point_map.is_null()) {
        if (rank0) {
          cerr << "Failed to load row map from file "
            "\"" << args.mapFilename << "\"!" << endl;
        }
        return EXIT_FAILURE;
      }

      // Read matrix
      if(rank0) std::cout<<"Reading matrix (as point)..."<<std::endl;
      RCP<const map_type> dummy_col_map;
      if (comm->getSize() == 0)
        A = reader_type::readSparseFile(args.matrixFilename, point_map, dummy_col_map, point_map, point_map);
      else {
        if(rank0) std::cout<<"Using per-rank reader..."<<std::endl;
        A = reader_type::readSparsePerRank(args.matrixFilename, ".mm", point_map, dummy_col_map, point_map, point_map,true,false,8,true);
      }
      if (A.is_null()) {
        if (rank0) {
          cerr << "Failed to load sparse matrix A from file "
            "\"" << args.matrixFilename << "\"!" << endl;
        }
        return EXIT_FAILURE;
      }
      if(rank0) std::cout<<"Matrix read complete..."<<std::endl;


      // Read right-hand side vector(s) B from Matrix Market file.
      if(rank0) std::cout<<"Reading rhs file..."<<std::endl;
      B = reader_type::readDenseFile(args.rhsFilename, comm, point_map);
      if (B.is_null()) {
        if (rank0) {
          cerr << "Failed to load right-hand side vector(s) from file \""
               << args.rhsFilename << "\"!" << endl;
        }
        return EXIT_FAILURE;
      }
      
      // Convert Matrix to Block
      if(rank0) std::cout<<"Converting A from point to block..."<<std::endl;
      Ablock = Tpetra::convertToBlockCrsMatrix<SC,LO,GO,NO>(*A, args.blockSize);


      // Read line information vector
      // We assume the vector contains the local line ids for each node
      if(rank0) std::cout<<"Reading line info file..."<<std::endl;
      RCP<const map_type> block_map = Ablock->getRowMap();
      line_info = LO_reader_type::readVectorFile(args.lineFilename, comm, block_map);
      if (line_info.is_null ()) {
        if (rank0) {
          cerr << "Failed to load line_info from file \""
               << args.lineFilename << "\"!" << endl;
        }
        return EXIT_FAILURE;
      }

    }


  
  if(rank0) {
    size_t numDomains = Ablock->getDomainMap()->getGlobalNumElements();
    size_t numRows = Ablock->getRowMap()->getGlobalNumElements();
    std::cout<<"Block Matrix has "<<numDomains<<" domains and "<<numRows
             << " rows with an implied block size of "<< ((double)numDomains / (double)numRows)<<std::endl;
  }

  if(args.useStackedTimer)
  {
    stackedTimer = rcp(new StackedTimer("BlockTriDiagonalSolver"));
    Teuchos::TimeMonitor::setStackedTimer(stackedTimer);
  }

  // Initial Guess
  if(rank0) std::cout<<"Allocating initial guess..."<<std::endl;
  X = rcp(new MV(Ablock->getRangeMap(),1));
  X->putScalar(Teuchos::ScalarTraits<SC>::zero());
    
  // Initial diagnostics
  Teuchos::Array<MT> normx(1),normb(1); 
  X->norm2(normx);
  B->norm2(normb);
  if(rank0) {
    std::cout<<"Initial norm X = "<<normx[0]<<" norm B = "<<normb[0]<<std::endl;
  }


  // Convert line_info vector to parts arrays
  // NOTE: Both of these needs to be nodes-leve guys, not parts-level.
  Teuchos::Array<Teuchos::Array<LO> > parts;
  {
    if(rank0) std::cout<<"Converting line info to parts..."<<std::endl;
    // Number of lines will vary per proc, so we need to count these
    auto line_ids = line_info->get1dView();
    LO max_line_id = 0;
    for(LO i=0; i<(LO)line_ids.size(); i++)
      max_line_id = std::max(max_line_id,line_ids[i]);

    LO num_local_lines = max_line_id + 1;

    for(LO i=0; i<num_local_lines; i++)
      parts.push_back(Teuchos::Array<LO>());


    // Assume contiguous blocks here
    for(LO i=0; i<(LO)line_ids.size(); i++) {
      LO block_lid = i;
      LO block_num = line_ids[i];
      parts[block_num].push_back(block_lid);     
    }      
    //    std::cout<<"On "<<line_ids.size()<<" local DOFs, detected "<<num_local_lines<<" lines"<<std::endl;
  }

  // Preposition the matrix on device by letting a matvec ensure a transfer
  {
    RCP<MV> temp = rcp(new MV(Ablock->getRangeMap(),1));
    Ablock->apply(*X,*temp);
  }

  // Create Ifpack2 preconditioner.
  if(rank0) std::cout<<"Creating preconditioner..."<<std::endl;
  RCP<BTDC> precond;

  {
    Teuchos::TimeMonitor precSetupTimeMon (*precSetupTime);
    precond = rcp(new BTDC(Ablock,parts));

    if(rank0) std::cout<<"Initializing preconditioner..."<<std::endl;
    precond->initialize ();

    if(rank0) std::cout<<"Computing preconditioner..."<<std::endl;
    precond->compute ();
  }

  // Solver Parameters
  auto ap                 = precond->createDefaultApplyParameters();
  ap.zeroStartingSolution = true;
  ap.tolerance            = args.tol;
  ap.maxNumSweeps         = args.numIters;
  ap.checkToleranceEvery  = 1;
 

  // Solve
  if(rank0) std::cout<<"Running solve..."<<std::endl;
  int nits;
  {
    Teuchos::TimeMonitor solveTimeMon (*solveTime);
    nits = precond->applyInverseJacobi(*B,*X,ap); 
  }

  auto norm0 = precond->getNorms0();
  auto normF = precond->getNormsFinal();

  if(rank0) {
    std::cout<<"Solver run for "<<nits<<" iterations (asked for "<<args.numIters<<") with residual reduction "<<normF/norm0<<std::endl;
    std::cout<<"  Norm0 = "<<norm0<<" NormF = "<<normF<<std::endl;
  }


  X->norm2(normx);
  B->norm2(normb);
  if(rank0) {
    std::cout<<"Final norm X = "<<normx[0]<<" norm B = "<<normb[0]<<std::endl;
  }


  // Report timings.
  if(args.useStackedTimer)
  {
    stackedTimer->stopBaseTimer();
    StackedTimer::OutputOptions options;
    options.num_histogram=3;
    options.print_warnings = false;
    options.output_histogram = true;
    options.output_fraction=true;
    options.output_minmax = true;
    stackedTimer->report(std::cout, comm, options);
    auto xmlOut = stackedTimer->reportWatchrXML(args.problemName, comm);
    if(rank0)
    {
      if(xmlOut.length())
        std::cout << "\nAlso created Watchr performance report " << xmlOut << '\n';
    }
  }
  else
  {
    // Stop the "Total" timer.
    if(!args.useStackedTimer)
      totalTimeMon = Teuchos::null;
    Teuchos::TimeMonitor::report (comm.ptr (), std::cout);
  }

  // Test output if this doesn't crash
  bool success=true;
  /*

  bool verbose=true;
  try {
    ;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);
  */

  comm->barrier();
  if(rank0) std::cout<<"End Result: TEST PASSED"<<std::endl;

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
