// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
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
  CmdLineArgs ():blockSize(-1),numIters(10),numRepeats(1),tol(1e-12),nx(172),ny(-1),nz(-1),mx(1),my(1),mz(1),sublinesPerLine(1),sublinesPerLineSchur(1),useStackedTimer(false),usePointMatrix(false),overlapCommAndComp(false),useSingleFile(false),skipLineFile(false){}

  std::string mapFilename;
  std::string matrixFilename;
  std::string rhsFilename;
  std::string lineFilename;
  int blockSize;
  int numIters;
  int numRepeats;
  double tol;
  int nx;
  int ny;
  int nz;
  int mx;
  int my;
  int mz;
  int sublinesPerLine;
  int sublinesPerLineSchur;
  bool useStackedTimer;
  bool usePointMatrix;
  bool overlapCommAndComp;
  bool useSingleFile;
  bool skipLineFile;
  std::string problemName;
  std::string matrixType;
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
  cmdp.setOption ("numIters", &args.numIters, "Number of iterations per Solve call");
  cmdp.setOption ("numRepeats", &args.numRepeats, "Number of times to run preconditioner compute & solve.");
  cmdp.setOption ("tol", &args.tol, "Solver tolerance");
  cmdp.setOption ("nx", &args.nx, "If using inline meshing, number of nodes in the x direction");
  cmdp.setOption ("ny", &args.ny, "If using inline meshing, number of nodes in the y direction");
  cmdp.setOption ("nz", &args.nz, "If using inline meshing, number of nodes in the z direction");
  cmdp.setOption ("mx", &args.mx, "If using inline meshing, number of procs in the x direction");
  cmdp.setOption ("my", &args.my, "If using inline meshing, number of procs in the y direction");
  cmdp.setOption ("mz", &args.mz, "If using inline meshing, number of procs in the z direction");
  cmdp.setOption ("sublinesPerLine", &args.sublinesPerLine, "If using inline meshing, number of sublines per mesh x line. If set to -1 the block Jacobi algorithm is used.");
  cmdp.setOption ("withStackedTimer", "withoutStackedTimer", &args.useStackedTimer,
      "Whether to run with a StackedTimer and print the timer tree at the end (and try to output Watchr report)");
  cmdp.setOption ("withPointMatrix", "withoutPointMatrix", &args.usePointMatrix,
      "Whether to run with a point matrix");
  cmdp.setOption ("withOverlapCommAndComp", "withoutOverlapCommAndComp", &args.overlapCommAndComp,
		  "Whether to run with overlapCommAndComp");
  cmdp.setOption ("useSingeFile", "useOneFilePerRank", &args.useSingleFile,
		  "Whether to read the matrix from one file or one file per rank");
  cmdp.setOption ("skipLineFile", "readLineFile", &args.skipLineFile,
		  "Whether to skip the lineFile and use a block Jacobi");
  cmdp.setOption("problemName", &args.problemName, "Human-readable problem name for Watchr plot");
  cmdp.setOption("matrixType", &args.matrixType, "matrixType");
  cmdp.setOption("sublinesPerLineSchur", &args.sublinesPerLineSchur, "sublinesPerLineSchur");
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
  RCP<const Map> map;
  if (matrixType == "Laplace1D")
    map = Galeri::Xpetra::CreateMap<LO,GO,NO>(lib, "Cartesian1D", comm, matrixList);
  else if (matrixType == "Laplace2D" || matrixType == "Star2D" || matrixType == "Cross2D")
    map = Galeri::Xpetra::CreateMap<LO,GO,NO>(lib, "Cartesian2D", comm, matrixList);
  else if(matrixType == "Elasticity2D") {
    GO numGlobalElements = 2*nx*ny;
    map = MapFactory::Build(lib, numGlobalElements, 0, comm);
  }
  else if (matrixType == "Laplace3D" || matrixType == "Brick3D")
    map = Galeri::Xpetra::CreateMap<LO,GO,NO>(lib, "Cartesian3D", comm, matrixList);
  else if  (matrixType == "Elasticity3D") {
    GO numGlobalElements = 3*nx*ny*nz;
    map = MapFactory::Build(lib, numGlobalElements, 0, comm);
  }
  else {
    std::string msg = matrixType + " is unsupported (in unit testing)";
    throw std::runtime_error(msg);
  }
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

  // Make the graph
  RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > FirstMatrix = BuildMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>(matrixList,comm);
  RCP<const Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node> > FGraph = FirstMatrix->getCrsGraph();
  
  const int blocksize = matrixList.get("blockSize", 3);
  RCP<const Xpetra::TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node> > TGraph = rcp_dynamic_cast<const Xpetra::TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node> >(FGraph);
  RCP<const Tpetra::CrsGraph<LocalOrdinal,GlobalOrdinal,Node> > TTGraph = TGraph->getTpetra_CrsGraph();
  
  using BCRS = Tpetra::BlockCrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  RCP<BCRS> bcrsmatrix = rcp(new BCRS (*TTGraph, blocksize));
  
  const Scalar zero   = Teuchos::ScalarTraits<Scalar>::zero();
  const Scalar one   = Teuchos::ScalarTraits<Scalar>::one();
  const Scalar two   = one+one;
  const Scalar three = two+one;
  
  Teuchos::Array<Scalar> basematrix(blocksize*blocksize, zero);
  Teuchos::Array<Scalar> offmatrix(blocksize*blocksize, zero);
  int entryIdx = 0;
  for(int rowIdx = 0; rowIdx < blocksize; ++rowIdx) {
    for(int colIdx = 0; colIdx < blocksize; ++colIdx) {
      entryIdx = rowIdx*blocksize + colIdx;
      basematrix[entryIdx] = ((entryIdx % 2) == 0 ? one : two);
    }
    // We enforce that the block matrix is SPD by setting the diagonal value
    // to a value which is greater than the sum of the absolute value of all the
    // entries on the corresponding row knowing that the off diagonal blocks are
    // set to minus identity.
    basematrix[rowIdx*blocksize + rowIdx] = three + 26 * one + blocksize * two;

    offmatrix[rowIdx*blocksize + rowIdx] = -one;
  }

  const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>& meshRowMap = *bcrsmatrix->getRowMap();
  const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>& meshColMap = *bcrsmatrix->getColMap();

  LocalOrdinal localNumRows = bcrsmatrix->getLocalNumRows();

#ifdef VERBOSE_OUTPUT
  std::stringstream streamI, streamJ, streamV;

  streamI << "I" << comm->getRank() << " = array([ ";
  streamJ << "J" << comm->getRank() << " = array([ ";
  streamV << "V" << comm->getRank() << " = array([ ";
  bool firstPrint = true;

#endif
  for (LocalOrdinal localRowInd = 0; localRowInd < localNumRows; ++localRowInd) {
    // Get a view of the current row.
    // You may modify the values, but not the column indices.
    typename BCRS::local_inds_host_view_type localColInds;
    typename BCRS::nonconst_values_host_view_type vals;

    LocalOrdinal numEntries = bcrsmatrix->getNumEntriesInLocalRow(localRowInd);
    if (numEntries == 0) {
      continue;
    }
    bcrsmatrix->getLocalRowViewNonConst (localRowInd, localColInds, vals);
    // Modify the entries in the current row.
    for (LocalOrdinal k = 0; k < numEntries; ++k) {
      LocalOrdinal offset = blocksize * blocksize * k;
#ifdef VERBOSE_OUTPUT
      if (bcrsmatrix->getGlobalNumRows() <= 100) {
        if (!firstPrint) {
          streamI << ", ";
          streamJ << ", ";
          streamV << ", ";
        } else {
          firstPrint = false;
        }
        streamI << meshRowMap.getGlobalElement(localRowInd);
        streamJ << meshColMap.getGlobalElement(localColInds(k));
        streamV << "1.";
      }
#endif
      if ( meshRowMap.getGlobalElement(localRowInd) == meshColMap.getGlobalElement(localColInds(k)) ) {
        // Blocks are stored in row-major format.
        for (LocalOrdinal j = 0; j < blocksize; ++j) {
          for (LocalOrdinal i = 0; i < blocksize; ++i) {
            vals(offset + i + j * blocksize) = basematrix[i + j * blocksize];
          }
        }
      } else {
        // Blocks are stored in row-major format.
        for (LocalOrdinal j = 0; j < blocksize; ++j) {
          for (LocalOrdinal i = 0; i < blocksize; ++i) {
            vals(offset + i + j * blocksize) = offmatrix[i + j * blocksize];
          }
        }
      }
    }
  }

#ifdef VERBOSE_OUTPUT
  if (bcrsmatrix->getGlobalNumRows() <= 100) {
    streamI << "])";
    streamJ << "])";
    streamV << "])";

    std::ofstream graph_file( "log_graph_"+std::to_string(comm->getRank())+".txt" );
    
    graph_file << streamI.str() << std::endl;
    graph_file << streamJ.str() << std::endl;
    graph_file << streamV.str() << std::endl;

    graph_file << "A" << comm->getRank() << " = sparse.coo_matrix((V" << comm->getRank() << ",(I" << comm->getRank() << ",J" << comm->getRank() << ")),shape=("<< bcrsmatrix->getGlobalNumRows() <<","<< bcrsmatrix->getGlobalNumCols() <<")).tocsr()" << std::endl;

    graph_file.close();
  }
#endif

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
  RCP<Time> precComputeTime = Teuchos::TimeMonitor::getNewTimer ("Preconditioner compute");
  RCP<Time> solveTime = Teuchos::TimeMonitor::getNewTimer ("Solve");
  RCP<Time> normTime = Teuchos::TimeMonitor::getNewTimer ("Norm");
  RCP<Time> warmupMatrixApplyTime = Teuchos::TimeMonitor::getNewTimer ("Preposition of the matrix on device");
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
    if (!args.skipLineFile && args.lineFilename == "") {
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
  RCP<row_matrix_type> Ablock;
  RCP<MV> B,X;
  RCP<IV> line_info;
#if defined(HAVE_IFPACK2_XPETRA)
  if(args.matrixFilename == "") {
    RCP<Time> matrixCreationTime = Teuchos::TimeMonitor::getNewTimer ("Create inline matrix");
    Teuchos::TimeMonitor matrixCreationTimeMon (*matrixCreationTime);
    if (args.usePointMatrix) {
      std::string msg = "usePointMatrix with inline matrix is not yet implemented";
      throw std::runtime_error(msg);      
    }
    // matrix
    Teuchos::ParameterList plist;
    if(args.matrixType == "") {
      plist.set("matrixType","Laplace1D");
    } else {
      plist.set("matrixType", args.matrixType);
      plist.set("nx", (GO)args.nx);
      plist.set("mx", (GO)args.mx);
      if(args.mx == 1 && args.ny == -1)
        plist.set("mx", (GO)comm->getSize());
      else if(args.mx != comm->getSize() && args.ny == -1) {
        std::string msg = "mx is not consistent with comm->getSize";
        throw std::runtime_error(msg);
      }
      else if(args.mx != 1) {
        std::string msg = "mx != 1 is not yet supported";
        throw std::runtime_error(msg);
      }
      if(args.ny > -1) {
        plist.set("ny", (GO)args.ny);
        plist.set("my", (GO)args.my);
        if(args.mx*args.my == 1 && args.nz == -1)
          plist.set("my", (GO)comm->getSize());
        else if(args.mx*args.my != comm->getSize() && args.nz == -1) {
          std::string msg = "mx and my are not consistent with comm->getSize";
          throw std::runtime_error(msg);
        }
      }
      if(args.nz > -1) {
        plist.set("nz", (GO)args.nz);
        plist.set("mz", (GO)args.mz);
        if(args.mx*args.my*args.mz == 1)
          plist.set("mz", (GO)comm->getSize());
        else if(args.mx*args.my*args.mz != comm->getSize()) {
          std::string msg = "mx, my, and mz are not consistent with comm->getSize";
          throw std::runtime_error(msg);
        }
      }
    }
    if(0 < args.blockSize) {
      plist.set("blockSize", args.blockSize);
    }
    Ablock = BuildBlockMatrix<SC,LO,GO,NO>(plist,comm);
    std::cout << "p=" << comm->getRank() << " | Ablock, local size: " << Ablock->getLocalNumRows() << std::endl;

    //rhs 
    B = rcp(new MV(Ablock->getRangeMap(),1));
    B->putScalar(Teuchos::ScalarTraits<SC>::one());

    // line info (sublinesPerLine lines per proc along direction x)
    if ( args.sublinesPerLine < 1 && args.sublinesPerLine != -1) {
      std::string msg = "the value of sublinesPerLine = " + std::to_string(args.sublinesPerLine) + " is not supported";
      throw std::runtime_error(msg);
    }

    int line_length = std::max(1, (int) std::ceil(args.nx  / args.sublinesPerLine));
    // We compute the number of lines oriented along the x direction of the mesh.
    // This number is called line_per_x_fiber where a fiber refers to an initial
    // x line in the mesh before dividing it in sublines.
    int line_per_x_fiber = std::ceil(args.nx  / line_length);
    line_info = rcp(new IV(Ablock->getRowMap()));
    auto line_ids = line_info->get1dViewNonConst();
    std::cout << "line_ids.size()=" << line_ids.size()
	      << ", args.nx*args.ny*args.nz/comm->getSize()=" << args.nx*args.ny*args.nz/comm->getSize() << std::endl;
    for(LO i=0; i<(LO)line_ids.size(); i++) {
      LO fiber_id = std::floor(i / args.nx);
      line_ids[i] = line_per_x_fiber * fiber_id + std::floor( (i % args.nx) / line_length );
    }

    if(rank0) {
      std::cout << "Using matrixType = " << plist.get<std::string>("matrixType")
		<< " nx = " << plist.get<GO>("nx")
		<< " ny = " << plist.get<GO>("ny")
		<< " nz = " << plist.get<GO>("nz")
		<< std::endl;
      std::cout<< "Using block_size = " << args.blockSize
	       << " # lines per proc (input provided by the user) = " << args.sublinesPerLine //*args.ny*args.nz
         << " # lines per fiber = " << line_per_x_fiber
	       << " and average line length = " << line_length<<std::endl;
    }
  }
  else
#endif 
    {
      RCP<Time> matrixReadingTime = Teuchos::TimeMonitor::getNewTimer ("Reading matrix input files");
      Teuchos::TimeMonitor matrixReadingTimeMon (*matrixReadingTime);
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
      if (args.useSingleFile || comm->getSize() == 1)
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
      {
        RCP<Time> matrixConversionTime = Teuchos::TimeMonitor::getNewTimer ("Matrix conversion");
        Teuchos::TimeMonitor matrixConversionTimeMon (*matrixConversionTime);
        Ablock = Tpetra::convertToBlockCrsMatrix<SC,LO,GO,NO>(*A, args.blockSize);
      }


      if(args.skipLineFile) {
        line_info = rcp(new IV(Ablock->getRowMap()));
        auto line_ids = line_info->get1dViewNonConst();
        for(LO i=0; i<(LO)line_ids.size(); i++) {
          line_ids[i] = i;
        }
      }
      else {
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

    }


  
  if(rank0) {
    size_t numDomains = Ablock->getDomainMap()->getGlobalNumElements();
    size_t numRows = Ablock->getRowMap()->getGlobalNumElements();
    std::cout<<"Block Matrix has "<<numDomains<<" domains and "<<numRows
             << " rows with an implied block size of "<< ((double)numDomains / (double)numRows)<<std::endl;
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
    Teuchos::TimeMonitor warmupMatrixApplyTimeMon (*warmupMatrixApplyTime);

    RCP<MV> temp = rcp(new MV(Ablock->getRangeMap(),1));
    Ablock->apply(*X,*temp);
  }

  if(args.useStackedTimer)
  {
    stackedTimer = rcp(new StackedTimer("BlockTriDiagonalSolver"));
    Teuchos::TimeMonitor::setStackedTimer(stackedTimer);
  }

  // Create Ifpack2 preconditioner.
  if(rank0) std::cout<<"Creating preconditioner..."<<std::endl;
  RCP<BTDC> precond;

  {
    Teuchos::TimeMonitor precSetupTimeMon (*precSetupTime);
    if(args.usePointMatrix)
      precond = rcp(new BTDC(A,parts,args.sublinesPerLineSchur,args.overlapCommAndComp, false, args.blockSize));
    else
      precond = rcp(new BTDC(Ablock,parts,args.sublinesPerLineSchur,args.overlapCommAndComp));

    if(rank0) std::cout<<"Initializing preconditioner..."<<std::endl;
    precond->initialize ();
    Kokkos::DefaultExecutionSpace().fence();
  }

  // Solver Parameters
  auto ap                 = precond->createDefaultApplyParameters();
  ap.zeroStartingSolution = true;
  ap.tolerance            = args.tol;
  ap.maxNumSweeps         = args.numIters;
  ap.checkToleranceEvery  = 10;
 

  // Solve
  for(int repeat=0; repeat < args.numRepeats; ++repeat)
  {
    if(rank0) std::cout<<"Computing preconditioner..."<<std::endl;
    {
      Teuchos::TimeMonitor precComputeTimeMon (*precComputeTime);
      precond->compute ();
      Kokkos::DefaultExecutionSpace().fence();
    }

    if(rank0) std::cout<<"Running solve..."<<std::endl;
    int nits;
    {
      Teuchos::TimeMonitor solveTimeMon (*solveTime);
      nits = precond->applyInverseJacobi(*B,*X,ap);
      Kokkos::DefaultExecutionSpace().fence();
    }

    auto norm0 = precond->getNorms0();
    auto normF = precond->getNormsFinal();

    if(rank0) {
      std::cout<<"Solver run for "<<nits<<" iterations (asked for "<<args.numIters<<") with residual reduction "<<normF/norm0<<std::endl;
      std::cout<<"  Norm0 = "<<norm0<<" NormF = "<<normF<<std::endl;
    }

    {
      Teuchos::TimeMonitor normTimeMon (*normTime);
      X->norm2(normx);
      B->norm2(normb);
      Kokkos::DefaultExecutionSpace().fence();
    }
    if(rank0) {
      std::cout<<"Final norm X = "<<normx[0]<<" norm B = "<<normb[0]<<std::endl;
    }
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
