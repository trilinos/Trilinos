//@HEADER
// ************************************************************************
//
//               Isorropia: Partitioning and Load Balancing Package
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Alan Williams (william@sandia.gov)
//                 or Erik Boman    (egboman@sandia.gov)
//
// ************************************************************************
//@HEADER

// Read in a matrix market file.  Use isorropia to do graph
// or hypergraph partitioning.  Compute the graph or hypergraph
// balance and cut metrics both before and after partitioning.
//
// This tests Isorropia::Epetra::create_partitioner followed by
// redistribution with a Isorropia::Epetra::Redistributor.
//
// For graph partitioning:
//
// The nonzeros of this matrix represent graph edges.  The row
// or column IDs represent the graph vertices.  Only square
// matrices will be processed with graph partitioning.
//
// Isorropia will ignore the self-edges (the nonzero diagonal entries).
//
// For hypergraph partitioning:
//
// By convention, the columns of this matrix are hyperedges, and we
// wish to balance the vertices, represented by the rows.
//
// If run with --v option, prints out partitioning before and after.
//
// This is what simple.mtx looks like.  25 rows, 25 cols, 105 non-zeroes.
//
//  0123456789012345678901234
// 0xx   x                   0
// 1xxx   x                  1
// 2 xxx   x                 2
// 3  xxx   x                3
// 4   xx    x               4
// 5x    xx   x              5
// 6 x   xxx   x             6
// 7  x   xxx   x            7
// 8   x   xxx   x           8
// 9    x   xx    x          9
// 0     x    xx   x         0
// 1      x   xxx   x        1
// 2       x   xxx   x       2
// 3        x   xxx   x      3
// 4         x   xx    x     4
// 5          x    xx   x    5
// 6           x   xxx   x   6
// 7            x   xxx   x  7
// 8             x   xxx   x 8
// 9              x   xx    x9
// 0               x    xx   0
// 1                x   xxx  1
// 2                 x   xxx 2
// 3                  x   xxx3
// 4                   x   xx4
//  0123456789012345678901234
//
// If run with --f={filename} a matrix market file other than simple.mtx
// will be processed.
//

#include <Isorropia_ConfigDefs.hpp>
#include <Isorropia_Epetra.hpp>
#include <Isorropia_EpetraColorer.hpp>
#include <Isorropia_EpetraCostDescriber.hpp>

#include <ispatest_lbeval_utils.hpp>
#include <ispatest_epetra_utils.hpp>

#ifdef HAVE_EPETRA
#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#ifdef HAVE_EPETRAEXT
#include <EpetraExt_CrsMatrixIn.h>
#include <Epetra_MapColoring.h>
#endif
#endif

// Only run the first test
//#define SHORT_TEST 1

#include <Teuchos_CommandLineProcessor.hpp>

#define GRAPH_PARTITIONING            1
#define HYPERGRAPH_PARTITIONING       2
#define NO_ZOLTAN                     3

#define SUPPLY_EQUAL_WEIGHTS               1
#define SUPPLY_UNEQUAL_WEIGHTS             2
#define NO_APPLICATION_SUPPLIED_WEIGHTS    3

#define EPETRA_CRSGRAPH               1
#define EPETRA_CRSMATRIX              2

#define FAILED() do { failures++; if (!runAll) goto Report; fail = 0;} while (0)

#define ERROREXIT(v, s) \
  do { if (v){               \
    test_type(objectType); \
    std::cout << s << std::endl << "FAIL" << std::endl; \
  }                     \
  exit(1); } while (0)

#define ERRORRETURN(v, s) \
  do { if (v){                 \
    test_type(objectType); \
    std::cout << s << std::endl << "FAIL" << std::endl; \
  }                       \
  return 1; } while(0)

static void test_type(int objectType)
{
  std::cout << "TEST: ";
  if (objectType == EPETRA_CRSGRAPH)
    std::cout << "using Epetra_CrsGraph interface";
  else
    std::cout << "using Epetra_CrsMatrix interface";

  std::cout << std::endl;

}

static int run_test(Teuchos::RCP<Epetra_CrsMatrix> matrix,
	  bool verbose,           // display the graph before & after
	  int objectType)         // use isorropia's CrsMatrix or CrsGraph
{
  int rc=0, fail = 0;
#ifdef HAVE_EPETRAEXT
  int localProc = 0;
  int numProcs = 1;
  int keepDenseEdges = 0;
  bool valid;

#ifdef HAVE_MPI
  const Epetra_MpiComm &Comm = dynamic_cast<const Epetra_MpiComm &>(matrix->Comm());
  localProc = Comm.MyPID();
  numProcs = Comm.NumProc();
#else
  const Epetra_SerialComm &Comm = dynamic_cast<const Epetra_SerialComm &>(matrix->Comm());
#endif

  int numRows = matrix->NumGlobalRows();

  if (numRows < (numProcs * 100)){
    // If dense edges are thrown out of a small
    // matrix, there may be nothing left.
    keepDenseEdges = 1;
  }

  double myShareBefore = 1.0 / numProcs;
  double myShare = myShareBefore;


  // Check that input matrix is valid.  This test constructs an "x"
  // with the matrix->DomainMap() and a "y" with matrix->RangeMap()
  // and then calculates y = Ax.

  valid = ispatest::test_matrix_vector_multiply(*matrix);

  if (!valid){
    ERRORRETURN((localProc==0), "Test matrix is not a valid Epetra matrix");
  }

  // Compute vertex and edge weights

  Teuchos::RCP<Isorropia::Epetra::CostDescriber> costs =
    Teuchos::rcp(new Isorropia::Epetra::CostDescriber);

  Teuchos::RCP<Epetra_Vector> vptr;

  Teuchos::RCP<Epetra_CrsMatrix> eptr;

  Teuchos::RCP<Epetra_Vector> hyperEdgeWeights;

  Teuchos::ParameterList params;


#ifdef HAVE_ISORROPIA_ZOLTAN

  // Set the Zoltan parameters for this problem

  Teuchos::ParameterList &sublist = params.sublist("ZOLTAN");

    //sublist.set("DEBUG_LEVEL", "1"); // Zoltan will print out parameters
    //sublist.set("DEBUG_LEVEL", "5");   // proc 0 will trace Zoltan calls
    //sublist.set("DEBUG_MEMORY", "2");  // Zoltan will trace alloc & free

  sublist.set("LB_METHOD", "GRAPH");

  // Perform partitioning with Zoltan (if we have it)

  Teuchos::RCP<Isorropia::Epetra::Colorer> colorer;

  if (objectType == EPETRA_CRSGRAPH){
    // Test the Epetra_CrsGraph interface of Isorropia
    Teuchos::RCP<const Epetra_CrsGraph> graph =
      Teuchos::rcp(new Epetra_CrsGraph(matrix->Graph()));
    colorer = Teuchos::rcp(new Isorropia::Epetra::Colorer(graph, params));
  }
  else{
      // Test the Epetra_CrsMatrix interface of Isorropia
    Teuchos::RCP<const Epetra_RowMatrix> rm = matrix;
    colorer = Teuchos::rcp(new Isorropia::Epetra::Colorer(rm, params));
  }

  colorer->color();

  int numberColors;
  numberColors = colorer->numColors();

  if (verbose && (localProc == 0)){
    std::cout << "Max number of colors :" << numberColors  << std::endl;
  }

  if ((numberColors < 0) || (numberColors > matrix->NumGlobalRows()))
    ERRORRETURN(verbose, "Inconsistant number of colors : " + numberColors);

#ifdef HAVE_EPETRAEXT
  Teuchos::RefCountPtr<Epetra_MapColoring> colorMap = colorer->generateMapColoring();
  int numberColorsExt;

  numberColorsExt = colorMap->MaxNumColors();

  if (numberColorsExt >= 10)
    ERRORRETURN(verbose, "Too many colors");

  if (numberColorsExt != numberColors)
    ERRORRETURN(verbose, "Inconsistant number of colors");
#endif /* HAVE_EPETRAEXT */

  for (int i = 1 ; i <= numberColors ; i++ ) {
    int numberElems;

    numberElems = colorer->numElemsWithColor(i);
    if (verbose && (localProc == 0)){
      std::cout << "Elems with color " << i << " : " << numberElems  << std::endl;
    }
    if ((numberElems < 0) || (numberElems > matrix->NumMyRows()))
      ERRORRETURN(verbose, "Inconsistant number of elements for color " + i);

    int *currentColor = new int[numberElems];
    colorer->elemsWithColor(i, currentColor, numberElems);
    for (int j =0 ; j < numberElems ; j++) {
      if ((currentColor[j]<0) || (currentColor[j]>= matrix->NumMyRows()) || 
	  (*colorer)[currentColor[j]] != i) {
	ERRORRETURN(verbose, "Inconsistant elements" << currentColor[j] << " for color " <<  i);
	delete[] currentColor;
      }
    }
    delete[] currentColor;
  }

  if (colorer->numElemsWithColor(numberColors + 1) != 0)
      ERRORRETURN(verbose, "Inconsistant number of elements for non existant color ");


  for (int i = 0 ; i < matrix->NumMyRows() ; i++ ) {
    if (verbose && (localProc == 0)){
      std::cout << " color[" << i << "] = " << (*colorer)[i]  << std::endl;
    }
  }
#endif /* HAVE_ISORROPIA_ZOLTAN */

#else
  std::<< "test_simple : currently can only test "
	 << "with Epetra and EpetraExt enabled." << std::endl;
  fail =  1;
#endif

  return (0);
}

int main(int argc, char** argv) {

  int rc=0, fail = 0;
#ifdef HAVE_EPETRAEXT
  bool verbose = false;
  int numProcs = 1;
  int localProc = 0;

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &localProc);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
  const Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  const Epetra_SerialComm Comm;
#endif

  if (getenv("DEBUGME")){
    std::cerr << localProc << " gdb test_simple.exe " << getpid() << std::endl;
    sleep(15);
  }

  Teuchos::CommandLineProcessor clp(false,true);

  // --f=fileName provides a different matrix market file for input
  // --v will print out the partitioning (small files only)
  // --run-all will continue to run all tests even if there is a failure

  std::string *inputFile = new std::string("simple.mtx");
  bool runAll = false;

  clp.setOption( "f", inputFile,
		"Name of input matrix market file");
  clp.setOption( "run-all", "abort", &runAll,
		"Don't abort if one test fails, run all of them.");
  clp.setOption( "v", "q", &verbose,
		"Display matrix before and after partitioning.");

  Teuchos::CommandLineProcessor::EParseCommandLineReturn parse_return =
    clp.parse(argc,argv);

  if( parse_return == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED){
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 0;
  }
  if( parse_return != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL ) {
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    return 1;
  }

  const char *fname = inputFile->c_str();

  // Read in the matrix market file and distribute its rows across the
  // processes.
  //
  // This reader uses the default Epetra_Map for number of rows for the
  // RowMap() and for the RangeMap().  For non-square matrices it uses
  // the default Epetra_Map for the number of columns for the DomainMap(),
  // otherwise it uses the RowMap().
  //
  // The maps can be specified with other versions of MMFtoCrsMatrix().

  Epetra_CrsMatrix *matrixPtr;
  rc = EpetraExt::MatrixMarketFileToCrsMatrix(fname, Comm, matrixPtr);
  if (rc < 0){
    if (localProc==0) std::cerr << "error reading input file"<< std::cout;
    exit(1);
  }

  bool square = (matrixPtr->NumGlobalRows() == matrixPtr->NumGlobalCols());

  // If matrix is square, determine if it's symmetric  TODO

  // Run some partitioning tests
  //   Test graph and hypergraph partitioning
  //   Test with and without application supplied weights
  //   Test the Epetra_CrsMatrix interface and also the Epetra_CrsGraph interface
  //   Do tests where the vertex or edge weights vary widely

  Teuchos::RCP<Epetra_CrsMatrix> testm = Teuchos::rcp(matrixPtr);

  int failures = 0;

  if (square){
#ifdef HAVE_ISORROPIA_ZOLTAN
    fail = run_test(testm,
	       verbose,            // draw graph before and after partitioning
	       EPETRA_CRSMATRIX);       // use the Epetra_CrsMatrix interface

    if (fail) FAILED();



  fail = run_test(testm,
	     verbose,
	     EPETRA_CRSMATRIX);

  if (fail) FAILED();
#endif

  }

#else
  fail = 0;
  if (localProc == 0){
    std::cout << "Test not run because it requires EPETRA_EXT" << std::endl;
  }
#endif

Report:

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  if (localProc == 0){
    if (failures){
      if (failures > 1)
	std::cout << std::endl << failures << " FAILURES" << std::endl;
      else
	std::cout << std::endl << "1 FAILURE" << std::endl;

      if (!runAll){
	std::cout <<
       "(Use option --run-all if you do not want this test to abort on failure)" << std::endl;
      }

    }
    else
      std::cout << std::endl << "PASS" << std::endl;
  }

  return fail;
}
