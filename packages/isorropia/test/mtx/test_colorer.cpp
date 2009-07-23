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
//
// ************************************************************************
//@HEADER

// Read in a matrix market file.  Use Isorropia to do graph coloring.
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
#ifdef HAVE_EPETRAEXT
  int localProc = 0;
  int numProcs = 1;
  int keepDenseEdges = 0;
  bool valid;

#ifdef HAVE_MPI
  const Epetra_MpiComm &Comm = dynamic_cast<const Epetra_MpiComm &>(matrix->Comm());
  localProc = Comm.MyPID();
  numProcs = Comm.NumProc();
#endif

  int numRows = matrix->NumGlobalRows();

  if (numRows < (numProcs * 100)){
    // If dense edges are thrown out of a small
    // matrix, there may be nothing left.
    keepDenseEdges = 1;
  }


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
    ERRORRETURN(verbose, "Inconsistent number of colors : " + numberColors);


#ifdef HAVE_EPETRAEXT
  Teuchos::RCP<Epetra_MapColoring> colorMap = colorer->generateRowMapColoring();
  int numberColorsExt;

  numberColorsExt = colorMap->MaxNumColors();


  colorMap = colorer->generateColMapColoring();
#endif /* HAVE_EPETRAEXT */

  for (int i = 1 ; i <= numberColors ; i++ ) {
    int numberElems;

    numberElems = colorer->numElemsWithColor(i);
    if (verbose && (localProc == 0)){
      std::cout << "(" << localProc << ") Elems with color " << i << " : " << numberElems  << std::endl;
    }
    if ((numberElems < 0) || (numberElems > matrix->NumMyRows()))
      ERRORRETURN(verbose, "Inconsistent number of elements for color " + i);

    int *currentColor = new int[numberElems];
    colorer->elemsWithColor(i, currentColor, numberElems);
    for (int j =0 ; j < numberElems ; j++) {
      if ((currentColor[j]<0) || (currentColor[j]>= matrix->NumMyRows()) || 
	  (*colorer)[currentColor[j]] != i) {
	ERRORRETURN(verbose, "Inconsistent elements " << currentColor[j] << " for color " <<  i);
	delete[] currentColor;
      }
    }
    delete[] currentColor;
  }

  if (colorer->numElemsWithColor(numberColors + 1) != 0)
      ERRORRETURN(verbose, "Inconsistent number of elements for non existant color ");


  int* colortab;
  int colorsize;
  colorer->extractColorsView(colorsize, (const int*&)colortab);
  if (colorsize != matrix->NumMyRows())
    ERRORRETURN(verbose, "Inconsistent colortab size (1)");

  for (int i = 0 ; i < matrix->NumMyRows() ; i++ ) {
    if (colortab[i] != (*colorer)[i])
      ERRORRETURN(verbose, "Inconsistent colortab (1)");
  }

  colortab = new int[colorsize];
  int len = (colorsize>2)?(colorsize-2):colorsize;
  colorer->extractColorsCopy(len, colorsize, colortab);
  if (colorsize != len) {
    delete[] colortab;
    ERRORRETURN(verbose, "Inconsistent colortab size (2)");
  }

  for (int i = 0 ; i < len ; i++ ) {
    if (colortab[i] != (*colorer)[i]) {
      delete[] colortab;
      ERRORRETURN(verbose, "Inconsistent colortab (2)");
    }
  }
  delete[] colortab;

  for (int i = 0 ; i < matrix->NumMyRows() ; i++ ) {
    if (verbose && (localProc == 0)){
      std::cout << " color[" << i << "] = " << (*colorer)[i]  << std::endl;
    }
  }
#endif /* HAVE_ISORROPIA_ZOLTAN */

#else
  std::cerr<< "test_simple : currently can only test "
	 << "with Epetra and EpetraExt enabled." << std::endl;
  fail =  1;
#endif

  return (0);
}

int main(int argc, char** argv) {

  int rc=0, fail = 0;
  int localProc = 0;
  int failures = 0;
  bool verbose = false;

#ifdef HAVE_MPI
  int numProcs;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &localProc);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
  const Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else /* HAVE_MPI */
  const Epetra_SerialComm Comm;
#endif /* HAVE_MPI */

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


  if (square){
#ifdef HAVE_ISORROPIA_ZOLTAN
    fail = run_test(testm,
	       verbose,            // draw graph before and after partitioning
	       EPETRA_CRSMATRIX);       // use the Epetra_CrsMatrix interface

    if (fail) FAILED();



  fail = run_test(testm,
	     verbose,
	     EPETRA_CRSGRAPH);

  if (fail) FAILED();
#endif

  }

Report:

#ifdef HAVE_MPI
  MPI_Finalize();
#endif /* HAVE_MPI */

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
