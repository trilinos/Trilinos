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
//
// Read in a matrix market file.  Use isorropia to do graph
// or hypergraph partitioning.  Compute the graph or hypergraph
// balance and cut metrics both before and after partitioning.
//
// This tests all variants of Isorropia::Epetra::create_balanced_copy().
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
#include <Isorropia_EpetraPartitioner.hpp>
#include <Isorropia_EpetraRedistributor.hpp>
#include <Isorropia_EpetraCostDescriber.hpp>

#include <ispatest_lbeval_utils.hpp>
#include <ispatest_epetra_utils.hpp>

static int tmp=0;
#define CHECK_FAILED() {      \
  Comm.SumAll(&fail, &tmp, 1);      \
  if (tmp){     \
    failures++;      \
    if (!runAll) goto Report;      \
  }     \
  fail = 0;     \
  }

#ifdef HAVE_EPETRA
#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#include <Epetra_Import.h>
#ifdef HAVE_EPETRAEXT
#include <EpetraExt_CrsMatrixIn.h>
#endif

#endif

#include <Teuchos_CommandLineProcessor.hpp>

#include <unistd.h>

//#define SHORT_TEST

#define GRAPH_PARTITIONING            1
#define HYPERGRAPH_PARTITIONING       2
#define NO_ZOLTAN                     3

#define SUPPLY_EQUAL_WEIGHTS               1
#define SUPPLY_UNEQUAL_WEIGHTS             2
#define NO_APPLICATION_SUPPLIED_WEIGHTS    3

#define EPETRA_CRSGRAPH               1
#define EPETRA_CRSMATRIX              2
#define EPETRA_ROWMATRIX              3
#define EPETRA_LINEARPROBLEM          4

#define ERROREXIT(v, s) \
  if (v){               \
    test_type(numPartitions, partitioningType, vertexWeightType, edgeWeightType, objectType); \
    std::cout << s << std::endl << "FAIL" << std::endl; \
  }                     \
  exit(1);

#define ERRORRETURN(v, s) \
  if (v){                 \
    test_type(numPartitions, partitioningType, vertexWeightType, edgeWeightType, objectType); \
    std::cout << s << std::endl << "FAIL" << std::endl; \
  }                       \
  return 1;


void test_type(int numPartitions, int partitioningType, int vertexWeightType,
	  int edgeWeightType, int objectType)
{
  std::cout << "TEST: ";
  if (partitioningType == NO_ZOLTAN)
    std::cout << "isorropia simple linear row partitioner, ";
  else if (partitioningType == GRAPH_PARTITIONING)
    std::cout << "graph partitioning, ";
  else
    std::cout << "hypergraph partitioning, ";

  if (vertexWeightType == SUPPLY_EQUAL_WEIGHTS)
    std::cout << "created equal vertex (row) weights, ";
  else if (vertexWeightType == SUPPLY_UNEQUAL_WEIGHTS)
    std::cout << "created unequal vertex (row) weights, ";
  else
    if (partitioningType == NO_ZOLTAN){
      std::cout << std::endl << "      ";
      std::cout << "did not supply row weights, default is number of row nonzeros,";
    }
    else
      std::cout << "did not supply vertex (row) weights, ";

  std::cout << std::endl << "      ";

  if (partitioningType != NO_ZOLTAN){

    if (partitioningType == GRAPH_PARTITIONING){
      if (edgeWeightType == SUPPLY_EQUAL_WEIGHTS)
	std::cout << "created equal edge weights, ";
      else if (edgeWeightType == SUPPLY_UNEQUAL_WEIGHTS)
	std::cout << "created unequal edge weights, ";
      else
	std::cout << "did not supply edge weights, ";
    }
    else{
      if (edgeWeightType == SUPPLY_EQUAL_WEIGHTS)
	std::cout << "created equal hyperedge (column) weights, ";
      else if (edgeWeightType == SUPPLY_UNEQUAL_WEIGHTS)
	std::cout << "created unequal hyperedge (column) weights,";
      else
	std::cout << "did not supply hyperedge (column) weights, ";
    }
  }

  if (objectType == EPETRA_CRSGRAPH)
    std::cout << "using Epetra_CrsGraph interface";
  else if (objectType == EPETRA_CRSMATRIX)
    std::cout << "using Epetra_CrsMatrix interface";
  else if (objectType == EPETRA_LINEARPROBLEM)
    std::cout << "using Epetra_LinearProblem interface";
  else if (objectType == EPETRA_ROWMATRIX)
    std::cout << "using Epetra_RowMatrix interface";

  if (numPartitions > 0){
    std::cout << std::endl << "      ";
    std::cout << "NUM_GLOBAL_PARTS is " << numPartitions;
  }
  std::cout << std::endl;

}

static int run_test(Teuchos::RCP<Epetra_CrsMatrix> matrix,
	  bool verbose,           // display the graph before & after
	  bool contract,          // set global number of partitions to 1/2 num procs
	  int partitioningType,   // hypergraph or graph partitioning, or simple
	  int vertexWeightType,   // use vertex weights?
	  int edgeWeightType,     // use edge/hyperedge weights?
	  int objectType)         // use isorropia's CrsMatrix or CrsGraph
{
  int rc=0, fail = 0;
#ifdef HAVE_EPETRAEXT
  int localProc = 0;
  double balance1, balance2, cutn1, cutn2, cutl1, cutl2;
  double balance3, cutn3, cutl3;
  double cutWgt1, cutWgt2, cutWgt3;;
  int numCuts1, numCuts2, numCuts3, valid;
  int numPartitions = 0;
  int keepDenseEdges = 0;
  int numProcs = 1;

#ifdef HAVE_MPI
  const Epetra_MpiComm &Comm = dynamic_cast<const Epetra_MpiComm &>(matrix->Comm());
  localProc = Comm.MyPID();
  numProcs = Comm.NumProc();
#else
  const Epetra_SerialComm &Comm = dynamic_cast<const Epetra_SerialComm &>(matrix->Comm());
#endif

  int numRows = matrix->NumGlobalRows();

  if (numRows < (numProcs * 100)){
    // By default Zoltan throws out dense edges, defined as those
    // whose number of non-zeros exceeds 25% of the number of vertices.
    //
    // If dense edges are thrown out of a small matrix, there may be nothing left.
    keepDenseEdges = 1;
  }

  double myShareBefore = 1.0 / numProcs;
  double myShare = myShareBefore;

  if (contract){
    numPartitions = numProcs / 2;

    if (numPartitions > numRows)
      numPartitions = numRows;

    if (numPartitions > 0){
      if (localProc < numPartitions){
	myShare = 1.0 / numPartitions;
      }
      else{
	myShare = 0.0;
      }
    }
    else{
      contract = 0;
    }
  }

  if (contract && (partitioningType == NO_ZOLTAN)){
    ERRORRETURN((localProc==0), "#Partitions < #Processes only works on Zoltan");
  }

  // If we want Zoltan's or Isorropia's default weights, then we don't
  // need to supply a CostDescriber object to create_balanced_copy,
  // so we get to test the API functions that don't take a CostDescriber.

  bool noCosts = ((vertexWeightType == NO_APPLICATION_SUPPLIED_WEIGHTS) &&
		   (edgeWeightType == NO_APPLICATION_SUPPLIED_WEIGHTS));

  // Test the interface that has no parameters, if possible

  bool noParams =
    ((partitioningType == HYPERGRAPH_PARTITIONING) && // default, so requires no params
     (numPartitions == 0) &&                          // >0 would require a parameter
     (keepDenseEdges == 0));                          // >0 would require a parameter

  // Maps for original object
  const Epetra_Map &sourceRowMap = matrix->RowMap();
  const Epetra_Map &sourceRangeMap = matrix->RangeMap();
//   const Epetra_Map &sourceColMap = matrix->ColMap();
  const Epetra_Map &sourceDomainMap = matrix->DomainMap();

  int numCols = matrix->NumGlobalCols();
  int nMyRows = sourceRowMap.NumMyElements();
  int base = sourceRowMap.IndexBase();

  // Compute vertex and edge weights

  Isorropia::Epetra::CostDescriber costs;

  Teuchos::RCP<Epetra_Vector> vptr;

  Teuchos::RCP<Epetra_CrsMatrix> eptr;

  Teuchos::RCP<Epetra_Vector> hyperEdgeWeights;

  if (edgeWeightType != NO_APPLICATION_SUPPLIED_WEIGHTS){

    if (partitioningType == GRAPH_PARTITIONING){

      // Create graph edge weights.

      eptr = Teuchos::rcp(new Epetra_CrsMatrix(*matrix));

      if (vertexWeightType == SUPPLY_EQUAL_WEIGHTS){
	eptr->PutScalar(1.0);   // set all nonzeros to 1.0
      }
      else{
	int maxRowSize = eptr->MaxNumEntries();
	double *newVal = NULL;
	if (maxRowSize > 0){
	  newVal = new double [maxRowSize];
	  for (int j=0; j<maxRowSize; j++){
	    newVal[j] = localProc + 1 + j;
	  }
	}
	int numEntries;
	int *idx;
	double *val;
	for (int i=0; i<nMyRows; i++){
	  rc = eptr->ExtractMyRowView(i, numEntries, val, idx);
	  for (int j=0; j<numEntries; j++){
	    val[j] = newVal[j];
	  }
	}
	if (newVal) delete [] newVal;
      }

      eptr->FillComplete(sourceDomainMap, sourceRangeMap);

      costs.setGraphEdgeWeights(eptr);
    }
    else{
      // Create hyperedge weights.  (Note that the list of hyperedges that a
      // process provides weights for has no relation to the columns
      // that it has non-zeroes for, or the rows that is has.  Hypergraphs
      // in general are not square.  Also more than one process can provide
      // a weight for the same edge.  Zoltan combines the weights according
      // to the value of the PHG_EDGE_WEIGHT_OPERATION parameter.  The default
      // for this parameter is to use the maximum edge weight provided by any
      // process for a given hyperedge.)

      Epetra_Map hyperEdgeMap(numCols, base, Comm);

      hyperEdgeWeights = Teuchos::rcp(new Epetra_Vector(hyperEdgeMap));

      int *edgeGIDs = NULL;
      double *weights = NULL;
      int numHEweights = hyperEdgeMap.NumMyElements();

      if (numHEweights){
	edgeGIDs = new int [numHEweights];
	weights = new double [numHEweights];

	if (edgeWeightType == SUPPLY_EQUAL_WEIGHTS){
	  for (int i=0; i<numHEweights; i++){
	    edgeGIDs[i] = hyperEdgeMap.GID(i);
	    weights[i] = 1.0;
	  }
	}
	else{
	  int hiVolumeStart = matrix->NumGlobalCols() / 3;
	  int hiVolumeEnd = hiVolumeStart * 2;
	  for (int i=0; i<numHEweights; i++){
	    edgeGIDs[i] = hyperEdgeMap.GID(i);
	    if ((edgeGIDs[i] < hiVolumeStart) || (edgeGIDs[i] >= hiVolumeEnd)){
	      weights[i] = 1.0;
	    }
	    else{
	      weights[i] = 3.0;
	    }
	  }
	}
	hyperEdgeWeights->ReplaceGlobalValues(numHEweights, weights, edgeGIDs);
      }

      if (weights){
	delete [] weights;
	delete [] edgeGIDs;
      }

      costs.setHypergraphEdgeWeights(hyperEdgeWeights);
    }
  }

  bool need_importer = false;

  if ((vertexWeightType != NO_APPLICATION_SUPPLIED_WEIGHTS) ||
      (partitioningType == NO_ZOLTAN)){

    need_importer = true;  // to redistribute row weights

    double *val = NULL;

    if (nMyRows){
      val = new double [nMyRows];

      if (vertexWeightType == SUPPLY_EQUAL_WEIGHTS){
	for (int i=0; i<nMyRows; i++){
	  val[i] = 1.0;
	}
      }
      else if (vertexWeightType == SUPPLY_UNEQUAL_WEIGHTS){
	for (int i=0; i<nMyRows; i++){
	  val[i] = 1.0 + ((localProc+1) / 2);
	}
      }
      else{
	// partitioningType is NO_ZOLTAN
	//
	// Isorropia will use row weights equal to the count of
	// non zeros in the row when doing its own row balancing.
	// We need to create that weight vector for the compute_*_metric
	// calls even though we don't give it to create_partition().

	for (int i=0; i<nMyRows; i++){
	  val[i] = matrix->NumMyEntries(i);
	}
      }
    }

    vptr = Teuchos::rcp(new Epetra_Vector(Copy, sourceRowMap, val));

    if (val) delete [] val;

    costs.setVertexWeights(vptr);
  }

  // Calculate partition quality metrics before calling Zoltan

  if (partitioningType == GRAPH_PARTITIONING){
    rc = ispatest::compute_graph_metrics(matrix->Graph(), costs,
	     myShare, balance1, numCuts1, cutWgt1, cutn1, cutl1);
    if (contract){
      // balance wrt target of balancing weight over *all* procs
      rc = ispatest::compute_graph_metrics(matrix->Graph(), costs,
	     myShareBefore, balance3, numCuts3, cutWgt3, cutn3, cutl3);
    }
  }
  else{
    rc = ispatest::compute_hypergraph_metrics(matrix->Graph(), costs,
	     myShare, balance1, cutn1, cutl1);
    if (contract){
      // balance wrt target of balancing weight over *all* procs
      rc = ispatest::compute_hypergraph_metrics(matrix->Graph(), costs,
	     myShareBefore, balance3, cutn3, cutl3);
    }
  }

  if (rc){
    ERROREXIT((localProc==0), "Error in computing partitioning metrics")
  }

  Teuchos::ParameterList params;

  if (partitioningType == NO_ZOLTAN){
    params.set("PARTITIONING_METHOD", "SIMPLE_LINEAR");
  }

#ifdef HAVE_ISORROPIA_ZOLTAN

  if ((partitioningType != NO_ZOLTAN) && !noParams){

    // We're using Zoltan for partitioning and supplying
    // parameters, overriding defaults.

    Teuchos::ParameterList &sublist = params.sublist("Zoltan");

    if (partitioningType == GRAPH_PARTITIONING){
      sublist.set("LB_METHOD", "GRAPH");
      sublist.set("GRAPH_PACKAGE", "PHG");
    }
    else{
      sublist.set("LB_METHOD", "HYPERGRAPH");
      sublist.set("LB_APPROACH", "PARTITION");
      sublist.set("PHG_CUT_OBJECTIVE", "CONNECTIVITY");  // "cutl"
    }

    if (keepDenseEdges){
      // only throw out rows that have no zeroes, default is to
      // throw out if .25 or more of the columns are non-zero
      sublist.set("PHG_EDGE_SIZE_THRESHOLD", "1.0");
    }
     if (numPartitions > 0){
	// test #Partitions < #Processes
	std::ostringstream os;
	os << numPartitions;
	std::string s = os.str();
	sublist.set("NUM_GLOBAL_PARTS", s);
      }

      //sublist.set("DEBUG_LEVEL", "1"); // Zoltan will print out parameters
      //sublist.set("DEBUG_LEVEL", "5");   // proc 0 will trace Zoltan calls
      //sublist.set("DEBUG_MEMORY", "2");  // Zoltan will trace alloc & free
  }

#else
  if (partitioningType != NO_ZOLTAN){
    ERROREXIT((localProc==0),
      "Zoltan partitioning required but Zoltan not available.")
  }
#endif

  // Function scope values

  Teuchos::RCP<Epetra_Vector> newvwgts;
  Teuchos::RCP<Epetra_CrsMatrix> newewgts;

  // Function scope values required for LinearProblem

  Epetra_LinearProblem *problem = NULL;
  Epetra_Map *LHSmap = NULL;
  Epetra_MultiVector *RHS = NULL;
  Epetra_MultiVector *LHS = NULL;

  // Reference counted pointer to balanced object

  Teuchos::RCP<Epetra_CrsMatrix> matrixPtr;
  Teuchos::RCP<Epetra_CrsGraph> graphPtr;
  Teuchos::RCP<Epetra_RowMatrix> rowMatrixPtr;
  Teuchos::RCP<Epetra_LinearProblem> problemPtr;

  // Row map for balanced object
  const Epetra_BlockMap *targetBlockRowMap=NULL;  // for input CrsGraph
  const Epetra_Map *targetRowMap=NULL;            // for all other inputs

  // Column map for balanced object
  const Epetra_BlockMap *targetBlockColMap=NULL;  // for input CrsGraph
  const Epetra_Map *targetColMap=NULL;            // for all other inputs

  if (objectType == EPETRA_CRSMATRIX){
    if (noParams && noCosts){
      matrixPtr = Isorropia::Epetra::create_balanced_copy(*matrix);
    }
    else if (noCosts){
      matrixPtr = Isorropia::Epetra::create_balanced_copy(*matrix, params);
    }
    else{
      // if noParams==true then param list is empty
      matrixPtr = Isorropia::Epetra::create_balanced_copy(*matrix, costs, params);

      if ( noParams && (vptr.get() != 0) &&
	   (eptr.get() == 0) && (hyperEdgeWeights.get() == 0)){

	// Test the interface where vector of row weights are provided,
	// but all else defaults.  This doesn't give the same answer
	// because random values are different the second time we
	// solve the same problem.

	double bal1, bal2, cn1, cn2, cl1, cl2;

	rc = ispatest::compute_hypergraph_metrics(*matrix, costs,
	     myShare, bal1, cn1, cl1);

	Teuchos::RCP<Epetra_CrsMatrix> comparePtr =
	  Isorropia::Epetra::create_balanced_copy(*matrix, *vptr);

	rc = ispatest::compute_hypergraph_metrics(*matrix, costs,
	     myShare, bal2, cn2, cl2);

	if (cl2 > cl1){
	  ERROREXIT((localProc==0),
	    "Matrix, weight vector interface failed");
	}
      }
    }
    targetRowMap = &(matrixPtr->RowMap());
    targetColMap = &(matrixPtr->ColMap());
  }
  else if (objectType == EPETRA_CRSGRAPH){
    const Epetra_CrsGraph graph = matrix->Graph();
    if (noParams && noCosts){
      graphPtr = Isorropia::Epetra::create_balanced_copy(graph);
    }
    else if (noCosts){
      graphPtr = Isorropia::Epetra::create_balanced_copy(graph, params);
    }
    else{
      // if noParams==true then param list is empty
      graphPtr = Isorropia::Epetra::create_balanced_copy(graph, costs, params);
    }
    targetBlockRowMap = &(graphPtr->RowMap());
    targetBlockColMap = &(graphPtr->ColMap());
  }
  else if (objectType == EPETRA_ROWMATRIX){
    if (noParams && noCosts){
      rowMatrixPtr = Isorropia::Epetra::create_balanced_copy(*matrix);
    }
    else if (noCosts){
      rowMatrixPtr = Isorropia::Epetra::create_balanced_copy(*matrix, params);
    }
    else{
      // if noParams==true then param list is empty
      rowMatrixPtr = Isorropia::Epetra::create_balanced_copy(*matrix, costs, params);
    }
    targetRowMap = &(rowMatrixPtr->RowMatrixRowMap());
    targetColMap = &(rowMatrixPtr->RowMatrixColMap());
  }
  else if (objectType == EPETRA_LINEARPROBLEM){

    // Create a linear problem with this matrix.

    LHSmap = new Epetra_Map(numCols, base, Comm);

    int myRHSsize = sourceRowMap.NumMyElements();
    int myLHSsize = LHSmap->NumMyElements();

    int valSize = ((myRHSsize > myLHSsize) ? myRHSsize : myLHSsize);

    double *vals = NULL;

    if (valSize){
      vals = new double [valSize];
    }

    if (valSize){
      for (int i=0; i < valSize; i++){
	// put my rank in my portion of LHS and my portion of RHS
	vals[i] = localProc;
      }
    }

    RHS = new Epetra_MultiVector(Copy, sourceRowMap, vals, 1, 1);

    LHS = new Epetra_MultiVector(Copy, *LHSmap, vals, 1, 1);

    if (valSize){
      delete [] vals;
    }

    problem = new Epetra_LinearProblem(matrix.get(), LHS, RHS);

    Epetra_LinearProblem lp = *problem;

    if (lp.CheckInput()){
      ERROREXIT((localProc==0), "Error creating a LinearProblem");
    }
    if (noParams && noCosts){
      problemPtr = Isorropia::Epetra::create_balanced_copy(lp);
    }
    else if (noCosts){
      problemPtr = Isorropia::Epetra::create_balanced_copy(lp, params);
    }
    else{
      // if noParams==true then param list is empty
      problemPtr = Isorropia::Epetra::create_balanced_copy(lp, costs, params);
    }

    targetRowMap = &(problemPtr->GetMatrix()->RowMatrixRowMap());
    targetColMap = &(problemPtr->GetMatrix()->RowMatrixColMap());
  }

  // Redistribute the edge weights
  // Comment this out since we don't redistribute columns

  if (edgeWeightType != NO_APPLICATION_SUPPLIED_WEIGHTS){

    if (partitioningType == GRAPH_PARTITIONING){

      Epetra_Import *importer = NULL;

      if (objectType == EPETRA_CRSGRAPH){
	newewgts = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *graphPtr));
	targetRowMap = &(newewgts->RowMap());
	targetColMap = &(newewgts->ColMap());
      }
      else{
	newewgts = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *targetRowMap, *targetColMap, 0));
      }

      importer = new Epetra_Import(*targetRowMap, sourceRowMap);
      newewgts->Import(*eptr, *importer, Insert);
      newewgts->FillComplete(*targetColMap, *targetRowMap);

      costs.setGraphEdgeWeights(newewgts);
    }
  }

  // Redistribute the vertex weights

  if ((vertexWeightType != NO_APPLICATION_SUPPLIED_WEIGHTS) ||
      (partitioningType == NO_ZOLTAN)){

    Epetra_Import *importer = NULL;

    if (objectType == EPETRA_CRSGRAPH){
      newvwgts = Teuchos::rcp(new Epetra_Vector(*targetBlockRowMap));
      importer = new Epetra_Import(*targetBlockRowMap, sourceRowMap);
    }
    else{
      newvwgts = Teuchos::rcp(new Epetra_Vector(*targetRowMap));
      importer = new Epetra_Import(*targetRowMap, sourceRowMap);
    }

    newvwgts->Import(*vptr, *importer, Insert);
    costs.setVertexWeights(newvwgts);
  }

  if (localProc == 0){
    test_type(numPartitions, partitioningType, vertexWeightType, edgeWeightType, objectType);
  }

  if (verbose){

    // Picture of problem before balancing

    if (objectType == EPETRA_LINEARPROBLEM){

      ispatest::show_matrix("Before load balancing", *problem, Comm);
    }
    else{
      ispatest::show_matrix("Before load balancing", matrix->Graph(), Comm);
    }

    // Picture of problem after balancing

    if (objectType == EPETRA_LINEARPROBLEM){
      ispatest::show_matrix("After load balancing (x in Ax=b is not redistributed)", *problemPtr, Comm);
    }
    else if (objectType == EPETRA_ROWMATRIX){
      ispatest::show_matrix("After load balancing", *rowMatrixPtr, Comm);
    }
    else if (objectType == EPETRA_CRSMATRIX){
      ispatest::show_matrix("After load balancing", matrixPtr->Graph(), Comm);
    }
    else if (objectType == EPETRA_CRSGRAPH){
      ispatest::show_matrix("After load balancing", *graphPtr, Comm);
    }
  }

  // After partitioning, recompute the metrics

  if (partitioningType == GRAPH_PARTITIONING){
    if (objectType == EPETRA_LINEARPROBLEM){
      rc = ispatest::compute_graph_metrics(*(problemPtr->GetMatrix()), costs,
	     myShare, balance2, numCuts2, cutWgt2, cutn2, cutl2);
    }
    else if (objectType == EPETRA_ROWMATRIX){
      rc = ispatest::compute_graph_metrics(*rowMatrixPtr, costs,
	     myShare, balance2, numCuts2, cutWgt2, cutn2, cutl2);
    }
    else if (objectType == EPETRA_CRSMATRIX){
      rc = ispatest::compute_graph_metrics(matrixPtr->Graph(), costs,
	     myShare, balance2, numCuts2, cutWgt2, cutn2, cutl2);
    }
    else {
      rc = ispatest::compute_graph_metrics(*graphPtr, costs,
	     myShare, balance2, numCuts2, cutWgt2, cutn2, cutl2);
    }
  }
  else{
    if (objectType == EPETRA_LINEARPROBLEM){
      rc = ispatest::compute_hypergraph_metrics(*(problemPtr->GetMatrix()), costs,
	     myShare, balance2, cutn2, cutl2);
    }
    else if (objectType == EPETRA_ROWMATRIX){
      rc = ispatest::compute_hypergraph_metrics(*rowMatrixPtr, costs,
	     myShare, balance2, cutn2, cutl2);
    }
    else if (objectType == EPETRA_CRSMATRIX){
      rc = ispatest::compute_hypergraph_metrics(matrixPtr->Graph(), costs,
	     myShare, balance2, cutn2, cutl2);
    }
    else{
      rc = ispatest::compute_hypergraph_metrics(*graphPtr, costs,
	     myShare, balance2, cutn2, cutl2);
    }
  }

  if (rc){
    ERROREXIT((localProc==0), "Error in computing partitioning metrics")
  }

  std::string why;

  if (partitioningType == GRAPH_PARTITIONING){
    fail = (cutWgt2 > cutWgt1);
    why = "New weighted edge cuts are worse";

    if (localProc == 0){
      std::cout << "Before partitioning: Balance " << balance1 ;
      std::cout << " cutn " << cutn1 ;
      std::cout << " cutl " << cutl1 ;

      if (contract){
	std::cout << "  (wrt balancing over " << numPartitions << " partitions)" << std::endl;
	std::cout << "Before partitioning: Balance " << balance3 ;
	std::cout << " cutn " << cutn3 ;
	std::cout << " cutl " << cutl3 ;
	std::cout << "  (wrt balancing over " << numProcs << " partitions)" ;
      }
      std::cout << std::endl;

      std::cout << " Total edge cuts: " << numCuts1;
      std::cout << " Total weighted edge cuts: " << cutWgt1 << std::endl;
      std::cout << "After partitioning: Balance " << balance2 ;
      std::cout << " cutn " << cutn2 ;
      std::cout << " cutl " << cutl2 << std::endl;
      std::cout << " Total edge cuts: " << numCuts2;
      std::cout << " Total weighted edge cuts: " << cutWgt2 << std::endl;
    }
  }
  else{
    if (partitioningType == NO_ZOLTAN){
      fail = (balance2 > balance1);
      why = "New balance is worse";
    }
    else{
      fail = (cutl2 > cutl1);
      why = "New cutl is worse";
    }

    if (localProc == 0){
      std::cout << "Before partitioning: Balance " << balance1 ;
      std::cout << " cutn " << cutn1 ;
      std::cout << " cutl " << cutl1 ;
      if (contract){
	std::cout << "  (wrt balancing over " << numPartitions << " partitions)" << std::endl;
	std::cout << "Before partitioning: Balance " << balance3 ;
	std::cout << " cutn " << cutn3 ;
	std::cout << " cutl " << cutl3 ;
	std::cout << "  (wrt balancing over " << numProcs << " partitions)" ;
      }
      std::cout << std::endl;
      std::cout << "After partitioning: Balance " << balance2 ;
      std::cout << " cutn " << cutn2 ;
      std::cout << " cutl " << cutl2 << std::endl;
    }
  }

  if (fail){
    if (localProc == 0) std::cout << "ERROR: "+why << std::endl;
  }

  // Check that input matrix is valid.  This test constructs an "x"
  // with the matrix->DomainMap() and a "y" with matrix->RangeMap()
  // and then calculates y = Ax.

  if (objectType == EPETRA_LINEARPROBLEM){
    valid = ispatest::test_matrix_vector_multiply(*problemPtr);
  }
  else if (objectType == EPETRA_ROWMATRIX){
    valid = ispatest::test_row_matrix_vector_multiply(*rowMatrixPtr);
  }
  else if (objectType == EPETRA_CRSMATRIX){
    valid = ispatest::test_matrix_vector_multiply(*matrixPtr);
  }
  else{
    valid = ispatest::test_matrix_vector_multiply(*graphPtr);
  }

  if (!valid){
    if (localProc == 0) std::cout << "Rebalanced matrix is not a valid Epetra matrix" << std::endl;
    fail = 1;
  }
  else{
    if (localProc == 0) std::cout << "Rebalanced matrix is a valid Epetra matrix" << std::endl;
  }

  if (localProc == 0)
    std::cout << std::endl;



#else
  std::cout << "test_simple main: currently can only test "
	 << "with Epetra and EpetraExt enabled." << std::endl;
  rc = -1;
#endif

  return fail;
}

int main(int argc, char** argv) {

  int rc=0, fail = 0;
#ifdef HAVE_EPETRAEXT
  bool verbose = false;
  int localProc = 0;
//   std::string *fstr;

#ifdef HAVE_MPI
  int numProcs;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &localProc);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
  const Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  const Epetra_SerialComm Comm;
#endif

   //if (getenv("DEBUGME")){
   // std::cout << localProc << " gdb test_create_balanced_copy.exe " << getpid() << std::endl;
   // sleep(15);
   //}

  Teuchos::CommandLineProcessor clp(false,true);

  // --f=fileName provides a different matrix market file for input
  // --v will print out the partitioning (small files only)

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
    if (localProc==0){
      std::cout << "error reading input file" << std::endl << "FAIL" << std::endl;
    }
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

#ifdef SHORT_TEST
  fail = run_test(testm,
	     verbose,
	     false,                 // do not test #partitions < #processes
	     NO_ZOLTAN,
	     NO_APPLICATION_SUPPLIED_WEIGHTS,
	     NO_APPLICATION_SUPPLIED_WEIGHTS,
	     EPETRA_CRSGRAPH);

  CHECK_FAILED();
  goto Report;

#else

  if (square){
#ifdef HAVE_ISORROPIA_ZOLTAN
    fail = run_test(testm,            // test matrix
	       verbose,               // display matrix before and after?
	       false,                 // do not test #partitions < #processes
	       GRAPH_PARTITIONING,    // perform zoltan graph partitioning
	       SUPPLY_EQUAL_WEIGHTS,  // supply equal vertex weights
	       SUPPLY_EQUAL_WEIGHTS,  // supply equal edge weights
	       EPETRA_LINEARPROBLEM); // use linear problem interface of isorropia

    CHECK_FAILED();

    fail = run_test(testm,
	       verbose,            // draw graph before and after partitioning?
	       false,                 // do not test #partitions < #processes
	       HYPERGRAPH_PARTITIONING,      // do graph partitioning
	       SUPPLY_EQUAL_WEIGHTS,    // supply vertex weights, all the same
	       NO_APPLICATION_SUPPLIED_WEIGHTS,  // go for default weights
	       EPETRA_CRSMATRIX);       // use the Epetra_CrsMatrix interface

    CHECK_FAILED();

    fail = run_test(testm,
	       verbose,
	       true,                 // test #partitions < #processes
	       GRAPH_PARTITIONING,
	       SUPPLY_UNEQUAL_WEIGHTS,
	       SUPPLY_EQUAL_WEIGHTS,
	       EPETRA_CRSMATRIX);

    CHECK_FAILED();

    fail = run_test(testm,
	       verbose,
	       false,                 // do not test #partitions < #processes
	       GRAPH_PARTITIONING,
	       SUPPLY_EQUAL_WEIGHTS,
	       SUPPLY_EQUAL_WEIGHTS,
	       EPETRA_LINEARPROBLEM);

    CHECK_FAILED();

    fail = run_test(testm,
	       verbose,
	       false,                 // do not test #partitions < #processes
	       GRAPH_PARTITIONING,
	       NO_APPLICATION_SUPPLIED_WEIGHTS,
	       NO_APPLICATION_SUPPLIED_WEIGHTS,
	       EPETRA_ROWMATRIX);

    CHECK_FAILED();
#else
  fail = 0;
  if (localProc == 0){
    std::cout << "Test not run because it requires EPETRA_EXT" << std::endl;
  }
#endif

    fail = run_test(testm,
	       verbose,
	       false,                 // do not test #partitions < #processes
	       NO_ZOLTAN,
	       SUPPLY_UNEQUAL_WEIGHTS,
	       SUPPLY_EQUAL_WEIGHTS,
	       EPETRA_CRSMATRIX);

    CHECK_FAILED();

    fail = run_test(testm,
	       verbose,
	       false,                 // do not test #partitions < #processes
	       NO_ZOLTAN,
	       SUPPLY_EQUAL_WEIGHTS,
	       SUPPLY_EQUAL_WEIGHTS,
	       EPETRA_LINEARPROBLEM);

    CHECK_FAILED();
  }

  fail = run_test(testm,
	     verbose,
	     false,                 // do not test #partitions < #processes
	     NO_ZOLTAN,
	     SUPPLY_EQUAL_WEIGHTS,
	     SUPPLY_EQUAL_WEIGHTS,
	     EPETRA_CRSMATRIX);

   CHECK_FAILED();

  fail = run_test(testm,
	     verbose,
	     false,                 // do not test #partitions < #processes
	     NO_ZOLTAN,
	     NO_APPLICATION_SUPPLIED_WEIGHTS,
	     NO_APPLICATION_SUPPLIED_WEIGHTS,
	     EPETRA_CRSGRAPH);

   CHECK_FAILED();

  fail = run_test(testm,
	     verbose,
	     false,                 // do not test #partitions < #processes
	     NO_ZOLTAN,
	     SUPPLY_UNEQUAL_WEIGHTS,
	     NO_APPLICATION_SUPPLIED_WEIGHTS,
	     EPETRA_CRSMATRIX);

   CHECK_FAILED();


#ifdef HAVE_ISORROPIA_ZOLTAN

  fail = run_test(testm,
	     verbose,
	     true,                 // test #partitions < #processes
	     HYPERGRAPH_PARTITIONING,
	     SUPPLY_EQUAL_WEIGHTS,
	     SUPPLY_UNEQUAL_WEIGHTS,
	     EPETRA_CRSGRAPH);

   CHECK_FAILED();

  fail = run_test(testm,
	     verbose,
	     false,                 // do not test #partitions < #processes
	     HYPERGRAPH_PARTITIONING,
	     SUPPLY_UNEQUAL_WEIGHTS,
	     SUPPLY_EQUAL_WEIGHTS,
	     EPETRA_ROWMATRIX);

   CHECK_FAILED();

  fail = run_test(testm,
	     verbose,
	     false,                 // do not test #partitions < #processes
	     HYPERGRAPH_PARTITIONING,
	     NO_APPLICATION_SUPPLIED_WEIGHTS,
	     NO_APPLICATION_SUPPLIED_WEIGHTS,
	     EPETRA_LINEARPROBLEM);

   CHECK_FAILED();

#endif

#endif // SHORT_TEST

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
