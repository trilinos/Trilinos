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
#include <Isorropia_EpetraPartitioner.hpp>
#include <Isorropia_EpetraRedistributor.hpp>
#include <Isorropia_EpetraCostDescriber.hpp>

#include <ispatest_lbeval_utils.hpp>

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
#endif

#endif

#include <Teuchos_CommandLineProcessor.hpp>

#define GRAPH_PARTITIONING            1
#define HYPERGRAPH_PARTITIONING       2
#define NO_ZOLTAN                     3

#define SUPPLY_EQUAL_WEIGHTS               1
#define SUPPLY_UNEQUAL_WEIGHTS             2
#define NO_APPLICATION_SUPPLIED_WEIGHTS    3

#define EPETRA_CRSGRAPH               1
#define EPETRA_CRSMATRIX              2  

#define ERROREXIT(v, s) \
  if (v) std::cout << s << std::endl << "FAIL" << std::endl; \
  exit(1);

static void test_type(int partitioningType, int vertexWeightType, 
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
  else
    std::cout << "using Epetra_CrsMatrix interface";

  std::cout << std::endl;
}

static int run_test(Teuchos::RCP<Epetra_CrsMatrix> matrix,
          bool verbose,           // display the graph before & after
          int partitioningType,   // hypergraph or graph partitioning, or simple
          int vertexWeightType,   // use vertex weights?
          int edgeWeightType,     // use edge/hyperedge weights?
          int objectType)         // use isorropia's CrsMatrix or CrsGraph
{
  int rc=0, fail = 0;  
#ifdef HAVE_EPETRAEXT
  int numProcs = 1;
  int localProc = 0;
  double balance1, balance2, cutn1, cutn2, cutl1, cutl2;
  double cutWgt1, cutWgt2; 
  int numCuts1, numCuts2;

#ifdef HAVE_MPI
  const Epetra_MpiComm &Comm = dynamic_cast<const Epetra_MpiComm &>(matrix->Comm());
  localProc = Comm.MyPID();
  numProcs = Comm.NumProc();
#else
  const Epetra_SerialComm &Comm = dynamic_cast<const Epetra_SerialComm &>(matrix->Comm());
#endif

  if (verbose)
    ispatest::show_matrix("Before load balancing", matrix->Graph(), Comm);


  // Compute vertex and edge weights

  Teuchos::RCP<Isorropia::Epetra::CostDescriber> costs =
    Teuchos::rcp(new Isorropia::Epetra::CostDescriber);

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
        int nrows = eptr->NumMyRows();
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
        for (int i=0; i<nrows; i++){
          rc = eptr->ExtractMyRowView(i, numEntries, val, idx);
          for (int j=0; j<numEntries; j++){
            val[j] = newVal[j];
          }
        }
        if (newVal) delete [] newVal;
      }
    
      eptr->FillComplete();
      eptr->OptimizeStorage();
    
      costs->setGraphEdgeWeights(eptr);
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
    
      Epetra_Map hyperEdgeMap(matrix->NumGlobalCols(), matrix->IndexBase(), Comm);
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

      costs->setHypergraphEdgeWeights(hyperEdgeWeights);
    }
  }
  
  if ((vertexWeightType != NO_APPLICATION_SUPPLIED_WEIGHTS) ||
      (partitioningType == NO_ZOLTAN)){
  
    const Epetra_Map &rowmap = matrix->RowMap();
    int nrows = rowmap.NumMyElements();
    double *val = NULL;
  
    if (nrows){
      val = new double [nrows];

      if (vertexWeightType == SUPPLY_EQUAL_WEIGHTS){
        for (int i=0; i<nrows; i++){
          val[i] = 1.0;
        }
      }
      else if (vertexWeightType == SUPPLY_UNEQUAL_WEIGHTS){
        for (int i=0; i<nrows; i++){
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

        for (int i=0; i<nrows; i++){
          val[i] = matrix->NumMyEntries(i);
        }
      }
    }
  
    vptr = Teuchos::rcp(new Epetra_Vector(Copy, rowmap, val));
  
    if (val) delete [] val;
  
    costs->setVertexWeights(vptr);
  }

  // Calculate partition quality metrics before calling Zoltan

  if (partitioningType == GRAPH_PARTITIONING){
    rc = ispatest::compute_graph_metrics(matrix->Graph(), *costs, 
             balance1, numCuts1, cutWgt1, cutn1, cutl1);
  }
  else{
    rc = ispatest::compute_hypergraph_metrics(matrix->Graph(), *costs,
             balance1, cutn1, cutl1);
  }

  if (rc){ 
    ERROREXIT((localProc==0), "Error in computing partitioning metrics")
  }

  Teuchos::ParameterList params;

  if (partitioningType == NO_ZOLTAN){
    params.set("PARTITIONING_METHOD", "SIMPLE_LINEAR");
  }

#ifdef HAVE_ISORROPIA_ZOLTAN

  // Set the Zoltan parameters for this problem

  if (partitioningType != NO_ZOLTAN){
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
  }
#else
  if (partitioningType != NO_ZOLTAN){
    ERROREXIT((localProc==0), 
      "Zoltan partitioning required but Zoltan not available.")
  }
#endif

  // Perform hyperedge partitioning with Zoltan (if we have it)


  Teuchos::RCP<Isorropia::Epetra::Partitioner> partitioner;

  if ((edgeWeightType == NO_APPLICATION_SUPPLIED_WEIGHTS) &&
      (vertexWeightType == NO_APPLICATION_SUPPLIED_WEIGHTS)){

    if (objectType == EPETRA_CRSGRAPH){
      // Test the Epetra_CrsGraph interface of Isorropia
      Teuchos::RCP<const Epetra_CrsGraph> graph = 
        Teuchos::rcp(new Epetra_CrsGraph(matrix->Graph()));
      partitioner = Isorropia::Epetra::create_partitioner(graph, params);
    }
    else{
      // Test the Epetra_CrsMatrix interface of Isorropia
      Teuchos::RCP<const Epetra_RowMatrix> rm = matrix;
      partitioner = Isorropia::Epetra::create_partitioner(rm, params);
    }

  }
  else{

    if (objectType == EPETRA_CRSGRAPH){
      // Test the Epetra_CrsGraph interface of Isorropia
      Teuchos::RCP<const Epetra_CrsGraph> graph = 
        Teuchos::rcp(new Epetra_CrsGraph(matrix->Graph()));
      partitioner = Isorropia::Epetra::create_partitioner(graph, costs, params);
    }
    else{
      // Test the Epetra_CrsMatrix interface of Isorropia
      Teuchos::RCP<const Epetra_RowMatrix> rm = matrix;
      partitioner = Isorropia::Epetra::create_partitioner(rm, costs, params);
    }
  }

  // Create a Redistributor based on the partitioning

  Isorropia::Epetra::Redistributor rd(partitioner);

  // Redistribute the matrix

  Teuchos::RCP<Epetra_CrsMatrix> newMatrix = rd.redistribute(*matrix);

  // Redistribute the vertex weights

  if ((vertexWeightType != NO_APPLICATION_SUPPLIED_WEIGHTS) ||
      (partitioningType == NO_ZOLTAN)){

    Teuchos::RCP<Epetra_Vector> newvwgts = rd.redistribute(*vptr);
    costs->setVertexWeights(newvwgts);
  }

  // Redistribute the edge weights

  if (edgeWeightType != NO_APPLICATION_SUPPLIED_WEIGHTS){

    if (partitioningType == GRAPH_PARTITIONING){
      Teuchos::RCP<Epetra_CrsMatrix> newewgts = rd.redistribute(*eptr);
      costs->setGraphEdgeWeights(newewgts);
    }
  }

  if (verbose)
    ispatest::show_matrix("After load balancing", newMatrix.get()->Graph(), Comm);

  // After partitioning, recompute the metrics

  if (partitioningType == GRAPH_PARTITIONING){
    rc = ispatest::compute_graph_metrics(newMatrix->Graph(), *costs, 
             balance2, numCuts2, cutWgt2, cutn2, cutl2);
  }
  else{
    rc = ispatest::compute_hypergraph_metrics(newMatrix->Graph(), *costs,
             balance2, cutn2, cutl2);
  }

  if (rc){ 
    ERROREXIT((localProc==0), "Error in computing partitioning metrics")
  }

  std::cout << std::endl;
  if (localProc == 0){
    test_type(partitioningType, vertexWeightType, edgeWeightType, objectType);
  }
  
  if (partitioningType == GRAPH_PARTITIONING){
    fail = (cutWgt2 > cutWgt1);
  
    if (localProc == 0){
      std::cout << "Before partitioning: Balance " << balance1 ;
      std::cout << " cutn " << cutn1 ;
      std::cout << " cutl " << cutl1 << std::endl;
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
    }
    else{
      fail = (cutl2 > cutl1);         // Zoltan hypergraph partitioning
    }
  
    if (localProc == 0){
      std::cout << "Before partitioning: Balance " << balance1 ;
      std::cout << " cutn " << cutn1 ;
      std::cout << " cutl " << cutl1 << std::endl;
      std::cout << "After partitioning: Balance " << balance2 ;
      std::cout << " cutn " << cutn2 ;
      std::cout << " cutl " << cutl2 << std::endl;
    }
  }


#else
  std::cout << "test_simple : currently can only test "
         << "with Epetra and EpetraExt enabled." << std::endl;
  rc = -1;
#endif

  return fail;
}

int main(int argc, char** argv) {

  int rc=0, fail = 0;  
#ifdef HAVE_EPETRAEXT
  bool verbose = false;
  int numProcs = 1;
  int localProc = 0;
  std::string *fstr;

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

  std::string *inputFile = new std::string("simple.mtx");

  clp.setOption( "f", inputFile, 
                "Name of input matrix market file");
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

  Epetra_CrsMatrix *matrixPtr;
  rc = EpetraExt::MatrixMarketFileToCrsMatrix(fname, Comm, matrixPtr);
  if (rc < 0){
    ERROREXIT((localProc==0), "error reading input file");
  }

  bool square = (matrixPtr->NumGlobalRows() == matrixPtr->NumGlobalCols());

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
               GRAPH_PARTITIONING,      // do graph partitioning
               SUPPLY_EQUAL_WEIGHTS,    // supply vertex weights, all the same
               SUPPLY_EQUAL_WEIGHTS,    // supply edge weights, all the same
               EPETRA_CRSMATRIX);       // use the Epetra_CrsMatrix interface

    if (fail){
      goto Report;
    }

    fail = run_test(testm,
               verbose,
               GRAPH_PARTITIONING, 
               SUPPLY_UNEQUAL_WEIGHTS,
               SUPPLY_EQUAL_WEIGHTS,
               EPETRA_CRSMATRIX);
  
    if (fail){
      goto Report;
      return 1;   
    }
  
    fail = run_test(testm,
               verbose,
               GRAPH_PARTITIONING,
               SUPPLY_EQUAL_WEIGHTS,
               SUPPLY_EQUAL_WEIGHTS,
               EPETRA_CRSGRAPH);
  
    if (fail){
      goto Report;
      return 1;   
    }

    fail = run_test(testm,
               verbose,
               GRAPH_PARTITIONING,
               NO_APPLICATION_SUPPLIED_WEIGHTS,
               NO_APPLICATION_SUPPLIED_WEIGHTS,
               EPETRA_CRSGRAPH);

    if (fail){
      goto Report;
    }
#endif
  }

#ifdef HAVE_ISORROPIA_ZOLTAN

  fail = run_test(testm,
             verbose, 
             HYPERGRAPH_PARTITIONING,
             SUPPLY_EQUAL_WEIGHTS,
             SUPPLY_EQUAL_WEIGHTS,
             EPETRA_CRSMATRIX);

  if (fail){
    goto Report;
  }


  fail = run_test(testm,
             verbose,
             HYPERGRAPH_PARTITIONING,
             SUPPLY_EQUAL_WEIGHTS,
             SUPPLY_UNEQUAL_WEIGHTS,
             EPETRA_CRSGRAPH);

  if (fail){
    goto Report;
  }

  fail = run_test(testm,
             verbose, 
             HYPERGRAPH_PARTITIONING,
             SUPPLY_UNEQUAL_WEIGHTS,
             SUPPLY_EQUAL_WEIGHTS,
             EPETRA_CRSMATRIX);

  if (fail){
    goto Report;
  }

  fail = run_test(testm,
             verbose,
             HYPERGRAPH_PARTITIONING,
             NO_APPLICATION_SUPPLIED_WEIGHTS,
             NO_APPLICATION_SUPPLIED_WEIGHTS,
             EPETRA_CRSGRAPH);

  if (fail){
    goto Report;
  }
#endif

  // Default row weight is number of non zeros in the row
  fail = run_test(testm,
             verbose,
             NO_ZOLTAN,
             NO_APPLICATION_SUPPLIED_WEIGHTS,
             NO_APPLICATION_SUPPLIED_WEIGHTS,
             EPETRA_CRSGRAPH);

  if (fail){
    goto Report;
  }

  fail = run_test(testm,
             verbose,
             NO_ZOLTAN,
             SUPPLY_EQUAL_WEIGHTS,
             NO_APPLICATION_SUPPLIED_WEIGHTS,
             EPETRA_CRSGRAPH);

  if (fail){
    goto Report;
  }

  fail = run_test(testm,
             verbose,
             NO_ZOLTAN,
             SUPPLY_UNEQUAL_WEIGHTS,
             NO_APPLICATION_SUPPLIED_WEIGHTS,
             EPETRA_CRSGRAPH);

  if (fail){
    goto Report;
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
    if (fail)
      std::cout << "FAIL" << std::endl;
    else
      std::cout << "PASS" << std::endl;
  }

  return fail;
}
