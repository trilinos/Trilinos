#include <iostream>

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>

/**********************************************************************************/
/* CREATE INITAL MATRIX                                                           */
/**********************************************************************************/
#define CTHULHU_USE_EPETRA
#include <Cthulhu_Map.hpp>
#include <Cthulhu_CrsMatrix.hpp>
#include <Cthulhu_EpetraCrsMatrix.hpp>
#include <Cthulhu_CrsOperator.hpp>
#include <Cthulhu.hpp>

#define CTHULHU_ENABLED //TODO
#include <MueLu_MatrixFactory.hpp>
/**********************************************************************************/

#include <MueLu_AggAlgorithm.hpp>

// For the moment, this file is just a modified version of ML_Linker.hpp

/**********************************************************************************/
/* Function to build MueLu_Graph (take an Epetra_CrsMatrix and                    */
/* extract out the Epetra_CrsGraph).                                              */
/**********************************************************************************/
MueLu_Graph *MueLu_BuildGraph(const Epetra_CrsMatrix & A, const char *name="")
{
  MueLu_Graph *graph;
  double *dtmp = NULL;

  graph = (MueLu_Graph *) malloc(sizeof(MueLu_Graph));
  graph->eGraph = NULL;
  graph->name = NULL;
  graph->name = (char *) malloc(sizeof(char)*80); strcpy(graph->name,name);
  graph->nVertices = A.NumMyRows(); // Q:local or global ? Amatrix->invec_leng;

  if ( A.NumMyRows() == 0) { // Q: ML_Operator* Amatrix->getrow->Nrows is local or global ?
     graph->vertexNeighbors    = NULL;
     graph->vertexNeighborsPtr = NULL;
     graph->nEdges             = 0;
  }
  else {
    A.ExtractCrsDataPointers(graph->vertexNeighborsPtr, graph->vertexNeighbors, dtmp);

    graph->nEdges = (graph->vertexNeighborsPtr)[A.NumMyRows()]; // Q:
    graph->eGraph = &(A.Graph());
  }
  if (graph->eGraph == NULL) graph->nGhost = 0;
  else {
    graph->nGhost = A.RowMatrixColMap().NumMyElements() - A.OperatorDomainMap().NumMyElements();
    if (graph->nGhost < 0) graph->nGhost = 0;
  }
  return graph;
}

int MueLu_DestroyGraph(MueLu_Graph *graph)
{
   if ( graph != NULL) {
      if (graph->name != NULL) free(graph->name);
      free(graph);
   }
   return 0;
}

int main(int argc, char *argv[]) {
  
  std::cout << "Hello World !" << std::endl;

  /**********************************************************************************/
  /* CREATE INITAL MATRIX                                                           */
  /**********************************************************************************/
  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);
  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  LO numThreads=1;
  GO nx=4;
  GO ny=4;
  GO nz=4;
  Teuchos::CommandLineProcessor cmdp(false,true);
  std::string matrixType("Laplace1D");
  cmdp.setOption("nt",&numThreads,"number of threads.");
  cmdp.setOption("nx",&nx,"mesh points in x-direction.");
  cmdp.setOption("ny",&ny,"mesh points in y-direction.");
  cmdp.setOption("nz",&nz,"mesh points in z-direction.");
  cmdp.setOption("matrixType",&matrixType,"matrix type: Laplace1D, Laplace2D, Star2D, Laplace3D");
  if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return EXIT_FAILURE;
  }

  std::cout << "#threads = " << numThreads << std::endl;
  std::cout << "problem size = " << nx*ny << std::endl;
  std::cout << "matrix type = " << matrixType << std::endl;

  Teuchos::ParameterList pl;
  pl.set("Num Threads",numThreads);

  GO numGlobalElements = nx*ny;
  if (matrixType == "Laplace3D")
    numGlobalElements *= nz;
  LO indexBase = 0;

  RCP<const Map > map;
  map = rcp( new MyMap(numGlobalElements, indexBase, comm) );

  Teuchos::ParameterList matrixList;
  matrixList.set("nx",nx);
  matrixList.set("ny",ny);
  matrixList.set("nz",nz);

  RCP<CrsOperator> Op = MueLu::Gallery::CreateCrsMatrix<SC,LO,GO, Map, CrsOperator>(matrixType,map,matrixList); //TODO: Operator vs. CrsOperator

  RCP<const Epetra_CrsMatrix> A;

  { // Get the underlying Epetra Mtx (Wow ! It's paintful ! => I should create a function to do that)
    RCP<const CrsMatrix> tmp_CrsMtx = Op->get_CrsMatrix();
    const RCP<const Cthulhu::EpetraCrsMatrix> &tmp_ECrsMtx = Teuchos::rcp_dynamic_cast<const Cthulhu::EpetraCrsMatrix>(tmp_CrsMtx);
    if (tmp_ECrsMtx == Teuchos::null) { std::cout << "Error !" << std::endl; return 1; }
    A = tmp_ECrsMtx->getEpetra_CrsMatrix();
  }
  
  MueLu_Graph *graph;
  std::string name = "Uncoupled";
  graph = MueLu_BuildGraph(*A, name.c_str());
  
  int printFlag=6;
  if (graph->eGraph->Comm().MyPID() == 0 && printFlag < MueLu_PrintLevel())
    printf("main() Aggregate_CoarsenUncoupled : \n");
  
  MueLu_AggOptions AggregateOptions;
  
  AggregateOptions.printFlag               = printFlag;      
  AggregateOptions.minNodesPerAggregate    = 2;  
  AggregateOptions.maxNeighAlreadySelected = 5;
  AggregateOptions.ordering                = 1;
  AggregateOptions.phase3AggCreation       = 0.5;
  
  Aggregates *aggregates = NULL;
  
  aggregates = MueLu_Aggregate_CoarsenUncoupled(&AggregateOptions,graph);

  name = "UC_CleanUp";
  MueLu_AggregateLeftOvers(&AggregateOptions, aggregates, name.c_str(), graph);
  
  Epetra_IntVector Final( aggregates->vertex2AggId->Map() );
  for (int i = 0; i < aggregates->vertex2AggId->Map().NumMyElements(); i++) 
    Final[i] = (*(aggregates->vertex2AggId))[i] + (*(aggregates->procWinner))[i]*1000;
  printf("finals\n");
  cout << Final << endl; sleep(2);
  
  MueLu_AggregateDestroy(aggregates); 
  MueLu_DestroyGraph(graph);
  
  return 0;

}
