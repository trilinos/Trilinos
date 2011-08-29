// @HEADER
// ***********************************************************************
//
//         Zoltan2: Sandia Partitioning Ordering & Coloring Library
//
//                Copyright message goes here.   TODO
//
// Questions? Contact Lee Ann Riesen (lriesen@sandia.gov)
//
// ***********************************************************************

// Testing graph input adapters.
//
//  Find a few small mtx files for graphs.
//  for each graph:
//     Read it in with EpetraExt
//     Create an Epetra::CrsGraph
//     Instaniate a Zoltan2::GraphInputAdapter
//     Test all queries, and verify.
//     Test copy and assignment.
//     Test reset - assigning a new graph to adapter.
//     Ditto for Tpetra
//
#include <EpetraExt_CrsMatrixIn.h>
#include <EpetraCrsGraph.h>
#include <TeuchosRCP.h>

static int numGraphs = 7;
static char **mtxFiles={
"ewgt.mtx", 
"nograph.mtx", 
"simple.mtx", 
"cage10.mtx", 
"diag500_4.mtx", 
"simple3d.mtx", 
"vwgt.mtx"
};

// Read an EpetraCrsMatrix from a file.  Diagonal weights are vertex weight.
//  off diagonal weights are edge weights.

template <typename Scalar>
Teuchos::RCP<Epetra::CrsMatrix> &getEpetraCrsGraphFromFile(std::string fname)
{
  Epetra::CrsMatrix *M;
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  int rc = EpetraExt::MatrixMarketFileToCrsMatrix(fname.c_str(), comm, M);

  if (rc < 0){
    std::cerr << "can't read file " << fname.c_str() << std::endl;
    exit(1);
  }

  Teuchos::RCP<Epetra::CrsMatrix> A = Teuchos::rcp(M, false);
  return A;
}


int main(int argc, char *argv[])
{

  std::string outcome("PASS");

  // We're testing the Graph input adapters.  We'll use the graph from a CrsMatrix
  // as our graph, and the non-zeros of the matrix as our weights.

  for (int i = 0; i < numGraph; is++){
    Teuchos::RCP<EpetraCrsMatrix> &mtx= getEpetraCrsGraphFromFile(mtxFiles[i]);
    int numMyVertices = mtx->NumMyRows();
    int numMyEdges = mtx->NumMyNonZeros();
    int maxNumEntries = mtx->MaxNumEntries();

    double *values = new double [maxNumEntries];
    int *indices= new int [maxNumEntries];

    Teuchos::ArrayRCP<float> vweights(NumMyVertices);
    Teuchos::ArrayRCP<float> eweights(NumMyEdges);

    // The Graph adapter can take the EpetraCrsGraph as its graph.  But we need
    // to supply the weights.  Vertex weights should be in order of the local vertex
    // ids of the graph.  Edge weights should appear in the order that edges appear in the
    // graph calls.

    int actualNum;
    float *ewgt=eweights.get();
    float *vwgt=vweights.get();

    for (i=0; i <  numMyVertices; i++){
      int gid = mtx->GRID(i);
      mtx->ExtractGlobalRowCopy(gid, maxNumEntries, actualNum, values, indices);
      for (j=0; j < actualNum; j++){
        if (indices[j] == gid)
          *vwgt++ = values[j];
        else
          *ewgt++ = values[j];
      }
   }

   // Create Zoltan2 graph adapter:

   Teuchos::RCP<Epetra::CrsGraph> graph = Teuchos::rcp(&mtx->Graph(), false);

   Zoltan::EpetraCrsGraphInput<float> (graph, vweights, eweights);
    
   
  }

}


