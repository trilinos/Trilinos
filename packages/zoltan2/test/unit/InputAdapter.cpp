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
//
// This test requires MPI.
//
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
#include <Zoltan2_Config.h>

#ifndef HAVE_MPI
int main(int argc, char *argv[])
{
  std::cerr << "Test " << argv[0] << "requires MPI." << std::endl;
  std::cout << "PASS" << std::endl;
}
#else

#include <EpetraExt_CrsMatrixIn.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_MpiComm.h>
#include <Teuchos_RCP.hpp>
#include <Teuchos_OpaqueWrapper.hpp>
#include <Teuchos_MPIComm.hpp>

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

template <typename Ordinal>
  Teuchos::RCP<Teuchos::MpiComm<Ordinal> > 
    &getTeuchosCommunicator(MPI_Comm &comm)
{
  Teuchos::Opaque<MPI_Comm> > handle = Teuchos::opaqueWrapper<MPI_Comm>(comm);

  Teuchos::RCP<Teuchos::Opaque<MPI_Comm> > commPtr = Teuchos::rcp(&handle)

  return Teuchos::rcp(&Teuchos::MpiComm(commPtr));
}

int rank, nprocs;

#if 0
// Read an EpetraCrsMatrix from a file.  Return an Epetra::CrsGraph
// and weights.  Diagonal non-zeros are vertex weights and
// off-diagonal weights are edge weights.

template <typename Scalar>
Teuchos::RCP<Epetra::CrsMatrix> &getEpetraCrsGraphFromFile(
     std::string fname,
     Teuchos::RCP<Epetra::CrsMatrix> &m,
     Scalar *vwgts,
     Scalar *ewgts)
{
  Epetra::CrsMatrix *M;
  Epetra::MpiComm comm(MPI_COMM_WORLD);

  int rc = EpetraExt::MatrixMarketFileToCrsMatrix(fname.c_str(), comm, M);

  if (rc < 0){
    std::cerr << "can't read file " << fname.c_str() << std::endl;
    if 
    exit(1);
  }

  Teuchos::RCP<Epetra::CrsMatrix> A = Teuchos::rcp(M, false);
  return A;
}
#endif


int main(int argc, char *argv[])
{
  // Set up MPI

  Teuchos::GlobalMPISession session(&argc, &argv);
  nprocs = session.getNProc();
  rank = session.getRank();

  Teuchos::RCP<Teuchos::Comm<int> > tcomm = getTeuchosCommunicator(MPI_COMM_WORLD);

  // Complication due to fact that MPI_Comms are implemented
  //  differently in different MPI implementations.

  Teuchos::RCP<Teuchos::OpaqueWrapper<MPI_Comm> > commDef = 
    Teuchos::rcp(Teuchos::opaqueWrapper(MPI_COMM_WORLD));

  Teuchos::MPIComm<int> comm(commDef);

Teuchos::RCP<Teuchos::OpaqueWrapper<MPI_Comm> > rawComm = comm.getRawMpiComm();
Epetra::MpiComm equivComm((*rawComm)());

#if 0
  Epetra::MpiComm comm(MPI_COMM_WORLD);
  comm->myPID(&rank);
  comm->NumProc(&nprocs);
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
#endif

}


#endif // HAVE_MPI
