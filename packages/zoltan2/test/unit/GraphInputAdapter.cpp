// @HEADER
// ***********************************************************************
//
//         Zoltan2: Sandia Partitioning Ordering & Coloring Library
//
//                Copyright message goes here.   TODO
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
#include <Zoltan2_config.h>
#include <string>
#include <vector>
#include <iostream>
#include <Zoltan2_Util.hpp>
//#include <Zoltan2_TpetraCrsGraphInput.hpp>
#include <Zoltan2_EpetraCrsGraphInput.hpp>
#include <EpetraExt_CrsMatrixIn.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_SerialComm.h>
#include <Epetra_MpiComm.h>

using namespace std;

int rank, nprocs;

int main(int argc, char *argv[])
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::SerialComm;
  using Teuchos::MpiComm;

#ifdef HAVE_MPI
  MPI_Comm c = MPI_COMM_WORLD;
  RCP<MpiComm<int> > comm = Zoltan2::getTeuchosMpiComm<int>(c);
  Epetra_MpiComm ecomm(MPI_COMM_WORLD);
#else
  RCP<SerialComm<int> > comm = rcp(new SerialComm<int>);
  Epetra_SerialComm ecomm;
#endif

  rank = comm->getRank();
  nprocs = comm->getSize();

  std::vector<std::string> mtxFiles;
  mtxFiles.push_back("ewgt.mtx");
  mtxFiles.push_back("simple.mtx");
  mtxFiles.push_back("cage10.mtx");
  mtxFiles.push_back("diag500_4.mtx");

  for (unsigned int i=0; i < mtxFiles.size(); i++){

    Epetra_CrsMatrix *M;

    int fail = EpetraExt::MatrixMarketFileToCrsMatrix(
      mtxFiles[i].c_str(), ecomm, M, 0, 0);

    if (fail){
      if (rank == 0) {
       std::cerr << "Error with " << mtxFiles[i] << std::cout;
       std::cerr << "EpetraExt::MatrixMarketFileToCrsMatrix" << std::endl;
      }
      comm->barrier();
      exit(1);
    }

    // I need a non-const copy of the graph.  A const reference
    // can not be used in EpetraCrsGraphInput because it can
    // not be used in Xpetra::EpetraCrsGraph.
    // TODO: why can't a const graph be used with Xpetra if we
    // don't call any Xpetra methods that change the graph.

    Epetra_CrsGraph g = M->Graph();

    RCP<Epetra_CrsGraph> graph = rcp(&g);

    Zoltan2::EpetraCrsGraphInput<float> adapter(graph);
  }
}

