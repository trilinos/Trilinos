#include <iostream>

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>
#include "Teuchos_VerboseObject.hpp"
#include <Teuchos_FancyOStream.hpp>

// Epetra
#include <Epetra_MpiComm.h>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_Import.h>

// Tpetra
#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_Import.hpp>

using Teuchos::RCP;
using Teuchos::ArrayRCP;

typedef int LO;
typedef int GO;

/**********************************************************************************/
// Maps of the import operation:
// -----------------------------
// SRC Map  Processor 0: Global IDs = 0 1 2 3 4 5 6  
//          Processor 1: Global IDs =                    9 10 11 12 13 14 15
//
// DEST Map Processor 0: Global IDs =  0 1 2 3 4 5 6 7 8
//          Processor 1: Global IDs =                7 8 9 10 11 12 13 14 15
//
//
// Vectors before import operation:
// --------------------------------
// srcVector  =  0  0  0  0  0 ...  0
// destVector = -1 -1 -1 -1 -1 ... -1
//
// Expected result:
// ----------------
// destVector Processor 0: Values = 0 0 0 0 0 0 0 -1 -1
//            Processor 1: Values =               -1 -1 0 0 0 0 0 0 0
ArrayRCP<GO> TEST1_srcGID(int MyPid) {
  ArrayRCP<GO> srcGID(7);
  if (MyPid == 0) {
    for(int i=0;i<7; i++)
      srcGID[i] = i;
  } else {
    for(int i=0;i<7; i++)
      srcGID[i] = i+9;
  }
  return srcGID;
}
ArrayRCP<GO> TEST1_destGID(int MyPid) {
  ArrayRCP<GO> destGID(9);
  if (MyPid == 0) {
    for(int i=0;i<9; i++)
      destGID[i] = i;
  } else {
    for(int i=0;i<9; i++)
      destGID[i] = i+7;
  }
  return destGID;
}

/**********************************************************************************/
// SRC Map  Processor 0: Global IDs = 
//          Processor 1: Global IDs = 0 1
//
// DEST Map Processor 0: Global IDs = 0 1 2
//          Processor 1: Global IDs =   1 2 
//
// Expected result:
// destVector Processor 0: Values = 0 0 -1
//            Processor 1: Values =   0 -1
ArrayRCP<GO> TEST2_srcGID(int MyPid) {
  if (MyPid == 0) {
    return ArrayRCP<GO>(0);
  } else {
    ArrayRCP<GO> srcGID(2);
    for(int i=0;i<2; i++)
      srcGID[i] = i;
    return srcGID;
  }
}
ArrayRCP<GO> TEST2_destGID(int MyPid) {
  if (MyPid == 0) {
    ArrayRCP<GO> destGID(3);
    for(int i=0;i<3; i++)
      destGID[i] = i;
    return destGID;
  } else {
    ArrayRCP<GO> destGID(2);
    for(int i=0;i<2; i++)
      destGID[i] = i+1;
    return destGID;
  }
}

/**********************************************************************************/
// SRC Map  Processor 0: Global IDs = 
//          Processor 1: Global IDs = 
//
// DEST Map Processor 0: Global IDs = 0 1
//          Processor 1: Global IDs = 0 1 
//
// Expected result:
// destVector Processor 0: Values = -1 -1
//            Processor 1: Values = -1 -1
ArrayRCP<GO> TEST3_srcGID(int MyPid) {
  ArrayRCP<GO> srcGID(0);
  return srcGID;
}
ArrayRCP<GO> TEST3_destGID(int MyPid) {
  ArrayRCP<GO> destGID(2);
  for(int i=0;i<2; i++)
    destGID[i] = i;
  return destGID;
}

/**********************************************************************************/

void TestEpetra(ArrayRCP<GO> srcGID, ArrayRCP<GO> destGID) {
  Epetra_MpiComm comm(MPI_COMM_WORLD);

  Epetra_Map srcMap (-1,  srcGID.size(),  srcGID.getRawPtr(), 0, comm);
  Epetra_Map destMap(-1, destGID.size(), destGID.getRawPtr(), 0, comm);
      
  Epetra_Vector srcVector (srcMap);
  Epetra_Vector destVector(destMap);
  destVector.PutScalar(-1);

  Epetra_Import importer(destMap, srcMap);
  int err = destVector.Import(srcVector, importer, Insert);

  TEST_FOR_EXCEPTION(err != 0,  std::runtime_error, "Epetra_Vector.Import() returned "+err);

  destVector.Print(std::cout);
}

/**********************************************************************************/

void TestTpetra(ArrayRCP<GO> srcGID, ArrayRCP<GO> destGID) {
  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  const Tpetra::global_size_t GSTI = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();

  RCP<Tpetra::Map<LO> >  srcMap = rcp( new Tpetra::Map<LO>(GSTI,  srcGID(), 0, comm));
  RCP<Tpetra::Map<LO> > destMap = rcp( new Tpetra::Map<LO>(GSTI, destGID(), 0, comm));
      
  Tpetra::Vector<int>  srcVector(srcMap);
  Tpetra::Vector<int> destVector(destMap);
  destVector.putScalar(-1);

  Tpetra::Import<LO> importer(srcMap, destMap);
  destVector.doImport(srcVector, importer, Tpetra::INSERT);

  Teuchos::FancyOStream out(Teuchos::rcp(&std::cout,false));
  destVector.describe(out, Teuchos::VERB_EXTREME);
}

/**********************************************************************************/

int main(int argc, char *argv[]) {
  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);
  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  int MyPid = comm->getRank();

  TestEpetra(TEST1_srcGID(MyPid), TEST1_destGID(MyPid));
  TestEpetra(TEST2_srcGID(MyPid), TEST2_destGID(MyPid));
  TestEpetra(TEST3_srcGID(MyPid), TEST3_destGID(MyPid));

  TestTpetra(TEST1_srcGID(MyPid), TEST1_destGID(MyPid)); // Exception: Tpetra::Distributor::doPosts(): imports must be large enough to store the imported data.
  TestTpetra(TEST2_srcGID(MyPid), TEST2_destGID(MyPid)); // Exception: Tpetra::Import::setupExport(): Target has GIDs not found in Source.
  TestTpetra(TEST3_srcGID(MyPid), TEST3_destGID(MyPid)); // Exception: Tpetra::Import<int, int, Kokkos::SerialNode>::setupSamePermuteRemote(): Target has remote LIDs but Source is not distributed globally.

  return EXIT_SUCCESS;
}
