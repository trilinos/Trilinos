#include <Teuchos_Comm.hpp>

#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include <Epetra_SerialComm.h>

#include <Cthulhu_EpetraComm.hpp>

#include "Tpetra_DefaultPlatform.hpp"

// This driver simply tests Teuchos2Epetra_Comm

int main(int argc, char** argv) 
{
  typedef int                                                       Ordinal;
  typedef double                                                    Scalar;

  typedef Tpetra::MpiPlatform<Kokkos::SerialNode>                   MpiPlatform;
  typedef Tpetra::SerialPlatform<Kokkos::SerialNode>                SerialPlatform;

  typedef MpiPlatform::NodeType                                     MpiNodeType;
  typedef SerialPlatform::NodeType                                  SerialNodeType;

  using namespace Teuchos;

  oblackholestream blackhole;
  GlobalMPISession mpiSession(&argc,&argv,&blackhole);

  ParameterList pl;

  Ordinal numThreads=1;
  pl.set("Num Threads",numThreads);
  RCP<MpiNodeType> mpiNode = rcp(new MpiNodeType(pl));
  RCP<SerialNodeType> serialNode = rcp(new SerialNodeType(pl));

  MpiPlatform    myMpiPlat(mpiNode);
  SerialPlatform mySerialPlat(serialNode);

  {
    RCP<const Comm<int> > teuchosComm = mySerialPlat.getComm();
    RCP<const Epetra_Comm> epetraComm = Teuchos2Epetra_Comm(teuchosComm);

    assert(epetraComm != Teuchos::null);
  }

  {
    RCP<const Comm<int> > teuchosComm = myMpiPlat.getComm();
    RCP<const Epetra_Comm> epetraComm = Teuchos2Epetra_Comm(teuchosComm);

    assert(epetraComm != Teuchos::null);
  }

  return(0);
} //main
