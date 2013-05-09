
//#include <stdio.h>
//#include <mpi.h>
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayView.hpp"
#include "Zoltan2_config.h"
#include "Zoltan2_Util.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include <Zoltan2_Environment.hpp>
#include <Zoltan2_Parameters.hpp>
#include <Zoltan2_TestHelpers.hpp>


#include <string>
#include <sstream>
#include <iostream>

typedef int LO;
typedef int GO;
typedef double Scalar;

int main(int narg, char **arg)
{
  Teuchos::GlobalMPISession mpiSession(&narg,&arg);

  Teuchos::RCP<const Teuchos::Comm<int> >
    comm = Teuchos::DefaultComm<int>::getComm();
  int me = comm->getRank();
  int nprocs = comm->getSize();

  if (nprocs != 4)
      std::cout << "Run with 4 MPI ranks " << std::endl;

  typedef Tpetra::Map<LO, GO> map_t;
  GO numGlobalCoords = 4000000;
  LO numLocalCoords = 1000000;
  Teuchos::ParameterList myParams("testParameterList");
  myParams.set("memory_procs", "0");
  myParams.set("memory_output_stream", "std::cout");

  LO newnumLocalCoords = 1000000;
  if (me == 0)
      newnumLocalCoords = 999999;
  else if (me == 1)
      newnumLocalCoords = 1000001;
  else
      newnumLocalCoords = 1000000;


  Zoltan2::Environment *defEnv = NULL;

  try{
    defEnv = new Zoltan2::Environment(myParams, comm);
  }
  catch(std::exception &e){
    std::cerr << e.what() << std::endl;
  }

  typedef Tpetra::MultiVector<Scalar, LO, GO> mvector_t;

  defEnv->memory("Before map construction");
  for (int i = 0 ; i < 1000; i++)
  {
      defEnv->memory("Inside the loop");
      Teuchos::RCP<const map_t> tmap = rcp(new map_t(numGlobalCoords, 
        numLocalCoords, 0, comm));
      Teuchos::RCP<const map_t> newTmap = rcp(new map_t(numGlobalCoords, 
        newnumLocalCoords, 0, comm));
      Teuchos::RCP<mvector_t> newMvector = rcp(new mvector_t(tmap, 3, true));
      RCP<Tpetra::Import<LO, GO> > importer = rcp(
        new Tpetra::Import<LO, GO>(tmap, newTmap));
      //defEnv->memory("Inside the loop after i = 0");
  }
  defEnv->memory("After map construction");


}
