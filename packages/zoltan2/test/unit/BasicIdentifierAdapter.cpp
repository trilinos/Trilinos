// @HEADER
// ***********************************************************************
//
//         Zoltan2: Sandia Partitioning Ordering & Coloring Library
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
//
// Basic testing of Zoltan2::BasicIdentifierInput 

#include <Zoltan2_BasicIdentifierInput.hpp>
#include <ErrorHandlingForTests.hpp>

#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_CommHelpers.hpp>

using Teuchos::RCP;
using Teuchos::Comm;
using Teuchos::DefaultComm;

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  RCP<const Comm<int> > comm = DefaultComm<int>::getComm();
  int rank = comm->getRank();
  int nprocs = comm->getSize();
  int fail = 0, gfail=0;

  typedef double scalar_t;
  typedef long   gid_t;
  typedef int    lno_t;
  typedef long   gno_t;

  // Create global identifiers with weights

  int numLocalIds = 10;
  int weightDim = 2;

  gid_t *myIds = new gid_t [numLocalIds];
  scalar_t *weights = new scalar_t [numLocalIds*weightDim];
  gid_t base = rank * numLocalIds * numLocalIds;

  for (lno_t i=0; i < numLocalIds; i++){
    myIds[i] = base+i;
    weights[i*weightDim] = 1.0;
    weights[i*weightDim + 1] = (nprocs-rank) / (i+1);
  }

  // Create a Zoltan2::BasicIdentifierInput object
  // and verify that it is correct

  typedef Zoltan2::BasicUserTypes<scalar_t, gid_t, lno_t, gno_t> userTypes_t;

  Zoltan2::BasicIdentifierInput<userTypes_t> ia(
    numLocalIds, weightDim, myIds, weights);

  if (!fail && ia.getLocalNumIds() != numLocalIds){
    fail = 4;
  }

  if (!fail && ia.getNumWeights() != weightDim)
    fail = 5;

  const gid_t *globalIdsIn;
  const scalar_t *weightsIn;

  if (!fail && ia.getIdList(globalIdsIn, weightsIn) != numLocalIds)
    fail = 6;

  if (!fail){
    for (lno_t i=0; i < numLocalIds; i++){
      if (globalIdsIn[i] != base+i){
        fail = 8;
        break;    
      }
      if (weightsIn[i*weightDim] != 1.0){
        fail = 9;
        break;    
      }
      if (weightsIn[i*weightDim + 1] != weights[i*weightDim+1]){
        fail = 10;
        break;    
      }
    }
  }

  delete [] myIds;
  delete [] weights;

  gfail = globalFail(comm, fail);
  if (gfail)
    printFailureCode(comm, fail);   // will exit(1)

  if (rank == 0)
    std::cout << "PASS" << std::endl;
}

