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
#include <Zoltan2_TestHelpers.hpp>

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

  // Create global identifiers with weights

  lno_t numLocalIds = 10;
  int weightDim = 2;

  gno_t *myIds = new gno_t [numLocalIds];
  scalar_t *weights = new scalar_t [numLocalIds*weightDim];
  gno_t base = rank * numLocalIds * numLocalIds;

  for (lno_t i=0; i < numLocalIds; i++){
    myIds[i] = base+i;
    weights[i*weightDim] = 1.0;
    weights[i*weightDim + 1] = (nprocs-rank) / (i+1);
  }

  // Create a Zoltan2::BasicIdentifierInput object
  // and verify that it is correct

  typedef Zoltan2::BasicUserTypes<scalar_t, gno_t, lno_t, gno_t> userTypes_t;
  scalar_t *weightPtrs[2] = {weights, weights+1};
  int strides[2] = {2,2};

  Zoltan2::BasicIdentifierInput<userTypes_t> ia( numLocalIds, weightDim, myIds,
    &weightPtrs[0], &strides[0]);

  if (rank == 0)
    std::cout << "Testing " << ia.inputAdapterName() << std::endl;

  if (!fail && ia.getLocalNumberOfIdentifiers() != numLocalIds){
    fail = 4;
  }

  if (!fail && ia.getNumberOfWeights() != weightDim)
    fail = 5;

  const gno_t *globalIdsIn;
  scalar_t const *weightsIn[2];
  int weightStridesIn[2];

  if (!fail && ia.getIdentifierList(globalIdsIn) != numLocalIds)
    fail = 6;

  for (int w=0; !fail && w < weightDim; w++){
    if (ia.getIdentifierWeights(w, weightsIn[w], weightStridesIn[w]) <
        numLocalIds * weightStridesIn[w])
      fail = 20;
  }

  const scalar_t *w1 = weightsIn[0];
  const scalar_t *w2 = weightsIn[1];
  int incr1 = weightStridesIn[0];
  int incr2 = weightStridesIn[1];

  for (lno_t i=0; !fail && i < numLocalIds; i++){

    if (globalIdsIn[i] != base+i)
      fail = 8;
    
    if (!fail && w1[i*incr1] != 1.0)
      fail = 9;
    
    if (!fail && w2[i*incr2] != weights[i*weightDim+1])
      fail = 10;
  }

  delete [] myIds;
  delete [] weights;

  gfail = globalFail(comm, fail);
  if (gfail)
    printFailureCode(comm, fail);   // will exit(1)

  if (rank == 0)
    std::cout << "PASS" << std::endl;
}

