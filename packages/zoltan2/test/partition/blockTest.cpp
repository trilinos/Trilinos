#include <Zoltan2_BasicIdentifierInput.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_TestHelpers.hpp>

#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>

using namespace std;
using Teuchos::Comm;
using Teuchos::RCP;

int main(int argc, char **argv)
{
  int fail=0, gfail=0;
  Teuchos::GlobalMPISession session(&argc, &argv);
  RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  int rank = comm->getRank();
  int nprocs = comm->getSize();

  gno_t numGlobalIdentifiers = 100;
  lno_t numMyIdentifiers = numGlobalIdentifiers / nprocs;
  if (rank < numGlobalIdentifiers % nprocs)
    numMyIdentifiers += 1;

  gno_t myBaseId = numGlobalIdentifiers * rank;

  int weightDim = 1;
  gno_t *myIds = new gno_t [numMyIdentifiers];
  scalar_t *myWeights = new scalar_t [numMyIdentifiers*weightDim];

  if (!myIds || !myWeights){
    fail = 1;
  }

  gfail = globalFail(comm, fail);

  if (gfail){
    if (rank==0){
      std::cout << "Memory allocation failure" << std::endl;
      std::cout << "FAIL" << std:: endl;
    }
    return 1;
  }

  for (lno_t i=0; i < numMyIdentifiers; i++){
    myIds[i] = myBaseId+i;
    myWeights[i] = rank%3 + 1;
  }

  std::vector<const scalar_t *> weightValues;
  std::vector<int> weightStrides;   // default is one
  weightValues.push_back(const_cast<const scalar_t *>(myWeights));

  typedef Zoltan2::BasicUserTypes<scalar_t, gno_t, lno_t, gno_t> mydata_t;
  typedef Zoltan2::BasicIdentifierInput<mydata_t> adapter_t;

  adapter_t adapter(numMyIdentifiers, myIds, weightValues, weightStrides);

  Teuchos::ParameterList params("test parameters");
  params.set("compute_metrics", "true");
  Teuchos::ParameterList &partitioningParams = params.sublist("partitioning");
  partitioningParams.set("num_local_parts", 1);
  partitioningParams.set("algorithm", "block");
  partitioningParams.set("approach", "partition");
  
  Zoltan2::PartitioningProblem<adapter_t> problem(&adapter, &params);

  problem.solve();

  Zoltan2::PartitioningSolution<adapter_t> solution = problem.getSolution();

  scalar_t *totalWeight = new scalar_t [nprocs];
  scalar_t *sumWeight = new scalar_t [nprocs];
  memset(totalWeight, 0, nprocs * sizeof(scalar_t));

  const gno_t *idList = solution.getIdList();
  const zoltan2_partId_t *partList = solution.getPartList();
  const scalar_t libImbalance = problem.getImbalance();

  for (lno_t i=0; !fail && i < numMyIdentifiers; i++){
    if (idList[i] != myIds[i])
      fail = 1;

    if (!fail)
      totalWeight[partList[i]] += myWeights[i];
  }

  gfail = globalFail(comm, fail);

  if (gfail){
    if (rank==0){
      std::cout << "failure in solution data" << std::endl;
      std::cout << "FAIL" << std:: endl;
    }
    return 1;
  }

  Teuchos::reduceAll<int, scalar_t>(*comm, Teuchos::REDUCE_SUM, nprocs, 
    totalWeight, sumWeight);

  double epsilon = 10e-6;

  if (rank == 0){
    std::cout << "Part weights: ";
    scalar_t total = 0;
    for (int i=0; i < nprocs; i++){
      std::cout << sumWeight[i] << " ";
      total += sumWeight[i];
    }
    std::cout << std::endl;

    scalar_t avg = total / scalar_t(nprocs);

    scalar_t imbalance = -1.0;

    for (int i=0; i < nprocs; i++){
      scalar_t imb = 0;
      if (sumWeight[i] > avg)
        imb = (sumWeight[i] - avg) / avg;
      else
        imb = (avg - sumWeight[i]) / avg;

      if (imb > imbalance)
        imbalance = imb;
    }
    imbalance += 1.0;

    std::cout << "Computed imbalance: " << imbalance << std::endl;
    std::cout << "Library's imbalance: " << libImbalance << std::endl;

    double err;
    if (imbalance > libImbalance)
      err = imbalance - libImbalance;
    else
      err = libImbalance - imbalance;

    if (err > epsilon)
      fail = 1;
  }
  else{
      fail = 0;
  }

  gfail = globalFail(comm, fail);

  if (gfail){
    if (rank==0){
      std::cout << "failure in solution's imbalance data" << std::endl;
      std::cout << "FAIL" << std:: endl;
    }
    return 1;
  }

  if (rank==0)
    std::cout << "PASS" << std:: endl;

  return 0;
}
