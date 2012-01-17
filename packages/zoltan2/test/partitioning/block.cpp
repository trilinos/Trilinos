#include <Zoltan2_BasicIdentifierInput.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_TestHelpers.hpp>

#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>

using namespace std;
using Teuchos::Comm;
using Teuchos::RCP;
using Teuchos::ParameterList;

int main(int argc, char **argv)
{
  int fail=0, gfail=0;
  Teuchos::GlobalMPISession session(&argc, &argv);
  RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  int rank = comm->getRank();
  int nprocs = comm->getSize();

  int numGlobalIdentifiers = 100;
  int numMyIdentifiers = numGlobalIdentifiers / nprocs;
  if (rank < numGlobalIdentifiers % nprocs)
    numMyIdentifiers += 1;

  int myBaseId = numGlobalIdentifiers * rank;

  int weightDim = 1;
  int *myIds = new int [numMyIdentifiers];
  float *myWeights = new float [numMyIdentifiers*weightDim];

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

  for (int i=0; i < numMyIdentifiers; i++){
    myIds[i] = myBaseId+i;
    myWeights[i] = rank%3 + 1;
  }

  float **weightPointers = new float * [weightDim];
  weightPointers[0] = myWeights;

  typedef Zoltan2::BasicUserTypes<float, int, int, int> mydata_t;
  typedef Zoltan2::BasicIdentifierInput<mydata_t> adapter_t;

  adapter_t adapter(numMyIdentifiers, weightDim, myIds, weightPointers, NULL);

  ParameterList params("test parameters");
  ParameterList &partitioningParams = params.sublist("partitioning");
  partitioningParams.set("num_local_parts", size_t(1));
  partitioningParams.set("algorithm", "block");
  partitioningParams.set("approach", "partition");
  
  Zoltan2::PartitioningProblem<adapter_t> problem(&adapter, &params);

  problem.solve();

  Zoltan2::PartitioningSolution<mydata_t> solution = problem.getSolution();

  float *totalWeight = new float [nprocs];
  float *sumWeight = new float [nprocs];
  memset(totalWeight, 0, nprocs * sizeof(float));

  const int *idList = solution.getGlobalIdList();
  const size_t *partList = solution.getPartList();
  const float *metrics = solution.getImbalance();

  for (int i=0; !fail && i < numMyIdentifiers; i++){
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

  Teuchos::reduceAll<int, float>(*comm, Teuchos::REDUCE_SUM, nprocs, 
    totalWeight, sumWeight);

  double epsilon = 10e-6;

  if (rank == 0){
    std::cout << "Part weights: ";
    float total = 0;
    for (int i=0; i < nprocs; i++){
      std::cout << sumWeight[i] << " ";
      total += sumWeight[i];
    }
    std::cout << std::endl;

    float avg = total / float(nprocs);

    float imbalance = -1.0;

    for (int i=0; i < nprocs; i++){
      float imb = 0;
      if (sumWeight[i] > avg)
        imb = (sumWeight[i] - avg) / avg;
      else
        imb = (avg - sumWeight[i]) / avg;

      if (imb > imbalance)
        imbalance = imb;
    }
    imbalance += 1.0;

    std::cout << "Computed imbalance: " << imbalance << std::endl;
    std::cout << "Returned imbalance: " << metrics[0] << std::endl;

    double err;
    if (imbalance > metrics[0])
      err = imbalance - metrics[0];
    else
      err = metrics[0] - imbalance;

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
