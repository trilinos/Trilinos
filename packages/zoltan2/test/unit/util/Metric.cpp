// @HEADER
// ***********************************************************************
//
//         Zoltan2: Sandia Partitioning Ordering & Coloring Library
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
//
// Basic testing of the namespace methods in Zoltan2_Metric.hpp.
//

#include <Zoltan2_Metric.hpp>
#include <Zoltan2_TestHelpers.hpp>

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::Comm;
using Teuchos::DefaultComm;
using Teuchos::Array;
using Teuchos::ArrayView;
using Zoltan2::imbalances;

int main(int argc, char *argv[])
{
  typedef float scalar_t;

  Teuchos::GlobalMPISession session(&argc, &argv);
  RCP<const Comm<int> > comm = DefaultComm<int>::getComm();
  int rank = comm->getRank();
  int nprocs = comm->getSize();
  int fail = 0;
  float epsilon = 1e-5;

  RCP<const Zoltan2::Environment> env = rcp(new Zoltan2::Environment);

  // The PART NUMBERS.  

  size_t numGlobalParts = nprocs;
  Array<size_t> partNums(numGlobalParts);
  for (size_t i=0; i < numGlobalParts; i++)
    partNums[i] = i;

  // The WEIGHTS per part under the current partitioning.

  int weightDim = 3;

  Array<float> passMin(weightDim, 1.0 - epsilon);
  Array<float> passMax(weightDim, 1.0 + epsilon);

  scalar_t *weight = new scalar_t [numGlobalParts * weightDim];

  for (size_t i=0; i < numGlobalParts; i++){
    weight[i] = 1.0;
  }

  ArrayView<scalar_t> w0(weight, numGlobalParts);

  size_t offset = numGlobalParts;

  for (size_t i=0; i < numGlobalParts; i++){
    if (i < numGlobalParts/2)
      weight[offset + i] = 1.0;
    else
      weight[offset + i] = 2.0;
  }

  ArrayView<scalar_t> w1(weight + offset, numGlobalParts);
  if (nprocs > 1){
    passMin[1] = 1.2;
    passMax[1] = 1.4;
  }

  offset += numGlobalParts;
  
  for (size_t i=0; i < numGlobalParts; i++){
    weight[offset + i] = nprocs - rank;
  }

  ArrayView<scalar_t> w2(weight + offset, numGlobalParts);

  Array<ArrayView<scalar_t> > weights;
  weights.push_back(w0);
  weights.push_back(w1);
  weights.push_back(w2);

  // The PART SIZES default to equal parts.

  ArrayView<float> p0(Teuchos::null);   // 3 empty part size arrays
  ArrayView<float> p1(Teuchos::null);   //   implies uniform parts
  ArrayView<float> p2(Teuchos::null);
  
  Array<ArrayView<float> > partSizes;
  partSizes.push_back(p0);
  partSizes.push_back(p1);
  partSizes.push_back(p2);

  // Answer

  Array<float> answer(weightDim, 0.0);

  /////////////////////////////////////////////////////////////
  // Test: multiple weights, uniform part sizes
  /////////////////////////////////////////////////////////////

  imbalances<scalar_t>(env, comm, numGlobalParts, 
      partSizes, partNums, weights, answer.view(0, weightDim));

  if (rank == 0){
    std::cout << "Test 1: " << answer[0] << " ";
    for (int i=1; i < weightDim; i++){
      std::cout << answer[i] << " ";
    }
    std::cout << std::endl;
  }

  fail = 0;

  for (int i=0; !fail && (i < weightDim); i++){
    if (answer[i] < passMin[i] || answer[i] > passMax[i])
      fail=1;
  }

  TEST_FAIL_AND_EXIT(*comm, fail==0, "invalid imbalance", 1);

  /////////////////////////////////////////////////////////////
  // Test: one weight, uniform part sizes
  /////////////////////////////////////////////////////////////

  float result;

  imbalances(env, comm, numGlobalParts, p1, partNums, w1, result);

  if (rank == 0)
    std::cout << "Test 2: " << result << std::endl;

  if (result < passMin[1] || result > passMax[1])
    fail=1;

  TEST_FAIL_AND_EXIT(*comm, fail==0, "invalid imbalance", 1);

  /////////////////////////////////////////////////////////////
  // Test: multiple weights, varying part sizes
  /////////////////////////////////////////////////////////////

  Array<float> varyingPartSizes(numGlobalParts, 1.0);
  float partSizeTotal = 0;
  for (size_t i=0; i < numGlobalParts/3; i++){
    varyingPartSizes[i] = .5;
    partSizeTotal += .5;
  }
  for (size_t i=numGlobalParts/3; i < 2*numGlobalParts/3; i++){
    partSizeTotal += 1.0;
  }
  for (size_t i=2*numGlobalParts/3; i < numGlobalParts; i++){
    varyingPartSizes[i] = 1.5;
    partSizeTotal += 1.5;
  }

  for (size_t i=0; i < numGlobalParts; i++)
    varyingPartSizes[i] /= partSizeTotal;

  partSizes[0] = varyingPartSizes;
  if (nprocs==2){
    passMin[0] = 1.25 - epsilon;
    passMax[0] = 1.25 + epsilon;
  }
  else if (nprocs > 2){
    passMin[0] = 2.0 - epsilon;
    passMax[0] = 2.5;
  }

  imbalances<scalar_t>(env, comm, numGlobalParts, 
      partSizes, partNums, weights, answer.view(0, weightDim));

  if (rank == 0){
    std::cout << "Test 3: " << answer[0] << " ";
    for (int i=1; i < weightDim; i++){
      std::cout << answer[i] << " ";
    }
    std::cout << std::endl;
  }

  fail = 0;

  for (int i=0; !fail && (i < weightDim); i++){
    if (answer[i] < passMin[i] || answer[i] > passMax[i])
      fail=1;
  }

  TEST_FAIL_AND_EXIT(*comm, fail==0, "invalid imbalance", 1);

  /////////////////////////////////////////////////////////////
  // Test: single weight, varying part sizes
  /////////////////////////////////////////////////////////////

  imbalances(env, comm, numGlobalParts, partSizes[0], partNums, w0, result);

  if (rank == 0)
    std::cout << "Test 4: " << result << std::endl;

  if (result < passMin[0] || result > passMax[0])
    fail=1;

  TEST_FAIL_AND_EXIT(*comm, fail==0, "invalid imbalance", 1);

  /////////////////////////////////////////////////////////////
  // Done
  /////////////////////////////////////////////////////////////
  
  if (rank == 0)
    std::cout << "PASS" << std::endl;
}
  
  
