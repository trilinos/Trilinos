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

#include <Zoltan2_Metric.hpp>
#include <ErrorHandlingForTests.hpp>

using Teuchos::RCP;
using Teuchos::Comm;
using Teuchos::DefaultComm;
using Teuchos::Array;
using Teuchos::ArrayView;

int main(int argc, char *argv[])
{
  typedef float scalar_t;

  Teuchos::GlobalMPISession session(&argc, &argv);
  RCP<const Comm<int> > comm = DefaultComm<int>::getComm();
  int rank = comm->getRank();
  int nprocs = comm->getSize();
  int fail = 0, gfail=0;

  RCP<const Zoltan2::Environment> env = Zoltan2::getDefaultEnvironment();

  // The part numbers.  

  size_t numGlobalParts = nprocs;
  Array<size_t> partNums(numGlobalParts);
  for (size_t i=0; i < numGlobalParts; i++)
    partNums[i] = i;

  // The weights under the current partitioning.

  int weightDim = 3;

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

  offset += numGlobalParts;
  
  for (size_t i=0; i < numGlobalParts; i++){
    weight[offset + i] = nprocs - rank;
  }

  ArrayView<scalar_t> w2(weight + offset, numGlobalParts);

  Array<ArrayView<scalar_t> > weights;
  weights.push_back(w0);
  weights.push_back(w1);
  weights.push_back(w2);

  // PART SIZES default to equal parts.

  ArrayView<scalar_t> p0(Teuchos::null);   // 3 empty part size arrays
  ArrayView<scalar_t> p1(Teuchos::null);   //   implies uniform parts
  ArrayView<scalar_t> p2(Teuchos::null);
  
  Array<ArrayView<scalar_t> > partSizes;
  partSizes.push_back(p0);
  partSizes.push_back(p1);
  partSizes.push_back(p2);

  // Answer

  Array<float> answer(weightDim, 0.0);

  // Compute imbalances

  Zoltan2::imbalances<scalar_t>(env, comm, numGlobalParts, 
      partSizes, partNums, weights, answer.view(0, weightDim));

}

  
