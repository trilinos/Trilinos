// @HEADER
// ***********************************************************************
//                Copyright message goes here.
// ***********************************************************************
//
// Test for Zoltan2::BasicCoordinateInput

/*! \file Metric.cpp
 *
 *  \brief Tests the Zoltan2_Metric.hpp functions.
 *  \todo test for exceptions
 *  \todo We need to check the correctness of the functions.
 */

#include <Zoltan2_TestHelpers.hpp>
#include <Zoltan2_Metric.hpp>

#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>

using Teuchos::Array;
using Teuchos::ArrayRCP;
using Teuchos::RCP;
using Teuchos::ArrayView;
using Teuchos::rcp;

using Zoltan2::StridedData;
using Zoltan2::Environment;
using Zoltan2::MetricValues;
using Zoltan2::objectMetrics;
using Zoltan2::printMetrics;

using namespace std;

void MetricTest(RCP<const Teuchos::Comm<int> > &comm)
{
  int nprocs = comm->getSize();
  RCP<const Environment> env = rcp(new Environment);
  
  int numLocalObj = 10;
  int numGlobalParts = nprocs;

  // Create three different weight arrays.

  Array<scalar_t> weights(3*numLocalObj);
  scalar_t *wgt = weights.getRawPtr();

  for (int i=0; i < numLocalObj; i++)
    wgt[i] = (i%2) + 1.0;

  scalar_t *upDownWeights = wgt;

  wgt += numLocalObj;
  for (int i=0; i < numLocalObj; i++)
    wgt[i] = 1.0;

  scalar_t *uniformWeights = wgt;

  wgt += numLocalObj;
  for (int i=0; i < numLocalObj; i++)
    wgt[i] = (i < numLocalObj/2)? 1.0 : 2.0;

  scalar_t *heavyHalfWeights = wgt;

  typedef StridedData<lno_t, scalar_t> strided_t;

  strided_t w1(ArrayView<scalar_t>(upDownWeights, numLocalObj), 1);
  strided_t w2(ArrayView<scalar_t>(uniformWeights, numLocalObj), 1);
  strided_t w3(ArrayView<scalar_t>(heavyHalfWeights, numLocalObj), 1);

  // Create three different part size arrays.

  Array<scalar_t> partSizes(3*numGlobalParts);
  scalar_t *psize = partSizes.getRawPtr();
  scalar_t psizeTotal = 0;

  for (int i=0; i < numGlobalParts; i++)
    psize[i] = 1.0 / numGlobalParts;

  scalar_t *uniformParts = psize;

  psize += numGlobalParts;
  for (int i=0; i < numGlobalParts; i++){
    psize[i] = (i ? 1.0 : 2.0);
    psizeTotal += psize[i];
  }

  for (int i=0; i < numGlobalParts; i++)
    psize[i] /= psizeTotal;

  scalar_t *oneBigPart = psize;

  psize += numGlobalParts;
  psizeTotal = 0;
  for (int i=0; i < numGlobalParts; i++){
    psize[i] = (i%2 ? 1.0 : 2.0);
    psizeTotal += psize[i];
  }

  for (int i=0; i < numGlobalParts; i++)
    psize[i] /= psizeTotal;

  scalar_t *alternatingSizeParts = psize;

  ArrayView<scalar_t> p1(uniformParts, numGlobalParts);
  ArrayView<scalar_t> p2(oneBigPart, numGlobalParts);
  ArrayView<scalar_t> p3(alternatingSizeParts, numGlobalParts);

  // Create a local part assignment array.

  Array<zoltan2_partId_t> parts(numLocalObj);
  for (int i=0; i < numLocalObj; i++){
    parts[i] = i % numGlobalParts;
  }

  // A multicriteria norm:
  // normMinimizeTotalWeight, normBalanceTotalMaximum or 
  // normMinimizeMaximumWeight

  Zoltan2::multiCriteriaNorm mcnorm = Zoltan2::normBalanceTotalMaximum;

  /*! \test Multiple non-uniform weights and part sizes
   */

  Array<strided_t> threeWeights(3);
  threeWeights[0] = w1;
  threeWeights[1] = w2;
  threeWeights[2] = w2;


  Array<ArrayView< scalar_t> > threePartSizes(3);
  threePartSizes[0] = p1;
  threePartSizes[1] = p2;
  threePartSizes[2] = p3;

  zoltan2_partId_t numParts, numNonemptyParts;
  ArrayRCP<MetricValues<scalar_t> > metrics;

  objectMetrics<scalar_t, lno_t>(env, comm, 
    numGlobalParts, threePartSizes.view(0,3), parts.view(0,numLocalObj),
    threeWeights.view(0,3), mcnorm,
    numParts, numNonemptyParts, metrics);

  printMetrics(std::cout, 
    numGlobalParts, numParts, numNonemptyParts, 
    metrics.view(0,metrics.size()));

  /*! \test Multiple non-uniform weights but uniform part sizes
   */

  mcnorm = Zoltan2::normMinimizeTotalWeight;

  objectMetrics<scalar_t, lno_t>(env, comm, 
    numGlobalParts, parts.view(0,numLocalObj),
    threeWeights.view(0,3), mcnorm,
    numParts, numNonemptyParts, metrics);

  printMetrics(std::cout, 
    numGlobalParts, numParts, numNonemptyParts, 
    metrics.view(0,metrics.size()));


  /*! \test One weight with non-uniform part sizes
   */

  objectMetrics<scalar_t, lno_t>(env, comm, 
    numGlobalParts, parts.view(0,numLocalObj), w1, 
    numParts, numNonemptyParts, metrics);

  printMetrics(std::cout, 
    numGlobalParts, numParts, numNonemptyParts, 
    metrics.view(0,metrics.size()));

  /*! \test One weight with uniform part sizes
   */

  objectMetrics<scalar_t, lno_t>(env, comm, 
    numGlobalParts, parts.view(0,numLocalObj), w1, 
    numParts, numNonemptyParts, metrics);

  printMetrics(std::cout, 
    numGlobalParts, numParts, numNonemptyParts, 
    metrics.view(0,metrics.size()));

  /*! \test Weights and part sizes are uniform.
   */

  objectMetrics<scalar_t, lno_t>(env, comm, numGlobalParts, 
    parts.view(0,numLocalObj),
    numParts, numNonemptyParts, metrics);

  printMetrics(std::cout, 
    numGlobalParts, numParts, numNonemptyParts, 
    metrics.view(0,metrics.size()));


  /*! \test Target number of parts is greater than current, non-uniform
   *            weights and part sizes.
   */

  int targetNumParts = numGlobalParts + 2;
  ArrayView<scalar_t> newPartSizeList(partSizes.getRawPtr(), targetNumParts);
  psizeTotal = 0;
  for (int i=0; i < targetNumParts; i++)
    psizeTotal += newPartSizeList[i];
  for (int i=0; i < targetNumParts; i++)
    newPartSizeList[i] /= psizeTotal;

  objectMetrics<scalar_t, lno_t>(env, comm, 
    targetNumParts, newPartSizeList, 
    parts.view(0,numLocalObj), w1, 
    numParts, numNonemptyParts, metrics);

  printMetrics(std::cout, 
    targetNumParts, numParts, numNonemptyParts, 
    metrics.view(0,metrics.size()));


  /*! \test Target number of parts is greater than current, uniform
   *            weights and part sizes.
   */

  objectMetrics<scalar_t, lno_t>(env, comm, targetNumParts, 
    parts.view(0,numLocalObj),
    numParts, numNonemptyParts, metrics);

  printMetrics(std::cout, 
    targetNumParts, numParts, numNonemptyParts, 
    metrics.view(0,metrics.size()));

  if (nprocs > 2){

    /*! \test Target number of parts is less than current, non-uniform
     *            weights and part sizes.
     */
  
    int targetNumParts = numGlobalParts - 1;
    newPartSizeList = ArrayView<scalar_t>(
      partSizes.getRawPtr(), targetNumParts);
    psizeTotal = 0;
    for (int i=0; i < targetNumParts; i++)
      psizeTotal += newPartSizeList[i];
    for (int i=0; i < targetNumParts; i++)
      newPartSizeList[i] /= psizeTotal;
  
    objectMetrics<scalar_t, lno_t>(env, comm, targetNumParts, newPartSizeList,
      parts.view(0,numLocalObj), w1, 
      numParts, numNonemptyParts, metrics);

    printMetrics(std::cout, 
      targetNumParts, numParts, numNonemptyParts, 
      metrics.view(0,metrics.size()));
  
    /*! \test Target number of parts is less than current, uniform
     *            weights and part sizes.
     */
  
    objectMetrics<scalar_t, lno_t>(env, comm, targetNumParts,
      parts.view(0,numLocalObj),
      numParts, numNonemptyParts, metrics);

    printMetrics(std::cout, 
      targetNumParts, numParts, numNonemptyParts, 
      metrics.view(0,metrics.size()));
  }

  return;
}

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  Teuchos::RCP<const Teuchos::Comm<int> > comm =
    Teuchos::DefaultComm<int>::getComm();

  MetricTest(comm);

  if (comm->getRank() > 0)
    std::cout << "PASS" << std::endl;

  return 0;
}

