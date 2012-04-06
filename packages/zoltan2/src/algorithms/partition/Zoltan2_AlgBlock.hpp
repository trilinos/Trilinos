#ifndef _ZOLTAN2_ALGBLOCK_HPP
#define _ZOLTAN2_ALGBLOCK_HPP_

#include <Zoltan2_IdentifierModel.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_Metric.hpp>

#include <Teuchos_ParameterList.hpp>

#include <sstream>
#include <string>
#include <bitset>

/*! \file Zoltan2_AlgBlock.hpp
 *  \brief The algorithm for block partitioning.
 */

typedef zoltan2_partId_t partId_t;

namespace Zoltan2{

/*! \brief The boolean parameters of interest to the Block algorithm.
 */
enum blockParams{
  block_balanceCount,            /*!< objective = balance_object_count */
  block_balanceWeight,          /*!< objective = balance_object_weight */
  block_minTotalWeight,      /*!< objective = mc_minimize_total_weight */
  block_minMaximumWeight,  /*!< objective = mc_minimize_maximum_weight */
  block_balanceTotalMaximum, /*!< objective = mc_balance_total_maximum */
  NUM_BLOCK_PARAMS
};

/*! Block partitioning method.
 *
 *  \param env   library configuration and problem parameters
 *  \param problemComm  the communicator for the problem
 *  \param ids    an Identifier model
 *  \param solution  a Solution object, containing part information
 *
 *  Preconditions: The parameters in the environment have been
 *    processed (committed).  No special requirements on the
 *    identifiers.
 *
 *   \todo Block partitioning uses one weight only
 *   \todo Block partitioning assumes one part per process.
 *   \todo check for memory allocation failures
 */


template <typename Adapter>
void AlgPTBlock(
  const RCP<const Environment> &env,
  const RCP<const Comm<int> > &problemComm,
  const RCP<const IdentifierModel<Adapter> > &ids, 
  RCP<PartitioningSolution<typename Adapter::user_t> > &solution
) 
{
  using std::string;
  using std::ostringstream;
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::scalar_t scalar_t;

  if (env->doStatus())
    env->debug(DETAILED_STATUS, string("Entering AlgBlock"));

  int rank = env->myRank_;
  int nprocs = env->numProcs_;

  ////////////////////////////////////////////////////////
  // Partitioning problem parameters of interest:
  //    objective
  //    imbalance_tolerance

  std::bitset<NUM_BLOCK_PARAMS> params;
  double imbalanceTolerance=0.0;
  bool isSet;
  string strChoice;

  env->getValue<string>(
     env->getList(env->getParameters(), "partitioning"),
    "objective", isSet, strChoice);

  if (isSet && strChoice == string("balance_object_count"))
    params.set(block_balanceCount);
  else if (isSet && strChoice ==
    string("multicriteria_minimize_total_weight"))
    params.set(block_minTotalWeight);
  else if (isSet && strChoice ==
    string("multicriteria_minimize_maximum_weight"))
    params.set(block_minMaximumWeight);
  else if (isSet && strChoice ==
    string("multicriteria_balance_total_maximum"))
    params.set(block_balanceTotalMaximum);
  else
    params.set(block_balanceWeight);

  env->getValue<double>(
     env->getList(env->getParameters(), "partitioning"),
    "imbalance_tolerance", isSet, imbalanceTolerance);

  if (!isSet)
    imbalanceTolerance = 1.1;

  ////////////////////////////////////////////////////////
  // From the IdentifierModel we need:
  //    the number of gnos
  //    number of weights per gno
  //    the weights

  size_t numGnos = ids->getLocalNumIdentifiers();
  int wtflag = ids->getIdentifierWeightDim();

  int weightDim = (wtflag ? wtflag : 1);

  ArrayView<const gno_t> idList;
  ArrayView<StridedData<lno_t, scalar_t> > wgtList;
  
  ids->getIdentifierList(idList, wgtList);

  ////////////////////////////////////////////////////////
  // From the Solution we get part information.

  size_t numGlobalParts = solution->getGlobalNumberOfParts();

  size_t numLocalParts = solution->getLocalNumberOfParts();
  if (numLocalParts != 1){
  }
  
  ////////////////////////////////////////////////////////
  // The algorithm
  //
  // Block partitioning algorithm lifted from zoltan/src/simple/block.c
  // The solution is:
  //    a list of part numbers in gno order
  //    an imbalance for each weight 

  scalar_t wtsum(0);

  if (wtflag){
    for (size_t i=0; i<numGnos; i++)
      wtsum += wgtList[0][i];          // [] operator knows stride
  }
  else
    wtsum = static_cast<scalar_t>(numGnos);

  Array<scalar_t> scansum(nprocs+1, 0);

  Teuchos::gatherAll<int, scalar_t>(*problemComm, 1, &wtsum, nprocs,
    scansum.getRawPtr()+1);

  /* scansum = sum of weights on lower processors, excluding self. */

  for (int i=2; i<=nprocs; i++)
    scansum[i] += scansum[i-1];

  scalar_t globalTotalWeight = scansum[nprocs];

  /* part_sizes - inclusive sum */

  Array<scalar_t> part_sizes(numGlobalParts);

  if (!solution->criteriaHasUniformPartSizes(0)){
    part_sizes[0] = solution->getCriteriaPartSize(0, 0);
    for (int i=1; i<numGlobalParts; i++)
      part_sizes[i] = part_sizes[i-1] + solution->getCriteriaPartSize(i, 0);
  }
  else{
    scalar_t onePart = 1.0/numGlobalParts;
    part_sizes[0] = onePart;
    for (int i=1; i<numGlobalParts; i++)
      part_sizes[i] += onePart;
  }

  // TODO assertion that last part sizes is about equal to 1.0

  if (env->doStatus()){
    ostringstream oss("Part sizes: ");
    for (int i=0; i < numGlobalParts; i++)
      oss << part_sizes[i] << " ";
    oss << "\n";
    oss << std::endl << "Weights : ";
    for (int i=0; i <= nprocs; i++)
      oss << scansum[i] << " ";
    oss << "\n";
    env->debug(VERBOSE_DETAILED_STATUS, oss.str());
  }

  /* Loop over objects and assign partition. */
  partId_t part = 0;
  wtsum = scansum[rank];
  Array<scalar_t> partTotal(numGlobalParts, 0);
  ArrayRCP<partId_t> gnoPart= arcp(new partId_t [numGnos], 0, numGnos);

  for (size_t i=0; i<numGnos; i++){
    scalar_t gnoWeight = (wtflag? wgtList[0][i] : 1.0);
    /* wtsum is now sum of all lower-ordered object */
    /* determine new partition number for this object,
       using the "center of gravity" */
    while (part<numGlobalParts-1 && 
           (wtsum+0.5*gnoWeight) > part_sizes[part]*globalTotalWeight)
      part++;
    gnoPart[i] = part;
    partTotal[part] += gnoWeight;
    wtsum += gnoWeight;
  }

  ////////////////////////////////////////////////////////////
  // Compute the imbalance.

  ArrayRCP<MetricValues<scalar_t> > metrics;
  partId_t numParts;
  partId_t numNonemptyParts;

  try{
    objectMetrics<scalar_t, partId_t, lno_t, gno_t>(
      env, problemComm, numGlobalParts, 
      gnoPart.view(0, numGnos),  wgtList[0],
      numParts, numNonemptyParts, metrics);
  }
  Z2_FORWARD_EXCEPTIONS;

  if (metrics.size() == 2)
    wgtImbalance = metrics[1].getMaxImbalance();
  else
    wgtImbalance = metrics[0].getMaxImbalance();

  ////////////////////////////////////////////////////////////
  // Done
  
  if (env->doStatus()){
#if 0
    if (wgtImbalance > Teuchos::as<scalar_t>(imbalanceTolerance)){
      ostringstream oss("Warning: imbalance is ");
      oss << wgtImbalance << std::endl;
      env->debug(BASIC_STATUS, oss.str());
    }
    else{
#endif
      ostringstream oss("Imbalance: ");
      oss << wgtImbalance << std::endl;
      env->debug(DETAILED_STATUS, oss.str());
#if 0
    }
#endif
  }

  // Done, update the solution TODO - compute the metrics
  //    with objectMetrics() call

  ArrayRCP<MetricValues<scalar_t> > emptyMetrics =
    arcp(new MetricValues<scalar_t> [2], 0, 2);

  ArrayRCP<const gno_t> gnos = arcpFromArrayView(idList);

  solution->setParts(gnos, gnoPart, emptyMetrics);

  if (env->doStatus())
    env->debug(DETAILED_STATUS, string("Exiting AlgBlock"));
}

}   // namespace Zoltan2

#endif
