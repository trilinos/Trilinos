#ifndef _ZOLTAN2_ALGBLOCK_HPP_
#define _ZOLTAN2_ALGBLOCK_HPP_

#include <Zoltan2_IdentifierModel.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_Metric.hpp>

#include <Teuchos_ParameterList.hpp>

#include <sstream>
#include <string>

/*! \file Zoltan2_AlgBlock.hpp
 *  \brief The algorithm for block partitioning.
 */

namespace Zoltan2{

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

  ////////////////////////////////////////////////////////
  // Library parameters of interest:
  //
  //    are we printing out debug messages
  //    are we timing
  //    are we computing memory used

  bool debug = env->doStatus();
  bool timing = env->doTiming();
  bool memstats = env->doMemoryProfiling();

  if (debug)
    env->debugOut_->print(DETAILED_STATUS, string("Entering AlgBlock"));

  Z2_GLOBAL_BUG_ASSERTION(*env, "parameters are not committed",
    env->parametersAreCommitted(), DEBUG_MODE_ASSERTION);

  int rank = env->myRank_;
  int nprocs = env->numProcs_;

  // Parameters that may drive algorithm choices
  //    speed_versus_quality
  //    memory_footprint_versus_runtime

  const Teuchos::ParameterList &pl = env->getParams();
  string mvr = pl.get<string>(string("memory_footprint_versus_runtime"));
  string svq = pl.get<string>(string("speed_versus_quality"));

  bool fastSolution = (svq==string("speed"));
  bool goodSolution = (svq==string("quality"));
  bool balancedSolution = (svq==string("balance"));
  
  bool lowMemory = (mvr==string("memory")) || (mvr==string("memory_footprint"));
  bool lowRunTime = (mvr==string("runtime"));
  bool balanceMemoryRunTime = (svq==string("balance"));

  ////////////////////////////////////////////////////////
  // Problem parameters of interest:
  //    objective
  //    imbalance_tolerance

  string objective("balance_object_weight");
  double imbalanceTolerance = 1.1;

  if (env->hasPartitioningParameters()){
    const Teuchos::ParameterList &pl = env->getPartitioningParams();
    objective = pl.get<string>(string("objective"));
    imbalanceTolerance = pl.get<double>(string("imbalance_tolerance"));
  }

  ////////////////////////////////////////////////////////
  // From the IdentifierModel we need:
  //    the number of gnos
  //    number of weights per gno
  //    the weights
  // TODO: modify algorithm for weight dimension greater than 1.

  size_t numGnos = ids->getLocalNumIdentifiers();
  int wtflag = ids->getIdentifierWeightDim();

  int weightDim = (wtflag ? wtflag : 1);

  ArrayView<const gno_t> idList;
  ArrayView<const StridedInput<lno_t, scalar_t> > wgtList;
  
  ids->getIdentifierList(idList, wgtList);

  ////////////////////////////////////////////////////////
  // From the Solution we get part information.
  //
  //   TODO: for now, we have 1 part per proc and all
  //   part sizes are the same.

  size_t numGlobalParts = solution->getGlobalNumberOfParts();
  size_t numLocalParts = solution->getLocalNumberOfParts();

  //const int *partDist = solution->getPartDistribution();
  //const size_t *procDist = solution->getProcsParts();
  //double *partSizes0 = getCriteriaPartSizes(0);
  
  ////////////////////////////////////////////////////////
  // What is the objective.
  // TODO: consider the objective in the algorithm

  enum {object_weight,     // balance on first weight
        object_count,      // ignore weight, balance on object count
        total_weight,      // balance on total of weights
        maximum_weight};   // balance on maximum of weights

  int balanceObjective = (wtflag ? object_weight : object_count);

  if (wtflag > 1){
    if (objective == string("balance_object_count"))
      balanceObjective = object_count;
    else if (objective == string("multicriteria_minimize_total_weight"))
      balanceObjective = total_weight;
    else if (objective == string("multicriteria_minimize_maximum_weight"))
      balanceObjective = maximum_weight;
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

  /* Overwrite part_sizes with cumulative sum (inclusive) part_sizes. */
  /* A cleaner way is to make a copy, but this works. */

  Array<scalar_t> part_sizes(numGlobalParts, 1.0/numGlobalParts); 
  for (int i=1; i<numGlobalParts; i++)
    part_sizes[i] += part_sizes[i-1];

  if (debug){
    ostringstream oss("Part sizes: ");
    for (int i=0; i <= nprocs; i++)
      oss << part_sizes[i];
    oss << std::endl << "Weights : ";
    for (int i=0; i <= nprocs; i++)
      oss << scansum[i];
    env->debugOut_->print(VERBOSE_DETAILED_STATUS, oss.str());
  }

  /* Loop over objects and assign partition. */
  size_t part = 0;
  wtsum = scansum[rank];
  Array<scalar_t> partTotal(numGlobalParts, 0);
  ArrayRCP<size_t> gnoPart= arcp(new size_t [numGnos], 0, numGnos);

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

  ArrayRCP<float> imbalance = arcp(new float[weightDim], 0, weightDim);

  // TODO - get part sizes from the solution object.  For now, 
  //    an empty part size array means uniform parts.

  ArrayView<float> defaultPartSizes(Teuchos::null);
  Array<ArrayView<float> > partSizes(weightDim, defaultPartSizes);

  // TODO have partNums default to 0 through numGlobalParts-1 in
  //    imbalances() call.
  Array<size_t> partNums(numGlobalParts);
  for (size_t i=0; i < numGlobalParts; i++) partNums[i] = i;

  Array<ArrayView<scalar_t> > partWeights(1);
  partWeights[0] = partTotal.view(0, numGlobalParts);

  try{
    imbalances<scalar_t>(env, problemComm, numGlobalParts, 
      partSizes, partNums.view(0, numGlobalParts),
         partWeights, imbalance.view(0, weightDim));
  }
  Z2_FORWARD_EXCEPTIONS;
  

  ////////////////////////////////////////////////////////////
  // Done
  
  if (debug){
    if (imbalance[0] > Teuchos::as<scalar_t>(imbalanceTolerance)){
      ostringstream oss("Warning: imbalance is ");
      oss << imbalance[0] << std::endl;
      env->debugOut_->print(BASIC_STATUS, oss.str());
    }
    else{
      ostringstream oss("Imbalance: ");
      oss << imbalance[0] << std::endl;
      env->debugOut_->print(DETAILED_STATUS, oss.str());
    }
  }

  // Done, update the solution

  solution->setParts(idList, gnoPart, imbalance);

  if (debug)
    env->debugOut_->print(DETAILED_STATUS, string("Exiting AlgBlock"));
}

}   // namespace Zoltan2

#endif
