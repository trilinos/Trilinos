#ifndef _ZOLTAN2_ALGBLOCK_HPP_
#define _ZOLTAN2_ALGBLOCK_HPP_

#include <Zoltan2_IdentifierModel.hpp>
#include <Zoltan2_PartitioningSolution.hpp>

/*! \file Zoltan2_AlgBlock.hpp
 *  \brief The algorithm for block partitioning.
 */

namespace Zoltan2{

/*! Block partitioning method.
 *
 *  \param env  parameters for the problem and library configuration
 *  \param problemComm  the communicator for the problem
 *  \param ids    an Identifier model
 *  \param solution is the Solution object
 *
 *  Preconditions: The parameters in the environment have been
 *    processed (committed).
 */


template <typename Adapter>
void AlgBlock(
  const RCP<const Environment> &env,
  const RCP<const Comm<int> > &problemComm,
  const RCP<const IdentifierModel<Adapter> > &ids, 
  RCP<PartitioningSolution<Adapter::gid_t, Adapter::lno_t, Adapter::gno_t> > &solution
) 
{
  if (debug){
    ostringstream oss << "Entering AlgBlock" << std::endl;
    env->debugOut_->print(DETAILED_STATUS, oss.str());
  }

  using std::string;
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::scalar_t scalar_t;

  Z2_GLOBAL_BUG_ASSERTION(env, "parameters are not committed",
    env->parametersAreCommitted(), DEBUG_MODE_ASSERTION);

  int rank = env->myRank_;
  int nprocs = env->numProcs_;

  ////////////////////////////////////////////////////////
  // Library parameters of interest:
  //
  //    are we printing out debug messages
  //    are we timing
  //    are we computing memory used
  //
  //    Parameters that may drive algorithm choices:
  //      speed_versus_quality
  //      memory_footprint_versus_runtime

  bool debug = env->doStatus();
  bool timing = env->doTiming();
  bool memstats = env->doMemoryProfiling();

  const Teuchos::ParameterList &pl = env->getParams();
  string &mvr = pl.get(string("memory_footprint_versus_runtime"));
  string &svq = pl.get(string("speed_versus_quality"));

  bool fastSolution = (svq==string("speed"));
  bool goodSolution = (svq==string("quality"));
  bool balancedSolution = (svq==string("balance"));
  
  bool lowMemory = (mvr==string("memory") || (mvr==string("memory_footprint"));
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
    objective = *(pl.getPtr(string("objective")));
    imbalanceTolerance = *(pl.getPtr(string("imbalance_tolerance")));
  }

  ////////////////////////////////////////////////////////
  // From the IdentifierModel we need:
  //    the number of gnos
  //    number of weights per gno
  //    the weights
  // TODO: modify algorithm for weight dimension greater than 1.

  size_t numGnos = getLocalNumIdentifiers();
  int wtflag = ids->getIdentifierWeightDim();

  int numWeights = (wtflag ? wtflag : 1);

  ArrayView<const gno_t> idList;
  ArrayView <const StridedInput> wgtList;
  
  ids->getIdentifierList(idList, wgtList);

  const StridedInput weight0 = wgtList[0];  // the first weight

  ////////////////////////////////////////////////////////
  // From the Solution we get part information.
  //
  //   TODO: for now, we have 1 part per proc and all
  //   part sizes are the same.

  size_t &numGlobalParts = solution->getGlobalNumberOfParts();
  size_t &numLocalParts = solution->getLocalNumberOfParts();

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
      wtsum += weight0[i];          // [] operator knows stride
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
  for (i=1; i<numGlobalParts; i++)
    part_sizes[i] += part_sizes[i-1];

  if (debug){
    ostringstream oss << "Part sizes: ";
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

  for (size_t i=0; i<numGnos; i++){
    scalar_t gnoWeight = (wtflag? weight0[i] : 1.0);
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

  // Compute the imbalance - taken from Zoltan's object_metrics function.

  ArrayRCP<scalar_t> globalPartTotal(nprocs); 

  Teuchos::reduceAll<int, scalar_t>(*problemComm, Teuchos::REDUCE_SUM, 
    numGlobalParts, partTotal.getRawPtr(), globalPartTotal.getRawPtr());

  bool emptyParts = false;
  for (int p=0; p < numGlobalParts; p++){
    if (globalPartTotal[p] == 0){
      emptyParts = true;
      break;
    }
  }

  scalar_t imbal = 0.0;
  if (!emptyParts){
    if (globalTotalWeight > 0){
      for (int i=0; i < numGlobalParts; i++){
        if (part_sizes[i] > 0){
          scalar_t tmp = 
            globalPartTotal[i] / (globalTotalWeight * part_sizes[]);
          if (tmp > imbal) imbal = tmp;
        }
      }
    }
    imbal = (imbal > 0 ? imbal : 1.0);
  }
  else{
    imbal= -1;  /* flag some part_sizes are zero */
  }

  Array<double> imbalance(weightDim);

  imbalance[0] = imbal;

  if (debug){
    if (imbalance[0] > Teuchos::as<scalar_t>(imbalanceTolerance)){
      ostringstream oss << "Warning: imbalance is " << imbal << std::endl;
      env->debugOut_->print(BASIC_STATUS, oss.str());
    }
    else{
      ostringstream oss << "Imbalance: " << imbal << std::endl;
      env->debugOut_->print(DETAILED_STATUS, oss.str());
    }
  }

  // Done, update the solution

  solution.setParts(idList, gnoPart, imbalance.view(0, weightDim));

  if (debug){
    ostringstream oss << "Exiting AlgBlock" << std::endl;
    env->debugOut_->print(DETAILED_STATUS, oss.str());
  }
}

}   // namespace Zoltan2

#endif
