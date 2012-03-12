// @HEADER
// ***********************************************************************
//                Copyright message goes here. 
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_AlgRCB.hpp
    \brief Contains the recursive coordinate bisection algorthm.
*/

#ifndef _ZOLTAN2_ALGRCB_HPP_
#define _ZOLTAN2_ALGRCB_HPP_

#include <Zoltan2_CoordinateModel.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_Metric.hpp>

#include <Teuchos_ParameterList.hpp>

#include <sstream>
#include <string>

namespace Zoltan2{

// The boolean parameters
enum rcbParams{
  doStatus,
  doTiming,
  doMemstats,
  fastSolution,
  goodSolution,
  balancedSolution,
  lowMemory,
  lowRunTime,
  balanceMemoryRunTime,
  balanceCount,
  balanceWeight,
  minTotalWeight,
  minMaximumWeight,
  balanceTotalMaximum,
  averageCuts,
  rectilinearBlocks,
  NUM_RCB_PARAMS
};

/*! \brief Recursive coordinate bisection partitioning.
 *
 *  \param env   library configuration and problem parameters
 *  \param problemComm  the communicator for the problem
 *  \param coords    a CoordinateModel with user data
 *  \param solution  a PartitioningSolution, on input it 
 *      contains part information, on return it also contains 
 *      the solution and quality metrics.
 *                    
 *   \todo timing and memory usage profiling
 *   \todo numParts != numProcs (requires work in PartitioningSolution)
 *   \todo  catch errors and pass back
 *   \todo write the rcb tree back to the solution
 *   \todo for "repartition", start with the tree in the solution
 */

template <typename Adapter>
void AlgRCB(
  const RCB<const Environment> &env,
  const RCB<const Comm<int> > &problemComm,
  const RCB<const CoordinateModel<Adapter> > &coords, 
  RCB<PartitioningSolution<typename Adapter::user_t> > &solution
) 
{
  using std::string;
  using std::ostringstream;
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::scalar_t scalar_t;

  int rank = env->myRank_;
  int nprocs = env->numProcs_;

  std::bitset<NUM_RCB_PARAMS> params;

  ////////////////////////////////////////////////////////
  // Library parameters of interest:
  //
  //    are we printing out debug messages
  //    are we timing
  //    are we computing memory used
  //    speed_versus_quality
  //    memory_versus_speed

  params.set(doStatus, env->doStatus());
  params.set(doTiming, env->doTiming());
  params.set(doMemstats, env->doMemoryProfiling());

  if (params.test(doStatus))
    env->debug(DETAILED_STATUS, string("Entering AlgRCB"));

  const Teuchos::ParameterList &pl = env->getParameters();
  const string defaultVal("balance");

  const string *mvr = pl.getPtr<string>("memory_versus_speed");
  if (!mvr)
    mvr = &defaultVal;
  
  const string *svq = pl.getPtr<string>("speed_versus_quality");
  if (!svq)
    svq = &defaultVal;

  params.set(fastSolution, (*svq==string("speed")));
  params.set(goodSolution, (*svq==string("quality")));
  params.set(balancedSolution, (*svq==string("balance")));
 
  params.set(lowMemory, (*mvr==string("memory")));
  params.set(lowRunTime, (*mvr==string("speed")));
  params.set(balancedMemoryRunTime, (*mvr==string("balance")));

  ////////////////////////////////////////////////////////
  // Partitioning problem parameters of interest:
  //    objective
  //    imbalance_tolerance

  const string *obj=NULL;
  const double *tol=NULL;

  if (env->hasPartitioningParameters()){
    const Teuchos::ParameterList &plPart = pl.sublist("partitioning");

    obj = plPart.getPtr<string>("objective");
    tol = plPart.getPtr<double>("imbalance_tolerance");
  }

  double imbalanceTolerance = (tol ? *tol : 1.1);
  string objective = (obj ? *obj : string("balance_object_weight"));

  params.set(balanceCount, (objective == string("balance_object_count")));
  params.set(balanceWeight, (objective == string("balance_object_weight")));
  params.set(minTotalWeight, 
    (objective == string("multicriteria_minimize_total_weight")));
  params.set(minMaximumWeight, 
    (objective == string("multicriteria_minimize_maximum_weight")));
  params.set(balanceTotalMaximum, 
    (objective == string("multicriteria_balance_total_maximum")));

  ////////////////////////////////////////////////////////
  // Geometric partitioning problem parameters of interest:
  //    average_cuts
  //    rectilinear_blocks
  //    bisection_num_test_cuts (experimental)

  params.set(averageCuts,false);
  params.set(rectilinearBlocks,false);
  int  numTestCuts=3; 

  if (env->hasPartitioningParameters()){
    const Teuchos::ParameterList &plPart = pl.sublist("partitioning");
    if (env->hasSublist(plPart, std::string("geometry"))){
      const Teuchos::ParameterList &geom = plPart.sublist("geometry");

      const int *zeroOne = geom.getPtr<int>("average_cuts");
      if (zeroOne && (*zeroOne==1)) params.set(averageCuts, true);

      const int *zeroOne = geom.getPtr<int>("rectilinear_blocks");
      if (zeroOne && (*zeroOne==1)) params.set(rectilinearBlocks, true);

      const int *numCuts = geom.getPtr<int>("bisection_num_test_cuts");
      if (numCuts) numTestCuts = *numCuts;
    }
  }

  ////////////////////////////////////////////////////////
  // From the CoordinateModel we need:
  //    coordinate values
  //    coordinate weights, if any
  //    coordinate global Ids

  typedef StridedInput<lno_t, scalar_t> input_t;

  int coordDim = coords->getCoordinateDim();
  int weightDim = coords->getCoordinateWeightDim();
  size_t numLocalCoords = coords->getLocalNumCoordinates();
  global_size_t numGlobalCoords = coords->getGlobalNumCoordinates();

  ArrayView<const gno_t> gnos;
  ArrayView<input_t>     xyz;
  ArrayView<input_t>     wgts;

  coords->getCoordinates(gnos, xyz, wgts);

  ////////////////////////////////////////////////////////
  // From the Solution we get part information.
  //

  size_t numGlobalParts = solution->getGlobalNumberOfParts();
  Array<bool> uniformParts(weightDim);
  for (int wdim = 0; wdim < weightDim; wdim++)
    uniformParts[wdim] = solution->criteriaHasUniformPartSizes(wdim);
  
  ////////////////////////////////////////////////////////
  // The algorithm
  //
  // A solution is:
  //    a list of part numbers in gno order
  //    an imbalance for each weight 
  //    a tree representing the cuts

}

template <typename scalar_t, typename lno_t, typename gno_t, typename node_t>
  static int findCut(
    const RCP<const Teuchos::Comm<int> &comm,
    Tpetra::Vector<scalar_t, lno_t, gno_t, node_t> &coords,
    Array<Tpetra::Vector<scalar_t> > &weights, 
    scalar_t coordGlobalMin,
    scalar_t coordGlobalMax,
    int numTestCuts,
    float imbalanceTolerance,
    std::bitset<NUM_RCB_PARAMS> &params)
{
  int rank = comm->getRank();
  int nprocs = comm->getSize();

  // For each weight dimension, are weights uniform or not?

  int weightDim = weights.size();
  Array<int> uniform(weightDim, 0);
  for (int i=0; i < weightDim; i++)
    if (weights[i].getGlobalLength() == 0)
      uniform[i] = 1;

  // Where do we make the first test cuts?

  int numRegions = numTestCuts + 1;
  Array<scalar_t> testCuts(numTestCuts);
  scalar_t diff = (coordGlobalMax - coordGlobalMin) / numRegions;
  for (int i=0; i < numTestCuts; i++){
    testCuts[i] = coordGlobalMin + diff;
    diff += diff;
  }

  // Sum weights or counts for each weight dimension

  Array<Array<gno_t> > counts[weightDim];
  Array<Array<scalar_t> > weightSums[weightDim];

  if (params.test(balanceCount)){
    // We ignore the weights and balance the number of objects
  }
  else if params.test(balanceWeight)){
    // If weightDim is one, balance the weights.  If not, it's minTotalWeight.
  }
  else if params.test(minTotalWeight)){
    // Balance the total of the weightDim weights
    // Zoltan RCB_MULTICRITERIA_NORM 1
  }
  else if params.test(minMaximumWeight)){
    // Check all weightDim weights, and minimize the maximum imbalance.
    // Zoltan RCB_MULTICRITERIA_NORM 3
  }
  else if params.test(balanceTotalMaximum)){
    // Zoltan RCB_MULTICRITERIA_NORM 2
  }
  else{
    // A logic error this should never happen
  }
  


}

}   // namespace Zoltan2
#endif
