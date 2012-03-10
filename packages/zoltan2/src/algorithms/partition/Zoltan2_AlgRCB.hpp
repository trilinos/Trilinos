// @HEADER
// ***********************************************************************
//                Copyright message goes here. 
// ***********************************************************************
// @HEADER

/*! \file Zoltan2_AlgRCB.hpp
    \brief Contains the recursive coordinate bisection algorthm.
*/

#ifndef _ZOLTAN2_ALGRCP_HPP_
#define _ZOLTAN2_ALGRCP_HPP_

#include <Zoltan2_CoordinateModel.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_Metric.hpp>

#include <Teuchos_ParameterList.hpp>

#include <sstream>
#include <string>

namespace Zoltan2{

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
void AlgRCP(
  const RCP<const Environment> &env,
  const RCP<const Comm<int> > &problemComm,
  const RCP<const CoordinateModel<Adapter> > &coords, 
  RCP<PartitioningSolution<typename Adapter::user_t> > &solution
) 
{
  using std::string;
  using std::ostringstream;
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::scalar_t scalar_t;

  int rank = env->myRank_;
  int nprocs = env->numProcs_;

  ////////////////////////////////////////////////////////
  // Library parameters of interest:
  //
  //    are we printing out debug messages
  //    are we timing
  //    are we computing memory used
  //    speed_versus_quality
  //    memory_versus_speed

  bool debug = env->doStatus();
  bool timing = env->doTiming();
  bool memstats = env->doMemoryProfiling();

  if (debug)
    env->debug(DETAILED_STATUS, string("Entering AlgRCB"));

  const Teuchos::ParameterList &pl = env->getParameters();
  const string defaultVal("balance");

  const string *mvr = pl.getPtr<string>("memory_versus_speed");
  if (!mvr)
    mvr = &defaultVal;
  
  const string *svq = pl.getPtr<string>("speed_versus_quality");
  if (!svq)
    svq = &defaultVal;

  bool fastSolution = (*svq==string("speed"));
  bool goodSolution = (*svq==string("quality"));
  bool balancedSolution = (*svq==string("balance"));
 
  bool lowMemory = (*mvr==string("memory"));
  bool lowRunTime = (*mvr==string("speed"));
  bool balanceMemoryRunTime = (*mvr==string("balance"));

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

  bool balanceCount = (objective == string("balance_object_count"));
  bool balanceWeight = (objective == string("balance_object_weight"));
  bool minTotalWeight = 
    (objective == string("multicriteria_minimize_total_weight"));
  bool minMaximumWeight = 
    (objective == string("multicriteria_minimize_maximum_weight"));
  bool balanceTotalMaximum = 
    (objective == string("multicriteria_balance_total_maximum"));

  ////////////////////////////////////////////////////////
  // Geometric partitioning problem parameters of interest:
  //    average_cuts
  //    rectilinear_blocks
  //    bisection_num_test_cuts (experimental)

  bool averageCuts=false;
  bool rectilinearBlocks=false;
  int  numTestCuts=3; 

  if (env->hasPartitioningParameters()){
    const Teuchos::ParameterList &plPart = pl.sublist("partitioning");
    if (env->hasSublist(plPart, std::string("geometry"))){
      const Teuchos::ParameterList &geom = plPart.sublist("geometry");

      const int *zeroOne = geom.getPtr<int>("average_cuts");
      if (zeroOne && (*zeroOne==1)) averageCuts = true;

      const int *zeroOne = geom.getPtr<int>("rectilinear_blocks");
      if (zeroOne && (*zeroOne==1)) rectilinearBlocks = true;

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
  size_t numLocalParts = solution->getLocalNumberOfParts();

  const int *partDist = solution->getPartDistribution();
  const size_t *procDist = solution->getProcsParts();
  double *partSizes0 = getCriteriaPartSizes(0);
  
  ////////////////////////////////////////////////////////
  // The algorithm
  //
  // A solution is:
  //    a list of part numbers in gno order
  //    an imbalance for each weight 
  //    a tree representing the cuts

}   // namespace Zoltan2

#endif
