#ifndef _ZOLTAN2_ALGBLOCK_HPP_
#define _ZOLTAN2_ALGBLOCK_HPP_

#include <Zoltan2_Standards.hpp>
#include <Zoltan2_Environment.hpp>
#include <Zoltan2_IdentifierModel.hpp>
#include <Zoltan2_PartitioningSolution.hpp>

/*! \file Zoltan2_AlgBlock.hpp
 *  \brief The algorithm for block partitioning.
 */

namespace Zoltan2{

/*! Block partitioning method.
 *
 *  \param env  parameters for the problem and library configuration
 *  \param model the application Ids and weights
 *  \param gnoPart on return, the resulting assignment of gno to part number
 *  \param imbalance on return, the imbalance achieved by this partitioning
 *
 *  Preconditions: The parameters in the environment have been
 *    processed (committed).  The model supplies consecutive identifiers
 *    increasing with process rank and beginning at some base.
 */


template <typename Adapter>
void AlgBlock(
  const RCP<const Environment> &env,
  const RCP<const IdentifierModel<Adapter> > &ids, 
  ArrayView<int> &gnoPart,
  ArrayView<double> &imbalance
) 
{
  using Teuchos::ParameterEntry;
  using std::string;
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::gid_t gid_t;
  typedef typename Adapter::lid_t lid_t;
  typedef typename Adapter::scalar_t scalar_t;

  Z2_GLOBAL_BUG_ASSERTION(env, "parameters are not committed",
    env.parametersAreCommitted(), DEBUG_MODE_ASSERTION);

  // Library parameters of interest:
  //    are we printing out debug messages
  //    are we timing
  //    are we computing memory used

  bool debug = env.doStatus();
  bool timing = env.doTiming();
  bool memstats = env.doMemoryProfiling();

  // Problem parameters of interest:
  //    objective
  //    imbalance_tolerance
  //    num_global_parts
  //    num_local_parts

  const Teuchos::ParameterList &pl = env.getParams();

  string &objective = pl.get(string("objective"));
  double &imbalanceTolerance = pl.get(string("imbalance_tolerance"));
  int &numGlobalParts = pl.get(string("num_global_parts"));
  int &numLocalParts = pl.get(string("num_local_parts"));

  // Parameters that may drive algorithm choices
  //    speed_versus_quality
  //    memory_footprint_versus_runtime

  string &svq = pl.get(string("speed_versus_quality"));
  bool fastSolution = (svq==string("speed"));
  bool goodSolution = (svq==string("quality"));
  bool balancedSolution = (svq==string("balance"));
  
  string &mvr = pl.get(string("memory_footprint_versus_runtime"));
  bool lowMemory = (mvr==string("memory") || (mvr==string("memory_footprint"));
  bool lowRunTime = (mvr==string("runtime"));
  bool balanceMemoryRunTime = (svq==string("balance"));

  // From the IdentifierModel we need:
  //    the number of gnos
  //    number of weights per gno
  //    the weights

  size_t numGnos = getLocalNumIdentifiers();
  int wtflag = ids->getIdentifierWeightDim();

  ArrayView<const gno_t> idList;
  ArrayView <const scalar_t> wgtList;
  
  ids->getIdentifierList(idList, wgtList);

  // Block partitioning algorithm lifted from zoltan/src/simple/block.c
  // The solution is:
  //    a list of part numbers in gno order
  //    an imbalance for each weight 

  int rank = env->myRank_;
  int nprocs = env->numProcs_;
  Array<float> part_sizes(numGlobalParts, 1.0);
  Array<int> newparts(numGnos);

  int i, part;
  scalar_t wtsum;

  Array<scalar_t> scansum(nprocs+1);

  if (wtflag){ /* Sum up local object weights. */
    wtsum = 0.0;
    for (i=0; i<numGnos; i++)
      wtsum += wgtList[i];
  }
  else
    wtsum = numGnos;



  /* Cumulative global wtsum FIXME */
  MPI_Allgather(&wtsum, 1, MPI_DOUBLE, &scansum[1], 1, MPI_DOUBLE,
                zz->Communicator);
  /* scansum = sum of weights on lower processors, excluding self. */
  scansum[0] = 0.;
  for (i=1; i<=nprocs; i++)
    scansum[i] += scansum[i-1];

  /* Overwrite part_sizes with cumulative sum (inclusive) part_sizes. */
  /* A cleaner way is to make a copy, but this works. */
  for (i=1; i<numGlobalParts; i++)
    part_sizes[i] += part_sizes[i-1];

  /* Loop over objects and assign partition. */
  part = 0;
  wtsum = scansum[rank];
  for (i=0; i<numGnos; i++){
    /* wtsum is now sum of all lower-ordered object */
    /* determine new partition number for this object,
       using the "center of gravity" */
    while (part<numGlobalParts-1 && 
           (wtsum+0.5*(wtflag? wgtList[i]: 1.0))
           > part_sizes[part]*scansum[nprocs])
      part++;
    gnoPart[i] = part;
    wtsum += (wtflag? wgtList[i] : 1.0);
  }
}

}

}
#endif
