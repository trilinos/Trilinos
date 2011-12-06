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
  const RCP<const Comm<int> > &problemComm,
  const RCP<const IdentifierModel<Adapter> > &ids, 
  ArrayView<size_t> &gnoPart,
  ArrayView<scalar_t> &imbalance
) 
{
  using std::string;
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::gno_t gno_t;
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

  // TODO - need a method to tell me local number of parts based
  //       on the value of those last two parameters.

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

  Z2_GLOBAL_BUG_ASSERTION(env, "partition array must be pre-allocated",
    gnoPart.size() >= numGnods, BASIC_ASSERTION);

  ArrayView<const gno_t> idList;
  ArrayView <const scalar_t> wgtList;
  
  ids->getIdentifierList(idList, wgtList);

  // Block partitioning algorithm lifted from zoltan/src/simple/block.c
  // The solution is:
  //    a list of part numbers in gno order
  //    an imbalance for each weight 

  int rank = env->myRank_;
  int nprocs = env->numProcs_;

  // TODO We are only using the first weight.

  scalar_t wtsum(0);

  if (wtflag){ /* Sum up first weight */
    for (size_t i=0; i<numGnos; i+=wtflag)
      wtsum += wgtList[i];
  }
  else
    wtsum = static_cast<scalar_t>(numGnos);

  Array<scalar_t> scansum(nprocs+1, 0);

  Teuchos::gatherAll<int, scalar_t>(*comm, 1, &wtsum, nprocs,
    scansum.getRawPtr()+1);

  /* scansum = sum of weights on lower processors, excluding self. */

  for (int i=2; i<=nprocs; i++)
    scansum[i] += scansum[i-1];

  scalar_t globalTotalWeight = scansum[nprocs];

  /* Overwrite part_sizes with cumulative sum (inclusive) part_sizes. */
  /* A cleaner way is to make a copy, but this works. */

  // TODO: we need an interface for part sizes - paramter list?
  Array<scalar_t> part_sizes(numGlobalParts, 1.0/numGlobalParts); 
  for (i=1; i<numGlobalParts; i++)
    part_sizes[i] += part_sizes[i-1];

  /* Loop over objects and assign partition. */
  size_t part = 0;
  wtsum = scansum[rank];
  Array<scalar_t> partTotal(numGlobalParts, 0);

  for (size_t i=0,j=0; i<numGnos; i++,j+=wtflag){
    /* wtsum is now sum of all lower-ordered object */
    /* determine new partition number for this object,
       using the "center of gravity" */
    while (part<numGlobalParts-1 && 
           (wtsum+0.5*(wtflag? wgtList[j]: 1.0))
           > part_sizes[part]*globalTotalWeight)
      part++;
    gnoPart[i] = part;
    partTotal[part] += wgtList[j];
    wtsum += (wtflag? wgtList[j] : 1.0);
  }

  // Compute the imbalance - taken from Zoltan's object_metrics function.
  //    TODO check this   copied with out thinking about it
  //   TODO we only get imbalance for one weight

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
  imbalance[0] = imbal;
}

}

}
#endif
