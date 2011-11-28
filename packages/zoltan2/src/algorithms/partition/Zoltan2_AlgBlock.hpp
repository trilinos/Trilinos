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
 *  \param gnoPart the resulting assignment of gno to part number
 *  \param imbalance the imbalance achieved by this partitioning
 *
 *  Preconditions: The parameters in the environment have been
 *    processed (committed).  The model supplies consecutive identifiers
 *    increasing with process rank and beginning at some base.
 */


template <typename Adapter>
void AlgBlock(
  const RCP<Environment> &env,
  const RCP<IdentifierModel<Adapter> > &model, 
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
  //    the gnos
  //    number of weights per gno
  //    the weights

  // The solution is:
  //    a list of part numbers in gno order
  //    an imbalance for each weight 
}

}
#endif
