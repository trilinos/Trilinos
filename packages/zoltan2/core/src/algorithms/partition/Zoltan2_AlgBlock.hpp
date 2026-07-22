// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _ZOLTAN2_ALGBLOCK_HPP_
#define _ZOLTAN2_ALGBLOCK_HPP_

#include <Zoltan2_Algorithm.hpp>
#include <Zoltan2_IdentifierModel.hpp>
#include <Zoltan2_PartitioningSolution.hpp>

#include <bitset>
#include <sstream>
#include <string>

/*! \file Zoltan2_AlgBlock.hpp
 *  \brief The algorithm for block partitioning.
 */

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
 *  \param adapter    adapter used to create IdentifierModel
 *
 *  Preconditions: The parameters in the environment have been
 *    processed (committed).  No special requirements on the
 *    identifiers.
 *
 *   \todo Block partitioning uses one weight only
 *   \todo check for memory allocation failures
 *   \todo The metrics come out really bad.  Is it an error in
 *                algorithm or in metrics.
 */


template <typename Adapter>
class AlgBlock : public Algorithm<Adapter>
{

private:
  const RCP<const Environment> env;
  const RCP<const Comm<int> > problemComm;
  const RCP<const typename Adapter::base_adapter_t > adapter;

public:
  typedef typename Adapter::lno_t lno_t;       // local ids
  typedef typename Adapter::gno_t gno_t;       // global ids
  typedef typename Adapter::scalar_t scalar_t; // scalars
  typedef typename Adapter::part_t part_t;     // part numbers

  // Constructor
  AlgBlock(
    const RCP<const Environment> &env_,
    const RCP<const Comm<int> > &problemComm_,
    const RCP<const typename Adapter::base_adapter_t > &adapter_
  ) :
    env(env_), problemComm(problemComm_), adapter(adapter_)
  {}

  // Partitioning method
  void partition(const RCP<PartitioningSolution<Adapter> > &solution)
  {
    env->debug(DETAILED_STATUS, std::string("Entering AlgBlock"));

    int rank = env->myRank_;
    int nprocs = env->numProcs_;

    ////////////////////////////////////////////////////////
    // From the IdentifierModel we need:
    //    the number of gnos
    //    number of weights per gno
    //    the weights

    modelFlag_t modelFlag;
    IdentifierModel<typename Adapter::base_adapter_t> ids(adapter, env, problemComm, modelFlag);
    size_t numGnos = ids.getLocalNumIdentifiers();

    ArrayView<const gno_t> idList;
    typedef StridedData<lno_t, scalar_t> input_t;
    ArrayView<input_t> wgtList;

    ids.getIdentifierList(idList, wgtList);

    // If user supplied no weights, we use uniform weights.
    bool uniformWeights = (wgtList.size() == 0);

    ////////////////////////////////////////////////////////
    // Partitioning problem parameters of interest:
    //    objective
    //    imbalance_tolerance

    const Teuchos::ParameterList &pl = env->getParameters();
    const Teuchos::ParameterEntry *pe;

    pe = pl.getEntryPtr("partitioning_objective");
    if (pe) {
      std::string po = pe->getValue<std::string>(&po);
      if (po == std::string("balance_object_count"))
        uniformWeights = true;    // User requests that we ignore weights
    }

    double imbalanceTolerance=1.1;
    pe = pl.getEntryPtr("imbalance_tolerance");
    if (pe) imbalanceTolerance = pe->getValue<double>(&imbalanceTolerance);

    ////////////////////////////////////////////////////////
    // From the Solution we get part information:
    // number of parts and part sizes

    size_t numGlobalParts = solution->getTargetGlobalNumberOfParts();

    Array<scalar_t> part_sizes(numGlobalParts);

    if (solution->criteriaHasUniformPartSizes(0))
      for (unsigned int i=0; i<numGlobalParts; i++)
        part_sizes[i] = 1.0 / numGlobalParts;
    else
      for (unsigned int i=0; i<numGlobalParts; i++)
        part_sizes[i] = solution->getCriteriaPartSize(0, i);

    for (unsigned int i=1; i<numGlobalParts; i++)
      part_sizes[i] += part_sizes[i-1];

    // TODO assertion that last part sizes is about equal to 1.0


    ////////////////////////////////////////////////////////
    // The algorithm
    //
    // Block partitioning algorithm lifted from zoltan/src/simple/block.c
    // The solution is:
    //    a list of part numbers in gno order
    //    an imbalance for each weight

    scalar_t wtsum(0);

    if (!uniformWeights) {
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

    if (env->getDebugLevel() >= VERBOSE_DETAILED_STATUS) {
      std::ostringstream oss("Part sizes: ");
      for (unsigned int i=0; i < numGlobalParts; i++)
        oss << part_sizes[i] << " ";
      oss << std::endl << std::endl << "Weights : ";
      for (int i=0; i <= nprocs; i++)
        oss << scansum[i] << " ";
      oss << std::endl;
      env->debug(VERBOSE_DETAILED_STATUS, oss.str());
    }

    /* Loop over objects and assign part. */
    part_t part = 0;
    wtsum = scansum[rank];
    Array<scalar_t> partTotal(numGlobalParts, 0);
    ArrayRCP<part_t> gnoPart= arcp(new part_t[numGnos], 0, numGnos);

    env->memory("Block algorithm memory");

    for (size_t i=0; i<numGnos; i++){
      scalar_t gnoWeight = (uniformWeights ? 1.0 : wgtList[0][i]);
      /* wtsum is now sum of all lower-ordered object */
      /* determine new part number for this object,
         using the "center of gravity" */
      while (unsigned(part)<numGlobalParts-1 &&
             (wtsum+0.5*gnoWeight) > part_sizes[part]*globalTotalWeight)
        part++;
      gnoPart[i] = part;
      partTotal[part] += gnoWeight;
      wtsum += gnoWeight;
    }

    ////////////////////////////////////////////////////////////
    // Done

    solution->setParts(gnoPart);

    env->debug(DETAILED_STATUS, std::string("Exiting AlgBlock"));
  }
};

}   // namespace Zoltan2

#endif
