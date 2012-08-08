// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef _ZOLTAN2_ALGBLOCK_HPP
#define _ZOLTAN2_ALGBLOCK_HPP_

#include <Zoltan2_IdentifierModel.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_GetParameter.hpp>

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
 *   \todo The metrics come out really bad.  Is it an error in
 *                algorithm or in metrics.
 */


template <typename Adapter>
void AlgBlock(
  const RCP<const Environment> &env,
  const RCP<Comm<int> > &problemComm,
  const RCP<const IdentifierModel<typename Adapter::base_adapter_t> > &ids, 
  RCP<PartitioningSolution<Adapter> > &solution
) 
{
  using std::string;
  using std::ostringstream;
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::scalar_t scalar_t;

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

  const Teuchos::ParameterList &pl = env->getParameters();

  getParameterValue<string>(pl, "partitioning",
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

  getParameterValue<double>(pl, "partitioning",
    "imbalance_tolerance", isSet, imbalanceTolerance);

  if (!isSet)
    imbalanceTolerance = 1.1;

  ////////////////////////////////////////////////////////
  // From the IdentifierModel we need:
  //    the number of gnos
  //    number of weights per gno
  //    the weights

  size_t numGnos = ids->getLocalNumIdentifiers();

  ArrayView<const gno_t> idList;
  typedef StridedData<lno_t, scalar_t> input_t;
  ArrayView<input_t> wgtList;
  
  ids->getIdentifierList(idList, wgtList);

  // If user supplied no weights, we use uniform weights.
  Array<input_t> uwArray(1);
  if (wgtList.size() == 0)
    wgtList = uwArray.view(0,1);

  bool uniformWeights = (wgtList[0].size() == 0);

  ////////////////////////////////////////////////////////
  // From the Solution we get part information.

  size_t numGlobalParts = solution->getTargetGlobalNumberOfParts();

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

  if (!uniformWeights){
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

  if (!solution->criteriaHasUniformPartSizes(0))
    for (unsigned int i=0; i<numGlobalParts; i++)
      part_sizes[i] = solution->getCriteriaPartSize(0, i);
  else
    for (unsigned int i=0; i<numGlobalParts; i++)
      part_sizes[i] = 1.0 / numGlobalParts;

  for (unsigned int i=1; i<numGlobalParts; i++)
    part_sizes[i] += part_sizes[i-1];

  // TODO assertion that last part sizes is about equal to 1.0

  if (env->getDebugLevel() >= VERBOSE_DETAILED_STATUS) {
    ostringstream oss("Part sizes: ");
    for (unsigned int i=0; i < numGlobalParts; i++)
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

  env->memory("Block algorithm memory");

  for (size_t i=0; i<numGnos; i++){
    scalar_t gnoWeight = (uniformWeights ? 1.0 : wgtList[0][i]);
    /* wtsum is now sum of all lower-ordered object */
    /* determine new partition number for this object,
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

  ArrayRCP<const gno_t> gnos = arcpFromArrayView(idList);

  solution->setParts(gnos, gnoPart);

  env->debug(DETAILED_STATUS, string("Exiting AlgBlock"));
}

}   // namespace Zoltan2

#endif
