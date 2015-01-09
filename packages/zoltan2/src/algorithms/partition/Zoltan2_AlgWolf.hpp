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
#ifndef _ZOLTAN2_ALGWOLF_HPP_
#define _ZOLTAN2_ALGWOLF_HPP_

#include <Zoltan2_IdentifierModel.hpp>
#include <Zoltan2_PartitioningSolution.hpp>
#include <Zoltan2_Algorithm.hpp>
#include <Zoltan2_AlgRCB.hpp>

#include <sstream>
#include <string>
#include <bitset>

/*! \file Zoltan2_AlgWolf.hpp
 *  \brief The algorithm for Wolf partitioning.
 */

namespace Zoltan2
{

// /*! \brief The boolean parameters of interest to the Block algorithm.
//  */
// enum blockParams{
//   block_balanceCount,            /*!< objective = balance_object_count */
//   block_balanceWeight,          /*!< objective = balance_object_weight */
//   block_minTotalWeight,      /*!< objective = mc_minimize_total_weight */
//   block_minMaximumWeight,  /*!< objective = mc_minimize_maximum_weight */
//   block_balanceTotalMaximum, /*!< objective = mc_balance_total_maximum */
//   NUM_BLOCK_PARAMS
// };

////////////////////////////////////////////////////////////////////////////////
/*! Wolf partitioning method.
 *
 *  \param env   library configuration and problem parameters
 *  \param problemComm  the communicator for the problem
 *  \param ids    an Identifier model
 *
 *  Preconditions: The parameters in the environment have been
 *    processed (committed).  No special requirements on the
 *    identifiers.
 *
 */
////////////////////////////////////////////////////////////////////////////////
template <typename Adapter>
class AlgWolf : public Algorithm<Adapter>
{

private:
  const RCP<const Environment> mEnv;
  const RCP<Comm<int> > mProblemComm;

  const RCP<const GraphModel<typename Adapter::base_adapter_t> > &mGraphModel;
  const RCP<const CoordinateModel<typename Adapter::base_adapter_t> > mIds;

  RCP<PartitioningSolution<Adapter> > mSolution; //Not sure if this should be saved

public:
  // Constructor
  AlgWolf(const RCP<const Environment> &env_,
	  const RCP<Comm<int> > &problemComm_,
	  const RCP<const GraphModel<typename Adapter::base_adapter_t> > &gModel_,
	  const RCP<const CoordinateModel<typename Adapter::base_adapter_t> > &cModel_)
    :mEnv(env_), mProblemComm(problemComm_), mGraphModel(gModel_), mIds(cModel_)
  {
#ifndef INCLUDE_ZOLTAN2_EXPERIMENTAL
    Z2_THROW_EXPERIMENTAL("Zoltan2 Wolf is strictly experimental software ")
#endif

#ifndef INCLUDE_ZOLTAN2_EXPERIMENTAL_WOLF
    Z2_THROW_EXPERIMENTAL_WOLF("Zoltan2 Wolf is strictly experimental software ")
#endif

  }

  // Partitioning method
  void partition(const RCP<PartitioningSolution<Adapter> > &solution_);

};
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
template <typename Adapter>
void AlgWolf<Adapter>::partition(
  const RCP<PartitioningSolution<Adapter> > &solution_
)
{
    // using std::string;
    // using std::ostringstream;

    // typedef typename Adapter::lno_t lno_t;     // local ids
    // typedef typename Adapter::gno_t gno_t;     // global ids
    // typedef typename Adapter::scalar_t scalar_t;   // scalars

    mEnv->debug(DETAILED_STATUS, std::string("Entering AlgWolf"));

    // int rank = env->myRank_;
    // int nprocs = env->numProcs_;

    // ////////////////////////////////////////////////////////
    // // From the CoordinateModel we need:
    // //    the number of gnos
    // //    number of weights per gno
    // //    the weights

    // size_t numGnos = ids->getLocalNumIdentifiers();

    // ArrayView<const gno_t> idList;
    // typedef StridedData<lno_t, scalar_t> input_t;
    // ArrayView<input_t> wgtList;
  
    // ids->getIdentifierList(idList, wgtList);

    // // If user supplied no weights, we use uniform weights.
    // bool uniformWeights = (wgtList.size() == 0);


   

    // First, let's partition with RCB

    // Q: can I use solution passed into alg or do I need to create a different one?
    //    For now using the one passed into alg

    {
    AlgRCB<Adapter> algrcb(this->mEnv, mProblemComm, this->mIds);
    algrcb.partition(solution_);
    }


    // ////////////////////////////////////////////////////////
    // // Partitioning problem parameters of interest:
    // //    objective
    // //    imbalance_tolerance

    // const Teuchos::ParameterList &pl = env->getParameters();
    // const Teuchos::ParameterEntry *pe;

    // pe = pl.getEntryPtr("partitioning_objective");
    // if (pe) {
    //   string po = pe->getValue<string>(&po);
    //   if (po == string("balance_object_count"))
    //     uniformWeights = true;    // User requests that we ignore weights
    // }

    // double imbalanceTolerance=1.1;
    // pe = pl.getEntryPtr("imbalance_tolerance");
    // if (pe) imbalanceTolerance = pe->getValue<double>(&imbalanceTolerance);

    // ////////////////////////////////////////////////////////
    // // From the Solution we get part information:
    // // number of parts and part sizes

    // size_t numGlobalParts = solution->getTargetGlobalNumberOfParts();

    // Array<scalar_t> part_sizes(numGlobalParts);

    // if (solution->criteriaHasUniformPartSizes(0))
    //   for (unsigned int i=0; i<numGlobalParts; i++)
    //     part_sizes[i] = 1.0 / numGlobalParts;
    // else
    //   for (unsigned int i=0; i<numGlobalParts; i++)
    //     part_sizes[i] = solution->getCriteriaPartSize(0, i);

    // for (unsigned int i=1; i<numGlobalParts; i++)
    //   part_sizes[i] += part_sizes[i-1];

    // // TODO assertion that last part sizes is about equal to 1.0


    // ////////////////////////////////////////////////////////
    // // The algorithm
    // //
    // // Wolf partitioning algorithm lifted from zoltan/src/simple/block.c
    // // The solution is:
    // //    a list of part numbers in gno order
    // //    an imbalance for each weight 

    // scalar_t wtsum(0);

    // if (!uniformWeights) {
    //   for (size_t i=0; i<numGnos; i++)
    //     wtsum += wgtList[0][i];          // [] operator knows stride
    // }
    // else
    //   wtsum = static_cast<scalar_t>(numGnos);

    // Array<scalar_t> scansum(nprocs+1, 0);

    // Teuchos::gatherAll<int, scalar_t>(*problemComm, 1, &wtsum, nprocs,
    //   scansum.getRawPtr()+1);

    // /* scansum = sum of weights on lower processors, excluding self. */

    // for (int i=2; i<=nprocs; i++)
    //   scansum[i] += scansum[i-1];

    // scalar_t globalTotalWeight = scansum[nprocs];

    // if (env->getDebugLevel() >= VERBOSE_DETAILED_STATUS) {
    //   ostringstream oss("Part sizes: ");
    //   for (unsigned int i=0; i < numGlobalParts; i++)
    //     oss << part_sizes[i] << " ";
    //   oss << std::endl << std::endl << "Weights : ";
    //   for (int i=0; i <= nprocs; i++)
    //     oss << scansum[i] << " ";
    //   oss << std::endl;
    //   env->debug(VERBOSE_DETAILED_STATUS, oss.str());
    // }

    // /* Loop over objects and assign part. */
    // partId_t part = 0;
    // wtsum = scansum[rank];
    // Array<scalar_t> partTotal(numGlobalParts, 0);
    // ArrayRCP<partId_t> gnoPart= arcp(new partId_t [numGnos], 0, numGnos);

    // env->memory("Wolf algorithm memory");

    // for (size_t i=0; i<numGnos; i++){
    //   scalar_t gnoWeight = (uniformWeights ? 1.0 : wgtList[0][i]);
    //   /* wtsum is now sum of all lower-ordered object */
    //   /* determine new part number for this object,
    //      using the "center of gravity" */
    //   while (unsigned(part)<numGlobalParts-1 && 
    //          (wtsum+0.5*gnoWeight) > part_sizes[part]*globalTotalWeight)
    //     part++;
    //   gnoPart[i] = part;
    //   partTotal[part] += gnoWeight;
    //   wtsum += gnoWeight;
    // }

    // ////////////////////////////////////////////////////////////
    // // Done

    // ArrayRCP<const gno_t> gnos = arcpFromArrayView(idList);
    // solution->setParts(gnos, gnoPart, true);

    mEnv->debug(DETAILED_STATUS, std::string("Exiting AlgWolf"));
}






}   // namespace Zoltan2

#endif
