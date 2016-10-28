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

/*! \file Zoltan2_OrderingMetricsUtility.hpp
 */

#ifndef ZOLTAN2_ORDERINGMETRICSUTILITY_HPP
#define ZOLTAN2_ORDERINGMETRICSUTILITY_HPP

#include <Zoltan2_OrderingMetrics.hpp>
#include <Zoltan2_MetricUtility.hpp>

namespace Zoltan2{

/*! \brief Compute ordering metrics.
 *
 *   \param env   The problem environment.
 *   \param comm  The problem communicator.
 *   \param ia the InputAdapter object which corresponds to the Solution.
 *   \param localSoln the local LocalOrderingSolution to be evaluated.
 *   \param metrics on return points to a list of named MetricValues objects
 */

template <typename Adapter>
  void localOrderingMetrics(
    const RCP<const Environment> &env,
    const RCP<const Comm<int> > &comm,
    const Adapter *ia,
    const LocalOrderingSolution<typename Adapter::lno_t> *localSoln,
    ArrayRCP<RCP<BaseClassMetrics<typename Adapter::scalar_t> > > &metrics)
{
  env->debug(DETAILED_STATUS, "Entering orderingMetrics"); // begin

  // set up some typedefs
  typedef typename Adapter::scalar_t scalar_t;
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::base_adapter_t base_adapter_t;
  typedef StridedData<lno_t, scalar_t> input_t;

  // get graph
  std::bitset<NUM_MODEL_FLAGS> modelFlags;
	RCP<GraphModel<base_adapter_t> > graph;
	const RCP<const base_adapter_t> bia =
	  rcp(dynamic_cast<const base_adapter_t *>(ia), false);
	graph = rcp(new GraphModel<base_adapter_t>(bia,env,comm,modelFlags));
  ArrayView<const gno_t> Ids;
  ArrayView<input_t> vwgts;
  ArrayView<const gno_t> edgeIds;
  ArrayView<const lno_t> offsets;
  ArrayView<input_t> wgts;
  ArrayView<input_t> vtx;
  graph->getEdgeList(edgeIds, offsets, wgts);
  lno_t numVertex = graph->getVertexList(Ids, vwgts);

  lno_t * perm = localSoln->getPermutationView();

  // print as matrix - this was debugging code which can be deleted later
  #define MDM
  #ifdef MDM
  for( int checkRank = 0; checkRank < comm->getSize(); ++checkRank ) {
    comm->barrier();
    if( checkRank == comm->getRank() ) {
      std::cout << "-----------------------------------------" << std::endl;
      std::cout << "Inspect rank: " << checkRank << std::endl;
      std::cout << std::endl;
      if(numVertex < 30) { // don't spam if it's too many...
        Array<lno_t> oldMatrix(numVertex*numVertex);
        Array<lno_t> newMatrix(numVertex*numVertex);

        // print the solution permutation
        std::cout << std::endl << "perm:  ";
        for(lno_t n = 0; n < numVertex; ++n) {
          std::cout << " " << perm[n] << " ";
        }

        lno_t * iperm = localSoln->getPermutationView(true);
        std::cout << std::endl << "iperm: ";
        for(lno_t n = 0; n < numVertex; ++n) {
          std::cout << " " << iperm[n] << " ";
        }
        std::cout << std::endl;
        // write 1's to old matrix (original form) and new matrix (using solution)
        for (lno_t y = 0; y < numVertex; y++) {
          for (lno_t n = offsets[y]; n < offsets[y+1]; ++n) {
            lno_t x = static_cast<lno_t>(edgeIds[n]); // to resolve
            if (x < numVertex && y < numVertex) { // to develop - for MPI this may not be local
              oldMatrix[x + y*numVertex] = 1;
              newMatrix[perm[x] + perm[y]*numVertex] = 1;
            }
          }
        }

        // print oldMatrix
        std::cout << std::endl << "unsolved graph in matrix form:" << std::endl;
        for(lno_t y = 0; y < numVertex; ++y) {
          for(lno_t x = 0; x < numVertex; ++x) {
            std::cout << " " << oldMatrix[x + y*numVertex];
          }
          std::cout << std::endl;
        }

        // print newMatrix
        std::cout << std::endl << "solved graph in matrix form:" << std::endl;
        for(lno_t y = 0; y < numVertex; ++y) {
          for(lno_t x = 0; x < numVertex; ++x) {
            std::cout << " " << newMatrix[x + y*numVertex];
          }
          std::cout << std::endl;
        }
        std::cout << std::endl;
      }
    }

    comm->barrier();
  }
  #endif // Ends temporary logging which can be deleted later

  // calculate bandwidth and envelope for unsolved and solved case
  lno_t bw_right_unsolved = 0;
  lno_t bw_left_unsolved = 0;
  lno_t bw_right_solved = 0;
  lno_t bw_left_solved = 0;

  lno_t envelope_unsolved = 0;
  lno_t envelope_solved = 0;

  for (lno_t j = 0; j < numVertex; j++) {
    lno_t y = Ids[j];
    for (auto n = offsets[j]; n < offsets[j+1]; ++n) {
      lno_t x = static_cast<lno_t>(edgeIds[n]); // to resolve
      if(x < numVertex) {
        lno_t x2 = perm[x];
        lno_t y2 = perm[y];

        // unsolved bandwidth calculation
        lno_t delta_right_unsolved = y - x;
        if (delta_right_unsolved > bw_right_unsolved) {
          bw_right_unsolved = delta_right_unsolved;
        }
        lno_t delta_left_unsolved = x - y;
        if (delta_left_unsolved > bw_left_unsolved) {
          bw_left_unsolved = delta_left_unsolved;
        }
        // solved bandwidth calculation
        lno_t delta_right_solved = y2 - x2;
        if (delta_right_solved > bw_right_solved) {
          bw_right_solved = delta_right_solved;
        }
        lno_t delta_left_solved = y2 - x2;
        if (delta_left_solved > bw_left_solved) {
          bw_left_solved = delta_left_solved;
        }

        // unsolved envelope calculation
        if(delta_right_unsolved > 0) {
          envelope_unsolved += delta_right_unsolved;
        }
        if(delta_left_unsolved > 0) {
          envelope_unsolved += delta_left_unsolved;
        }
        envelope_unsolved += 1; // need to check this

        // solved envelope calculation
        if(delta_right_solved > 0) {
          envelope_solved += delta_right_solved;
        }
        if(delta_left_solved > 0) {
          envelope_solved += delta_left_solved;
        }
        envelope_solved += 1; // need to check this
      }
    }
  }

  lno_t bw_solved = (bw_left_solved + bw_right_solved + 1);
  lno_t bw_unsolved = (bw_left_unsolved + bw_right_unsolved + 1);

  // add the new metrics
  typedef OrderingMetrics<scalar_t> om_t;

  // set up the bandwidth metric
  RCP<om_t> bandWidthMetric = addNewMetric<om_t, scalar_t>(env, metrics);
  bandWidthMetric->setName("bandwidth");
  bandWidthMetric->setMetricValue("unsolved", bw_unsolved);
  bandWidthMetric->setMetricValue("solved", bw_solved);

  // add the envelope metric
  RCP<om_t> envelopeMetric = addNewMetric<om_t, scalar_t>(env, metrics);
  envelopeMetric->setName("envelope");
  envelopeMetric->setMetricValue("unsolved", envelope_unsolved);
  envelopeMetric->setMetricValue("solved", envelope_solved);

  // add the separator size metric
  RCP<om_t> separatorSizeMetric = addNewMetric<om_t, scalar_t>(env, metrics);
  separatorSizeMetric->setName("separator size");
  separatorSizeMetric->setMetricValue("unsolved", -111); // placeholder
  separatorSizeMetric->setMetricValue("solved", -111); // placeholder

  env->debug(DETAILED_STATUS, "Exiting orderingMetrics"); // end
}

/*! \brief Print out header info for ordering metrics.
 */
template <typename scalar_t>
void printOrderingMetricsHeader(std::ostream &os)
{
  os << "Ordering Metrics" << std::endl;
  OrderingMetrics<scalar_t>::printHeader(os);
}

/*! \brief Print out list of ordering metrics.
 */
template <typename scalar_t>
void printOrderingMetrics(std::ostream &os,
  const ArrayView<RCP<BaseClassMetrics<scalar_t>>> &infoList)
{
  printOrderingMetricsHeader<scalar_t>(os);
  for (int i=0; i < infoList.size(); i++) {
    if (infoList[i]->getName() != METRICS_UNSET_STRING) {
      infoList[i]->printLine(os);
    }
  }
  os << std::endl;
}

/*! \brief Print out header and a single ordering metric.
 */
template <typename scalar_t>
void printOrderingMetrics(std::ostream &os,
  RCP<BaseClassMetrics<scalar_t>> metricValue)
{
  printOrderingMetricsHeader<scalar_t>(os);
  metricValue->printLine(os);
}

} //namespace Zoltan2


#endif
