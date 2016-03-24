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

/*! \file Zoltan2_EvaluatePartition.hpp
 *  \brief Defines the EvaluatePartition class.
 */

#ifndef ZOLTAN2_EVALUATEPARTITION_HPP
#define ZOLTAN2_EVALUATEPARTITION_HPP

#include <Zoltan2_Metric.hpp>
#include <Zoltan2_PartitioningSolution.hpp>

namespace Zoltan2{

/*! \brief A class that computes and returns quality metrics.
 *  \todo For some problems it will be necessary to build the
 *          Model again in order to compute metrics.  For now
 *          we don't have any problems like that.
    \todo write a unit test for this class
 */

  template <typename Adapter>
  class EvaluatePartition {

private:

  typedef typename Adapter::base_adapter_t base_adapter_t;
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::part_t part_t;
  typedef typename Adapter::scalar_t scalar_t;
  typedef StridedData<lno_t, scalar_t> input_t;

  const RCP<const Environment> env_;

  part_t numGlobalParts_;           // desired
  part_t targetGlobalParts_;        // actual
  part_t numNonEmpty_;              // of actual

  ArrayRCP<MetricValues<scalar_t> > metrics_;
  ArrayRCP<const MetricValues<scalar_t> > metricsConst_;
  ArrayRCP<GraphMetricValues<scalar_t> > graphMetrics_;
  ArrayRCP<const GraphMetricValues<scalar_t> > graphMetricsConst_;

public:

  /*! \brief Constructor
      \param env   the problem environment
      \param problemComm  the problem communicator
      \param ia the problem input adapter
      \param soln  the solution
      \param graphModel the graph model

      The constructor does global communication to compute the metrics.
      The rest of the  methods are local.
   */
  EvaluatePartition(const RCP<const Environment> &env,
    const RCP<const Comm<int> > &problemComm,
		    const Adapter *ia, 
    const PartitioningSolution<Adapter> *soln,
    const RCP<const GraphModel<typename Adapter::base_adapter_t> > &graphModel=
		    Teuchos::null);

  /*! \brief Return the metric values.
   *  \param values on return is the array of values.
   */
  ArrayRCP<const MetricValues<scalar_t> > getMetrics() const{
    //BDD return metricsConst_;
      if(metricsConst_.is_null()) return metrics_;
      return metricsConst_;
  }

  /*! \brief Return the object count imbalance.
   *  \param imbalance on return is the object count imbalance.
   */
  void getObjectCountImbalance(scalar_t &imbalance) const{
    imbalance = metrics_[0].getMaxImbalance();
  }

  /*! \brief Return the object normed weight imbalance.
   *  \param imbalance on return is the object normed weight imbalance.
   *  If there were no weights, this is the object count imbalance.
   *  If there was one weight, it is the imbalance with respect to that weight.
   */
  void getNormedImbalance(scalar_t &imbalance) const{
    if (metrics_.size() > 1)
      imbalance = metrics_[1].getMaxImbalance();
    else 
      imbalance = metrics_[0].getMaxImbalance();
  }

  /*! \brief Return the imbalance for the requested weight.
   *  \param imbalance on return is the requested value.
   *  \param idx is the weight index requested, ranging from zero
   *     to one less than the number of weights provided in the input.
   *  If there were no weights, this is the object count imbalance.
   */
  void getWeightImbalance(scalar_t &imbalance, int idx=0) const{
    imbalance = 0;
    if (metrics_.size() > 2)  // idx of multiple weights
      imbalance = metrics_[idx+2].getMaxImbalance();
    else if (metrics_.size() == 2)   //  only one weight
      imbalance = metrics_[1].getMaxImbalance();
    else                       // no weights, return object count imbalance
      imbalance = metrics_[0].getMaxImbalance();
  }  

  /*! \brief Print all the metrics
   */
  void printMetrics(std::ostream &os) const {
    Zoltan2::printMetrics<scalar_t, part_t>(os, 
      targetGlobalParts_, numGlobalParts_, numNonEmpty_, 
      metrics_.view(0, metrics_.size()));
  }

  /*! \brief Return the graph metric values.
   *  \param values on return is the array of values.
   */
  ArrayRCP<const GraphMetricValues<scalar_t> > getGraphMetrics() const{
      if(graphMetricsConst_.is_null()) return graphMetrics_;
      return graphMetricsConst_;
  }

  /*! \brief Return the max cut for the requested weight.
   *  \param cut on return is the requested value.
   *  \param idx is the weight index reqested, ranging from zero
   *     to one less than the number of weights provided in the input.
   *  If there were no weights, this is the cut count.
   */
  void getWeightCut(scalar_t &cut, int idx=0) const{
    if (graphMetrics_.size() < idx)  // idx too high
      cut = graphMetrics_[graphMetrics_.size()-1].getGlobalMax();
    else if (idx < 0)   //  idx too low
      cut = graphMetrics_[0].getGlobalMax();
    else                       // idx weight
      cut = graphMetrics_[idx].getGlobalMax();
  }

  /*! \brief Print all the metrics
   */
  void printGraphMetrics(std::ostream &os) const {
    Zoltan2::printMetrics<scalar_t, part_t>(os, 
      targetGlobalParts_, numGlobalParts_, 
      graphMetrics_.view(0, graphMetrics_.size()));
  }
};

  template <typename Adapter>
  EvaluatePartition<Adapter>::EvaluatePartition(
  const RCP<const Environment> &env,
  const RCP<const Comm<int> > &problemComm,
  const Adapter *ia, 
  const PartitioningSolution<Adapter> *soln,
  const RCP<const GraphModel<typename Adapter::base_adapter_t> > &graphModel):
    env_(env), numGlobalParts_(0), targetGlobalParts_(0), numNonEmpty_(0),
    metrics_(),  metricsConst_(), graphMetrics_(), graphMetricsConst_()
{

  env->debug(DETAILED_STATUS, std::string("Entering EvaluatePartition"));
  env->timerStart(MACRO_TIMERS, "Computing metrics");

  // When we add parameters for which weights to use, we
  // should check those here.  For now we compute metrics
  // using all weights.

  const Teuchos::ParameterList &pl = env->getParameters();

  multiCriteriaNorm mcnorm = normBalanceTotalMaximum;
  const Teuchos::ParameterEntry *pe = pl.getEntryPtr("partitioning_objective");
  if (pe){
    std::string strChoice = pe->getValue<std::string>(&strChoice);
    if (strChoice == std::string("multicriteria_minimize_total_weight"))
      mcnorm = normMinimizeTotalWeight;
    else if (strChoice == std::string("multicriteria_minimize_maximum_weight"))
      mcnorm = normMinimizeMaximumWeight;
  } 

  const RCP<const base_adapter_t> bia =
    rcp(dynamic_cast<const base_adapter_t *>(ia), false);

  try{
    objectMetrics<Adapter>(env, problemComm, mcnorm, ia, soln, graphModel,
			   numGlobalParts_, numNonEmpty_,metrics_);
  }
  Z2_FORWARD_EXCEPTIONS;

  if (soln)
  targetGlobalParts_ = soln->getTargetGlobalNumberOfParts();
  else
    targetGlobalParts_ = problemComm->getSize();

  env->timerStop(MACRO_TIMERS, "Computing metrics");

  BaseAdapterType inputType = bia->adapterType();

  if (inputType == GraphAdapterType ||
      inputType == MatrixAdapterType ||
      inputType == MeshAdapterType){
    env->timerStart(MACRO_TIMERS, "Computing graph metrics");
    // When we add parameters for which weights to use, we
    // should check those here.  For now we compute graph metrics
    // using all weights.

    std::bitset<NUM_MODEL_FLAGS> modelFlags;

    // Create a GraphModel based on input data.

    RCP<GraphModel<base_adapter_t> > graph;
    if (graphModel == Teuchos::null)
      graph = rcp(new GraphModel<base_adapter_t>(bia, env, problemComm, 
						 modelFlags));

    // Local number of objects.

    size_t numLocalObjects = bia->getLocalNumIDs();

    // Parts to which objects are assigned.

    const part_t *parts;
    if (soln) {
      parts = soln->getPartListView();
      env->localInputAssertion(__FILE__, __LINE__, "parts not set",
        ((numLocalObjects == 0) || parts), BASIC_ASSERTION);
    } else {
      part_t *procs = new part_t [numLocalObjects];
      for (size_t i=0; i<numLocalObjects; i++) procs[i]=problemComm->getRank();
      parts = procs;
    }
    ArrayView<const part_t> partArray(parts, numLocalObjects);

    ArrayRCP<scalar_t> globalSums;

    if (graphModel == Teuchos::null) {
      try{
	globalWeightedCutsByPart<Adapter>(env,
					  problemComm, graph, partArray,
					  numGlobalParts_, graphMetrics_,
					  globalSums);
      }
      Z2_FORWARD_EXCEPTIONS;
    } else {
      try{
	globalWeightedCutsByPart<Adapter>(env,
					  problemComm, graphModel, partArray,
					  numGlobalParts_, graphMetrics_,
					  globalSums);
      }
      Z2_FORWARD_EXCEPTIONS;
    }

    env->timerStop(MACRO_TIMERS, "Computing graph metrics");
  }
  env->debug(DETAILED_STATUS,
	     std::string("Exiting EvaluatePartition"));
}

}   // namespace Zoltan2

#endif
