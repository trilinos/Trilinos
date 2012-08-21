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

/*! \file Zoltan2_PartitioningSolutionQuality.hpp
 *  \brief Defines the PartitioningSolutionQuality class.
 */

#ifndef ZOLTAN2_SOLUTIONQUALITY_HPP
#define ZOLTAN2_SOLUTIONQUALITY_HPP

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
  class PartitioningSolutionQuality {

private:

  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::gid_t gid_t;
  typedef typename Adapter::scalar_t scalar_t;

  const RCP<const Environment> env_;

  partId_t numGlobalParts_;           // desired
  partId_t targetGlobalParts_;        // actual
  partId_t numNonEmpty_;              // of actual

  ArrayRCP<MetricValues<scalar_t> > metrics_;
  ArrayRCP<const MetricValues<scalar_t> > metricsConst_;

public:

  /*! \brief Constructor
      \param env   the problem environment
      \param problemComm  the problem communicator
      \param ia the problem input adapter
      \param soln  the solution

      The constructor does global communication to compute the metrics.
      The rest of the  methods are local.
   */
  PartitioningSolutionQuality(const RCP<const Environment> &env,
    const RCP<const Comm<int> > &problemComm,
    const RCP<const Adapter> &ia, 
    const RCP<const PartitioningSolution<Adapter> > &soln);

  /*! \brief Return the metric values.
   *  \param values on return is the array of values.
   */
  void getMetrics(ArrayRCP<const MetricValues<scalar_t> > &values) const{
    values = metricsConst_;
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

  /*! \brief Return the imbalance for the requested weight dimension.
   *  \param imbalance on return is the requested value.
   *  \param dim is the weight dimension requested, ranging from zero
   *     to one less than the number of weights provided in the input.
   *  If there were no weights, this is the object count imbalance.
   */
  void getWeightImbalance(scalar_t &imbalance, int dim=0) const{
    imbalance = 0;
    if (metrics_.size() > 2)  // dimension dim of multiple weights
      imbalance = metrics_[dim+2].getMaxImbalance();
    else if (metrics_.size() == 2)   //  only one weight
      imbalance = metrics_[1].getMaxImbalance();
    else                       // no weights, return object count imbalance
      imbalance = metrics_[0].getMaxImbalance();
  }  

  /*! \brief Print all the metrics
   */
  void printMetrics(ostream &os) const {
    Zoltan2::printMetrics<scalar_t>(os, 
      targetGlobalParts_, numGlobalParts_, numNonEmpty_, 
      metrics_.view(0, metrics_.size()));
  }
};

template <typename Adapter>
  PartitioningSolutionQuality<Adapter>::PartitioningSolutionQuality(
  const RCP<const Environment> &env,
  const RCP<const Comm<int> > &problemComm,
  const RCP<const Adapter> &ia, 
  const RCP<const PartitioningSolution<Adapter> > &soln):
    env_(env), numGlobalParts_(0), targetGlobalParts_(0), numNonEmpty_(0),
    metrics_(),  metricsConst_()
{

  env->debug(DETAILED_STATUS, string("Entering PartitioningSolutionQuality"));
  env->timerStart(MACRO_TIMERS, "Computing metrics");

  // When we add parameters for which weights to use, we
  // should check those here.  For now we compute metrics
  // using all weights.

  const Teuchos::ParameterList &pl = env->getParameters();
  multiCriteriaNorm mcnorm = normBalanceTotalMaximum;
  bool balanceCountOnly=false;

  bool isSet;
  string strChoice;

  getParameterValue<string>(pl, "partitioning", "objective", 
    isSet, strChoice);

  if (isSet){
    if (strChoice == string("balance_object_count"))
      balanceCountOnly=true;
    else if (strChoice == string("multicriteria_minimize_total_weight"))
      mcnorm = normMinimizeTotalWeight;
    else if (strChoice == string("multicriteria_minimize_maximum_weight"))
      mcnorm = normMinimizeMaximumWeight;
  } 

  try{
    objectMetrics<Adapter>(env, problemComm, mcnorm, ia, soln,
      numGlobalParts_, numNonEmpty_, metrics_);
  }
  Z2_FORWARD_EXCEPTIONS;

  targetGlobalParts_ = soln->getTargetGlobalNumberOfParts();

  env->timerStop(MACRO_TIMERS, "Computing metrics");
  env->debug(DETAILED_STATUS, string("Exiting PartitioningSolutionQuality"));
}

}   // namespace Zoltan2

#endif
