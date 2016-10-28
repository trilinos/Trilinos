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

/*! \file Zoltan2_EvaluateOrdering.hpp
 *  \brief Defines the Zoltan2_EvaluateOrdering.hpp class.
 */

#ifndef ZOLTAN2_EVALUATEORDERING_HPP
#define ZOLTAN2_EVALUATEORDERING_HPP

#include <Zoltan2_EvaluateBaseClass.hpp>
#include <Zoltan2_OrderingSolution.hpp>

namespace Zoltan2{

/*! \brief A class that computes and returns quality metrics.
 *  \A base class for the local and global ordering versions.
 */
template <typename Adapter>
class EvaluateOrdering : public EvaluateBaseClass<Adapter> {

  // these define key names which convert to an API call
  // the could go into getMetricResult except they are also being used for
  // the api calls - however that may change - right now I have the API calls
  // looking strictly for a single matching metric which is different than
  // the way EvaluatePartition handles things (has normed and weights)
  // so just need to find out exactly how we want to throw errors and organize
  #define API_STRING_getBandwidth "bandwidth"
  #define API_STRING_getEnvelope "envelope"
  #define API_STRING_getSeparatorSize "separator size"

private:
  using EvaluateBaseClass<Adapter>::getAllMetricsOfType;
  using EvaluateBaseClass<Adapter>::getAllMetrics;

  typedef typename Adapter::base_adapter_t base_adapter_t;
  typedef typename Adapter::lno_t lno_t;
  typedef typename Adapter::gno_t gno_t;
  typedef typename Adapter::part_t part_t;
  typedef typename Adapter::scalar_t scalar_t;
  typedef BaseClassMetrics<scalar_t> base_metric_t;

  void sharedConstructor(const Adapter *ia,
                         ParameterList *p,
                         const RCP<const Comm<int> > &problemComm,
                         const LocalOrderingSolution<lno_t> *localSoln,
                         const GlobalOrderingSolution<gno_t> *globalSoln);
public:

  /*! \brief Constructor where communicator is Teuchos default.
      \param ia the problem input adapter
      \param p the parameter list
      \param localSoln the local solution
      \param globalSoln the global solution
      The constructor does global communication to compute the metrics.
      The rest of the  methods are local.
   */
  EvaluateOrdering(
    const Adapter *ia,
    ParameterList *p,
    const LocalOrderingSolution<lno_t> *localSoln,
    const GlobalOrderingSolution<gno_t> *globalSoln)
  {
    RCP<const Comm<int> > problemComm = DefaultComm<int>::getComm();
    sharedConstructor(ia, p, problemComm, localSoln, globalSoln);
  }

  /*! \brief Constructor where Teuchos communicator is specified
      \param ia the problem input adapter
      \param p the parameter list
      \param problemComm  the problem communicator
      \param localSoln the local solution
      \param globalSoln the global solution
      The constructor does global communication to compute the metrics.
      The rest of the  methods are local.
   */
  EvaluateOrdering(
    const Adapter *ia,
    ParameterList *p,
    const RCP<const Comm<int> > &problemComm,
    const LocalOrderingSolution<lno_t> *localSoln,
    const GlobalOrderingSolution<gno_t> *globalSoln)
  {
    sharedConstructor(ia, p, problemComm, localSoln, globalSoln);
  }

#ifdef HAVE_ZOLTAN2_MPI
  /*! \brief Constructor for MPI builds
      \param ia the problem input adapter
      \param p the parameter list
      \param comm  the problem communicator
      \param localSoln the local solution
      \param globalSoln the global solution
      The constructor does global communication to compute the metrics.
      The rest of the  methods are local.
   */
  EvaluateOrdering(
    const Adapter *ia,
    ParameterList *p,
    MPI_Comm comm,
    const LocalOrderingSolution<lno_t> *localSoln,
    const GlobalOrderingSolution<gno_t> *globalSoln)
  {
    RCP<Teuchos::OpaqueWrapper<MPI_Comm> > wrapper =
        Teuchos::opaqueWrapper(comm);
    RCP<const Comm<int> > problemComm =
        rcp<const Comm<int> >(new Teuchos::MpiComm<int>(wrapper));
    sharedConstructor(ia, p, problemComm, localSoln, globalSoln);
  }
#endif

  // To do - resolve how to handle local and global
  lno_t getBandwidth() const {
    auto bandwidthMetric = EvaluateBaseClass<Adapter>::getMetric(
      ORDERING_METRICS_TYPE_NAME, "bandwidth"); // can throw
    return bandwidthMetric->getMetricValue("solved");
  }

  // To do - resolve how to handle local and global
  lno_t getEnvelops() const {
    auto envelopeMetric = EvaluateBaseClass<Adapter>::getMetric(
      ORDERING_METRICS_TYPE_NAME, "envelope"); // can throw
    return envelopeMetric->getMetricValue("solved");
  }

  // To do - resolve how to handle local and global
  lno_t getSeparationSize() const {
    auto sepSizeMetric = EvaluateBaseClass<Adapter>::getMetric(
      ORDERING_METRICS_TYPE_NAME, "separator size");  // can throw
    return sepSizeMetric->getMetricValue("solved");
  }

  /*! \brief Print all metrics of type metricType based on the metric object type
   *  Note that parent class currently suppresses this if the list is empty.
   */
  virtual void callStaticPrintMetrics(std::ostream &os,
    ArrayView<RCP<base_metric_t>> metrics, std::string metricType) const {
      if( metricType == ORDERING_METRICS_TYPE_NAME ) {
        Zoltan2::printOrderingMetrics<scalar_t>(os, metrics);
      }
  }

  /*! \brief Return true for any names we accept.
   */
  virtual bool isMetricCheckNameValid(std::string metricCheckName) const {
    // currently EvaluateOrdering does not implement any special anems
    return false;
  }

  /*! \brief Reads a metric value for bounds checking.
   * Handle any special optional parameters.
   */
  virtual MetricAnalyzerInfo<typename Adapter::scalar_t> getMetricResult(
    const ParameterList & metricCheckParameters, std::string keyWord) const {

    MetricAnalyzerInfo<typename Adapter::scalar_t> result;

    if (keyWord == API_STRING_getBandwidth) {
      result.theValue = getBandwidth();
    }
    else if (keyWord == API_STRING_getEnvelope) {
      result.theValue = getEnvelops();
    }
    else if (keyWord == API_STRING_getSeparatorSize) {
      result.theValue = getSeparationSize();
    }
    else {
      // we have found an invalid key word - throw an error
      throw std::logic_error( "The parameter '" +
        std::string(KEYWORD_PARAMETER_NAME) + "' was specified as '" +
        keyWord + "' which is not understood." );
    }

    result.parameterDescription = keyWord; // just give default name for now

    return result;
  }
};

// sharedConstructor
template <typename Adapter>
void EvaluateOrdering<Adapter>::sharedConstructor(
  const Adapter *ia,
  ParameterList *p,
  const RCP<const Comm<int> > &comm,
  const LocalOrderingSolution<lno_t> *localSoln,
  const GlobalOrderingSolution<gno_t> *globalSoln)
{
  RCP<const Comm<int> > problemComm = (comm == Teuchos::null) ?
    DefaultComm<int>::getComm() : comm;

  RCP<Environment> env;
  try{
    env = rcp(new Environment(*p, problemComm));
  }
  Z2_FORWARD_EXCEPTIONS

  env->debug(DETAILED_STATUS, std::string("Entering EvaluateOrdering"));
  env->timerStart(MACRO_TIMERS, "Computing ordering metrics");

  try{
    // May want to move these into the specific derived classes
    // But it depends on whether we eventually may have both types and perhaps
    // want to combine the metrics
    if(localSoln) {
      localOrderingMetrics<Adapter>(env, problemComm, ia, localSoln,
        getAllMetrics());
    }

    if(globalSoln) {
      throw std::logic_error("EvaluateOrdering not set up for global ordering.");
    }
  }
  Z2_FORWARD_EXCEPTIONS;
  env->timerStop(MACRO_TIMERS, "Computing ordering metrics");

  env->debug(DETAILED_STATUS, std::string("Exiting EvaluateOrdering"));
}

template <typename Adapter>
class EvaluateLocalOrdering : public EvaluateOrdering<Adapter> {
private:
  typedef typename Adapter::lno_t lno_t;

public:
  /*! \brief Constructor where communicator is Teuchos default.
      \param ia the problem input adapter
      \param p the parameter list
      \param localSoln the local solution
      The constructor does global communication to compute the metrics.
      The rest of the  methods are local.
   */
  EvaluateLocalOrdering(
    const Adapter *ia,
    ParameterList *p,
    const LocalOrderingSolution<lno_t> *localSoln) :
    EvaluateOrdering<Adapter>(ia, p, localSoln, nullptr) {}

  /*! \brief Constructor where Teuchos communicator is specified
      \param ia the problem input adapter
      \param p the parameter list
      \param problemComm  the problem communicator
      \param localSoln the local solution
      The constructor does global communication to compute the metrics.
      The rest of the  methods are local.
   */
  EvaluateLocalOrdering(
    const Adapter *ia,
    ParameterList *p,
    const RCP<const Comm<int> > &problemComm,
    const LocalOrderingSolution<lno_t> *localSoln) :
    EvaluateOrdering<Adapter>(ia, p, problemComm, localSoln, nullptr) {}

#ifdef HAVE_ZOLTAN2_MPI
  /*! \brief Constructor for MPI builds
      \param ia the problem input adapter
      \param p the parameter list
      \param comm  the problem communicator
      \param localSoln the local solution
      \param globalSoln the global solution
      The constructor does global communication to compute the metrics.
      The rest of the  methods are local.
   */
  EvaluateLocalOrdering(
    const Adapter *ia,
    ParameterList *p,
    MPI_Comm comm,
    const LocalOrderingSolution<lno_t> *localSoln) :
    EvaluateOrdering<Adapter>(ia, p, comm, localSoln, nullptr) {}
#endif
};

template <typename Adapter>
class EvaluateGlobalOrdering : public EvaluateOrdering<Adapter> {
private:
  typedef typename Adapter::gno_t gno_t;

public:
  /*! \brief Constructor where communicator is Teuchos default.
      \param ia the problem input adapter
      \param p the parameter list
      \param localSoln the local solution
      The constructor does global communication to compute the metrics.
      The rest of the  methods are local.
   */
  EvaluateGlobalOrdering(
    const Adapter *ia,
    ParameterList *p,
    const GlobalOrderingSolution<gno_t> *globalSoln) :
    EvaluateOrdering<Adapter>(ia, p, nullptr, globalSoln) {}

  /*! \brief Constructor where Teuchos communicator is specified
      \param ia the problem input adapter
      \param p the parameter list
      \param problemComm  the problem communicator
      \param localSoln the local solution
      The constructor does global communication to compute the metrics.
      The rest of the  methods are local.
   */
  EvaluateGlobalOrdering(
    const Adapter *ia,
    ParameterList *p,
    const RCP<const Comm<int> > &problemComm,
    const GlobalOrderingSolution<gno_t> *globalSoln) :
    EvaluateOrdering<Adapter>(ia, p, problemComm, nullptr, globalSoln) {}

#ifdef HAVE_ZOLTAN2_MPI
  /*! \brief Constructor for MPI builds
      \param ia the problem input adapter
      \param p the parameter list
      \param comm  the problem communicator
      \param localSoln the local solution
      \param globalSoln the global solution
      The constructor does global communication to compute the metrics.
      The rest of the  methods are local.
   */
  EvaluateGlobalOrdering(
    const Adapter *ia,
    ParameterList *p,
    MPI_Comm comm,
    const GlobalOrderingSolution<gno_t> *globalSoln) :
    EvaluateOrdering<Adapter>(ia, p, comm, nullptr, globalSoln) {}
#endif
};

}   // namespace Zoltan2

#endif
