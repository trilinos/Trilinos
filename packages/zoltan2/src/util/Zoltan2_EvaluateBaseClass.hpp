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

/*! \file Zoltan2_EvaluateBaseClass.hpp
 *  \brief Base class for the EvaluatePartition and EvaluateOrdering classes.
 */

#ifndef ZOLTAN2_EVALUATEBASECLASS_HPP
#define ZOLTAN2_EVALUATEBASECLASS_HPP

#include <Zoltan2_Metric.hpp>

namespace Zoltan2{

/*! \brief A base class for EvaluatePartition, EvaluateOrdering, ...
 */

template<typename scalar_t>
struct MetricAnalyzerInfo
{
  scalar_t theValue;
  double upperValue;
  double lowerValue;
  bool bFoundUpperBound;
  bool bFoundLowerBound;
  std::string parameterDescription;
};

template <typename Adapter>
class EvaluateBaseClass {

public:
  typedef typename Adapter::scalar_t scalar_t;
  typedef BaseClassMetrics<scalar_t> base_metric_t;
  typedef ArrayRCP<RCP<base_metric_t>> base_metric_array_t;

  EvaluateBaseClass(): metricsBase_() {}

  /*! \brief Return the full metric list which will be mixed types.
   */
  base_metric_array_t & getAllMetrics() {
    return metricsBase_;
  }

  /*! \get the unique class names which are actually used
   */
  std::vector<std::string> getMetricTypes() const {
    std::set<std::string> temporarySet;
    std::vector<std::string> returnVector; // Michel did this because he wants
                                           // set behavior but preserve the
                                           // ordering and a vector is
                                           // convenient later since we may
                                           // also be looking at the vector
                                           // set of all possible names from
                                           // BaseClassMetrics
                                           // other option is to use only types
                                           // that are currently being used
    for( int n = 0; n < metricsBase_.size(); ++n ) {
      std::string checkName = metricsBase_[n]->getMetricType();
      if (temporarySet.find(checkName) == temporarySet.end()) {
        temporarySet.insert(checkName);
        returnVector.push_back(checkName);
      }
    }
    return returnVector;
  }


  /*! \brief Get a single metric result. This should only be used if only 1
   * metric with the type and name is expected
   */
  RCP<base_metric_t> getMetric(std::string metricType, std::string name) const {
    RCP<base_metric_t> result;
    for(auto n = 0; n < metricsBase_.size(); ++n) {
      if( metricsBase_[n]->getMetricType() == metricType &&
        metricsBase_[n]->getName() == name ) {
        if (result != Teuchos::null) { // check if more than 1?
          throw std::logic_error( "getMetric called for metric class: '" +
            metricType + "' and metric name: '" + name + "' but more than" +
            " 1 match was found which is not expected." );
        }
        result = metricsBase_[n];
      }
    }
    if (result == Teuchos::null) { // check if didn't find any?
      throw std::logic_error( "getMetric called for metric class: '" +
        metricType + "' and metric name: '" + name + "' but no matching" +
        " metric was found which is not expected." );
    }
    return result;
  }

  /*! \brief Return the metric list for types matching the given metric type.
   */
  ArrayView<RCP<base_metric_t>> getAllMetricsOfType(std::string metricType) const {
    // find the beginning and the end of the contiguous block
    // the list is an ArrayRCP and must preserve any ordering
    int beginIndex = -1;
    int sizeOfArrayView = 0;
    for(auto n = 0; n < metricsBase_.size(); ++n) {
      if( metricsBase_[n]->getMetricType() == metricType ) {
        if (beginIndex == -1) {
          beginIndex = int(n);
        }
        ++sizeOfArrayView;
      }
    }
    if (sizeOfArrayView == 0) {
      return ArrayView<RCP<base_metric_t> >(); // empty array view
    }
    return metricsBase_.view(beginIndex, sizeOfArrayView);
  }

  /*! \brief Print all the metrics based on the  metric object type
   */
  void printMetrics(std::ostream &os, bool bIncludeUnusedTypes = true) const {
    std::vector<std::string> types =
                    (bIncludeUnusedTypes ?
                     BaseClassMetrics<scalar_t>::static_allMetricNames_ :
                     getMetricTypes());
    for( auto metricType : types ) {
      printMetrics (os, metricType);
    }
  }

  /*! \brief Print all metrics of type metricType based on the metric object type
   */
  void printMetrics(std::ostream &os, std::string metricType) const {
    // Might need changing. The issue here is that they each have unique
    // parameters to pass.
    ArrayView<RCP<base_metric_t>> metrics = getAllMetricsOfType(metricType);

    // this could be a critical decision - do we want a blank table with
    // headers when the list is empty - for debugging that is probably better
    // but it's very messy to have lots of empty tables in the logs
    if (metrics.size() != 0) {
      callStaticPrintMetrics(os, metrics, metricType);
    }
  }

  /*! \brief print metric functions with each base class
   *  Abstract but may relax that for future classes.
   */
  virtual void callStaticPrintMetrics(std::ostream &os,
    ArrayView<RCP<base_metric_t>>, std::string metricType) const  = 0;

  /*! \brief getMetricAnalyzerInfo is responsible for reading a metric value
   * and then checking it against upper and lower bounds. Any fomratting errors
   * should throw.
   */
  MetricAnalyzerInfo<typename Adapter::scalar_t> getMetricAnalyzerInfo(
    const ParameterList & metricCheckParameters) const {

    #define KEYWORD_PARAMETER_NAME "check" // should be the first entry
    #define UPPER_PARAMETER_NAME "upper"
    #define LOWER_PARAMETER_NAME "lower"

    // first validate that all the string names in the metric check are correct
    for (auto iterateAllKeys = metricCheckParameters.begin();
         iterateAllKeys != metricCheckParameters.end(); ++iterateAllKeys) {
      auto checkName = metricCheckParameters.name(iterateAllKeys);

      bool bIsGeneralName = (checkName == KEYWORD_PARAMETER_NAME ||
                             checkName == UPPER_PARAMETER_NAME ||
                             checkName == LOWER_PARAMETER_NAME );

      if (!bIsGeneralName && !isMetricCheckNameValid(checkName)) {
        throw std::logic_error(
          "Key name: '" + checkName + "' is not understood" );
      }
    }

    if( !metricCheckParameters.isParameter(KEYWORD_PARAMETER_NAME)) {
     throw std::logic_error( "Matric check must specify a key name using "
       "the keyword " + std::string(KEYWORD_PARAMETER_NAME) );
    }

    std::string keyWord =
      metricCheckParameters.get<std::string>(KEYWORD_PARAMETER_NAME);

    // one of the parameters called "check" should define a string which is a
    // keyword which correlates to an API call forthis class.
    MetricAnalyzerInfo<typename Adapter::scalar_t> result
      = getMetricResult(metricCheckParameters, keyWord);

    // now we can obtain the upper and lower bounds for this test
    result.bFoundUpperBound =
      metricCheckParameters.isParameter(UPPER_PARAMETER_NAME);
    result.bFoundLowerBound =
      metricCheckParameters.isParameter(LOWER_PARAMETER_NAME);

    if (result.bFoundUpperBound) {
      result.upperValue =
        metricCheckParameters.get<double>(UPPER_PARAMETER_NAME);
    }
    if (result.bFoundLowerBound) {
      result.lowerValue =
        metricCheckParameters.get<double>(LOWER_PARAMETER_NAME);
    }

    return result;
  }

  /*! \brief getMetricValue is abstract and the derived class must define
   * the proper method to check optional values and determine the resulting
   * scalar value. The derived class will also throw if formatting is incorrect.
   */
  virtual MetricAnalyzerInfo<typename Adapter::scalar_t> getMetricResult(
    const ParameterList & metricCheckParameters, std::string keyWord) const = 0;

  /*! \brief isMetricCheckNameValid should return true for any special
   * metric names used by the derived class. This is abstract and must be
   * defined by the base class.
   */
  virtual bool isMetricCheckNameValid(std::string metricCheckName) const = 0;

  private:
    base_metric_array_t metricsBase_;
};

}   // namespace Zoltan2

#endif
