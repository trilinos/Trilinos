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

/* \file Zoltan2_MetricAnalyzer.hpp
 * \brief Used by the Zoltan2 test driver for running \
          simple pass fail tests based on ranges of problem metrics.
 */
#ifndef ZOLTAN2_METRIC_ANALYZER_HPP
#define ZOLTAN2_METRIC_ANALYZER_HPP

#include <Zoltan2_TestHelpers.hpp>
#include <Zoltan2_Typedefs.hpp>
#include <Zoltan2_EvaluateOrdering.hpp>
#include <Zoltan2_EvaluatePartition.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_XMLObject.hpp>
#include <Teuchos_FileInputSource.hpp>

#include <sstream>
#include <string>
#include <iostream>

using Teuchos::ParameterList;
using Teuchos::Comm;
using Teuchos::RCP;
using Teuchos::ArrayRCP;
using namespace Zoltan2_TestingFramework;

struct MetricAnalyzerInfo
{
  double theValue;
  double upperValue;
  double lowerValue;
  bool bFoundUpperBound;
  bool bFoundLowerBound;
  std::string parameterDescription;
};

template<class Adapter>
class MetricAnalyzer {

protected:
  RCP<Zoltan2::EvaluateBaseClass<Adapter>> evaluate_;

public:
  #define KEYWORD_PARAMETER_NAME "check" // should be the first entry
  #define UPPER_PARAMETER_NAME "upper"
  #define LOWER_PARAMETER_NAME "lower"

  /*! \brief MetricAnalyzer constructor takes an EvaluateBaseClass such
       as EvaluateOrdering or EvaluatePartition.
   */
  MetricAnalyzer(RCP<Zoltan2::EvaluateBaseClass<Adapter>> evaluate)
    : evaluate_(evaluate) {
  }

  /*! \brief analyzeMetrics for a problem based on a range of tolerances
       @param metricsPlist parameter list defining tolerances
       @param[out] msg_stream to return information from the analysis
       @return returns a boolean value indicated pass/failure.
   */
  bool analyzeMetrics(const ParameterList &metricsParameters,
    std::ostringstream & msg_stream )
  {
    if (metricsParameters.numParams() == 0) {
      // specification is that we do nothing - we may just be testing our status
      return true;
    }

    bool bAllPassed = true;

    std::vector<MetricAnalyzerInfo> metricInfoSet;
    LoadMetricInfo(metricInfoSet, metricsParameters);

    int countFailedMetricChecks = 0;
    for (auto metricInfo = metricInfoSet.begin();
         metricInfo != metricInfoSet.end(); ++metricInfo) {
      if (!MetricAnalyzer::executeMetricCheck(*metricInfo, msg_stream)) {
        ++countFailedMetricChecks;
      }
    }

    // this code prints a summary of all metric checks and indicates how many
    // failed, if any did fail
    if(countFailedMetricChecks == 0) {
      msg_stream << metricsParameters.numParams() << " out of " <<
        metricsParameters.numParams() << " metric checks"  << " PASSED."
        << std::endl;
    }
    else {
      msg_stream << countFailedMetricChecks << " out of " <<
        metricsParameters.numParams() << " metric checks " << " FAILED."
        << std::endl;
      bAllPassed = false;
    }
    msg_stream << std::endl; // cosmetic spacer
    return bAllPassed;
  }

  /*! \brief getMetricValue is abstract and the derived class must define
   * the proper method to check optional values and determine the resulting
   * scalar value. The derived class will also throw if formatting is incorrect.
   */
  virtual MetricAnalyzerInfo getMetricResult(
    const ParameterList & metricCheckParameters, std::string keyWord) const = 0;

  void LoadMetricInfo(std::vector<MetricAnalyzerInfo> & metricInfoSet,
    const ParameterList &metricsParameters) {

    int headingIndex = 1;

    for (auto iterateArbitraryHeadingNames = metricsParameters.begin(); 
         iterateArbitraryHeadingNames != metricsParameters.end();
         ++iterateArbitraryHeadingNames) {
      auto headingName = metricsParameters.name(iterateArbitraryHeadingNames);

      // we could be flexible on these headers but for now let's enforce it to
      // get any convention inconsistencies cleaned up
      std::string expectedHeadingName = "metriccheck" + std::to_string(headingIndex);
      if( expectedHeadingName != headingName) {
        throw std::logic_error(
          "The parameter list expected to find a heading with name '"
          + expectedHeadingName + "' but instead found '" + headingName );
      }

      // get the parameters specific to the check we want to run
      MetricAnalyzerInfo metricInfo = getMetricAnalyzerInfo(
        metricsParameters.sublist(headingName));
      metricInfoSet.push_back(metricInfo);
      ++headingIndex;
    }
  }

  /*! \brief getMetricAnalyzerInfo is responsible for reading a metric value
   * and then checking it against upper and lower bounds. Any fomratting errors
   * should throw.
   */
  MetricAnalyzerInfo getMetricAnalyzerInfo(
    const ParameterList & metricCheckParameters) const {

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
    MetricAnalyzerInfo result
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

  /*! \brief Return true for any names we accept.
   */
  virtual bool isMetricCheckNameValid(std::string metricCheckName) const {
    // EvaluatePartition will have special key words 'weight' and 'normed'
    return false;
  }

private:
   bool executeMetricCheck(const MetricAnalyzerInfo & metricInfo,
    std::ostringstream &msg_stream)
  {
    bool bDoesThisTestPass = true; // will set this false if a test fails
    if (metricInfo.bFoundUpperBound && metricInfo.bFoundLowerBound) {
      if (metricInfo.theValue < metricInfo.lowerValue ||
          metricInfo.theValue > metricInfo.upperValue) {
        msg_stream << "FAILED: " << metricInfo.parameterDescription
                   << " value: " << metricInfo.theValue << " is not in range: "
                   << metricInfo.lowerValue << " to "
                   << metricInfo.upperValue << std::endl;
        bDoesThisTestPass = false;
      }
      else {
        msg_stream << "Success: " << metricInfo.parameterDescription
                   << " value: " << metricInfo.theValue << " is in range: "
                   << metricInfo.lowerValue << " to "
                   << metricInfo.upperValue << std::endl;
      }
    }
    else if (metricInfo.bFoundUpperBound) {
      if (metricInfo.theValue > metricInfo.upperValue) {
        msg_stream << "FAILED: " << metricInfo.parameterDescription
                   << " value: " << metricInfo.theValue << " is not below "
                   << metricInfo.upperValue << std::endl;
        bDoesThisTestPass = false;
      }
      else {
        msg_stream << "Success: " << metricInfo.parameterDescription
                   << " value: " << metricInfo.theValue << " is below: "
                   << metricInfo.upperValue << std::endl;
      }
    }
    else if (metricInfo.bFoundLowerBound) {
      if (metricInfo.theValue < metricInfo.lowerValue) {
        msg_stream << "FAILED: " << metricInfo.parameterDescription
                   << " value: " << metricInfo.theValue << " is not above "
                   << metricInfo.lowerValue << std::endl;
        bDoesThisTestPass = false;
      }
      else {
        msg_stream << "Success: " << metricInfo.parameterDescription  
                   << " value: " << metricInfo.theValue << " is above: "
                   << metricInfo.lowerValue << std::endl;
      }
    }
    return bDoesThisTestPass;
  }
};

template<class Adapter>
class MetricAnalyzerEvaluatePartition : public MetricAnalyzer<Adapter> {
public:
  // defines for metric checks - these are key words used in xml
  #define WEIGHT_PARAMETER_NAME "weight"
  #define NORMED_PARAMETER_NAME "normed"

  /*! \brief MetricAnalyzerEvaluatePartition constructor.
   */
  MetricAnalyzerEvaluatePartition(
    RCP<Zoltan2::EvaluateBaseClass<Adapter>> evaluate)
    : MetricAnalyzer<Adapter>(evaluate) {
  }

  /*! \brief Reads a metric value for bounds checking.
   * Handle any special optional parameters.
   */
  virtual MetricAnalyzerInfo getMetricResult(
    const ParameterList & metricCheckParameters, std::string keyWord) const {

    RCP<Zoltan2::EvaluatePartition<Adapter>> pEvaluatePartition =
      Teuchos::rcp_dynamic_cast<Zoltan2::EvaluatePartition<Adapter>>(this->evaluate_);

    MetricAnalyzerInfo result;

    // didn't want to duplicate this value - a weight index should be 0 or
    // larger but it's optional to specify it
    #define UNDEFINED_PARAMETER_INT_INDEX -1

    // Read the weight index parameter and throw if not a good format
    // This is an optional parameter for EvaluatePartition
    int weightIndex = UNDEFINED_PARAMETER_INT_INDEX;
    if( metricCheckParameters.isParameter(WEIGHT_PARAMETER_NAME)) {
      weightIndex = metricCheckParameters.get<int>(WEIGHT_PARAMETER_NAME);
      if( weightIndex < 0 ) {
        throw std::logic_error( "Optional weight index was specified as: " +
          std::to_string(weightIndex) +
          "   Weight index must be 0 or positive." );
      }
    }

    // Read the norm index and throw if not a good format
    // This is an optional parameter for EvaluatePartition
    int normedSetting = UNDEFINED_PARAMETER_INT_INDEX;
    if( metricCheckParameters.isParameter(NORMED_PARAMETER_NAME)) {
      bool bNormSetting = metricCheckParameters.get<bool>(NORMED_PARAMETER_NAME);
      normedSetting = bNormSetting ? 1 : 0;
      if( normedSetting != 0 && normedSetting != 1 ) {
        throw std::logic_error( "Optional normed parameter was specified as: "
          + std::to_string(normedSetting) +
          "   Normed parameter must be true or false." );
      }
    }

    if( weightIndex != UNDEFINED_PARAMETER_INT_INDEX &&
      normedSetting != UNDEFINED_PARAMETER_INT_INDEX ) {
      throw std::logic_error( "Both parameters 'normed' and 'weight' were "
        " specified. They should never appear together." );
    }

    // these define key names which convert to an API call
    #define API_STRING_getWeightImbalance "imbalance"
    #define API_STRING_getTotalEdgeCuts "total edge cuts"
    #define API_STRING_getMaxEdgeCuts "max edge cuts"

    // throw if normed set and weight is set
    if( keyWord != API_STRING_getWeightImbalance &&
      normedSetting != UNDEFINED_PARAMETER_INT_INDEX ) {
      throw std::logic_error( "'normed' was specified but this only has meaning"
       " for the 'imbalance' parameter." );
    }

    // Enforcing parallel usage to the API calls exist in EvaluatePartition
    if (keyWord == API_STRING_getWeightImbalance) {
      if( weightIndex == UNDEFINED_PARAMETER_INT_INDEX ) {
        if( normedSetting == 1 ) {
          result.theValue = pEvaluatePartition->getNormedImbalance();
        }
        else {
          result.theValue = pEvaluatePartition->getObjectCountImbalance();
        }
      }
      else {
        // this will get the proper index specified
        result.theValue = pEvaluatePartition->getWeightImbalance(weightIndex);
      }
    }
    else if (keyWord == API_STRING_getTotalEdgeCuts) {
      if( weightIndex == UNDEFINED_PARAMETER_INT_INDEX ) {
        result.theValue = pEvaluatePartition->getTotalEdgeCut();
      }
      else {
        result.theValue = pEvaluatePartition->getTotalWeightEdgeCut(weightIndex);
      }
    }
    else if (keyWord == API_STRING_getMaxEdgeCuts) {
      if( weightIndex == UNDEFINED_PARAMETER_INT_INDEX ) {
        result.theValue = pEvaluatePartition->getMaxEdgeCut();
      }
      else {
        result.theValue = pEvaluatePartition->getMaxWeightEdgeCut(weightIndex);
      }
    }
    else {
      // we have found an invalid key word - throw an error
      throw std::logic_error( "The parameter '" +
        std::string(KEYWORD_PARAMETER_NAME) + "' was specified as '" +
        keyWord + "' which is not understood." );
    }

    result.parameterDescription = keyWord;
    if( weightIndex != UNDEFINED_PARAMETER_INT_INDEX ) {
      result.parameterDescription = result.parameterDescription +
        " (weight: " + std::to_string(weightIndex) + ")";
    }
    else if( normedSetting != UNDEFINED_PARAMETER_INT_INDEX ) {
      // throw above would catch the case where both of these were set
      result.parameterDescription = result.parameterDescription + " (normed: "
        + ( ( normedSetting == 0 ) ? "false" : "true" ) + ")";
    }

    return result;
  }

  /*! \brief Return true for any names we accept.
   */
  virtual bool isMetricCheckNameValid(std::string metricCheckName) const {
    return (metricCheckName == WEIGHT_PARAMETER_NAME ||
            metricCheckName == NORMED_PARAMETER_NAME);
  }
};

template<class Adapter>
class MetricAnalyzerEvaluateOrdering : public MetricAnalyzer<Adapter> {
public:

  // these define key names which convert to an API call
  #define API_STRING_getBandwidth "bandwidth"
  #define API_STRING_getEnvelope "envelope"
  #define API_STRING_getSeparatorSize "separator size"

  /*! \brief MetricAnalyzerEvaluatePartition constructor.
   */
  MetricAnalyzerEvaluateOrdering(
    RCP<Zoltan2::EvaluateBaseClass<Adapter>> evaluate)
    : MetricAnalyzer<Adapter>(evaluate) {
  }

  /*! \brief Reads a metric value for bounds checking.
   * Handle any special optional parameters.
   */
  virtual MetricAnalyzerInfo getMetricResult(
    const ParameterList & metricCheckParameters, std::string keyWord) const {

    RCP<Zoltan2::EvaluateOrdering<Adapter>> pEvaluateOrdering =
      Teuchos::rcp_dynamic_cast<Zoltan2::EvaluateOrdering<Adapter>>(this->evaluate_);

    MetricAnalyzerInfo result;

    if (keyWord == API_STRING_getBandwidth) {
      result.theValue = pEvaluateOrdering->getBandwidth();
    }
    else if (keyWord == API_STRING_getEnvelope) {
      result.theValue = pEvaluateOrdering->getEnvelope();
    }
    else if (keyWord == API_STRING_getSeparatorSize) {
      result.theValue = pEvaluateOrdering->getSeparatorSize();
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

#endif //ZOLTAN2_METRIC_ANALYZER_HPP
