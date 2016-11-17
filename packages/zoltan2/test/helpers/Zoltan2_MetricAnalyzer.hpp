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
#include <Zoltan2_EvaluateBaseClass.hpp>
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

template<class Adapter>
class MetricAnalyzer {
public:
  typedef Zoltan2::MetricAnalyzerInfo<typename Adapter::scalar_t> metric_analyzer_info_t;

  /// \brief Analyze metrics for a problem based on a range of tolerances
  ///
  /// @param metricsPlist parameter list defining tolerances
  /// @param problem the problem whose metrics are to be analyzed
  /// @param[out] msg_stream to return information from the analysis
  ///
  /// @return returns a boolean value indicated pass/failure.
  static bool analyzeMetrics(
    RCP<Zoltan2::EvaluateBaseClass<Adapter>> pEvaluate,
    const ParameterList &metricsParameters, 
    std::ostringstream & msg_stream )
  {
    if (metricsParameters.numParams() == 0) {
      // specification is that we do nothing - we may just be testing our status
      return true;
    }

    bool bAllPassed = true;

    std::vector<metric_analyzer_info_t> metricInfoSet;
    LoadMetricInfo(metricInfoSet, pEvaluate, metricsParameters);

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

  static void LoadMetricInfo(
    std::vector<metric_analyzer_info_t> & metricInfoSet,
    RCP<Zoltan2::EvaluateBaseClass<Adapter>> pEvaluate,
    const ParameterList &metricsParameters) {

    // at this point we should be looking at a metricsPlist with the following 
    // format - note that weight is optional

    //      <ParameterList name="metriccheck1">
    //        <Parameter name="check" type="string" value="imbalance"/>
    //        <Parameter name="lower" type="double" value="0.99"/>
    //        <Parameter name="upper" type="double" value="1.4"/>
    //      </ParameterList>
    //      <ParameterList name="metriccheck2">
    //        <Parameter name="check" type="string" value="imbalance"/>
    //        <Parameter name="weight" type="int" value="0"/>
    //        <Parameter name="lower" type="double" value="0.99"/>
    //        <Parameter name="upper" type="double" value="1.4"/>
    //      </ParameterList>

    // first let's get a list of all the headings, so "metriccheck1",
    // "metriccheck2" in this case. I've currently got this enforcing those
    // names strictly to make sure formatting is correct. But really the
    // headings could just be any unique names and are arbitrary
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
      metric_analyzer_info_t metricInfo = pEvaluate->getMetricAnalyzerInfo(
        metricsParameters.sublist(headingName));
      metricInfoSet.push_back(metricInfo);
      ++headingIndex;
    }
  }

private:

  static bool executeMetricCheck(const metric_analyzer_info_t & metricInfo,
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


#endif //ZOLTAN2_METRIC_ANALYZER_HPP
