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

template <class Problem>
class MetricAnalyzer {
public:
  /// \brief analyze metrics for a problem based on a range of tolerances
  ///
  /// @param metricsPlist parameter list defining tolerances
  /// @param problem the problem whose metrics are to be analyzed
  /// @param comm an RCP for a communicator
  /// @param[out] msg_stream a std::ostringstream stream to return information from the analysis
  ///
  /// @return returns a boolean value indicated pass/failure.
  static bool analyzeMetrics( const ParameterList &metricsPlist,
                              const Problem * problem,
                              const RCP<const Comm<int>> &comm,
                              std::ostringstream &msg_stream) {
    msg_stream << "Metric analysis is not supported for problems of this type." << std::endl;
    return false;
  }

private:

};

template<>
class MetricAnalyzer<partitioning_problem_t> {

public:
  /// \brief Analyze metrics for a problem based on a range of tolerances
  ///
  /// @param metricsPlist parameter list defining tolerances
  /// @param problem the problem whose metrics are to be analyzed
  /// @param comm an RCP for a communicator
  /// @param[out] msg_stream a std::ostringstream stream to return information from the analysis
  ///
  /// @return returns a boolean value indicated pass/failure.
  static bool analyzeMetrics( const ParameterList &metricsPlist,
                              const partitioning_problem_t * problem,
			      const RCP<const Zoltan2::EvaluatePartition
			      <basic_id_t> > &metricObject,
                              const RCP<const Comm<int>> &comm,
                              std::ostringstream &msg_stream) {

      ArrayRCP<const metric_t> metrics = (problem)->getMetrics();
      
      string test_name;
      bool all_tests_pass = true;
      for(int i = 0; i < metrics.size(); i++)
      {
        // print their names...
        test_name = metrics[i].getName();
        if(metricsPlist.isSublist(test_name))
        {
          if(!MetricAnalyzer<partitioning_problem_t>::MetricBoundsTest(comm,
                metrics[i],
                metricsPlist.sublist(test_name),
                msg_stream)) {
            all_tests_pass = false;
          }
        }
      }

      return all_tests_pass;
  }

private:

  /// \brief Preforms the analysis
  ///
  /// @param comm RCP for a communicator
  /// @param metric the metrics for the problem
  /// @param metricsPlist parameter list defining tolerances
  /// @param[out] msg_stream a std::ostringstream stream to return information from the analysis
  ///
  /// @return returns a boolean value indicated pass/failure.
  static bool MetricBoundsTest(const RCP<const Comm<int>> &comm,
                        const metric_t & metric,
                        const Teuchos::ParameterList & metricPlist,
                        ostringstream &msg)
  {
    // run a comparison of min and max agains a given metric
    // return an error message on failure
    bool pass = true;
    string test_name = metric.getName() + " test";
    double value = metric.getMaxImbalance();
    
    // Perfom tests
    if (metricPlist.isParameter("lower"))
    {
      double min = metricPlist.get<double>("lower");
      
      if(value < min)
      {
        if (comm->getRank() == 0)
          msg << test_name << " FAILED: imbalance per part, "
              << value << ", less than specified allowable minimum, " 
              << min << ".\n";
        pass = false;
      }else{
        if (comm->getRank() == 0)
          msg << test_name << " PASSED: imbalance per part, "
              << value << ", greater than specified allowable minimum, " 
              << min << ".\n";
      }
    }
    
    if(metricPlist.isParameter("upper" ) && pass != false) {
      double max = metricPlist.get<double>("upper");
      if (value > max)
      {
        if (comm->getRank() == 0)
          msg << test_name << " FAILED: imbalance per part, "
              << value << ", greater than specified allowable maximum, " 
              << max << ".\n";
        pass = false;
      }else{
        if (comm->getRank() == 0)
          msg << test_name << " PASSED: imbalance per part, "
              << value << ", less than specified allowable maximum, " 
              << max << ".\n";
      }
      
    }
    
    return pass;
  }
};

#endif //ZOLTAN2_METRIC_ANALYZER_HPP
