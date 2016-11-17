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

/*! \file Zoltan2_EvaluateFactory.hpp
    \brief Returns a pointer to new test classes.
    Is not responsible for memory management!
*/
#ifndef ZOLTAN2_EVALUATE_FACTORY_HPP
#define ZOLTAN2_EVALUATE_FACTORY_HPP

#include <Zoltan2_Typedefs.hpp>
#include <Zoltan2_EvaluatePartition.hpp>
#include <Zoltan2_EvaluateOrdering.hpp>
#include <Zoltan2_OrderingProblem.hpp>
#include <Zoltan2_ProblemFactory.hpp>

using namespace Zoltan2_TestingFramework;
using namespace Zoltan2;

namespace Zoltan2_TestingFramework {
/// \brief ProblemFactory class contains 1 static factory method
  class EvaluateFactory {
  public:

    /// \brif Zoltan2::EvaluateBaseClass factory method
    ///
    /// @param kind A string equal to the type of the problem (paritioning, ordering, coloring)
    /// @param kind A string equal to the type of the adapter
    /// @param input Input adapter used to construct the problem
    /// @param params Zolta2 parameter list
    /// @param (MPI) MPI world communicator
    ///
    /// @return returns a pointer to new Zoltan2::Problem or a nullptr if kind was not known.

    EvaluateFactory(const std::string &problem_name,
                    RCP<AdapterFactory> adapterFactory,
                    ParameterList *params,
                    RCP<ProblemFactory> problem) {

      adapter_template_name = adapterFactory->adaptersSet.main.template_name;
      this->problem_name = problem_name;

      if (problem_name == "partitioning") {
        #define PARTITIONING_PROBLEM(adapterClass) rcp_dynamic_cast<           \
          PartitioningProblem<adapterClass>> (problem->getProblem())

        #define EVALUATE_PARTITION(adapterClass)                               \
        if (adapter_template_name == #adapterClass) {                          \
          const adapterClass * pAdapterClassUpCast = dynamic_cast<             \
            const adapterClass *>(adapterFactory->adaptersSet.main.adapter);   \
          if(!pAdapterClassUpCast) throw std::logic_error( "Bad adapter class cast!"  );   \
          evaluate = rcp(new EvaluatePartition<adapterClass>(                  \
             pAdapterClassUpCast,                                              \
              params, problem->getComm(),                                      \
              (&PARTITIONING_PROBLEM(adapterClass)->getSolution())));          \
        }

        TEMPLATE_CONVERSION(EVALUATE_PARTITION)
      }
      else if(problem_name == "ordering") {
        #define ORDERING_PROBLEM(adapterClass) rcp_dynamic_cast<               \
          OrderingProblem<adapterClass>> (problem->getProblem())

        #define LOCAL_ORDERING(adapterClass)                                   \
        if (adapter_template_name == #adapterClass) {                          \
          const adapterClass * pAdapterClassUpCast = dynamic_cast<             \
            const adapterClass *>(adapterFactory->adaptersSet.main.adapter);   \
          if(!pAdapterClassUpCast) throw std::logic_error( "Bad adapter class cast!"  );   \
          evaluate = rcp(new EvaluateLocalOrdering<adapterClass>(              \
              pAdapterClassUpCast,                                             \
              params, problem->getComm(),                                      \
              ORDERING_PROBLEM(adapterClass)->getLocalOrderingSolution()));    \
        }

        // EvaluateGlobalOrdering not tested/implemented yet
        #define GLOBAL_ORDERING(adapterClass)                                  \
        if (adapter_template_name == #adapterClass) {                          \
          const adapterClass * pAdapterClassUpCast = dynamic_cast<             \
            const adapterClass *>(adapterFactory->adaptersSet.main);           \
          if(!pAdapterClassUpCast) throw std::logic_error( "Bad adapter class cast!"  );   \
          evaluate = rcp(new EvaluateGlobalOrdering<adapterClass>(             \
              pAdapterClassUpCast,                                             \
              params, ORDERING_PROBLEM(adapterClass)->getComm(),               \
              ORDERING_PROBLEM(adapterClass)->getGlobalOrderingSolution()));   \
        }

        TEMPLATE_CONVERSION(LOCAL_ORDERING)
      }
      else if(problem_name == "coloring") {
        // Coloring code here... EvaluateColoringFactory not created yet
        // return EvaluateColoringFactory::newEvaluatColoring(
        //           dynamic_cast<coloring_problem_t*> (problem),
        //           adapter_name, input, params);
      }
    }

    void printMetrics(std::ostringstream & msg) {

      #define PRINT_METRICS(adapterClass)                          \
        if (adapter_template_name == #adapterClass) {                                      \
          RCP<EvaluateBaseClass<adapterClass>> pCast = rcp_dynamic_cast<EvaluateBaseClass<adapterClass>>(evaluate); \
          if(pCast == Teuchos::null) throw std::logic_error( "Bad evaluate class cast in printMetrics!"  );   \
          pCast->printMetrics(msg);                                         \
        }
        TEMPLATE_CONVERSION(PRINT_METRICS)
    }

    bool analyzeMetrics(std::ostringstream & msg,
      const ParameterList &problem_parameters) {

      #define ANALYZE_METRICS(adapterClass)                        \
        if (adapter_template_name == #adapterClass) {                                      \
          RCP<EvaluateBaseClass<adapterClass>> pCast = rcp_dynamic_cast<EvaluateBaseClass<adapterClass>>(evaluate); \
          if(pCast == Teuchos::null) throw std::logic_error( "Bad evaluate class cast in analyzeMetrics!"  );   \
          return MetricAnalyzer<adapterClass>::analyzeMetrics(pCast, problem_parameters.sublist("Metrics"), msg);                                         \
        }
        TEMPLATE_CONVERSION(ANALYZE_METRICS)
        return false;
    }

    void loadMetricInfo(std::vector<MetricAnalyzerInfo<zscalar_t>> & metricInfo,
      const ParameterList &metricsPlist) {

      #define LOAD_METRIC_INFO(adapterClass)                       \
        if (adapter_template_name == #adapterClass) {                                      \
        RCP<EvaluateBaseClass<adapterClass>> pCast = rcp_dynamic_cast<EvaluateBaseClass<adapterClass>>(evaluate); \
        if(pCast == Teuchos::null) throw std::logic_error( "Bad evaluate class cast in loadMetricInfo!"  );   \
          MetricAnalyzer<adapterClass>::LoadMetricInfo(metricInfo,             \
          pCast,          \
          metricsPlist.sublist("Metrics"));                                    \
        }
        TEMPLATE_CONVERSION(LOAD_METRIC_INFO)
    }

    bool isValid() const {
      return ( evaluate != Teuchos::null );
    }

    private:
      std::string problem_name;  // string converts to a problem type
      std::string adapter_template_name;  // string converts to an adapter type
      RCP<EvaluateBaseClassRoot> evaluate;
  };
}
#endif // ZOLTAN2_EVALUATE_FACTORY_HPP

