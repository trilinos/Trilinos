// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
    /// @param params Zoltan2 parameter list
    /// @param (MPI) MPI world communicator
    ///
    /// @return returns a pointer to new Zoltan2::Problem or a nullptr if kind was not known.

    EvaluateFactory(const std::string & problemName,
                    RCP<AdapterFactory> adapterFactory,
                    ParameterList *params,
                    RCP<ProblemFactory> problemFactory) {

      adapterType = adapterFactory->getMainAdapterType();
      problem_name = problemName;

      if (problem_name == "partitioning") {
        #define PARTITIONING_PROBLEM(adapterClass) rcp_dynamic_cast<           \
          PartitioningProblem<adapterClass>> (problemFactory->getProblem())

        #define EVALUATE_PARTITION(adapterClass)                               \
          const adapterClass * pAdapterClassUpCast = dynamic_cast<             \
            const adapterClass *>(adapterFactory->getMainAdapter());   \
          if(!pAdapterClassUpCast) throw std::logic_error(                     \
            "Bad adapter class cast!"  );                                      \
          evaluate = rcp(new EvaluatePartition<adapterClass>(                  \
             pAdapterClassUpCast, params,                                      \
             problemFactory->getProblem()->getComm(),                          \
              (&PARTITIONING_PROBLEM(adapterClass)->getSolution())));

        Z2_TEST_UPCAST(adapterType, EVALUATE_PARTITION)
      }
      else if(problem_name == "ordering") {
        #define ORDERING_PROBLEM(adapterClass) rcp_dynamic_cast<               \
          OrderingProblem<adapterClass>> (problemFactory->getProblem())

        #define LOCAL_ORDERING(adapterClass)                                   \
          const adapterClass * pAdapterClassUpCast = dynamic_cast<             \
            const adapterClass *>(adapterFactory->getMainAdapter());           \
          if(!pAdapterClassUpCast) throw std::logic_error(                     \
            "Bad adapter class cast!");                                        \
          evaluate = rcp(new EvaluateLocalOrdering<adapterClass>(              \
              pAdapterClassUpCast, params,                                     \
              problemFactory->getProblem()->getComm(),                         \
              ORDERING_PROBLEM(adapterClass)->getLocalOrderingSolution()));

        // EvaluateGlobalOrdering not tested/implemented yet
        #define GLOBAL_ORDERING(adapterClass)                                  \
          const adapterClass * pAdapterClassUpCast = dynamic_cast<             \
            const adapterClass *>(adapterFactory->getMainAdapter());           \
          if(!pAdapterClassUpCast) throw std::logic_error(                     \
            "Bad adapter class cast!"  );                                      \
          evaluate = rcp(new EvaluateGlobalOrdering<adapterClass>(             \
              pAdapterClassUpCast,                                             \
              params, ORDERING_PROBLEM(adapterClass)->getComm(),               \
              ORDERING_PROBLEM(adapterClass)->getGlobalOrderingSolution()));

        Z2_TEST_UPCAST(adapterType, LOCAL_ORDERING)
      }
      else if(problem_name == "coloring") {
        // Coloring code here... EvaluateColoringFactory not created yet
        // return EvaluateColoringFactory::newEvaluatColoring(
        //           dynamic_cast<coloring_problem_t*> (problem),
        //           adapter_name, input, params);
      }

      if(evaluate == Teuchos::null) {
        throw std::logic_error("EvaluateFactory failed to create!");
      }
    }

    RCP<EvaluateBaseClassRoot>  getEvaluateClass() { return evaluate; }
    const std::string & getProblemName() const { return problem_name; }
    EAdapterType getAdapterType() const { return adapterType; }

    private:
      std::string problem_name;  // string converts to a problem type
      EAdapterType adapterType;  // converts to an adapter type
      RCP<EvaluateBaseClassRoot> evaluate;
  };
}
#endif // ZOLTAN2_EVALUATE_FACTORY_HPP

