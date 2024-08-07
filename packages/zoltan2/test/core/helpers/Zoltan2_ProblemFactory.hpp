// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file Zoltan2_TestFactory.hpp
    \brief Returns a pointer to new test classes.
    Is not responsible for memory management!
*/
#ifndef ZOLTAN2_PROBLEM_FACTORY_HPP
#define ZOLTAN2_PROBLEM_FACTORY_HPP
#include <Zoltan2_TestHelpers.hpp>
#include <Zoltan2_Problem.hpp>
#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_OrderingProblem.hpp>
#include <Zoltan2_ColoringProblem.hpp>
#include <Zoltan2_Typedefs.hpp>

using namespace Zoltan2_TestingFramework;
using namespace Zoltan2;

namespace Zoltan2_TestingFramework {
/// \brief ProblemFactory class contains 1 static factory method
  class ProblemFactory {
  public:
    /// \brif Zoltan2::Problem factory method
    ///
    /// @param kind A string equal to the type of the problem (paritioning, ordering, coloring)
    /// @param kind A string equal to the type of the adapter
    /// @param params Zolta2 parameter list
    /// @param (MPI) MPI world communicator
    ///
    /// @return returns a pointer to new Zoltan2::Problem or a nullptr if
    /// problem_name was not known.
    ProblemFactory(const std::string & problemName,
                   RCP<AdapterFactory> adapterFactory,
                   ParameterList *params
                   #ifdef HAVE_ZOLTAN2_MPI
                    , MPI_Comm comm
                   #endif
                   ) {

      problem_name = problemName;
      adapterType = adapterFactory->getMainAdapterType();

      #ifdef HAVE_ZOLTAN2_MPI
        #define CREATE_PRBLM(problemClass, adapterClass)                       \
          adapterClass * pCast = dynamic_cast<adapterClass *>                  \
            (adapterFactory->getMainAdapter());                                \
          if(!pCast) { throw std::logic_error(                                 \
            "ProblemFactory adapter dynamic_cast failed for problem name "     \
              + problem_name + " and adapterClass " + #adapterClass ); }       \
          problem = rcp(new problemClass<adapterClass>(pCast, params, comm));
      #else
        #define CREATE_PRBLM(problemClass, adapterClass)                       \
          adapterClass * pCast = dynamic_cast<adapterClass *>                  \
            (adapterFactory->getMainAdapter());                                \
          if(!pCast) { throw std::logic_error(                                 \
            "ProblemFactory adapter dynamic_cast failed for problem name "     \
              + problem_name + " and adapterClass " + #adapterClass ); }       \
          problem = rcp(new problemClass<adapterClass>(pCast, params));
      #endif

      #define MAKE_PARTITION_PROBLEM(adapterClass)  \
        CREATE_PRBLM(PartitioningProblem, adapterClass);

      #define MAKE_ORDERING_PROBLEM(adapterClass)  \
        CREATE_PRBLM(OrderingProblem, adapterClass);

      // PartitioningProblem
      if(problem_name == "partitioning") {
        Z2_TEST_UPCAST(adapterType, MAKE_PARTITION_PROBLEM)
      }
      else if(problem_name == "ordering") {
        Z2_TEST_UPCAST(adapterType, MAKE_ORDERING_PROBLEM)
      }

      if(problem == Teuchos::null) {
        throw std::logic_error("ProblemFactory failed to create Problem!");
      }
    }

    RCP<ProblemRoot> getProblem() { return problem; }
    const std::string & getProblemName() const { return problem_name; }
    EAdapterType getAdapterType() const { return adapterType; }

    private:
      std::string problem_name; // string converts to a problem type
      EAdapterType adapterType; // converts to an adapter type
      RCP<ProblemRoot> problem;
  };
}
#endif // ZOLTAN2_PROBLEM_FACTORY_HPP

