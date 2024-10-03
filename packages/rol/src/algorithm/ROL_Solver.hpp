// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_SOLVER_HPP
#define ROL_SOLVER_HPP

#include "ROL_TypeU_AlgorithmFactory.hpp"
#include "ROL_TypeB_AlgorithmFactory.hpp"
#include "ROL_TypeE_AlgorithmFactory.hpp"
#include "ROL_TypeG_AlgorithmFactory.hpp"
#include "ROL_Problem.hpp"

#include "ROL_ParameterList.hpp"
#include "ROL_Stream.hpp"

/** \class ROL::Solver
    \brief Provides a simplified interface for solving a wide range of
           optimization problems
 */

namespace ROL {

template<typename Real>
class Solver {
private:

  const Ptr<Problem<Real>> opt_;
  const EProblem problemType_;

  Ptr<TypeU::Algorithm<Real>> algoU_;
  Ptr<TypeB::Algorithm<Real>> algoB_;
  Ptr<TypeE::Algorithm<Real>> algoE_;
  Ptr<TypeG::Algorithm<Real>> algoG_;

public:

  /** \brief Constructor.
  
       @param[in] opt       the OptimizationProblem to be solved
       @param[in] parlist   algorithm and step input parameters

      ---
  */
  Solver( const Ptr<Problem<Real>> &opt,
          ParameterList            &list,
          const Ptr<Secant<Real>>  &secant = nullPtr );

  /** \brief Solve optimization problem with no iteration output.

      @param[in] status          is a user-defined StatusTest
      @param[in] combineStatus   if true, the user-defined StatusTest will be combined with the default StatusTest

      ---
  */
  int solve( const Ptr<StatusTest<Real>> &status = nullPtr,
             bool combineStatus = true);

  /** \brief Solve optimization problem.

      @param[in] outStream       is the output stream to collect iteration history
      @param[in] status          is a user-defined StatusTest
      @param[in] combineStatus   if true, the user-defined StatusTest will be combined with the default StatusTest

      ---
  */
  int solve( std::ostream &outStream,
             const Ptr<StatusTest<Real>> &status = nullPtr,
             bool combineStatus = true );

  /** \brief Return the AlgorithmState.

      ---
  */
  //Ptr<const AlgorithmState<Real>>& getAlgorithmState() const;
  Ptr<const AlgorithmState<Real>> getAlgorithmState() const;

  /** \brief Reset both Algorithm and Step.

      This function will reset the AlgorithmState and reinitialize the
      Step.  This function does not permit changing the Step specified
      upon construction.  To change the Step, reinitialize the
      OptimizationSolver.

      ---
  */
  void reset();

}; // class Solver

} // namespace ROL

#include "ROL_Solver_Def.hpp"

#endif // ROL_SOLVER_HPP


