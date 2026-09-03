// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#if !defined(ROL_MiniTensor_MiniSolver_hpp)
#define ROL_MiniTensor_MiniSolver_hpp

#include "MiniTensor_Solvers.h"

namespace ROL {

using Index = minitensor::Index;

///
/// Minimizer Struct
///
template<typename T, Index N>
struct MiniTensor_Minimizer
{
public:
  MiniTensor_Minimizer();

  template<typename FN>
  void
  solve(
      std::string const & algoname,
      Teuchos::ParameterList & params,
      FN & fn,
      minitensor::Vector<T, N> & soln);

  template<typename FN, typename BC>
  void
  solve(
      std::string const & algoname,
      Teuchos::ParameterList & params,
      FN & fn,
      BC & bc,
      minitensor::Vector<T, N> & soln);

  template<typename FN, typename EIC, Index NC>
  void
  solve(
      std::string const & algoname,
      Teuchos::ParameterList & params,
      FN & fn,
      EIC & eic,
      minitensor::Vector<T, N> & soln,
      minitensor::Vector<T, NC> & cv);

  template<typename FN, typename EIC, typename BC, Index NC>
  void
  solve(
      std::string const & algoname,
      Teuchos::ParameterList & params,
      FN & fn,
      EIC & eic,
      BC & bc,
      minitensor::Vector<T, N> & soln,
      minitensor::Vector<T, NC> & cv);

  void
  printReport(std::ostream & os);

private:
  void
  updateConvergenceCriterion(T const abs_error);

  void
  updateDivergenceCriterion(T const fn_value);

  template<typename FN>
  void
  recordFinals(FN & fn, minitensor::Vector<T, N> const & x);

public:
  T
  initial_norm{1.0};

  T
  rel_tol{1.0e-12};

  T
  rel_error{1.0};

  T
  abs_tol{1.0e-12};

  T
  abs_error{1.0};

  T
  growth_limit{1.0};

  T
  initial_value{0.0};

  T
  previous_value{0.0};

  T
  final_value{0.0};

  bool
  verbose{false};

  bool
  failed{false};

  bool
  converged{false};

  bool
  monotonic{true};

  bool
  bounded{true};

  bool
  enforce_monotonicity{false};

  bool
  enforce_boundedness{false};

  minitensor::Vector<T, N>
  initial_guess;

  minitensor::Vector<T, N>
  final_soln;

  minitensor::Vector<T, N>
  final_gradient;

  minitensor::Tensor<T, N>
  final_hessian;

  char const *
  step_method_name{nullptr};

  char const *
  function_name{nullptr};

  char const *
  failure_message{"No failure detected"};
};

} // namespace ROL

#include "ROL_MiniTensor_MiniSolver_Def.hpp"

#endif // ROL_MiniTensor_MiniSolver_hpp
