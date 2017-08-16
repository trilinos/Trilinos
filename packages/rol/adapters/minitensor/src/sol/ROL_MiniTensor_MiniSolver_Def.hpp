// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions: Alejandro Mota (amota@sandia.gov)
//
// ************************************************************************
// @HEADER
#include <type_traits>

#include "ROL_Algorithm.hpp"
#include "ROL_MiniTensor_BoundConstraint.hpp"
#include "ROL_MiniTensor_EqualityConstraint.hpp"
#include "ROL_MiniTensor_Function.hpp"
#include "ROL_MiniTensor_InequalityConstraint.hpp"
#include "ROL_MiniTensor_Vector.hpp"

namespace ROL {

using Index = minitensor::Index;

//
//
//
template<typename T, Index N>
MiniTensor_Minimizer<T, N>::
MiniTensor_Minimizer()
{
  constexpr bool
  is_fad = Sacado::IsADType<T>::value == true;

  static_assert(is_fad == false, "AD types not allowed for type T");

  return;
}

//
//
//
template<typename T, Index N>
template<typename FN>
void
MiniTensor_Minimizer<T, N>::
solve(
    std::string const & algoname,
    Teuchos::ParameterList & params,
    FN & fn,
    minitensor::Vector<T, N> & soln)
{
  step_method_name = algoname.c_str();
  function_name = FN::NAME;
  initial_guess = soln;

  minitensor::Vector<T, N>
  resi = fn.gradient(soln);

  initial_value = fn.value(soln);
  previous_value = initial_value;
  failed = failed || fn.get_failed();
  if (fn.get_failed() == true) failure_message = fn.get_failure_message();
  initial_norm = minitensor::norm(resi);

  // Define algorithm.
  ROL::MiniTensor_Objective<FN, T, N>
  obj(fn);

  ROL::Algorithm<T>
  algo(algoname, params);

  // Set Initial Guess
  ROL::MiniTensorVector<T, N>
  x(soln);

  // Run Algorithm
  algo.run(x, obj, verbose);

  soln = ROL::MTfromROL<T, N>(x);

  resi = fn.gradient(soln);

  T const
  norm_resi = minitensor::norm(resi);

  updateConvergenceCriterion(norm_resi);

  ROL::AlgorithmState<T> &
  state = const_cast<ROL::AlgorithmState<T> &>(*(algo.getState()));

  ROL::StatusTest<T>
  status(params);

  converged = status.check(state) == false;
  recordFinals(fn, soln);

  return;
}

//
//
//
template<typename T, Index N>
template<typename FN, typename BC>
void
MiniTensor_Minimizer<T, N>::
solve(
    std::string const & algoname,
    Teuchos::ParameterList & params,
    FN & fn,
    BC & bc,
    minitensor::Vector<T, N> & soln)
{
  step_method_name = algoname.c_str();
  function_name = FN::NAME;
  initial_guess = soln;

  minitensor::Vector<T, N>
  resi = fn.gradient(soln);

  initial_value = fn.value(soln);
  previous_value = initial_value;
  failed = failed || fn.get_failed();
  if (fn.get_failed() == true) failure_message = fn.get_failure_message();
  initial_norm = minitensor::norm(resi);

  // Define algorithm.
  ROL::MiniTensor_Objective<FN, T, N>
  obj(fn);

  MiniTensorVector<T, N>
  lo(bc.lower);

  MiniTensorVector<T, N>
  hi(bc.upper);

  // Define bound constraint
  ROL::MiniTensor_BoundConstraint<T, N>
  bound_constr(lo, hi);

  ROL::Algorithm<T>
  algo(algoname, params);

  // Set initial guess
  ROL::MiniTensorVector<T, N>
  x(soln);

  // Run Algorithm
  algo.run(x, obj, bound_constr, verbose);

  soln = ROL::MTfromROL<T, N>(x);

  resi = fn.gradient(soln);

  T const
  norm_resi = minitensor::norm(resi);

  updateConvergenceCriterion(norm_resi);

  ROL::AlgorithmState<T> &
  state = const_cast<ROL::AlgorithmState<T> &>(*(algo.getState()));

  ROL::StatusTest<T>
  status(params);

  converged = status.check(state) == false;
  recordFinals(fn, soln);

  return;
}

//
//
//
template<typename T, Index N>
template<typename FN, typename EIC, Index NC>
void
MiniTensor_Minimizer<T, N>::
solve(
    std::string const & algoname,
    Teuchos::ParameterList & params,
    FN & fn,
    EIC & eic,
    minitensor::Vector<T, N> & soln,
    minitensor::Vector<T, NC> & cv)
{
  step_method_name = algoname.c_str();
  function_name = FN::NAME;
  initial_guess = soln;

  minitensor::Vector<T, N>
  resi = fn.gradient(soln);

  initial_value = fn.value(soln);
  previous_value = initial_value;
  failed = failed || fn.get_failed();
  if (fn.get_failed() == true) failure_message = fn.get_failure_message();
  initial_norm = minitensor::norm(resi);

  // Define algorithm.
  ROL::MiniTensor_Objective<FN, T, N>
  obj(fn);

  // Define constraint
  using MTCONSTR = typename std::conditional<EIC::IS_EQUALITY,
      ROL::MiniTensor_EqualityConstraint<EIC, T, NC, N>,
      ROL::MiniTensor_InequalityConstraint<EIC, T, NC, N>>::type;

  MTCONSTR
  constr(eic);

  ROL::Algorithm<T>
  algo(algoname, params);

  // Set Initial Guess
  ROL::MiniTensorVector<T, N>
  x(soln);

  // Set constraint vector
  ROL::MiniTensorVector<T, NC>
  c(cv);

  // Run Algorithm
  algo.run(x, c, obj, constr, verbose);

  soln = ROL::MTfromROL<T, N>(x);

  resi = fn.gradient(soln);

  T const
  norm_resi = minitensor::norm(resi);

  updateConvergenceCriterion(norm_resi);

  ROL::AlgorithmState<T> &
  state = const_cast<ROL::AlgorithmState<T> &>(*(algo.getState()));

  ROL::StatusTest<T>
  status(params);

  converged = status.check(state) == false;

  recordFinals(fn, soln);

  return;
}

//
//
//
template<typename T, Index N>
template<typename FN, typename EIC, typename BC, Index NC>
void
MiniTensor_Minimizer<T, N>::
solve(
    std::string const & algoname,
    Teuchos::ParameterList & params,
    FN & fn,
    EIC & eic,
    BC & bc,
    minitensor::Vector<T, N> & soln,
    minitensor::Vector<T, NC> & cv)
{
  step_method_name = algoname.c_str();
  function_name = FN::NAME;
  initial_guess = soln;

  minitensor::Vector<T, N>
  resi = fn.gradient(soln);

  initial_value = fn.value(soln);
  previous_value = initial_value;
  failed = failed || fn.get_failed();
  if (fn.get_failed() == true) failure_message = fn.get_failure_message();
  initial_norm = minitensor::norm(resi);

  // Define algorithm.
  ROL::MiniTensor_Objective<FN, T, N>
  obj(fn);

  MiniTensorVector<T, N>
  lo(bc.lower);

  MiniTensorVector<T, N>
  hi(bc.upper);

  // Define bound constraint
  ROL::MiniTensor_BoundConstraint<T, N>
  bound_constr(lo, hi);

  // Define equality/inequality constraint
  using MTCONSTR = typename std::conditional<EIC::IS_EQUALITY,
      ROL::MiniTensor_EqualityConstraint<EIC, T, NC, N>,
      ROL::MiniTensor_InequalityConstraint<EIC, T, NC, N>>::type;

  MTCONSTR
  eqineq_constr(eic);

  ROL::Algorithm<T>
  algo(algoname, params);

  // Set Initial Guess
  ROL::MiniTensorVector<T, N>
  x(soln);

  // Set constraint vector
  ROL::MiniTensorVector<T, NC>
  c(cv);

  // Run Algorithm
  algo.run(x, c, obj, eqineq_constr, bound_constr, verbose);

  soln = ROL::MTfromROL<T, N>(x);

  resi = fn.gradient(soln);

  T const
  norm_resi = minitensor::norm(resi);

  updateConvergenceCriterion(norm_resi);

  ROL::AlgorithmState<T> &
  state = const_cast<ROL::AlgorithmState<T> &>(*(algo.getState()));

  ROL::StatusTest<T>
  status(params);

  converged = status.check(state) == false;

  recordFinals(fn, soln);

  return;
}
//
//
//
template<typename T, Index N>
void
MiniTensor_Minimizer<T, N>::
printReport(std::ostream & os)
{
  char const * const
  converged_string = converged == true ? "YES" : "NO";

  // Happy / frowny face
  //char const * const
  //converged_string = converged == true ? "\U0001F60A" : "\U0001F623";

  os << "\n\n";
  os << "Method       : " << step_method_name << '\n';
  os << "Function     : " << function_name << '\n';
  os << "Converged    : " << converged_string << '\n';

  os << std::scientific << std::setprecision(16);

  os << "Initial |R|  : " << std::setw(24) << initial_norm << '\n';
  os << "Abs Tol      : " << std::setw(24) << abs_tol << '\n';
  os << "Abs Error    : " << std::setw(24) << abs_error << '\n';
  os << "Rel Tol      : " << std::setw(24) << rel_tol << '\n';
  os << "Rel Error    : " << std::setw(24) << rel_error << '\n';
  os << "Initial X    : " << initial_guess << '\n';
  os << "Initial f(X) : " << std::setw(24) << initial_value << '\n';
  os << "Final X      : " << final_soln << '\n';
  os << "FInal f(X)   : " << std::setw(24) << final_value << '\n';
  os << "Final Df(X)  : " << final_gradient << '\n';
  os << "FInal DDf(X) : " << final_hessian << '\n';
  os << '\n';

  return;
}

//
//
//
template<typename T, Index N>
void
MiniTensor_Minimizer<T, N>::
updateConvergenceCriterion(T const ae)
{
  abs_error = ae;
  rel_error = initial_norm > 0.0 ? abs_error / initial_norm : T(0.0);

  bool const
  converged_absolute = abs_error <= abs_tol;

  bool const
  converged_relative = rel_error <= rel_tol;

  converged = converged_absolute || converged_relative;

  return;
}

//
//
//
template<typename T, Index N>
void
MiniTensor_Minimizer<T, N>::
updateDivergenceCriterion(T const fn_value)
{
  monotonic = fn_value <= previous_value;

  if (enforce_monotonicity == true && monotonic == false) {
    failed = true;
  }

  previous_value = fn_value;

  bounded = fn_value <= growth_limit * initial_value;

  if (enforce_boundedness == true && bounded == false) {
    failed = true;
  }

  return;
}

//
//
//
template<typename T, Index N>
template<typename FN>
void
MiniTensor_Minimizer<T, N>::
recordFinals(FN & fn, minitensor::Vector<T, N> const & x)
{
  final_soln = x;
  final_value = fn.value(x);
  final_gradient = fn.gradient(x);
  final_hessian = fn.hessian(x);
}

} // namespace ROL
