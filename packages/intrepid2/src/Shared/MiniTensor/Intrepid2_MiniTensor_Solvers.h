// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
//                 Copyright (2007) Sandia Corporation
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

#if !defined(Intrepid2_MiniTensor_Solvers_h)
#define Intrepid2_MiniTensor_Solvers_h

#include <memory>
#include <utility>

#include <Intrepid2_MiniTensor.h>

namespace Intrepid2
{
/// The Fad type to use.
template<typename T, int N>
using FAD = Sacado::Fad::SLFad<T, N>;

///
/// Function base class that defines the interface to Mini Solvers.
///
template<typename FunctionDerived, typename S>
struct Function_Base
{
public:
  Function_Base()
  {
  }

  ///
  /// By default use merit function 0.5 dot(gradient, gradient)
  /// as the target to optimize if only the gradient is provided.
  ///
  template<typename T, Index N>
  T
  value(FunctionDerived & f, Vector<T, N> const & x);

  ///
  /// By default compute gradient with AD from value().
  ///
  template<typename T, Index N>
  Vector<T, N>
  gradient(FunctionDerived & f, Vector<T, N> const & x);

  ///
  /// By default compute Hessian with AD from gradient().
  ///
  template<typename T, Index N>
  Tensor<T, N>
  hessian(FunctionDerived & f, Vector<T, N> const & x);

  ///
  /// Signal that something has gone horribly wrong.
  ///
  bool
  failed{false};

};

///
/// Minimizer Struct
///
template<typename T, Index N>
struct Minimizer
{
public:
  Minimizer()
  {
    constexpr bool
    is_fad = Sacado::IsADType<T>::value == true;

    static_assert(is_fad == false, "AD types not allowed for type T");
  }

  template<typename STEP, typename FN>
  void
  solve(STEP & step_method, FN & fn, Vector<T, N> & x);

  void
  printReport(std::ostream & os);

private:
  void
  updateConvergenceCriterion(T const abs_error);

  void
  updateDivergenceCriterion(T const fn_value);

  bool
  continueSolve() const;

  template<typename FN>
  void
  recordFinals(FN & fn, Vector<T, N> const & x)
  {
    final_soln = x;
    final_value = fn.value(x);
    final_gradient = fn.gradient(x);
    final_hessian = fn.hessian(x);
  }

public:
  Index
  max_num_iter{256};

  Index
  min_num_iter{0};

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

  T
  initial_value{0.0};

  T
  previous_value{0.0};

  T
  final_value{0.0};

  Vector<T, N>
  final_gradient;

  Tensor<T, N>
  final_hessian;

private:
  T
  initial_norm{1.0};

  Index
  num_iter{0};

  Vector<T, N>
  initial_guess;

  Vector<T, N>
  final_soln;

  char const *
  step_method_name{nullptr};

  char const *
  function_name{nullptr};
};

///
/// Newton line search
///
template<typename T, Index N>
struct NewtonLineSearch
{
  template<typename FN>
  Vector<T, N>
  step(FN & fn, Vector<T, N> const & direction, Vector<T, N> const & soln);

  Index
  max_num_iter{16};

  T
  tolerance{1.0e-6};
};

///
/// Trust region subproblem. Exact algorithm, Nocedal 2nd Ed 4.3
///
template<typename T, Index N>
struct TrustRegionExact
{
  Vector<T, N>
  step(Tensor<T, N> const & Hessian, Vector<T, N> const & gradient);

  Index
  max_num_iter{4};

  T
  region_size{1.0};

  T
  initial_lambda{0.0};
};

///
/// Step Base
///
template<typename FN, typename T, Index N>
struct StepBase
{
  StepBase()
  {
    constexpr bool
    is_fad = Sacado::IsADType<T>::value == true;

    static_assert(is_fad == false, "AD types not allowed for type T");
  }

  virtual
  char const * const
  name() = 0;

  virtual
  void
  initialize(FN & fn, Vector<T, N> const & x, Vector<T, N> const & r) = 0;

  virtual
  Vector<T, N>
  step(FN & fn, Vector<T, N> const & x, Vector<T, N> const & r) = 0;

  virtual
  ~StepBase() {}
};

///
/// The step types
///
enum class StepType
{
  UNDEFINED = 0, NEWTON = 1, TRUST_REGION = 2, CG = 3, LINE_SEARCH_REG = 4
};

///
///
///
template<typename FN, typename T, Index N>
std::unique_ptr<StepBase<FN, T, N>>
stepFactory(StepType step_type);

///
/// Plain Newton Step
///
template<typename FN, typename T, Index N>
struct NewtonStep final : public StepBase<FN, T, N>
{
  static constexpr
  char const * const
  NAME{"Newton"};

  virtual
  char const * const
  name()
  {
    return NAME;
  }

  virtual
  void
  initialize(FN & fn, Vector<T, N> const & x, Vector<T, N> const & r);

  virtual
  Vector<T, N>
  step(FN & fn, Vector<T, N> const & x, Vector<T, N> const & r);

  virtual
  ~NewtonStep() {}
};

///
/// Trust Region Step
///
template<typename FN, typename T, Index N>
struct TrustRegionStep final : public StepBase<FN, T, N>
{
  static constexpr
  char const * const
  NAME{"Trust Region"};

  virtual
  char const * const
  name()
  {
    return NAME;
  }

  virtual
  void
  initialize(FN & fn, Vector<T, N> const & x, Vector<T, N> const & r);

  virtual
  Vector<T, N>
  step(FN & fn, Vector<T, N> const & x, Vector<T, N> const & r);

  virtual
  ~TrustRegionStep() {}

  T
  max_region_size{10.0};

  T
  initial_region_size{10.0};

  T
  min_reduction{0.0};

private:
  T
  region_size{0.0};
};

///
/// Conjugate Gradient Step
///
template<typename FN, typename T, Index N>
struct ConjugateGradientStep final : public StepBase<FN, T, N>
{
  static constexpr
  char const * const
  NAME{"Preconditioned Conjugate Gradient"};

  virtual
  char const * const
  name()
  {
    return NAME;
  }

  virtual
  void
  initialize(FN & fn, Vector<T, N> const & x, Vector<T, N> const & r);

  virtual
  Vector<T, N>
  step(FN & fn, Vector<T, N> const & x, Vector<T, N> const & r);

  virtual
  ~ConjugateGradientStep() {}

  Index
  restart_directions_interval{32};

private:
  Vector<T, N>
  search_direction;

  Vector<T, N>
  precon_resi;

  T
  projection_new{0.0};

  Index
  restart_directions_counter{0};
};

///
/// Line Search Regularized Step
///
template<typename FN, typename T, Index N>
struct LineSearchRegularizedStep final : public StepBase<FN, T, N>
{
  static constexpr
  char const * const
  NAME{"Line Search Regularized"};

  virtual
  char const * const
  name()
  {
    return NAME;
  }

  virtual
  void
  initialize(FN & fn, Vector<T, N> const & x, Vector<T, N> const & r);

  virtual
  Vector<T, N>
  step(FN & fn, Vector<T, N> const & x, Vector<T, N> const & r);

  virtual
  ~LineSearchRegularizedStep() {}

  T
  step_length{1.0};

  T
  hessian_cond_tol{1.0e+08};
};

} // namespace Intrepid2

#include "Intrepid2_MiniTensor_Solvers.t.h"

#endif // Intrepid2_MiniTensor_Solvers_h
