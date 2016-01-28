//*****************************************************************//
//    Albany 2.0:  Copyright 2012 Sandia Corporation               //
//    This Software is released under the BSD license detailed     //
//    in the file "license.txt" in the top-level Albany directory  //
//*****************************************************************//

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
class Function_Base
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

};

///
/// Minimizer Struct
///
template<typename T, Index N>
struct Minimizer
{
public:
  Minimizer():
  max_num_iter(256),
  min_num_iter(0),
  rel_tol(1.0e-12),
  rel_error(1.0),
  abs_tol(1.0e-12),
  abs_error(1.0),
  converged(false),
  initial_norm(1.0),
  final_value(0.0),
  num_iter(0),
  step_method_name(nullptr),
  function_name(nullptr)
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

  bool
  converged{false};

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

  T
  hessian_singular_tol{1.0e-12};
};

} // namespace Intrepid2

#include "Intrepid2_MiniTensor_Solvers.t.h"

#endif // Intrepid2_MiniTensor_Solvers_h
