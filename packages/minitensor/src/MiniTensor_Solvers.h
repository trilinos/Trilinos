// @HEADER
// *****************************************************************************
//                           MiniTensor Package
//
// Copyright 2016 NTESS and the MiniTensor contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#if !defined(MiniTensor_Solvers_h)
#define MiniTensor_Solvers_h

#include <memory>
#include <utility>

#include <MiniTensor.h>

namespace minitensor
{
/// The Fad type to use.
template<typename T, int N>
using FAD = Sacado::Fad::SLFad<T, N>;

///
/// Function base class that defines the interface to Mini Solvers.
/// The dimensions M and N are different in case vectors are passed
/// to the methods that have different capacity. They must have
/// dimension N.
///
template<typename FunctionDerived, typename S, Index M>
struct Function_Base
{
public:
  static constexpr
  Index
  DIMENSION{M};

  Function_Base()
  {
  }

  ///
  /// By default use merit function 0.5 dot(residual,residual)
  /// as the target to optimize if only the residual is provided.
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
  /// Defined explicitly.
  ///
  template<typename T, Index N>
  Vector<T, N>
  residual(FunctionDerived & f, Vector<T, N> const & x);

  ///
  /// By default compute Hessian with AD from gradient().
  ///
  template<typename T, Index N>
  Tensor<T, N>
  hessian(FunctionDerived & f, Vector<T, N> const & x);

  void
  set_failed(char const * const msg = nullptr);

  bool
  get_failed();

  void
  clear_failed();

  void
  set_failure_message(char const * const msg = nullptr);

  char const *
  get_failure_message();

protected:
  ///
  /// Signal that something has gone horribly wrong.
  ///
  bool
  failed{false};

  ///
  /// Keep a message to inform what went wrong above.
  char const *
  failure_message{nullptr};  
};

///
/// Equality constraint base class that defines the interface to Mini Solvers.
///
template<typename ConstraintDerived, typename S, Index NC, Index NV>
struct Equality_Constraint
{
public:
  Equality_Constraint()
  {
  }

  ///
  ///
  template<typename T, Index N>
  Vector<T, NC>
  value(ConstraintDerived & c, Vector<T, N> const & x);

  ///
  /// By default compute gradient with AD from value().
  ///
  template<typename T, Index N>
  Matrix<T, NC, NV>
  gradient(ConstraintDerived & c, Vector<T, N> const & x);

  static constexpr
  bool
  IS_EQUALITY{true};

  ///
  /// Signal that something has gone horribly wrong.
  ///
  bool
  failed{false};

  static constexpr
  Index
  NUM_CONSTR{NC};

  static constexpr
  Index
  NUM_VAR{NV};
};

///
/// Inequality constraint base class that defines the interface to Mini Solvers.
///
template<typename ConstraintDerived, typename S, Index NC, Index NV>
struct Inequality_Constraint :
    public Equality_Constraint<ConstraintDerived, S, NC, NV>
{
  static constexpr
  bool
  IS_EQUALITY{false};
};

///
/// Bounds constraint
///
template<typename T, Index N>
struct Bounds
{
  Bounds(Vector<T, N> const & l, Vector<T, N> const & u);

  Vector<T, N>
  lower;

  Vector<T, N>
  upper;
};

///
/// Minimizer Struct
///
template<typename T, Index N>
struct Minimizer
{
public:
  Minimizer();

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
  recordFinals(FN & fn, Vector<T, N> const & x);

public:
  Index
  max_num_iter{256};

  Index
  min_num_iter{0};

  Index
  num_iter{0};

  Index
  num_stagnation_iter{0};

  Index
  max_stagnation_iter{0};

  T
  initial_norm{1.0};

  T
  rel_tol{1.0e-12};

  T
  rel_error{1.0};

  T
  abs_tol{1.0e-12};

  T
  acc_tol{1.0e-12};  

  T
  stagnation_tol{1.0};

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
  failed{false};

  bool
  warning{false};
  
  bool
  converged{false};

  bool
  monotonic{true};

  bool
  bounded{true};

  bool
  non_stagnant{true};

  bool
  enforce_monotonicity{false};

  bool
  enforce_boundedness{false};

  bool
  enforce_non_stagnation{false};

  Vector<T, N>
  initial_guess;

  Vector<T, N>
  final_soln;

  Vector<T, N>
  final_gradient;

  Tensor<T, N>
  final_hessian;

  char const *
  step_method_name{nullptr};

  char const *
  function_name{nullptr};

  char const *
  failure_message{"No failure detected"};

  char const *
  warning_message{"No warning detected"};  
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
/// Back-track line search
///
template<typename T, Index N>
struct BacktrackingLineSearch
{
  template<typename FN>
  Vector<T, N>
  step(FN & fn, Vector<T, N> const & direction, Vector<T, N> const & soln);

  Index
  max_num_iter{100};

  Index
  max_line_iter{10};

  T
  search_parameter{0.5};

  T
  search_increment{0.1};

  T
  alpha{1.0};

  T
  tolerance{1.0e-6};
};

///
/// Trust region subproblem base
///
template<typename T, Index N>
struct TrustRegionSubproblemBase
{
  PreconditionerType
  preconditioner_type{PreconditionerType::IDENTITY};

  Vector<T, N>
  lin_solve(Tensor<T, N> const & A, Vector<T, N> const & b);
};

///
/// Trust region subproblem with a given objective function.
/// Exact algorithm, Nocedal 2nd Ed 4.3
///
template<typename T, Index N>
struct TrustRegionExactValue : public TrustRegionSubproblemBase<T, N>
{
  Vector<T, N>
  step(Tensor<T, N> const & Hessian, Vector<T, N> const & gradient);

  Index
  max_num_iter{4};

  T
  region_size{1.0};
};

///
/// Trust region subproblem with a given gradient/residual.
/// Exact algorithm, Nocedal 2nd Ed 4.3
///
template<typename T, Index N>
struct TrustRegionExactGradient : public TrustRegionSubproblemBase<T, N>
{
  Vector<T, N>
  step(Tensor<T, N> const & Hessian, Vector<T, N> const & gradient);

  Index
  max_num_iter{4};

  T
  region_size{1.0};
}; 

///
/// Trust region subproblem with a given objective function.
/// Dog leg algorithm.
///
template<typename T, Index N>
struct TrustRegionDogLegValue : public TrustRegionSubproblemBase<T, N>
{
  Vector<T, N>
  step(Tensor<T, N> const & Hessian, Vector<T, N> const & gradient);

  T
  region_size{1.0};
};

///
/// Trust region subproblem with a given gradient/residual.
/// Dog leg algorithm.
///
template<typename T, Index N>
struct TrustRegionDogLegGradient : public TrustRegionSubproblemBase<T, N>
{
  Vector<T, N>
  step(Tensor<T, N> const & Hessian, Vector<T, N> const & gradient);

  T
  region_size{1.0};
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
  char const *
  name() = 0;

  virtual
  void
  initialize(FN & fn, Vector<T, N> const & x, Vector<T, N> const & r) = 0;

  virtual
  Vector<T, N>
  step(FN & fn, Vector<T, N> const & x, Vector<T, N> const & r) = 0;

  virtual
  ~StepBase() {}

  PreconditionerType
  preconditioner_type{PreconditionerType::IDENTITY};

  Vector<T, N>
  lin_solve(Tensor<T, N> const & A, Vector<T, N> const & b);
};

///
/// The step types
///
enum class StepType
{
  UNDEFINED = 0,
  NEWTON = 1,
  NEWTON_LS = 2,
  TRUST_REGION = 3,
  CG = 4,
  LINE_SEARCH_REG = 5
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
struct NewtonStep : public StepBase<FN, T, N>
{
  static constexpr
  char const * const
  NAME{"Newton"};

  virtual
  char const *
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
/// Newton Step with line search
///
template<typename FN, typename T, Index N>
struct NewtonWithLineSearchStep : public StepBase<FN, T, N>
{
  static constexpr
  char const * const
  NAME{"Newton with Line Search"};

  virtual
  char const *
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
  ~NewtonWithLineSearchStep() {}
};

///
/// Trust Region Step
///
template<typename FN, typename T, Index N>
struct TrustRegionStep : public StepBase<FN, T, N>
{
  static constexpr
  char const * const
  NAME{"Trust Region"};

  virtual
  char const *
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
struct ConjugateGradientStep : public StepBase<FN, T, N>
{
  static constexpr
  char const * const
  NAME{"Preconditioned Conjugate Gradient"};

  virtual
  char const *
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
struct LineSearchRegularizedStep : public StepBase<FN, T, N>
{
  static constexpr
  char const * const
  NAME{"Line Search Regularized"};

  virtual
  char const *
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

} // namespace minitensor

#include "MiniTensor_Solvers.t.h"

#endif // MiniTensor_Solvers_h
