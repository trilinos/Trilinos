// @HEADER
// *****************************************************************************
//                           MiniTensor Package
//
// Copyright 2016 NTESS and the MiniTensor contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

namespace minitensor
{

//
//
//
template<typename FunctionDerived, typename S, Index M>
template<typename T, Index N>
T
Function_Base<FunctionDerived, S, M>::
value(FunctionDerived & f, Vector<T, N> const & x)
{
  assert(x.get_dimension() <= DIMENSION);

  Vector<T, N> const
  r = residual(f, x);

  return 0.5 * dot(r, r);
}

//
//
//
template<typename FunctionDerived, typename S, Index M>
template<typename T, Index N>
Vector<T, N>
Function_Base<FunctionDerived, S, M>::
gradient(FunctionDerived & f, Vector<T, N> const & x)
{
  using AD = FAD<T, N>;

  Index const
  dimension = x.get_dimension();

  assert(dimension <= DIMENSION);

  Vector<AD, N>
  x_ad(dimension);

  for (Index i{0}; i < dimension; ++i) {
    x_ad(i) = AD(dimension, i, x(i));
  }

  AD const
  f_ad = f.value(x_ad);

  Vector<T, N>
  gradient(dimension);

  for (Index i{0}; i < dimension; ++i) {
    gradient(i) = f_ad.dx(i);
  }

  return gradient;
}

//
//
//
template<typename FunctionDerived, typename S, Index M>
template<typename T, Index N>
Vector<T, N>
Function_Base<FunctionDerived, S, M>::
residual(FunctionDerived & f, Vector<T, N> const & x)
{
  return f.gradient(x);
}

//
//
//
template<typename FunctionDerived, typename S, Index M>
template<typename T, Index N>
Tensor<T, N>
Function_Base<FunctionDerived, S, M>::
hessian(FunctionDerived & f, Vector<T, N> const & x)
{
  using AD = FAD<T, N>;

  Index const
  dimension = x.get_dimension();

  assert(dimension <= DIMENSION);

  Vector<AD, N>
  x_ad(dimension);

  for (Index i{0}; i < dimension; ++i) {
    x_ad(i) = AD(dimension, i, x(i));
  }

  Vector<AD, N> const
  r_ad = f.gradient(x_ad);

  Tensor<T, N>
  Hessian(dimension);

  for (Index i{0}; i < dimension; ++i) {
    for (Index j{0}; j < dimension; ++j) {
      Hessian(i, j) = r_ad(i).dx(j);
    }
  }

  return Hessian;
}

//
//
//
template<typename FunctionDerived, typename S, Index M>
void
Function_Base<FunctionDerived, S, M>::
set_failed(char const * const msg)
{
  failed = true;
  failure_message = msg;
  return;
}

//
//
//
template<typename FunctionDerived, typename S, Index M>
bool
Function_Base<FunctionDerived, S, M>::
get_failed()
{
  return failed;
}

//
//
//
template<typename FunctionDerived, typename S, Index M>
void
Function_Base<FunctionDerived, S, M>::
clear_failed()
{
  failed = false;
  return;
}

//
//
//
template<typename FunctionDerived, typename S, Index M>
void
Function_Base<FunctionDerived, S, M>::
set_failure_message(char const * const msg)
{
  failure_message = msg;
  return;
}

//
//
//
template<typename FunctionDerived, typename S, Index M>
char const *
Function_Base<FunctionDerived, S, M>::
get_failure_message()
{
  return failure_message;
}

//
//
//
template<typename ConstraintDerived, typename S, Index NC, Index NV>
template<typename T, Index N>
Vector<T, NC>
Equality_Constraint<ConstraintDerived, S, NC, NV>::
value(ConstraintDerived & c, Vector<T, N> const & x)
{
  assert(x.get_dimension() <= NUM_VAR);
  return c.value(x);
}

//
//
//
template<typename ConstraintDerived, typename S, Index NC, Index NV>
template<typename T, Index N>
Matrix<T, NC, NV>
Equality_Constraint<ConstraintDerived, S, NC, NV>::
gradient(ConstraintDerived & c, Vector<T, N> const & x)
{
  using AD = FAD<T, N>;

  Index const
  num_var = x.get_dimension();

  assert(num_var <= NUM_VAR);

  Vector<AD, N>
  x_ad(num_var);

  for (Index i{0}; i < num_var; ++i) {
    x_ad(i) = AD(num_var, i, x(i));
  }

  Vector<AD, NC> const
  r_ad = c.value(x_ad);

  Index const
  num_constr = r_ad.get_dimension();

  Matrix<T, NC, NV>
  Jacobian(num_constr, num_var);

  for (Index i{0}; i < num_constr; ++i) {
    for (Index j{0}; j < num_var; ++j) {
      Jacobian(i, j) = r_ad(i).dx(j);
    }
  }

  return Jacobian;
}

//
//
//
template<typename T, Index N>
Bounds<T, N>::
Bounds(Vector<T, N> const & l, Vector<T, N> const & u) : lower(l), upper(u)
{
  return;
}

//
//
//
template<typename T, Index N>
Minimizer<T, N>::
Minimizer()
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
template<typename STEP, typename FN>
void
Minimizer<T, N>::
solve(STEP & step_method, FN & fn, Vector<T, N> & soln)
{
  step_method_name = step_method.name();
  function_name = FN::NAME;
  initial_guess = soln;

  Vector<T, N>
  resi = fn.gradient(soln);

  initial_value = fn.value(soln);
  previous_value = initial_value;
  failed = failed || fn.get_failed();
  if (fn.get_failed() == true) failure_message = fn.get_failure_message();
  initial_norm = norm(resi);

  updateConvergenceCriterion(initial_norm);

  step_method.initialize(fn, soln, resi);

  while (continueSolve() == true) {

    Vector<T, N> const
    step = step_method.step(fn, soln, resi);

    soln += step;

    resi = fn.gradient(soln);

    failed = failed || fn.get_failed();
    if (fn.get_failed() == true) failure_message = fn.get_failure_message();

    T const
    norm_resi = norm(resi);

    updateConvergenceCriterion(norm_resi);

    T const
    value = fn.value(soln);

    failed = failed || fn.get_failed();
    if (fn.get_failed() == true) failure_message = fn.get_failure_message();

    updateDivergenceCriterion(value);

    ++num_iter;
  }

  recordFinals(fn, soln);
  return;
}

//
//
//
template<typename T, Index N>
void
Minimizer<T, N>::
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
  os << "Max Iters    : " << max_num_iter << '\n';
  os << "Iters Taken  : " << num_iter << '\n';

  os << std::scientific << std::setprecision(17);

  os << "Initial |R|  : " << std::setw(24) << initial_norm << '\n';
  os << "Abs Tol      : " << std::setw(24) << abs_tol << '\n';
  os << "Abs Error    : " << std::setw(24) << abs_error << '\n';
  os << "Rel Tol      : " << std::setw(24) << rel_tol << '\n';
  os << "Rel Error    : " << std::setw(24) << rel_error << '\n';
  os << "Initial X    : " << initial_guess << '\n';
  os << "Initial f(X) : " << std::setw(24) << initial_value << '\n';
  os << "Final X      : " << final_soln << '\n';
  os << "Final f(X)   : " << std::setw(24) << final_value << '\n';
  os << "Final Df(X)  : " << final_gradient << '\n';
  os << "Final DDf(X) : " << final_hessian << '\n';
  os << '\n';

  return;
}

//
//
//
template<typename T, Index N>
void
Minimizer<T, N>::
updateConvergenceCriterion(T const ae)
{
  abs_error = ae;
  rel_error = initial_norm > 0.0 ? abs_error / initial_norm : T(0.0);

  bool const
  converged_absolute = abs_error <= abs_tol;

  bool const
  converged_relative = rel_error <= rel_tol;

  converged = converged_absolute || converged_relative;

  bool const
  converged_acceptable = abs_error <= acc_tol && num_iter == max_num_iter - 1;  
  
  if (converged == false && converged_acceptable == true) {
    converged = true;
    warning = true;
    warning_message = "Reached acceptable tolerance";
  }  

  return;
}

//
//
//
template<typename T, Index N>
void
Minimizer<T, N>::
updateDivergenceCriterion(T const fn_value)
{
  monotonic = fn_value <= previous_value;

  if (enforce_monotonicity == true && monotonic == false) {
    failed = true;
    failure_message = "Non-monotonic";
  }

  T reduction_ratio = previous_value > 0.0 ? (fn_value / previous_value) : 0.0;

  if (reduction_ratio > stagnation_tol) {
    ++num_stagnation_iter;
  }
  else {
    num_stagnation_iter = 0;
  }

  non_stagnant = num_stagnation_iter < max_stagnation_iter;

  // only set warning for stagnant residual
  if (enforce_non_stagnation == true && non_stagnant == false) {
    warning = true;
    warning_message = "Stagnant residual";
  }

  previous_value = fn_value;

  bounded = fn_value <= growth_limit * initial_value;

  if (enforce_boundedness == true && bounded == false) {
    failed = true;
    failure_message = "Growing unbounded";
  }

  return;
}

//
//
//
template<typename T, Index N>
bool
Minimizer<T, N>::
continueSolve() const
{
  // If failure has occurred, stop immediately.
  if (failed == true) return false;

  // Regardless of other criteria, if the residual is zero stop solving.
  bool const
  zero_resi = ((abs_error > 0.0) == false);

  if (zero_resi == true) return false;

  // Minimum iterations takes precedence over maximum iterations and
  // convergence. Continue solving if not exceeded.
  bool const
  exceeds_min_iter = num_iter >= min_num_iter;

  if (exceeds_min_iter == false) return true;

  // Maximum iterations takes precedence over convergence.
  // Stop solving if exceeded.
  bool const
  exceeds_max_iter = num_iter >= max_num_iter;

  if (exceeds_max_iter == true) return false;

  // Lastly check for convergence.
  bool const
  continue_solve = (converged == false);

  return continue_solve;
}

//
//
//
template<typename T, Index N>
template<typename FN>
void
Minimizer<T, N>::
recordFinals(FN & fn, Vector<T, N> const & x)
{
  final_soln = x;
  final_value = fn.value(x);
  final_gradient = fn.gradient(x);
  final_hessian = fn.hessian(x);
}

//
// Linear solver for step methods
//
template<typename FN, typename T, Index N>
Vector<T, N>
StepBase<FN, T, N>::
lin_solve(Tensor<T, N> const & A, Vector<T, N> const & b)
{
  return solve(A, b, preconditioner_type);
}

//
// Linear solve for trust region subproblems
//
template<typename T, Index N>
Vector<T, N>
TrustRegionSubproblemBase<T, N>::
lin_solve(Tensor<T, N> const & A, Vector<T, N> const & b)
{
  return solve(A, b, preconditioner_type);
}

//
// Trust region subproblem with given objective function.
// Exact algorithm, Nocedal 2nd Ed 4.3
//
template<typename T, Index N>
Vector<T, N>
TrustRegionExactValue<T, N>::
step(Tensor<T, N> const & Hessian, Vector<T, N> const & gradient)
{
  Index const
  dimension = gradient.get_dimension();

  Tensor<T, N> const
  I = identity<T, N>(dimension);

  Vector<T, N>
  step(dimension);

  // set to Hessian norm to ensure that K is positive definite
  T
  lambda = norm(Hessian);

  for (Index i{0}; i < max_num_iter; ++i) {

    Tensor<T, N> const
    K = Hessian + lambda * I;

    Tensor<T, N>
    L(Filler::ZEROS);

    bool
    is_posdef{false};

    std::tie(L, is_posdef) = cholesky(K);

    if (is_posdef == false) {
      MT_ERROR_EXIT("Trust region subproblem encountered singular Hessian.");
    }

    step = - this->lin_solve(K, gradient);

    Vector<T, N> const
    q = this->lin_solve(L, step);

    T const
    np = norm(step);

    T const
    nps = np * np;

    T const
    nqs = norm_square(q);

    T const
    lambda_incr = nps * (np - region_size) / nqs / region_size;

    lambda += std::max(lambda_incr, 0.0);
      
  }

  return step;
}
 

//
// Trust region subproblem with given gradient/residual.
// Exact algorithm, Nocedal 2nd Ed 4.3
//
template<typename T, Index N>
Vector<T, N>
TrustRegionExactGradient<T, N>::
step(Tensor<T, N> const & Hessian, Vector<T, N> const & gradient)
{
  Index const
  dimension = gradient.get_dimension();

  Tensor<T, N> const
  I = identity<T, N>(dimension);

  Vector<T, N>
  step(dimension);

  // set to Hessian norm to ensure that K is positive definite
  T
  lambda = norm(Hessian);

  for (Index i{0}; i < max_num_iter; ++i) {

    
    Tensor<T, N> const
    HTH = dot(transpose(Hessian), Hessian);

    Vector<T, N> const
    HTr = dot(transpose(Hessian), gradient);
      
    Tensor<T, N> const
    K = HTH + lambda * I;

    Tensor<T, N>
    L(Filler::ZEROS);

    bool
    is_posdef{false};

    std::tie(L, is_posdef) = cholesky(K);

    if (is_posdef == false) {
      MT_ERROR_EXIT("Trust region subproblem encountered singular Hessian.");
    }

    step = - this->lin_solve(K, HTr);

    Vector<T, N> const
    q = this->lin_solve(L, step);

    T const
    np = norm(step);

    T const
    nps = np * np;

    T const
    nqs = norm_square(q);

    T const
    lambda_incr = nps * (np - region_size) / nqs / region_size;

    lambda += std::max(lambda_incr, 0.0);
      
  }

  return step;
}


//
// Trust region subproblem. Dog-leg algorithm with given gradient/residual.
// See pp. 73 - 74, Nocedal 2nd Ed 4.3
//
template<typename T, Index N>
Vector<T, N>
TrustRegionDogLegValue<T, N>::
step(Tensor<T, N> const & Hessian, Vector<T, N> const & gradient)
{    
  T const
  normg_H = dot(gradient, dot(Hessian, gradient));

  if (normg_H < 0.0) {
    MT_ERROR_EXIT("Trust region subproblem encountered singular Hessian.");
  }

  if (normg_H == 0.0) return Vector<T, N>(gradient.get_dimension(), Filler::ZEROS);

  T const
  normg_squared = dot(gradient, gradient);

  Vector<T, N> const
  step_minimizer = - normg_squared / normg_H * gradient;

  Vector<T, N> const
  step_unconstrained = - this->lin_solve(Hessian, gradient);

  T const
  normg_cubed = norm(gradient) * normg_squared;

  T const
  tau = std::min(1.0, normg_cubed / (region_size * normg_H));
    
  Vector<T, N>
  step{tau * step_minimizer};

  if (tau > 1.0) {
    step = step_minimizer + (tau - 1.0) * (step_unconstrained - step_minimizer);
  }

  return step;
}
 

//
// Trust region subproblem. Dog-leg algorithm with given gradient/residual.
// See pp. 73 - 74, Nocedal 2nd Ed 4.3
//
template<typename T, Index N>
Vector<T, N>
TrustRegionDogLegGradient<T, N>::
step(Tensor<T, N> const & Hessian, Vector<T, N> const & gradient)
{    
  Vector<T, N> const
  HTr = dot(transpose(Hessian), gradient);
  
  T const
  normHTr_squared = dot(HTr, HTr);

  Tensor<T, N> const
  HTH = dot(transpose(Hessian), Hessian);

  T const
  normHTr_HTH = dot(HTr, dot(HTH, HTr));

  if (normHTr_HTH <  0.0) {
    MT_ERROR_EXIT("Trust region subproblem encountered singular Hessian.");
  }

  if (normHTr_HTH == 0.0) return Vector<T, N>(gradient.get_dimension(), Filler::ZEROS);

  T const
  normHTr_cubed = norm(HTr) * normHTr_squared;

  T const
  tau = std::min(1.0, normHTr_cubed / (region_size * normHTr_HTH));
    
  Vector<T, N> const
  step_minimizer = - normHTr_squared / normHTr_HTH * HTr;

  Vector<T, N> const
  step_unconstrained = - this->lin_solve(HTH, HTr);

  Vector<T, N>
  step{tau * step_minimizer};

  if (tau > 1.0) {
    step = step_minimizer + (tau - 1.0) * (step_unconstrained - step_minimizer);
  }

  return step;
}
 
//
// Newton line search
//
template<typename T, Index N>
template<typename FN>
Vector<T, N>
NewtonLineSearch<T, N>::
step(FN & fn, Vector<T, N> const & direction, Vector<T, N> const & soln)
{
  Index const
  dimension = soln.get_dimension();

  Vector<T, N>
  step(dimension, Filler::ZEROS);

  T const
  projection_direction = dot(direction, direction);

  for (Index i{0}; i < max_num_iter; ++i) {

    Vector<T, N> const
    soln_next = soln + step;

    Vector<T, N> const
    gradient_next = fn.gradient(soln_next);

    Tensor<T, N> const
    Hessian_next = fn.hessian(soln_next);

    T const
    projection = dot(gradient_next, direction);

    T const
    contraction = dot(direction, dot(Hessian_next, direction));

    T const
    step_length = - projection / contraction;

    step += step_length * direction;

    T const
    ls_length2 = step_length * step_length * projection_direction;

    bool const
    line_search_converged = ls_length2 <= tolerance * tolerance;

    if (line_search_converged == true) break;

  }

  return step;
}

//
// Back-tracking line search
// Taken from box I.1 in the appendix of Armero and Perez-Foguet (Jan. 2000)
//
template<typename T, Index N>
template<typename FN>
Vector<T, N>
BacktrackingLineSearch<T, N>::
step(FN & fn, Vector<T, N> const & direction, Vector<T, N> const & soln)
{
  Index const
  dimension = soln.get_dimension();

  Vector <T, N>
  step = direction;

  Vector<T, N>
  step_line_search(dimension, Filler::ZEROS);

  Vector<T, N> const
  resid = fn.gradient(soln);

  Vector<T, N> const
  resid_newton = fn.gradient(soln + step);

  Index
  line_iter{0};

  T const
  resid_norm = dot(resid, resid);

  T const
  resid_newton_norm = dot(resid_newton, resid_newton);

  T
  resid_line_norm = resid_newton_norm;

  for (Index i{0}; i < max_num_iter; i++) {

    if (line_iter == max_line_iter) {

      alpha = 1.0;
      line_iter = 0;
      search_parameter += search_increment;
    }

    if (search_parameter >= 1.0) break;

    step_line_search = alpha * step;

    Vector<T, N> const
    soln_line_search = soln + step_line_search;

    Vector<T, N> const
    gradient_line_search = fn.gradient(soln_line_search);

    T const
    resid_line_norm_old = resid_line_norm;

    resid_line_norm = dot(gradient_line_search, gradient_line_search);

    if (resid_line_norm <= (search_parameter * resid_norm)) {
      step = step_line_search;
      break;
    }

    T const
    num = 0.25 * alpha * alpha * resid_line_norm_old;

    T const
    den = resid_line_norm + resid_line_norm_old * (0.5 * alpha - 1.0);

    T const
    den_tol = 1.0e-10;

    bool const
    is_zero_den = (std::fabs(den) < den_tol);

    alpha = is_zero_den == true ? 0.5 * alpha : max(0.5 * alpha, num / den);

    line_iter++;

  } //Index i

  return step;
}

//
//
//
template<typename FN, typename T, Index N>
void
NewtonStep<FN, T, N>::
initialize(FN &, Vector<T, N> const &, Vector<T, N> const &)
{
  return;
}

//
// Plain Newton step.
//
template<typename FN, typename T, Index N>
Vector<T, N>
NewtonStep<FN, T, N>::
step(FN & fn, Vector<T, N> const & soln, Vector<T, N> const & resi)
{
  Tensor<T, N> const
  Hessian = fn.hessian(soln);

  Vector<T, N> const
  step = - this->lin_solve(Hessian, resi);

  return step;
}


//
//
//
template<typename FN, typename T, Index N>
void
NewtonWithLineSearchStep<FN, T, N>::
initialize(FN &, Vector<T, N> const &, Vector<T, N> const &)
{
  return;
}

//
// Plain Newton step with line search
//
template<typename FN, typename T, Index N>
Vector<T, N>
NewtonWithLineSearchStep<FN, T, N>::
step(FN & fn, Vector<T, N> const & soln, Vector<T, N> const & resi)
{
  Tensor<T, N> const
  Hessian = fn.hessian(soln);

  Vector<T, N> const
  step = - this->lin_solve(Hessian, resi);

  // Newton back-tracking line search.
  BacktrackingLineSearch<T, N>
  newton_backtrack_ls;

  Vector<T, N> const
  ls_step = newton_backtrack_ls.step(fn, step, soln);

  return ls_step;
}

//
//
//
template<typename FN, typename T, Index N>
void
TrustRegionStep<FN, T, N>::
initialize(FN &, Vector<T, N> const &, Vector<T, N> const &)
{
  region_size = initial_region_size;

  return;
}

//
// Trust Region method.  See Nocedal's algorithm 11.5.
//
template<typename FN, typename T, Index N>
Vector<T, N>
TrustRegionStep<FN, T, N>::
step(FN & fn, Vector<T, N> const & soln, Vector<T, N> const & resi)
{
  Tensor<T, N> const
  Hessian = fn.hessian(soln);

  // Compute full Newton step first.
  Vector<T, N>
  step = - this->lin_solve(Hessian, resi);

  T const
  norm_step = minitensor::norm(step);

  // Take full Newton step if inside trust region.
  if (norm_step < region_size) return step;

  // Trust region subproblem. Exact algorithm, Nocedal 2nd Ed 4.3
  TrustRegionExactValue<T, N>
  tr_exact;

  tr_exact.region_size = region_size;

  step = tr_exact.step(Hessian, resi);

  Vector<T, N> const
  soln_next = soln + step;

  Vector<T, N> const
  resi_next = fn.gradient(soln_next);

  // Compute reduction factor \rho_k in Nocedal's algorithm 11.5
  T const
  nr = norm_square(resi);

  T const
  nrp = norm_square(resi_next);

  T const
  nrKp = norm_square(resi + dot(Hessian, step));

  T const
  reduction = (nr - nrp) / (nr - nrKp);

  // Determine whether the trust region should be increased, decreased
  // or left the same.
  T const
  computed_size = norm(step);

  if (reduction < 0.25) {

    region_size = 0.25 * computed_size;

  } else {

    bool const
    increase_region_size = reduction > 0.75;

    if (increase_region_size == true) {
      region_size = std::min(2.0 * region_size, max_region_size);
    }

  }

  if (reduction <= min_reduction) {
    step.fill(Filler::ZEROS);
  }

  return step;
}

//
//
//
template<typename FN, typename T, Index N>
void
ConjugateGradientStep<FN, T, N>::
initialize(FN & fn, Vector<T, N> const & soln, Vector<T, N> const & gradient)
{
  Tensor<T, N> const
  Hessian = fn.hessian(soln);

  precon_resi = - this->lin_solve(Hessian, gradient);

  search_direction = precon_resi;

  projection_new = - dot(gradient, search_direction);

  restart_directions_counter = 0;

  return;
}

//
// Conjugate Gradient Method step.
// For now the Gram-Schmidt method is fixed to Polak-Ribiere
// and preconditioning with the Hessian.
// This is taken from J.R. Shewchuck "painless" conjugate gradient
// manuscript that is all over the place on the net.
//
template<typename FN, typename T, Index N>
Vector<T, N>
ConjugateGradientStep<FN, T, N>::
step(FN & fn, Vector<T, N> const & soln, Vector<T, N> const &)
{
  // Newton line search.
  NewtonLineSearch<T, N>
  newton_ls;

  Vector<T, N> const
  step = newton_ls.step(fn, search_direction, soln);

  Vector<T, N> const
  soln_next = soln + step;

  Vector<T, N> const
  gradient_next = fn.gradient(soln_next);

  T const
  projection_old = projection_new;

  T const
  projection_mid = - dot(gradient_next, precon_resi);

  Tensor<T, N> const
  Hessian = fn.hessian(soln_next);

  precon_resi = - this->lin_solve(Hessian, gradient_next);

  projection_new = - dot(gradient_next, precon_resi);

  T const
  gram_schmidt_factor = (projection_new - projection_mid) / projection_old;

  ++restart_directions_counter;

  bool const
  rewind = restart_directions_counter == restart_directions_interval;

  bool const
  bad_directions = gram_schmidt_factor <= 0.0;

  bool const
  restart_directions = rewind || bad_directions;

  if (restart_directions == true) {

    search_direction = precon_resi;
    restart_directions_counter = 0;

  } else {

    search_direction = precon_resi + gram_schmidt_factor * search_direction;

  }

  return step;
}

//
//
//
template<typename FN, typename T, Index N>
void
LineSearchRegularizedStep<FN, T, N>::
initialize(FN &, Vector<T, N> const &, Vector<T, N> const &)
{
  return;
}

//
// Line Search Newton-like method.  See Nocedal's algorithm 11.4.
//
template<typename FN, typename T, Index N>
Vector<T, N>
LineSearchRegularizedStep<FN, T, N>::
step(FN & fn, Vector<T, N> const & soln, Vector<T, N> const & gradient)
{
  Index const
  dimension = soln.get_dimension();

  Tensor<T, N> const
  Hessian = fn.hessian(soln);

  Vector<T, N>
  step(dimension);

  bool const
  bad_hessian = inv_cond(Hessian) * hessian_cond_tol < 1.0;

  // Regularize Hessian if it is bad.
  if (bad_hessian == true) {

    // Trust region subproblem. Exact algorithm, Nocedal 2nd Ed 4.3
    TrustRegionExactValue<T, N>
    tr_exact;

    tr_exact.region_size = step_length;

    step = tr_exact.step(Hessian, gradient);
  } else {

    // Standard Newton step
    step = - this->lin_solve(Hessian, gradient);

  }

  NewtonLineSearch<T, N>
  newton_ls;

  Vector<T, N> const
  ls_step = newton_ls.step(fn, step, soln);

  return ls_step;
}

//
//
//
template<typename FN, typename T, Index N>
std::unique_ptr<StepBase<FN, T, N>>
stepFactory(StepType step_type)
{
  using STUP = std::unique_ptr<StepBase<FN, T, N>>;

  switch (step_type) {

  default:
    MT_ERROR_EXIT("ERROR: Unknown step method type.");
    break;

  case StepType::NEWTON:
    return STUP(new NewtonStep<FN, T, N>());
    break;

  case StepType::TRUST_REGION:
    return STUP(new TrustRegionStep<FN, T, N>());
    break;

  case StepType::CG:
    return STUP(new ConjugateGradientStep<FN, T, N>());
    break;

  case StepType::LINE_SEARCH_REG:
    return STUP(new LineSearchRegularizedStep<FN, T, N>());
    break;

  case StepType::NEWTON_LS:
    return STUP(new NewtonWithLineSearchStep<FN, T, N>());
    break;    
  }

  return STUP(nullptr);
}

} // namespace minitensor
