// @HEADER
// ************************************************************************
//
//                           Intrepid Package
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

namespace Intrepid
{

//
//
//
template<typename FunctionDerived, typename S>
template<typename T, Index N>
T
Function_Base<FunctionDerived, S>::
value(FunctionDerived & f, Vector<T, N> const & x)
{
  Vector<T, N> const
  r = f.gradient(x);

  return 0.5 * dot(r, r);
}

//
//
//
template<typename FunctionDerived, typename S>
template<typename T, Index N>
Vector<T, N>
Function_Base<FunctionDerived, S>::
gradient(FunctionDerived & f, Vector<T, N> const & x)
{
  using AD = FAD<T, N>;

  Index const
  dimension = x.get_dimension();

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
template<typename FunctionDerived, typename S>
template<typename T, Index N>
Tensor<T, N>
Function_Base<FunctionDerived, S>::
hessian(FunctionDerived & f, Vector<T, N> const & x)
{
  using AD = FAD<T, N>;

  Index const
  dimension = x.get_dimension();

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
  failed = failed || fn.failed;
  initial_norm = norm(resi);

  updateConvergenceCriterion(initial_norm);

  step_method.initialize(fn, soln, resi);

  while (continueSolve() == true) {

    Vector<T, N> const
    step = step_method.step(fn, soln, resi);

    soln += step;

    resi = fn.gradient(soln);

    failed = failed || fn.failed;

    T const
    norm_resi = norm(resi);

    updateConvergenceCriterion(norm_resi);

    T const
    value = fn.value(soln);

    failed = failed || fn.failed;

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
Minimizer<T, N>::
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

  // Last check for convergence.
  bool const
  continue_solve = (converged == false);

  return continue_solve;
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
// Trust region subproblem. Exact algorithm, Nocedal 2nd Ed 4.3
//
template<typename T, Index N>
Vector<T, N>
TrustRegionExact<T, N>::
step(Tensor<T, N> const & Hessian, Vector<T, N> const & gradient)
{
  Index const
  dimension = gradient.get_dimension();

  Tensor<T, N> const
  I = identity<T, N>(dimension);

  Vector<T, N>
  step(dimension);

  T
  lambda = initial_lambda;

  for (Index i{0}; i < max_num_iter; ++i) {

    Tensor<T, N> const
    K = Hessian + lambda * I;

    Tensor<T, N> const
    L = cholesky(K).first;

    step = - Intrepid::solve(K, gradient);

    Vector<T, N> const
    q = Intrepid::solve(L, step);

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
  step(dimension, ZEROS);

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
  step_line_search(dimension, ZEROS);

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
    quad_approx = num / den;

    alpha = std::max(0.5 * alpha, quad_approx);

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
  step = - Intrepid::solve(Hessian, resi);

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
  step = - Intrepid::solve(Hessian, resi);

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
  step = - Intrepid::solve(Hessian, resi);

  T const
  norm_step = Intrepid::norm(step);

  // Take full Newton step if inside trust region.
  if (norm_step < region_size) return step;

  // Trust region subproblem. Exact algorithm, Nocedal 2nd Ed 4.3
  TrustRegionExact<T, N>
  tr_exact;

  tr_exact.initial_lambda = 0.0;
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
    step.fill(ZEROS);
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

  precon_resi = - Intrepid::solve(Hessian, gradient);

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

  precon_resi = - Intrepid::solve(Hessian, gradient_next);

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
    TrustRegionExact<T, N>
    tr_exact;

    tr_exact.initial_lambda = 1.0;
    tr_exact.region_size = step_length;

    step = tr_exact.step(Hessian, gradient);

  } else {

    // Standard Newton step
    step = - Intrepid::solve(Hessian, gradient);

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
    std::cerr << __PRETTY_FUNCTION__ << '\n';
    std::cerr << "ERROR: Unknown step method type\n";
    exit(1);
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

} // namespace Intrepid
