#include <Akri_DiagWriter.hpp>
#include <Akri_DistributedVector.hpp>
#include <Akri_Optimize.hpp>
#include <Akri_OptimizeTypes.hpp>
#include <krino/mesh_utils/Akri_AllReduce.hpp>
#include <Akri_ObjectiveInterface.hpp>
#include <stk_math/StkVector.hpp>
#include <stk_util/util/ReportHandler.hpp>
#include <deque>
#include <functional>
#include <iostream>
#include <tuple>
#include <vector>

namespace krino {


stk::math::Vector3d xpby(const stk::math::Vector3d & x, const double b, const stk::math::Vector3d & y)
{
  return x + b*y;
}

stk::math::Vector3d scalar_times_vector(const double a, const stk::math::Vector3d & x)
{
  return a*x;
}

stk::math::Vector3d vectorSubtract(const stk::math::Vector3d& a, const stk::math::Vector3d& b)
{
  return a-b;
}

template<typename Objective, typename VEC>
double compute_scaled_value(const Objective & objFn, const double scale, const VEC &x)
{
  return objFn.compute_value(x)/scale;
}

template<typename Objective, typename VEC>
void fill_scaled_gradient(const Objective & objFn, const double scale, const VEC &x,  VEC & grad)
{
  const double invScale = 1./scale;
  objFn.fill_gradient(x,grad);
  for (double & g : grad)
    g *= invScale;
}

template<typename VEC>
void rescale(double & scale, double & f, VEC & grad, double & gradMag)
{
  const double invScale = 1./f;

  scale *= f;
  f = 1.;

  for (double & g : grad)
    g *= invScale;

  gradMag *= invScale;
}

template<typename Objective, typename VEC>
void project_solution(const Objective & objFn, VEC &x)
{
  objFn.project_solution(x);
}

//----------------------------------------
// Cubic minimizer in [xLo, xHi]
//   matches phi & dphi at both ends.
// If the root is outside [lo+eps, hi-eps], falls back to midpoint.
//----------------------------------------
double cubic_minimizer(
    const double xLo,  const double xHi,
    const double fLo,  const double fHi,
    const double dfLo, const double dfHi,
    const bool doAvoidSmallSteps)
{
    // compute coefficients of the interpolating cubic’s derivative
    const double d1 = dfLo + dfHi - 3.0*(fLo - fHi)/(xLo - xHi);
    const double d2sq = d1*d1 - dfLo*dfHi;
    const double d2   = (d2sq > 0.0 ? std::sqrt(d2sq) : 0.0);

    const double numerator   = (xHi - xLo)*(dfHi + d2 - d1);
    const double denominator = (dfHi - dfLo + 2.0*d2);
    double xCubic = xHi - (denominator != 0.0 ? numerator/denominator : 0.5*(xLo + xHi));
    if (!doAvoidSmallSteps)
      return std::max(xLo,std::min(xHi,xCubic));

    // clamp to [lo+tol, hi-tol] to avoid too-small steps
    const double tol = 0.1;
    const double loBound = std::min(xLo, xHi) + tol*std::abs(xHi - xLo);
    const double hiBound = std::max(xLo, xHi) - tol*std::abs(xHi - xLo);
    if (xCubic < loBound || xCubic > hiBound)
      xCubic = 0.5*(xLo + xHi);
    return xCubic;
}

//----------------------------------------
// Zoom phase: given bracket [alphaLo,alphaHi],
// shrink it via cubic interpolation until
// Wolfe condition is satisfied (or max_iters).
//----------------------------------------
template<typename Objective, typename VEC>
void zoom_cubic(
    const Objective & objFn,
    const double scale,
    const VEC &x0, const VEC &p,
    double alphaLo, double alphaHi,
    double fLo,   double fHi,
    double dfLo,  double dfHi,
    VEC & x,
    double & f,
    VEC & grad,
    const double minAlpha,
    const double c1, const double c2)
{
    const double f0 = fLo;
    const double df0 = dfLo;
    bool doAvoidSmallSteps = false;

    unsigned iter = 0;
    while(std::abs(alphaHi - alphaLo) > minAlpha)
    {
        if (iter++ > 2)
          doAvoidSmallSteps = true;

        // trial point by cubic interpolation
        double alpha_j = cubic_minimizer( alphaLo, alphaHi, fLo,   fHi, dfLo,  dfHi, doAvoidSmallSteps);
        x = xpby(x0, alpha_j, p);
        f = compute_scaled_value(objFn, scale, x);
        fill_scaled_gradient(objFn, scale, x, grad);
        const double df = Dot(grad, p);

        // check Armijo:
        if (f > f0 + c1*alpha_j*df0 || f >= fLo)
        {
            // insufficient decrease → shrink right end
            alphaHi = alpha_j;
            fHi   = f;
            dfHi  = df;
        }
        else
        {
            // Armijo ok, check curvature:
            if (df >= c2*df0)
            {
                return;  // Wolfe satisfied
            }

            // curvature failed
            if (df*(alphaHi - alphaLo) >= 0.0)
            {
                // force bracket to have opposite slope
                alphaHi = alphaLo;
                fHi   = fLo;
                dfHi  = dfLo;
            }
            // shrink left end
            alphaLo = alpha_j;
            fLo   = f;
            dfLo  = df;
        }
    }

    if (f >= f0)
    {
      // fallback: go back to x0 solution
      x = x0;
      f = f0;
      fill_scaled_gradient(objFn, scale, x, grad);
    }
}

template<typename Objective, typename VEC>
void line_search_wolfe_cubic(const Objective & objFn,
    const double scale,
    const VEC &x0,
    const double f0,
    const VEC &grad0,
    const VEC &p,       // descent dir: grad(x)^T p < 0
    VEC & x,
    double & f,
    VEC & grad,
    const double minAlpha,
    const double c1 = 1e-4,  // Armijo param
    const double c2 = 0.9)   // curvature param (Wolfe)
{
    // evaluate at alpha=0
    const double df0 = Dot(grad0, p);
    if (df0 >= 0)
    {
      x = x0;
      f = f0;
      grad = grad0;
      return;
    }

    // evaluate at alpha
    double alpha = 1.0;
    x = xpby(x0, alpha, p);
    f = compute_scaled_value(objFn, scale, x);

    while (!std::isfinite(f) || (f > f0 + c1*alpha*df0)) // Armijo
    {
      alpha *= 0.5;

      if (alpha < minAlpha)
      {
        x = x0;
        f = f0;
        grad = grad0;
        return;
      }

      x = xpby(x0, alpha, p);
      f = compute_scaled_value(objFn, scale, x);
    }

    fill_scaled_gradient(objFn, scale, x, grad);
    const double df = Dot(grad, p);

    // unless alpha already satisfies Armijo and strong Wolfe, and is still decreasing, use bracketed cubic minimizer
    if (df > 0 || std::abs(df) > -c2*df0)
    {
      zoom_cubic(objFn, scale,
          x0, p,
          /*alphaLo=*/0.0,/*alphaHi=*/alpha,
          /*fLo=*/f0, /*fHi=*/f,
          /*dfLo=*/df0, /*dfHi=*/df,
          x, f, grad,
          minAlpha,
          c1, c2);
    }
}

template<typename Objective, typename VEC>
void line_search_armijo(const Objective & objFn,
    const double scale,
    const VEC & x0,
    const double f0,
    const VEC & grad0,
    const VEC & p,
    VEC & x,
    double & f,
    const double c1 = 1.e-4,
    const unsigned maxIters = 45)
{
  const double gradDotP = Dot(grad0, p); // Directional derivative
  if (gradDotP >= 0)
  {
    x = x0;
    f = f0;
    return;
  }

  double alpha = 1.0; // Start with maximum step size

  for(unsigned iter=0; iter<maxIters; ++iter)
  {
    x = xpby(x0, alpha, p);
    f = compute_scaled_value(objFn, scale, x);

    // Check Armijo condition (sufficient decrease)
    if (f <= f0 + c1 * alpha * gradDotP) // won't eval to true if f is NaN
      return;

    alpha *= 0.5; // Reduce step size
  }

  if (f >= f0)
  {
    x = x0; // Fallback to 0 step
    f = f0;
  }
}

// -----------------------------------------------------------------------------
// Memory-update helper for L-BFGS
// -----------------------------------------------------------------------------
template<typename VEC>
void update_LBFGS_Memory(
    std::deque<VEC>& sList,
    std::deque<VEC>& yList,
    std::deque<double>& rhoList,
    const VEC& s,
    const VEC& y,
    const unsigned maxLevels)
{
    const double ys = Dot(y, s);
    if (ys <= 0.0)
    {
        // curvature condition failed; skip update
        return;
    }

    const double rho = 1.0 / ys;

    // if we're already storing m pairs, drop the oldest
    if (sList.size() == maxLevels)
    {
        sList.pop_front();
        yList.pop_front();
        rhoList.pop_front();
    }

    sList.push_back(s);
    yList.push_back(y);
    rhoList.push_back(rho);
}

template<typename VEC>
void reset_LBFGS_Memory(std::deque<VEC>& sList, std::deque<VEC>& yList, std::deque<double>& rhoList)
{
    sList.clear();
    yList.clear();
    rhoList.clear();
}

// -----------------------------------------------------------------------------
// Two-loop recursion to compute H_k * grad (search direction)
// -----------------------------------------------------------------------------
template<typename VEC>
VEC compute_LBFGS_Hgrad(
    const VEC& grad,
    const std::deque<VEC>& sList,
    const std::deque<VEC>& yList,
    const std::deque<double>& rhoList)
{
    const unsigned numLevels = sList.size();
    VEC q = grad;
    std::vector<double> alpha(numLevels);
    std::vector<double> beta(numLevels);

    // backward loop
    for (int i = numLevels - 1; i >= 0; --i)
    {
        alpha[i] = rhoList[i] * Dot(sList[i], q);
        q = xpby(q, -alpha[i], yList[i]);
    }

    // initial Hessian-scaling: gamma * I
    double gamma = 1.0;
    if (numLevels > 0)
    {
        const auto& sLast = sList.back();
        const auto& yLast = yList.back();
        const double sy = Dot(sLast, yLast);
        const double yy = Dot(yLast, yLast);
        if (yy > 1e-20)
            gamma = sy / yy;
    }
    VEC r = scalar_times_vector(gamma, q);

    // forward loop
    for (unsigned i = 0; i < numLevels; ++i)
    {
        beta[i] = rhoList[i] * Dot(yList[i], r);
        r = xpby(r, alpha[i] - beta[i], sList[i]);
    }

    return r;  // this is H_k * grad
}

// -----------------------------------------------------------------------------
// L-BFGS optimizer
// -----------------------------------------------------------------------------
template<typename Objective, typename VEC>
void lbfgs(const Objective & objFn,
    VEC& x,
    const double xTol,
    const double gradTol,
    const double scale0,
    const unsigned maxIter,
    const unsigned maxLevels,
    const bool doOutput)
{
    double scale = scale0;
    std::deque<VEC> sList, yList;
    std::deque<double> rhoList;

    double f = compute_scaled_value(objFn, scale, x);
    VEC grad;
    fill_scaled_gradient(objFn, scale, x, grad);

    double gradMag = std::sqrt(Dot(grad, grad));
    double dxMag = 0.;

    for (unsigned iter = 0; iter < maxIter; ++iter)
    {
        if (gradMag < gradTol)
        {
          if (doOutput)
            krinolog << "Gradient converged at iteration " << iter << stk::diag::dendl;
          return;
        }

        // compute search direction p = - H_k * grad
        const VEC Hgrad = compute_LBFGS_Hgrad(grad, sList, yList, rhoList);
        const VEC p = scalar_times_vector(-1., Hgrad);

        const double minAlpha = std::max(1.0e-14, std::min(1., xTol / std::sqrt(Dot(p, p))));

        const auto xOld = x;
        const auto gradOld = grad;
        const double gradMagOld = gradMag;
        const double fOld = f;
        line_search_wolfe_cubic(objFn, scale, xOld, fOld, gradOld, p, x, f, grad, minAlpha);
        project_solution(objFn, x);

        const VEC s = vectorSubtract(x, xOld);
        dxMag = std::sqrt(Dot(s,s));
        if (doOutput)
          krinolog << "L-BFGS iteration " << iter << ", |f|= " << std::abs(f*scale/scale0) << ", |grad|= " << gradMag*scale/scale0 << ", |dx|= " << dxMag << stk::diag::dendl;

        if (dxMag < xTol)
        {
          if (doOutput)
            krinolog << "Solution converged at iteration " << stk::diag::dendl;
          return;
        }

        gradMag = std::sqrt(Dot(grad, grad));

        if (gradMag < xTol*gradMagOld)
        {
          if (doOutput)
            krinolog << "Restarting/rescaling due to large drop in gradient magnitude " << gradMag/gradMagOld << "\n";
          rescale(scale, f, grad, gradMag);
          reset_LBFGS_Memory(sList, yList, rhoList);
        }
        else
        {
          const auto y = vectorSubtract(grad, gradOld);
          update_LBFGS_Memory(sList, yList, rhoList, s, y, maxLevels);
        }
    }

    krinolog << "Reached max iterations " << maxIter<< ", |f|= " << std::abs(f) << ", |grad|= " << gradMag << ", |dx|= " << dxMag << stk::diag::dendl;
    return;
}

template<typename Objective, typename VEC>
void steepest_descent(const Objective & objFn,
    VEC& x,
    const double xTol,
    const double gradTol,
    const double scale,
    const unsigned maxIter,
    const bool doOutput)
{
  VEC grad;
  fill_scaled_gradient(objFn, scale, x, grad);
  double f = compute_scaled_value(objFn, scale, x);

  double gradMag = std::sqrt(Dot(grad, grad));
  if (gradMag == 0.)
    return;
  double scaling = 1.0/gradMag;

  double dxSqr = 0.;

  for (unsigned iter = 0; iter < maxIter; ++iter)
  {
    gradMag = std::sqrt(Dot(grad, grad));
    if (gradMag < gradTol)
    {
      if (doOutput)
        krinolog << "Gradient converged at iteration " << iter << stk::diag::dendl;
      return;
    }

    // Compute the search direction
    const VEC p = scalar_times_vector(-scaling, grad);

    const auto xOld = x;
    const auto gradOld = grad;
    const double fOld = f;
    line_search_armijo(objFn, scale, xOld, fOld, gradOld, p, x, f);
    project_solution(objFn, x);

    const VEC s = vectorSubtract(x, xOld);
    dxSqr = Dot(s,s);
    if (doOutput)
      krinolog << "Steepest descent iteration " << iter << ", |f|= " << std::abs(f) << ", |grad|= " << gradMag << ", |dx|= " << std::sqrt(dxSqr) << stk::diag::dendl;
    if (std::sqrt(dxSqr) < xTol)
    {
      krinolog << "Solution converged at iteration " << iter << stk::diag::dendl;
      return;
    }

    fill_scaled_gradient(objFn, scale, x, grad);
    const VEC y = vectorSubtract(grad, gradOld);
    const double ys = Dot(y, s);
    if (ys > 0.)
      scaling = std::min(1.e6,std::max(1.e-6,dxSqr/ys)); //Barzilai–Borwein
  }

  krinolog << "Reached max iterations " << maxIter << ", |grad|= " << std::sqrt(Dot(grad,grad)) << ", |dx|= " << std::sqrt(dxSqr) << stk::diag::dendl;
}

template
void steepest_descent(const Objective3DInterface & objFn,
    stk::math::Vector3d & x,
    const double xTol,
    const double gradTol,
    const double scale,
    const unsigned maxIter,
    const bool doOutput);

template
void steepest_descent(const ObjectiveInterface & objFn,
    DistributedVector & x,
    const double xTol,
    const double gradTol,
    const double scale,
    const unsigned maxIter,
    const bool doOutput);

template
void lbfgs(const Objective3DInterface & objFn,
    stk::math::Vector3d & x,
    const double xTol,
    const double gradTol,
    const double scale,
    const unsigned maxIter,
    const unsigned maxLevels,
    const bool doOutput);

template
void lbfgs(const ObjectiveInterface & objFn,
    DistributedVector & x,
    const double xTol,
    const double gradTol,
    const double scale,
    const unsigned maxIter,
    const unsigned maxLevels,
    const bool doOutput);

}

