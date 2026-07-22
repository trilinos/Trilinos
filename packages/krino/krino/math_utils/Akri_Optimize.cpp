#include <Akri_DiagWriter.hpp>
#include <Akri_DistributedVector.hpp>
#include <Akri_Optimize.hpp>
#include <krino/mesh_utils/Akri_AllReduce.hpp>
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

//----------------------------------------
// Cubic minimizer in [xLo, xHi]
//   matches phi & dphi at both ends.
// If the root is outside [lo+eps, hi-eps], falls back to midpoint.
//----------------------------------------
double cubic_minimizer(
    const double xLo,  const double xHi,
    const double fLo,  const double fHi,
    const double dfLo, const double dfHi)
{
    // compute coefficients of the interpolating cubic’s derivative
    const double d1 = dfLo + dfHi - 3.0*(fLo - fHi)/(xLo - xHi);
    const double d2sq = d1*d1 - dfLo*dfHi;
    const double d2   = (d2sq > 0.0 ? std::sqrt(d2sq) : 0.0);

    const double numerator   = (xHi - xLo)*(dfHi + d2 - d1);
    const double denominator = (dfHi - dfLo + 2.0*d2);
    double xCubic = xHi - (denominator != 0.0 ? numerator/denominator : 0.5*(xLo + xHi));

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
template<typename VEC>
void zoom(
    const std::function<double(const VEC&)> &objFn,
    const std::function<void(const VEC&, VEC&)> &gradObjFn,
    const VEC &x0, const VEC &p,
    const double f0, const double dphi0,     // phi(0), phi'(0)
    double alphaLo, double alphaHi,
    double phiLo,   double phiHi,
    double dphiLo,  double dphiHi,
    VEC & x,
    double & f,
    VEC & grad,
    const double c1, const double c2,
    const unsigned maxIters = 20)
{
    for (unsigned iter = 0; iter < maxIters; ++iter)
    {
        // trial point by cubic interpolation
        const double alpha_j = cubic_minimizer( alphaLo, alphaHi, phiLo,   phiHi, dphiLo,  dphiHi);

        x = xpby(x0, alpha_j, p);
        f = objFn(x);
        gradObjFn(x, grad);
        const double df = Dot(grad, p);

        // check Armijo:
        if (f > f0 + c1*alpha_j*dphi0 || f >= phiLo)
        {
            // insufficient decrease → shrink right end
            alphaHi = alpha_j;
            phiHi   = f;
            dphiHi  = df;
        }
        else
        {
            // Armijo ok, check curvature:
            if (df >= c2*dphi0)
            {
                return;  // Wolfe satisfied
            }
            // curvature fails
            if (df*(alphaHi - alphaLo) >= 0.0)
            {
                // force bracket to have opposite slope
                alphaHi = alphaLo;
                phiHi   = phiLo;
                dphiHi  = dphiLo;
            }
            // shrink left end
            alphaLo = alpha_j;
            phiLo   = f;
            dphiLo  = df;
        }
    }

    // fallback: midpoint
    const double alphaMid = 0.5*(alphaLo + alphaHi);
    x = xpby(x0, alphaMid, p);
    f = objFn(x);
    gradObjFn(x, grad);
}

template<typename VEC>
void line_search_wolfe_cubic(
    std::function<double(const VEC&)> const &objFn,
    std::function<void(const VEC&, VEC&)> const &gradObjFn,
    const VEC &x0,
    const double f0,
    const VEC &grad0,
    const VEC &p,       // descent dir: grad(x)^T p < 0
    VEC & x,
    double & f,
    VEC & grad,
    const double c1     = 1e-4,  // Armijo param
    const double c2     = 0.9,   // curvature param (Wolfe)
    const int maxZoom = 20)
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

    // evaluate at alpha0
    const double alpha0 = 1.0;
    x = xpby(x0, alpha0, p);
    f = objFn(x);
    gradObjFn(x, grad);
    const double df = Dot(grad, p);

    // if alpha0 already satisfies Armijo and strong Wolfe, accept it
    if (f <= f0 + c1*alpha0*df0 && std::abs(df) <= -c2*df0)
    {
        return;
    }

    // otherwise zoom in [0, alpha0]
    zoom(objFn, gradObjFn, x0, p,
        f0, df0,
        /*alphaLo=*/0.0,/*alphaHi=*/alpha0,
        /*fLo=*/f0,     /*fHi=*/f,
        /*dfLo=*/df0,   /*dfHi=*/df,
        x, f, grad,
        c1, c2, maxZoom);
}

template<typename VEC>
void line_search_armijo(
    const std::function<double(const VEC &)>& objectiveFunction,
    const VEC & x0,
    const double f0,
    const VEC & grad0,
    const VEC & p,
    VEC & x,
    double & f,
    const double c1 = 1.e-4,
    const unsigned maxIters = 25)
{
  const double gradDotP = Dot(grad0, p); // Directional derivative

  double alpha = 1.0; // Start with maximum step size

  for(unsigned iter=0; iter<maxIters; ++iter)
  {
    x = xpby(x0, alpha, p);
    f = objectiveFunction(x);

    // Check Armijo condition (sufficient decrease)
    if (f <= f0 + c1 * alpha * gradDotP)
      return;

    //krinolog << " " << alpha << " " << f << " " << f-(f0 + c1 * alpha * gradDotP) << stk::diag::dendl;

    alpha *= 0.5; // Reduce step size
  }

  x = x0; // Fallback to 0 step
  f = f0;
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
template<typename VEC>
void lbfgs(const std::function<double(const VEC&)> & calc_objective,
    const std::function<void(const VEC&, VEC&)> & fill_gradient,
    VEC& x,
    const double xTol,
    const double gradTol,
    const unsigned maxIter,
    const unsigned maxLevels)
{
    std::deque<VEC> sList, yList;
    std::deque<double> rhoList;

    double f = calc_objective(x);
    VEC grad;
    fill_gradient(x, grad);

    double gradMag = 0.;
    double dxMag = 0.;

    for (unsigned iter = 0; iter < maxIter; ++iter)
    {
        gradMag = std::sqrt(Dot(grad, grad));
        if (gradMag < gradTol)
        {
            //krinolog << "Gradient converged at iteration " << iter << ", |grad|= " << gradMag << stk::diag::dendl;
            return;
        }

        // compute search direction p = - H_k * grad
        const VEC Hgrad = compute_LBFGS_Hgrad(grad, sList, yList, rhoList);
        const VEC p = scalar_times_vector(-1., Hgrad);

        const auto xOld = x;
        const auto gradOld = grad;
        const double fOld = f;
        line_search_wolfe_cubic(calc_objective, fill_gradient, xOld, fOld, gradOld, p, x, f, grad);

        const VEC s = vectorSubtract(x, xOld);
        dxMag = std::sqrt(Dot(s,s));
        //krinolog << "L-BFGS iteration " << iter << ", |grad|= " << gradMag << ", |dx|= " << dxMag << stk::diag::dendl;
        if (dxMag < xTol)
        {
          //krinolog << "Solution converged at iteration " << iter << ", |dx|= " << dxMag << stk::diag::dendl;
          return;
        }

        const auto y = vectorSubtract(grad, gradOld);
        update_LBFGS_Memory(sList, yList, rhoList, s, y, maxLevels);
    }

    krinolog << "Reached max iterations " << maxIter << ", |grad|= " << gradMag << ", |dx|= " << dxMag << stk::diag::dendl;
    return;
}

template<typename VEC>
void steepest_descent(const std::function<double(const VEC&)> & calc_objective,
    const std::function<void(const VEC&, VEC&)> & fill_gradient,
    VEC& x,
    const double xTol,
    const double gradTol,
    const unsigned maxIter)
{
  VEC grad;
  fill_gradient(x, grad);
  double f = calc_objective(x);

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
      //krinolog << "Gradient converged at iteration " << iter << ", |grad|= " << gradMag << stk::diag::dendl;
      return;
    }

    // Compute the search direction
    const VEC p = scalar_times_vector(-scaling, grad);

    const auto xOld = x;
    const auto gradOld = grad;
    const double fOld = f;
    line_search_armijo(calc_objective, xOld, fOld, gradOld, p, x, f);

    const VEC s = vectorSubtract(x, xOld);
    dxSqr = Dot(s,s);
    //krinolog << "Steepest descent iteration " << iter << ", |grad|= " << gradMag << ", |dx|= " << std::sqrt(dxSqr) << stk::diag::dendl;
    if (std::sqrt(dxSqr) < xTol)
    {
      //krinolog << "Solution converged at iteration " << iter << ", |dx|= " << std::sqrt(dxSqr) << stk::diag::dendl;
      return;
    }

    fill_gradient(x, grad);
    const VEC y = vectorSubtract(grad, gradOld);
    const double ys = Dot(y, s);
    if (ys > 0.)
      scaling = std::min(1.e6,std::max(1.e-6,dxSqr/ys)); //Barzilai–Borwein
  }

  krinolog << "Reached max iterations " << maxIter << ", |grad|= " << std::sqrt(Dot(grad,grad)) << ", |dx|= " << std::sqrt(dxSqr) << stk::diag::dendl;
}

template
void steepest_descent(const Vector3dObjectiveFn & calc_objective,
    const Vector3dObjectiveSensFn & fill_gradient,
    stk::math::Vector3d& x,
    const double xTol,
    const double gradTol,
    const unsigned maxIter);

template
void steepest_descent(const DistributedVectorObjectiveFn & calc_objective,
    const DistributedVectorObjectiveSensFn & fill_gradient,
    DistributedVector& x,
    const double xTol,
    const double gradTol,
    const unsigned maxIter);

template
void lbfgs(const Vector3dObjectiveFn & calc_objective,
    const Vector3dObjectiveSensFn & fill_gradient,
    stk::math::Vector3d& x,
    const double xTol,
    const double gradTol,
    const unsigned maxIter,
    const unsigned maxLevels);

template
void lbfgs(const DistributedVectorObjectiveFn & calc_objective,
    const DistributedVectorObjectiveSensFn & fill_gradient,
    DistributedVector& x,
    const double xTol,
    const double gradTol,
    const unsigned maxIter,
    const unsigned maxLevels);

}

