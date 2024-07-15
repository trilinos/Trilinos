// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_NONLINEARLEASTSQUARESOBJECTIVE_H
#define ROL_NONLINEARLEASTSQUARESOBJECTIVE_H

#include "ROL_Objective.hpp"
#include "ROL_Constraint.hpp"
#include "ROL_Types.hpp"

/** @ingroup func_group
    \class ROL::NonlinearLeastSquaresObjective
    \brief Provides the interface to evaluate nonlinear least squares objective
           functions.

    ROL's nonlinear least squares objective function interface constructs the
    the nonlinear least squares objective function associated with the equality
    constraint \f$c(x)=0\f$.  That is,
    \f[
       J(x) = \langle \mathfrak{R} c(x),c(x) \rangle_{\mathcal{C}^*,\mathcal{C}}
    \f]
    where \f$c:\mathcal{X}\to\mathcal{C}\f$ and \f$\mathfrak{R}\in\mathcal{L}(
    \mathcal{C},\mathcal{C}^*)\f$ denotes the Riesz map from \f$\mathcal{C}\f$
    into \f$\mathcal{C}^*\f$.

    ---
*/


namespace ROL {

template<typename Real>
class NonlinearLeastSquaresObjective : public Objective<Real> {
private:
  const Ptr<Constraint<Real> > con_;
  const bool GaussNewtonHessian_;

  Ptr<Vector<Real> > c1_, c2_, c1dual_, x_;

public:
  /** \brief Constructor. 

      This function constructs a nonlinear least squares objective function. 
      @param[in]          con   is the nonlinear equation to be solved.
      @param[in]          vec   is a constraint space vector used for cloning.
      @param[in]          GHN   is a flag dictating whether or not to use the Gauss-Newton Hessian.
  */
  NonlinearLeastSquaresObjective(const Ptr<Constraint<Real> > &con,
                                 const Vector<Real> &optvec,
                                 const Vector<Real> &convec,
                                 const bool GNH = false);

  void update( const Vector<Real> &x, UpdateType type, int iter = -1 ) override;
  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) override;
  Real value( const Vector<Real> &x, Real &tol ) override;
  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) override;
  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) override;
  void precond( Vector<Real> &Pv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) override;

// Definitions for parametrized (stochastic) equality constraints
public:
  void setParameter(const std::vector<Real> &param) override;
};

} // namespace ROL

#include "ROL_NonlinearLeastSquaresObjective_Def.hpp"

#endif
