// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#ifndef ROL_TYPEE_COMPOSITESTEPALGORITHM_H
#define ROL_TYPEE_COMPOSITESTEPALGORITHM_H

#include "ROL_TypeE_Algorithm.hpp"

/** \class ROL::TypeE::CompositeStepAlgorithm
    \brief Provides an interface to run equality constrained
           optimization algorithms using the Composite-Step
           Trust-Region Sequential Quadratic Programming (SQP) method.
*/

namespace ROL {
namespace TypeE {

template<typename Real>
class CompositeStepAlgorithm : public TypeE::Algorithm<Real> {
private:
  // Parameter list.
  ParameterList list_;

  // Vectors used for cloning.
  ROL::Ptr<Vector<Real> > xvec_;
  ROL::Ptr<Vector<Real> > gvec_;
  ROL::Ptr<Vector<Real> > cvec_;
  ROL::Ptr<Vector<Real> > lvec_;

  // Diagnostic return flags for subalgorithms.
  int flagCG_;
  int flagAC_;
  int iterCG_;

  // Stopping conditions.
  int maxiterCG_;
  int maxiterOSS_;
  Real tolCG_;
  Real tolOSS_;
  bool tolOSSfixed_;

  // Tolerances and stopping conditions for subalgorithms.
  Real lmhtol_;
  Real qntol_;
  Real pgtol_;
  Real projtol_;
  Real tangtol_;
  Real tntmax_;

  // Trust-region parameters.
  Real zeta_;
  Real Delta_;
  Real penalty_;
  Real eta_;
  bool useConHess_;

  Real ared_;
  Real pred_;
  Real snorm_;
  Real nnorm_;
  Real tnorm_;

  // Output flags.
  bool infoQN_;
  bool infoLM_;
  bool infoTS_;
  bool infoAC_;
  bool infoLS_;
  bool infoALL_;

  // Performance summary.
  int totalIterCG_;
  int totalProj_;
  int totalNegCurv_;
  int totalRef_;
  int totalCallLS_;
  int totalIterLS_;

  // Verbosity flags.
  int verbosity_;
  bool printHeader_;

  using TypeE::Algorithm<Real>::state_;
  using TypeE::Algorithm<Real>::status_;

  /** \brief Initialize algorithm by computing a few quantities.
  */
  void initialize(Vector<Real>        &x,
                  const Vector<Real>  &g,
                  Vector<Real>        &l,
                  const Vector<Real>  &c,
                  Objective<Real>     &obj,
                  Constraint<Real>    &con,
                  std::ostream        &outStream = std::cout);

  /** \brief Compute trial step.
  */
  void computeTrial(Vector<Real> &s,
                    const Vector<Real> &x,
                    const Vector<Real> &l,
                    Objective<Real>    &obj,
                    Constraint<Real>   &con,
                    std::ostream       &os);

  /** \brief Update trust-region radius, take step, etc.
  */
  void updateRadius(Vector<Real>       &x,
                    Vector<Real>       &l,
                    const Vector<Real> &s,
                    Objective<Real>    &obj,
                    Constraint<Real>   &con,
                    std::ostream       &os);

  /** \brief Compute Lagrange multipliers by solving the least-squares
             problem minimizing the gradient of the Lagrangian, via the
             augmented system formulation.

             @param[out]      l   is the updated Lagrange multiplier; a dual constraint-space vector
             @param[in]       x   is the current iterate; an optimization-space vector
             @param[in]       gf  is the gradient of the objective function; a dual optimization-space vector
             @param[in]       con is the equality constraint object
  */
  void computeLagrangeMultiplier(Vector<Real>       &l,
                                 const Vector<Real> &x,
                                 const Vector<Real> &gf,
                                 Constraint<Real>   &con,
                                 std::ostream       &os);

  /** \brief Compute quasi-normal step by minimizing the norm of
             the linearized constraint.

             Compute an approximate solution of the problem
             \f[
               \begin{array}{rl}
                 \min_{n} & \|c'(x_k)n + c(x_k)\|^2_{\mathcal{X}} \\
                 \mbox{subject to} & \|n\|_{\mathcal{X}} \le \delta
               \end{array}
             \f]
             The approximate solution is computed using Powell's dogleg
             method. The dogleg path is computed using the Cauchy point and
             a full Newton step.  The path's intersection with the trust-region
             constraint gives the quasi-normal step.

             @param[out]      n     is the quasi-normal step; an optimization-space vector
             @param[in]       c     is the value of equality constraints; a constraint-space vector
             @param[in]       x     is the current iterate; an optimization-space vector
             @param[in]       delta is the trust-region radius for the quasi-normal step
             @param[in]       con   is the equality constraint object

  */
  void computeQuasinormalStep(Vector<Real>       &n,
                              const Vector<Real> &c,
                              const Vector<Real> &x,
                              Real               delta,
                              Constraint<Real>   &con,
                              std::ostream       &os);

  /** \brief Solve tangential subproblem.

             @param[out]      t     is the solution of the tangential subproblem; an optimization-space vector
             @param[out]      tCP   is the Cauchy point for the tangential subproblem; an optimization-space vector
             @param[out]      Wg    is the dual of the projected gradient of the Lagrangian; an optimization-space vector
             @param[in]       x     is the current iterate; an optimization-space vector
             @param[in]       g     is the gradient of the Lagrangian; a dual optimization-space vector
             @param[in]       n     is the quasi-normal step; an optimization-space vector
             @param[in]       l     is the Lagrange multiplier; a dual constraint-space vector
             @param[in]       delta is the trust-region radius for the tangential subproblem
             @param[in]       obj   is the objective function
             @param[in]       con   is the equality constraint object

  */
  void solveTangentialSubproblem(Vector<Real>       &t,
                                 Vector<Real>       &tCP,
                                 Vector<Real>       &Wg,
                                 const Vector<Real> &x,
                                 const Vector<Real> &g,
                                 const Vector<Real> &n,
                                 const Vector<Real> &l,
                                 Real               delta,
                                 Objective<Real>    &obj,
                                 Constraint<Real>   &con,
                                 std::ostream       &os);

  /** \brief Check acceptance of subproblem solutions, adjust merit function penalty parameter, ensure global convergence.
  */
  void accept(Vector<Real> &s, Vector<Real> &n, Vector<Real> &t, Real f_new, Vector<Real> &c_new,
              Vector<Real> &gf_new, Vector<Real> &l_new, Vector<Real> &g_new,
              const Vector<Real> &x, const Vector<Real> &l, Real f, const Vector<Real> &gf, const Vector<Real> &c,
              const Vector<Real> &g, Vector<Real> &tCP, Vector<Real> &Wg,
              Objective<Real> &obj, Constraint<Real> &con, std::ostream &os);

  template<typename T> int sgn(T val) const;

  void printInfoLS(const std::vector<Real> &res, std::ostream& os) const;

  Real setTolOSS(const Real intol) const;


public:

  CompositeStepAlgorithm(ParameterList &list);

  using TypeE::Algorithm<Real>::run;
  virtual void run(Vector<Real>       &x,
                   const Vector<Real> &g,
                   Objective<Real>    &obj,
                   Constraint<Real>   &econ,
                   Vector<Real>       &emul,
                   const Vector<Real> &eres,
                   std::ostream       &outStream = std::cout) override;

  virtual void writeHeader(std::ostream& os) const override;

  virtual void writeName(std::ostream& os) const override;

  virtual void writeOutput(std::ostream& os, const bool print_header = false) const override;

}; // class ROL::TypeE::CompositeStepAlgorithm

} // namespace TypeE
} // namespace ROL

#include "ROL_TypeE_CompositeStepAlgorithm_Def.hpp"

#endif
