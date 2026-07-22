// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_COMPOSITESTEP_H
#define ROL_COMPOSITESTEP_H

#include "ROL_Types.hpp"
#include "ROL_Step.hpp"
#include "ROL_LAPACK.hpp"
#include "ROL_LinearAlgebra.hpp"
#include <sstream>
#include <iomanip>

/** \class ROL::CompositeStep
    \brief Implements the computation of optimization steps
           with composite-step trust-region methods.
*/


namespace ROL {

template <class Real>
class CompositeStep : public Step<Real> {
private:

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

  template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
  }

  void printInfoLS(const std::vector<Real> &res) const {
    if (infoLS_) {
      std::stringstream hist;
      hist << std::scientific << std::setprecision(8);
      hist << "\n    Augmented System Solver:\n";
      hist << "    True Residual\n";
      for (unsigned j=0; j<res.size(); j++) {
        hist << "    " << std::left << std::setw(14) << res[j] << "\n";
      }
      hist << "\n";
      std::cout << hist.str();
    }
  }

  Real setTolOSS(const Real intol) const {
    return tolOSSfixed_ ? tolOSS_ : intol;
  }

public:
  using Step<Real>::initialize;
  using Step<Real>::compute;
  using Step<Real>::update;  

  virtual ~CompositeStep() {}

  CompositeStep( ROL::ParameterList & parlist ) : Step<Real>() {
    //ROL::Ptr<StepState<Real> > step_state = Step<Real>::getState();
    flagCG_ = 0;
    flagAC_ = 0;
    iterCG_ = 0;

    ROL::ParameterList& steplist = parlist.sublist("Step").sublist("Composite Step");

    //maxiterOSS_  = steplist.sublist("Optimality System Solver").get("Iteration Limit", 50);
    tolOSS_      = steplist.sublist("Optimality System Solver").get("Nominal Relative Tolerance", 1e-8);
    tolOSSfixed_ = steplist.sublist("Optimality System Solver").get("Fix Tolerance", true);
 
    maxiterCG_   = steplist.sublist("Tangential Subproblem Solver").get("Iteration Limit", 20);
    tolCG_       = steplist.sublist("Tangential Subproblem Solver").get("Relative Tolerance", 1e-2);
    Delta_       = steplist.get("Initial Radius", 1e2);
    useConHess_  = steplist.get("Use Constraint Hessian", true);

    int outLvl   = steplist.get("Output Level", 0);

    lmhtol_  = tolOSS_;
    qntol_   = tolOSS_;
    pgtol_   = tolOSS_;
    projtol_ = tolOSS_;
    tangtol_ = tolOSS_;
    tntmax_  = 2.0;

    zeta_    = 0.8;
    penalty_ = 1.0;
    eta_     = 1e-8;

    snorm_   = 0.0;
    nnorm_   = 0.0;
    tnorm_   = 0.0;

    infoALL_  = false;
    if (outLvl > 0) { 
      infoALL_ = true;
    }
    infoQN_  = false;
    infoLM_  = false;
    infoTS_  = false;
    infoAC_  = false;
    infoLS_  = false;
    infoQN_  = infoQN_ || infoALL_;
    infoLM_  = infoLM_ || infoALL_;
    infoTS_  = infoTS_ || infoALL_;
    infoAC_  = infoAC_ || infoALL_;
    infoLS_  = infoLS_ || infoALL_;

    totalIterCG_  = 0;
    totalProj_    = 0;
    totalNegCurv_ = 0;
    totalRef_     = 0;
    totalCallLS_  = 0;
    totalIterLS_  = 0;
  }

  /** \brief Initialize step.
  */
  void initialize( Vector<Real> &x, const Vector<Real> &g, Vector<Real> &l, const Vector<Real> &c,
                   Objective<Real> &obj, Constraint<Real> &con,
                   AlgorithmState<Real> &algo_state ) {
    //ROL::Ptr<StepState<Real> > step_state = Step<Real>::getState();
    ROL::Ptr<StepState<Real> > state = Step<Real>::getState();
    state->descentVec    = x.clone();
    state->gradientVec   = g.clone();
    state->constraintVec = c.clone();

    xvec_ = x.clone();
    gvec_ = g.clone();
    lvec_ = l.clone();
    cvec_ = c.clone();

    ROL::Ptr<Vector<Real> > ajl = gvec_->clone();
    ROL::Ptr<Vector<Real> > gl  = gvec_->clone();

    algo_state.nfval = 0;
    algo_state.ncval = 0;
    algo_state.ngrad = 0;

    Real zerotol = std::sqrt(ROL_EPSILON<Real>());

    // Update objective and constraint.
    obj.update(x,true,algo_state.iter);
    algo_state.value = obj.value(x, zerotol);
    algo_state.nfval++;
    con.update(x,true,algo_state.iter);
    con.value(*cvec_, x, zerotol);
    algo_state.cnorm = cvec_->norm();
    algo_state.ncval++;
    obj.gradient(*gvec_, x, zerotol);

    // Compute gradient of Lagrangian at new multiplier guess.
    computeLagrangeMultiplier(l, x, *gvec_, con);
    con.applyAdjointJacobian(*ajl, l, x, zerotol);
    gl->set(*gvec_); gl->plus(*ajl);
    algo_state.ngrad++;
    algo_state.gnorm = gl->norm();
  }

  /** \brief Compute step.
  */
  void compute( Vector<Real> &s, const Vector<Real> &x, const Vector<Real> &l,
                Objective<Real> &obj, Constraint<Real> &con,
                AlgorithmState<Real> &algo_state ) {
    //ROL::Ptr<StepState<Real> > step_state = Step<Real>::getState();
    Real zerotol = std::sqrt(ROL_EPSILON<Real>());
    Real f = 0.0;
    ROL::Ptr<Vector<Real> > n   = xvec_->clone();
    ROL::Ptr<Vector<Real> > c   = cvec_->clone();
    ROL::Ptr<Vector<Real> > t   = xvec_->clone();
    ROL::Ptr<Vector<Real> > tCP = xvec_->clone();
    ROL::Ptr<Vector<Real> > g   = gvec_->clone();
    ROL::Ptr<Vector<Real> > gf  = gvec_->clone();
    ROL::Ptr<Vector<Real> > Wg  = xvec_->clone();
    ROL::Ptr<Vector<Real> > ajl = gvec_->clone();

    Real f_new = 0.0;
    ROL::Ptr<Vector<Real> > l_new  = lvec_->clone();
    ROL::Ptr<Vector<Real> > c_new  = cvec_->clone();
    ROL::Ptr<Vector<Real> > g_new  = gvec_->clone();
    ROL::Ptr<Vector<Real> > gf_new = gvec_->clone();

    // Evaluate objective ... should have been stored.
    f = obj.value(x, zerotol);
    algo_state.nfval++;
    // Compute gradient of objective ... should have been stored.
    obj.gradient(*gf, x, zerotol);
    // Evaluate constraint ... should have been stored.
    con.value(*c, x, zerotol);

    // Compute quasi-normal step.
    computeQuasinormalStep(*n, *c, x, zeta_*Delta_, con);

    // Compute gradient of Lagrangian ... should have been stored.
    con.applyAdjointJacobian(*ajl, l, x, zerotol);
    g->set(*gf);
    g->plus(*ajl);
    algo_state.ngrad++;

    // Solve tangential subproblem.
    solveTangentialSubproblem(*t, *tCP, *Wg, x, *g, *n, l, Delta_, obj, con);
    totalIterCG_ += iterCG_;

    // Check acceptance of subproblem solutions, adjust merit function penalty parameter, ensure global convergence.
    accept(s, *n, *t, f_new, *c_new, *gf_new, *l_new, *g_new, x, l, f, *gf, *c, *g, *tCP, *Wg, obj, con, algo_state);
  }

  /** \brief Update step, if successful.
  */
  void update( Vector<Real> &x, Vector<Real> &l, const Vector<Real> &s,
               Objective<Real> &obj, Constraint<Real> &con,
               AlgorithmState<Real> &algo_state ) {
    //ROL::Ptr<StepState<Real> > state = Step<Real>::getState();

    Real zero(0);
    Real one(1);
    Real two(2);
    Real seven(7);
    Real half(0.5);
    Real zp9(0.9);
    Real zp8(0.8);
    Real em12(1e-12);
    Real zerotol = std::sqrt(ROL_EPSILON<Real>());//zero;
    Real ratio(zero);

    ROL::Ptr<Vector<Real> > g   = gvec_->clone();
    ROL::Ptr<Vector<Real> > ajl = gvec_->clone();
    ROL::Ptr<Vector<Real> > gl  = gvec_->clone();
    ROL::Ptr<Vector<Real> > c   = cvec_->clone();

    // Determine if the step gives sufficient reduction in the merit function,
    // update the trust-region radius.
    ratio = ared_/pred_;
    if ((std::abs(ared_) < em12) && std::abs(pred_) < em12) {
      ratio = one;
    }
    if (ratio >= eta_) {
      x.plus(s);
      if (ratio >= zp9) {
          Delta_ = std::max(seven*snorm_, Delta_);
      }
      else if (ratio >= zp8) {
          Delta_ = std::max(two*snorm_, Delta_);
      }
      obj.update(x,true,algo_state.iter);
      con.update(x,true,algo_state.iter);
      flagAC_ = 1;
    }
    else {
      Delta_ = half*std::max(nnorm_, tnorm_);
      obj.update(x,false,algo_state.iter);
      con.update(x,false,algo_state.iter);
      flagAC_ = 0;
    } // if (ratio >= eta)

    Real val = obj.value(x, zerotol);
    algo_state.nfval++;
    obj.gradient(*g, x, zerotol);
    computeLagrangeMultiplier(l, x, *g, con);
    con.applyAdjointJacobian(*ajl, l, x, zerotol);
    gl->set(*g); gl->plus(*ajl);
    algo_state.ngrad++;
    con.value(*c, x, zerotol);

    ROL::Ptr<StepState<Real> > state = Step<Real>::getState();
    state->gradientVec->set(*gl);
    state->constraintVec->set(*c);

    algo_state.value = val;
    algo_state.gnorm = gl->norm();
    algo_state.cnorm = c->norm();
    algo_state.iter++;
    algo_state.snorm = snorm_;

    // Update algorithm state
    //(algo_state.iterateVec)->set(x);
  }

  /** \brief Compute step for bound constraints; here only to satisfy the
             interface requirements, does nothing, needs refactoring.
  */
  void compute( Vector<Real> &s, const Vector<Real> &x, Objective<Real> &obj,
                        BoundConstraint<Real> &con,
                        AlgorithmState<Real> &algo_state ) {}

  /** \brief Update step, for bound constraints; here only to satisfy the
             interface requirements, does nothing, needs refactoring.
  */
  void update( Vector<Real> &x, const Vector<Real> &s, Objective<Real> &obj,
                       BoundConstraint<Real> &con,
                       AlgorithmState<Real> &algo_state ) {}

  /** \brief Print iterate header.
  */
  std::string printHeader( void ) const  {
    std::stringstream hist;
    hist << "  ";
    hist << std::setw(6)  << std::left << "iter";
    hist << std::setw(15) << std::left << "fval";
    hist << std::setw(15) << std::left << "cnorm";
    hist << std::setw(15) << std::left << "gLnorm";
    hist << std::setw(15) << std::left << "snorm";
    hist << std::setw(10) << std::left << "delta";
    hist << std::setw(10) << std::left << "nnorm";
    hist << std::setw(10) << std::left << "tnorm";
    hist << std::setw(8) << std::left << "#fval";
    hist << std::setw(8) << std::left << "#grad";
    hist << std::setw(8) << std::left << "iterCG";
    hist << std::setw(8) << std::left << "flagCG";
    hist << std::setw(8) << std::left << "accept";
    hist << std::setw(8) << std::left << "linsys";
    hist << "\n";
    return hist.str();
  }

  std::string printName( void ) const {
    std::stringstream hist;
    hist << "\n" << " Composite-step trust-region solver";
    hist << "\n";
    return hist.str();
  }

  /** \brief Print iterate status.
  */
  std::string print( AlgorithmState<Real> & algo_state, bool pHeader = false ) const  {
    //const ROL::Ptr<const StepState<Real> >& step_state = Step<Real>::getStepState();

    std::stringstream hist;
    hist << std::scientific << std::setprecision(6);
    if ( algo_state.iter == 0 ) {
      hist << printName();
    }
    if ( pHeader ) {
      hist << printHeader();
    }
    if ( algo_state.iter == 0 ) {
      hist << "  ";
      hist << std::setw(6)  << std::left << algo_state.iter;
      hist << std::setw(15) << std::left << algo_state.value;
      hist << std::setw(15) << std::left << algo_state.cnorm;
      hist << std::setw(15) << std::left << algo_state.gnorm;
      hist << "\n";
    }
    else {
      hist << "  ";
      hist << std::setw(6)  << std::left << algo_state.iter;
      hist << std::setw(15) << std::left << algo_state.value;
      hist << std::setw(15) << std::left << algo_state.cnorm;
      hist << std::setw(15) << std::left << algo_state.gnorm;
      hist << std::setw(15) << std::left << algo_state.snorm;
      hist << std::scientific << std::setprecision(2);
      hist << std::setw(10) << std::left << Delta_;
      hist << std::setw(10) << std::left << nnorm_;
      hist << std::setw(10) << std::left << tnorm_;
      hist << std::scientific << std::setprecision(6);
      hist << std::setw(8) << std::left << algo_state.nfval;
      hist << std::setw(8) << std::left << algo_state.ngrad;
      hist << std::setw(8) << std::left << iterCG_;
      hist << std::setw(8) << std::left << flagCG_;
      hist << std::setw(8) << std::left << flagAC_;
      hist << std::left << totalCallLS_ << "/" << totalIterLS_;
      hist << "\n";
    }
    return hist.str();
  }

  /** \brief Compute Lagrange multipliers by solving the least-squares
             problem minimizing the gradient of the Lagrangian, via the
             augmented system formulation.

             @param[out]      l   is the updated Lagrange multiplier; a dual constraint-space vector
             @param[in]       x   is the current iterate; an optimization-space vector
             @param[in]       gf  is the gradient of the objective function; a dual optimization-space vector
             @param[in]       con is the equality constraint object

             On return ... fill the blanks.
  */
  void computeLagrangeMultiplier(Vector<Real> &l, const Vector<Real> &x, const Vector<Real> &gf, Constraint<Real> &con) {

    Real one(1);
    Real zerotol = std::sqrt(ROL_EPSILON<Real>());
    std::vector<Real> augiters;

    if (infoLM_) {
      std::stringstream hist;
      hist << "\n  Lagrange multiplier step\n";
      std::cout << hist.str();
    }

    /* Apply adjoint of constraint Jacobian to current multiplier. */
    ROL::Ptr<Vector<Real> > ajl = gvec_->clone();
    con.applyAdjointJacobian(*ajl, l, x, zerotol);

    /* Form right-hand side of the augmented system. */
    ROL::Ptr<Vector<Real> > b1 = gvec_->clone();
    ROL::Ptr<Vector<Real> > b2 = cvec_->clone();
    // b1 is the negative gradient of the Lagrangian
    b1->set(gf); b1->plus(*ajl); b1->scale(-one);
    // b2 is zero
    b2->zero();

    /* Declare left-hand side of augmented system. */
    ROL::Ptr<Vector<Real> > v1 = xvec_->clone();
    ROL::Ptr<Vector<Real> > v2 = lvec_->clone();

    /* Compute linear solver tolerance. */
    Real b1norm  = b1->norm();
    Real tol = setTolOSS(lmhtol_*b1norm);

    /* Solve augmented system. */
    augiters = con.solveAugmentedSystem(*v1, *v2, *b1, *b2, x, tol);
    totalCallLS_++;
    totalIterLS_ = totalIterLS_ + augiters.size();
    printInfoLS(augiters);

    /* Return updated Lagrange multiplier. */
    // v2 is the multiplier update
    l.plus(*v2);

  }  // computeLagrangeMultiplier


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
  void computeQuasinormalStep(Vector<Real> &n, const Vector<Real> &c, const Vector<Real> &x, Real delta, Constraint<Real> &con) {

    if (infoQN_) {
      std::stringstream hist;
      hist << "\n  Quasi-normal step\n";
      std::cout << hist.str();
    }

    Real zero(0);
    Real one(1);
    Real zerotol = std::sqrt(ROL_EPSILON<Real>()); //zero;
    std::vector<Real> augiters;

    /* Compute Cauchy step nCP. */
    ROL::Ptr<Vector<Real> > nCP     = xvec_->clone();
    ROL::Ptr<Vector<Real> > nCPdual = gvec_->clone();
    ROL::Ptr<Vector<Real> > nN      = xvec_->clone();
    ROL::Ptr<Vector<Real> > ctemp   = cvec_->clone();
    ROL::Ptr<Vector<Real> > dualc0  = lvec_->clone();
    dualc0->set(c.dual());
    con.applyAdjointJacobian(*nCPdual, *dualc0, x, zerotol);
    nCP->set(nCPdual->dual());
    con.applyJacobian(*ctemp, *nCP, x, zerotol);

    Real normsquare_ctemp = ctemp->dot(*ctemp);
    if (normsquare_ctemp != zero) {
      nCP->scale( -(nCP->dot(*nCP))/normsquare_ctemp );
    }

    /* If the  Cauchy step nCP is outside the trust region,
       return the scaled Cauchy step. */
    Real norm_nCP = nCP->norm();
    if (norm_nCP >= delta) {
      n.set(*nCP);
      n.scale( delta/norm_nCP );
      if (infoQN_) {
        std::stringstream hist;
        hist << "  taking partial Cauchy step\n";
        std::cout << hist.str();
      }
      return;
    }

    /* Compute 'Newton' step, for example, by solving a problem
       related to finding the minimum norm solution of min || c(x_k)*s + c ||^2. */
    // Compute tolerance for linear solver.
    con.applyJacobian(*ctemp, *nCP, x, zerotol);
    ctemp->plus(c);
    Real tol = setTolOSS(qntol_*ctemp->norm());
    // Form right-hand side.
    ctemp->scale(-one);
    nCPdual->set(nCP->dual());
    nCPdual->scale(-one);
    // Declare left-hand side of augmented system.
    ROL::Ptr<Vector<Real> > dn = xvec_->clone();
    ROL::Ptr<Vector<Real> > y  = lvec_->clone();
    // Solve augmented system.
    augiters = con.solveAugmentedSystem(*dn, *y, *nCPdual, *ctemp, x, tol);
    totalCallLS_++;
    totalIterLS_ = totalIterLS_ + augiters.size();
    printInfoLS(augiters);

    nN->set(*dn);
    nN->plus(*nCP);

    /* Either take full or partial Newton step, depending on
       the trust-region constraint. */
    Real norm_nN = nN->norm();
    if (norm_nN <= delta) {
      // Take full feasibility step.
      n.set(*nN);
      if (infoQN_) {
        std::stringstream hist;
        hist << "  taking full Newton step\n";
        std::cout << hist.str();
      }
    }
    else {
      // Take convex combination n = nCP+tau*(nN-nCP),
      // so that ||n|| = delta.  In other words, solve
      // scalar quadratic equation: ||nCP+tau*(nN-nCP)||^2 = delta^2.
      Real aa  = dn->dot(*dn);
      Real bb  = dn->dot(*nCP);
      Real cc  = norm_nCP*norm_nCP - delta*delta;
      Real tau = (-bb+sqrt(bb*bb-aa*cc))/aa;
      n.set(*nCP);
      n.axpy(tau, *dn);
      if (infoQN_) {
        std::stringstream hist;
        hist << "  taking dogleg step\n";
        std::cout << hist.str();
      }
    }

  } // computeQuasinormalStep


  /** \brief Solve tangential subproblem.

             @param[out]      t     is the solution of the tangential subproblem; an optimization-space vector
             @param[out]      tCP   is the Cauchy point for the tangential subproblem; an optimization-space vector
             @param[out]      Wg    is the dual of the projected gradient of the Lagrangian; an optimization-space vector
             @param[in]       x     is the current iterate; an optimization-space vector
             @param[in]       g     is the gradient of the Lagrangian; a dual optimization-space vector
             @param[in]       n     is the quasi-normal step; an optimization-space vector
             @param[in]       l     is the Lagrange multiplier; a dual constraint-space vector
             @param[in]       delta is the trust-region radius for the tangential subproblem
             @param[in]       con   is the equality constraint object

  */
  void solveTangentialSubproblem(Vector<Real> &t, Vector<Real> &tCP, Vector<Real> &Wg,
                                 const Vector<Real> &x, const Vector<Real> &g, const Vector<Real> &n, const Vector<Real> &l,
                                 Real delta, Objective<Real> &obj, Constraint<Real> &con) {

    /* Initialization of the CG step. */
    bool orthocheck = true;  // set to true if want to check orthogonality
                             // of Wr and r, otherwise set to false
    Real tol_ortho = 0.5;    // orthogonality measure; represets a bound on norm( \hat{S}, 2), where
                             // \hat{S} is defined in Heinkenschloss/Ridzal., "A Matrix-Free Trust-Region SQP Method"
    Real S_max = 1.0;        // another orthogonality measure; norm(S) needs to be bounded by
                             // a modest constant; norm(S) is small if the approximation of
                             // the null space projector is good
    Real zero(0);
    Real one(1);
    Real zerotol =  std::sqrt(ROL_EPSILON<Real>());
    std::vector<Real> augiters;
    iterCG_ = 1;
    flagCG_ = 0;
    t.zero();
    tCP.zero();
    ROL::Ptr<Vector<Real> > r     = gvec_->clone();
    ROL::Ptr<Vector<Real> > pdesc = xvec_->clone();
    ROL::Ptr<Vector<Real> > tprev = xvec_->clone();
    ROL::Ptr<Vector<Real> > Wr    = xvec_->clone();
    ROL::Ptr<Vector<Real> > Hp    = gvec_->clone();
    ROL::Ptr<Vector<Real> > xtemp = xvec_->clone();
    ROL::Ptr<Vector<Real> > gtemp = gvec_->clone();
    ROL::Ptr<Vector<Real> > ltemp = lvec_->clone();
    ROL::Ptr<Vector<Real> > czero = cvec_->clone();
    czero->zero();
    r->set(g);
    obj.hessVec(*gtemp, n, x, zerotol);
    r->plus(*gtemp);
    if (useConHess_) {
      con.applyAdjointHessian(*gtemp, l, n, x, zerotol);
      r->plus(*gtemp);
    }
    Real normg  = r->norm();
    Real normWg = zero;
    Real pHp    = zero;
    Real rp     = zero;
    Real alpha  = zero;
    Real normp  = zero;
    Real normr  = zero;
    Real normt  = zero;
    std::vector<Real> normWr(maxiterCG_+1, zero);

    std::vector<ROL::Ptr<Vector<Real > > >  p;    // stores search directions
    std::vector<ROL::Ptr<Vector<Real > > >  Hps;  // stores duals of hessvec's applied to p's
    std::vector<ROL::Ptr<Vector<Real > > >  rs;   // stores duals of residuals
    std::vector<ROL::Ptr<Vector<Real > > >  Wrs;  // stores duals of projected residuals

    Real rptol(1e-12);

    if (infoTS_) {
      std::stringstream hist;
      hist << "\n  Tangential subproblem\n";
      hist << std::setw(6)  << std::right << "iter" << std::setw(18) << "||Wr||/||Wr0||" << std::setw(15) << "||s||";
      hist << std::setw(15) << "delta" << std::setw(15) << "||c'(x)s||" << "\n";
      std::cout << hist.str();
    }

    if (normg == 0) {
      if (infoTS_) {
        std::stringstream hist;
        hist << "    >>> Tangential subproblem: Initial gradient is zero! \n";
        std::cout << hist.str();
      }
      iterCG_ = 0; Wg.zero(); flagCG_ = 0;
      return;
    }

    /* Start CG loop. */
    while (iterCG_ < maxiterCG_) {

      // Store tangential Cauchy point (which is the current iterate in the second iteration).
      if (iterCG_ == 2) {
        tCP.set(t);
      }

      // Compute (inexact) projection W*r.
      if (iterCG_ == 1) {
        // Solve augmented system.
        Real tol = setTolOSS(pgtol_);
        augiters = con.solveAugmentedSystem(*Wr, *ltemp, *r, *czero, x, tol);
        totalCallLS_++;
        totalIterLS_ = totalIterLS_ + augiters.size();
        printInfoLS(augiters);

        Wg.set(*Wr);
        normWg = Wg.norm();
        if (orthocheck) {
          Wrs.push_back(xvec_->clone());
          (Wrs[iterCG_-1])->set(*Wr);
        }
        // Check if done (small initial projected residual).
        if (normWg == zero) {
          flagCG_ = 0;
          iterCG_--;
          if (infoTS_) {
            std::stringstream hist;
            hist << "  Initial projected residual is close to zero! \n";
            std::cout << hist.str();
          }
          return;
        }
        // Set first residual to projected gradient.
        // change r->set(Wg);
        r->set(Wg.dual());
        if (orthocheck) {
          rs.push_back(xvec_->clone());
          // change (rs[0])->set(*r);
          (rs[0])->set(r->dual());
        }
      }
      else {
        // Solve augmented system.
        Real tol = setTolOSS(projtol_);
        augiters = con.solveAugmentedSystem(*Wr, *ltemp, *r, *czero, x, tol);
        totalCallLS_++;
        totalIterLS_ = totalIterLS_ + augiters.size();
        printInfoLS(augiters);

        if (orthocheck) {
          Wrs.push_back(xvec_->clone());
          (Wrs[iterCG_-1])->set(*Wr);
        }
      }

      normWr[iterCG_-1] = Wr->norm();

      if (infoTS_) {
        ROL::Ptr<Vector<Real> > ct = cvec_->clone();
        con.applyJacobian(*ct, t, x, zerotol);
        Real linc = ct->norm();
        std::stringstream hist;
        hist << std::scientific << std::setprecision(6);
        hist << std::setw(6)  << std::right << iterCG_-1 << std::setw(18) << normWr[iterCG_-1]/normWg << std::setw(15) << t.norm();
        hist << std::setw(15) << delta << std::setw(15) << linc << "\n";
        std::cout << hist.str();
      }

      // Check if done (small relative residual).
      if (normWr[iterCG_-1]/normWg < tolCG_) {
        flagCG_ = 0;
        iterCG_ = iterCG_-1;
        if (infoTS_) {
          std::stringstream hist;
          hist << "  || W(g + H*(n+s)) || <= cgtol*|| W(g + H*n)|| \n";
          std::cout << hist.str();
        }
        return;
      }

      // Check nonorthogonality, one-norm of (WR*R/diag^2 - I)
      if (orthocheck) {
        ROL::LA::Matrix<Real> Wrr(iterCG_,iterCG_);  // holds matrix Wrs'*rs
        ROL::LA::Matrix<Real> T(iterCG_,iterCG_);    // holds matrix T=(1/diag)*Wrs'*rs*(1/diag)
        ROL::LA::Matrix<Real> Tm1(iterCG_,iterCG_);  // holds matrix Tm1=T-I
        for (int i=0; i<iterCG_; i++) {
          for (int j=0; j<iterCG_; j++) {
            Wrr(i,j)  = (Wrs[i])->dot(*rs[j]);
            T(i,j)    = Wrr(i,j)/(normWr[i]*normWr[j]);
            Tm1(i,j)  = T(i,j);
            if (i==j) {
              Tm1(i,j) = Tm1(i,j) - one;
            }
          }
        }
        if (Tm1.normOne() >= tol_ortho) {
          ROL::LAPACK<int,Real> lapack;
          std::vector<int>          ipiv(iterCG_);
          int                       info;
          std::vector<Real>         work(3*iterCG_);
          // compute inverse of T
          lapack.GETRF(iterCG_, iterCG_, T.values(), T.stride(), &ipiv[0], &info);
          lapack.GETRI(iterCG_, T.values(), T.stride(), &ipiv[0], &work[0], 3*iterCG_, &info);
          Tm1 = T;
          for (int i=0; i<iterCG_; i++) {
            Tm1(i,i) = Tm1(i,i) - one;
          }
          if (Tm1.normOne() > S_max) {
            flagCG_ = 4;
            if (infoTS_) {
              std::stringstream hist;
              hist << "  large nonorthogonality in W(R)'*R detected \n";
              std::cout << hist.str();
            }
            return;
          }
        }
      }

      // Full orthogonalization.
      p.push_back(xvec_->clone());
      (p[iterCG_-1])->set(*Wr);
      (p[iterCG_-1])->scale(-one);
      for (int j=1; j<iterCG_; j++) {
        Real scal = (p[iterCG_-1])->dot(*(Hps[j-1])) / (p[j-1])->dot(*(Hps[j-1]));
        ROL::Ptr<Vector<Real> > pj = xvec_->clone();
        pj->set(*p[j-1]);
        pj->scale(-scal);
        (p[iterCG_-1])->plus(*pj);
      }

      // change Hps.push_back(gvec_->clone());
      Hps.push_back(xvec_->clone());
      // change obj.hessVec(*(Hps[iterCG_-1]), *(p[iterCG_-1]), x, zerotol);
      obj.hessVec(*Hp, *(p[iterCG_-1]), x, zerotol);
      if (useConHess_) {
        con.applyAdjointHessian(*gtemp, l, *(p[iterCG_-1]), x, zerotol);
        // change (Hps[iterCG_-1])->plus(*gtemp);
        Hp->plus(*gtemp);
      }
      // "Preconditioning" step.
      (Hps[iterCG_-1])->set(Hp->dual());

      pHp = (p[iterCG_-1])->dot(*(Hps[iterCG_-1]));
      // change rp  = (p[iterCG_-1])->dot(*r);
      rp  = (p[iterCG_-1])->dot(*(rs[iterCG_-1]));

      normp = (p[iterCG_-1])->norm();
      normr = r->norm();

      // Negative curvature stopping condition.
      if (pHp <= 0) {
        pdesc->set(*(p[iterCG_-1])); // p is the descent direction
        if ((std::abs(rp) >= rptol*normp*normr) && (sgn(rp) == 1)) {
          pdesc->scale(-one); // -p is the descent direction
        }
	flagCG_ = 2;
        Real a = pdesc->dot(*pdesc);
        Real b = pdesc->dot(t);
        Real c = t.dot(t) - delta*delta;
        // Positive root of a*theta^2 + 2*b*theta + c = 0.
        Real theta = (-b + std::sqrt(b*b - a*c)) / a;
        xtemp->set(*(p[iterCG_-1]));
        xtemp->scale(theta);
        t.plus(*xtemp);
        // Store as tangential Cauchy point if terminating in first iteration.
        if (iterCG_ == 1) {
          tCP.set(t);
        }
	if (infoTS_) {
          std::stringstream hist;
          hist << "  negative curvature detected \n";
          std::cout << hist.str();
        }
	return;
      }

      // Want to enforce nonzero alpha's.
      if (std::abs(rp) < rptol*normp*normr) {
        flagCG_ = 5;
        if (infoTS_) {
          std::stringstream hist;
          hist << "  Zero alpha due to inexactness. \n";
          std::cout << hist.str();
        }
	return;
      }

      alpha = - rp/pHp;

      // Iterate update.
      tprev->set(t);
      xtemp->set(*(p[iterCG_-1]));
      xtemp->scale(alpha);
      t.plus(*xtemp);

      // Trust-region stopping condition.
      normt = t.norm();
      if (normt >= delta) {
        pdesc->set(*(p[iterCG_-1])); // p is the descent direction
        if (sgn(rp) == 1) {
          pdesc->scale(-one); // -p is the descent direction
        }
        Real a = pdesc->dot(*pdesc);
        Real b = pdesc->dot(*tprev);
        Real c = tprev->dot(*tprev) - delta*delta;
        // Positive root of a*theta^2 + 2*b*theta + c = 0.
        Real theta = (-b + std::sqrt(b*b - a*c)) / a;
        xtemp->set(*(p[iterCG_-1]));
        xtemp->scale(theta);
        t.set(*tprev);
        t.plus(*xtemp);
        // Store as tangential Cauchy point if terminating in first iteration.
        if (iterCG_ == 1) {
          tCP.set(t);
        }
	flagCG_ = 3;
	if (infoTS_) {
           std::stringstream hist;
           hist << "  trust-region condition active \n";
           std::cout << hist.str();
        }
	return;
      }

      // Residual update.
      xtemp->set(*(Hps[iterCG_-1]));
      xtemp->scale(alpha);
      // change r->plus(*gtemp);
      r->plus(xtemp->dual());
      if (orthocheck) {
        // change rs.push_back(gvec_->clone());
        rs.push_back(xvec_->clone());
        // change (rs[iterCG_])->set(*r);
        (rs[iterCG_])->set(r->dual());
      }

      iterCG_++;

    } // while (iterCG_ < maxiterCG_)

    flagCG_ = 1;
    if (infoTS_) {
      std::stringstream hist;
      hist << "  maximum number of iterations reached \n";
      std::cout << hist.str();
    }

  } // solveTangentialSubproblem


  /** \brief Check acceptance of subproblem solutions, adjust merit function penalty parameter, ensure global convergence.
  */
  void accept(Vector<Real> &s, Vector<Real> &n, Vector<Real> &t, Real f_new, Vector<Real> &c_new,
              Vector<Real> &gf_new, Vector<Real> &l_new, Vector<Real> &g_new,
              const Vector<Real> &x, const Vector<Real> &l, Real f, const Vector<Real> &gf, const Vector<Real> &c,
              const Vector<Real> &g, Vector<Real> &tCP, Vector<Real> &Wg,
              Objective<Real> &obj, Constraint<Real> &con, AlgorithmState<Real> &algo_state) {

    Real beta         = 1e-8;              // predicted reduction parameter
    Real tol_red_tang = 1e-3;              // internal reduction factor for tangtol
    Real tol_red_all  = 1e-1;              // internal reduction factor for qntol, lmhtol, pgtol, projtol, tangtol
    //bool glob_refine  = true;              // true  - if subsolver tolerances are adjusted in this routine, keep adjusted values globally
                                           // false - if subsolver tolerances are adjusted in this routine, discard adjusted values
    Real tol_fdiff    = 1e-12;             // relative objective function difference for ared computation
    int ct_max        = 10;                // maximum number of globalization tries
    Real mintol       = 1e-16;             // smallest projection tolerance value

    // Determines max value of |rpred|/pred.
    Real rpred_over_pred = 0.5*(1-eta_);

    if (infoAC_) {
      std::stringstream hist;
      hist << "\n  Composite step acceptance\n";
      std::cout << hist.str();
    }

    Real zero      =  0.0;
    Real one       =  1.0;
    Real two       =  2.0;
    Real half      =  one/two;
    Real zerotol   =  std::sqrt(ROL_EPSILON<Real>());
    std::vector<Real> augiters;

    Real pred          = zero;
    Real ared          = zero;
    Real rpred         = zero;
    Real part_pred     = zero;
    Real linc_preproj  = zero;
    Real linc_postproj = zero;
    Real tangtol_start = zero;
    Real tangtol = tangtol_;
    //Real projtol = projtol_;
    bool flag = false;
    int num_proj = 0;
    bool try_tCP = false;
    Real fdiff = zero;

    ROL::Ptr<Vector<Real> > xtrial  = xvec_->clone();
    ROL::Ptr<Vector<Real> > Jl      = gvec_->clone();
    ROL::Ptr<Vector<Real> > gfJl    = gvec_->clone();
    ROL::Ptr<Vector<Real> > Jnc     = cvec_->clone();
    ROL::Ptr<Vector<Real> > t_orig  = xvec_->clone();
    ROL::Ptr<Vector<Real> > t_dual  = gvec_->clone();
    ROL::Ptr<Vector<Real> > Jt_orig = cvec_->clone();
    ROL::Ptr<Vector<Real> > t_m_tCP = xvec_->clone();
    ROL::Ptr<Vector<Real> > ltemp   = lvec_->clone();
    ROL::Ptr<Vector<Real> > xtemp   = xvec_->clone();
    ROL::Ptr<Vector<Real> > rt      = cvec_->clone();
    ROL::Ptr<Vector<Real> > Hn      = gvec_->clone();
    ROL::Ptr<Vector<Real> > Hto     = gvec_->clone();
    ROL::Ptr<Vector<Real> > cxxvec  = gvec_->clone();
    ROL::Ptr<Vector<Real> > czero   = cvec_->clone();
    czero->zero();
    Real Jnc_normsquared = zero;
    Real c_normsquared = zero;

    // Compute and store some quantities for later use. Necessary
    // because of the function and constraint updates below.
    con.applyAdjointJacobian(*Jl, l, x, zerotol);
    con.applyJacobian(*Jnc, n, x, zerotol);
    Jnc->plus(c);
    Jnc_normsquared = Jnc->dot(*Jnc);
    c_normsquared = c.dot(c);

    for (int ct=0; ct<ct_max; ct++) {

      try_tCP = true;
      t_m_tCP->set(t);
      t_m_tCP->scale(-one);
      t_m_tCP->plus(tCP);
      if (t_m_tCP->norm() == zero) {
        try_tCP = false;
      }

      t_orig->set(t);
      con.applyJacobian(*Jt_orig, *t_orig, x, zerotol);
      linc_preproj = Jt_orig->norm();
      pred  = one;
      rpred = two*rpred_over_pred*pred;
      flag = false;
      num_proj = 1;
      tangtol_start = tangtol;

      while (std::abs(rpred)/pred > rpred_over_pred) {
        // Compute projected tangential step.
        if (flag) {
          tangtol  = tol_red_tang*tangtol;
          num_proj++;
          if (tangtol < mintol) {
            if (infoAC_) {
              std::stringstream hist;
              hist << "\n The projection of the tangential step cannot be done with sufficient precision.\n";
              hist << " Is the quasi-normal step very small? Continuing with no global convergence guarantees.\n";
              std::cout << hist.str();
            }
            break;
          }
        }
        // Solve augmented system.
        Real tol = setTolOSS(tangtol);
        // change augiters = con.solveAugmentedSystem(t, *ltemp, *t_orig, *czero, x, tol);
        t_dual->set(t_orig->dual());
        augiters = con.solveAugmentedSystem(t, *ltemp, *t_dual, *czero, x, tol);
        totalCallLS_++;
        totalIterLS_ = totalIterLS_ + augiters.size();
        printInfoLS(augiters);
        totalProj_++;
	con.applyJacobian(*rt, t, x, zerotol);
        linc_postproj = rt->norm();

        // Compute composite step.
        s.set(t);
        s.plus(n);

        // Compute some quantities before updating the objective and the constraint.
        obj.hessVec(*Hn, n, x, zerotol);
        if (useConHess_) {
          con.applyAdjointHessian(*cxxvec, l, n, x, zerotol);
          Hn->plus(*cxxvec);
        }
        obj.hessVec(*Hto, *t_orig, x, zerotol);
        if (useConHess_) {
          con.applyAdjointHessian(*cxxvec, l, *t_orig, x, zerotol);
          Hto->plus(*cxxvec);
        }

        // Compute objective, constraint, etc. values at the trial point.
        xtrial->set(x);
        xtrial->plus(s);
        obj.update(*xtrial,false,algo_state.iter);
        con.update(*xtrial,false,algo_state.iter);
        f_new = obj.value(*xtrial, zerotol);
        obj.gradient(gf_new, *xtrial, zerotol);
        con.value(c_new, *xtrial, zerotol);
        l_new.set(l);
        computeLagrangeMultiplier(l_new, *xtrial, gf_new, con);

        // Penalty parameter update.
        part_pred = - Wg.dot(*t_orig);
        gfJl->set(gf);
        gfJl->plus(*Jl);
        // change part_pred -= gfJl->dot(n);
        part_pred -= n.dot(gfJl->dual());
        // change part_pred -= half*Hn->dot(n);
        part_pred -= half*n.dot(Hn->dual());
        // change part_pred -= half*Hto->dot(*t_orig);
        part_pred -= half*t_orig->dot(Hto->dual());
        ltemp->set(l_new);
        ltemp->axpy(-one, l);
        // change part_pred -= Jnc->dot(*ltemp);
        part_pred -= Jnc->dot(ltemp->dual());

        if ( part_pred < -half*penalty_*(c_normsquared-Jnc_normsquared) ) {
          penalty_ = ( -two * part_pred / (c_normsquared-Jnc_normsquared) ) + beta;
        }

        pred = part_pred + penalty_*(c_normsquared-Jnc_normsquared);

        // Computation of rpred.
        // change rpred = - ltemp->dot(*rt) - penalty_ * rt->dot(*rt) - two * penalty_ * rt->dot(*Jnc);
        rpred = - rt->dot(ltemp->dual()) - penalty_ * rt->dot(*rt) - two * penalty_ * rt->dot(*Jnc);
        // change ROL::Ptr<Vector<Real> > lrt   = lvec_->clone();
        //lrt->set(*rt);
        //rpred = - ltemp->dot(*rt) - penalty_ * std::pow(rt->norm(), 2) - two * penalty_ * lrt->dot(*Jnc);
        flag = 1;

      } // while (std::abs(rpred)/pred > rpred_over_pred)

      tangtol = tangtol_start;

      // Check if the solution of the tangential subproblem is
      // disproportionally large compared to total trial step.
      xtemp->set(n);
      xtemp->plus(t);
      if ( t_orig->norm()/xtemp->norm() < tntmax_ ) {
        break;
      }
      else {
        t_m_tCP->set(*t_orig);
        t_m_tCP->scale(-one);
        t_m_tCP->plus(tCP);
        if ((t_m_tCP->norm() > 0) && try_tCP) {
          if (infoAC_) {
            std::stringstream hist;
            hist << "       ---> now trying tangential Cauchy point\n";
            std::cout << hist.str();
          }
          t.set(tCP);
        }
        else {
          if (infoAC_) {
            std::stringstream hist;
            hist << "       ---> recomputing quasi-normal step and re-solving tangential subproblem\n";
            std::cout << hist.str();
          }
          totalRef_++;
          // Reset global quantities.
          obj.update(x);
          con.update(x);
          /*lmhtol  = tol_red_all*lmhtol;
          qntol   = tol_red_all*qntol;
          pgtol   = tol_red_all*pgtol;
          projtol = tol_red_all*projtol;
          tangtol = tol_red_all*tangtol;
          if (glob_refine) {
            lmhtol_  = lmhtol;
            qntol_   = qntol;
            pgtol_   = pgtol;
            projtol_ = projtol;
            tangtol_ = tangtol;
          }*/
          if (!tolOSSfixed_) {
            lmhtol_  *= tol_red_all;
            qntol_   *= tol_red_all;
            pgtol_   *= tol_red_all;
            projtol_ *= tol_red_all;
            tangtol_ *= tol_red_all;
          }
          // Recompute the quasi-normal step.
          computeQuasinormalStep(n, c, x, zeta_*Delta_, con);
          // Solve tangential subproblem.
          solveTangentialSubproblem(t, tCP, Wg, x, g, n, l, Delta_, obj, con);
          totalIterCG_ += iterCG_;
          if (flagCG_ == 1) {
            totalNegCurv_++;
          }
        }
      } // else w.r.t. ( t_orig->norm()/xtemp->norm() < tntmax )

    } // for (int ct=0; ct<ct_max; ct++)

    // Compute actual reduction;
    fdiff = f - f_new;
    // Heuristic 1: If fdiff is very small compared to f, set it to 0,
    // in order to prevent machine precision issues.
    Real em24(1e-24);
    Real em14(1e-14);
    if (std::abs(fdiff / (f+em24)) < tol_fdiff) {
      fdiff = em14;
    }
    // change ared = fdiff  + (l.dot(c) - l_new.dot(c_new)) + penalty_*(c.dot(c) - c_new.dot(c_new));
    // change ared = fdiff  + (l.dot(c) - l_new.dot(c_new)) + penalty_*(std::pow(c.norm(),2) - std::pow(c_new.norm(),2));
    ared = fdiff  + (c.dot(l.dual()) - c_new.dot(l_new.dual())) + penalty_*(c.dot(c) - c_new.dot(c_new));

    // Store actual and predicted reduction.
    ared_ = ared;
    pred_ = pred;

    // Store step and vector norms.
    snorm_ = s.norm();
    nnorm_ = n.norm();
    tnorm_ = t.norm();

    // Print diagnostics.
    if (infoAC_) {
        std::stringstream hist;
        hist << "\n         Trial step info ...\n";
        hist <<   "         n_norm              = " << nnorm_ << "\n";
        hist <<   "         t_norm              = " << tnorm_ << "\n";
        hist <<   "         s_norm              = " << snorm_ << "\n";
        hist <<   "         xtrial_norm         = " << xtrial->norm() << "\n";
        hist <<   "         f_old               = " << f << "\n";
        hist <<   "         f_trial             = " << f_new << "\n";
        hist <<   "         f_old-f_trial       = " << f-f_new << "\n";
        hist <<   "         ||c_old||           = " << c.norm() << "\n";
        hist <<   "         ||c_trial||         = " << c_new.norm() << "\n";
        hist <<   "         ||Jac*t_preproj||   = " << linc_preproj << "\n";
        hist <<   "         ||Jac*t_postproj||  = " << linc_postproj << "\n";
        hist <<   "         ||t_tilde||/||t||   = " << t_orig->norm() / t.norm() << "\n";
        hist <<   "         ||t_tilde||/||n+t|| = " << t_orig->norm() / snorm_ << "\n";
        hist <<   "         # projections       = " << num_proj << "\n";
        hist <<   "         penalty param       = " << penalty_ << "\n";
        hist <<   "         ared                = " << ared_ << "\n";
        hist <<   "         pred                = " << pred_ << "\n";
        hist <<   "         ared/pred           = " << ared_/pred_ << "\n";
        std::cout << hist.str();
    }

  } // accept

}; // class CompositeStep

} // namespace ROL

#endif
