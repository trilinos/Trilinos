// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef ROL_COMPOSITESTEPSQP_H
#define ROL_COMPOSITESTEPSQP_H

#include "ROL_Types.hpp"
#include "ROL_Step.hpp"
#include <sstream>
#include <iomanip>
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_LAPACK.hpp"

/** \class ROL::CompositeStepSQP
    \brief Implements the computation of optimization steps
           with composite-step trust-region SQP methods.
*/


namespace ROL {

template <class Real>
class CompositeStepSQP : public Step<Real> {
private:

  // Diagnostic return flags for subalgorithms. 
  int flagTR_;
  int flagCG_;
  int iterCG_;

  // Stopping conditions.
  int maxiterCG_;
  Real tolCG_;

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
  bool infoALL_;

  // Performance summary.
  int totalIterCG_;
  int totalProj_;
  int totalNegCurv_;
  int totalRef_;

  template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
  }

public:

  virtual ~CompositeStepSQP() {}

  CompositeStepSQP( Teuchos::ParameterList & parlist ) : Step<Real>() {
    //Teuchos::RCP<StepState<Real> > step_state = Step<Real>::getState();
    flagTR_ = 0;
    flagCG_ = 0;
    iterCG_ = 0;

    maxiterCG_ = 20;
    tolCG_ = 1e-4;

    Real nominal_tol = parlist.get("Nominal SQP Optimality Solver Tolerance", 1e-3);
    lmhtol_  = nominal_tol;
    qntol_   = nominal_tol;
    pgtol_   = nominal_tol;       
    projtol_ = nominal_tol;     
    tangtol_ = nominal_tol;
    tntmax_  = 2.0;
    
    zeta_    = 0.9;
    Delta_   = 1e2;
    penalty_ = 1.0;
    eta_     = 1e-8;

    snorm_   = 0.0;
    nnorm_   = 0.0;
    tnorm_   = 0.0;

    infoQN_  = false;
    infoLM_  = false;
    infoTS_  = false;
    infoAC_  = false;
    infoALL_ = true;
    infoQN_  = infoQN_ || infoALL_;
    infoLM_  = infoLM_ || infoALL_;
    infoTS_  = infoTS_ || infoALL_;
    infoAC_  = infoAC_ || infoALL_;

    totalIterCG_  = 0;
    totalProj_    = 0;
    totalNegCurv_ = 0;
    totalRef_ = 0;
  }

  /** \brief Initialize step.
  */
  void initialize( Vector<Real> &x, Vector<Real> &l,
                   Objective<Real> &obj, EqualityConstraint<Real> &con, 
                   AlgorithmState<Real> &algo_state ) {
    //Teuchos::RCP<StepState<Real> > step_state = Step<Real>::getState();

    algo_state.nfval = 0;
    algo_state.ncval = 0;
    algo_state.ngrad = 0;

    Real zerotol = 0.0;

    // Update objective and constraint.
    obj.update(x,true,algo_state.iter);    
    algo_state.value = obj.value(x, zerotol); 
    algo_state.nfval++;
    con.update(x,true,algo_state.iter);
    Teuchos::RCP<Vector<Real> > c = l.clone();
    con.value(*c, x, zerotol);    
    algo_state.cnorm = c->norm(); 
    algo_state.ncval++;
    Teuchos::RCP<Vector<Real> > g = x.clone();
    obj.gradient(*g, x, zerotol);

    // Compute gradient of Lagrangian at new multiplier guess.
    computeLagrangeMultiplier(l, x, *g, con);
    Teuchos::RCP<Vector<Real> > ajl = x.clone();
    con.applyAdjointJacobian(*ajl, l, x, zerotol);
    Teuchos::RCP<Vector<Real> > gl = x.clone();
    gl->set(*g); gl->plus(*ajl);
    algo_state.ngrad++;
    algo_state.gnorm = gl->norm();
  }

  /** \brief Compute step.
  */
  void compute( Vector<Real> &s, const Vector<Real> &x, const Vector<Real> &l,
                Objective<Real> &obj, EqualityConstraint<Real> &con, 
                AlgorithmState<Real> &algo_state ) {
    //Teuchos::RCP<StepState<Real> > step_state = Step<Real>::getState();
    Real zerotol = 0.0;
    Real f = 0.0;
    Teuchos::RCP<Vector<Real> > n   = s.clone();
    Teuchos::RCP<Vector<Real> > c   = l.clone();
    Teuchos::RCP<Vector<Real> > t   = s.clone();
    Teuchos::RCP<Vector<Real> > tCP = s.clone();
    Teuchos::RCP<Vector<Real> > g   = x.clone();
    Teuchos::RCP<Vector<Real> > gf  = x.clone();
    Teuchos::RCP<Vector<Real> > Wg  = x.clone();
    Teuchos::RCP<Vector<Real> > ajl = x.clone();

    Real f_new = 0.0;
    Teuchos::RCP<Vector<Real> > l_new  = l.clone();
    Teuchos::RCP<Vector<Real> > c_new  = l.clone();
    Teuchos::RCP<Vector<Real> > g_new  = x.clone();
    Teuchos::RCP<Vector<Real> > gf_new = x.clone();

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

    //s.set(*n);
    //s.plus(*t);
  }

  /** \brief Update step, if successful.
  */
  void update( Vector<Real> &x, Vector<Real> &l, const Vector<Real> &s,
               Objective<Real> &obj, EqualityConstraint<Real> &con, 
               AlgorithmState<Real> &algo_state ) {
    //Teuchos::RCP<StepState<Real> > state = Step<Real>::getState();

    Real zero = 0.0;
    Real half = 0.5;
    Real zerotol = zero;
    Real ratio = zero;

    Teuchos::RCP<Vector<Real> > g   = x.clone();
    Teuchos::RCP<Vector<Real> > ajl = x.clone();
    Teuchos::RCP<Vector<Real> > gl = x.clone();
    Teuchos::RCP<Vector<Real> > c = l.clone();

    // Determine if the step gives sufficient reduction in the merit function,
    // update the trust-region radius.
    ratio = ared_/pred_;
    if (ratio >= eta_) {
      x.plus(s);
      if (ratio >= 0.9) {
          Delta_ = std::max(7.0*snorm_, Delta_);
      }
      else if (ratio >= 0.8) {
          Delta_ = std::max(2.0*snorm_, Delta_);
      }
      obj.update(x,true,algo_state.iter);
      con.update(x,true,algo_state.iter);
    }
    else {
      Delta_ = half*std::max(nnorm_, tnorm_);
      obj.update(x,false,algo_state.iter);
      con.update(x,false,algo_state.iter);
    } // if (ratio >= eta) 

    Real val = obj.value(x, zerotol);
    algo_state.nfval++;
    obj.gradient(*g, x, zerotol);
    computeLagrangeMultiplier(l, x, *g, con);
    con.applyAdjointJacobian(*ajl, l, x, zerotol);
    gl->set(*g); gl->plus(*ajl);
    algo_state.ngrad++;
    con.value(*c, x, zerotol);

    algo_state.value = val;
    algo_state.gnorm = gl->norm();
    algo_state.cnorm = c->norm();
    algo_state.iter++;
    algo_state.snorm = snorm_;

    //Real tol = std::sqrt(ROL_EPSILON);

    //Real eps = 0.0;
  
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
    hist << std::setw(15) << std::left << "delta";
    hist << std::setw(10) << std::left << "#fval";
    hist << std::setw(10) << std::left << "#grad";
    hist << std::setw(10) << std::left << "flagTR";
    hist << std::setw(10) << std::left << "iterCG";
    hist << std::setw(10) << std::left << "flagCG";
    hist << "\n";
    return hist.str();
  }

  std::string printName( void ) const {
    std::stringstream hist;
    hist << "\n" << " Composite-step trust-region SQP solver";
    hist << "\n";
    return hist.str();
  }

  /** \brief Print iterate status.
  */
  std::string print( AlgorithmState<Real> & algo_state, bool pHeader = false ) const  {
    //const Teuchos::RCP<const StepState<Real> >& step_state = Step<Real>::getStepState();

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
      hist << std::setw(15) << std::left << Delta_; 
      hist << std::setw(10) << std::left << algo_state.nfval;              
      hist << std::setw(10) << std::left << algo_state.ngrad;              
      hist << std::setw(10) << std::left << flagTR_;              
      hist << std::setw(10) << std::left << iterCG_;
      hist << std::setw(10) << std::left << flagCG_;
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
  void computeLagrangeMultiplier(Vector<Real> &l, const Vector<Real> &x, const Vector<Real> &gf, EqualityConstraint<Real> &con) {

    Real zerotol = 0.0;

    if (infoLM_) {
      std::stringstream hist;
      hist << "\n  SQP_lagrange_multiplier\n";
      std::cout << hist.str();
    }

    /* Apply adjoint of constraint Jacobian to current multiplier. */
    Teuchos::RCP<Vector<Real> > ajl = x.clone();
    con.applyAdjointJacobian(*ajl, l, x, zerotol);

    /* Form right-hand side of the augmented system. */
    Teuchos::RCP<Vector<Real> > b1 = x.clone();
    Teuchos::RCP<Vector<Real> > b2 = l.clone();
    // b1 is the negative gradient of the Lagrangian
    b1->set(gf); b1->plus(*ajl); b1->scale(-1.0);
    // b2 is zero
    b2->zero();

    /* Declare left-hand side of augmented system. */
    Teuchos::RCP<Vector<Real> > v1 = x.clone();
    Teuchos::RCP<Vector<Real> > v2 = l.clone();
    
    /* Compute linear solver tolerance. */
    Real b1norm  = b1->norm();
    Real tol = lmhtol_*b1norm;

    /* Solve augmented system. */
    con.solveAugmentedSystem(*v1, *v2, *b1, *b2, x, tol);

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
  void computeQuasinormalStep(Vector<Real> &n, const Vector<Real> &c, const Vector<Real> &x, Real delta, EqualityConstraint<Real> &con) {

    if (infoQN_) {
      std::stringstream hist;
      hist << "\n  SQP_quasi-normal_step\n";
      std::cout << hist.str();
    }

    Real zero    = 0.0;
    Real one     = 1.0;
    Real zerotol = zero;

    /* Compute Cauchy step nCP. */
    Teuchos::RCP<Vector<Real> > nCP     = n.clone();
    Teuchos::RCP<Vector<Real> > nCPtemp = n.clone();
    Teuchos::RCP<Vector<Real> > nN      = n.clone();
    Teuchos::RCP<Vector<Real> > ctemp   = c.clone();
    con.applyAdjointJacobian(*nCP, c, x, zerotol);
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
        std::cout << "  taking partial Cauchy step\n";
      }
      return;
    }

    /* Compute 'Newton' step, for example, by solving a problem
       related to finding the minimum norm solution of min || c(x_k)*s + c ||^2. */
    // Compute tolerance for linear solver.
    con.applyJacobian(*ctemp, *nCP, x, zerotol);
    ctemp->plus(c);
    Real tol = qntol_*ctemp->norm();
    // Form right-hand side.
    ctemp->scale(-one);
    nCPtemp->set(*nCP);
    nCPtemp->scale(-one);
    // Declare left-hand side of augmented system.
    Teuchos::RCP<Vector<Real> > dn = n.clone();
    Teuchos::RCP<Vector<Real> > y  = c.clone();
    // Solve augmented system.
    con.solveAugmentedSystem(*dn, *y, *nCPtemp, *ctemp, x, tol);
    nN->set(*dn);
    nN->plus(*nCP);

    /* Either take full or partial Newton step, depending on
       the trust-region constraint. */
    Real norm_nN = nN->norm();
    if (norm_nN <= delta) {
      // Take full feasibility step.
      n.set(*nN);
      if (infoQN_) {
        std::cout << "  taking full Newton step\n";
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
      n.set(*dn);
      n.axpy(tau, *nCP);
      if (infoQN_) {
        std::cout << "  taking dogleg step\n";
      }
    }

  } // computeQuasinormalStep


  /** \brief Solve tangential subproblem.

             @param[out]      t     is the solution of the tangential subproblem; an optimization-space vector
             @param[out]      tCP   is the Cauchy point for the tangential subproblem; an optimization-space vector
             @param[out]      Wg    is the projected gradient of the Lagrangian; a dual optimization-space vector
             @param[in]       x     is the current iterate; an optimization-space vector
             @param[in]       g     is the gradient of the Lagrangian; a dual optimization-space vector
             @param[in]       n     is the quasi-normal step; an optimization-space vector
             @param[in]       l     is the Lagrange multiplier; a dual constraint-space vector
             @param[in]       delta is the trust-region radius for the tangential subproblem
             @param[in]       con   is the equality constraint object

  */
  void solveTangentialSubproblem(Vector<Real> &t, Vector<Real> &tCP, Vector<Real> &Wg,
                                 const Vector<Real> &x, const Vector<Real> &g, const Vector<Real> &n, const Vector<Real> &l,
                                 Real delta, Objective<Real> &obj, EqualityConstraint<Real> &con) {

    /* Initialization of the CG step. */
    bool orthocheck = true;  // set to true if want to check orthogonality
                             // of Wr and r, otherwise set to false
    Real tol_ortho  = 0.5;   // orthogonality measure; represets a bound on norm( \hat{S}, 2), where
                             // \hat{S} is defined in Heinkenschloss/Ridzal., "A Matrix-Free Trust-Region SQP Method"
    Real S_max      = 1.0;   // another orthogonality measure; norm(S) needs to be bounded by
                             // a modest constant; norm(S) is small if the approximation of
                             // the null space projector is good
    Real zero    =  0.0;
    Real one     =  1.0;
    Real zerotol =  zero;
    iterCG_ = 1;
    flagCG_ = 0;
    t.zero();
    tCP.zero();
    Teuchos::RCP<Vector<Real> > r = g.clone();
    Teuchos::RCP<Vector<Real> > pdesc = g.clone();
    Teuchos::RCP<Vector<Real> > tprev = t.clone();
    Teuchos::RCP<Vector<Real> > Wr = g.clone();
    Teuchos::RCP<Vector<Real> > vtemp = g.clone();
    Teuchos::RCP<Vector<Real> > ltemp = l.clone();
    Teuchos::RCP<Vector<Real> > czero = l.clone();
    czero->zero();
    r->set(g);
    obj.hessVec(*vtemp, n, x, zerotol);
    r->plus(*vtemp);
    con.applyAdjointHessian(*vtemp, l, n, x, zerotol);
    r->plus(*vtemp);
    Real normg  = r->norm();
    Real normWg = zero;
    Real pHp    = zero;
    Real rp     = zero;
    Real alpha  = zero;
    Real normp  = zero;
    Real normr  = zero;
    Real normt  = zero;
    std::vector<Real> normWr(maxiterCG_+1, zero);

    std::vector<Teuchos::RCP<Vector<Real > > >  p;   // stores search directions
    std::vector<Teuchos::RCP<Vector<Real > > >  Hp;  // stores hessvec's applied to p's
    std::vector<Teuchos::RCP<Vector<Real > > >  rs;  // stores residuals
    std::vector<Teuchos::RCP<Vector<Real > > >  Wrs; // stores projected residuals

    Real rptol = 1e-12;

    if (infoTS_) {
      std::stringstream hist;
      hist << "\n  SQP_tangential_subproblem\n";
      hist << std::setw(6)  << std::right << "iter" << std::setw(18) << "||Wr||/||Wr0||" << std::setw(15) << "||s||";
      hist << std::setw(15) << "delta" << std::setw(15) << "||c'(x)s||" << "\n";
      std::cout << hist.str();
    }

    if (normg == 0) {
      if (infoTS_) {
        std::stringstream hist;
        hist << "    >>> Tangential subproblem: Initial gradient is zero! \n";
        std::cout << hist;
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
        Real tol = pgtol_;
        con.solveAugmentedSystem(*Wr, *ltemp, *r, *czero, x, tol);
        Wg.set(*Wr);
        normWg = Wg.norm();
        if (orthocheck) {
          Wrs.push_back(Wr->clone());
          (Wrs[iterCG_-1])->set(*Wr);
        }
        // Check if done (small initial projected residual).
        if (normWg < 1e-14) {
          flagCG_ = 0;
          iterCG_ = iterCG_-1;
          if (infoTS_) {
            std::stringstream hist;
            hist << "  Initial projected residual is close to zero! \n";
            std::cout << hist;
          }
          return;
        }
        // Set first residual to projected gradient.
        r->set(Wg);
        if (orthocheck) {
          rs.push_back(r->clone());
          (rs[0])->set(*r);
        }
      }
      else {
        // Solve augmented system.
        Real tol = projtol_;
        con.solveAugmentedSystem(*Wr, *ltemp, *r, *czero, x, tol);
        if (orthocheck) {
          Wrs.push_back(Wr->clone());
          (Wrs[iterCG_-1])->set(*Wr);
        }
      }

      normWr[iterCG_-1] = Wr->norm();

      if (infoTS_) {
        Teuchos::RCP<Vector<Real> > ct = l.clone();
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
          std::cout << "  || W(g + H*(n+s)) || <= cgtol*|| W(g + H*n)|| \n";
        }
        return;
      }

      // Check nonorthogonality, one-norm of (WR*R/diag^2 - I)
      if (orthocheck) {
        Teuchos::SerialDenseMatrix<int,Real> Wrr(iterCG_,iterCG_);  // holds matrix Wrs'*rs
        Teuchos::SerialDenseMatrix<int,Real> T(iterCG_,iterCG_);    // holds matrix T=(1/diag)*Wrs'*rs*(1/diag)
        Teuchos::SerialDenseMatrix<int,Real> Tm1(iterCG_,iterCG_);  // holds matrix Tm1=T-I
        for (int i=0; i<iterCG_; i++) {
          for (int j=0; j<iterCG_; j++) {
            Wrr(i,j)  = (Wrs[i])->dot(*rs[j]);
            T(i,j)    = Wrr(i,j)/(normWr[i]*normWr[j]);
            Tm1(i,j)  = T(i,j);
            if (i==j) {
              Tm1(i,j) = Tm1(i,j) - 1.0;
            }
          }
        }
        if (Tm1.normOne() >= tol_ortho) {
          Teuchos::LAPACK<int,Real> lapack;
          std::vector<int>          ipiv(iterCG_);
          int                       info;
          std::vector<Real>         work(3*iterCG_);
          // compute inverse of T
          lapack.GETRF(iterCG_, iterCG_, T.values(), T.stride(), &ipiv[0], &info);
          lapack.GETRI(iterCG_, T.values(), T.stride(), &ipiv[0], &work[0], 3*iterCG_, &info);
          Tm1 = T;
          for (int i=0; i<iterCG_; i++) {
            Tm1(i,i) = Tm1(i,i) - 1.0;
          }
          if (Tm1.normOne() > S_max) {
            if (infoTS_) {
              std::cout << "  large nonorthogonality in W(R)'*R detected \n";
            }
            return;
          }
        }
      }

      // Full orthogonalization.
      p.push_back(Wr->clone());
      (p[iterCG_-1])->set(*Wr);
      (p[iterCG_-1])->scale(-one);
      for (int j=1; j<iterCG_; j++) {
        Real scal = (p[iterCG_-1])->dot(*(Hp[j-1])) / (p[j-1])->dot(*(Hp[j-1]));
        Teuchos::RCP<Vector<Real> > pj = (p[j-1])->clone();
        pj->set(*p[j-1]);
        pj->scale(-scal);
        (p[iterCG_-1])->plus(*pj);
      }

      Hp.push_back(x.clone());
      obj.hessVec(*(Hp[iterCG_-1]), *(p[iterCG_-1]), x, zerotol);
      con.applyAdjointHessian(*vtemp, l, *(p[iterCG_-1]), x, zerotol);
      (Hp[iterCG_-1])->plus(*vtemp);
      pHp = (p[iterCG_-1])->dot(*(Hp[iterCG_-1]));
      rp  = (p[iterCG_-1])->dot(*r);

      normp = (p[iterCG_-1])->norm();
      normr = r->norm();

      // Negative curvature stopping condition.
      if (pHp <= 0) {
        pdesc->set(*(p[iterCG_-1])); // p is the descent direction
        if ((std::abs(rp) >= rptol*normp*normr) && (sgn(rp) == 1)) {
          pdesc->scale(-one); // -p is the descent direction
        }
	flagCG_ = 1;
        Real a = pdesc->dot(*pdesc);
        Real b = pdesc->dot(t);
        Real c = t.dot(t) - delta*delta;
        // Positive root of a*theta^2 + 2*b*theta + c = 0.
        Real theta = (-b + std::sqrt(b*b - a*c)) / a;
        vtemp->set(*(p[iterCG_-1]));
        vtemp->scale(theta);
        t.plus(*vtemp);
        // Store as tangential Cauchy point if terminating in first iteration.
        if (iterCG_ == 1) {
          tCP.set(t);
        }
	if (infoTS_) {
           std::cout << "  negative curvature detected \n";
        }
	return;
      }

      // Want to enforce nonzero alpha's.
      if (std::abs(rp) < rptol*normp*normr) {
        flagCG_ = 4;
        if (infoTS_) {
          std::cout << "  Zero alpha due to inexactness. \n";
        }
	return;
      }

      alpha = - rp/pHp;

      // Iterate update.
      tprev->set(t);
      vtemp->set(*(p[iterCG_-1]));
      vtemp->scale(alpha);
      t.plus(*vtemp);

      // Trust-region stopping condition.
      normt = t.norm();
      if (normt >= delta) {
        pdesc->set(*(p[iterCG_-1])); // p is the descent direction
        if (sgn(rp) == 1) {
          pdesc->scale(-one); // -p is the descent direction
        }
	flagCG_ = 1;
        Real a = pdesc->dot(*pdesc);
        Real b = pdesc->dot(t);
        Real c = t.dot(t) - delta*delta;
        // Positive root of a*theta^2 + 2*b*theta + c = 0.
        Real theta = (-b + std::sqrt(b*b - a*c)) / a;
        vtemp->set(*(p[iterCG_-1]));
        vtemp->scale(theta);
        t.set(*tprev);
        t.plus(*vtemp);
        flagCG_ = 2;
        // Store as tangential Cauchy point if terminating in first iteration.
        if (iterCG_ == 1) {
          tCP.set(t);
        }
	if (infoTS_) {
           std::cout << "  trust-region condition active \n";
        }
	return;
      }

      // Residual update.
      vtemp->set(*(Hp[iterCG_-1]));
      vtemp->scale(alpha);
      r->plus(*vtemp);
      if (orthocheck) {
        rs.push_back(r->clone());
        (rs[iterCG_])->set(*r);
      }

      iterCG_++;

    } // while (iterCG_ < maxiterCG_)

    flagCG_ = 3;
    if (infoTS_) {
      std::cout << "  maximum number of iterations reached \n";
    }

  } // solveTangentialSubproblem

  
  /** \brief Check acceptance of subproblem solutions, adjust trust-region radius and parameter, ensure global convergence.
  */
  void accept(Vector<Real> &s, Vector<Real> &n, Vector<Real> &t, Real f_new, Vector<Real> &c_new,
              Vector<Real> &gf_new, Vector<Real> &l_new, Vector<Real> &g_new,
              const Vector<Real> &x, const Vector<Real> &l, Real f, const Vector<Real> &gf, const Vector<Real> &c,
              const Vector<Real> &g, Vector<Real> &tCP, Vector<Real> &Wg,
              Objective<Real> &obj, EqualityConstraint<Real> &con, AlgorithmState<Real> &algo_state) {

    Real beta         = 1e-8;              // predicted reduction parameter
    Real tol_red_tang = 1e-3;              // internal reduction factor for tangtol
    Real tol_red_all  = 1e-1;              // internal reduction factor for qntol, lmhtol, pgtol, projtol, tangtol
    //bool glob_refine  = true;              // true  - if subsolver tolerances are adjusted in this routine, keep adjusted values globally
                                           // false - if subsolver tolerances are adjusted in this routine, discard adjusted values
    int ct_max        = 10;                // maximum number of globalization tries
    int mintol        = 1e-16;             // smallest tolerance value

    // Determines max value of |rpred|/pred. 
    Real rpred_over_pred = 0.5*(1-eta_);

    if (infoAC_) {
      std::stringstream hist;
      hist << "\n  SQP_accept\n";
      std::cout << hist.str();
    }

    Real zero      =  0.0;
    Real one       =  1.0;
    Real two       =  2.0;
    Real half      =  one/two;
    Real zerotol   =  zero;

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

    Teuchos::RCP<Vector<Real> > xtrial = x.clone();
    Teuchos::RCP<Vector<Real> > Jl = x.clone();
    Teuchos::RCP<Vector<Real> > gfJl = x.clone();
    Teuchos::RCP<Vector<Real> > Jnc = c.clone();
    Teuchos::RCP<Vector<Real> > t_orig = t.clone();
    Teuchos::RCP<Vector<Real> > Jt_orig = c.clone();
    Teuchos::RCP<Vector<Real> > t_m_tCP = t.clone();
    Teuchos::RCP<Vector<Real> > ltemp = l.clone();
    Teuchos::RCP<Vector<Real> > xtemp = x.clone();
    Teuchos::RCP<Vector<Real> > rt = c.clone();
    Teuchos::RCP<Vector<Real> > Hn = x.clone();
    Teuchos::RCP<Vector<Real> > Hto = x.clone();
    Teuchos::RCP<Vector<Real> > cxxvec = x.clone();
    Teuchos::RCP<Vector<Real> > czero = c.clone();
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
        Real tol = tangtol;
        con.solveAugmentedSystem(t, *ltemp, *t_orig, *czero, x, tol);         
        totalProj_++;
	con.applyJacobian(*rt, t, x, zerotol);
        linc_postproj = rt->norm();

        // Compute composite step.
        s.set(t);
        s.plus(n);

        // Compute some quantities before updating the objective and the constraint.
        obj.hessVec(*Hn, n, x, zerotol);
        con.applyAdjointHessian(*cxxvec, l, n, x, zerotol);
        Hn->plus(*cxxvec);
        obj.hessVec(*Hto, *t_orig, x, zerotol);
        con.applyAdjointHessian(*cxxvec, l, *t_orig, x, zerotol);
        Hto->plus(*cxxvec);

        // Compute objective, constraint, etc. values at the trial point.
        xtrial->set(x);
        xtrial->plus(s);
        obj.update(*xtrial,false,algo_state.iter);
        con.update(*xtrial,false,algo_state.iter);
        f_new = obj.value(*xtrial, zerotol);
std::cout << xtrial->dot(*(xtrial->basis(0))) << "\n";
std::cout << xtrial->dot(*(xtrial->basis(1))) << "\n";
std::cout << xtrial->dot(*(xtrial->basis(2))) << "\n";
std::cout << xtrial->dot(*(xtrial->basis(3))) << "\n";
std::cout << xtrial->dot(*(xtrial->basis(4))) << "\n";
std::cout << xtrial->norm() << "  " << f_new << "\n";
std::cout << t.dot(*(t.basis(0))) << "\n";
std::cout << t.dot(*(t.basis(1))) << "\n";
std::cout << t.dot(*(t.basis(2))) << "\n";
std::cout << t.dot(*(t.basis(3))) << "\n";
std::cout << t.dot(*(t.basis(4))) << "\n";
        obj.gradient(gf_new, *xtrial, zerotol);
        con.value(c_new, *xtrial, zerotol);
        l_new.set(l);
        computeLagrangeMultiplier(l_new, *xtrial, gf_new, con);

        // Penalty parameter update.
        part_pred = - Wg.dot(*t_orig);
        gfJl->set(gf);
        gfJl->plus(*Jl);
        part_pred -= gfJl->dot(n);
        part_pred -= half*Hn->dot(n);
        part_pred -= half*Hto->dot(*t_orig);
        ltemp->set(l);
        ltemp->axpy(-one, l_new);
        part_pred -= Jnc->dot(*ltemp);

        if ( part_pred < -half*penalty_*(c_normsquared-Jnc_normsquared) ) {
          penalty_ = ( -two * part_pred / (c_normsquared-Jnc_normsquared) ) + beta;
        }

        pred = part_pred + penalty_*(c_normsquared-Jnc_normsquared);

        // Computation of rpred.
        rpred = - ltemp->dot(*rt) - penalty_ * rt->dot(*rt) - two * penalty_ * rt->dot(*Jnc);
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
          lmhtol_  *= tol_red_all;
          qntol_   *= tol_red_all;
          pgtol_   *= tol_red_all;
          projtol_ *= tol_red_all;
          tangtol_ *= tol_red_all;
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

    // Compute actual reduction.
    ared = (f - f_new)  + (l.dot(c) - l_new.dot(c_new)) + penalty_*(c.dot(c) - c_new.dot(c_new));

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
        //hist <<   "         ||grad(L)_trial||   = %12.5e\n', normgl);
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

}; // class CompositeStepSQP

} // namespace ROL

#endif
