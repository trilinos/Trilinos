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

/** \class ROL::CompositeStepSQP
    \brief Implements the computation of optimization steps
           with composite-step trust-region SQP methods.
*/


namespace ROL {

template <class Real>
class CompositeStepSQP : public Step<Real> {
private:

  int TRflag_;
  int CGflag_;
  int CGiter_;

  int maxiterCG_;
  Real tolCG_;

  Real lmhtol_;
  Real qntol_;
  Real pgtol_;
  Real projtol_;
  Real tangtol_;
  Real tntmax_;

  Real zeta_;
  Real Delta_;

  bool debugQN_;
  bool debugTANGSUB_;

public:

  virtual ~CompositeStepSQP() {}

  CompositeStepSQP( Teuchos::ParameterList & parlist ) : Step<Real>() {
    //Teuchos::RCP<StepState<Real> > step_state = Step<Real>::getState();
    TRflag_ = 0;
    CGflag_ = 0;
    CGiter_ = 0;

    maxiterCG_ = 20;
    tolCG_ = 1e-4;

    Real nominal_tol = parlist.get("Nominal SQP Optimality Solver Tolerance", 1e-3);
    lmhtol_  = nominal_tol;
    qntol_   = nominal_tol;
    pgtol_   = nominal_tol;       
    projtol_ = nominal_tol;     
    tangtol_ = nominal_tol;     
    
    zeta_  = 0.9;
    Delta_ = 1e2;

    debugQN_ = true;
    debugTANGSUB_ = true;
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
    Teuchos::RCP<Vector<Real> > n = s.clone();
    Teuchos::RCP<Vector<Real> > c = l.clone();
    con.value(*c, x, zerotol);
    computeQuasinormalStep(*n, *c, x, zeta_*Delta_, con);
  }

  /** \brief Update step, if successful.
  */
  void update( Vector<Real> &x, Vector<Real> &l, const Vector<Real> &s,
               Objective<Real> &obj, EqualityConstraint<Real> &con, 
               AlgorithmState<Real> &algo_state ) {
    //Teuchos::RCP<StepState<Real> > state = Step<Real>::getState();

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
      hist << std::setw(10) << std::left << algo_state.nfval;              
      hist << std::setw(10) << std::left << algo_state.ngrad;              
      hist << std::setw(10) << std::left << TRflag_;              
      hist << std::setw(10) << std::left << CGiter_;
      hist << std::setw(10) << std::left << CGflag_;
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

    if (debugQN_) {
      std::cout << "\n  SQP_quasi-normal_step\n";
    }

    Real zero      = 0.0;
    Real zerotol   = 0.0;
    Real negative1 = -1.0;

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
      if (debugQN_) {
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
    ctemp->scale(negative1);
    nCPtemp->set(*nCP);
    nCPtemp->scale(negative1);
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
      if (debugQN_) {
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
      if (debugQN_) {
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
  void solveTangentialSubproblem(Vector<Real> &t, Vector<Real> &tCP, Vector<Real> Wg,
                                 const Vector<Real> &x, const Vector<Real> &g, const Vector<Real> &n, const Vector<Real> &l,
                                 Real delta, Objective<Real> &obj, EqualityConstraint<Real> &con) {

    /* Initialization of the CG step. */
    Real zero = 0.0;
    Real zerotol = 0.0;
    int iter = 1;
    int flag = 0;
    t.zero();
    tCP.zero();
    Teuchos::RCP<Vector<Real> > r = g.clone();
    Teuchos::RCP<Vector<Real> > Wr = g.clone();
    Teuchos::RCP<Vector<Real> > vtemp = g.clone();
    Teuchos::RCP<Vector<Real> > ltemp = l.clone();
    Teuchos::RCP<Vector<Real> > czero = l.clone();
    czero->zero();
    obj.hessVec(*r, n, x, zerotol);
    con.applyAdjointHessian(*vtemp, l, n, x, zerotol);
    r->plus(*vtemp);
    Real normg = r->norm();
    Real normWg = zero;
    Real normWr = zero;

    std::vector<Teuchos::RCP<Vector<Real > > >  p;   // stores search directions
    std::vector<Teuchos::RCP<Vector<Real > > >  Hp;  // stores hessvec's applied to p's
    std::vector<Teuchos::RCP<Vector<Real > > >  rs;  // stores residuals
    std::vector<Teuchos::RCP<Vector<Real > > >  Wrs; // stores projected residuals

    Real rptol = 1e-12;

    if (debugTANGSUB_) {
      std::cout << "\n  SQP_tangential_subproblem\n";
      std::cout << "  iter    ||W*r||/||W*r0||     ||D*s||        delta       ||Jac*s||\n";
    }

    if (normg == 0) {
      if (debugTANGSUB_) {
        std::cout << "    >>> Tangential subproblem: Initial gradient is zero! \n";
      }
      iter = 0; Wg.zero(); flag = 0;
      return;
    }

    /* Start CG loop. */
    while (iter < maxiterCG_)

      // Store tangential Cauchy point (which is the current iterate in the second iteration).
      if (iter == 2) {
        tCP.set(t);
      }

      // Compute (inexact) projection W*r.
      if (iter == 1) {
        // Solve augmented system.
        Real tol = pgtol_;
        con.solveAugmentedSystem(*Wr, *ltemp, *r, *czero, x, tol);
        Wg.set(*Wr);
        Real normWg = Wg.norm();
        // Check if done (small initial projected residual).
        if (normWg < 1e-14) {
          flag = 0;
          iter = iter-1;
          if (debugTANGSUB_) {
            std::cout << "  Initial projected residual is close to zero! \n";
          }
          return;
        }
        // Set first residual to projected gradient.
        r->set(Wg);
      }
      else {
        // Solve augmented system.
        Real tol = projtol_;
        con.solveAugmentedSystem(*Wr, *ltemp, *r, *czero, x, tol);
      }
      normWr = Wr->norm();

      if (debugTANGSUB_) {
         Teuchos::RCP<Vector<Real> > ct = l.clone();
         con.applyJacobian(*ct, t, x, zerotol);
         Real linc = ct->norm();
         printf("%5d     %12.5e     %12.5e  %12.5e  %12.5e \n", iter-1, normWr/normWg, t->norm(), delta, linc);
      }

      // Check if done (small relative residual).
      if (normWr/normWg < tolCG_) {
        flag = 0;
        iter = iter-1;
        if (debugTANGSUB_) {
          std::cout << "  || W(g + H*(n+s)) || <= cgtol*|| W(g + H*n)|| \n";
        }
        return;
      }


  } // solveTangentialSubproblem



}; // class CompositeStepSQP

} // namespace ROL

#endif
