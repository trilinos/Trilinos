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

#ifndef ROL_PRIMALDUALACTIVESETSTEP_H
#define ROL_PRIMALDUALACTIVESETSTEP_H

#include "ROL_Step.hpp"
#include "ROL_Vector.hpp"
#include "ROL_Objective.hpp"
#include "ROL_BoundConstraint.hpp"
#include "ROL_Types.hpp"
#include "ROL_Secant.hpp"
#include "Teuchos_ParameterList.hpp"

/** @ingroup step_group
    \class ROL::PrimalDualActiveSetStep
    \brief Implements the computation of optimization steps 
           with the Newton primal-dual active set method.

    To describe primal-dual active set (PDAS), we consider the following 
    abstract setting.  Suppose \f$\mathcal{X}\f$ is a Hilbert space of 
    functions mapping \f$\Xi\f$ to \f$\mathbb{R}\f$.  For example, 
    \f$\Xi\subset\mathbb{R}^n\f$ and \f$\mathcal{X}=L^2(\Xi)\f$ or 
    \f$\Xi = \{1,\ldots,n\}\f$ and \f$\mathcal{X}=\mathbb{R}^n\f$. We 
    assume \f$ f:\mathcal{X}\to\mathbb{R}\f$ is twice-continuously Fr&eacute;chet 
    differentiable and \f$a,\,b\in\mathcal{X}\f$ with \f$a\le b\f$ almost 
    everywhere in \f$\Xi\f$.  Note that the PDAS algorithm will also work 
    with secant approximations of the Hessian. 

    Traditionally, PDAS is an algorithm for the minimizing quadratic objective 
    functions subject to bound constraints.  ROL implements a Newton PDAS which 
    extends PDAS to general bound-constrained nonlinear programs, i.e., 
    \f[
        \min_x \quad f(x) \quad \text{s.t.} \quad a \le x \le b.
    \f] 
    Given the \f$k\f$-th iterate \f$x_k\f$, the Newton PDAS algorithm computes 
    steps by applying PDAS to the quadratic subproblem 
    \f[
        \min_s \quad \langle \nabla^2 f(x_k)s + \nabla f(x_k),s \rangle_{\mathcal{X}}
        \quad \text{s.t.} \quad a \le x_k + s \le b.
    \f]
    For the \f$k\f$-th quadratic subproblem, PDAS builds an approximation of the 
    active set \f$\mathcal{A}_k\f$ using the dual variable \f$\lambda_k\f$ as 
    \f[
       \mathcal{A}^+_k = \{\,\xi\in\Xi\,:\,(\lambda_k + c(x_k-b))(\xi) > 0\,\}, \quad
       \mathcal{A}^-_k = \{\,\xi\in\Xi\,:\,(\lambda_k + c(x_k-a))(\xi) < 0\,\}, \quad\text{and}\quad
       \mathcal{A}_k = \mathcal{A}^-_k\cup\mathcal{A}^+_k.
    \f] 
    We define the inactive set \f$\mathcal{I}_k=\Xi\setminus\mathcal{A}_k\f$.
    The solution to the quadratic subproblem is then computed iteratively by solving 
    \f[
       \nabla^2 f(x_k) s_k + \lambda_{k+1} = -\nabla f(x_k), \quad
       x_k+s_k = a \;\text{on}\;\mathcal{A}^-_k,\quad x_k+s_k = b\;\text{on}\;\mathcal{A}^+_k,
       \quad\text{and}\quad
       \lambda_{k+1} = 0\;\text{on}\;\mathcal{I}_k
    \f]
    and updating the active and inactive sets. 
 
    One can rewrite this system by consolidating active and inactive parts, i.e., 
    \f[
       \begin{pmatrix}
           \nabla^2 f(x_k)_{\mathcal{A}_k,\mathcal{A}_k}  & \nabla^2 f(x_k)_{\mathcal{A}_k,\mathcal{I}_k} \\
           \nabla^2 f(x_k)_{\mathcal{I}_k,\mathcal{A}_k}  & \nabla^2 f(x_k)_{\mathcal{I}_k,\mathcal{I}_k} 
       \end{pmatrix}
       \begin{pmatrix}
         (s_k)_{\mathcal{A}_k} \\
         (s_k)_{\mathcal{I}_k}
       \end{pmatrix}
       +
       \begin{pmatrix}
         (\lambda_{k+1})_{\mathcal{A}_k} \\
         0
       \end{pmatrix}
       = - 
       \begin{pmatrix}
         \nabla f(x_k)_{\mathcal{A}_k}\\
         \nabla f(x_k)_{\mathcal{I}_k}
       \end{pmatrix}.
    \f]
    Here the subscripts \f$\mathcal{A}_k\f$ and \f$\mathcal{I}_k\f$ denote the active and inactive 
    components, respectively.  Moreover, the active components of \f$s_k\f$ are 
    \f$s_k(\xi) = a(\xi)-x_k(\xi)\f$ if \f$\xi\in\mathcal{A}^-_k\f$ and \f$s_k(\xi) = b(\xi)-x_k(\xi)\f$
    if \f$\xi\in\mathcal{A}^+_k\f$.  Since \f$(s_k)_{\mathcal{A}_k}\f$ is fixed, we only need to solve 
    for the inactive components of \f$s_k\f$ which we can do this using conjugate residuals (CR) (i.e., the 
    Hessian operator corresponding to the inactive indices may not be positive definite).  Once 
    \f$(s_k)_{\mathcal{I}_k}\f$ is computed, it is straight forward to update the dual variables.
*/

namespace ROL {

template <class Real>
class PrimalDualActiveSetStep : public Step<Real> {
private:
  // Krylov Parameters
  int maxitCR_; ///< Maximum number of CR iterations to solve quadratic subproblem
  int iterCR_;  ///< CR iteration counter
  int flagCR_;  ///< CR termination flag
  Real tol1_;   ///< CR absolute tolerance
  Real tol2_;   ///< CR relative tolerance
  Real itol_;   ///< Inexact CR tolerance

  // PDAS Parameters
  int maxit_;      ///< Maximum number of PDAS iterations 
  int iter_;       ///< PDAS iteration counter
  int flag_;       ///< PDAS termination flag
  Real stol_;      ///< PDAS minimum step size stopping tolerance
  Real gtol_;      ///< PDAS gradient stopping tolerance
  Real scale_;     ///< Scale for dual variables in the active set, \f$c\f$
  Real neps_;      ///< \f$\epsilon\f$-active set parameter 
  bool feasible_;  ///< Flag whether the current iterate is feasible or not

  // Dual Variable
  Teuchos::RCP<Vector<Real> > lambda_; ///< Container for dual variables

  // Secant Information
  ESecant esec_;                       ///< Enum for secant type
  Teuchos::RCP<Secant<Real> > secant_; ///< Secant object

  /** \brief Compute the gradient-based criticality measure.

             The criticality measure is 
             \f$\|x_k - P_{[a,b]}(x_k-\nabla f(x_k))\|_{\mathcal{X}}\f$.
             Here, \f$P_{[a,b]}\f$ denotes the projection onto the
             bound constraints.
 
             @param[in]    x     is the current iteration
             @param[in]    obj   is the objective function
             @param[in]    con   are the bound constraints
             @param[in]    tol   is a tolerance for inexact evaluations of the objective function
  */ 
  Real computeCriticalityMeasure(Vector<Real> &x, Objective<Real> &obj, BoundConstraint<Real> &con, Real tol) {
    Teuchos::RCP<StepState<Real> > step_state = Step<Real>::getState();
    obj.gradient(*(step_state->gradientVec),x,tol);
    Teuchos::RCP<Vector<Real> > xnew = x.clone();
    xnew->set(x);
    xnew->axpy(-1.0,*(step_state->gradientVec));
    con.project(*xnew);
    xnew->axpy(-1.0,x);
    return xnew->norm();
  }

  /** \brief Apply the inactive components of the Hessian operator.
 
             I.e., the components corresponding to \f$\mathcal{I}_k\f$.

             @param[out]       hv     is the result of applying the Hessian at @b x to 
                                      @b v
             @param[in]        v      is the direction in which we apply the Hessian
             @param[in]        x      is the current iteration vector \f$x_k\f$
             @param[in]        xlam   is the vector \f$x_k + c\lambda_k\f$
             @param[in]        obj    is the objective function
             @param[in]        con    are the bound constraints
  */
  void applyInactiveHessian(Vector<Real> &hv, Vector<Real> &v, const Vector<Real> &x, 
                      const Vector<Real> &xlam, Objective<Real> &obj, BoundConstraint<Real> &con) {
    con.pruneActive(v,xlam,this->neps_);
    if ( this->secant_ == Teuchos::null ) {
      obj.hessVec(hv,v,x,this->itol_);
    }
    else {
      this->secant_->applyB(hv,v,x);
    }
    con.pruneActive(hv,xlam,this->neps_);
  }

  /** \brief Apply the inactive components of the preconditioner operator.

             I.e., the components corresponding to \f$\mathcal{I}_k\f$.

             @param[out]       hv     is the result of applying the preconditioner at @b x to 
                                      @b v
             @param[in]        v      is the direction in which we apply the preconditioner
             @param[in]        x      is the current iteration vector \f$x_k\f$
             @param[in]        xlam   is the vector \f$x_k + c\lambda_k\f$
             @param[in]        obj    is the objective function
             @param[in]        con    are the bound constraints
  */
  void applyInactivePrecond(Vector<Real> &pv, Vector<Real> &v, const Vector<Real> &x,
                      const Vector<Real> &xlam, Objective<Real> &obj, BoundConstraint<Real> &con) {
    con.pruneActive(v,xlam,this->neps_);
    obj.precond(pv,v,x,this->itol_);
    con.pruneActive(pv,xlam,this->neps_);
  }

  /** \brief Solve the inactive part of the PDAS optimality system.  

             The inactive PDAS optimality system is 
             \f[
                 \nabla^2 f(x_k)_{\mathcal{I}_k,\mathcal{I}_k}s  = 
                     -\nabla f(x_k)_{\mathcal{I}_k}
                     -\nabla^2 f(x_k)_{\mathcal{I}_k,\mathcal{A}_k} (s_k)_{\mathcal{A}_k}.
             \f]
             Since the inactive part of the Hessian may not be positive definite, we solve 
             using CR.
   
             @param[out]       sol    is the vector containing the solution
             @param[in]        rhs    is the right-hand side vector
             @param[in]        xlam   is the vector \f$x_k + c\lambda_k\f$
             @param[in]        x      is the current iteration vector \f$x_k\f$
             @param[in]        obj    is the objective function
             @param[in]        con    are the bound constraints
  */
  // Solve the inactive part of the optimality system using conjugate residuals
  void solve(Vector<Real> &sol, const Vector<Real> &rhs, const Vector<Real> &xlam, const Vector<Real> &x, 
             Objective<Real> &obj, BoundConstraint<Real> &con) {
    // Initialize Residual
    Teuchos::RCP<Vector<Real> > res = rhs.clone();
    res->set(rhs);
    Real rnorm  = res->norm(); 
    Real rtol   = std::min(this->tol1_,this->tol2_*rnorm);
    if ( false ) { this->itol_ = rtol/(this->maxitCR_*rnorm); }
    sol.zero();

    // Apply preconditioner to residual r = Mres
    Teuchos::RCP<Vector<Real> > r = x.clone();
    this->applyInactivePrecond(*r,*res,x,xlam,obj,con);

    // Initialize direction p = v
    Teuchos::RCP<Vector<Real> > p = x.clone();
    p->set(*r);

    // Apply Hessian to v
    Teuchos::RCP<Vector<Real> > Hr = x.clone();
    this->applyInactiveHessian(*Hr,*r,x,xlam,obj,con);

    // Apply Hessian to p
    Teuchos::RCP<Vector<Real> > Hp  = x.clone();
    Teuchos::RCP<Vector<Real> > MHp = x.clone();
    Hp->set(*Hr);

    this->iterCR_ = 0;
    this->flagCR_ = 0;

    Real kappa = 0.0, beta  = 0.0, alpha = 0.0, tmp = 0.0, rHr = Hr->dot(*r);

    for (this->iterCR_ = 0; this->iterCR_ < this->maxitCR_; this->iterCR_++) {
      // Precondition Hp
      this->applyInactivePrecond(*MHp,*Hp,x,xlam,obj,con);

      kappa = Hp->dot(*MHp);  // p' H M H p
      alpha = rHr/kappa;      // r' M H M r
      sol.axpy(alpha,*p);     // update step
      res->axpy(-alpha,*Hp);  // residual
      r->axpy(-alpha,*MHp);   // preconditioned residual
      
      // recompute rnorm and decide whether or not to exit
      rnorm = res->norm();
      if ( rnorm < rtol ) { break; }

      // Apply Hessian to v
      this->itol_ = rtol/(this->maxitCR_*rnorm);
      this->applyInactiveHessian(*Hr,*r,x,xlam,obj,con);

      tmp  = rHr;
      rHr  = Hr->dot(*r);
      beta = rHr/tmp;
      p->scale(beta);
      p->axpy(1.0,*r);
      Hp->scale(beta);
      Hp->axpy(1.0,*Hr);
    }
    if ( this->iterCR_ == this->maxitCR_ ) {
      this->flagCR_ = 1;
    }
    else {
      this->iterCR_++;
    }
  }

public:
  /** \brief Constructor.
     
             @param[in]     parlist   is a parameter list containing relevent algorithmic information
             @param[in]     useSecant is a bool which determines whether or not the algorithm uses 
                                      a secant approximation of the Hessian
  */
  PrimalDualActiveSetStep( Teuchos::ParameterList &parlist, bool useSecant = false ) 
    : Step<Real>::Step(), iterCR_(0), flagCR_(0), iter_(0), flag_(0), neps_(-ROL_EPSILON), feasible_(false) {
    maxitCR_ = parlist.get("Maximum Number of Krylov Iterations", 50);
    tol1_    = parlist.get("Absolute Krylov Tolerance", 1.e-4);
    tol2_    = parlist.get("Relative Krylov Tolerance", 1.e-2);
    esec_    = StringToESecant(parlist.get("Secant Type","Limited-Memory BFGS"));

    maxit_   = parlist.get("PDAS Maximum Number of Iterations",10);
    stol_    = parlist.get("PDAS Relative Step Tolerance",1.e-8);
    gtol_    = parlist.get("PDAS Relative Gradient Tolerance",1.e-6);
    scale_   = parlist.get("PDAS Dual Scaling", 1.0);
  
    secant_  = Teuchos::null;
    if ( useSecant ) {
      int L   = parlist.get("Maximum Secant Storage",10);
      int BB  = parlist.get("Barzilai-Borwein",1);
      secant_ = getSecant<Real>(esec_,L,BB); 
    }
  }

  /** \brief Initialize step.  

             This includes projecting the initial guess onto the constraints, 
             computing the initial objective function value and gradient, 
             and initializing the dual variables.

             @param[in,out]    x           is the initial guess 
             @param[in]        obj         is the objective function
             @param[in]        con         are the bound constraints
             @param[in]        algo_state  is the current state of the algorithm
  */
  void initialize( Vector<Real> &x, Objective<Real> &obj, BoundConstraint<Real> &con, 
                   AlgorithmState<Real> &algo_state ) {
    Teuchos::RCP<StepState<Real> > step_state = Step<Real>::getState();
    // Initialize state descent direction and gradient storage
    step_state->descentVec  = x.clone();
    step_state->gradientVec = x.clone();
    step_state->searchSize  = 0.0;
    // Project x onto constraint set
    con.project(x);
    // Update objective function, get value, and get gradient
    Real tol = std::sqrt(ROL_EPSILON);
    obj.update(x,true,algo_state.iter);
    algo_state.value = obj.value(x,tol);
    algo_state.nfval++;
    algo_state.gnorm = this->computeCriticalityMeasure(x,obj,con,tol);
    algo_state.ngrad++;
    // Initialize dual variable
    this->lambda_ = x.clone(); this->lambda_->zero();
    this->lambda_->set(*(step_state->gradientVec));
    this->lambda_->scale(-1.0);
    //con.setVectorToLowerBound(*this->lambda_);
  }

  /** \brief Compute step.

             Given \f$x_k\f$, this function first builds the 
             primal-dual active sets
             \f$\mathcal{A}_k^-\f$ and \f$\mathcal{A}_k^+\f$.  
             Next, it uses CR to compute the inactive 
             components of the step by solving 
             \f[
                 \nabla^2 f(x_k)_{\mathcal{I}_k,\mathcal{I}_k}(s_k)_{\mathcal{I}_k}  = 
                     -\nabla f(x_k)_{\mathcal{I}_k}
                     -\nabla^2 f(x_k)_{\mathcal{I}_k,\mathcal{A}_k} (s_k)_{\mathcal{A}_k}.
             \f]
             Finally, it updates the active components of the 
             dual variables as 
             \f[
                \lambda_{k+1} = -\nabla f(x_k)_{\mathcal{A}_k} 
                                -(\nabla^2 f(x_k) s_k)_{\mathcal{A}_k}.
             \f]

             @param[out]       s           is the step computed via PDAS
             @param[in]        x           is the current iterate
             @param[in]        obj         is the objective function
             @param[in]        con         are the bound constraints
             @param[in]        algo_state  is the current state of the algorithm
  */
  void compute( Vector<Real> &s, const Vector<Real> &x, Objective<Real> &obj, BoundConstraint<Real> &con, 
                AlgorithmState<Real> &algo_state ) {
    // Initialize storage
    Teuchos::RCP<Vector<Real> > x0   = x.clone();
    Teuchos::RCP<Vector<Real> > xlam = x.clone();
    Teuchos::RCP<Vector<Real> > xbnd = x.clone();
    Teuchos::RCP<Vector<Real> > As   = x.clone();
    Teuchos::RCP<Vector<Real> > Ag   = x.clone();
    Teuchos::RCP<Vector<Real> > rhs  = x.clone();
    Teuchos::RCP<Vector<Real> > Hs   = x.clone();
    Teuchos::RCP<Vector<Real> > IHs  = x.clone();
    Teuchos::RCP<Vector<Real> > IHAs = x.clone();
    Teuchos::RCP<Vector<Real> > tmp  = x.clone();
    Teuchos::RCP<Vector<Real> > res  = x.clone();

    Teuchos::RCP<StepState<Real> > step_state = Step<Real>::getState();
    s.zero();
    x0->set(x);
    res->set(*(step_state->gradientVec));
    for ( this->iter_ = 0; this->iter_ < this->maxit_; this->iter_++ ) {
      /********************************************************************/
      // MODIFY ITERATE VECTOR TO CHECK ACTIVE SET
      /********************************************************************/
      xlam->set(*x0);                               // xlam = x0
      xlam->axpy(this->scale_,*(this->lambda_));    // xlam = x0 + c*lambda
      /********************************************************************/
      // PROJECT x ONTO PRIMAL DUAL FEASIBLE SET
      /********************************************************************/
      As->zero();                                   // As   = 0
  
      con.setVectorToUpperBound(*xbnd);             // xbnd = u        
      xbnd->axpy(-1.0,x);                           // xbnd = u - x    
      tmp->set(*xbnd);                              // tmp  = u - x    
      con.pruneUpperActive(*tmp,*xlam,this->neps_); // tmp  = I(u - x) 
      xbnd->axpy(-1.0,*tmp);                        // xbnd = A(u - x)  
      As->plus(*xbnd);                              // As  += A(u - x)

      con.setVectorToLowerBound(*xbnd);             // xbnd = l
      xbnd->axpy(-1.0,x);                           // xbnd = l - x
      tmp->set(*xbnd);                              // tmp  = l - x
      con.pruneLowerActive(*tmp,*xlam,this->neps_); // tmp  = I(l - x)
      xbnd->axpy(-1.0,*tmp);                        // xbnd = A(l - x)
      As->plus(*xbnd);                              // As  += A(l - x)
      /********************************************************************/
      // APPLY HESSIAN TO ACTIVE COMPONENTS OF s AND REMOVE INACTIVE 
      /********************************************************************/
      this->itol_ = std::sqrt(ROL_EPSILON);
      if ( this->secant_ == Teuchos::null ) {       // IHAs = H*As
        obj.hessVec(*IHAs,*As,x,this->itol_);
      }
      else {
        this->secant_->applyB(*IHAs,*As,x);
      }
      con.pruneActive(*IHAs,*xlam,this->neps_);     // IHAs = I(H*As)
      /********************************************************************/
      // SEPARATE ACTIVE AND INACTIVE COMPONENTS OF THE GRADIENT
      /********************************************************************/
      rhs->set(*(step_state->gradientVec));    // Inactive components
      con.pruneActive(*rhs,*xlam,this->neps_);

      Ag->set(*(step_state->gradientVec));     // Active components
      Ag->axpy(-1.0,*rhs);
      /********************************************************************/
      // SOLVE REDUCED NEWTON SYSTEM 
      /********************************************************************/
      rhs->plus(*IHAs);
      rhs->scale(-1.0);                        // rhs = -Ig - I(H*As)
      s.zero();
      if ( rhs->norm() > 0.0 ) {             
        this->solve(s,*rhs,*xlam,x,obj,con);   // Call conjugate residuals
        con.pruneActive(s,*xlam,this->neps_);  // s <- Is
      }
      s.plus(*As);                             // s = Is + As
      /********************************************************************/
      // UPDATE MULTIPLIER 
      /********************************************************************/
      if ( this->secant_ == Teuchos::null ) {
        obj.hessVec(*Hs,s,x,this->itol_);
      }
      else {
        this->secant_->applyB(*Hs,s,x);
      }
      IHs->set(*Hs);
      con.pruneActive(*IHs,*xlam,this->neps_);
      this->lambda_->set(*Hs);
      this->lambda_->axpy(-1.0,*IHs);
      this->lambda_->plus(*Ag);
      this->lambda_->scale(-1.0);
      /********************************************************************/
      // UPDATE STEP 
      /********************************************************************/
      x0->set(x);
      x0->plus(s);
      res->set(*(step_state->gradientVec));
      res->plus(*Hs);
      // Compute criticality measure  
      tmp->set(*x0);
      tmp->axpy(-1.0,*res);
      con.project(*tmp);
      tmp->axpy(-1.0,*x0);
//      std::cout << s.norm()              << "  " 
//                << tmp->norm()           << "  " 
//                << res->norm()           << "  " 
//                << this->lambda_->norm() << "  " 
//                << this->flagCR_         << "  " 
//                << this->iterCR_         << "\n";
      if ( tmp->norm() < this->gtol_*algo_state.gnorm ) {
        this->flag_ = 0;
        break;
      }
      if ( s.norm() < this->stol_*x.norm() ) {
        this->flag_ = 2;
        break;
      } 
    }
    if ( this->iter_ == this->maxit_ ) {
      this->flag_ = 1;
    }
    else {
      this->iter_++;
    }
  }

  /** \brief Update step, if successful.

             This function returns \f$x_{k+1} = x_k + s_k\f$.
             It also updates secant information if being used.

             @param[in]        x           is the new iterate
             @param[out]       s           is the step computed via PDAS
             @param[in]        obj         is the objective function
             @param[in]        con         are the bound constraints
             @param[in]        algo_state  is the current state of the algorithm
  */
  void update( Vector<Real> &x, const Vector<Real> &s, Objective<Real> &obj, BoundConstraint<Real> &con,
               AlgorithmState<Real> &algo_state ) {
    Teuchos::RCP<StepState<Real> > step_state = Step<Real>::getState();

    x.plus(s);
    this->feasible_ = con.isFeasible(x);
    algo_state.snorm = s.norm();
    algo_state.iter++;
    Real tol = std::sqrt(ROL_EPSILON);
    obj.update(x,true,algo_state.iter);
    algo_state.value = obj.value(x,tol);
    algo_state.nfval++;
    
    Teuchos::RCP<Vector<Real> > gp;
    if ( this->secant_ != Teuchos::null ) {
      gp = x.clone();
      gp->set(*(step_state->gradientVec));
    }
    algo_state.gnorm = this->computeCriticalityMeasure(x,obj,con,tol);
    algo_state.ngrad++;

    if ( this->secant_ != Teuchos::null ) {
      this->secant_->update(*(step_state->gradientVec),*gp,s,algo_state.snorm,algo_state.iter+1);
    }
  }

  /** \brief Print iterate header.

             This function produces a string containing 
             header information.  
  */
  std::string printHeader( void ) const {
    std::stringstream hist;
    hist << "  ";
    hist << std::setw(6) << std::left << "iter";
    hist << std::setw(15) << std::left << "value";
    hist << std::setw(15) << std::left << "gnorm";
    hist << std::setw(15) << std::left << "snorm";
    hist << std::setw(10) << std::left << "#fval";
    hist << std::setw(10) << std::left << "#grad";
    if ( this->maxit_ > 1 ) {
      hist << std::setw(10) << std::left << "iterPDAS";
      hist << std::setw(10) << std::left << "flagPDAS";
    }
    else {
      hist << std::setw(10) << std::left << "iterCR";
      hist << std::setw(10) << std::left << "flagCR";
    }
    hist << std::setw(10) << std::left << "feasible";
    hist << "\n";
    return hist.str();
  }

  /** \brief Print step name.

             This function produces a string containing 
             the algorithmic step information.  
  */
  std::string printName( void ) const {
    std::stringstream hist;
    hist << "\nPrimal Dual Active Set Newton's Method\n";
    return hist.str();
  }

  /** \brief Print iterate status.
    
             This function prints the iteration status.

             @param[in]        algo_state  is the current state of the algorithm
             @param[in]        printHeader if set to true will print the header at each iteration
  */
  virtual std::string print( AlgorithmState<Real> &algo_state, bool printHeader = false ) const {
    std::stringstream hist;
    hist << std::scientific << std::setprecision(6);
    if ( algo_state.iter == 0 ) {
      hist << this->printName();
    }
    if ( printHeader ) {
      hist << this->printHeader();
    }
    if ( algo_state.iter == 0 ) {
      hist << "  ";
      hist << std::setw(6) << std::left << algo_state.iter;
      hist << std::setw(15) << std::left << algo_state.value;
      hist << std::setw(15) << std::left << algo_state.gnorm;
      hist << "\n";
    }
    else {
      hist << "  ";
      hist << std::setw(6) << std::left << algo_state.iter;
      hist << std::setw(15) << std::left << algo_state.value;
      hist << std::setw(15) << std::left << algo_state.gnorm;
      hist << std::setw(15) << std::left << algo_state.snorm;
      hist << std::setw(10) << std::left << algo_state.nfval;
      hist << std::setw(10) << std::left << algo_state.ngrad;
      if ( this->maxit_ > 1 ) {
        hist << std::setw(10) << std::left << this->iter_;
        hist << std::setw(10) << std::left << this->flag_;
      }
      else {
        hist << std::setw(10) << std::left << this->iterCR_;
        hist << std::setw(10) << std::left << this->flagCR_;
      }
      if ( this->feasible_ ) {
        hist << std::setw(10) << std::left << "YES";
      }
      else {
        hist << std::setw(10) << std::left << "NO";
      }
      hist << "\n";
    }
    return hist.str();
  }

}; // class PrimalDualActiveSetStep

} // namespace ROL

#endif

//  void solve(Vector<Real> &sol, const Vector<Real> &rhs, const Vector<Real> &xlam, const Vector<Real> &x, 
//             Objective<Real> &obj, BoundConstraint<Real> &con) {
//    Real rnorm  = rhs.norm(); 
//    Real rtol   = std::min(this->tol1_,this->tol2_*rnorm);
//    this->itol_ = std::sqrt(ROL_EPSILON);
//    sol.zero();
//
//    Teuchos::RCP<Vector<Real> > res = rhs.clone();
//    res->set(rhs);
//
//    Teuchos::RCP<Vector<Real> > v = x.clone();
//    con.pruneActive(*res,xlam,this->neps_);
//    obj.precond(*v,*res,x,this->itol_);
//    con.pruneActive(*v,xlam,this->neps_);
//
//    Teuchos::RCP<Vector<Real> > p = x.clone();
//    p->set(*v);
//
//    Teuchos::RCP<Vector<Real> > Hp = x.clone();
//
//    this->iterCR_ = 0;
//    this->flagCR_ = 0;
//
//    Real kappa = 0.0, beta  = 0.0, alpha = 0.0, tmp = 0.0, rv = v->dot(*res);
//
//    for (this->iterCR_ = 0; this->iterCR_ < this->maxitCR_; this->iterCR_++) {
//      if ( false ) {
//        this->itol_ = rtol/(this->maxitCR_*rnorm);
//      }
//      con.pruneActive(*p,xlam,this->neps_);
//      if ( this->secant_ == Teuchos::null ) {
//        obj.hessVec(*Hp, *p, x, this->itol_);
//      }
//      else {
//        this->secant_->applyB( *Hp, *p, x );
//      }
//      con.pruneActive(*Hp,xlam,this->neps_);
//
//      kappa = p->dot(*Hp);
//      if ( kappa <= 0.0 ) { this->flagCR_ = 2; break; }
//      alpha = rv/kappa;
//      sol.axpy(alpha,*p);
//
//      res->axpy(-alpha,*Hp);
//      rnorm = res->norm();
//      if ( rnorm < rtol ) { break; }
//
//      con.pruneActive(*res,xlam,this->neps_);
//      obj.precond(*v,*res,x,this->itol_);
//      con.pruneActive(*v,xlam,this->neps_);
//      tmp  = rv;
//      rv   = v->dot(*res);
//      beta = rv/tmp;
//
//      p->scale(beta);
//      p->axpy(1.0,*v);
//    }
//    if ( this->iterCR_ == this->maxitCR_ ) {
//      this->flagCR_ = 1;
//    }
//    else {
//      this->iterCR_++;
//    }
//  }
