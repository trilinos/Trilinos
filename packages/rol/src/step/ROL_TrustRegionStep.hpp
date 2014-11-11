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

#ifndef ROL_TRUSTREGIONSTEP_H
#define ROL_TRUSTREGIONSTEP_H

#include "ROL_Types.hpp"
#include "ROL_Step.hpp"
#include "ROL_Secant.hpp"
#include "ROL_TrustRegion.hpp"
#include <sstream>
#include <iomanip>

/** @ingroup step_group
    \class ROL::TrustRegionStep
    \brief Provides the interface to compute optimization steps
           with trust regions.

    Suppose \f$\mathcal{X}\f$ is a Hilbert space of 
    functions mapping \f$\Xi\f$ to \f$\mathbb{R}\f$.  For example, 
    \f$\Xi\subset\mathbb{R}^n\f$ and \f$\mathcal{X}=L^2(\Xi)\f$ or 
    \f$\Xi = \{1,\ldots,n\}\f$ and \f$\mathcal{X}=\mathbb{R}^n\f$. We 
    assume \f$f:\mathcal{X}\to\mathbb{R}\f$ is twice-continuously Fr&eacute;chet 
    differentiable and \f$a,\,b\in\mathcal{X}\f$ with \f$a\le b\f$ almost 
    everywhere in \f$\Xi\f$.  Note that these trust-region algorithms will also work 
    with secant approximations of the Hessian. 
    This step applies to unconstrained and bound constrained optimization problems,
    \f[
        \min_x\quad f(x) \qquad\text{and}\qquad \min_x\quad f(x)\quad\text{s.t.}\quad a\le x\le b,
    \f]
    respectively.  

    For unconstrained problems, given the \f$k\f$-th iterate \f$x_k\f$ the trial step
    \f$s_k\f$ is computed by approximately solving the trust-region subproblem 
    \f[
       \min_{s} \frac{1}{2}\langle B_k s, s\rangle_{\mathcal{X}} + \langle g_k,s\rangle_{\mathcal{X}}
           \quad\text{s.t.}\quad \|s\|_{\mathcal{X}} \le \Delta_k
    \f]
    where \f$B_k\in L(\mathcal{X},\mathcal{X})\f$, \f$g_k\approx\nabla f(x_k)\f$, and \f$\Delta_k > 0\f$.
    The approximate minimizer \f$s_k\f$ must satisfy the fraction of Cauchy decrease condition
    \f[
       -\frac{1}{2}\langle B_k s, s\rangle_{\mathcal{X}} - \langle g_k,s\rangle_{\mathcal{X}}
          \ge \kappa_0 \|g_k\|_{\mathcal{X}}
          \min\left\{\,\Delta_k,\,
          \frac{\|g_k\|_{\mathcal{X}}}{1+\|B_k\|_{L(\mathcal{X},\mathcal{X}})}\,\right\}
    \f]
    for some \f$\kappa_0>0\f$ independent of \f$k\f$.
    ROL's trust-region algorithm allows for both inexact objective function and gradient evaluation.  
    The user must ensure that the inexact objective function, \f$f_k\f$ satisfies
    \f[
       |(f(x_k+s_k)-f_k(x_k+s_k)) - (f(x_k)-f_k(x_k))| \le 
        \eta_1 \min\{\,-\frac{1}{2}\langle B_k s, s\rangle_{\mathcal{X}} - \langle g_k,s\rangle_{\mathcal{X}},
                       \,r_k\,\}
    \f]
    where \f$\eta_1\f$ is the step acceptance threshold and \f$r_k\f$ is a user-defined forcing sequence of 
    positive numbers converging to zero.  The inexact gradient, \f$g_k\f$, must satisfy
    \f[
       \|g_k-\nabla J(x_k)\|_{\mathcal{X}} \le \kappa_1\min\{\,\|g_k\|_{\mathcal{X}},\,\Delta_k\,\}
    \f]
    where \f$\kappa_1 > 0\f$ is independent of \f$k\f$.

    For bound constrained problems, ROL employs projected Newton-type methods.  
    For these methods, ROL requires the notion of an active set of an iterate \f$x_k\f$, 
    \f[
       \mathcal{A}_k = \{\, \xi\in\Xi\,:\,x_k(\xi) = a(\xi)\,\}\cap
                       \{\, \xi\in\Xi\,:\,x_k(\xi) = b(\xi)\,\}.
    \f]
    Given \f$\mathcal{A}_k\f$ and a gradient approximation \f$g_k\f$, we define the binding set as
    \f[
       \mathcal{B}_k = \{\, \xi\in\Xi\,:\,x_k(\xi) = a(\xi) \;\text{and}\; -g_k(\xi) < 0 \,\}\cap
                       \{\, \xi\in\Xi\,:\,x_k(\xi) = b(\xi) \;\text{and}\; -g_k(\xi) > 0 \,\}.
    \f]
    The binding set contains the values of \f$\xi\in\Xi\f$ such that if \f$x_k(\xi)\f$ is on a 
    bound, then \f$(x_k+s_k)(\xi)\f$ will violate bound.  Using these definitions, ROL 
    prunes the variables in the binding set and runs a standard trust-region subproblem solver on the 
    free variables.  ROL then must perform a projected search to ensure the fraction of Cauchy decrease 
    condition is satisfied.

    TrustRegionStep implements a number of algorithms for both bound constrained and unconstrained 
    optimization.  These algorithms are: Cauchy Point, Dogleg, Double Dogleg, and Truncated CG.  
    Each of these methods can be run using a secant approximation of the Hessian. 
    These methods are chosen through the ETrustRegion enum.
*/


namespace ROL {

template <class Real>
class TrustRegionStep : public Step<Real> {
private:

  Teuchos::RCP<Secant<Real> >       secant_;      ///< Container for secant approximation.
  Teuchos::RCP<TrustRegion<Real> >  trustRegion_; ///< Container for trust-region object.

  Teuchos::RCP<Vector<Real> >       xnew_;        ///< Container for updated iteration vector.
  Teuchos::RCP<Vector<Real> >       xold_;        ///< Container for previous iteration vector.
  Teuchos::RCP<Vector<Real> >       gp_;          ///< Container for previous gradient vector.

  ETrustRegion      etr_;        ///< Trust-region subproblem solver type.
  ESecant           esec_;       ///< Secant type.

  bool useSecantHessVec_;        ///< Flag whether to use a secant Hessian.
  bool useSecantPrecond_;        ///< Flag whether to use a secant preconditioner. 

  bool useProjectedGrad_;        ///< Flag whether to use the projected gradient criticality measure.

  std::vector<bool> useInexact_; ///< Contains flags for inexact (0) objective function, (1) gradient, (2) Hessian.
  int               TRflag_  ;   ///< Trust-region exit flag.
  int               TR_nfval_;   ///< Trust-region function evaluation counter.
  int               TR_ngrad_;   ///< Trust-region gradient evaluation counter.
  int               CGflag_;     ///< Truncated CG termination flag.
  int               CGiter_;     ///< Truncated CG iteration count.

  Real              alpha_init_; ///< Initial line-search parameter for projected methods.
  int               max_fval_;   ///< Maximum function evaluations in line-search for projected methods.

  Real              scale0_; ///< Scale for inexact gradient computation.
  Real              scale1_; ///< Scale for inexact gradient computation.

  /** \brief Update gradient to iteratively satisfy inexactness condition.

      This function attempts to ensure that the inexact gradient condition,
      \f[
         \|g_k-\nabla J(x_k)\|_{\mathcal{X}} \le \kappa_1\min\{\,\|g_k\|_{\mathcal{X}},\,\Delta_k\,\},
      \f]
      is satisfied.  This function works under the assumption that the gradient function returns 
      a gradient approximation which satisfies the error tolerance prescribed by the tol input 
      parameter.  
      @param[in]      x          is the current optimization variable.
      @param[in]      obj        is the objective function.
      @param[in]      con        is the bound constraint.
      @param[in,out]  algo_state is the algorithm state.
  */
  void updateGradient( Vector<Real> &x, Objective<Real> &obj, BoundConstraint<Real> &con, 
                       AlgorithmState<Real> &algo_state ) {
    Teuchos::RCP<StepState<Real> > state = Step<Real>::getState();
    if ( this->useInexact_[1] ) {
      Real c = this->scale0_*std::max(1.e-2,std::min(1.0,1.e4*algo_state.gnorm));
      Real gtol1  = c*(state->searchSize);
      Real gtol0  = this->scale1_*gtol1 + 1.0;
      while ( gtol0 > gtol1*this->scale1_ ) {
        obj.gradient(*(state->gradientVec),x,gtol1);
        algo_state.gnorm = this->computeCriticalityMeasure(*(state->gradientVec),x,con);
        gtol0 = gtol1;
        c = this->scale0_*std::max(1.e-2,std::min(1.0,1.e4*algo_state.gnorm));
        gtol1 = c*std::min(algo_state.gnorm,state->searchSize);
      }
      algo_state.ngrad++;
    }
    else {
      Real gtol = std::sqrt(ROL_EPSILON);
      obj.gradient(*(state->gradientVec),x,gtol);
      algo_state.ngrad++;
      algo_state.gnorm = this->computeCriticalityMeasure(*(state->gradientVec),x,con);
    }
  }

  /** \brief Compute the criticality measure.

      This function computes either the norm of the gradient projected onto the tangent cone or 
      the norm of \f$x_k - P_{[a,b]}(x_k-g_k)\f$.
       @param[in]       g     is the current gradient.
       @param[in]       x     is the current iterate.
       @param[in]       con   is the bound constraint.
  */
  Real computeCriticalityMeasure( const Vector<Real> &g, const Vector<Real> &x, BoundConstraint<Real> &con ) {
    if ( con.isActivated() ) {
      if ( this->useProjectedGrad_ ) {
        xnew_->set(g);
        con.computeProjectedGradient( *xnew_, x );
      }
      else {
        xnew_->set(x);
        xnew_->axpy(-1.0,g);
        con.project(*xnew_);
        xnew_->axpy(-1.0,x);
      }
      return xnew_->norm();
    }
    else {
      return g.norm();
    }
  }

public:

  virtual ~TrustRegionStep() {}

  /** \brief Constructor.

      Standard constructor to build a TrustRegionStep object.  Algorithmic 
      specifications are passed in through a Teuchos::ParameterList.

      @param[in]     parlist    is a parameter list containing algorithmic specifications
  */
  TrustRegionStep( Teuchos::ParameterList & parlist ) : Step<Real>() {
    Teuchos::RCP<StepState<Real> > step_state = Step<Real>::getState();

    // Enumerations
    etr_   = StringToETrustRegion(parlist.get("Trust-Region Subproblem Solver Type","Cauchy Point"));  
    esec_  = StringToESecant(parlist.get("Secant Type","Limited-Memory BFGS"));
    // Secant Information
    useSecantPrecond_ = parlist.get("Use Secant Preconditioning", false);
    useSecantHessVec_ = parlist.get("Use Secant Hessian-Times-A-Vector", false);
    // Trust-Region Parameters
    step_state->searchSize = parlist.get("Initial Trust-Region Radius", -1.0);
    // Inexactness Information
    useInexact_.clear();
    useInexact_.push_back(parlist.get("Use Inexact Objective Function", false));
    useInexact_.push_back(parlist.get("Use Inexact Gradient", false));
    useInexact_.push_back(parlist.get("Use Inexact Hessian-Times-A-Vector", false));
    this->scale0_ = parlist.get("Gradient Update Tolerance Scaling",1.e-1);
    this->scale1_ = parlist.get("Gradient Update Relative Tolerance",2.0);     
     
    // Initialize Trust Region Subproblem Solver Object
    useProjectedGrad_ = parlist.get("Use Projected Gradient Criticality Measure", false);
    max_fval_         = parlist.get("Maximum Number of Function Evaluations", 20);
    alpha_init_       = parlist.get("Initial Linesearch Parameter", 1.0);
    trustRegion_      = Teuchos::rcp( new TrustRegion<Real>(parlist) );

    // Secant Parameters
    secant_ = Teuchos::null;
    if ( useSecantPrecond_ || useSecantHessVec_ ) {
      int L      = parlist.get("Maximum Secant Storage",10);
      int BBtype = parlist.get("Barzilai-Borwein Type",1);
      secant_ = getSecant<Real>(esec_,L,BBtype);
    }
  }

  /** \brief Constructor.

      Constructor to build a TrustRegionStep object with a user-defined 
      secant object.  Algorithmic specifications are passed in through 
      a Teuchos::ParameterList.

      @param[in]     secant     is a user-defined secant object
      @param[in]     parlist    is a parameter list containing algorithmic specifications
  */
  TrustRegionStep( Teuchos::RCP<Secant<Real> > &secant, Teuchos::ParameterList &parlist ) 
    : Step<Real>(), secant_(secant) {
    Teuchos::RCP<StepState<Real> > step_state = Step<Real>::getState();

    // Enumerations
    etr_   = StringToETrustRegion(parlist.get("Trust-Region Subproblem Solver Type","Cauchy Point"));  
    esec_  = StringToESecant(parlist.get("Secant Type","Limited-Memory BFGS"));
    // Secant Information
    useSecantPrecond_ = parlist.get("Use Secant Preconditioning", false);
    useSecantHessVec_ = parlist.get("Use Secant Hessian-Times-A-Vector", false);
    // Trust-Region Parameters
    step_state->searchSize = parlist.get("Initial Trust-Region Radius", -1.0);
    // Inexactness Information
    useInexact_.clear();
    useInexact_.push_back(parlist.get("Use Inexact Objective Function", false));
    useInexact_.push_back(parlist.get("Use Inexact Gradient", false));
    useInexact_.push_back(parlist.get("Use Inexact Hessian-Times-A-Vector", false));
    this->scale0_ = parlist.get("Gradient Update Tolerance Scaling",1.e-1);
    this->scale1_ = parlist.get("Gradient Update Relative Tolerance",2.0);     

    // Initialize Trust Region Subproblem Solver Object
    useProjectedGrad_ = parlist.get("Use Projected Gradient Criticality Measure", false);
    max_fval_         = parlist.get("Maximum Number of Function Evaluations", 20);
    alpha_init_       = parlist.get("Initial Linesearch Parameter", 1.0);
    trustRegion_      = Teuchos::rcp( new TrustRegion<Real>(parlist) );
  }

  /** \brief Initialize step.

      This function initializes the information necessary to run the trust-region algorithm.
      @param[in]     x           is the initial guess for the optimization vector.
      @param[in]     obj         is the objective function.
      @param[in]     con         is the bound constraint.
      @param[in]     algo_state  is the algorithm state.
  */
  void initialize( Vector<Real> &x, Objective<Real> &obj, BoundConstraint<Real> &con, 
                   AlgorithmState<Real> &algo_state ) {
    Teuchos::RCP<StepState<Real> > step_state = Step<Real>::getState();

    trustRegion_->initialize(x);

    algo_state.nfval = 0;
    algo_state.ngrad = 0;

    Real htol = std::sqrt(ROL_EPSILON);
    Real ftol = ROL_OVERFLOW; 

    step_state->descentVec  = x.clone();
    step_state->gradientVec = x.clone();

    if ( con.isActivated() ) {
      con.project(x);
      xnew_ = x.clone();
      xold_ = x.clone();
    }
 
    if ( con.isActivated() || secant_ != Teuchos::null ) {
      gp_ = x.clone();
    }

    // Update approximate gradient and approximate objective function.
    obj.update(x,true,algo_state.iter);    
    this->updateGradient(x,obj,con,algo_state);
    algo_state.snorm = 1.e10;
    algo_state.value = obj.value(x,ftol); 
    algo_state.nfval++;

    // Evaluate Objective Function at Cauchy Point
    if ( step_state->searchSize <= 0.0 ) {
      Teuchos::RCP<Vector<Real> > Bg = x.clone();
      if ( this->useSecantHessVec_ ) {
        this->secant_->applyB(*Bg,*(step_state->gradientVec),x);
      }
      else {
        obj.hessVec(*Bg,*(step_state->gradientVec),x,htol);
      }
      Real gBg = Bg->dot(*(step_state->gradientVec));
      Real alpha = 1.0;
      if ( gBg > ROL_EPSILON ) {
        alpha = algo_state.gnorm*algo_state.gnorm/gBg;
      }
      // Evaluate the objective function at the Cauchy point
      Teuchos::RCP<Vector<Real> > cp = x.clone();
      cp->set(*(step_state->gradientVec)); 
      cp->scale(-alpha);
      Teuchos::RCP<Vector<Real> > xcp = x.clone();
      xcp->set(x);
      xcp->plus(*cp);
      obj.update(*xcp);
      Real fnew = obj.value(*xcp,ftol); // MUST DO SOMETHING HERE WITH FTOL
      algo_state.nfval++;
      // Perform cubic interpolation to determine initial trust region radius
      Real gs = (step_state->gradientVec)->dot(*cp);
      Real a  = fnew - algo_state.value - gs - 0.5*alpha*alpha*gBg;
      if ( std::abs(a) < ROL_EPSILON ) { 
        // a = 0 implies the objective is quadratic in the negative gradient direction
        step_state->searchSize = alpha*algo_state.gnorm;
      }
      else {
        Real b  = 0.5*alpha*alpha*gBg;
        Real c  = gs;
        if ( b*b-3.0*a*c > ROL_EPSILON ) {
          // There is at least one critical point
          Real t1 = (-b-std::sqrt(b*b-3.0*a*c))/(3.0*a);
          Real t2 = (-b+std::sqrt(b*b-3.0*a*c))/(3.0*a);
          if ( 6.0*a*t1 + 2.0*b > 0.0 ) {
            // t1 is the minimizer
            step_state->searchSize = t1*alpha*algo_state.gnorm;          
          }
          else {
            // t2 is the minimizer
            step_state->searchSize = t2*alpha*algo_state.gnorm;
          }
        }
        else {
          step_state->searchSize = alpha*algo_state.gnorm;
        }
      }
    }
  }

  /** \brief Compute step.

      Computes a trial step, \f$s_k\f$ by solving the trust-region subproblem.  
      The trust-region subproblem solver is defined by the enum ETrustRegion.  
      @param[out]      s          is the computed trial step
      @param[in]       x          is the current iterate
      @param[in]       obj        is the objective function
      @param[in]       con        are the bound constraints
      @param[in]       algo_state contains the current state of the algorithm
  */
  void compute( Vector<Real> &s, const Vector<Real> &x, Objective<Real> &obj, BoundConstraint<Real> &con, 
                AlgorithmState<Real> &algo_state ) {
    Teuchos::RCP<StepState<Real> > step_state = Step<Real>::getState();

    Real eps = 0.0;
    if ( con.isActivated() ) {
      eps = algo_state.gnorm;
    }
    ProjectedObjective<Real> pObj(obj,con,this->secant_,this->useSecantPrecond_,this->useSecantHessVec_,eps);

    this->CGflag_ = 0;
    this->CGiter_ = 0;
    this->trustRegion_->run(s,algo_state.snorm,step_state->searchSize,this->CGflag_,this->CGiter_,
                            x,*(step_state->gradientVec),algo_state.gnorm,pObj);
  }

  /** \brief Update step, if successful.

      Given a trial step, \f$s_k\f$, this function updates \f$x_{k+1}=x_k+s_k\f$. 
      This function also updates the secant approximation.

      @param[in,out]   x          is the updated iterate
      @param[in]       s          is the computed trial step
      @param[in]       obj        is the objective function
      @param[in]       con        are the bound constraints
      @param[in]       algo_state contains the current state of the algorithm
  */
  void update( Vector<Real> &x, const Vector<Real> &s, Objective<Real> &obj, BoundConstraint<Real> &con, 
               AlgorithmState<Real> &algo_state ) {
    Teuchos::RCP<StepState<Real> > state = Step<Real>::getState();

    Real tol = std::sqrt(ROL_EPSILON);

    Real eps = 0.0;
    if ( con.isActivated() ) {
      eps = algo_state.gnorm;
    }
    ProjectedObjective<Real> pObj(obj,con,this->secant_,this->useSecantPrecond_,this->useSecantHessVec_,eps);

    // Store previous step for constraint computations
    if ( con.isActivated() ) {
      xold_->set(x);
    }

    // Update trust-region information
    this->TRflag_   = 0;
    this->TR_nfval_ = 0;
    this->TR_ngrad_ = 0;
    Real fold = algo_state.value;
    Real fnew = 0.0;
    algo_state.iter++;
    this->trustRegion_->update(x,fnew,state->searchSize,this->TR_nfval_,this->TR_ngrad_,this->TRflag_,
                               s,algo_state.snorm,fold,*(state->gradientVec),algo_state.iter,pObj);
    algo_state.value = fnew;
    algo_state.nfval += this->TR_nfval_;
    algo_state.ngrad += this->TR_ngrad_;

    // Compute new gradient and update secant storage
    if ( this->TRflag_ == 0 || this->TRflag_ == 1 ) {  
      // Perform line search (smoothing) to ensure decrease 
      if ( con.isActivated() ) {
        // Compute new gradient
        obj.gradient(*gp_,x,tol); // MUST DO SOMETHING HERE WITH TOL
        algo_state.ngrad++;
        // Compute smoothed step
        Real alpha = 1.0;
        xnew_->set(x);
        xnew_->axpy(-alpha*this->alpha_init_,*gp_);
        con.project(*xnew_);
        // Compute new objective value
        obj.update(*xnew_,true,algo_state.iter);
        Real ftmp = obj.value(*xnew_,tol); // MUST DO SOMETHING HERE WITH TOL
        algo_state.nfval++;
        // Perform smoothing
        int cnt = 0;
        alpha = 1.0/this->alpha_init_;
        while ( (fnew-ftmp) <= 1.e-4*(fnew-fold) ) { 
          xnew_->set(x);
          xnew_->axpy(-alpha*this->alpha_init_,*gp_);
          con.project(*xnew_);
          obj.update(*xnew_,true,algo_state.iter);
          ftmp = obj.value(*xnew_,tol); // MUST DO SOMETHING HERE WITH TOL
          algo_state.nfval++;
          if ( cnt >= this->max_fval_ ) {
            break;
          }
          alpha *= 0.5;
          cnt++;
        }
        // Store objective function and iteration information
        fnew = ftmp;
        x.set(*xnew_);
      }

      // Store previous gradient for secant update
      if ( this->secant_ != Teuchos::null ) {
        gp_->set(*(state->gradientVec));
      }

      // Update objective function and approximate model
      //obj.update(x,true,algo_state.iter);
      this->updateGradient(x,obj,con,algo_state);

      // Update secant information
      if ( this->secant_ != Teuchos::null ) {
        if ( con.isActivated() ) { // Compute new constrained step
          xnew_->set(x);
          xnew_->axpy(-1.0,*xold_);
          secant_->update(*(state->gradientVec),*gp_,*xnew_,algo_state.snorm,algo_state.iter+1);
        }
        else {
          secant_->update(*(state->gradientVec),*gp_,s,algo_state.snorm,algo_state.iter+1);
        }
      }

      // Update algorithm state
      (algo_state.iterateVec)->set(x);
    }    
  }

  /** \brief Print iterate header.

      This function produces a string containing header information.
  */
  std::string printHeader( void ) const  {
    std::stringstream hist;
    hist << "  ";
    hist << std::setw(6)  << std::left << "iter";
    hist << std::setw(15) << std::left << "value";
    hist << std::setw(15) << std::left << "gnorm";
    hist << std::setw(15) << std::left << "snorm";
    hist << std::setw(15) << std::left << "delta";
    hist << std::setw(10) << std::left << "#fval";
    hist << std::setw(10) << std::left << "#grad";
    hist << std::setw(10) << std::left << "tr_flag";
    if ( this->etr_ == TRUSTREGION_TRUNCATEDCG ) {
      hist << std::setw(10) << std::left << "iterCG";
      hist << std::setw(10) << std::left << "flagCG";
    }
    hist << "\n";
    return hist.str();
  }

  /** \brief Print step name.

      This function produces a string containing the algorithmic step information.
  */
  std::string printName( void ) const {
    std::stringstream hist;
    hist << "\n" << ETrustRegionToString(this->etr_) << " Trust-Region solver";
    if ( this->useSecantPrecond_ || this->useSecantHessVec_ ) {
      if ( this->useSecantPrecond_ && !this->useSecantHessVec_ ) {
        hist << " with " << ESecantToString(this->esec_) << " preconditioning\n";
      }
      else if ( !this->useSecantPrecond_ && this->useSecantHessVec_ ) {
        hist << " with " << ESecantToString(this->esec_) << " Hessian approximation\n";
      }
      else {
        hist << " with " << ESecantToString(this->esec_) << " preconditioning and Hessian approximation\n";
      }
    }
    else {
      hist << "\n";
    }
    return hist.str();
  }

  /** \brief Print iterate status.

      This function prints the iteration status.

      @param[in]     algo_state    is the current state of the algorithm
      @param[in]     printHeader   if ste to true will print the header at each iteration
  */
  std::string print( AlgorithmState<Real> & algo_state, bool printHeader = false ) const  {
    const Teuchos::RCP<const StepState<Real> >& step_state = Step<Real>::getStepState();

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
      hist << std::setw(6)  << std::left << algo_state.iter;
      hist << std::setw(15) << std::left << algo_state.value;
      hist << std::setw(15) << std::left << algo_state.gnorm;
      hist << std::setw(15) << std::left << " "; 
      hist << std::setw(15) << std::left << step_state->searchSize; 
      hist << "\n";
    }
    else {
      hist << "  "; 
      hist << std::setw(6)  << std::left << algo_state.iter;  
      hist << std::setw(15) << std::left << algo_state.value; 
      hist << std::setw(15) << std::left << algo_state.gnorm; 
      hist << std::setw(15) << std::left << algo_state.snorm; 
      hist << std::setw(15) << std::left << step_state->searchSize; 
      hist << std::setw(10) << std::left << algo_state.nfval;              
      hist << std::setw(10) << std::left << algo_state.ngrad;              
      hist << std::setw(10) << std::left << this->TRflag_;              
      if ( this->etr_ == TRUSTREGION_TRUNCATEDCG ) {
        hist << std::setw(10) << std::left << this->CGiter_;
        hist << std::setw(10) << std::left << this->CGflag_;
      }
      hist << "\n";
    }
    return hist.str();
  }

}; // class Step

} // namespace ROL

#endif
