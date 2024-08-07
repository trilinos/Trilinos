// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#ifndef ROL_PRIMALDUALACTIVESETSTEP_H
#define ROL_PRIMALDUALACTIVESETSTEP_H

#include "ROL_Step.hpp"
#include "ROL_Vector.hpp"
#include "ROL_KrylovFactory.hpp"
#include "ROL_Objective.hpp"
#include "ROL_BoundConstraint.hpp"
#include "ROL_Types.hpp"
#include "ROL_Secant.hpp"
#include "ROL_ParameterList.hpp"

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

  ROL::Ptr<Krylov<Real> > krylov_;

  // Krylov Parameters
  int iterCR_;  ///< CR iteration counter
  int flagCR_;  ///< CR termination flag
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
  ROL::Ptr<Vector<Real> > lambda_; ///< Container for dual variables
  ROL::Ptr<Vector<Real> > xlam_;   ///< Container for primal plus dual variables
  ROL::Ptr<Vector<Real> > x0_;     ///< Container for initial priaml variables
  ROL::Ptr<Vector<Real> > xbnd_;   ///< Container for primal variable bounds
  ROL::Ptr<Vector<Real> > As_;     ///< Container for step projected onto active set
  ROL::Ptr<Vector<Real> > xtmp_;   ///< Container for temporary primal storage
  ROL::Ptr<Vector<Real> > res_;    ///< Container for optimality system residual for quadratic model
  ROL::Ptr<Vector<Real> > Ag_;     ///< Container for gradient projected onto active set
  ROL::Ptr<Vector<Real> > rtmp_;   ///< Container for temporary right hand side storage
  ROL::Ptr<Vector<Real> > gtmp_;   ///< Container for temporary gradient storage
 
  // Secant Information
  ESecant esec_;                       ///< Enum for secant type
  ROL::Ptr<Secant<Real> > secant_; ///< Secant object
  bool useSecantPrecond_; 
  bool useSecantHessVec_;

  class HessianPD : public LinearOperator<Real> {
  private:
    const ROL::Ptr<Objective<Real> > obj_;
    const ROL::Ptr<BoundConstraint<Real> > bnd_;
    const ROL::Ptr<Vector<Real> > x_;
    const ROL::Ptr<Vector<Real> > xlam_;
    ROL::Ptr<Vector<Real> > v_;
    Real eps_;
    const ROL::Ptr<Secant<Real> > secant_;
    bool useSecant_;
  public:
    HessianPD(const ROL::Ptr<Objective<Real> > &obj,
              const ROL::Ptr<BoundConstraint<Real> > &bnd,
              const ROL::Ptr<Vector<Real> > &x,
              const ROL::Ptr<Vector<Real> > &xlam,
              const Real eps = 0,
              const ROL::Ptr<Secant<Real> > &secant = ROL::nullPtr,
              const bool useSecant = false )
      : obj_(obj), bnd_(bnd), x_(x), xlam_(xlam),
        eps_(eps), secant_(secant), useSecant_(useSecant) {
      v_ = x_->clone();
      if ( !useSecant || secant == ROL::nullPtr ) {
        useSecant_ = false;
      }
    }
    void apply( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
      v_->set(v);
      bnd_->pruneActive(*v_,*xlam_,eps_);
      if ( useSecant_ ) {
        secant_->applyB(Hv,*v_);
      }
      else {
        obj_->hessVec(Hv,*v_,*x_,tol);
      }
      bnd_->pruneActive(Hv,*xlam_,eps_);
    }
  };

  class PrecondPD : public LinearOperator<Real> {
  private:
    const ROL::Ptr<Objective<Real> > obj_;
    const ROL::Ptr<BoundConstraint<Real> > bnd_;
    const ROL::Ptr<Vector<Real> > x_;
    const ROL::Ptr<Vector<Real> > xlam_;
    ROL::Ptr<Vector<Real> > v_;
    Real eps_;
    const ROL::Ptr<Secant<Real> > secant_;
    bool useSecant_;
  public:
    PrecondPD(const ROL::Ptr<Objective<Real> > &obj,
              const ROL::Ptr<BoundConstraint<Real> > &bnd,
              const ROL::Ptr<Vector<Real> > &x,
              const ROL::Ptr<Vector<Real> > &xlam,
              const Real eps = 0,
              const ROL::Ptr<Secant<Real> > &secant = ROL::nullPtr,
              const bool useSecant = false )
      : obj_(obj), bnd_(bnd), x_(x), xlam_(xlam),
        eps_(eps), secant_(secant), useSecant_(useSecant) {
      v_ = x_->dual().clone();
      if ( !useSecant || secant == ROL::nullPtr ) {
        useSecant_ = false;
      }
    }
    void apply( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
      Hv.set(v.dual());
    }
    void applyInverse( Vector<Real> &Hv, const Vector<Real> &v, Real &tol ) const {
      v_->set(v);
      bnd_->pruneActive(*v_,*xlam_,eps_);
      if ( useSecant_ ) {
        secant_->applyH(Hv,*v_);
      }
      else {
        obj_->precond(Hv,*v_,*x_,tol);
      }
      bnd_->pruneActive(Hv,*xlam_,eps_);
    }
  };

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
    Real one(1);
    ROL::Ptr<StepState<Real> > step_state = Step<Real>::getState();
    obj.gradient(*(step_state->gradientVec),x,tol);
    xtmp_->set(x);
    xtmp_->axpy(-one,(step_state->gradientVec)->dual());
    con.project(*xtmp_);
    xtmp_->axpy(-one,x);
    return xtmp_->norm();
  }

public:
  /** \brief Constructor.
     
             @param[in]     parlist   is a parameter list containing relevent algorithmic information
             @param[in]     useSecant is a bool which determines whether or not the algorithm uses 
                                      a secant approximation of the Hessian
  */
  PrimalDualActiveSetStep( ROL::ParameterList &parlist ) 
    : Step<Real>::Step(), krylov_(ROL::nullPtr),
      iterCR_(0), flagCR_(0), itol_(0),
      maxit_(0), iter_(0), flag_(0), stol_(0), gtol_(0), scale_(0),
      neps_(-ROL_EPSILON<Real>()), feasible_(false),
      lambda_(ROL::nullPtr), xlam_(ROL::nullPtr), x0_(ROL::nullPtr),
      xbnd_(ROL::nullPtr), As_(ROL::nullPtr), xtmp_(ROL::nullPtr),
      res_(ROL::nullPtr), Ag_(ROL::nullPtr), rtmp_(ROL::nullPtr),
      gtmp_(ROL::nullPtr),
      esec_(SECANT_LBFGS), secant_(ROL::nullPtr), useSecantPrecond_(false),
      useSecantHessVec_(false) {
    Real one(1), oem6(1.e-6), oem8(1.e-8);
    // Algorithmic parameters
    maxit_ = parlist.sublist("Step").sublist("Primal Dual Active Set").get("Iteration Limit",10);
    stol_ = parlist.sublist("Step").sublist("Primal Dual Active Set").get("Relative Step Tolerance",oem8);
    gtol_ = parlist.sublist("Step").sublist("Primal Dual Active Set").get("Relative Gradient Tolerance",oem6);
    scale_ = parlist.sublist("Step").sublist("Primal Dual Active Set").get("Dual Scaling", one);
    // Build secant object
    esec_ = StringToESecant(parlist.sublist("General").sublist("Secant").get("Type","Limited-Memory BFGS"));
    useSecantHessVec_ = parlist.sublist("General").sublist("Secant").get("Use as Hessian", false); 
    useSecantPrecond_ = parlist.sublist("General").sublist("Secant").get("Use as Preconditioner", false);
    if ( useSecantHessVec_ || useSecantPrecond_ ) {
      secant_ = SecantFactory<Real>(parlist);
    }
    // Build Krylov object
    krylov_ = KrylovFactory<Real>(parlist);
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
  using Step<Real>::initialize;
  void initialize( Vector<Real> &x, const Vector<Real> &s, const Vector<Real> &g, 
                   Objective<Real> &obj, BoundConstraint<Real> &con, 
                   AlgorithmState<Real> &algo_state ) {
    ROL::Ptr<StepState<Real> > step_state = Step<Real>::getState();
    Real zero(0), one(1);
    // Initialize state descent direction and gradient storage
    step_state->descentVec  = s.clone();
    step_state->gradientVec = g.clone();
    step_state->searchSize  = zero;
    // Initialize additional storage
    xlam_ = x.clone(); 
    x0_   = x.clone();
    xbnd_ = x.clone();
    As_   = s.clone(); 
    xtmp_ = x.clone(); 
    res_  = g.clone();
    Ag_   = g.clone(); 
    rtmp_ = g.clone(); 
    gtmp_ = g.clone(); 
    // Project x onto constraint set
    con.project(x);
    // Update objective function, get value, and get gradient
    Real tol = std::sqrt(ROL_EPSILON<Real>());
    obj.update(x,true,algo_state.iter);
    algo_state.value = obj.value(x,tol);
    algo_state.nfval++;
    algo_state.gnorm = computeCriticalityMeasure(x,obj,con,tol);
    algo_state.ngrad++;
    // Initialize dual variable
    lambda_ = s.clone(); 
    lambda_->set((step_state->gradientVec)->dual());
    lambda_->scale(-one);
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
  using Step<Real>::compute;
  void compute( Vector<Real> &s, const Vector<Real> &x, Objective<Real> &obj, BoundConstraint<Real> &con, 
                AlgorithmState<Real> &algo_state ) {
    ROL::Ptr<StepState<Real> > step_state = Step<Real>::getState();
    Real zero(0), one(1);
    s.zero();
    x0_->set(x);
    res_->set(*(step_state->gradientVec));
    for ( iter_ = 0; iter_ < maxit_; iter_++ ) {
      /********************************************************************/
      // MODIFY ITERATE VECTOR TO CHECK ACTIVE SET
      /********************************************************************/
      xlam_->set(*x0_);                          // xlam = x0
      xlam_->axpy(scale_,*(lambda_));            // xlam = x0 + c*lambda
      /********************************************************************/
      // PROJECT x ONTO PRIMAL DUAL FEASIBLE SET
      /********************************************************************/
      As_->zero();                               // As   = 0
   
      xbnd_->set(*con.getUpperBound());          // xbnd = u
      xbnd_->axpy(-one,x);                       // xbnd = u - x
      xtmp_->set(*xbnd_);                        // tmp  = u - x
      con.pruneUpperActive(*xtmp_,*xlam_,neps_); // tmp  = I(u - x)
      xbnd_->axpy(-one,*xtmp_);                  // xbnd = A(u - x)
      As_->plus(*xbnd_);                         // As  += A(u - x)

      xbnd_->set(*con.getLowerBound());          // xbnd = l
      xbnd_->axpy(-one,x);                       // xbnd = l - x
      xtmp_->set(*xbnd_);                        // tmp  = l - x
      con.pruneLowerActive(*xtmp_,*xlam_,neps_); // tmp  = I(l - x)
      xbnd_->axpy(-one,*xtmp_);                  // xbnd = A(l - x)
      As_->plus(*xbnd_);                         // As  += A(l - x)
      /********************************************************************/
      // APPLY HESSIAN TO ACTIVE COMPONENTS OF s AND REMOVE INACTIVE
      /********************************************************************/
      itol_ = std::sqrt(ROL_EPSILON<Real>());
      if ( useSecantHessVec_ && secant_ != ROL::nullPtr ) {        // IHAs = H*As
        secant_->applyB(*gtmp_,*As_);
      }
      else {
        obj.hessVec(*gtmp_,*As_,x,itol_);
      }
      con.pruneActive(*gtmp_,*xlam_,neps_);     // IHAs = I(H*As)
      /********************************************************************/
      // SEPARATE ACTIVE AND INACTIVE COMPONENTS OF THE GRADIENT
      /********************************************************************/
      rtmp_->set(*(step_state->gradientVec));    // Inactive components
      con.pruneActive(*rtmp_,*xlam_,neps_);

      Ag_->set(*(step_state->gradientVec));     // Active components
      Ag_->axpy(-one,*rtmp_);
      /********************************************************************/
      // SOLVE REDUCED NEWTON SYSTEM 
      /********************************************************************/
      rtmp_->plus(*gtmp_);
      rtmp_->scale(-one);                        // rhs = -Ig - I(H*As)
      s.zero();
      if ( rtmp_->norm() > zero ) {             
        // Initialize Hessian and preconditioner
        ROL::Ptr<Objective<Real> > obj_ptr = ROL::makePtrFromRef(obj);
        ROL::Ptr<BoundConstraint<Real> > con_ptr = ROL::makePtrFromRef(con);
        ROL::Ptr<LinearOperator<Real> > hessian
          = ROL::makePtr<HessianPD>(obj_ptr,con_ptr,
              algo_state.iterateVec,xlam_,neps_,secant_,useSecantHessVec_);
        ROL::Ptr<LinearOperator<Real> > precond
          = ROL::makePtr<PrecondPD>(obj_ptr,con_ptr,
              algo_state.iterateVec,xlam_,neps_,secant_,useSecantPrecond_);
        //solve(s,*rtmp_,*xlam_,x,obj,con);   // Call conjugate residuals
        krylov_->run(s,*hessian,*rtmp_,*precond,iterCR_,flagCR_);
        con.pruneActive(s,*xlam_,neps_);        // s <- Is
      }
      s.plus(*As_);                             // s = Is + As
      /********************************************************************/
      // UPDATE MULTIPLIER 
      /********************************************************************/
      if ( useSecantHessVec_ && secant_ != ROL::nullPtr ) {
        secant_->applyB(*rtmp_,s);
      }
      else {
        obj.hessVec(*rtmp_,s,x,itol_);
      }
      gtmp_->set(*rtmp_);
      con.pruneActive(*gtmp_,*xlam_,neps_);
      lambda_->set(*rtmp_);
      lambda_->axpy(-one,*gtmp_);
      lambda_->plus(*Ag_);
      lambda_->scale(-one);
      /********************************************************************/
      // UPDATE STEP 
      /********************************************************************/
      x0_->set(x);
      x0_->plus(s);
      res_->set(*(step_state->gradientVec));
      res_->plus(*rtmp_);
      // Compute criticality measure  
      xtmp_->set(*x0_);
      xtmp_->axpy(-one,res_->dual());
      con.project(*xtmp_);
      xtmp_->axpy(-one,*x0_);
//      std::cout << s.norm()               << "  " 
//                << tmp->norm()            << "  " 
//                << res_->norm()           << "  " 
//                << lambda_->norm()  << "  " 
//                << flagCR_          << "  " 
//                << iterCR_          << "\n";
      if ( xtmp_->norm() < gtol_*algo_state.gnorm ) {
        flag_ = 0;
        break;
      }
      if ( s.norm() < stol_*x.norm() ) {
        flag_ = 2;
        break;
      } 
    }
    if ( iter_ == maxit_ ) {
      flag_ = 1;
    }
    else {
      iter_++;
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
  using Step<Real>::update;
  void update( Vector<Real> &x, const Vector<Real> &s, Objective<Real> &obj, BoundConstraint<Real> &con,
               AlgorithmState<Real> &algo_state ) {
    ROL::Ptr<StepState<Real> > step_state = Step<Real>::getState();
    step_state->SPiter = (maxit_ > 1) ? iter_ : iterCR_;
    step_state->SPflag = (maxit_ > 1) ? flag_ : flagCR_;

    x.plus(s);
    feasible_ = con.isFeasible(x);
    algo_state.snorm = s.norm();
    algo_state.iter++;
    Real tol = std::sqrt(ROL_EPSILON<Real>());
    obj.update(x,true,algo_state.iter);
    algo_state.value = obj.value(x,tol);
    algo_state.nfval++;
    
    if ( secant_ != ROL::nullPtr ) {
      gtmp_->set(*(step_state->gradientVec));
    }
    algo_state.gnorm = computeCriticalityMeasure(x,obj,con,tol);
    algo_state.ngrad++;

    if ( secant_ != ROL::nullPtr ) {
      secant_->updateStorage(x,*(step_state->gradientVec),*gtmp_,s,algo_state.snorm,algo_state.iter+1);
    }
    (algo_state.iterateVec)->set(x);
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
    if ( maxit_ > 1 ) {
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
  virtual std::string print( AlgorithmState<Real> &algo_state, bool print_header = false ) const {
    std::stringstream hist;
    hist << std::scientific << std::setprecision(6);
    if ( algo_state.iter == 0 ) {
      hist << printName();
    }
    if ( print_header ) {
      hist << printHeader();
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
      if ( maxit_ > 1 ) {
        hist << std::setw(10) << std::left << iter_;
        hist << std::setw(10) << std::left << flag_;
      }
      else {
        hist << std::setw(10) << std::left << iterCR_;
        hist << std::setw(10) << std::left << flagCR_;
      }
      if ( feasible_ ) {
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
//    Real rtol   = std::min(tol1_,tol2_*rnorm);
//    itol_ = std::sqrt(ROL_EPSILON<Real>());
//    sol.zero();
//
//    ROL::Ptr<Vector<Real> > res = rhs.clone();
//    res->set(rhs);
//
//    ROL::Ptr<Vector<Real> > v = x.clone();
//    con.pruneActive(*res,xlam,neps_);
//    obj.precond(*v,*res,x,itol_);
//    con.pruneActive(*v,xlam,neps_);
//
//    ROL::Ptr<Vector<Real> > p = x.clone();
//    p->set(*v);
//
//    ROL::Ptr<Vector<Real> > Hp = x.clone();
//
//    iterCR_ = 0;
//    flagCR_ = 0;
//
//    Real kappa = 0.0, beta  = 0.0, alpha = 0.0, tmp = 0.0, rv = v->dot(*res);
//
//    for (iterCR_ = 0; iterCR_ < maxitCR_; iterCR_++) {
//      if ( false ) {
//        itol_ = rtol/(maxitCR_*rnorm);
//      }
//      con.pruneActive(*p,xlam,neps_);
//      if ( secant_ == ROL::nullPtr ) {
//        obj.hessVec(*Hp, *p, x, itol_);
//      }
//      else {
//        secant_->applyB( *Hp, *p, x );
//      }
//      con.pruneActive(*Hp,xlam,neps_);
//
//      kappa = p->dot(*Hp);
//      if ( kappa <= 0.0 ) { flagCR_ = 2; break; }
//      alpha = rv/kappa;
//      sol.axpy(alpha,*p);
//
//      res->axpy(-alpha,*Hp);
//      rnorm = res->norm();
//      if ( rnorm < rtol ) { break; }
//
//      con.pruneActive(*res,xlam,neps_);
//      obj.precond(*v,*res,x,itol_);
//      con.pruneActive(*v,xlam,neps_);
//      tmp  = rv;
//      rv   = v->dot(*res);
//      beta = rv/tmp;
//
//      p->scale(beta);
//      p->axpy(1.0,*v);
//    }
//    if ( iterCR_ == maxitCR_ ) {
//      flagCR_ = 1;
//    }
//    else {
//      iterCR_++;
//    }
//  }


//  /** \brief Apply the inactive components of the Hessian operator.
// 
//             I.e., the components corresponding to \f$\mathcal{I}_k\f$.
//
//             @param[out]       hv     is the result of applying the Hessian at @b x to 
//                                      @b v
//             @param[in]        v      is the direction in which we apply the Hessian
//             @param[in]        x      is the current iteration vector \f$x_k\f$
//             @param[in]        xlam   is the vector \f$x_k + c\lambda_k\f$
//             @param[in]        obj    is the objective function
//             @param[in]        con    are the bound constraints
//  */
//  void applyInactiveHessian(Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, 
//                      const Vector<Real> &xlam, Objective<Real> &obj, BoundConstraint<Real> &con) {
//    ROL::Ptr<Vector<Real> > tmp = v.clone();
//    tmp->set(v);
//    con.pruneActive(*tmp,xlam,neps_);
//    if ( secant_ == ROL::nullPtr ) {
//      obj.hessVec(hv,*tmp,x,itol_);
//    }
//    else {
//      secant_->applyB(hv,*tmp,x);
//    }
//    con.pruneActive(hv,xlam,neps_);
//  }
//
//  /** \brief Apply the inactive components of the preconditioner operator.
//
//             I.e., the components corresponding to \f$\mathcal{I}_k\f$.
//
//             @param[out]       hv     is the result of applying the preconditioner at @b x to 
//                                      @b v
//             @param[in]        v      is the direction in which we apply the preconditioner
//             @param[in]        x      is the current iteration vector \f$x_k\f$
//             @param[in]        xlam   is the vector \f$x_k + c\lambda_k\f$
//             @param[in]        obj    is the objective function
//             @param[in]        con    are the bound constraints
//  */
//  void applyInactivePrecond(Vector<Real> &pv, const Vector<Real> &v, const Vector<Real> &x,
//                      const Vector<Real> &xlam, Objective<Real> &obj, BoundConstraint<Real> &con) {
//    ROL::Ptr<Vector<Real> > tmp = v.clone();
//    tmp->set(v);
//    con.pruneActive(*tmp,xlam,neps_);
//    obj.precond(pv,*tmp,x,itol_);
//    con.pruneActive(pv,xlam,neps_);
//  }
//
//  /** \brief Solve the inactive part of the PDAS optimality system.  
//
//             The inactive PDAS optimality system is 
//             \f[
//                 \nabla^2 f(x_k)_{\mathcal{I}_k,\mathcal{I}_k}s  = 
//                     -\nabla f(x_k)_{\mathcal{I}_k}
//                     -\nabla^2 f(x_k)_{\mathcal{I}_k,\mathcal{A}_k} (s_k)_{\mathcal{A}_k}.
//             \f]
//             Since the inactive part of the Hessian may not be positive definite, we solve 
//             using CR.
//   
//             @param[out]       sol    is the vector containing the solution
//             @param[in]        rhs    is the right-hand side vector
//             @param[in]        xlam   is the vector \f$x_k + c\lambda_k\f$
//             @param[in]        x      is the current iteration vector \f$x_k\f$
//             @param[in]        obj    is the objective function
//             @param[in]        con    are the bound constraints
//  */
//  // Solve the inactive part of the optimality system using conjugate residuals
//  void solve(Vector<Real> &sol, const Vector<Real> &rhs, const Vector<Real> &xlam, const Vector<Real> &x, 
//             Objective<Real> &obj, BoundConstraint<Real> &con) {
//    // Initialize Residual
//    ROL::Ptr<Vector<Real> > res = rhs.clone();
//    res->set(rhs);
//    Real rnorm  = res->norm(); 
//    Real rtol   = std::min(tol1_,tol2_*rnorm);
//    if ( false ) { itol_ = rtol/(maxitCR_*rnorm); }
//    sol.zero();
//
//    // Apply preconditioner to residual r = Mres
//    ROL::Ptr<Vector<Real> > r = x.clone();
//    applyInactivePrecond(*r,*res,x,xlam,obj,con);
//
//    // Initialize direction p = v
//    ROL::Ptr<Vector<Real> > p = x.clone();
//    p->set(*r);
//
//    // Apply Hessian to v
//    ROL::Ptr<Vector<Real> > Hr = x.clone();
//    applyInactiveHessian(*Hr,*r,x,xlam,obj,con);
//
//    // Apply Hessian to p
//    ROL::Ptr<Vector<Real> > Hp  = x.clone();
//    ROL::Ptr<Vector<Real> > MHp = x.clone();
//    Hp->set(*Hr);
//
//    iterCR_ = 0;
//    flagCR_ = 0;
//
//    Real kappa = 0.0, beta  = 0.0, alpha = 0.0, tmp = 0.0, rHr = Hr->dot(*r);
//
//    for (iterCR_ = 0; iterCR_ < maxitCR_; iterCR_++) {
//      // Precondition Hp
//      applyInactivePrecond(*MHp,*Hp,x,xlam,obj,con);
//
//      kappa = Hp->dot(*MHp);  // p' H M H p
//      alpha = rHr/kappa;      // r' M H M r
//      sol.axpy(alpha,*p);     // update step
//      res->axpy(-alpha,*Hp);  // residual
//      r->axpy(-alpha,*MHp);   // preconditioned residual
//      
//      // recompute rnorm and decide whether or not to exit
//      rnorm = res->norm();
//      if ( rnorm < rtol ) { break; }
//
//      // Apply Hessian to v
//      itol_ = rtol/(maxitCR_*rnorm);
//      applyInactiveHessian(*Hr,*r,x,xlam,obj,con);
//
//      tmp  = rHr;
//      rHr  = Hr->dot(*r);
//      beta = rHr/tmp;
//      p->scale(beta);
//      p->axpy(1.0,*r);
//      Hp->scale(beta);
//      Hp->axpy(1.0,*Hr);
//    }
//    if ( iterCR_ == maxitCR_ ) {
//      flagCR_ = 1;
//    }
//    else {
//      iterCR_++;
//    }
//  }
