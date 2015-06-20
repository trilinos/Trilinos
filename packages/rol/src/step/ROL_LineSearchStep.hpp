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

#ifndef ROL_LINESEARCHSTEP_H
#define ROL_LINESEARCHSTEP_H

#include "ROL_Types.hpp"
#include "ROL_HelperFunctions.hpp"

#include "ROL_Step.hpp"
#include "ROL_Secant.hpp"
#include "ROL_Krylov.hpp"
#include "ROL_NonlinearCG.hpp"
#include "ROL_LineSearch.hpp"
#include "ROL_ProjectedHessian.hpp"
#include "ROL_ProjectedPreconditioner.hpp"

#include <sstream>
#include <iomanip>

/** @ingroup step_group
    \class ROL::LineSearchStep
    \brief Provides the interface to compute optimization steps
           with line search.

    Suppose \f$\mathcal{X}\f$ is a Hilbert space of 
    functions mapping \f$\Xi\f$ to \f$\mathbb{R}\f$.  For example, 
    \f$\Xi\subset\mathbb{R}^n\f$ and \f$\mathcal{X}=L^2(\Xi)\f$ or 
    \f$\Xi = \{1,\ldots,n\}\f$ and \f$\mathcal{X}=\mathbb{R}^n\f$. We 
    assume \f$f:\mathcal{X}\to\mathbb{R}\f$ is twice-continuously Fr&eacute;chet 
    differentiable and \f$a,\,b\in\mathcal{X}\f$ with \f$a\le b\f$ almost 
    everywhere in \f$\Xi\f$.  Note that these line-search algorithms will also work 
    with secant approximations of the Hessian. 
    This step applies to unconstrained and bound constrained optimization problems,
    \f[
        \min_x\quad f(x) \qquad\text{and}\qquad \min_x\quad f(x)\quad\text{s.t.}\quad a\le x\le b,
    \f]
    respectively.  

    For unconstrained problems, given the \f$k\f$-th iterate \f$x_k\f$ and a descent direction
    \f$s_k\f$, the line search approximately minimizes the 1D objective function 
    \f$\phi_k(t) = f(x_k + t s_k)\f$.  The approximate minimizer \f$t\f$ must satisfy 
    sufficient decrease and curvature conditions into guarantee global convergence.  The 
    sufficient decrease condition (often called the Armijo condition) is 
    \f[
       \phi_k(t) \le \phi_k(0) + c_1 t \phi_k'(0) \qquad\iff\qquad
       f(x_k+ts_k) \le f(x_k) + c_1 t \langle \nabla f(x_k),s_k\rangle_{\mathcal{X}}
    \f]
    where \f$0 < c_1 < 1\f$.  The curvature conditions implemented in ROL include:

    <CENTER>
    | Name              | Condition                                                     | Parameters |
    | :---------------- | :-----------------------------------------------------------: | :---------------------: |
    | Wolfe             | \f$\phi_k'(t) \ge c_2\phi_k'(0)\f$                            | \f$c_1<c_2<1\f$         |
    | Strong Wolfe      | \f$\left|\phi_k'(t)\right| \le c_2 \left|\phi_k'(0)\right|\f$ | \f$c_1<c_2<1\f$         |
    | Generalized Wolfe | \f$c_2\phi_k'(0)\le \phi_k'(t)\le-c_3\phi_k'(0)\f$            | \f$0<c_3<1\f$           |
    | Approximate Wolfe | \f$c_2\phi_k'(0)\le \phi_k'(t)\le(2c_1-1)\phi_k'(0)\f$        | \f$c_1<c_2<1\f$         |
    | Goldstein         | \f$\phi_k(0)+(1-c_1)t\phi_k'(0)\le \phi_k(t)\f$               | \f$0<c_1<\frac{1}{2}\f$ |
    </CENTER>
    
    Note that \f$\phi_k'(t) = \langle \nabla f(x_k+ts_k),s_k\rangle_{\mathcal{X}}\f$.

    For bound constrained problems, the line-search step performs a projected search.  That is, 
    the line search approximately minimizes the 1D objective function 
    \f$\phi_k(t) = f(P_{[a,b]}(x_k+ts_k))\f$ where \f$P_{[a,b]}\f$ denotes the projection onto 
    the upper and lower bounds.  Such line-search algorithms correspond to projected gradient 
    and projected Newton-type algorithms.  

    For projected methods, we require the notion of an active set of an iterate \f$x_k\f$, 
    \f[
       \mathcal{A}_k = \{\, \xi\in\Xi\,:\,x_k(\xi) = a(\xi)\,\}\cap
                       \{\, \xi\in\Xi\,:\,x_k(\xi) = b(\xi)\,\}.
    \f]
    Given \f$\mathcal{A}_k\f$ and a search direction \f$s_k\f$, we define the binding set as
    \f[
       \mathcal{B}_k = \{\, \xi\in\Xi\,:\,x_k(\xi) = a(\xi) \;\text{and}\; s_k(\xi) < 0 \,\}\cap
                       \{\, \xi\in\Xi\,:\,x_k(\xi) = b(\xi) \;\text{and}\; s_k(\xi) > 0 \,\}.
    \f]
    The binding set contains the values of \f$\xi\in\Xi\f$ such that if \f$x_k(\xi)\f$ is on a 
    bound, then \f$(x_k+s_k)(\xi)\f$ will violate bound.  Using these definitions, we can 
    redefine the sufficient decrease and curvature conditions from the unconstrained case to 
    the case of bound constraints.

    LineSearchStep implements a number of algorithms for both bound constrained and unconstrained 
    optimization.  These algorithms are: Steepest descent; Nonlinear CG; Quasi-Newton methods;
    Inexact Newton methods; Newton's method. These methods are chosen through the EDescent enum.
*/


namespace ROL {

template <class Real>
class LineSearchStep : public Step<Real> {
private:

  Teuchos::RCP<Secant<Real> >              secant_;      ///< Secant object (used for quasi-Newton)
  Teuchos::RCP<Krylov<Real> >              krylov_;      ///< Krylov solver object (used for inexact Newton)
  Teuchos::RCP<NonlinearCG<Real> >         nlcg_;        ///< Nonlinear CG object (used for nonlinear CG)
  Teuchos::RCP<LineSearch<Real> >          lineSearch_;  ///< Line-search object

  Teuchos::RCP<ProjectedHessian<Real> >             hessian_;
  Teuchos::RCP<ProjectedPreconditioner<Real> >      precond_;

  Teuchos::RCP<Vector<Real> > d_;
  Teuchos::RCP<Vector<Real> > gp_; 

  int iterKrylov_; ///< Number of Krylov iterations (used for inexact Newton)
  int flagKrylov_; ///< Termination flag for Krylov method (used for inexact Newton)

  ELineSearch         els_;   ///< enum determines type of line search
  ENonlinearCG        enlcg_; ///< enum determines type of nonlinear CG
  ECurvatureCondition econd_; ///< enum determines type of curvature condition
  EDescent            edesc_; ///< enum determines type of descent step
  ESecant             esec_;  ///< enum determines type of secant approximation
  EKrylov             ekv_;   ///< enum determines type of Krylov solver
 
  int ls_nfval_; ///< Number of function evaluations during line search
  int ls_ngrad_; ///< Number of gradient evaluations during line search
 
  bool useSecantHessVec_; ///< Whether or not a secant approximation is used for Hessian-times-a-vector
  bool useSecantPrecond_; ///< Whether or not a secant approximation is used for preconditioning inexact Newton

  bool useProjectedGrad_; ///< Whether or not to use to the projected gradient criticality measure

  std::vector<bool> useInexact_; ///< Flags for inexact objective function, gradient, and Hessian evaluation

  bool              softUp_;

  void LineSearchFactory(Teuchos::ParameterList &parlist) {
    switch(els_) {
      case LINESEARCH_ITERATIONSCALING: 
        lineSearch_ = Teuchos::rcp( new IterationScaling<Real>(parlist) );     break;
      case LINESEARCH_PATHBASEDTARGETLEVEL: 
        lineSearch_ = Teuchos::rcp( new PathBasedTargetLevel<Real>(parlist) ); break;
      case LINESEARCH_BACKTRACKING:
        lineSearch_ = Teuchos::rcp( new BackTracking<Real>(parlist) );         break;
      case LINESEARCH_BISECTION:
        lineSearch_ = Teuchos::rcp( new Bisection<Real>(parlist) );            break;
      case LINESEARCH_BRENTS:
        lineSearch_ = Teuchos::rcp( new Brents<Real>(parlist) );               break;
      case LINESEARCH_GOLDENSECTION:
        lineSearch_ = Teuchos::rcp( new GoldenSection<Real>(parlist) );        break;
      case LINESEARCH_CUBICINTERP:
      default:
        lineSearch_ = Teuchos::rcp( new CubicInterp<Real>(parlist) );          break;
    }
  }

  void KrylovFactory(Teuchos::ParameterList &parlist) {
    Real absTol = parlist.get("Absolute Krylov Tolerance", 1.e-4);
    Real relTol = parlist.get("Relative Krylov Tolerance", 1.e-2);
    int maxit   = parlist.get("Maximum Number of Krylov Iterations", 20);
    switch(ekv_) {
      case KRYLOV_CR:
        krylov_ = Teuchos::rcp( new ConjugateResiduals<Real>(absTol,relTol,maxit,useInexact_[2]) ); break;
      case KRYLOV_CG:
      default:
        krylov_ = Teuchos::rcp( new ConjugateGradients<Real>(absTol,relTol,maxit,useInexact_[2]) ); break;
    }
  }

public:

  virtual ~LineSearchStep() {}

  /** \brief Constructor.

      Standard constructor to build a LineSearchStep object.  Algorithmic 
      specifications are passed in through a Teuchos::ParameterList.

      @param[in]     parlist    is a parameter list containing algorithmic specifications
  */
  LineSearchStep( Teuchos::ParameterList &parlist ) : Step<Real>() {
    // Enumerations
    edesc_ = StringToEDescent(parlist.get("Descent Type","Quasi-Newton Method"));
    enlcg_ = StringToENonlinearCG(parlist.get("Nonlinear CG Type","Hagar-Zhang"));
    els_   = StringToELineSearch(parlist.get("Linesearch Type","Cubic Interpolation"));
    econd_ = StringToECurvatureCondition(parlist.get("Linesearch Curvature Condition","Strong Wolfe Conditions"));
    esec_  = StringToESecant(parlist.get("Secant Type","Limited-Memory BFGS"));
    ekv_   = StringToEKrylov(parlist.get("Krylov Type","Conjugate Gradients"));

    // Inexactness Information
    useInexact_.clear();
    useInexact_.push_back(parlist.get("Use Inexact Objective Function", false));
    useInexact_.push_back(parlist.get("Use Inexact Gradient", false));
    useInexact_.push_back(parlist.get("Use Inexact Hessian-Times-A-Vector", false));
     
    // Initialize Linesearch Object
    useProjectedGrad_ = parlist.get("Use Projected Gradient Criticality Measure", false);
    LineSearchFactory(parlist);

    // Changing Objective Functions
    softUp_ = parlist.get("Variable Objective Function",false);

    // Initialize Krylov Object
    useSecantHessVec_ = parlist.get("Use Secant Hessian-Times-A-Vector", false);
    useSecantPrecond_ = parlist.get("Use Secant Preconditioning", false);
    krylov_ = Teuchos::null;
    iterKrylov_ = 0;
    flagKrylov_ = 0;
    if ( edesc_ == DESCENT_NEWTONKRYLOV ) {
      KrylovFactory(parlist);
    }

    // Initialize Secant Object
    secant_ = Teuchos::null;
    if ( edesc_ == DESCENT_SECANT || (edesc_ == DESCENT_NEWTONKRYLOV && useSecantPrecond_) ) {
      int L      = parlist.get("Maximum Secant Storage",10);
      int BBtype = parlist.get("Barzilai-Borwein Type",1);
      secant_ = getSecant<Real>(esec_,L,BBtype);
    }
    if ( edesc_ == DESCENT_SECANT ) {
      useSecantHessVec_ = true;
    }

    // Initialize Nonlinear CG Object
    nlcg_ = Teuchos::null;
    if ( edesc_ == DESCENT_NONLINEARCG ) {
      nlcg_ = Teuchos::rcp( new NonlinearCG<Real>(enlcg_) );
    }

    ls_nfval_ = 0;
    ls_ngrad_ = 0;
  }

  /** \brief Constructor.

      Standard constructor to build a LineSearchStep object.  Algorithmic 
      specifications are passed in through a Teuchos::ParameterList.

      @param[in]     lineSearch is a user-defined line search object
      @param[in]     parlist    is a parameter list containing algorithmic specifications
  */
  LineSearchStep(Teuchos::RCP<LineSearch<Real> > &lineSearch, Teuchos::ParameterList &parlist) 
    : Step<Real>(), lineSearch_(lineSearch) {
    // Enumerations
    edesc_ = StringToEDescent(parlist.get("Descent Type","Quasi-Newton Method"));
    enlcg_ = StringToENonlinearCG(parlist.get("Nonlinear CG Type","Hagar-Zhang"));
    econd_ = StringToECurvatureCondition(parlist.get("Linesearch Curvature Condition","Strong Wolfe Conditions"));
    esec_  = StringToESecant(parlist.get("Secant Type","Limited-Memory BFGS"));
    ekv_   = StringToEKrylov(parlist.get("Krylov Type","Conjugate Gradients"));

    // Inexactness Information
    useInexact_.clear();
    useInexact_.push_back(parlist.get("Use Inexact Objective Function", false));
    useInexact_.push_back(parlist.get("Use Inexact Gradient", false));
    useInexact_.push_back(parlist.get("Use Inexact Hessian-Times-A-Vector", false));
     
    // Initialize Linesearch Object
    useProjectedGrad_ = parlist.get("Use Projected Gradient Criticality Measure", false);
    els_ = LINESEARCH_USERDEFINED;

    // Changing Objective Functions
    softUp_ = parlist.get("Variable Objective Function",false);

    // Initialize Krylov Object
    useSecantHessVec_ = parlist.get("Use Secant Hessian-Times-A-Vector", false);
    useSecantPrecond_ = parlist.get("Use Secant Preconditioning", false);
    krylov_ = Teuchos::null;
    iterKrylov_ = 0;
    flagKrylov_ = 0;
    if ( edesc_ == DESCENT_NEWTONKRYLOV ) {
      KrylovFactory(parlist);
    }

    // Initialize Secant Object
    secant_ = Teuchos::null;
    if ( edesc_ == DESCENT_SECANT || (edesc_ == DESCENT_NEWTONKRYLOV && useSecantPrecond_) ) {
      int L      = parlist.get("Maximum Secant Storage",10);
      int BBtype = parlist.get("Barzilai-Borwein Type",1);
      secant_ = getSecant<Real>(esec_,L,BBtype);
    }
    if ( edesc_ == DESCENT_SECANT ) {
      useSecantHessVec_ = true;
    }

    // Initialize Nonlinear CG Object
    nlcg_ = Teuchos::null;
    if ( edesc_ == DESCENT_NONLINEARCG ) {
      nlcg_ = Teuchos::rcp( new NonlinearCG<Real>(enlcg_) );
    }

    ls_nfval_ = 0;
    ls_ngrad_ = 0;
  }

  /** \brief Constructor.

      Constructor to build a LineSearchStep object with a user-defined 
      secant object.  Algorithmic specifications are passed in through 
      a Teuchos::ParameterList.

      @param[in]     secant     is a user-defined secant object
      @param[in]     parlist    is a parameter list containing algorithmic specifications
  */
  LineSearchStep( Teuchos::RCP<Secant<Real> > &secant, Teuchos::ParameterList &parlist ) 
    : Step<Real>(), secant_(secant) {
    // Enumerations
    edesc_ = StringToEDescent(parlist.get("Descent Type","Quasi-Newton Method"));
    enlcg_ = StringToENonlinearCG(parlist.get("Nonlinear CG Type","Hagar-Zhang"));
    els_   = StringToELineSearch(parlist.get("Linesearch Type","Cubic Interpolation"));
    econd_ = StringToECurvatureCondition(parlist.get("Linesearch Curvature Condition","Strong Wolfe Conditions"));
    esec_  = StringToESecant(parlist.get("Secant Type","Limited-Memory BFGS"));
    ekv_   = StringToEKrylov(parlist.get("Krylov Type","Conjugate Gradients"));

    // Inexactness Information
    useInexact_.clear();
    useInexact_.push_back(parlist.get("Use Inexact Objective Function", false));
    useInexact_.push_back(parlist.get("Use Inexact Gradient", false));
    useInexact_.push_back(parlist.get("Use Inexact Hessian-Times-A-Vector", false));

    // Initialize Linesearch Object
    useProjectedGrad_ = parlist.get("Use Projected Gradient Criticality Measure", false);
    LineSearchFactory(parlist);

    // Changing Objective Functions
    softUp_ = parlist.get("Variable Objective Function",false);

    // Initialize Krylov Object
    useSecantHessVec_ = parlist.get("Use Secant Hessian-Times-A-Vector", false);
    useSecantPrecond_ = parlist.get("Use Secant Preconditioner", false);
    krylov_ = Teuchos::null;
    iterKrylov_ = 0;
    flagKrylov_ = 0;
    if ( edesc_ == DESCENT_NEWTONKRYLOV ) {
      KrylovFactory(parlist);
    }

    // Secant Information
    if ( edesc_ == DESCENT_SECANT ) {
      useSecantHessVec_ = true;
    }
     
    // Initialize Nonlinear CG Object
    nlcg_ = Teuchos::null;
    if ( edesc_ == DESCENT_NONLINEARCG ) {
      nlcg_ = Teuchos::rcp( new NonlinearCG<Real>(enlcg_) );
    }

    ls_nfval_ = 0;
    ls_ngrad_ = 0;
  }

  /** \brief Constructor.

      Standard constructor to build a LineSearchStep object.  Algorithmic 
      specifications are passed in through a Teuchos::ParameterList.

      @param[in]     krylov     is a user-defined Krylov object
      @param[in]     parlist    is a parameter list containing algorithmic specifications
  */
  LineSearchStep( Teuchos::RCP<Krylov<Real> > &krylov, Teuchos::ParameterList &parlist ) 
    : Step<Real>(), krylov_(krylov) {
    // Enumerations
    edesc_ = StringToEDescent(parlist.get("Descent Type","Quasi-Newton Method"));
    enlcg_ = StringToENonlinearCG(parlist.get("Nonlinear CG Type","Hagar-Zhang"));
    els_   = StringToELineSearch(parlist.get("Linesearch Type","Cubic Interpolation"));
    econd_ = StringToECurvatureCondition(parlist.get("Linesearch Curvature Condition","Strong Wolfe Conditions"));
    esec_  = StringToESecant(parlist.get("Secant Type","Limited-Memory BFGS"));
    ekv_   = StringToEKrylov(parlist.get("Krylov Type","Conjugate Gradients"));

    // Inexactness Information
    useInexact_.clear();
    useInexact_.push_back(parlist.get("Use Inexact Objective Function", false));
    useInexact_.push_back(parlist.get("Use Inexact Gradient", false));
    useInexact_.push_back(parlist.get("Use Inexact Hessian-Times-A-Vector", false));
     
    // Initialize Linesearch Object
    useProjectedGrad_ = parlist.get("Use Projected Gradient Criticality Measure", false);
    LineSearchFactory(parlist);

    // Changing Objective Functions
    softUp_ = parlist.get("Variable Objective Function",false);

    // Initialize Krylov Object
    useSecantHessVec_ = parlist.get("Use Secant Hessian-Times-A-Vector", false);
    useSecantPrecond_ = parlist.get("Use Secant Preconditioning", false);
    iterKrylov_ = 0;
    flagKrylov_ = 0;

    // Initialize Secant Object
    secant_ = Teuchos::null;
    if ( edesc_ == DESCENT_SECANT || (edesc_ == DESCENT_NEWTONKRYLOV && useSecantPrecond_) ) {
      int L      = parlist.get("Maximum Secant Storage",10);
      int BBtype = parlist.get("Barzilai-Borwein Type",1);
      secant_ = getSecant<Real>(esec_,L,BBtype);
    }
    if ( edesc_ == DESCENT_SECANT ) {
      useSecantHessVec_ = true;
    }

    // Initialize Nonlinear CG Object
    nlcg_ = Teuchos::null;
    if ( edesc_ == DESCENT_NONLINEARCG ) {
      nlcg_ = Teuchos::rcp( new NonlinearCG<Real>(enlcg_) );
    }

    ls_nfval_ = 0;
    ls_ngrad_ = 0;
  }

  /** \brief Constructor.

      Constructor to build a LineSearchStep object with a user-defined 
      secant object.  Algorithmic specifications are passed in through 
      a Teuchos::ParameterList.

      @param[in]     lineSearch is a user-defined line search object
      @param[in]     secant     is a user-defined secant object
      @param[in]     parlist    is a parameter list containing algorithmic specifications
  */
  LineSearchStep( Teuchos::RCP<LineSearch<Real> > &lineSearch, Teuchos::RCP<Secant<Real> > &secant, 
                  Teuchos::ParameterList &parlist ) : Step<Real>(), lineSearch_(lineSearch), secant_(secant) {
    // Enumerations
    edesc_ = StringToEDescent(parlist.get("Descent Type","Quasi-Newton Method"));
    enlcg_ = StringToENonlinearCG(parlist.get("Nonlinear CG Type","Hagar-Zhang"));
    econd_ = StringToECurvatureCondition(parlist.get("Linesearch Curvature Condition","Strong Wolfe Conditions"));
    esec_  = StringToESecant(parlist.get("Secant Type","Limited-Memory BFGS"));
    ekv_   = StringToEKrylov(parlist.get("Krylov Type","Conjugate Gradients"));

    // Inexactness Information
    useInexact_.clear();
    useInexact_.push_back(parlist.get("Use Inexact Objective Function", false));
    useInexact_.push_back(parlist.get("Use Inexact Gradient", false));
    useInexact_.push_back(parlist.get("Use Inexact Hessian-Times-A-Vector", false));

    // Initialize Linesearch Object
    useProjectedGrad_ = parlist.get("Use Projected Gradient Criticality Measure", false);
    els_ = LINESEARCH_USERDEFINED;

    // Changing Objective Functions
    softUp_ = parlist.get("Variable Objective Function",false);

    // Initialize Krylov Object
    useSecantHessVec_ = parlist.get("Use Secant Hessian-Times-A-Vector", false);
    useSecantPrecond_ = parlist.get("Use Secant Preconditioner", false);
    krylov_ = Teuchos::null;
    iterKrylov_ = 0;
    flagKrylov_ = 0;
    if ( edesc_ == DESCENT_NEWTONKRYLOV ) {
      KrylovFactory(parlist);
    }

    // Secant Information
    if ( edesc_ == DESCENT_SECANT ) {
      useSecantHessVec_ = true;
    }
     
    // Initialize Nonlinear CG Object
    nlcg_ = Teuchos::null;
    if ( edesc_ == DESCENT_NONLINEARCG ) {
      nlcg_ = Teuchos::rcp( new NonlinearCG<Real>(enlcg_) );
    }
    
    ls_nfval_ = 0;
    ls_ngrad_ = 0;
  }

  void initialize( Vector<Real> &x, const Vector<Real> &s, const Vector<Real> &g, 
                   Objective<Real> &obj, BoundConstraint<Real> &con, 
                   AlgorithmState<Real> &algo_state ) {
    Step<Real>::initialize(x,s,g,obj,con,algo_state);
    Teuchos::RCP<StepState<Real> > step_state = Step<Real>::getState();
    lineSearch_->initialize(x, s, *(step_state->gradientVec),obj,con);
    if ( edesc_ == DESCENT_NEWTONKRYLOV || edesc_ == DESCENT_NEWTON || edesc_ == DESCENT_SECANT ) {
      Teuchos::RCP<Objective<Real> > obj_ptr = Teuchos::rcp(&obj, false);
      Teuchos::RCP<BoundConstraint<Real> > con_ptr = Teuchos::rcp(&con, false);
      hessian_ = Teuchos::rcp(
        new ProjectedHessian<Real>(secant_,obj_ptr,con_ptr,algo_state.iterateVec,step_state->gradientVec,
                                   useSecantHessVec_));
      precond_ = Teuchos::rcp(
        new ProjectedPreconditioner<Real>(secant_,obj_ptr,con_ptr,algo_state.iterateVec,
                                          step_state->gradientVec,useSecantPrecond_));
    }
    if ( con.isActivated() ) {
      d_ = s.clone();
    }
    if ( con.isActivated() || edesc_ == DESCENT_SECANT 
                           || (edesc_ == DESCENT_NEWTONKRYLOV && useSecantPrecond_) ) {
      gp_ = g.clone();
    }
  }


  /** \brief Compute step.

      Computes a trial step, \f$s_k\f$ as defined by the enum EDescent.  Once the 
      trial step is determined, this function determines an approximate minimizer 
      of the 1D function \f$\phi_k(t) = f(x_k+ts_k)\f$.  This approximate 
      minimizer must satisfy sufficient decrease and curvature conditions.

      @param[out]      s          is the computed trial step
      @param[in]       x          is the current iterate
      @param[in]       obj        is the objective function
      @param[in]       con        are the bound constraints
      @param[in]       algo_state contains the current state of the algorithm
  */
  void compute( Vector<Real> &s, const Vector<Real> &x, Objective<Real> &obj, BoundConstraint<Real> &con, 
                AlgorithmState<Real> &algo_state ) {
    Teuchos::RCP<StepState<Real> > step_state = Step<Real>::getState();

    Real tol = std::sqrt(ROL_EPSILON);

    // Set active set parameter
    Real eps = 0.0;
    if ( con.isActivated() ) {
      eps = algo_state.gnorm;
    }
    lineSearch_->setData(eps);
    if ( hessian_ != Teuchos::null ) {
      hessian_->setData(eps);
    }
    if ( precond_ != Teuchos::null ) {
      precond_->setData(eps);
    }

    // Compute step s
    switch(edesc_) {
      case DESCENT_NEWTONKRYLOV:
        flagKrylov_ = 0;
        krylov_->run(s,*hessian_,*(step_state->gradientVec),*precond_,iterKrylov_,flagKrylov_);
        break;
      case DESCENT_NEWTON:
      case DESCENT_SECANT:
        hessian_->applyInverse(s,*(step_state->gradientVec),tol);
        break;
      case DESCENT_NONLINEARCG:
        nlcg_->run(s,*(step_state->gradientVec),x,obj);
        break;
      case DESCENT_STEEPEST:
        s.set(step_state->gradientVec->dual());
        break;
      default: break;
    }

    // Compute g.dot(s)
    Real gs = 0.0;
    if ( !con.isActivated() ) {
      gs = -s.dot((step_state->gradientVec)->dual());
    }
    else {
      if ( edesc_ == DESCENT_STEEPEST ) {
        d_->set(x);
        d_->axpy(-1.0,s);
        con.project(*d_);
        d_->scale(-1.0);
        d_->plus(x);
        //d->set(s);
        //con.pruneActive(*d,s,x,eps);
        //con.pruneActive(*d,*(step_state->gradientVec),x,eps);
        gs = -d_->dot((step_state->gradientVec)->dual());
      }
      else {
        d_->set(s);
        con.pruneActive(*d_,*(step_state->gradientVec),x,eps);
        gs = -d_->dot((step_state->gradientVec)->dual());
        d_->set(x);
        d_->axpy(-1.0,(step_state->gradientVec)->dual());
        con.project(*d_);
        d_->scale(-1.0);
        d_->plus(x);
        con.pruneInactive(*d_,*(step_state->gradientVec),x,eps);
        gs -= d_->dot((step_state->gradientVec)->dual());
      }
    }

    // Check if s is a descent direction i.e., g.dot(s) < 0
    if ( gs >= 0.0 || (flagKrylov_ == 2 && iterKrylov_ <= 1) ) {
      s.set((step_state->gradientVec)->dual());
      if ( con.isActivated() ) {
        d_->set(s);
        con.pruneActive(*d_,s,x);
        gs = -d_->dot((step_state->gradientVec)->dual());
      }
      else {
        gs = -s.dot((step_state->gradientVec)->dual());
      }
    }
    s.scale(-1.0);

    // Perform line search
    Real fnew  = algo_state.value;
    ls_nfval_ = 0;
    ls_ngrad_ = 0;
    lineSearch_->run(step_state->searchSize,fnew,ls_nfval_,ls_ngrad_,gs,s,x,obj,con);
    algo_state.nfval += ls_nfval_;
    algo_state.ngrad += ls_ngrad_;

    // Compute get scaled descent direction
    s.scale(step_state->searchSize);
    if ( con.isActivated() ) {
      s.plus(x);
      con.project(s);
      s.axpy(-1.0,x);
    }

    // Update step state information
    (step_state->descentVec)->set(s);

    // Update algorithm state information
    algo_state.snorm = s.norm();
    algo_state.value = fnew;
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
    Real tol = std::sqrt(ROL_EPSILON);
    Teuchos::RCP<StepState<Real> > step_state = Step<Real>::getState();

    // Update iterate
    algo_state.iter++;
    x.axpy(1.0, s);
    if ( softUp_ ) {
      obj.update(x,true,algo_state.iter);
    }

    // Compute new gradient
    if ( edesc_ == DESCENT_SECANT || 
        (edesc_ == DESCENT_NEWTONKRYLOV && useSecantPrecond_) ) {
      gp_->set(*(step_state->gradientVec));
    }
    obj.gradient(*(step_state->gradientVec),x,tol);
    algo_state.ngrad++;

    // Update Secant Information
    if ( edesc_ == DESCENT_SECANT || 
        (edesc_ == DESCENT_NEWTONKRYLOV && useSecantPrecond_) ) {
      secant_->update(*(step_state->gradientVec),*gp_,s,algo_state.snorm,algo_state.iter+1);
    }

    // Update algorithm state
    (algo_state.iterateVec)->set(x);
    if ( con.isActivated() ) {
      if ( useProjectedGrad_ ) {
        gp_->set(*(step_state->gradientVec));
        con.computeProjectedGradient( *gp_, x );
        algo_state.gnorm = gp_->norm();
      }
      else {
        d_->set(x);
        d_->axpy(-1.0,(step_state->gradientVec)->dual());
        con.project(*d_);
        d_->axpy(-1.0,x);
        algo_state.gnorm = d_->norm();
      }
    }
    else {
      algo_state.gnorm = (step_state->gradientVec)->norm();
    }
  }

  /** \brief Print iterate header.

      This function produces a string containing header information.
  */
  std::string printHeader( void ) const  {
    std::stringstream hist;
    hist << "  ";
    hist << std::setw(6) << std::left << "iter";  
    hist << std::setw(15) << std::left << "value";
    hist << std::setw(15) << std::left << "gnorm"; 
    hist << std::setw(15) << std::left << "snorm";
    hist << std::setw(10) << std::left << "#fval";
    hist << std::setw(10) << std::left << "#grad";
    hist << std::setw(10) << std::left << "ls_#fval";
    hist << std::setw(10) << std::left << "ls_#grad";
    if ( edesc_ == DESCENT_NEWTONKRYLOV ) {
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
    hist << "\n" << EDescentToString(edesc_) 
         << " with " << ELineSearchToString(els_) 
         << " Linesearch satisfying " 
         << ECurvatureConditionToString(econd_) << "\n";
    if ( edesc_ == DESCENT_NEWTONKRYLOV ) {
      hist << "Krylov Type: " << EKrylovToString(ekv_) << "\n";
    }
    if ( edesc_ == DESCENT_SECANT || 
        (edesc_ == DESCENT_NEWTONKRYLOV && (useSecantPrecond_ || useSecantHessVec_)) ) {
      hist << "Secant Type: " << ESecantToString(esec_) << "\n";
    }
    if ( edesc_ == DESCENT_NONLINEARCG ) {
      hist << "Nonlinear CG Type: " << ENonlinearCGToString(enlcg_) << "\n";
    }
    return hist.str();
  }

  /** \brief Print iterate status.

      This function prints the iteration status.

      @param[in]     algo_state    is the current state of the algorithm
      @param[in]     printHeader   if ste to true will print the header at each iteration
  */
  std::string print( AlgorithmState<Real> & algo_state, bool print_header = false ) const  {
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
      hist << std::setw(10) << std::left << ls_nfval_;              
      hist << std::setw(10) << std::left << ls_ngrad_;              
      if ( edesc_ == DESCENT_NEWTONKRYLOV ) {
        hist << std::setw(10) << std::left << iterKrylov_;
        hist << std::setw(10) << std::left << flagKrylov_;
      }
      hist << "\n";
    }
    return hist.str();
  }

  // struct StepState (scalars, vectors) map?

  // getState

  // setState

}; // class Step

} // namespace ROL

#endif
