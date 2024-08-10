// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_LINESEARCHSTEP_H
#define ROL_LINESEARCHSTEP_H

#include "ROL_Types.hpp"
#include "ROL_Step.hpp"
#include "ROL_LineSearch.hpp"

// Unconstrained Methods
#include "ROL_GradientStep.hpp"
#include "ROL_NonlinearCGStep.hpp"
#include "ROL_SecantStep.hpp"
#include "ROL_NewtonStep.hpp"
#include "ROL_NewtonKrylovStep.hpp"

// Bound Constrained Methods
#include "ROL_ProjectedSecantStep.hpp"
#include "ROL_ProjectedNewtonStep.hpp"
#include "ROL_ProjectedNewtonKrylovStep.hpp"

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

  ROL::Ptr<Step<Real> >        desc_;       ///< Unglobalized step object
  ROL::Ptr<Secant<Real> >      secant_;     ///< Secant object (used for quasi-Newton)
  ROL::Ptr<Krylov<Real> >      krylov_;     ///< Krylov solver object (used for inexact Newton)
  ROL::Ptr<NonlinearCG<Real> > nlcg_;       ///< Nonlinear CG object (used for nonlinear CG)
  ROL::Ptr<LineSearch<Real> >  lineSearch_; ///< Line-search object

  ROL::Ptr<Vector<Real> > d_;

  ELineSearch         els_;   ///< enum determines type of line search
  ECurvatureCondition econd_; ///< enum determines type of curvature condition

  bool acceptLastAlpha_;  ///< For backwards compatibility. When max function evaluations are reached take last step

  bool usePreviousAlpha_; ///< If true, use the previously accepted step length (if any) as the new initial step length

  int verbosity_;
  bool computeObj_;
  Real fval_;

  ROL::ParameterList parlist_;

  std::string lineSearchName_;  

  Real GradDotStep(const Vector<Real> &g, const Vector<Real> &s,
                   const Vector<Real> &x,
                   BoundConstraint<Real> &bnd, Real eps = 0) {
    Real gs(0), one(1);
    if (!bnd.isActivated()) {
      gs = s.dot(g.dual());
    }
    else {
      d_->set(s);
      bnd.pruneActive(*d_,g,x,eps);
      gs = d_->dot(g.dual());
      d_->set(x);
      d_->axpy(-one,g.dual());
      bnd.project(*d_);
      d_->scale(-one);
      d_->plus(x);
      bnd.pruneInactive(*d_,g,x,eps);
      gs -= d_->dot(g.dual());
    }
    return gs;
  }

public:

  using Step<Real>::initialize;
  using Step<Real>::compute;
  using Step<Real>::update;

  /** \brief Constructor.

      Standard constructor to build a LineSearchStep object.  Algorithmic 
      specifications are passed in through a ROL::ParameterList.  The
      line-search type, secant type, Krylov type, or nonlinear CG type can
      be set using user-defined objects.

      @param[in]     parlist    is a parameter list containing algorithmic specifications
      @param[in]     lineSearch is a user-defined line search object
      @param[in]     secant     is a user-defined secant object
      @param[in]     krylov     is a user-defined Krylov object
      @param[in]     nlcg       is a user-defined Nonlinear CG object
  */
  LineSearchStep( ROL::ParameterList &parlist,
                  const ROL::Ptr<LineSearch<Real> > &lineSearch = ROL::nullPtr,
                  const ROL::Ptr<Secant<Real> > &secant = ROL::nullPtr,
                  const ROL::Ptr<Krylov<Real> > &krylov = ROL::nullPtr,
                  const ROL::Ptr<NonlinearCG<Real> > &nlcg = ROL::nullPtr )
    : Step<Real>(), desc_(ROL::nullPtr), secant_(secant),
      krylov_(krylov), nlcg_(nlcg), lineSearch_(lineSearch),
      els_(LINESEARCH_USERDEFINED), econd_(CURVATURECONDITION_WOLFE),
      verbosity_(0), computeObj_(true), fval_(0), parlist_(parlist) {
    // Parse parameter list
    ROL::ParameterList& Llist = parlist.sublist("Step").sublist("Line Search");
    ROL::ParameterList& Glist = parlist.sublist("General");
    econd_ = StringToECurvatureCondition(Llist.sublist("Curvature Condition").get("Type","Strong Wolfe Conditions") );
    acceptLastAlpha_ = Llist.get("Accept Last Alpha", false); 
    verbosity_ = Glist.get("Print Verbosity",0);
    computeObj_ = Glist.get("Recompute Objective Function",false);
    // Initialize Line Search
    if (lineSearch_ == ROL::nullPtr) {
      lineSearchName_ = Llist.sublist("Line-Search Method").get("Type","Cubic Interpolation"); 
      els_ = StringToELineSearch(lineSearchName_);
      lineSearch_ = LineSearchFactory<Real>(parlist);
    } 
    else { // User-defined linesearch provided
      lineSearchName_ = Llist.sublist("Line-Search Method").get("User Defined Line-Search Name",
                                                                "Unspecified User Defined Line-Search");
    }

  }

  void initialize( Vector<Real> &x, const Vector<Real> &s, const Vector<Real> &g, 
                   Objective<Real> &obj, BoundConstraint<Real> &bnd, 
                   AlgorithmState<Real> &algo_state ) {
    d_ = x.clone();

    // Initialize unglobalized step
    ROL::ParameterList& list
      = parlist_.sublist("Step").sublist("Line Search").sublist("Descent Method");
    EDescent edesc = StringToEDescent(list.get("Type","Quasi-Newton Method"));
    if (bnd.isActivated()) {
      switch(edesc) {
        case DESCENT_STEEPEST: {
          desc_ = ROL::makePtr<GradientStep<Real>>(parlist_,computeObj_);
          break;
        }
        case DESCENT_NONLINEARCG: {
          desc_ = ROL::makePtr<NonlinearCGStep<Real>>(parlist_,nlcg_,computeObj_);
          break;
        }
        case DESCENT_SECANT: {
          desc_ = ROL::makePtr<ProjectedSecantStep<Real>>(parlist_,secant_,computeObj_);
          break;
        }
        case DESCENT_NEWTON: {
          desc_ = ROL::makePtr<ProjectedNewtonStep<Real>>(parlist_,computeObj_);
          break;
        }
        case DESCENT_NEWTONKRYLOV: {
          desc_ = ROL::makePtr<ProjectedNewtonKrylovStep<Real>>(parlist_,krylov_,secant_,computeObj_);
          break;
        }
        default:
          ROL_TEST_FOR_EXCEPTION(true,std::invalid_argument,
            ">>> (LineSearchStep::Initialize): Undefined descent type!");
      }
    }
    else {
      switch(edesc) {
        case DESCENT_STEEPEST: {
          desc_ = ROL::makePtr<GradientStep<Real>>(parlist_,computeObj_);
          break;
        }
        case DESCENT_NONLINEARCG: {
          desc_ = ROL::makePtr<NonlinearCGStep<Real>>(parlist_,nlcg_,computeObj_);
          break;
        }
        case DESCENT_SECANT: {
          desc_ = ROL::makePtr<SecantStep<Real>>(parlist_,secant_,computeObj_);
          break;
        }
        case DESCENT_NEWTON: {
          desc_ = ROL::makePtr<NewtonStep<Real>>(parlist_,computeObj_);
          break;
        }
        case DESCENT_NEWTONKRYLOV: {
          desc_ = ROL::makePtr<NewtonKrylovStep<Real>>(parlist_,krylov_,secant_,computeObj_);
          break;
        }
        default:
          ROL_TEST_FOR_EXCEPTION(true,std::invalid_argument,
            ">>> (LineSearchStep::Initialize): Undefined descent type!");
      }
    }
    desc_->initialize(x,s,g,obj,bnd,algo_state);

    // Initialize line search
    lineSearch_->initialize(x,s,g,obj,bnd);
    //const ROL::Ptr<const StepState<Real> > desc_state = desc_->getStepState();
    //lineSearch_->initialize(x,s,*(desc_state->gradientVec),obj,bnd);
  }

  /** \brief Compute step.

      Computes a trial step, \f$s_k\f$ as defined by the enum EDescent.  Once the 
      trial step is determined, this function determines an approximate minimizer 
      of the 1D function \f$\phi_k(t) = f(x_k+ts_k)\f$.  This approximate 
      minimizer must satisfy sufficient decrease and curvature conditions.

      @param[out]      s          is the computed trial step
      @param[in]       x          is the current iterate
      @param[in]       obj        is the objective function
      @param[in]       bnd        are the bound constraints
      @param[in]       algo_state contains the current state of the algorithm
  */
  void compute( Vector<Real> &s, const Vector<Real> &x,
                Objective<Real> &obj, BoundConstraint<Real> &bnd, 
                AlgorithmState<Real> &algo_state ) {
    Real zero(0), one(1);
    // Compute unglobalized step
    desc_->compute(s,x,obj,bnd,algo_state);

    // Ensure that s is a descent direction
    // ---> If not, then default to steepest descent
    const ROL::Ptr<const StepState<Real> > desc_state = desc_->getStepState();
    Real gs = GradDotStep(*(desc_state->gradientVec),s,x,bnd,algo_state.gnorm);
    if (gs >= zero) {
      s.set((desc_state->gradientVec)->dual());
      s.scale(-one);
      gs = GradDotStep(*(desc_state->gradientVec),s,x,bnd,algo_state.gnorm);
    }

    // Perform line search
    ROL::Ptr<StepState<Real> > step_state = Step<Real>::getState();
    fval_ = algo_state.value;
    step_state->nfval = 0; step_state->ngrad = 0;
    lineSearch_->setData(algo_state.gnorm,*(desc_state->gradientVec));
    lineSearch_->run(step_state->searchSize,fval_,step_state->nfval,step_state->ngrad,gs,s,x,obj,bnd);

    // Make correction if maximum function evaluations reached
    if(!acceptLastAlpha_) {
      lineSearch_->setMaxitUpdate(step_state->searchSize,fval_,algo_state.value);
    }

    // Compute scaled descent direction
    s.scale(step_state->searchSize);
    if ( bnd.isActivated() ) {
      s.plus(x);
      bnd.project(s);
      s.axpy(static_cast<Real>(-1),x);
    }
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
  void update( Vector<Real> &x, const Vector<Real> &s,
               Objective<Real> &obj, BoundConstraint<Real> &bnd,
               AlgorithmState<Real> &algo_state ) {
    ROL::Ptr<StepState<Real> > step_state = Step<Real>::getState();
    algo_state.nfval += step_state->nfval;
    algo_state.ngrad += step_state->ngrad;
    desc_->update(x,s,obj,bnd,algo_state);
    step_state->flag = desc_->getStepState()->flag;
    step_state->SPiter = desc_->getStepState()->SPiter;
    step_state->SPflag = desc_->getStepState()->SPflag;
    if ( !computeObj_ ) {
      algo_state.value = fval_;
    }
  }

  /** \brief Print iterate header.

      This function produces a string containing header information.
  */
  std::string printHeader( void ) const {
    std::string head = desc_->printHeader();
    head.erase(std::remove(head.end()-3,head.end(),'\n'), head.end());
    std::stringstream hist;
    hist.write(head.c_str(),head.length());
    hist << std::setw(10) << std::left << "ls_#fval";
    hist << std::setw(10) << std::left << "ls_#grad";
    hist << "\n";
    return hist.str();
  }
  
  /** \brief Print step name.

      This function produces a string containing the algorithmic step information.
  */
  std::string printName( void ) const {
    std::string name = desc_->printName();
    std::stringstream hist;
    hist << name;
    hist << "Line Search: " << lineSearchName_;
    hist << " satisfying " << ECurvatureConditionToString(econd_) << "\n";
    return hist.str();
  }

  /** \brief Print iterate status.

      This function prints the iteration status.

      @param[in]     algo_state    is the current state of the algorithm
      @param[in]     printHeader   if set to true will print the header at each iteration
  */
  std::string print( AlgorithmState<Real> & algo_state, bool print_header = false ) const  {
    const ROL::Ptr<const StepState<Real> > step_state = Step<Real>::getStepState();
    std::string desc = desc_->print(algo_state,false);
    desc.erase(std::remove(desc.end()-3,desc.end(),'\n'), desc.end());
    std::string name = desc_->printName();
    size_t pos = desc.find(name);
    if ( pos != std::string::npos ) {
      desc.erase(pos, name.length());
    }

    std::stringstream hist;
    if ( algo_state.iter == 0 ) {
      hist << printName();
    }
    if ( print_header ) {
      hist << printHeader();
    }
    hist << desc;
    if ( algo_state.iter == 0 ) {
      hist << "\n";
    }
    else {
      hist << std::setw(10) << std::left << step_state->nfval;              
      hist << std::setw(10) << std::left << step_state->ngrad;
      hist << "\n";
    }
    return hist.str();
  }
}; // class LineSearchStep

} // namespace ROL

#endif
