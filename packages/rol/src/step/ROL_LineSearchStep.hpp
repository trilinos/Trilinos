//@HEADER
// ***********************************************************************
//
//                     Rapid Optimization Library
//
// Questions? Contact:    Drew Kouri (dpkouri@sandia.gov)
//                      Denis Ridzal (dridzal@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef ROL_LINESEARCHSTEP_H
#define ROL_LINESEARCHSTEP_H

#include "ROL_Types.hpp"
#include "ROL_Step.hpp"
#include "ROL_Secant.hpp"
#include "ROL_Krylov.hpp"
#include "ROL_NonlinearCG.hpp"
#include "ROL_LineSearch.hpp"
#include <sstream>
#include <iomanip>

/** \class ROL::LineSearchStep
    \brief Provides the interface to compute optimization steps
           with line search.
*/


namespace ROL {

template <class Real>
class LineSearchStep : public Step<Real> {
private:

  Teuchos::RCP<Secant<Real> >      secant_;
  Teuchos::RCP<Krylov<Real> >      krylov_;
  Teuchos::RCP<NonlinearCG<Real> > nlcg_;
  Teuchos::RCP<LineSearch<Real> >  lineSearch_;

  int iterKrylov_;
  int flagKrylov_;

  ELineSearch         els_; 
  ECurvatureCondition econd_;
  EDescent            edesc_;
  ESecant             esec_;
 
  int ls_nfval_;
  int ls_ngrad_;

  std::vector<bool> useInexact_;

public:

  virtual ~LineSearchStep() {}

  LineSearchStep( Teuchos::ParameterList &parlist ) {
    // Enumerations
    edesc_ = parlist.get("Descent Type",                   DESCENT_SECANT);
    els_   = parlist.get("Linesearch Type",                LINESEARCH_CUBICINTERP);
    econd_ = parlist.get("Linesearch Curvature Condition", CURVATURECONDITION_STRONGWOLFE);
    esec_  = parlist.get("Secant Type",                    SECANT_LBFGS);
    // Inexactness Information
    useInexact_.clear();
    useInexact_.push_back(parlist.get("Use Inexact Objective Function", false));
    useInexact_.push_back(parlist.get("Use Inexact Gradient", false));
    useInexact_.push_back(parlist.get("Use Inexact Hessian-Times-A-Vector", false));
     
    // Initialize Linesearch Object
    lineSearch_ = Teuchos::rcp( new LineSearch<Real>(parlist) );

    // Initialize Krylov Object
    krylov_ = Teuchos::null;
    if ( edesc_ == DESCENT_NEWTONKRYLOV || edesc_ == DESCENT_SECANTPRECOND ) {
      Real CGtol1 = parlist.get("Absolute Krylov Tolerance", 1.e-4);
      Real CGtol2 = parlist.get("Relative Krylov Tolerance", 1.e-2);
      int maxitCG = parlist.get("Maximum Number of Krylov Iterations", 20);
      krylov_ = Teuchos::rcp( new Krylov<Real>(CGtol1,CGtol2,maxitCG,useInexact_[2]) );
      iterKrylov_ = 0;
      flagKrylov_ = 0;
    }

    // Initialize Secant Object
    secant_ = Teuchos::null;
    if ( edesc_ == DESCENT_SECANT || edesc_ == DESCENT_SECANTPRECOND ) {
      int L      = parlist.get("Maximum Secant Storage",10);
      int BBtype = parlist.get("Barzilai-Borwein Type",1);
      secant_ = getSecant<Real>(esec_,L,BBtype);
    }

    // Initialize Nonlinear CG Object
    nlcg_ = Teuchos::null;
    if ( edesc_ == DESCENT_NONLINEARCG ) {
      nlcg_ = Teuchos::rcp( new NonlinearCG<Real>(NONLINEARCG_HAGAR_ZHANG) );
    }
  }

  LineSearchStep( Teuchos::RCP<Secant<Real> > &secant, Teuchos::ParameterList &parlist ) :
    secant_(secant) {
    // Enumerations
    edesc_ = parlist.get("Descent Type",                   DESCENT_SECANT);
    els_   = parlist.get("Linesearch Type",                LINESEARCH_CUBICINTERP);
    econd_ = parlist.get("Linesearch Curvature Condition", CURVATURECONDITION_STRONGWOLFE);
    esec_  = SECANT_USERDEFINED;
    // Inexactness Information
    useInexact_.clear();
    useInexact_.push_back(parlist.get("Use Inexact Objective Function", false));
    useInexact_.push_back(parlist.get("Use Inexact Gradient", false));
    useInexact_.push_back(parlist.get("Use Inexact Hessian-Times-A-Vector", false));
     
    // Initialize Linesearch Object
    lineSearch_ = Teuchos::rcp( new LineSearch<Real>(parlist) );

    // Initialize Krylov Object
    krylov_ = Teuchos::null;
    if ( edesc_ == DESCENT_NEWTONKRYLOV || edesc_ == DESCENT_SECANTPRECOND ) {
      Real CGtol1 = parlist.get("Absolute Krylov Tolerance", 1.e-4);
      Real CGtol2 = parlist.get("Relative Krylov Tolerance", 1.e-2);
      int maxitCG = parlist.get("Maximum Number of Krylov Iterations", 20);
      krylov_ = Teuchos::rcp( new Krylov<Real>(CGtol1,CGtol2,maxitCG,useInexact_) );
      iterKrylov_ = 0;
      flagKrylov_ = 0;
    }

    // Initialize Nonlinear CG Object
    nlcg_ = Teuchos::null;
    if ( edesc_ == DESCENT_NONLINEARCG ) {
      nlcg_ = Teuchos::rcp( new NonlinearCG<Real>(NONLINEARCG_HAGAR_ZHANG) );
    }
  }

  /** \brief Compute step.
  */
  void compute( Vector<Real> &s, const Vector<Real> &x, Objective<Real> &obj, AlgorithmState<Real> &algo_state ) {
    // Compute step s
    if ( this->edesc_ == DESCENT_NEWTONKRYLOV || this->edesc_ == DESCENT_SECANTPRECOND ) {
      this->flagKrylov_ = 0;
      this->krylov_->CG(s,this->iterKrylov_,this->flagKrylov_,*(Step<Real>::state_->gradientVec),x,obj,this->secant_);
    }
    else if ( this->edesc_ == DESCENT_NEWTON ) {
      Real tol = std::sqrt(ROL_EPSILON);
      obj.invHessVec(s,*(Step<Real>::state_->gradientVec),x,tol);
    }
    else if ( this->edesc_ == DESCENT_SECANT ) {
      this->secant_->applyH(s,*(Step<Real>::state_->gradientVec),x);
    }
    else if ( this->edesc_ == DESCENT_NONLINEARCG ) {
      this->nlcg_->run(s,*(Step<Real>::state_->gradientVec),x,obj);
    }

    // Check if s is a descent direction
    Real gs = -(Step<Real>::state_->gradientVec)->dot(s);
    if ( gs > 0.0 || this->flagKrylov_ == 2 || this->edesc_ == DESCENT_STEEPEST ) {
      s.set(*(Step<Real>::state_->gradientVec));
      gs = -(Step<Real>::state_->gradientVec)->dot(s);
    }
    s.scale(-1.0);

    // Perform line search
    Real alpha = 0.0;
    Real fnew  = algo_state.value;
    this->ls_nfval_ = 0;
    this->ls_ngrad_ = 0;
    this->lineSearch_->run(alpha,fnew,this->ls_nfval_,this->ls_ngrad_,gs,s,x,obj);
    algo_state.nfval += this->ls_nfval_;
    algo_state.ngrad += this->ls_ngrad_;
    s.scale(alpha);

    // Update step state information
    (Step<Real>::state_->descentVec)->set(s);

    // Update algorithm state information
    algo_state.snorm = s.norm();
    algo_state.value = fnew;
  }

  /** \brief Update step, if successful.
  */
  void update( Vector<Real> &x, const Vector<Real> &s, Objective<Real> &obj, AlgorithmState<Real> &algo_state ) {
    Real tol = std::sqrt(ROL_EPSILON);

    // Update iterate
    x.axpy(1.0, s);

    // Compute new gradient
    Teuchos::RCP<Vector<Real> > gp;
    if ( this->edesc_ == DESCENT_SECANT || this->edesc_ == DESCENT_SECANTPRECOND ) {
      gp = x.clone();
      gp->set(*(Step<Real>::state_->gradientVec));
    }
    obj.gradient(*(Step<Real>::state_->gradientVec),x,tol);
    algo_state.ngrad++;

    // Update Secant Information
    if ( this->edesc_ == DESCENT_SECANT || this->edesc_ == DESCENT_SECANTPRECOND ) {
      secant_->update(*(Step<Real>::state_->gradientVec),*gp,s,algo_state.snorm,algo_state.iter+1);
    }

    // Update algorithm state
    (algo_state.iterateVec)->set(x);
    algo_state.gnorm = (Step<Real>::state_->gradientVec)->norm();
    algo_state.iter++;
  }

  /** \brief Print iterate status.
  */
  std::string print( AlgorithmState<Real> & algo_state, bool printHeader = false ) const  {
    std::stringstream hist;
    if ( algo_state.iter == 0 ) {
      hist << "\n" << EDescentToString(this->edesc_) 
           << " with " << ELineSearchToString(this->els_) 
           << " Linesearch satisfying " 
           << ECurvatureConditionToString(this->econd_) << "\n";
      if ( this->edesc_ == DESCENT_SECANT || this->edesc_ == DESCENT_SECANTPRECOND ) {
        hist << "Secant Type: " << ESecantToString(this->esec_) << "\n";
      }
      hist << "  ";
      hist << std::setw(6) << std::left << "iter";  
      hist << std::setw(15) << std::left << "value";
      hist << std::setw(15) << std::left << "gnorm"; 
      hist << std::setw(15) << std::left << "snorm";
      hist << std::setw(10) << std::left << "#fval";
      hist << std::setw(10) << std::left << "#grad";
      hist << std::setw(10) << std::left << "ls_#fval";
      hist << std::setw(10) << std::left << "ls_#grad";
      if (    this->edesc_ == DESCENT_NEWTONKRYLOV 
           || this->edesc_ == DESCENT_SECANTPRECOND ) {
        hist << std::setw(10) << std::left << "iterCG";
        hist << std::setw(10) << std::left << "flagCG";
      }
      hist << "\n";
      hist << std::setfill('-') << std::setw(110) << "-" << "\n";
      hist << std::setfill(' ') << "  ";
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
      hist << std::setw(10) << std::left << this->ls_nfval_;              
      hist << std::setw(10) << std::left << this->ls_ngrad_;              
      if (    this->edesc_ == DESCENT_NEWTONKRYLOV 
           || this->edesc_ == DESCENT_SECANTPRECOND ) {
        hist << std::setw(10) << std::left << this->iterKrylov_;
        hist << std::setw(10) << std::left << this->flagKrylov_;
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
