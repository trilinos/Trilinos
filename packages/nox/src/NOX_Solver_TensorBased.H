//@HEADER
//@HEADER

#ifndef NOX_SOLVER_TENSORBASED_H
#define NOX_SOLVER_TENSORBASED_H

#include "NOX_Solver_Generic.H"	    // base class
#include "NOX_LineSearch_Manager.H" // class data element
#include "NOX_Direction_Manager.H"  // class data element
#include "NOX_Parameter_List.H"	    // class data element
#include "NOX_Abstract_Vector.H"    // for NOX::Abstract::Vector::NormType

namespace NOX {
namespace Solver {

/*!
  \brief Basic line search based solver
  
  The following parameters are valid for this solver:
  [bwb: These need to be corrected]

  
  - "Line Search" - Sublist of the line search parameters, passed to
    the NOX::LineSearch::Manager constructor. Defaults to an empty list.

  - "%Direction" - Sublist of the direction parameters, passed to
    the NOX::Direction::Manager constructor. Defaults to an empty list.

  - "Linear %Solver" - Sublist of the linear solver paramaters, passed
    to NOX::Abstract::Group::computeNewton(). Furthermore, the "Tolerance"
    within this list may be modified by the
    resetForcingTerm(). Defaults to an empty list.

    The following parameters can be specified for this line search:

<ul>
<li> "Default Step" - starting step length (defaults to 1.0)
<li> "Minimum Step" - minimum acceptable step length (defaults to 1.0e-12)
<li> "Recovery Step" - step to take when the line search fails
     (defaults to value for "Default Step")
<li> "Max Iters" - maximum number of iterations (i.e., RHS computations)
<li> "Decrease Condition" - choices are the following:
     <ul>
     <li> "Max Norm" - simple decrease using the infinity norm of the RHS (default)
     <li> "Two Norm" - simple decrease using the Euclidean norm of the RHS
     </ul>
</ul>

    \author Tammy Kolda (SNL 8950), Roger Pawlowski (SNL 9233)
*/

class TensorBased : public Generic {

public:

  //! %TensorBased-method constructor. 
  TensorBased(Abstract::Group& xgrp, StatusTest::Generic& t, const Parameter::List& p);

  //! Destructor
  virtual ~TensorBased();

  virtual bool reset(Abstract::Group& xgrp, StatusTest::Generic& t, const Parameter::List& p);
  virtual StatusTest::StatusType getStatus();
  virtual StatusTest::StatusType iterate();
  virtual StatusTest::StatusType solve();
  virtual const Abstract::Group& getSolutionGroup() const;
  virtual const Abstract::Group& getPreviousSolutionGroup() const;
  virtual int getNumIterations() const;
  virtual const Parameter::List& getParameterList() const;

protected:
  
  //! Print out initialization information and calcuation the RHS.
  virtual void init();

  //! Prints the current iteration information.
  virtual void printUpdate();

protected:
  
  //! Current solution.
  Abstract::Group* solnptr;		

  //! Previous solution pointer. 
  /*! We have both a pointer and a reference because we need to create
    a DERIVED object and then want to have a reference to it. */
  Abstract::Group* oldsolnptr;	
  //! Previous solution reference.
  Abstract::Group& oldsoln;	

  //! Current search direction.pointer.
  /*! We have both a pointer and a reference because we need to create
    a DERIVED object and then want to have a reference to it. */
  Abstract::Vector* dirptr;
  //! Current search direction.reference.
  Abstract::Vector& dir;

  //! Stopping test.
  StatusTest::Generic* testptr;		

  //! Input parameters.
  Parameter::List iparams;	

  //! Output paramters.
  mutable Parameter::List oparams; 

  //! %Search %Direction. 
  Direction::Manager direction; 

  //! Current step.
  double step;			

  //! Number of nonlinear iterations.
  int niter;                    

  //! %Status of nonlinear solver.
  StatusTest::StatusType status;

  //! "alpha" parameter in Armijo-Goldstein sufficient decrease condition.
  double alpha;
  
private:

  /*! \brief Compute the correct norm (according to the input
      parameters) of the given grp. 
  */
  double getNormF(const Abstract::Group& grp) const;

private:  

  //! Minimum step length (i.e., when we give up)
  double minStep;
 
  //! Default step
  double defaultStep;

  //! Default step for linesearch failure
  double recoveryStep;

  //! Maximum number of iterations
  int maxiters;
  
  //! Norm used in simple decrease condition
  NOX::Abstract::Vector::NormType normtype;

};
} // namespace Solver
} // namespace NOX

#endif

