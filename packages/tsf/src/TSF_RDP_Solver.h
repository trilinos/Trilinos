#ifndef _TSF_RDP_SOLVER_H_
#define _TSF_RDP_SOLVER_H_

//! TSF_RDP_Solver:  The Trilinos Real Double Precision Solver Class.
/*! The TSF_RDP_Solver class is a pure virtual class that specifies
    the required interfaces that any Trilinos-compliant real double precision 
    solver must implement.
    This class is most basic of all Trilinos RDP solver classes.
*/
#include "TSF_RDP_Operator.h"

class TSF_RDP_Solver: public virtual TSF_RDP_Operator {
    
  public:
  //! TSF_RDP_Solver Destructor.
  virtual ~TSF_RDPSolver(void);


  //! Performs any required setup that is needed by Apply(), must be called after SetParameters().
  virtual int SolverSetup(void) = 0;

  //! Finds solution for a given TSF_RDP_Problem.
  /*! Solves for X such that AX = B as defined by the TSF_RDP_Problem object.
    \param InOut
           Problem - A TSF_Problem object containing an operator, set of
	   right hand sides and corresponding solution vectors.
  */
  virtual int Solve(TSF_RDP_Problem & Problem) = 0;

};

#endif /* _TSF_RDP_SOLVER_H_ */
