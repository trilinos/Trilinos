#ifndef _TSF_RDP_ITERATIVESOLVER_H_
#define _TSF_RDP_ITERATIVESOLVER_H_

//! TSF_RDP_IterativeSolver:  The Trilinos Real Double Precision Iterative Solver Class.
/*! The TSF_RDP_IterativeSolver class is a pure virtual class that specifies
    the required interfaces that any Trilinos-compliant real double precision iterative
    solver must implement.
    This class is derives from the TSF_RDP_Solver class.
*/
#include "TSF_RDP_Solver.h"

class TSF_RDP_IterativeSolver: public virtual TSF_RDP_Solver {
    
  public:
  //! TSF_RDP_IterativeSolver Destructor.
  virtual ~TSF_RDP_IterativeSolver(void);

  //! Set Problem \f$Ax = b\f$.
  int SetProblem(TSF_RDP_Problem & Problem) = 0;

  //! Get Problem \f$Ax = b\f$.
  int GetProblem(TSF_RDP_Problem * & Problem) = 0;

  //!Iterate on the given TSF_RDP_Problem.
  /*! Performs up to MaxIters iteration to find X such that AX = B as defined by the TSF_RDP_Problem object.
      Uses the Converged() method from TSF_RDP_Problem object to test convergence.
    \param In
           MaxIters - Maximum number of iterations that will be performed.
    \param Out
           NumIters - Number of iterations actually performed.
    \return Returns zero if convergence reached.
  */
  virtual int Iterate(int MaxIters, int & NumIters) = 0;

};

#endif /* _TSF_RDP_ITERATIVESOLVER_H_ */
