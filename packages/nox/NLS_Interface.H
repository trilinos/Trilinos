#ifndef _NLS_INTERFACE_H
#define _NLS_INTERFACE_H

//! NLS_INTERFACE: Nonlinear Solver Function Interface
/*! The  interface between the nonlinear solver and the
  associated function. We wish to solve F(x) = 0.  */

#include "NLS_Vector.H"

class NLS_Interface {

protected:

  //! NLS_Interface Constructor.
  NLS_Interface();

public:

  //! NLS_Interface Deconstructor
  virtual ~NLS_Interface() = 0;

  /** @name Creating stuff */
  //@{
  virtual NLS_Vector& createVectorByDeepCopy(NLS_Vector& copyFrom) = 0;
  virtual NLS_Group& createGroupByDeepCopy(NLS_Group& copyFrom) = 0;
  //@}

  //! Return initial guess
  virtual NLS_Group& getInitialGuess() = 0;

};


#endif
