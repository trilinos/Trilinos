#ifndef _NLS_METHOD_H
#define _NLS_METHOD_H

#include <string>
#include "NLS_Parameter.H"
#include "NLS_ParameterList.H"
#include "NLS_Group.H"

//! Abstract nonlinear solver method interface.
class NLS_Method 
{

protected: 
  //! Constructor 
  NLS_Method();

  //! This is what the basic constructor for a derived class should look like
  // NLS_Method(NLS_Function& f, NLS_ParameterList& p) {};
  
public:

  //! Deconstructor
  virtual ~NLS_Method();

  //! Change parameters for the nonlinear method.
  virtual void resetParameters(NLS_ParameterList& p) = 0;

  //! Check for convergence of the nonlinear problem
  virtual bool isConverged() = 0;

  //! Do one nonlinear method iteration
  virtual int iterate() = 0;

  //! Solve the nonlinear problem
  virtual int solve() = 0;

  //! Return a reference to the current solution group
  virtual NLS_Group& getSolutionGroup() const = 0;

  //! Fill in the parameter with the named profile information
  virtual bool getProfileInfo(string& name, NLS_Parameter& p) const = 0;

};



#endif
