#ifndef NLS_METHOD_MANAGER_H
#define NLS_METHOD_MANAGER_H

#include "NLS_Method.H"		// Base class

//! Concrete class to control the instatiation of the derived 
//! NLS_Method class.  
/*! 
  This is the class that takes control from the user code and hands it over
  to the nonlinear solver component. 
*/

//-----------------------------------------------------------------------------
class NLS_MethodManager : public NLS_Method {

public:
  //! Generic interface to any of the various methods. 
  NLS_MethodManager(NLS_Group& initial_guess, NLS_ParameterList& p);

  //! Destructor
  ~NLS_MethodManager();


  void resetParameters(NLS_ParameterList& p);
  bool isConverged();
  int iterate();
  int solve();
  NLS_Group& getSolutionGroup() const;
  bool getProfileInfo(string& name, NLS_Parameter& p) const;


private:

  //! Pointer to the solver object
  NLS_Method *ptr;
};

#endif
