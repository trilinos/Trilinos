
#ifndef NLS_METHOD_MANAGER_H
#define NLS_METHOD_MANAGER_H

// Forward Declarations
class NLS_ParameterList;
class NLS_Method;

//! Concrete class to control the instatiation of the derived 
//! NLS_Method class.  
/*! 
  This is the class that takes control from the user code and hands it over
  to the nonlinear solver component. 
*/

//-----------------------------------------------------------------------------
class NLS_MethodManager
{
public:
  //! NLS_MethodManager constructor.  User must supply a NLS_ParameterList
  //! object.  The method manager then instantiates the correct 
  //! concrete class and returns a NLS_Method object. 
  NLS_MethodManager(NLS_ParameterList& params);

  //! Destructor
  ~NLS_MethodManager();

  //! Allocate memory for the NLS_Method object.
  void createMethod();

  //! Call one step of the nonlinear solver method
  int iterate();
  
  //! Solve the nonlinear problem
  int solve();

private:
  //!Reference to the ParameterList object
  NLS_ParameterList &ParameterList;

  //! Pointer to the Method object
  NLS_Method *Method;
};

#endif
