//-----------------------------------------------------------------------------
// ASCI_NLS_NonlinearSolver
//-----------------------------------------------------------------------------
#include <string>
//#include "NLS_Interface.H"
#include "NLS_ParameterList.H"
#include "NLS_Method.H"
#include "NLS_MethodManager.H"
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Constructors/Destructors
//-----------------------------------------------------------------------------
NLS_MethodManager::NLS_MethodManager(NLS_ParameterList& params)
  :  ParameterList(params)
{ Method=0; }

NLS_MethodManager::~NLS_MethodManager()
{ if (Method != 0) delete Method; }

//-----------------------------------------------------------------------------
// Create the Global Nonlinear Solver Method
//-----------------------------------------------------------------------------

void NLS_MethodManager::createMethod()
{
  string SearchMethod = ParameterList.getParameter("Search Method", "default");
  
  if (SearchMethod=="Newton") {
    //Method = new NLS_Newton(InitialGuess, ParameterList );
    cout << "Global Method: Newton" << endl; 

  } else if (SearchMethod=="LineSearchNewton") {
    //Method = new NLS_Method(InitialGuess, ParameterList );
    cout << "Global Method: LineSearchNewton" << endl; 

  } else if (SearchMethod=="default") {
    //Method = new NLS_Newton(InitialGuess, ParameterList );
    cout << "WARNING: In NLS_MethodManager - no global method found!" << endl;
    cout << "         Defaulting to Full Newton Method" << endl;
  } else {
    cout << "ERROR: In NLS_MethodManager - no global method found!" << endl; 
    exit (-1);
  }
}
//-----------------------------------------------------------------------------
// Solve one nonlinear iteration of the Nonlinear System
//-----------------------------------------------------------------------------
int NLS_MethodManager::iterate()
{
  int status=0;
  status = Method->iterate();
  return status;
}
//-----------------------------------------------------------------------------
// Solve the Nonlinear System
//-----------------------------------------------------------------------------
int NLS_MethodManager::solve()
{
  int status=0;
  status = Method->solve();
  return status;
}
//-----------------------------------------------------------------------------




