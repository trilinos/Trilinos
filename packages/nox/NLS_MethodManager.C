//-----------------------------------------------------------------------------
// ASCI_NLS_NonlinearSolver
//-----------------------------------------------------------------------------

#include "NLS_MethodManager.H"

NLS_MethodManager::NLS_MethodManager(NLS_Group& i, NLS_ParameterList& p) :
  method(NULL)
{

  string method = ParameterList.getParameter("Nonlinear Solver", "Newton");

  cout << "Nonlinear Solver: " << method << endl; 
  
  if (method == "Newton") {
    ptr = new NLS_Newton(i, p);
  } 
  else {
    cout << "ERROR: invalid choice for nonlinear solver "
	 << "in NLS_MethodManager constructor" << endl;
    throw 1;
  }

}

NLS_MethodManager::~NLS_MethodManager()
{
  delete ptr;
}


void NLS_MethodManager::resetParameters(NLS_ParameterList& p)
{
  ptr->resetParameters(p);
}

bool NLS_MethodManager::isConverged()
{
  return ptr->isConverged();
}

int NLS_MethodManager::iterate()
{
  return ptr->iterate();
}

int NLS_MethodManager::solve()
{
  return ptr->solve();
}

NLS_Group& NLS_MethodManager::getSolutionGroup() const
{
  return ptr->getSolutionGroup();
}

bool NLS_MethodManager::getProfileInfo(string& name, NLS_Parameter& p) const
{
  return ptr->getProfileInfo(name, p);
}




