// $Id$ 
// $Source$ 

// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#include "NLS_Newton.H"
#include "NLS_MethodManager.H"

NLS_MethodManager::NLS_MethodManager(NLS_Group& i, 
				     NLS_Group& s, 
				     NLS_ParameterList& p) :
  ptr(NULL)
{

  string method = p.getParameter("Nonlinear Solver", "Newton");

  cout << "Nonlinear Solver: " << method << endl; 
  
  if (method == "Newton") {
    ptr = new NLS_Newton(i, s, p);
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




