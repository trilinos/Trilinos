
// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#include <cmath>
#include <cstring>
#include "NLS_Vector.H"
#include "NLS_Utilities.H"

#include "NLS_ForcingTerm.H"

NLS_ForcingTerm::NLS_ForcingTerm() 
{
  eta_max = 0.01;
}

NLS_ForcingTerm::~NLS_ForcingTerm() 
{
  
}

double NLS_ForcingTerm::getForcingTerm(NLS_ParameterList& p)
{
  eta_min = p.getParameter("Minimum Linear Solver Tolerance", 1.0e-6);
  eta_km1 = p.getParameter("Linear Solver Tolerance", 0.0);

  string Method = p.getParameter("Forcing Term Method","Failed!");

  if (Method == "Constant") {    

    eta_k = eta_min;

    if (NLS_Utilities::doPrint(3)) 
      cout << "      Forcing Term: Constant = " << eta_k << endl << endl;

    return eta_k;
  }        

  else if (Method == "Type 1") {
    
    residual_norm_km1 = p.getParameter("Residual Norm km1",0.0);
    residual_norm_k = p.getParameter("Residual Norm k",0.0);
    residual_norm_linear_model = p.getParameter("Linearized Residual Norm",0.0);
    
    eta_k = fabs(residual_norm_k -
		 residual_norm_linear_model)/residual_norm_km1;
     
    if (NLS_Utilities::doPrint(4)) {
      cout << "      Forcing Term Type 1" << endl;
      cout << "      Residual Norm k-1 =            " 
	   << residual_norm_km1 << endl;
      cout << "      Residual Norm Linear Model k = " 
	   << residual_norm_linear_model << endl;
      cout << "      Residual Norm k =              " 
	   << residual_norm_k << endl;
      cout << "      Calculated eta_k =             " << eta_k << endl;
    }
  
    // Impose safeguard and constraints ...
    alpha = (1.0 + sqrt(5.0))/2.0;
    eta_k_alpha = pow(eta_km1,alpha);
    eta_k = max(eta_k, eta_k_alpha);
    eta_k = max(eta_k, eta_min);
    eta_k = min(eta_max, eta_k);
    
    if (NLS_Utilities::doPrint(4)) {
      cout << "      Applying bounds...             " << endl;
      cout << "      Final eta_k =                  " << eta_k << endl;
    }
  
    if (NLS_Utilities::doPrint(3)) 
      cout << "      Forcing Term: Type 1 = "<< eta_k << endl << endl;
     
    return eta_k;
  }
    
  else if (Method == "Type 2") {  
    
    alpha = p.getParameter("Forcing Term alpha", 1.5);
    gamma = p.getParameter("Forcing Term gamma", 0.9);
    residual_norm_km1 = p.getParameter("Residual Norm km1",0.0);
    residual_norm_k = p.getParameter("Residual Norm k",0.0);
    residual_ratio = residual_norm_k/residual_norm_km1;
    
    eta_k = gamma * pow(residual_ratio, alpha);
     
    if (NLS_Utilities::doPrint(4)) {
      cout << "      Forcing Term Type 2" << endl;
      cout << "      Residual Norm k-1 =            " 
	   << residual_norm_km1 << endl;
      cout << "      Residual Norm k =              " 
	   << residual_norm_k << endl;
      cout << "      Calculated eta_k =             " << eta_k << endl;
    }
  
    // Impose safeguard and constraints ... 
    eta_k_alpha = gamma * pow(eta_km1,alpha);
    eta_k = max(eta_k, eta_k_alpha);
    eta_k = max(eta_k, eta_min);
    eta_k = min(eta_max, eta_k);
    
    if (NLS_Utilities::doPrint(4)) {
      cout << "      Applying bounds...             " << endl;
      cout << "      Final eta_k =                  " << eta_k << endl;
    }
  
    if (NLS_Utilities::doPrint(3)) 
      cout << "      Forcing Term: Type 2 = "<< eta_k << endl << endl;
      
    return eta_k;
  }

  // Default 

  eta_k = eta_min;

  if (NLS_Utilities::doPrint(3)) 
    cout << "Warning: Forcing Term Method was not set."  << endl 
	 << "Returning default linear solver tolerance of " 
	 << eta_k << "." << endl << endl;

  return eta_k;
}
