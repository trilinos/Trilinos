
// Nonlinear Solver Package (NLSPACK)
// COPYRIGHT (2002) Sandia Corporation.
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// LICENSE & WARRANTY INFORMATION in README.txt and LICENSE.txt.
// CONTACT T. Kolda (tgkolda@sandia.gov) or R. Pawlowski (rppawlo@sandia.gov)

#include "NLS_AdaptiveTolerance.H" // struct definition

#include <cmath>		// for min, max
#include "NLS_Utilities.H"	// for static function doPrint

double NLS_AdaptiveTolerance::operator()(NLS_ParameterList& p, double normrhs, 
				   double normoldrhs, double normpredrhs)
{
  const double eta_min = p.getParameter("Minimum Linear Solver Tolerance", 1.0e-6);
  const double eta_km1 = p.getParameter("Linear Solver Tolerance", 0.0);
  const double eta_max = 0.01;

  const string Method = p.getParameter("Forcing Term Method", "Failed!");

  if (Method == "Constant") {    

    if (NLS_Utilities::doPrint(3)) 
      cout << "      Forcing Term: Constant = " << eta_min << "\n" << endl;

    return eta_min;
  }        

  else if (Method == "Type 1") {
    
    double eta_k = fabs(normrhs - normpredrhs) / normoldrhs;
     
    if (NLS_Utilities::doPrint(4)) {
      cout << "      Forcing Term Type 1" << "\n";
      cout << "      Residual Norm k-1 =            " 
	   << normoldrhs << "\n";
      cout << "      Residual Norm Linear Model k = " 
	   << normpredrhs << "\n";
      cout << "      Residual Norm k =              " 
	   << normrhs << "\n";
      cout << "      Calculated eta_k =             " << eta_k << endl;
    }
  
    // Impose safeguard and constraints ...
    const double alpha = (1.0 + sqrt(5.0)) / 2.0;
    const double eta_k_alpha = pow(eta_km1, alpha);
    eta_k = max(eta_k, eta_k_alpha);
    eta_k = max(eta_k, eta_min);
    eta_k = min(eta_max, eta_k);
    
    if (NLS_Utilities::doPrint(4)) {
      cout << "      Applying bounds...             " << "\n";
      cout << "      Final eta_k =                  " << eta_k << endl;
    }
  
    if (NLS_Utilities::doPrint(3)) 
      cout << "      Forcing Term: Type 1 = "<< eta_k << "\n" << endl;
     
    return eta_k;
  }
    
  else if (Method == "Type 2") {  
    
    const double alpha = p.getParameter("Forcing Term alpha", 1.5);
    const double gamma = p.getParameter("Forcing Term gamma", 0.9);
    const double residual_ratio = normrhs / normoldrhs;
    
    double eta_k = gamma * pow(residual_ratio, alpha);
     
    if (NLS_Utilities::doPrint(4)) {
      cout << "      Forcing Term Type 2" << "\n";
      cout << "      Residual Norm k-1 =            " 
	   << normoldrhs << "\n";
      cout << "      Residual Norm k =              " 
	   << normrhs << "\n";
      cout << "      Calculated eta_k =             " << eta_k << endl;
    }
  
    // Impose safeguard and constraints ... 
    const double eta_k_alpha = gamma * pow(eta_km1, alpha);
    eta_k = max(eta_k, eta_k_alpha);
    eta_k = max(eta_k, eta_min);
    eta_k = min(eta_max, eta_k);
    
    if (NLS_Utilities::doPrint(4)) {
      cout << "      Applying bounds...             " << "\n";
      cout << "      Final eta_k =                  " << eta_k << endl;
    }
  
    if (NLS_Utilities::doPrint(3)) 
      cout << "      Forcing Term: Type 2 = "<< eta_k << "\n" << endl;
      
    return eta_k;
  }

  // Default 

  if (NLS_Utilities::doPrint(3)) 
    cout << "Warning: Forcing Term Method was not set."  << "\n" 
	 << "Returning default linear solver tolerance of " 
	 << eta_min << "." << "\n" << endl;

  return eta_min;
}
