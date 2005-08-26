//
// @HEADER
// ***********************************************************************
// 
//                           Capo Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER


#ifndef Capo_EA_INTEGRATOR_H
#define Capo_EA_INTEGRATOR_H

/************************************************************ 
File:      EA_integrator.hpp
Purpose:   The Elezgaray & Arneodo model problem.
Date:      8-1-05
Author:    Joseph Simonis
**************************************************************/


/**** Include Files ****/


/**** EAIntegrator Class ****/
class EAIntegrator
{
public:
  // Constructor
  EAIntegrator(int n);
  
  // Destructor
  ~EAIntegrator();
  
  // Integration function
  bool Integrate(double *u_T, const double *u, double t, double lambda);
  
private:
  int numElements_;

  //! Tridiagonal system solver
  bool thomas(double* a, double* b, double* c, double* r, int n);
  
  //! Runga-Kutta integrator
  bool ea_rungakutta(double T, double* u, int n);

};

#endif //Capo_EA_INTEGRATOR_H
