// $Id$
// $Source$

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
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
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#include "Teuchos_LAPACK.hpp"
#include "Teuchos_ParameterList.hpp"
#include "LOCA_GlobalData.H"
#include "LOCA_EigenvalueSort_Strategies.H"

NOX::Abstract::Group::ReturnType 
LOCA::EigenvalueSort::LargestMagnitude::sort(int n, double* evals, 
					     std::vector<int>* perm) const
{
  int i, j, tempord;
  double temp, temp2;
  Teuchos::LAPACK<int,double> lapack;

  // Reset the index	
  if (perm) {
    for (i=0; i < n; i++) {
      (*perm)[i] = i;
    }
  }

  //---------------------------------------------------------------
  // Sort eigenvalues in decreasing order of magnitude
  //---------------------------------------------------------------
  for (j=1; j < n; ++j) {
    temp = evals[j]; 
    if (perm)
      tempord = (*perm)[j];
    temp2 = evals[j]*evals[j];
    for (i=j-1; i>=0 && (evals[i]*evals[i])<temp2; --i) {
      evals[i+1]=evals[i];
      if (perm)
	(*perm)[i+1]=(*perm)[i];
    }
    evals[i+1] = temp; 
    if (perm)
      (*perm)[i+1] = tempord;	
  }
  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType 
LOCA::EigenvalueSort::LargestMagnitude::sort(int n, double* r_evals, 
					     double* i_evals,
					     std::vector<int>* perm) const
{
  int i, j, tempord;
  double temp, tempr, tempi;
  Teuchos::LAPACK<int,double> lapack;

  //
  // Reset the index	
  if (perm) {
    for (i=0; i < n; i++) {
      (*perm)[i] = i;
    }
  }

  //---------------------------------------------------------------
  // Sort eigenvalues in decreasing order of magnitude
  //---------------------------------------------------------------
  for (j=1; j < n; ++j) {
    tempr = r_evals[j]; tempi = i_evals[j]; 
    if (perm)
      tempord = (*perm)[j];
    temp=lapack.LAPY2(r_evals[j],i_evals[j]);
    for (i=j-1; i>=0 && lapack.LAPY2(r_evals[i],i_evals[i])<temp; --i) {
      r_evals[i+1]=r_evals[i]; i_evals[i+1]=i_evals[i];
      if (perm)
	(*perm)[i+1]=(*perm)[i];
    }
    r_evals[i+1] = tempr; i_evals[i+1] = tempi; 
    if (perm)
      (*perm)[i+1] = tempord;	
  }	
  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType 
LOCA::EigenvalueSort::SmallestMagnitude::sort(int n, double* evals, 
					     std::vector<int>* perm) const
{
  int i, j, tempord;
  double temp, temp2;
  Teuchos::LAPACK<int,double> lapack;

  // Reset the index	
  if (perm) {
    for (i=0; i < n; i++) {
      (*perm)[i] = i;
    }
  }

  //---------------------------------------------------------------
  // Sort eigenvalues in increasing order of magnitude
  //---------------------------------------------------------------
  for (j=1; j < n; ++j) {
    temp = evals[j]; 
    if (perm)
      tempord = (*perm)[j];
    temp2 = evals[j]*evals[j];
    for (i=j-1; i>=0 && (evals[i]*evals[i])>temp2; --i) {
      evals[i+1]=evals[i];
      if (perm)
	(*perm)[i+1]=(*perm)[i];
    }
    evals[i+1] = temp; 
    if (perm) 
      (*perm)[i+1] = tempord;	
  }
  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType 
LOCA::EigenvalueSort::SmallestMagnitude::sort(int n, double* r_evals, 
					     double* i_evals,
					     std::vector<int>* perm) const
{
  int i, j, tempord;
  double temp, tempr, tempi;
  Teuchos::LAPACK<int,double> lapack;

  //
  // Reset the index	
  if (perm) {
    for (i=0; i < n; i++) {
      (*perm)[i] = i;
    }
  }

  //---------------------------------------------------------------
  // Sort eigenvalues in increasing order of magnitude
  //---------------------------------------------------------------
  for (j=1; j < n; ++j) {
    tempr = r_evals[j]; tempi = i_evals[j]; 
    if (perm)
      tempord = (*perm)[j];
    temp=lapack.LAPY2(r_evals[j],i_evals[j]);
    for (i=j-1; i>=0 && lapack.LAPY2(r_evals[i],i_evals[i])>temp; --i) {
      r_evals[i+1]=r_evals[i]; i_evals[i+1]=i_evals[i];
      if (perm)
	(*perm)[i+1]=(*perm)[i];
    }
    r_evals[i+1] = tempr; i_evals[i+1] = tempi; 
    if (perm)
      (*perm)[i+1] = tempord;	
  }	
  return NOX::Abstract::Group::Ok;
}


NOX::Abstract::Group::ReturnType 
LOCA::EigenvalueSort::LargestReal::sort(int n, double* evals, 
					std::vector<int>* perm) const
{
  int i, j, tempord;
  double temp;
  Teuchos::LAPACK<int,double> lapack;

  // Reset the index	
  if (perm) {
    for (i=0; i < n; i++) {
      (*perm)[i] = i;
    }
  }

  //---------------------------------------------------------------
  // Sort eigenvalues in decreasing order of real part
  //---------------------------------------------------------------
  for (j=1; j < n; ++j) {
    temp = evals[j]; 
    if (perm)
      tempord = (*perm)[j];
    for (i=j-1; i>=0 && evals[i]<temp; --i) {
      evals[i+1]=evals[i];
      if (perm)
	(*perm)[i+1]=(*perm)[i];
    }
    evals[i+1] = temp; 
    if (perm)
      (*perm)[i+1] = tempord;	
  }
  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType 
LOCA::EigenvalueSort::LargestReal::sort(int n, double* r_evals, 
					double* i_evals,
					std::vector<int>* perm) const
{
  int i, j, tempord;
  double tempr, tempi;
  Teuchos::LAPACK<int,double> lapack;

  // Reset the index	
  if (perm) {
    for (i=0; i < n; i++) {
      (*perm)[i] = i;
    }
  }

  //---------------------------------------------------------------
  // Sort eigenvalues in decreasing order of real part
  //---------------------------------------------------------------
  for (j=1; j < n; ++j) {
    tempr = r_evals[j]; tempi = i_evals[j]; 
    if (perm)
      tempord = (*perm)[j];
    for (i=j-1; i>=0 && r_evals[i]<tempr; --i) {
      r_evals[i+1]=r_evals[i]; i_evals[i+1]=i_evals[i];
      if (perm)
	(*perm)[i+1]=(*perm)[i];
    }
    r_evals[i+1] = tempr; i_evals[i+1] = tempi; 
    if (perm)
      (*perm)[i+1] = tempord;	
  }	
  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType 
LOCA::EigenvalueSort::SmallestReal::sort(int n, double* evals, 
					 std::vector<int>* perm) const
{
  int i, j, tempord;
  double temp;
  Teuchos::LAPACK<int,double> lapack;

  // Reset the index	
  if (perm) {
    for (i=0; i < n; i++) {
      (*perm)[i] = i;
    }
  }

  //---------------------------------------------------------------
  // Sort eigenvalues in increasing order of real part
  //---------------------------------------------------------------
  for (j=1; j < n; ++j) {
    temp = evals[j]; 
    if (perm)
      tempord = (*perm)[j];
    for (i=j-1; i>=0 && evals[i]>temp; --i) {
      evals[i+1]=evals[i];
      if (perm)
	(*perm)[i+1]=(*perm)[i];
    }
    evals[i+1] = temp; 
    if (perm)
      (*perm)[i+1] = tempord;	
  }
  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType 
LOCA::EigenvalueSort::SmallestReal::sort(int n, double* r_evals, 
					double* i_evals,
					std::vector<int>* perm) const
{
  int i, j, tempord;
  double tempr, tempi;
  Teuchos::LAPACK<int,double> lapack;

  // Reset the index	
  if (perm) {
    for (i=0; i < n; i++) {
      (*perm)[i] = i;
    }
  }

  //---------------------------------------------------------------
  // Sort eigenvalues in increasing order of real part
  //---------------------------------------------------------------
  for (j=1; j < n; ++j) {
    tempr = r_evals[j]; tempi = i_evals[j]; 
    if (perm)
      tempord = (*perm)[j];
    for (i=j-1; i>=0 && r_evals[i]>tempr; --i) {
      r_evals[i+1]=r_evals[i]; i_evals[i+1]=i_evals[i];
      if (perm)
	(*perm)[i+1]=(*perm)[i];
    }
    r_evals[i+1] = tempr; i_evals[i+1] = tempi; 
    if (perm)
      (*perm)[i+1] = tempord;	
  }	
  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType 
LOCA::EigenvalueSort::LargestImaginary::sort(int n, double* evals, 
					     std::vector<int>* perm) const
{
  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType 
LOCA::EigenvalueSort::LargestImaginary::sort(int n, double* r_evals, 
					     double* i_evals,
					     std::vector<int>* perm) const
{
  int i, j, tempord;
  double tempr, tempi;
  Teuchos::LAPACK<int,double> lapack;

  // Reset the index	
  if (perm) {
    for (i=0; i < n; i++) {
      (*perm)[i] = i;
    }
  }

  //---------------------------------------------------------------
  // Sort eigenvalues in decreasing order of imaginary part
  //---------------------------------------------------------------
  for (j=1; j < n; ++j) {
    tempr = r_evals[j]; tempi = i_evals[j]; 
    if (perm)
      tempord = (*perm)[j];
    for (i=j-1; i>=0 && i_evals[i]<tempi; --i) {
      r_evals[i+1]=r_evals[i]; i_evals[i+1]=i_evals[i];
      if (perm)
	(*perm)[i+1]=(*perm)[i];
    }
    r_evals[i+1] = tempr; i_evals[i+1] = tempi; 
    if (perm)
      (*perm)[i+1] = tempord;	
  }
  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType 
LOCA::EigenvalueSort::SmallestImaginary::sort(int n, double* evals, 
					      std::vector<int>* perm) const
{
  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType 
LOCA::EigenvalueSort::SmallestImaginary::sort(int n, double* r_evals, 
					     double* i_evals,
					     std::vector<int>* perm) const
{
  int i, j, tempord;
  double tempr, tempi;
  Teuchos::LAPACK<int,double> lapack;

  // Reset the index	
  if (perm) {
    for (i=0; i < n; i++) {
      (*perm)[i] = i;
    }
  }

  //---------------------------------------------------------------
  // Sort eigenvalues in increasing order of imaginary part
  //---------------------------------------------------------------
  for (j=1; j < n; ++j) {
    tempr = r_evals[j]; tempi = i_evals[j]; 
    if (perm)
      tempord = (*perm)[j];
    for (i=j-1; i>=0 && i_evals[i]>tempi; --i) {
      r_evals[i+1]=r_evals[i]; i_evals[i+1]=i_evals[i];
      if (perm)
	(*perm)[i+1]=(*perm)[i];
    }
    r_evals[i+1] = tempr; i_evals[i+1] = tempi; 
    if (perm)
      (*perm)[i+1] = tempord;	
  }
  return NOX::Abstract::Group::Ok;
}

LOCA::EigenvalueSort::LargestRealInverseCayley::LargestRealInverseCayley(
	     const Teuchos::RCP<LOCA::GlobalData>& global_data,
	     const Teuchos::RCP<Teuchos::ParameterList>& eigenParams) :
  sigma(0.0),
  mu(0.0)
{
  sigma = eigenParams->get("Cayley Pole",0.0);
  mu = eigenParams->get("Cayley Zero",0.0);
}

NOX::Abstract::Group::ReturnType 
LOCA::EigenvalueSort::LargestRealInverseCayley::sort(
						int n, double* evals, 
						std::vector<int>* perm) const
{
  int i, j, tempord;
  double temp, templambda;
  Teuchos::LAPACK<int,double> lapack;

  // Reset the index	
  if (perm) {
    for (i=0; i < n; i++) {
      (*perm)[i] = i;
    }
  }

  //---------------------------------------------------------------
  // Sort eigenvalues according to the real part of the inverse-Cayley
  // transformation
  //---------------------------------------------------------------
  for (j=1; j < n; ++j) {
    temp = evals[j]; 
    tempord = (*perm)[j];
    templambda=realLambda(evals[j],0);
    for (i=j-1; i>=0 && realLambda(evals[i],0)<templambda; --i) {
      evals[i+1]=evals[i]; 
      (*perm)[i+1]=(*perm)[i];
    }
    evals[i+1] = temp; (*perm)[i+1] = tempord;	
  }	
  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType 
LOCA::EigenvalueSort::LargestRealInverseCayley::sort(
					        int n, double* r_evals, 
						double* i_evals,
						std::vector<int>* perm) const
{
  int i, j, tempord;
  double temp, tempr, tempi;
  Teuchos::LAPACK<int,double> lapack;

  // Reset the index	
  if (perm) {
    for (i=0; i < n; i++) {
      (*perm)[i] = i;
    }
  }

  //---------------------------------------------------------------
  // Sort eigenvalues according to the real part of the inverse-Cayley
  // transformation
  //---------------------------------------------------------------
  for (j=1; j < n; ++j) {
    tempr = r_evals[j]; tempi = i_evals[j]; 
    tempord = (*perm)[j];
    temp=realLambda(r_evals[j],i_evals[j]);
    for (i=j-1; i>=0 && realLambda(r_evals[i],i_evals[i])<temp; --i) {
      r_evals[i+1]=r_evals[i]; i_evals[i+1]=i_evals[i];
      (*perm)[i+1]=(*perm)[i];
    }
    r_evals[i+1] = tempr; i_evals[i+1] = tempi; (*perm)[i+1] = tempord;	
  }
  return NOX::Abstract::Group::Ok;
}

double
LOCA::EigenvalueSort::LargestRealInverseCayley::realLambda(double er, 
							   double ei) const
{
  // Reject if it is to the right of sigma --- these are junk
  double reLambda = (sigma*(er*er+ei*ei) - (sigma+mu)*er + mu) / 
    ( (er-1.0)*(er-1.0) + ei*ei);
  if (reLambda > sigma) 
    return -1.0e6;
  else                   
    return reLambda;
}
