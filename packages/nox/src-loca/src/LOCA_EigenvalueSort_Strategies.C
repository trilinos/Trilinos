// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_LAPACK.hpp"
#include "Teuchos_ParameterList.hpp"
#include "LOCA_GlobalData.H"
#include "LOCA_EigenvalueSort_Strategies.H"

NOX::Abstract::Group::ReturnType
LOCA::EigenvalueSort::LargestMagnitude::sort(int n, double* evals,
                         std::vector<int>* perm) const
{
  int i, j;
  int tempord = 0;
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
  int i, j;
  int tempord = 0;
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
  int i, j;
  int tempord = 0;
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
  int i, j;
  int tempord = 0;
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
  int i, j;
  int tempord = 0;
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
  int i, j;
  int tempord = 0;
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
  int i, j;
  int tempord = 0;
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
  int i, j;
  int tempord = 0;
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
LOCA::EigenvalueSort::LargestImaginary::sort(int /* n */, double* /* evals */,
                         std::vector<int>* /* perm */) const
{
  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType
LOCA::EigenvalueSort::LargestImaginary::sort(int n, double* r_evals,
                         double* i_evals,
                         std::vector<int>* perm) const
{
  int i, j;
  int tempord = 0;
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
LOCA::EigenvalueSort::SmallestImaginary::sort(int /* n */, double* /* evals */,
                          std::vector<int>* /* perm */) const
{
  return NOX::Abstract::Group::NotDefined;
}

NOX::Abstract::Group::ReturnType
LOCA::EigenvalueSort::SmallestImaginary::sort(int n, double* r_evals,
                         double* i_evals,
                         std::vector<int>* perm) const
{
  int i, j;
  int tempord = 0;
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
         const Teuchos::RCP<LOCA::GlobalData>& /* global_data */,
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
  int i, j;
  int tempord = 0;
  double temp, templambda;
  Teuchos::LAPACK<int,double> lapack;

  TEUCHOS_TEST_FOR_EXCEPT(perm == NULL);

  // Reset the index
  for (i=0; i < n; i++)
    (*perm)[i] = i;

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
  int i, j;
  int tempord = 0;
  double temp, tempr, tempi;
  Teuchos::LAPACK<int,double> lapack;

  TEUCHOS_TEST_FOR_EXCEPT(perm == NULL);

  // Reset the index
  for (i=0; i < n; i++)
    (*perm)[i] = i;

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
