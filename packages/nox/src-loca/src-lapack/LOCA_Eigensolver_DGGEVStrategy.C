// $Id$
// $Source$

//@HEADER
// ************************************************************************
//
//                  LOCA Continuation Algorithm Package
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
// Questions? Contact Andy Salinger (agsalin@sandia.gov) or Eric Phipps
// (etphipp@sandia.gov), Sandia National Laboratories.
//
// ************************************************************************
//@HEADER

#include "LOCA_Eigensolver_DGGEVStrategy.H"
#include "NOX_Parameter_List.H"
#include "NOX_LAPACK_Wrappers.H"
#include "LOCA_LAPACK_Group.H"
#include "LOCA_Utils.H"

LOCA::Eigensolver::DGGEVStrategy::DGGEVStrategy(
		 const Teuchos::RefCountPtr<NOX::Parameter::List>& eigParams,
		 const Teuchos::RefCountPtr<NOX::Parameter::List>& solParams) 
{
}

LOCA::Eigensolver::DGGEVStrategy::~DGGEVStrategy() 
{
}

NOX::Abstract::Group::ReturnType
LOCA::Eigensolver::DGGEVStrategy::computeEigenvalues(
						 NOX::Abstract::Group* group)
{

  // Get LAPACK group
  LOCA::LAPACK::Group* grp = 
    dynamic_cast<LOCA::LAPACK::Group*>(group);
 
  // Check to make sure we have dggev available if we need generalized 
  // eigenvalues.
#ifndef HAVE_LAPACK_GENEV
  if (hasMassMatrix) {
    if (LOCA::Utils::doPrint(Utils::StepperIteration)) {
      errorCheck.printWarning(
	"LOCA::Eigensolver::DGGEVStrategy::computeEigenvalues",
	"LAPACK Generalized eigensolver (dggev) requested but not available!");
    }
    return LOCA::Abstract::Group::Ok;
  }

#endif

  bool hasMassMatrix = grp->hasMass();

  // Make sure Jacobian & mass matrices are fresh
  grp->computeJacobian();
  if (hasMassMatrix) 
    grp->computeMassMatrix();

  // Get data out of group
  NOX::LAPACK::Matrix& jacobianMatrix = grp->getJacobianMatrix();
  NOX::LAPACK::Matrix& massMatrix = grp->getMassMatrix();
  
  // Size of matrix
  int n = jacobianMatrix.numRows();
  int lda = jacobianMatrix.numRows();
  int ldb = massMatrix.numRows();

  // Space to hold right eigenvectors
  double *vr = new double[n*n];

  // Space to hold real and imaginary eigenvalues
  double *alphar = new double[n];
  double *alphai = new double[n];
  double *beta = new double[n];

  // Size of work array, set to -1 to do a workspace query
  int lwork = -1;

  // Initial work "array"
  double work0;

  // Actual work array
  double *work;

  // Return code
  int info;

  // Copy Jacobian matrix since lapack routines overwrite it
  NOX::LAPACK::Matrix J(jacobianMatrix);

  NOX::LAPACK::Matrix M;

  // First do a workspace query
  if (hasMassMatrix) {

    // Copy mass matrix since lapack routines overwrite it
    M = massMatrix;

#ifdef HAVE_LAPACK_GENEV
    DGGEV_F77("N", "V", &n, &J(0,0), &lda, &M(0,0), &ldb, alphar, alphai, beta,
	      vr, &n, vr, &n, &work0, &lwork, &info);
#endif
  }
  else {
    DGEEV_F77("N", "V", &n, &J(0,0), &lda, alphar, alphai, 
	      vr, &n, vr, &n, &work0, &lwork, &info);
  }

  // Allocate work array
  lwork = (int) work0;
  work = new double[lwork];

  // Calculate eigenvalues, eigenvectors
  if (hasMassMatrix) {
#ifdef HAVE_LAPACK_GENEV
    DGGEV_F77("N", "V", &n, &J(0,0), &lda, &M(0,0), &ldb, alphar, alphai, beta,
	      vr, &n, vr, &n, work, &lwork, &info);
#endif
  }
  else {
    DGEEV_F77("N", "V", &n, &J(0,0), &lda, alphar, alphai, 
	      vr, &n, vr, &n, work, &lwork, &info);
  }

  // Check for success
  if (info != 0)
    return NOX::Abstract::Group::Failed;

  // Print out eigenvalues
  if (LOCA::Utils::doPrint(LOCA::Utils::StepperIteration)) {
    if (hasMassMatrix) {
      cout << "Generalized eigenvalues: " << endl;
      for (int i=0; i<n; i++)
	cout << "\t" << LOCA::Utils::sci(alphar[i]/beta[i]) << " + i" << LOCA::Utils::sci(alphai[i]/beta[i]) << endl;
    }
    else {
      cout << "Eigenvalues: " << endl;
      for (int i=0; i<n; i++)
	cout << "\t" << LOCA::Utils::sci(alphar[i]) << " + i" << LOCA::Utils::sci(alphai[i]) << endl;
    }
  }

  delete [] alphar;
  delete [] alphai;
  delete [] beta;
  delete [] vr;
  delete [] work;

  return NOX::Abstract::Group::Ok;
}
