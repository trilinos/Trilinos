// $Id$
// $Source$

//@HEADER
// ************************************************************************
//
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

#include "Teuchos_ParameterList.hpp"

#include "LOCA_GlobalData.H"
#include "LOCA_Factory.H"
#include "LOCA_ErrorCheck.H"
#include "NOX_Utils.H"
#include "LOCA_Parameter_SublistParser.H"

#include "LOCA_Eigensolver_DGGEVStrategy.H"
#include "LOCA_EigenvalueSort_Strategies.H"

#include "Teuchos_LAPACK_wrappers.hpp"
#include "LOCA_LAPACK_Group.H"


LOCA::Eigensolver::DGGEVStrategy::DGGEVStrategy(
    const Teuchos::RCP<LOCA::GlobalData>& global_data,
    const Teuchos::RCP<LOCA::Parameter::SublistParser>& tpParams,
    const Teuchos::RCP<Teuchos::ParameterList>& eigParams) :
  globalData(global_data),
  topParams(tpParams),
  eigenParams(eigParams),
  nev(4),
  which("LM")
{
  nev = eigenParams->get("Num Eigenvalues", 4);
  which = eigenParams->get("Sorting Order","LM");
}

LOCA::Eigensolver::DGGEVStrategy::~DGGEVStrategy()
{
}

NOX::Abstract::Group::ReturnType
LOCA::Eigensolver::DGGEVStrategy::computeEigenvalues(
         NOX::Abstract::Group& group,
         Teuchos::RCP< std::vector<double> >& evals_r,
         Teuchos::RCP< std::vector<double> >& evals_i,
         Teuchos::RCP< NOX::Abstract::MultiVector >& evecs_r,
         Teuchos::RCP< NOX::Abstract::MultiVector >& evecs_i)
{

  // Get LAPACK group
  LOCA::LAPACK::Group & grp =
    dynamic_cast<LOCA::LAPACK::Group&>(group);

  const bool hasMassMatrix = true;

  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperIteration)) {
    globalData->locaUtils->out()
      << std::endl << globalData->locaUtils->fill(64,'=') << std::endl
      << "LAPACK ";
    if (hasMassMatrix)
      globalData->locaUtils->out() << "DGGEV ";
    else
      globalData->locaUtils->out() << "DGEEV ";
    globalData->locaUtils->out() << "Eigensolver starting."
                 << std::endl << std::endl;;
  }

  // Make sure Jacobian & mass matrices are fresh
  grp.computeJacobian();
  if (hasMassMatrix)
    grp.computeShiftedMatrix(0.0, 1.0);

  // Get data out of group
  NOX::LAPACK::Matrix<double>& jacobianMatrix = grp.getJacobianMatrix();
  NOX::LAPACK::Matrix<double>& massMatrix = grp.getShiftedMatrix();

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
  NOX::LAPACK::Matrix<double> J(jacobianMatrix);

  NOX::LAPACK::Matrix<double> M;

  // First do a workspace query
  if (hasMassMatrix) {

    // Copy mass matrix since lapack routines overwrite it
    M = massMatrix;

    DGGEV_F77("N", "V", &n, &J(0,0), &lda, &M(0,0), &ldb, alphar, alphai, beta,
          vr, &n, vr, &n, &work0, &lwork, &info);
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
    DGGEV_F77("N", "V", &n, &J(0,0), &lda, &M(0,0), &ldb, alphar, alphai, beta,
          vr, &n, vr, &n, work, &lwork, &info);
  }
  else {
    DGEEV_F77("N", "V", &n, &J(0,0), &lda, alphar, alphai,
          vr, &n, vr, &n, work, &lwork, &info);
  }

  // Check for success
  if (info != 0)
    return NOX::Abstract::Group::Failed;

  // Compute all of the eigenvalues and eigenvectors before sorting
  std::vector<double> evals_r_tmp(n);
  std::vector<double> evals_i_tmp(n);
  Teuchos::RCP<NOX::Abstract::MultiVector> evecs_r_tmp =
    group.getX().createMultiVector(n, NOX::ShapeCopy);
  Teuchos::RCP<NOX::Abstract::MultiVector>evecs_i_tmp =
    group.getX().createMultiVector(n, NOX::ShapeCopy);
  double rnext;
  double inext;
  bool isComplexEval = false;
  bool isPrevComplexEval = false;
  for (int j=0; j<n; j++) {

    // Compute eigenvalues
    if (hasMassMatrix) {
      evals_r_tmp[j] = alphar[j]/beta[j];
      evals_i_tmp[j] = alphai[j]/beta[j];
    }
    else {
      evals_r_tmp[j] = alphar[j];
      evals_i_tmp[j] = alphai[j];
    }

    // Compute next eigenvalue
    if (!isPrevComplexEval && j < n-1) {
      if (hasMassMatrix) {
    rnext = alphar[j+1]/beta[j+1];
    inext = alphai[j+1]/beta[j+1];
      }
      else {
    rnext = alphar[j+1];
    inext = alphai[j+1];
      }

      // Determine if this eigenvalue is a complex conjugate pair
      if (fabs(evals_r_tmp[j] - rnext) < 1.0e-14*fabs(1.0+evals_r_tmp[j]) &&
      fabs(evals_i_tmp[j] + inext) < 1.0e-14*fabs(1.0+evals_i_tmp[j]))
    isComplexEval = true;
      else
    isComplexEval = false;
    }
    else if (!isPrevComplexEval && j == n-1)
      isComplexEval = false;

    NOX::LAPACK::Vector & tmpr = dynamic_cast<NOX::LAPACK::Vector&>((*evecs_r_tmp)[j]);
    NOX::LAPACK::Vector & tmpi = dynamic_cast<NOX::LAPACK::Vector&>((*evecs_i_tmp)[j]);

    if (isComplexEval)
      for (int i=0; i<n; i++) {
        tmpr(i) =  vr[i+j*n];
        tmpi(i) =  vr[i+(j+1)*n];
      }
    else if (isPrevComplexEval)
      for (int i=0; i<n; i++) {
        tmpr(i) =  vr[i+(j-1)*n];
        tmpi(i) = -vr[i+j*n];
      }
    else
      for (int i=0; i<n; i++) {
        tmpr(i) = vr[i+j*n];
        tmpi(i) = 0.0;;
      }

    if (isPrevComplexEval) {
      isPrevComplexEval = false;
      isComplexEval = false;
    }
    if (isComplexEval) {
      isPrevComplexEval = true;
      isComplexEval = false;
    }

  }

  // Instantiate a sorting strategy
  Teuchos::RCP<LOCA::EigenvalueSort::AbstractStrategy> evalSort =
    globalData->locaFactory->createEigenvalueSortStrategy(topParams,
                              eigenParams);

  // Create permutation array
  std::vector<int> perm(n);

  // Sort eigenvalues
  evalSort->sort(n, &evals_r_tmp[0], &evals_i_tmp[0], &perm);

  // Get first nev entries of perm
  std::vector<int> perm_short(perm.begin(), perm.begin()+nev);

  // Get sorted eigenvalues and eigenvectors
  evals_r = Teuchos::rcp(new std::vector<double>(evals_r_tmp.begin(),
                         evals_r_tmp.begin()+nev));
  evals_i = Teuchos::rcp(new std::vector<double>(evals_i_tmp.begin(),
                         evals_i_tmp.begin()+nev));
  evecs_r = evecs_r_tmp->subCopy(perm_short);
  evecs_i = evecs_i_tmp->subCopy(perm_short);

  // Print out eigenvalues
  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperIteration)) {
    for (int i=0; i<nev; i++)
      globalData->locaUtils->out()
    << "Eigenvalue " << i << " : "
    << globalData->locaUtils->sciformat((*evals_r)[i]) << " "
    << globalData->locaUtils->sciformat((*evals_i)[i]) << " i"
    << std::endl;
  }

  if (globalData->locaUtils->isPrintType(NOX::Utils::StepperIteration))
    globalData->locaUtils->out()
      << "\nLAPACK Eigensolver finished.\n"
      << globalData->locaUtils->fill(64,'=') << "\n" << std::endl;

  delete [] alphar;
  delete [] alphai;
  delete [] beta;
  delete [] vr;
  delete [] work;

  return NOX::Abstract::Group::Ok;
}
