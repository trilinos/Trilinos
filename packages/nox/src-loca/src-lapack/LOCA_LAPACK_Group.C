// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#include "LOCA_LAPACK_Group.H"	// class definition
#include "NOX_BLAS_Wrappers.H"
#include "NOX_LAPACK_Wrappers.H"
#include "LOCA_Utils.H"
#include "Teuchos_LAPACK.hpp"

LOCA::LAPACK::Group::Group(LOCA::LAPACK::Interface& interface,
			   bool hasMassMat) : 
  NOX::LAPACK::Group(interface), 
  LOCA::Abstract::Group(),
  locaProblemInterface(interface), 
  params(),
  massMatrix(),
  hasMassMatrix(hasMassMat),
  isValidMass(false)
{
  if (hasMassMat)
    massMatrix = NOX::LAPACK::Matrix(jacobianMatrix.numRows(),
				     jacobianMatrix.numCols());
}

LOCA::LAPACK::Group::Group(NOX::Parameter::List& params,
			   LOCA::LAPACK::Interface& interface,
			   bool hasMassMat) : 
  NOX::LAPACK::Group(interface), 
  LOCA::Abstract::Group(params),
  locaProblemInterface(interface), 
  params(),
  massMatrix(),
  hasMassMatrix(hasMassMat),
  isValidMass(false)
{
  if (hasMassMat)
    massMatrix = NOX::LAPACK::Matrix(jacobianMatrix.numRows(),
				     jacobianMatrix.numCols());
}

LOCA::LAPACK::Group::Group(LOCA::LAPACK::Interface& interface,
			   int m, int n, bool hasMassMat) : 
  NOX::LAPACK::Group(interface, m, n), 
  LOCA::Abstract::Group(),
  locaProblemInterface(interface), 
  params(),
  massMatrix(),
  hasMassMatrix(hasMassMat),
  isValidMass(false)
{
  if (hasMassMat)
    massMatrix = NOX::LAPACK::Matrix(m, n);
}

LOCA::LAPACK::Group::Group(NOX::Parameter::List& params,
			   LOCA::LAPACK::Interface& interface,
			   int m, int n, bool hasMassMat) : 
  NOX::LAPACK::Group(interface, m, n), 
  LOCA::Abstract::Group(params),
  locaProblemInterface(interface), 
  params(),
  massMatrix(),
  hasMassMatrix(hasMassMat),
  isValidMass(false)
{
  if (hasMassMat)
    massMatrix = NOX::LAPACK::Matrix(m, n);
}

LOCA::LAPACK::Group::Group(const LOCA::LAPACK::Group& source, 
			   NOX::CopyType type) : 
  NOX::LAPACK::Group(source,type), 
  LOCA::Abstract::Group(source,type),
  locaProblemInterface(source.locaProblemInterface), 
  params(source.params),
  massMatrix(source.massMatrix),
  hasMassMatrix(source.hasMassMatrix),
  isValidMass(source.isValidMass)
{
}

LOCA::LAPACK::Group::~Group() 
{}

NOX::Abstract::Group& 
LOCA::LAPACK::Group::operator=(const NOX::Abstract::Group& source) {
  return operator=(dynamic_cast<const LOCA::LAPACK::Group&>(source));
}

LOCA::Abstract::Group& 
LOCA::LAPACK::Group::operator=(const LOCA::Abstract::Group& source) {
  return operator=(dynamic_cast<const LOCA::LAPACK::Group&>(source));
}

NOX::LAPACK::Group&
LOCA::LAPACK::Group::operator=(const NOX::LAPACK::Group& source) {
  return operator=(dynamic_cast<const LOCA::LAPACK::Group&>(source));
}

LOCA::LAPACK::Group& 
LOCA::LAPACK::Group::operator=(const LOCA::LAPACK::Group& source) {

  NOX::LAPACK::Group::operator=(source);
  LOCA::Abstract::Group::operator=(source);

  params = source.params;
  massMatrix = source.massMatrix;
  hasMassMatrix = source.hasMassMatrix;
  isValidMass = source.isValidMass;

  return *this;
}

NOX::Abstract::Group*
LOCA::LAPACK::Group::clone(NOX::CopyType type) const {
  return new LOCA::LAPACK::Group(*this, type);
}

NOX::Abstract::Group::ReturnType
LOCA::LAPACK::Group::computeF() {
  locaProblemInterface.setParams(params);
  return NOX::LAPACK::Group::computeF();
}

NOX::Abstract::Group::ReturnType
LOCA::LAPACK::Group::computeJacobian() {
  locaProblemInterface.setParams(params);
  return NOX::LAPACK::Group::computeJacobian();
}

NOX::Abstract::Group::ReturnType
LOCA::LAPACK::Group::applyJacobianInverseMulti(NOX::Parameter::List& params,
			    const NOX::Abstract::Vector* const* inputs,
			    NOX::Abstract::Vector** outputs, int nVecs) const
{
  if (nVecs < 1)
    return NOX::Abstract::Group::Failed;

  int m = jacobianMatrix.numRows();
  int n = jacobianMatrix.numCols();
  int lda = jacobianMatrix.numRows();
  int info;

  // Copy all input vectors into one matrix
  NOX::LAPACK::Matrix B(m,nVecs);
  const NOX::LAPACK::Vector* constVecPtr;
  for (int j=0; j<nVecs; j++) {
    constVecPtr = dynamic_cast<const NOX::LAPACK::Vector*>(inputs[j]);
    for (int i=0; i<m; i++)
      B(i,j) = (*constVecPtr)(i);
  }

  // Compute Jacobian LU factorization if invalid
  if (!NOX::LAPACK::Group::isValidJacobianLUFact) {
    NOX::LAPACK::Group::jacobianLUFact = NOX::LAPACK::Group::jacobianMatrix;
    DGETRF_F77(&m, &n, &jacobianLUFact(0,0), &lda, 
	       &NOX::LAPACK::Group::pivots[0], &info);

    if (info != 0)
      return NOX::Abstract::Group::Failed;

    NOX::LAPACK::Group::isValidJacobianLUFact = true;
  }

  // Backsolve using LU factorization
  DGETRS_F77("N", &n, &nVecs, &jacobianLUFact(0,0), &lda, &pivots[0], 
  	     &B(0,0), &m, &info);

  if (info != 0)
      return NOX::Abstract::Group::Failed;

  // Copy result from matrix
  NOX::LAPACK::Vector* vecPtr;
  for (int j=0; j<nVecs; j++) {
    vecPtr = dynamic_cast<NOX::LAPACK::Vector*>(outputs[j]);
    for (int i=0; i<m; i++)
      (*vecPtr)(i) = B(i,j);
  }

  return NOX::Abstract::Group::Ok;
}

void
LOCA::LAPACK::Group::setParams(const LOCA::ParameterVector& p) 
{
  resetIsValid();
  params = p;
}

const LOCA::ParameterVector& 
LOCA::LAPACK::Group::getParams() const
{
  return params;
}

void
LOCA::LAPACK::Group::setParam(int paramID, double val)
{
  resetIsValid();
  params.setValue(paramID, val);
}

double
LOCA::LAPACK::Group::getParam(int paramID) const
{
  return params.getValue(paramID);
}

void
LOCA::LAPACK::Group::setParam(string paramID, double val)
{
  resetIsValid();
  params.setValue(paramID, val);
}

double
LOCA::LAPACK::Group::getParam(string paramID) const
{
  return params.getValue(paramID);
}

double
LOCA::LAPACK::Group::computeScaledDotProduct(
				       const NOX::Abstract::Vector& a,
				       const NOX::Abstract::Vector& b) const
{
  return a.dot(b) / a.length();
}

void
LOCA::LAPACK::Group::scaleVector(NOX::Abstract::Vector& x) const
{
  x.scale(1.0 / sqrt(static_cast<double>(x.length())));
}

void 
LOCA::LAPACK::Group::printSolution(const double conParam) const
{
   printSolution(xVector, conParam);
}

void
LOCA::LAPACK::Group::printSolution(const NOX::LAPACK::Vector& x_,
                                   const double conParam) const
{
   locaProblemInterface.printSolution(x_, conParam);
}

void
LOCA::LAPACK::Group::printSolution(const NOX::Abstract::Vector& x_,
                                   const double conParam) const
{
   printSolution(dynamic_cast<const NOX::LAPACK::Vector&>(x_), conParam);
}

NOX::Abstract::Group::ReturnType 
LOCA::LAPACK::Group::computeEigenvalues(NOX::Parameter::List& params)
{

  // Check to make sure we have dggev available if we need generalized 
  // eigenvalues.
#ifndef HAVE_LAPACK_GENEV
  if (hasMassMatrix) {
    if (Utils::doPrint(Utils::StepperIteration)) {
      errorCheck.printWarning("LOCA::LAPACK::Group::computeEigenvalues",
				     "LAPACK Generalized eigensolver (dggev) requested but not available!");
    }
    return LOCA::Abstract::Group::Ok;
  }

#endif
  
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

  // Make sure Jacobian is fresh
  computeJacobian();

  // Copy Jacobian matrix since lapack routines overwrite it
  NOX::LAPACK::Matrix J(NOX::LAPACK::Group::jacobianMatrix);

  NOX::LAPACK::Matrix M;

  // First do a workspace query
  if (hasMassMatrix) {

    // Make sure mass matrix is fresh
    computeMassMatrix();

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

void
LOCA::LAPACK::Group::projectToDraw(const NOX::Abstract::Vector& x,
				   double *px) const
{
  const NOX::LAPACK::Vector& lx = 
    dynamic_cast<const NOX::LAPACK::Vector&>(x);
  locaProblemInterface.projectToDraw(lx, px);
}

int
LOCA::LAPACK::Group::projectToDrawDimension() const
{
  return locaProblemInterface.projectToDrawDimension();
}

NOX::Abstract::Group::ReturnType 
LOCA::LAPACK::Group::augmentJacobianForHomotopy(double conParamValue)
{
  int size = jacobianMatrix.numRows();

  // Scale the matrix by the value of the homotopy continuation param
  jacobianMatrix.scale(conParamValue);

  // Add the scaled identity matrix to the jacobian
  for (int i = 0; i < size; i++) 
    jacobianMatrix(i,i) += (1.0 - conParamValue);

  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType
LOCA::LAPACK::Group::computeMassMatrix()
{
  // Skip if the mass matrix is already valid
  if (isValidMass || !hasMassMatrix)
    return NOX::Abstract::Group::Ok;

  isValidMass = locaProblemInterface.computeMass(massMatrix, xVector);

  if (isValidMass)
    return NOX::Abstract::Group::Ok;
  else
    return NOX::Abstract::Group::Failed;
}

NOX::Abstract::Group::ReturnType
LOCA::LAPACK::Group::applyMassMatrix(const NOX::Abstract::Vector& input,
				     NOX::Abstract::Vector& result) const
{
  const NOX::LAPACK::Vector& lapackInput = 
    dynamic_cast<const NOX::LAPACK::Vector&>(input);
  NOX::LAPACK::Vector& lapackResult = 
    dynamic_cast<NOX::LAPACK::Vector&>(result);
  return applyMassMatrix(lapackInput, lapackResult);
}

NOX::Abstract::Group::ReturnType
LOCA::LAPACK::Group::applyMassMatrix(const NOX::LAPACK::Vector& input,
		      NOX::LAPACK::Vector& result) const
{
  // Check validity of the mass matrix
  if (!isMassMatrix()) 
    return NOX::Abstract::Group::BadDependency;

  // Compute result = M * input
  int m = massMatrix.numRows();
  int n = massMatrix.numCols();
  int lda = massMatrix.numRows();

  DGEMV_F77("N", &m, &n, &NOX::LAPACK::d_one, &massMatrix(0,0), &lda, 
	    &input(0), &NOX::LAPACK::i_one, &NOX::LAPACK::d_zero, &result(0), 
	    &NOX::LAPACK::i_one);

  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType
LOCA::LAPACK::Group::applyComplexInverse(
			       NOX::Parameter::List& params,
			       const NOX::Abstract::Vector& input_real,
			       const NOX::Abstract::Vector& input_imag,
			       double frequency,
			       NOX::Abstract::Vector& result_real,
			       NOX::Abstract::Vector& result_imag) const
{
   NOX::Abstract::Group::ReturnType res;
   const NOX::Abstract::Vector* inputs_real[1];
   const NOX::Abstract::Vector* inputs_imag[1];
   NOX::Abstract::Vector* results_real[1];
   NOX::Abstract::Vector* results_imag[1];

   inputs_real[0] = &input_real;
   inputs_imag[0] = &input_imag;
   results_real[0] = &result_real;
   results_imag[0] = &result_imag;
   
   res = applyComplexInverseMulti(params, inputs_real, inputs_imag, frequency, 
				  results_real, results_imag, 1);

   return res;
}

NOX::Abstract::Group::ReturnType
LOCA::LAPACK::Group::applyComplexInverseMulti(
			       NOX::Parameter::List& params,
			       const NOX::Abstract::Vector* const* inputs_real,
			       const NOX::Abstract::Vector* const* inputs_imag,
			       double frequency,
			       NOX::Abstract::Vector** results_real,
			       NOX::Abstract::Vector** results_imag,
			       int nVecs) const
{
  if (!isMassMatrix()) 
    return NOX::Abstract::Group::BadDependency;

  if (nVecs < 1)
    return NOX::Abstract::Group::Failed;

  int n = jacobianMatrix.numRows();
  int m = nVecs;
  int info;
  int *piv = new int[n];

  // Copy all input vectors into one (complex) matrix
  complex<double>* B = new complex<double>[n*m];
  const NOX::LAPACK::Vector* constVecPtrR;
  const NOX::LAPACK::Vector* constVecPtrI;
  for (int j=0; j<m; j++) {
    constVecPtrR = dynamic_cast<const NOX::LAPACK::Vector*>(inputs_real[j]);
    constVecPtrI = dynamic_cast<const NOX::LAPACK::Vector*>(inputs_imag[j]);
    for (int i=0; i<n; i++) {
      B[i+n*j] = complex<double>((*constVecPtrR)(i), (*constVecPtrI)(i));
    }
  }

  // Create complex matrix J+i*w*M
  complex<double>* A = new complex<double>[n*n];
  for (int j=0; j<n; j++) {
    for (int i=0; i<n; i++) {
      A[i+n*j] = 
	complex<double>(jacobianMatrix(i,j), frequency*massMatrix(i,j));
    }
  }

  // Solve A*X = B
  Teuchos::LAPACK< int,complex<double> > L;
  L.GESV(n, m, A, n, piv, B, n, &info);
  //ZGESV_F77(&n, &m, &A(0,0), &n, piv, &B(0,0), &n, &info);

  if (info != 0)
      return NOX::Abstract::Group::Failed;

  // Copy result from matrix
  NOX::LAPACK::Vector* vecPtrR;
  NOX::LAPACK::Vector* vecPtrI;
  for (int j=0; j<m; j++) {
    vecPtrR = dynamic_cast<NOX::LAPACK::Vector*>(results_real[j]);
    vecPtrI = dynamic_cast<NOX::LAPACK::Vector*>(results_imag[j]);
    for (int i=0; i<n; i++) {
      (*vecPtrR)(i) = B[i+n*j].real();
      (*vecPtrI)(i) = B[i+n*j].imag();
    }
  }

  delete [] piv;
  delete [] A;
  delete [] B;

  return NOX::Abstract::Group::Ok;
}

bool
LOCA::LAPACK::Group::hasMass() const 
{
  return hasMassMatrix;
}

bool
LOCA::LAPACK::Group::isMassMatrix() const 
{
  return isValidMass && hasMassMatrix;
}

void
LOCA::LAPACK::Group::resetIsValid()
{
  NOX::LAPACK::Group::resetIsValid();
  isValidMass = false;
}

NOX::Abstract::Group::ReturnType
LOCA::LAPACK::Group::applyBorderedJacobianInverse(bool useTrans,
				   NOX::Parameter::List& params,
				   const NOX::Abstract::Vector& a,
				   const NOX::Abstract::Vector& b,
				   const NOX::Abstract::Vector& vInput,
				   double sInput,
				   NOX::Abstract::Vector& vResult,
				   double& sResult) const
{
  if (!isJacobian()) {
    cerr << "ERROR: NOX::LAPACK::Group::applyJacobianInverse() - invalid Jacobian" << endl;
    throw "NOX Error";
  }

  // Cast input, results to lapack vectors
  const NOX::LAPACK::Vector& lapack_vInput = 
    dynamic_cast<const NOX::LAPACK::Vector&>(vInput);
  const NOX::LAPACK::Vector& lapack_a = 
    dynamic_cast<const NOX::LAPACK::Vector&>(a);
  const NOX::LAPACK::Vector& lapack_b = 
    dynamic_cast<const NOX::LAPACK::Vector&>(b);
  NOX::LAPACK::Vector& lapack_vResult = 
    dynamic_cast<NOX::LAPACK::Vector&>(vResult);

  int m = jacobianMatrix.numRows();
  int n = jacobianMatrix.numCols();
  int info;

  int M = m+1;
  int N = n+1;

  const char *trans = "N";
  if (useTrans)
    trans = "T";

  // Create and fill extended matrix
  NOX::LAPACK::Matrix extendedJac(m+1, n+1);
  for (int j=0; j<n; j++) {
    for (int i=0; i<m; i++) 
      extendedJac(i,j) = jacobianMatrix(i,j);
    extendedJac(m, j) = lapack_b(j);
  }
  for (int i=0; i<m; i++)
    extendedJac(i, n) = lapack_a(i);
  extendedJac(m,n) = 0.0;

  // Create and fill extended RHS
  double *extendedInput = new double[m+1];
  for (int i=0; i<m; i++)
    extendedInput[i] = lapack_vInput(i);
  extendedInput[m] = sInput;

  // Create extended pivot array
  int *extendedPivots = new int[m+1];
  for (int i=0; i<m+1; i++)
    extendedPivots[i] = 0;

  // Factor extended matrix
  DGETRF_F77(&M, &N, &extendedJac(0,0), &M, &extendedPivots[0], &info);

  if (info != 0)
    return NOX::Abstract::Group::Failed;

  // Backsolve using LU factorization
  DGETRS_F77(trans, &N, &NOX::LAPACK::i_one, &extendedJac(0,0), &M, 
	     &extendedPivots[0], extendedInput, &M, &info);

  // Copy out extended solution
  for (int i=0; i<m; i++)
    lapack_vResult(i) = extendedInput[i];
  sResult = extendedInput[m];

  delete [] extendedInput;
  delete [] extendedPivots;
    
  return (info == 0) ? (NOX::Abstract::Group::Ok) : (NOX::Abstract::Group::Failed);

}
