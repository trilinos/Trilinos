// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_LAPACK_Group.H"    // class definition
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"
#include "Teuchos_LAPACK.hpp"

LOCA::LAPACK::Group::Group(
            const Teuchos::RCP<LOCA::GlobalData>& global_data,
            LOCA::LAPACK::Interface& interface) :
  NOX::LAPACK::Group(interface),
  LOCA::Abstract::Group(global_data),
  globalData(global_data),
  locaProblemInterface(interface),
  params(),
  shiftedSolver(jacSolver.getMatrix().numRows()),
  freq(0.0),
  isValidComplex(false)
#ifdef HAVE_TEUCHOS_COMPLEX
  ,complexSolver(jacSolver.getMatrix().numRows())
#endif
{
}

LOCA::LAPACK::Group::Group(const LOCA::LAPACK::Group& source,
               NOX::CopyType type) :
  NOX::LAPACK::Group(source,type),
  LOCA::Abstract::Group(source,type),
  globalData(source.globalData),
  locaProblemInterface(source.locaProblemInterface),
  params(source.params),
  shiftedSolver(source.shiftedSolver),
  freq(source.freq),
  isValidComplex(source.isValidComplex)
#ifdef HAVE_TEUCHOS_COMPLEX
  ,complexSolver(source.complexSolver)
#endif
{
}

LOCA::LAPACK::Group::~Group()
{}

LOCA::LAPACK::Group&
LOCA::LAPACK::Group::operator=(const LOCA::LAPACK::Group& source) {

  NOX::LAPACK::Group::operator=(source);
  LOCA::Abstract::Group::copy(source);

  globalData = source.globalData;
  params = source.params;
  shiftedSolver = source.shiftedSolver;
  freq = source.freq;
  isValidComplex = source.isValidComplex;
#ifdef HAVE_TEUCHOS_COMPLEX
  complexSolver = source.complexSolver;
#endif

  return *this;
}

NOX::Abstract::Group&
LOCA::LAPACK::Group::operator=(const NOX::Abstract::Group& source) {
  operator=(dynamic_cast<const LOCA::LAPACK::Group&>(source));
  return *this;
}

NOX::LAPACK::Group&
LOCA::LAPACK::Group::operator=(const NOX::LAPACK::Group& source) {
  operator=(dynamic_cast<const LOCA::LAPACK::Group&>(source));
  return *this;
}

Teuchos::RCP<NOX::Abstract::Group>
LOCA::LAPACK::Group::clone(NOX::CopyType type) const {
  return Teuchos::rcp(new LOCA::LAPACK::Group(*this, type));
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
LOCA::LAPACK::Group::applyJacobianTransposeInverse(
                     Teuchos::ParameterList& /* p */,
                     const NOX::Abstract::Vector& input,
                     NOX::Abstract::Vector& result) const
{

  if (!isJacobian()) {
    std::cerr << "ERROR: "
     << "LOCA::LAPACK::Group::applyJacobianTransposeInverse()"
     << " - invalid Jacobian" << std::endl;
    throw std::runtime_error("NOX Error");
  }

  const NOX::LAPACK::Vector& lapack_input =
    dynamic_cast<const NOX::LAPACK::Vector&> (input);
  NOX::LAPACK::Vector& lapack_result =
    dynamic_cast<NOX::LAPACK::Vector&> (result);

  // Solve Jacobian transpose
  lapack_result = lapack_input;
  bool res = jacSolver.solve(true, 1, &lapack_result(0));

  return res ? (NOX::Abstract::Group::Ok) : (NOX::Abstract::Group::Failed);
}

NOX::Abstract::Group::ReturnType
LOCA::LAPACK::Group::applyJacobianTransposeInverseMultiVector(
                     Teuchos::ParameterList& /* p */,
                     const NOX::Abstract::MultiVector& input,
                     NOX::Abstract::MultiVector& result) const
{

  if (!isJacobian()) {
    std::cerr << "ERROR: "
     << "LOCA::LAPACK::Group::applyJacobianTransposeInverseMultiVector()"
     << " - invalid Jacobian" << std::endl;
    throw std::runtime_error("NOX Error");
  }

  // Number of RHS
  int nVecs = input.numVectors();

  int m = jacSolver.getMatrix().numRows();

  // Copy all input vectors into one matrix
  NOX::LAPACK::Matrix<double> B(m,nVecs);
  const NOX::LAPACK::Vector* constVecPtr;
  for (int j=0; j<nVecs; j++) {
    constVecPtr = dynamic_cast<const NOX::LAPACK::Vector*>(&(input[j]));
    TEUCHOS_ASSERT(constVecPtr != NULL);
    for (int i=0; i<m; i++)
      B(i,j) = (*constVecPtr)(i);
  }

  // Solve Jacobian transpose
  bool res = jacSolver.solve(true, nVecs, &B(0,0));

  if (!res)
    return NOX::Abstract::Group::Failed;

  // Copy result from matrix
  NOX::LAPACK::Vector* vecPtr;
  for (int j=0; j<nVecs; j++) {
    vecPtr = dynamic_cast<NOX::LAPACK::Vector*>(&(result[j]));
    TEUCHOS_ASSERT(vecPtr != NULL);
    for (int i=0; i<m; i++)
      (*vecPtr)(i) = B(i,j);
  }

  return NOX::Abstract::Group::Ok;
}

void
LOCA::LAPACK::Group::copy(const NOX::Abstract::Group& source) {
  *this = source;
}

void
LOCA::LAPACK::Group::setParams(const LOCA::ParameterVector& p)
{
  resetIsValid();
  params = p;
}

void
LOCA::LAPACK::Group::setParam(int paramID, double val)
{
  resetIsValid();
  params.setValue(paramID, val);
}

void
LOCA::LAPACK::Group::setParam(std::string paramID, double val)
{
  resetIsValid();
  params.setValue(paramID, val);
}

const LOCA::ParameterVector&
LOCA::LAPACK::Group::getParams() const
{
  return params;
}

double
LOCA::LAPACK::Group::getParam(int paramID) const
{
  return params.getValue(paramID);
}

double
LOCA::LAPACK::Group::getParam(std::string paramID) const
{
  return params.getValue(paramID);
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

double
LOCA::LAPACK::Group::computeScaledDotProduct(
                       const NOX::Abstract::Vector& a,
                       const NOX::Abstract::Vector& b) const
{
  return a.innerProduct(b) / a.length();
}

void
LOCA::LAPACK::Group::printSolution(const double conParam) const
{
  printSolution(xVector, conParam);
}

void
LOCA::LAPACK::Group::printSolution(const NOX::Abstract::Vector& x_,
                                   const double conParam) const
{
  locaProblemInterface.printSolution(dynamic_cast<const NOX::LAPACK::Vector&>(x_), conParam);
}

void
LOCA::LAPACK::Group::scaleVector(NOX::Abstract::Vector& x) const
{
  x.scale(1.0 / sqrt(static_cast<double>(x.length())));
}

NOX::Abstract::Group::ReturnType
LOCA::LAPACK::Group::computeShiftedMatrix(double alpha, double beta)
{
  // Compute alpha*J+beta*M
  bool res =
    locaProblemInterface.computeShiftedMatrix(alpha, beta, xVector,
                          shiftedSolver.getMatrix());

  if (res)
    return NOX::Abstract::Group::Ok;
  else
    return NOX::Abstract::Group::Failed;
}

NOX::Abstract::Group::ReturnType
LOCA::LAPACK::Group::applyShiftedMatrix(const NOX::Abstract::Vector& input,
                    NOX::Abstract::Vector& result) const
{
  // Cast inputs to LAPACK vectors
  const NOX::LAPACK::Vector& lapack_input =
    dynamic_cast<const NOX::LAPACK::Vector&>(input);
  NOX::LAPACK::Vector& lapack_result =
    dynamic_cast<NOX::LAPACK::Vector&>(result);

  // Apply shifted matrix
  shiftedSolver.apply(false, 1, &lapack_input(0), &lapack_result(0));

  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType
LOCA::LAPACK::Group::applyShiftedMatrixMultiVector(
                   const NOX::Abstract::MultiVector& input,
                   NOX::Abstract::MultiVector& result) const
{
  // Number of RHS
  int nVecs = input.numVectors();

  int m = shiftedSolver.getMatrix().numRows();

  // Copy all input vectors into one matrix
  NOX::LAPACK::Matrix<double> B(m,nVecs);
  NOX::LAPACK::Matrix<double> C(m,nVecs);
  const NOX::LAPACK::Vector* constVecPtr;
  for (int j=0; j<nVecs; j++) {
    constVecPtr = dynamic_cast<const NOX::LAPACK::Vector*>(&(input[j]));
    TEUCHOS_ASSERT(constVecPtr != NULL);
    for (int i=0; i<m; i++)
      B(i,j) = (*constVecPtr)(i);
  }

  // Apply shifted matrix
  shiftedSolver.apply(false, nVecs, &B(0,0), &C(0,0));

  // Copy result from matrix
  NOX::LAPACK::Vector* vecPtr;
  for (int j=0; j<nVecs; j++) {
    vecPtr = dynamic_cast<NOX::LAPACK::Vector*>(&(result[j]));
    TEUCHOS_ASSERT(vecPtr != NULL);
    for (int i=0; i<m; i++)
      (*vecPtr)(i) = C(i,j);
  }

  return NOX::Abstract::Group::Ok;
}

NOX::Abstract::Group::ReturnType
LOCA::LAPACK::Group::applyShiftedMatrixInverseMultiVector(
                         Teuchos::ParameterList& /* params */,
                     const NOX::Abstract::MultiVector& input,
                     NOX::Abstract::MultiVector& result) const
{
  // Number of RHS
  int nVecs = input.numVectors();

  int m = shiftedSolver.getMatrix().numRows();

  // Copy all input vectors into one matrix
  NOX::LAPACK::Matrix<double> B(m,nVecs);
  const NOX::LAPACK::Vector* constVecPtr;
  for (int j=0; j<nVecs; j++) {
    constVecPtr = dynamic_cast<const NOX::LAPACK::Vector*>(&(input[j]));
    TEUCHOS_ASSERT(constVecPtr != NULL);
    for (int i=0; i<m; i++)
      B(i,j) = (*constVecPtr)(i);
  }

  bool res = shiftedSolver.solve(false, nVecs, &B(0,0));

  if (!res)
    return NOX::Abstract::Group::Failed;

  // Copy result from matrix
  NOX::LAPACK::Vector* vecPtr;
  for (int j=0; j<nVecs; j++) {
    vecPtr = dynamic_cast<NOX::LAPACK::Vector*>(&(result[j]));
    for (int i=0; i<m; i++)
      (*vecPtr)(i) = B(i,j);
  }

  return NOX::Abstract::Group::Ok;
}

bool
LOCA::LAPACK::Group::isComplex() const
{
  return isValidComplex;
}

NOX::Abstract::Group::ReturnType
LOCA::LAPACK::Group::computeComplex(double frequency)
{
  (void)frequency;
  std::string callingFunction = "LOCA::LAPACK::computeComplex()";

#ifdef HAVE_TEUCHOS_COMPLEX
  NOX::Abstract::Group::ReturnType finalStatus;

  freq = frequency;

  // Compute Jacobian
  finalStatus = computeJacobian();
  globalData->locaErrorCheck->checkReturnType(finalStatus, callingFunction);

  // Compute Mass matrix
  bool res =
    locaProblemInterface.computeShiftedMatrix(0.0, 1.0, xVector,
                          shiftedSolver.getMatrix());

  // Compute complex matrix
  NOX::LAPACK::Matrix<double>& jacobianMatrix = jacSolver.getMatrix();
  NOX::LAPACK::Matrix<double>& massMatrix = shiftedSolver.getMatrix();
  NOX::LAPACK::Matrix< std::complex<double> >& complexMatrix =
    complexSolver.getMatrix();
  int n = jacobianMatrix.numRows();
  for (int j=0; j<n; j++) {
    for (int i=0; i<n; i++) {
      complexMatrix(i,j) =
    std::complex<double>(jacobianMatrix(i,j), frequency*massMatrix(i,j));
    }
  }

  if (finalStatus == NOX::Abstract::Group::Ok && res)
    isValidComplex = true;

  if (res)
    return finalStatus;
  else
    return NOX::Abstract::Group::Failed;
#else
  globalData->locaErrorCheck->throwError(
    callingFunction,
    "TEUCHOS_COMPLEX must be enabled for complex support!  Reconfigure with -D Teuchos_ENABLE_COMPLEX");
  return NOX::Abstract::Group::BadDependency;
#endif
}

NOX::Abstract::Group::ReturnType
LOCA::LAPACK::Group::applyComplex(const NOX::Abstract::Vector& input_real,
                  const NOX::Abstract::Vector& input_imag,
                  NOX::Abstract::Vector& result_real,
                  NOX::Abstract::Vector& result_imag) const
{
  (void)input_real;
  (void)input_imag;
  (void)result_real;
  (void)result_imag;
#ifdef HAVE_TEUCHOS_COMPLEX
   // Check validity of the Jacobian
  if (!isComplex())
    return NOX::Abstract::Group::BadDependency;

  int n = complexSolver.getMatrix().numCols();

  // Copy inputs into a complex vector
  std::vector< std::complex<double> > input(n);
  std::vector< std::complex<double> > result(n);
  const NOX::LAPACK::Vector& lapack_input_real =
    dynamic_cast<const NOX::LAPACK::Vector&>(input_real);
  const NOX::LAPACK::Vector& lapack_input_imag =
    dynamic_cast<const NOX::LAPACK::Vector&>(input_imag);
  for (int i=0; i<n; i++)
    input[i] = std::complex<double>(lapack_input_real(i),
                    lapack_input_imag(i));

  // Apply complex matrix
  complexSolver.apply(false, 1, &input[0], &result[0]);

  // Copy result into NOX vectors
  NOX::LAPACK::Vector& lapack_result_real =
    dynamic_cast<NOX::LAPACK::Vector&>(result_real);
  NOX::LAPACK::Vector& lapack_result_imag =
    dynamic_cast<NOX::LAPACK::Vector&>(result_imag);
  for (int i=0; i<n; i++) {
    lapack_result_real(i) = result[i].real();
    lapack_result_imag(i) = result[i].imag();
  }

  return NOX::Abstract::Group::Ok;
#else
  globalData->locaErrorCheck->throwError(
    "LOCA::LAPACK::Group::applyComplex()",
    "TEUCHOS_COMPLEX must be enabled for complex support!  Reconfigure with -D Teuchos_ENABLE_COMPLEX");
  return NOX::Abstract::Group::BadDependency;
#endif
}

NOX::Abstract::Group::ReturnType
LOCA::LAPACK::Group::applyComplexMultiVector(
                const NOX::Abstract::MultiVector& input_real,
                const NOX::Abstract::MultiVector& input_imag,
                NOX::Abstract::MultiVector& result_real,
                NOX::Abstract::MultiVector& result_imag) const
{
  (void)input_real;
  (void)input_imag;
  (void)result_real;
  (void)result_imag;
#ifdef HAVE_TEUCHOS_COMPLEX
   // Check validity of the Jacobian
  if (!isComplex())
    return NOX::Abstract::Group::BadDependency;

  int n = complexSolver.getMatrix().numRows();
  int p = input_real.numVectors();

  // Copy inputs into a complex vector
  std::vector< std::complex<double> > input(n*p);
  std::vector< std::complex<double> > result(n*p);
  const NOX::LAPACK::Vector* lapack_input_real;
  const NOX::LAPACK::Vector* lapack_input_imag;
  for (int j=0; j<p; j++) {
    lapack_input_real =
      dynamic_cast<const NOX::LAPACK::Vector*>(&(input_real[j]));
    TEUCHOS_ASSERT(lapack_input_real != NULL);
    lapack_input_imag =
      dynamic_cast<const NOX::LAPACK::Vector*>(&(input_imag[j]));
    TEUCHOS_ASSERT(lapack_input_imag != NULL);
    for (int i=0; i<n; i++)
      input[i+n*j] = std::complex<double>((*lapack_input_real)(i),
                      (*lapack_input_imag)(i));
  }

  // Apply complex matrix
  complexSolver.apply(false, p, &input[0], &result[0]);

  // Copy result into NOX vectors
  NOX::LAPACK::Vector* lapack_result_real;
  NOX::LAPACK::Vector* lapack_result_imag;
  for (int j=0; j<p; j++) {
    lapack_result_real =
      dynamic_cast<NOX::LAPACK::Vector*>(&(result_real[j]));
    TEUCHOS_ASSERT(lapack_result_real != NULL);
    lapack_result_imag =
      dynamic_cast<NOX::LAPACK::Vector*>(&(result_imag[j]));
    TEUCHOS_ASSERT(lapack_result_imag != NULL);
    for (int i=0; i<n; i++) {
      (*lapack_result_real)(i) = result[i+n*j].real();
      (*lapack_result_imag)(i) = result[i+n*j].imag();
    }
  }

  return NOX::Abstract::Group::Ok;
#else
  globalData->locaErrorCheck->throwError(
    "LOCA::LAPACK::Group::applyComplexMultiVector()",
    "TEUCHOS_COMPLEX must be enabled for complex support!  Reconfigure with -D Teuchos_ENABLE_COMPLEX");
  return NOX::Abstract::Group::BadDependency;
#endif
}

NOX::Abstract::Group::ReturnType
LOCA::LAPACK::Group::applyComplexInverseMultiVector(
                Teuchos::ParameterList& /* params */,
                const NOX::Abstract::MultiVector& input_real,
                const NOX::Abstract::MultiVector& input_imag,
                NOX::Abstract::MultiVector& result_real,
                NOX::Abstract::MultiVector& result_imag) const
{
  (void)input_real;
  (void)input_imag;
  (void)result_real;
  (void)result_imag;
#ifdef HAVE_TEUCHOS_COMPLEX
   // Check validity of the Jacobian
  if (!isComplex())
    return NOX::Abstract::Group::BadDependency;

  int n = complexSolver.getMatrix().numRows();
  int p = input_real.numVectors();

  // Copy inputs into a complex vector
  std::vector< std::complex<double> > input(n*p);
  const NOX::LAPACK::Vector* lapack_input_real;
  const NOX::LAPACK::Vector* lapack_input_imag;
  for (int j=0; j<p; j++) {
    lapack_input_real =
      dynamic_cast<const NOX::LAPACK::Vector*>(&(input_real[j]));
    TEUCHOS_ASSERT(lapack_input_real != NULL);
    lapack_input_imag =
      dynamic_cast<const NOX::LAPACK::Vector*>(&(input_imag[j]));
    TEUCHOS_ASSERT(lapack_input_imag != NULL);
    for (int i=0; i<n; i++)
      input[i+n*j] = std::complex<double>((*lapack_input_real)(i),
                      (*lapack_input_imag)(i));
  }

  // Solve complex matrix
  bool res = complexSolver.solve(false, p, &input[0]);

  // Copy result into NOX vectors
  NOX::LAPACK::Vector* lapack_result_real;
  NOX::LAPACK::Vector* lapack_result_imag;
  for (int j=0; j<p; j++) {
    lapack_result_real =
      dynamic_cast<NOX::LAPACK::Vector*>(&(result_real[j]));
    TEUCHOS_ASSERT(lapack_result_real != NULL);
    lapack_result_imag =
      dynamic_cast<NOX::LAPACK::Vector*>(&(result_imag[j]));
    TEUCHOS_ASSERT(lapack_result_imag != NULL);
    for (int i=0; i<n; i++) {
      (*lapack_result_real)(i) = input[i+n*j].real();
      (*lapack_result_imag)(i) = input[i+n*j].imag();
    }
  }

  if (res)
    return NOX::Abstract::Group::Ok;
  else
    return NOX::Abstract::Group::Failed;
#else
  globalData->locaErrorCheck->throwError(
    "LOCA::LAPACK::Group::applyComplexInverseMultiVector()",
    "TEUCHOS_COMPLEX must be enabled for complex support!  Reconfigure with -D Teuchos_ENABLE_COMPLEX");
  return NOX::Abstract::Group::BadDependency;
#endif
}

NOX::Abstract::Group::ReturnType
LOCA::LAPACK::Group::applyComplexTranspose(
                  const NOX::Abstract::Vector& input_real,
                  const NOX::Abstract::Vector& input_imag,
                  NOX::Abstract::Vector& result_real,
                  NOX::Abstract::Vector& result_imag) const
{
  (void)input_real;
  (void)input_imag;
  (void)result_real;
  (void)result_imag;
#ifdef HAVE_TEUCHOS_COMPLEX
   // Check validity of the Jacobian
  if (!isComplex())
    return NOX::Abstract::Group::BadDependency;

  int n = complexSolver.getMatrix().numCols();

  // Copy inputs into a complex vector
  std::vector< std::complex<double> > input(n);
  std::vector< std::complex<double> > result(n);
  const NOX::LAPACK::Vector& lapack_input_real =
    dynamic_cast<const NOX::LAPACK::Vector&>(input_real);
  const NOX::LAPACK::Vector& lapack_input_imag =
    dynamic_cast<const NOX::LAPACK::Vector&>(input_imag);
  for (int i=0; i<n; i++)
    input[i] = std::complex<double>(lapack_input_real(i),
                    lapack_input_imag(i));

  // Apply complex matrix
  complexSolver.apply(true, 1, &input[0], &result[0]);

  // Copy result into NOX vectors
  NOX::LAPACK::Vector& lapack_result_real =
    dynamic_cast<NOX::LAPACK::Vector&>(result_real);
  NOX::LAPACK::Vector& lapack_result_imag =
    dynamic_cast<NOX::LAPACK::Vector&>(result_imag);
  for (int i=0; i<n; i++) {
    lapack_result_real(i) = result[i].real();
    lapack_result_imag(i) = result[i].imag();
  }

  return NOX::Abstract::Group::Ok;
#else
  globalData->locaErrorCheck->throwError(
    "LOCA::LAPACK::Group::applyComplexTranspose()",
    "TEUCHOS_COMPLEX must be enabled for complex support!  Reconfigure with -D Teuchos_ENABLE_COMPLEX");
  return NOX::Abstract::Group::BadDependency;
#endif
}

NOX::Abstract::Group::ReturnType
LOCA::LAPACK::Group::applyComplexTransposeMultiVector(
                const NOX::Abstract::MultiVector& input_real,
                const NOX::Abstract::MultiVector& input_imag,
                NOX::Abstract::MultiVector& result_real,
                NOX::Abstract::MultiVector& result_imag) const
{
  (void)input_real;
  (void)input_imag;
  (void)result_real;
  (void)result_imag;
#ifdef HAVE_TEUCHOS_COMPLEX
   // Check validity of the Jacobian
  if (!isComplex())
    return NOX::Abstract::Group::BadDependency;

  int n = complexSolver.getMatrix().numRows();
  int p = input_real.numVectors();

  // Copy inputs into a complex vector
  std::vector< std::complex<double> > input(n*p);
  std::vector< std::complex<double> > result(n*p);
  const NOX::LAPACK::Vector* lapack_input_real;
  const NOX::LAPACK::Vector* lapack_input_imag;
  for (int j=0; j<p; j++) {
    lapack_input_real =
      dynamic_cast<const NOX::LAPACK::Vector*>(&(input_real[j]));
    TEUCHOS_ASSERT(lapack_input_real != NULL);
    lapack_input_imag =
      dynamic_cast<const NOX::LAPACK::Vector*>(&(input_imag[j]));
    TEUCHOS_ASSERT(lapack_input_imag != NULL);
    for (int i=0; i<n; i++)
      input[i+n*j] = std::complex<double>((*lapack_input_real)(i),
                      (*lapack_input_imag)(i));
  }

  // Apply complex matrix
  complexSolver.apply(true, p, &input[0], &result[0]);

  // Copy result into NOX vectors
  NOX::LAPACK::Vector* lapack_result_real;
  NOX::LAPACK::Vector* lapack_result_imag;
  for (int j=0; j<p; j++) {
    lapack_result_real =
      dynamic_cast<NOX::LAPACK::Vector*>(&(result_real[j]));
    TEUCHOS_ASSERT(lapack_result_real != NULL);
    lapack_result_imag =
      dynamic_cast<NOX::LAPACK::Vector*>(&(result_imag[j]));
    TEUCHOS_ASSERT(lapack_result_imag != NULL);
    for (int i=0; i<n; i++) {
      (*lapack_result_real)(i) = result[i+n*j].real();
      (*lapack_result_imag)(i) = result[i+n*j].imag();
    }
  }

  return NOX::Abstract::Group::Ok;
#else
  globalData->locaErrorCheck->throwError(
    "LOCA::LAPACK::Group::applyComplexTransposeMultiVector()",
    "TEUCHOS_COMPLEX must be enabled for complex support!  Reconfigure with -D Teuchos_ENABLE_COMPLEX");
  return NOX::Abstract::Group::BadDependency;
#endif
}

NOX::Abstract::Group::ReturnType
LOCA::LAPACK::Group::applyComplexTransposeInverseMultiVector(
                Teuchos::ParameterList& /* params */,
                const NOX::Abstract::MultiVector& input_real,
                const NOX::Abstract::MultiVector& input_imag,
                NOX::Abstract::MultiVector& result_real,
                NOX::Abstract::MultiVector& result_imag) const
{
  (void)input_real;
  (void)input_imag;
  (void)result_real;
  (void)result_imag;
#ifdef HAVE_TEUCHOS_COMPLEX
   // Check validity of the Jacobian
  if (!isComplex())
    return NOX::Abstract::Group::BadDependency;

  int n = complexSolver.getMatrix().numRows();
  int p = input_real.numVectors();

  // Copy inputs into a complex vector
  std::vector< std::complex<double> > input(n*p);
  const NOX::LAPACK::Vector* lapack_input_real;
  const NOX::LAPACK::Vector* lapack_input_imag;
  for (int j=0; j<p; j++) {
    lapack_input_real =
      dynamic_cast<const NOX::LAPACK::Vector*>(&(input_real[j]));
    TEUCHOS_ASSERT(lapack_input_real != NULL);
    lapack_input_imag =
      dynamic_cast<const NOX::LAPACK::Vector*>(&(input_imag[j]));
    TEUCHOS_ASSERT(lapack_input_imag != NULL);
    for (int i=0; i<n; i++)
      input[i+n*j] = std::complex<double>((*lapack_input_real)(i),
                      (*lapack_input_imag)(i));
  }

  // Solve complex matrix
  bool res = complexSolver.solve(true, p, &input[0]);

  // Copy result into NOX vectors
  NOX::LAPACK::Vector* lapack_result_real;
  NOX::LAPACK::Vector* lapack_result_imag;
  for (int j=0; j<p; j++) {
    lapack_result_real =
      dynamic_cast<NOX::LAPACK::Vector*>(&(result_real[j]));
    TEUCHOS_ASSERT(lapack_result_real != NULL);
    lapack_result_imag =
      dynamic_cast<NOX::LAPACK::Vector*>(&(result_imag[j]));
    TEUCHOS_ASSERT(lapack_result_imag != NULL);
    for (int i=0; i<n; i++) {
      (*lapack_result_real)(i) = input[i+n*j].real();
      (*lapack_result_imag)(i) = input[i+n*j].imag();
    }
  }

  if (res)
    return NOX::Abstract::Group::Ok;
  else
    return NOX::Abstract::Group::Failed;
#else
  globalData->locaErrorCheck->throwError(
    "LOCA::LAPACK::Group::applyComplexTransposeInverseMultiVector()",
    "TEUCHOS_COMPLEX must be enabled for complex support!  Reconfigure with -D Teuchos_ENABLE_COMPLEX");
  return NOX::Abstract::Group::BadDependency;
#endif
}

NOX::Abstract::Group::ReturnType
LOCA::LAPACK::Group::augmentJacobianForHomotopy(double a, double b)
{
  NOX::LAPACK::Matrix<double>& jacobianMatrix = jacSolver.getMatrix();
  int size = jacobianMatrix.numRows();

  // Scale the matrix by the value of the homotopy continuation param
  jacobianMatrix.scale(a);

  // Add the scaled identity matrix to the jacobian
  for (int i = 0; i < size; i++)
    jacobianMatrix(i,i) += b;

  return NOX::Abstract::Group::Ok;
}

void
LOCA::LAPACK::Group::resetIsValid()
{
  NOX::LAPACK::Group::resetIsValid();
  shiftedSolver.reset(); // Reset factorization
  isValidComplex = false;
#ifdef HAVE_TEUCHOS_COMPLEX
  complexSolver.reset(); // Reset factorization
#endif
}
