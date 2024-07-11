// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


// Class Header
#include "LOCA_DerivUtils.H"

#include "LOCA_MultiContinuation_AbstractGroup.H"
#include "LOCA_Hopf_MooreSpence_AbstractGroup.H"
#include "LOCA_Hopf_MinimallyAugmented_AbstractGroup.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Abstract_MultiVector.H"
#include "LOCA_Parameter_Vector.H"
#include "NOX_Common.H"  // For fabs function
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"

LOCA::DerivUtils::DerivUtils(
            const Teuchos::RCP<LOCA::GlobalData>& global_data,
            double _perturb) :
  globalData(global_data),
  perturb(_perturb)
{
  // assert (perturb>0.0);
}

LOCA::DerivUtils::DerivUtils(const DerivUtils& source) :
  globalData(source.globalData),
  perturb(source.perturb)
{

}

LOCA::DerivUtils::~DerivUtils()
{

}

Teuchos::RCP<LOCA::DerivUtils>
LOCA::DerivUtils::clone(NOX::CopyType /* type */) const
{
  return Teuchos::rcp(new DerivUtils(*this));  //Call Copy Constructor
}

NOX::Abstract::Group::ReturnType
LOCA::DerivUtils::computeDfDp(LOCA::MultiContinuation::AbstractGroup& grp,
                  const std::vector<int>& param_ids,
                  NOX::Abstract::MultiVector& result,
                  bool isValidF) const
{
  std::string callingFunction =
    "LOCA::DerivUtils::computeDfDp()";
  NOX::Abstract::Group::ReturnType status, finalStatus;

  // Views of f, df/dp
  NOX::Abstract::Vector *f = &result[0];
  NOX::Abstract::Vector *dfdp = NULL;

  // Compute base residual F
  if (!isValidF) {
    finalStatus = grp.computeF();
    globalData->locaErrorCheck->checkReturnType(finalStatus, callingFunction);
    *f = grp.getF();
  }
  else
    finalStatus = NOX::Abstract::Group::Ok;

  double param;
  double eps;

  // Loop over each parameter
  for (unsigned int i=0; i<param_ids.size(); i++) {

    // Perturb single parameter in this group, and return perturbation, eps
    eps = perturbParam(grp, param, param_ids[i]);

    // Compute perturbed residual
    status = grp.computeF();
    finalStatus =
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                                 finalStatus,
                                 callingFunction);

    // Difference perturbed and base vector
    dfdp = &result[i+1];
    dfdp->update(1.0, grp.getF(), -1.0, *f, 0.0);
    dfdp->scale(1.0/eps);

    // Restore original parameter value
    grp.setParam(param_ids[i], param);

  }

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::DerivUtils::computeDJnDp(LOCA::MultiContinuation::AbstractGroup& grp,
                   const std::vector<int>& paramIDs,
                   const NOX::Abstract::Vector& nullVector,
                   NOX::Abstract::MultiVector& result,
                   bool isValid) const
{
  std::string callingFunction =
    "LOCA::DerivUtils::computeDJnDp()";
  NOX::Abstract::Group::ReturnType status, finalStatus;

  // Views of Jn, d(Jn)/dp
  NOX::Abstract::Vector *Jn = &result[0];
  NOX::Abstract::Vector *dJndp = NULL;

  // Compute base residual
  if (!isValid) {
    finalStatus = grp.computeJacobian();
    globalData->locaErrorCheck->checkReturnType(finalStatus, callingFunction);

    status = grp.applyJacobian(nullVector, *Jn);
    finalStatus =
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                                 finalStatus,
                                 callingFunction);
  }
  else
    finalStatus = NOX::Abstract::Group::Ok;

  double param;
  double eps;

  // Loop over each parameter
  for (unsigned int i=0; i<paramIDs.size(); i++) {

    // Perturb single parameter in this group, and return perturbation, eps
    eps = perturbParam(grp, param, paramIDs[i]);

    // Fill perturbed Jn vector
    status = grp.computeJacobian();
    finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
                               callingFunction);

    dJndp = &result[i+1];
    status = grp.applyJacobian(nullVector, *dJndp);
    finalStatus =
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                                 finalStatus,
                                 callingFunction);

    // Difference perturbed and base vector
    dJndp->update(-1.0, *Jn, 1.0);
    dJndp->scale(1.0/eps);

    // Restore original parameter value
    grp.setParam(paramIDs[i], param);

  }

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::DerivUtils::computeDJnDxa(LOCA::MultiContinuation::AbstractGroup& grp,
                const NOX::Abstract::Vector& nullVector,
                const NOX::Abstract::MultiVector& aVector,
                NOX::Abstract::MultiVector& result) const
{
  std::string callingFunction =
    "LOCA::DerivUtils::computeDJnDxa()";
  NOX::Abstract::Group::ReturnType status, finalStatus;

  // Allocate base Jn vector and fill with J times n
  Teuchos::RCP<NOX::Abstract::Vector> baseJnVectorPtr =
    nullVector.clone(NOX::ShapeCopy);

  if (!grp.isJacobian()) {
    finalStatus = grp.computeJacobian();
    globalData->locaErrorCheck->checkReturnType(finalStatus, callingFunction);
  }
  else
    finalStatus = NOX::Abstract::Group::Ok;

  status = grp.applyJacobian(nullVector, *baseJnVectorPtr);
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
                               callingFunction);

  // Now that Jn is known, call other routine
  status = computeDJnDxa(grp, nullVector, aVector, *baseJnVectorPtr, result);
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
                               callingFunction);

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::DerivUtils::computeDJnDxa(LOCA::MultiContinuation::AbstractGroup& grp,
                const NOX::Abstract::Vector& nullVector,
                const NOX::Abstract::MultiVector& aVector,
                const NOX::Abstract::Vector& JnVector,
                NOX::Abstract::MultiVector& result) const
{
  std::string callingFunction =
    "LOCA::DerivUtils::computeDJnDxa()";
  NOX::Abstract::Group::ReturnType status;
  NOX::Abstract::Group::ReturnType finalStatus =
    NOX::Abstract::Group::NotDefined;

  // Copy original solution vector
  Teuchos::RCP<NOX::Abstract::Vector> Xvec =
    grp.getX().clone(NOX::DeepCopy);

  // Loop over each column of multivector
  for (int i=0; i<aVector.numVectors(); i++) {

    // Perturb solution vector in direction of aVector, return perturbation
    double eps = perturbXVec(grp, *Xvec, aVector[i]);

    // Fill perturbed Jn vector
    finalStatus = grp.computeJacobian();
    globalData->locaErrorCheck->checkReturnType(finalStatus, callingFunction);

    status = grp.applyJacobian(nullVector, result[i]);
    finalStatus =
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                                 finalStatus,
                                 callingFunction);

    // Difference perturbed and base vector
    result[i].update(-1.0, JnVector, 1.0);
    result[i].scale(1.0/eps);

  }

  // Restore original solution vector
  grp.setX(*Xvec);

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::DerivUtils::computeDwtJnDp(
                  LOCA::MultiContinuation::AbstractGroup& grp,
                  const std::vector<int>& paramIDs,
                  const NOX::Abstract::Vector& w,
                  const NOX::Abstract::Vector& nullVector,
                  NOX::Abstract::MultiVector::DenseMatrix& result,
                  bool isValid) const
{
  std::string callingFunction =
    "LOCA::DerivUtils::computeDwtJnDp()";
  NOX::Abstract::Group::ReturnType status, finalStatus;

  // Vector to store J*n
  Teuchos::RCP<NOX::Abstract::Vector> Jn =
    w.clone(NOX::ShapeCopy);
  double base_wtJn;

  // Compute base w^T*J*n
  if (!isValid) {

    // Compute J
    finalStatus = grp.computeJacobian();
    globalData->locaErrorCheck->checkReturnType(finalStatus, callingFunction);

    // Compute J*n
    status = grp.applyJacobian(nullVector, *Jn);
    finalStatus =
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                                 finalStatus,
                                 callingFunction);

    // Compute w^T*J*n
    base_wtJn = w.innerProduct(*Jn);
    result(0,0) = base_wtJn;
  }
  else {
    base_wtJn = result(0,0);
    finalStatus = NOX::Abstract::Group::Ok;
  }

  double param;
  double eps;
  double perturb_wtJn;

  // Loop over each parameter
  for (unsigned int i=0; i<paramIDs.size(); i++) {

    // Perturb single parameter in this group, and return perturbation, eps
    eps = perturbParam(grp, param, paramIDs[i]);

    // Compute perturbed J
    status = grp.computeJacobian();
    finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
                               callingFunction);
    // Compute perturbed J*n
    status = grp.applyJacobian(nullVector, *Jn);
    finalStatus =
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                                 finalStatus,
                                 callingFunction);

    // Compute perturbed w^T*J*n
    perturb_wtJn = w.innerProduct(*Jn);

    // Difference perturbed and base values
    result(0,i+1) = (perturb_wtJn - base_wtJn) / eps;

    // Restore original parameter value
    grp.setParam(paramIDs[i], param);

  }

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::DerivUtils::computeDwtJDp(LOCA::MultiContinuation::AbstractGroup& grp,
                const std::vector<int>& paramIDs,
                const NOX::Abstract::Vector& w,
                NOX::Abstract::MultiVector& result,
                bool isValid) const
{
  std::string callingFunction =
    "LOCA::DerivUtils::computeDwtJDp()";
  NOX::Abstract::Group::ReturnType status, finalStatus;

  // Views of w^T*J, d(w^T*J)/dp
  NOX::Abstract::Vector *wtJ = &result[0];
  NOX::Abstract::Vector *dwtJdp = NULL;

  // Compute base residual w^T*J
  if (!isValid) {
    finalStatus = grp.computeJacobian();
    globalData->locaErrorCheck->checkReturnType(finalStatus, callingFunction);

    status = grp.applyJacobianTranspose(w, *wtJ);
    finalStatus =
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                                 finalStatus,
                                 callingFunction);
  }
  else
    finalStatus = NOX::Abstract::Group::Ok;

  double param;
  double eps;

  // Loop over each parameter
  for (unsigned int i=0; i<paramIDs.size(); i++) {

    // Perturb single parameter in this group, and return perturbation, eps
    eps = perturbParam(grp, param, paramIDs[i]);

    // Fill perturbed w^T*J vector
    status = grp.computeJacobian();
    finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
                               callingFunction);

    dwtJdp = &result[i+1];
    status = grp.applyJacobianTranspose(w, *dwtJdp);
    finalStatus =
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                                 finalStatus,
                                 callingFunction);

    // Difference perturbed and base vector
    dwtJdp->update(-1.0, *wtJ, 1.0);
    dwtJdp->scale(1.0/eps);

    // Restore original parameter value
    grp.setParam(paramIDs[i], param);

  }

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::DerivUtils::computeDwtJnDx(LOCA::MultiContinuation::AbstractGroup& grp,
                 const NOX::Abstract::Vector& w,
                 const NOX::Abstract::Vector& nullVector,
                 NOX::Abstract::Vector& result) const
{
  std::string callingFunction =
    "LOCA::DerivUtils::computeDwtJnDx()";
  NOX::Abstract::Group::ReturnType status, finalStatus;

  // Vector to store w^T*J
  Teuchos::RCP<NOX::Abstract::Vector> wtJ =
    w.clone(NOX::ShapeCopy);

  // Compute base w^T*J
  finalStatus = grp.computeJacobian();
  globalData->locaErrorCheck->checkReturnType(finalStatus, callingFunction);

  status = grp.applyJacobianTranspose(w, *wtJ);
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
                               callingFunction);

  // Copy original solution vector
  Teuchos::RCP<NOX::Abstract::Vector> Xvec =
    grp.getX().clone(NOX::DeepCopy);

  // Perturb solution vector in direction of nullVector, return perturbation
  double eps = perturbXVec(grp, *Xvec, nullVector);

  // Fill perturbed w^T*J vector
  finalStatus = grp.computeJacobian();
  globalData->locaErrorCheck->checkReturnType(finalStatus, callingFunction);

  status = grp.applyJacobianTranspose(w, result);
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                               finalStatus,
                               callingFunction);

  // Difference perturbed and base vector
  result.update(-1.0, *wtJ, 1.0);
  result.scale(1.0/eps);

  // Restore original solution vector
  grp.setX(*Xvec);

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::DerivUtils::computeDwtJnDx(LOCA::MultiContinuation::AbstractGroup& grp,
                 const NOX::Abstract::MultiVector& w,
                 const NOX::Abstract::Vector& nullVector,
                 NOX::Abstract::MultiVector& result) const
{
  std::string callingFunction =
    "LOCA::DerivUtils::computeDwtJnDx()";
  NOX::Abstract::Group::ReturnType status, finalStatus;

  // Vector to store w^T*J
  Teuchos::RCP<NOX::Abstract::MultiVector> wtJ =
    w.clone(NOX::ShapeCopy);

  // Compute base w^T*J
  finalStatus = grp.computeJacobian();
  globalData->locaErrorCheck->checkReturnType(finalStatus, callingFunction);

  status = grp.applyJacobianTransposeMultiVector(w, *wtJ);
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
                               callingFunction);

  // Copy original solution vector
  Teuchos::RCP<NOX::Abstract::Vector> Xvec =
    grp.getX().clone(NOX::DeepCopy);

  // Perturb solution vector in direction of nullVector, return perturbation
  double eps = perturbXVec(grp, *Xvec, nullVector);

  // Fill perturbed w^T*J vector
  finalStatus = grp.computeJacobian();
  globalData->locaErrorCheck->checkReturnType(finalStatus, callingFunction);

  status = grp.applyJacobianTransposeMultiVector(w, result);
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                               finalStatus,
                               callingFunction);

  // Difference perturbed and base vector
  result.update(-1.0, *wtJ, 1.0);
  result.scale(1.0/eps);

  // Restore original solution vector
  grp.setX(*Xvec);

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::DerivUtils::computeDCeDp(LOCA::Hopf::MooreSpence::AbstractGroup& grp,
                   const std::vector<int>& paramIDs,
                   const NOX::Abstract::Vector& yVector,
                   const NOX::Abstract::Vector& zVector,
                   double w,
                   NOX::Abstract::MultiVector& result_real,
                   NOX::Abstract::MultiVector& result_imag,
                   bool isValid) const
{
  std::string callingFunction =
    "LOCA::DerivUtils::computeDCeDp()";
  NOX::Abstract::Group::ReturnType status, finalStatus;

  // Views of Ce, d(Ce)/dp
  NOX::Abstract::Vector& CeReal = result_real[0];
  NOX::Abstract::Vector& CeImag = result_imag[0];
  NOX::Abstract::Vector* dCedpReal;
  NOX::Abstract::Vector* dCedpImag;

  // Compute base residual
  if (!isValid) {
    finalStatus = grp.computeComplex(w);
    globalData->locaErrorCheck->checkReturnType(finalStatus, callingFunction);

    status = grp.applyComplex(yVector, zVector, CeReal, CeImag);
    finalStatus =
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                                 finalStatus,
                                 callingFunction);
  }
  else
    finalStatus = NOX::Abstract::Group::Ok;

  double param;
  double eps;

  // Loop over each parameter
  for (unsigned int i=0; i<paramIDs.size(); i++) {

    // Perturb single parameter in this group, and return perturbation, eps
    eps = perturbParam(grp, param, paramIDs[i]);

    // Fill perturbed Ce vectors
    status = grp.computeComplex(w);
    finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
                               callingFunction);

    dCedpReal = &(result_real[i+1]);
    dCedpImag = &(result_imag[i+1]);
    status = grp.applyComplex(yVector, zVector, *dCedpReal, *dCedpImag);
    finalStatus =
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                                 finalStatus,
                                 callingFunction);

    // Difference perturbed and base vectors
    dCedpReal->update(-1.0, CeReal, 1.0);
    dCedpReal->scale(1.0/eps);
    dCedpImag->update(-1.0, CeImag, 1.0);
    dCedpImag->scale(1.0/eps);

    // Restore original parameter value
    grp.setParam(paramIDs[i], param);

  }

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::DerivUtils::computeDCeDxa(
                LOCA::Hopf::MooreSpence::AbstractGroup& grp,
                const NOX::Abstract::Vector& yVector,
                const NOX::Abstract::Vector& zVector,
                double w,
                const NOX::Abstract::MultiVector& aVector,
                NOX::Abstract::MultiVector& result_real,
                NOX::Abstract::MultiVector& result_imag) const
{
  std::string callingFunction =
    "LOCA::DerivUtils::computeDCeDxa()";
  NOX::Abstract::Group::ReturnType status, finalStatus;

  // Allocate base Ce
  Teuchos::RCP<NOX::Abstract::Vector> CeReal =
    yVector.clone(NOX::ShapeCopy);
  Teuchos::RCP<NOX::Abstract::Vector> CeImag =
    zVector.clone(NOX::ShapeCopy);

  // Compute base Ce
  finalStatus = grp.computeComplex(w);
  globalData->locaErrorCheck->checkReturnType(finalStatus, callingFunction);

  status = grp.applyComplex(yVector, zVector, *CeReal, *CeImag);
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                               finalStatus,
                               callingFunction);

  // Now that Ce is known, call other routine
  status =
    computeDCeDxa(grp, yVector, zVector, w, aVector, *CeReal, *CeImag,
          result_real, result_imag);
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
                         callingFunction);

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::DerivUtils::computeDCeDxa(
                LOCA::Hopf::MooreSpence::AbstractGroup& grp,
                const NOX::Abstract::Vector& yVector,
                const NOX::Abstract::Vector& zVector,
                double w,
                const NOX::Abstract::MultiVector& aVector,
                const NOX::Abstract::Vector& Ce_real,
                const NOX::Abstract::Vector& Ce_imag,
                NOX::Abstract::MultiVector& result_real,
                NOX::Abstract::MultiVector& result_imag) const
{
  std::string callingFunction =
    "LOCA::DerivUtils::computeDCeDxa()";
  NOX::Abstract::Group::ReturnType status, finalStatus =
    NOX::Abstract::Group::Ok;

  // Copy original solution vector
  Teuchos::RCP<NOX::Abstract::Vector> Xvec =
    grp.getX().clone(NOX::DeepCopy);

  // Loop over each column of multivector
  for (int i=0; i<aVector.numVectors(); i++) {

    // Perturb solution vector in direction of aVector, return perturbation
    double eps = perturbXVec(grp, *Xvec, aVector[i]);

    // Compute perturbed Ce vectors
    status = grp.computeComplex(w);
    finalStatus =
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                                 finalStatus,
                                 callingFunction);

    status =
      grp.applyComplex(yVector, zVector, result_real[i], result_imag[i]);
    finalStatus =
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                                 finalStatus,
                                 callingFunction);

    // Difference perturbed and base vector and return approximate derivative
    result_real[i].update(-1.0, Ce_real, 1.0); result_real[i].scale(1.0/eps);
    result_imag[i].update(-1.0, Ce_imag, 1.0); result_imag[i].scale(1.0/eps);

  }

  // Restore original solution vector
  grp.setX(*Xvec);

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::DerivUtils::computeDwtCeDp(
              LOCA::Hopf::MinimallyAugmented::AbstractGroup& grp,
              const std::vector<int>& paramIDs,
              const NOX::Abstract::Vector& w1,
              const NOX::Abstract::Vector& w2,
              const NOX::Abstract::Vector& yVector,
              const NOX::Abstract::Vector& zVector,
              double omega,
              NOX::Abstract::MultiVector::DenseMatrix& result_real,
              NOX::Abstract::MultiVector::DenseMatrix& result_imag,
              bool isValid) const
{
  std::string callingFunction =
    "LOCA::DerivUtils::computeDwtCeDp()";
  NOX::Abstract::Group::ReturnType status, finalStatus;

  // Views of Ce, d(Ce)/dp
  Teuchos::RCP<NOX::Abstract::Vector> CeReal =
    w1.clone(NOX::ShapeCopy);
  Teuchos::RCP<NOX::Abstract::Vector> CeImag =
    w2.clone(NOX::ShapeCopy);

  // Compute base w^T*C*e
  if (!isValid) {
    finalStatus = grp.computeComplex(omega);
    globalData->locaErrorCheck->checkReturnType(finalStatus, callingFunction);

    status = grp.applyComplex(yVector, zVector, *CeReal, *CeImag);
    finalStatus =
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                                 finalStatus,
                                 callingFunction);

    // Conjugate transpose
    result_real(0,0) = w1.innerProduct(*CeReal) + w2.innerProduct(*CeImag);
    result_imag(0,0) = w1.innerProduct(*CeImag) - w2.innerProduct(*CeReal);
  }
  else
    finalStatus = NOX::Abstract::Group::Ok;

  double param;
  double eps;

  // Loop over each parameter
  for (unsigned int i=0; i<paramIDs.size(); i++) {

    // Perturb single parameter in this group, and return perturbation, eps
    eps = perturbParam(grp, param, paramIDs[i]);

    // Fill perturbed Ce vectors
    status = grp.computeComplex(omega);
    finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
                               callingFunction);

    status = grp.applyComplex(yVector, zVector, *CeReal, *CeImag);
    finalStatus =
      globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                                 finalStatus,
                                 callingFunction);

    // Difference perturbed and base values
    // Conjugate transpose
    result_real(0,i+1) = (w1.innerProduct(*CeReal) + w2.innerProduct(*CeImag) -
              result_real(0,0)) / eps;
    result_imag(0,i+1) = (w1.innerProduct(*CeImag) - w2.innerProduct(*CeReal) -
              result_imag(0,0)) /eps;

    // Restore original parameter value
    grp.setParam(paramIDs[i], param);

  }

  return finalStatus;
}

NOX::Abstract::Group::ReturnType
LOCA::DerivUtils::computeDwtCeDx(
                LOCA::Hopf::MinimallyAugmented::AbstractGroup& grp,
                const NOX::Abstract::Vector& w1,
                const NOX::Abstract::Vector& w2,
                const NOX::Abstract::Vector& yVector,
                const NOX::Abstract::Vector& zVector,
                double omega,
                NOX::Abstract::Vector& result_real,
                NOX::Abstract::Vector& result_imag) const
{
  std::string callingFunction =
    "LOCA::DerivUtils::computeDwtCeDxa()";
  NOX::Abstract::Group::ReturnType status, finalStatus =
    NOX::Abstract::Group::Ok;

  // Vectors to store w^T*C
  Teuchos::RCP<NOX::Abstract::Vector> wtC_real =
    w1.clone(NOX::ShapeCopy);
  Teuchos::RCP<NOX::Abstract::Vector> wtC_imag =
    w2.clone(NOX::ShapeCopy);

  // Compute base w^T*C
  finalStatus = grp.computeComplex(omega);
  globalData->locaErrorCheck->checkReturnType(finalStatus, callingFunction);

  status = grp.applyComplexTranspose(w1, w2, *wtC_real, *wtC_imag);
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
                               callingFunction);

  // Copy original solution vector
  Teuchos::RCP<NOX::Abstract::Vector> Xvec =
    grp.getX().clone(NOX::DeepCopy);

  // Perturb solution vector in direction of yVector, return perturbation
  double eps = perturbXVec(grp, *Xvec, yVector);

  // Compute perturbed wtC vectors
  status = grp.computeComplex(omega);
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                               finalStatus,
                               callingFunction);

  status =
    grp.applyComplexTranspose(w1, w2, result_real, result_imag);
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                               finalStatus,
                               callingFunction);

  // Difference perturbed and base vector and return approximate derivative
  result_real.update(-1.0, *wtC_real, 1.0); result_real.scale(1.0/eps);
  result_imag.update(-1.0, *wtC_imag, 1.0); result_imag.scale(1.0/eps);

  // Take conjugate
  result_imag.scale(-1.0);

  // Restore original solution vector
  grp.setX(*Xvec);

  // Perturb solution vector in direction of zVector, return perturbation
  eps = perturbXVec(grp, *Xvec, zVector);

  // Compute perturbed wtC vectors
  status = grp.computeComplex(omega);
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                               finalStatus,
                               callingFunction);

  Teuchos::RCP<NOX::Abstract::Vector> tmp_r =
    result_real.clone(NOX::ShapeCopy);
  Teuchos::RCP<NOX::Abstract::Vector> tmp_i =
    result_imag.clone(NOX::ShapeCopy);
  status =
    grp.applyComplexTranspose(w1, w2, *tmp_r, *tmp_i);
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                               finalStatus,
                               callingFunction);

  // Difference perturbed and base vector and return approximate derivative
  tmp_r->update(-1.0, *wtC_real, 1.0); tmp_r->scale(1.0/eps);
  tmp_i->update(-1.0, *wtC_imag, 1.0); tmp_i->scale(1.0/eps);

  // Take conjugate
  tmp_i->scale(-1.0);

  result_real.update(-1.0, *tmp_i, 1.0);
  result_imag.update( 1.0, *tmp_r, 1.0);

  // Restore original solution vector
  grp.setX(*Xvec);

  return finalStatus;
}

//
// Protected methods start here.
//

double
LOCA::DerivUtils::perturbParam(LOCA::MultiContinuation::AbstractGroup& grp,
                   double& paramOrig,
                   int param_id) const
{
  paramOrig = grp.getParam(param_id);

  // Find perturbation size and perturb parameter
  double eps = epsScalar(paramOrig);
  double param = paramOrig + eps;

  // Copy this perturbed parameter vector into group
  grp.setParam(param_id, param);

  // Return perturbation size
  return eps;
}

double
LOCA::DerivUtils::perturbXVec(LOCA::MultiContinuation::AbstractGroup& grp,
                  const NOX::Abstract::Vector& xVector,
                  const NOX::Abstract::Vector& aVector) const
{
  // Allocate tempertory xVector
  Teuchos::RCP<NOX::Abstract::Vector> tmpXVecPtr =
    xVector.clone(NOX::DeepCopy);

  // Get perturbation size for directional derivative
  double eps = epsVector(*tmpXVecPtr, aVector);

  // Perturb temp vector and copy into group's x vector
  grp.setX(tmpXVecPtr->update(eps, aVector, 1.0));

  // Return perturbation size
  return eps;
}

double
LOCA::DerivUtils::epsScalar(double p) const
{
   return perturb * (perturb + fabs(p));
}

double
LOCA::DerivUtils::epsVector(const NOX::Abstract::Vector& xVector,
                const NOX::Abstract::Vector& aVector) const
{
   return perturb * (perturb + xVector.norm(NOX::Abstract::Vector::TwoNorm)
                  / (aVector.norm(NOX::Abstract::Vector::TwoNorm) + perturb));
}
