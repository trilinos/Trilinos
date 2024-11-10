// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_AnasaziOperator_Cayley2Matrix.H"
#include "Teuchos_ParameterList.hpp"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"

LOCA::AnasaziOperator::Cayley2Matrix::Cayley2Matrix(
    const Teuchos::RCP<LOCA::GlobalData>& global_data,
    const Teuchos::RCP<LOCA::Parameter::SublistParser>& /* topParams */,
    const Teuchos::RCP<Teuchos::ParameterList>& eigenParams_,
    const Teuchos::RCP<Teuchos::ParameterList>& solverParams_,
    const Teuchos::RCP<LOCA::TimeDependent::AbstractGroup>& grp_)
  : globalData(global_data),
    myLabel("Cayley2Matrix Transformation"),
    eigenParams(eigenParams_),
    solverParams(solverParams_),
    grp(grp_),
    tmp_r(),
    tmp_i(),
    sigma(0.0),
    mu(0.0)
{
  sigma = eigenParams->get("Cayley Pole",0.0);
  mu = eigenParams->get("Cayley Zero",0.0);
}

LOCA::AnasaziOperator::Cayley2Matrix::~Cayley2Matrix()
{
}

const std::string&
LOCA::AnasaziOperator::Cayley2Matrix::label() const
{
  return myLabel;
}

void
LOCA::AnasaziOperator::Cayley2Matrix::apply(const NOX::Abstract::MultiVector& input,
                     NOX::Abstract::MultiVector& output) const
{
  std::string callingFunction =
    "LOCA::AnasaziOperator::Cayley2Matrix::apply()";

  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // Allocate temporary vector
  if (tmp_r == Teuchos::null || tmp_r->numVectors() != input.numVectors())
    tmp_r = input.clone(NOX::ShapeCopy);

  // Compute J-mu*M -- moved to preProcessSeedVector

  // Compute (J-mu*M)*input
  status = grp->applySecondShiftedMatrixMultiVector(input, *tmp_r);
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                               finalStatus,
                               callingFunction);

  // Compute J-sigma*M -- moved to preProcessSeedVector

  // Solve (J-sigma*M)*output = (J-mu*M)*input
  status = grp->applyShiftedMatrixInverseMultiVector(*solverParams, *tmp_r,
                             output);

  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                               finalStatus,
                               callingFunction);
}

void
LOCA::AnasaziOperator::Cayley2Matrix::preProcessSeedVector(NOX::Abstract::MultiVector& ivec)
{
  // Changes random seed vector ivec:   ivec = (J - sigma*M)^{-1}*M*ivec
  std::string callingFunction =
    "LOCA::AnasaziOperator::Cayley2Matrix::preProcessSeedVector()";

  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // Allocate temporary vector
  if (tmp_r == Teuchos::null || tmp_r->numVectors() != ivec.numVectors())
    tmp_r = ivec.clone(NOX::ShapeCopy);

  // Compute M
  status = grp->computeShiftedMatrix(0.0, 1.0);
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                               finalStatus,
                               callingFunction);

  // Compute M*ivec
  status = grp->applyShiftedMatrixMultiVector(ivec, *tmp_r);
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                               finalStatus,
                               callingFunction);

  // Compute J-sigma*M
  status = grp->computeShiftedMatrix(1.0, -sigma);
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                               finalStatus,
                               callingFunction);

  // Solve (J-sigma*M)*output = (M)*ivec
  status = grp->applyShiftedMatrixInverseMultiVector(*solverParams, *tmp_r,
                             ivec);

  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status,
                               finalStatus,
                               callingFunction);

  // Compute J-mu*M
  status = grp->computeSecondShiftedMatrix(1.0, -mu);
}

void
LOCA::AnasaziOperator::Cayley2Matrix::beginPostProcessing()
{
  // Make sure mass matrix is up-to-date
  grp->computeShiftedMatrix(0.0, 1.0);

  // Make sure Jacobian is up-to-date
  grp->computeJacobian();
}

void
LOCA::AnasaziOperator::Cayley2Matrix::transformEigenvalue(double& ev_r,
                           double& ev_i) const
{
  // compute inverse of eigenvalue, then shift
  double mag = (1.0 - ev_r)*(1.0 - ev_r) + ev_i*ev_i;
  ev_r = (sigma*(ev_r*ev_r + ev_i*ev_i) - (sigma+mu)*ev_r + mu) / mag;
  ev_i = (mu-sigma)*ev_i/mag;
}

NOX::Abstract::Group::ReturnType
LOCA::AnasaziOperator::Cayley2Matrix::rayleighQuotient(
                         NOX::Abstract::Vector& evec_r,
                     NOX::Abstract::Vector& evec_i,
                     double& rq_r, double& rq_i) const
{
  std::string callingFunction =
    "LOCA::AnasaziOperator::Cayley2Matrix::rayleighQuotient()";

  // Allocate temporary vectors
  if (tmp_r == Teuchos::null)
    tmp_r = evec_r.createMultiVector(1, NOX::ShapeCopy);
  if (tmp_i == Teuchos::null)
    tmp_i = evec_i.createMultiVector(1, NOX::ShapeCopy);

  NOX::Abstract::Group::ReturnType finalStatus = NOX::Abstract::Group::Ok;
  NOX::Abstract::Group::ReturnType status;

  // Compute z^h J z
  status = grp->applyJacobian(evec_r, (*tmp_r)[0]);
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
                               callingFunction);

  status = grp->applyJacobian(evec_i, (*tmp_i)[0]);
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
                               callingFunction);

  rq_r = evec_r.innerProduct((*tmp_r)[0]) + evec_i.innerProduct((*tmp_i)[0]);
  rq_i = evec_r.innerProduct((*tmp_i)[0]) - evec_i.innerProduct((*tmp_r)[0]);

  // Compute z^h M z
  status = grp->applyShiftedMatrix(evec_r, (*tmp_r)[0]);
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
                               callingFunction);

  status = grp->applyShiftedMatrix(evec_i, (*tmp_i)[0]);
  finalStatus =
    globalData->locaErrorCheck->combineAndCheckReturnTypes(status, finalStatus,
                               callingFunction);

  double m_r =
    evec_r.innerProduct((*tmp_r)[0]) + evec_i.innerProduct((*tmp_i)[0]);
  double m_i =
    evec_r.innerProduct((*tmp_i)[0]) - evec_i.innerProduct((*tmp_r)[0]);
  double m = m_r*m_r + m_i*m_i;

  // Compute z^h J z / z^h M z
  double tmp = (rq_r*m_r + rq_i*m_i) / m;
  rq_i = (rq_i*m_r - rq_r*m_i) / m;
  rq_r = tmp;

  if (eigenParams->get("Normalize Eigenvectors with Mass Matrix",false)) {
    double scl = 1.0 / sqrt(sqrt(m));
    evec_r.scale(scl);
    evec_i.scale(scl);
  }

  return finalStatus;
}
