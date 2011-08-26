// @HEADER
// ************************************************************************
// 
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
// 
// ************************************************************************
// @HEADER

#include "Piro_Epetra_SensitivityOperator.hpp"
#include "Teuchos_TestForException.hpp"
#include "Epetra_Map.h"
#include "Epetra_LocalMap.h"
#include "Epetra_Vector.h"
#include "NOX_Epetra_MultiVector.H"

#include "Stokhos_MatrixFreeOperator.hpp"
#include "Stokhos_EpetraMultiVectorOperator.hpp"
#include "Stokhos_EpetraOperatorOrthogPoly.hpp"

Piro::Epetra::SensitivityOperator::
SensitivityOperator(
  const Teuchos::RCP<const Epetra_Map>& g_map_,
  const Teuchos::RCP<const Epetra_Map>& p_map_) : 
  label("Piro Epetra Sensitivity Operator"),
  g_map(g_map_),
  p_map(p_map_),
  useTranspose(false),
  grp()
{
}

void
Piro::Epetra::SensitivityOperator::
setup(
  const EpetraExt::ModelEvaluator::Derivative& dfdp_,
  const EpetraExt::ModelEvaluator::Derivative& dgdx_,
  const EpetraExt::ModelEvaluator::Derivative& dgdp_,
  Teuchos::RCP<Teuchos::ParameterList>& piroParams_,
  const Teuchos::RCP<NOX::Epetra::Group>& grp_,
  const Teuchos::RCP<LOCA::Epetra::TransposeLinearSystem::AbstractStrategy>& tls_strategy_)
{
  dfdp = dfdp_;
  dgdx = dgdx_;
  dgdp = dgdp_;
  piroParams = piroParams_;
  grp = grp_;
  tls_strategy = tls_strategy_;
}

Piro::Epetra::SensitivityOperator::
~SensitivityOperator()
{
}

int 
Piro::Epetra::SensitivityOperator::
SetUseTranspose(bool UseTranspose) 
{
  if (useTranspose == false && UseTranspose == true) {
    Teuchos::RCP<NOX::Epetra::LinearSystem> linSys = grp->getLinearSystem();
    Teuchos::RCP<Epetra_Operator> jac = linSys->getJacobianOperator();
    jac->SetUseTranspose(false);
    linSys->setJacobianOperatorForSolve(jac);
    linSys->destroyPreconditioner();
  }

  useTranspose = UseTranspose;

  if (dfdp.getLinearOp() != Teuchos::null)
    dfdp.getLinearOp()->SetUseTranspose(UseTranspose);
  if (dgdx.getLinearOp() != Teuchos::null)
    dgdx.getLinearOp()->SetUseTranspose(UseTranspose);
  if (dgdp.getLinearOp() != Teuchos::null)
    dgdp.getLinearOp()->SetUseTranspose(UseTranspose);

  if (useTranspose) {
    tls_strategy->createJacobianTranspose();
    const NOX::Epetra::Vector& x_nox = 
      dynamic_cast<const NOX::Epetra::Vector&>(grp->getX());
    tls_strategy->createTransposePreconditioner(x_nox, *piroParams);
  }

  return 0;
}

int 
Piro::Epetra::SensitivityOperator::
Apply(const Epetra_MultiVector& Input, Epetra_MultiVector& Result) const
{
  int ret;
  if (useTranspose)
    ret = applyAdjoint(Input, Result);
  else
    ret = applyForward(Input, Result);
  return ret;
}

int 
Piro::Epetra::SensitivityOperator::
ApplyInverse(const Epetra_MultiVector& Input, Epetra_MultiVector& Result) const
{
  throw "SensitivityOperator::ApplyInverse not defined!";
  return -1;
}

double 
Piro::Epetra::SensitivityOperator::
NormInf() const
{
  return 1.0;
}


const char* 
Piro::Epetra::SensitivityOperator::
Label () const
{
  return const_cast<char*>(label.c_str());
}
  
bool 
Piro::Epetra::SensitivityOperator::
UseTranspose() const
{
  return useTranspose;
}

bool 
Piro::Epetra::SensitivityOperator::
HasNormInf() const
{
  return false;
}

const Epetra_Comm & 
Piro::Epetra::SensitivityOperator::
Comm() const
{
  return p_map->Comm();
}
const Epetra_Map& 
Piro::Epetra::SensitivityOperator::
OperatorDomainMap() const
{
  if (useTranspose)
    return *g_map;
  return *p_map;
}

const Epetra_Map& 
Piro::Epetra::SensitivityOperator::
OperatorRangeMap() const
{
  if (useTranspose)
    return *p_map;
  return *g_map;
}

int
Piro::Epetra::SensitivityOperator::
applyForward(const Epetra_MultiVector& Input, 
	     Epetra_MultiVector& Result) const
{
  int m = Input.NumVectors();

  // Calculate Result = dg/dp*Input
  Teuchos::RCP<Epetra_Operator> dgdp_lo = dgdp.getLinearOp();
  Teuchos::RCP<Epetra_MultiVector> dgdp_mv = dgdp.getMultiVector();
  if (dgdp_lo != Teuchos::null) {
    dgdp_lo->Apply(Input, Result);
  }
  else if (dgdp_mv != Teuchos::null) {
    EpetraExt::ModelEvaluator::EDerivativeMultiVectorOrientation dgdp_orient =
      dgdp.getMultiVectorOrientation();
    if (dgdp_orient == EpetraExt::ModelEvaluator::DERIV_MV_BY_COL)
      Result.Multiply('N', 'N', 1.0, *dgdp_mv, Input, 0.0);
    else
      Result.Multiply('T', 'N', 1.0, *dgdp_mv, Input, 0.0);
  }
  else 
    TEST_FOR_EXCEPTION(
      true, std::logic_error,
      std::endl <<"Piro::Epetra::SensitivityOperator::applyForward():  " <<
      "dgdp is NULL" << std::endl);

  // Compute w = df/dp*Input
  Teuchos::RCP<Epetra_Operator> dfdp_lo = dfdp.getLinearOp();
  Teuchos::RCP<Epetra_MultiVector> dfdp_mv = dfdp.getMultiVector();
  Teuchos::RCP<Epetra_MultiVector> w;
  if (dfdp_lo != Teuchos::null) {
    w = Teuchos::rcp(new Epetra_MultiVector(dfdp_lo->OperatorRangeMap(), m));
    dfdp_lo->Apply(Input, *w);
  }
  else if (dfdp_mv != Teuchos::null) {
    EpetraExt::ModelEvaluator::EDerivativeMultiVectorOrientation dfdp_orient =
      dfdp.getMultiVectorOrientation();
    if (dfdp_orient == EpetraExt::ModelEvaluator::DERIV_MV_BY_COL) {
      w = Teuchos::rcp(new Epetra_MultiVector(dfdp_mv->Map(), m));
      w->Multiply('N', 'N', 1.0, *dfdp_mv, Input, 0.0);
    }
    else {
      Epetra_LocalMap f_map(dfdp_mv->NumVectors(), 0, dfdp_mv->Comm());
      w = Teuchos::rcp(new Epetra_MultiVector(f_map, m));
      w->Multiply('T', 'N', 1.0, *dfdp_mv, Input, 0.0);
    }
  }
  else 
    TEST_FOR_EXCEPTION(
      true, std::logic_error,
      std::endl <<"Piro::Epetra::SensitivityOperator::applyForward():  " <<
      "dfdp is NULL" << std::endl);

  // Calculate y = -J^{-1}*w
  Teuchos::RCP<Epetra_MultiVector> y = 
    Teuchos::rcp(new Epetra_MultiVector(w->Map(), m));
  NOX::Epetra::MultiVector w_nox(w, NOX::DeepCopy,  
				 NOX::Epetra::MultiVector::CreateView);
  NOX::Epetra::MultiVector y_nox(y, NOX::DeepCopy,  
				 NOX::Epetra::MultiVector::CreateView);
  grp->applyJacobianInverseMultiVector(*piroParams, w_nox, y_nox);

  // Calculate Result = Result - dg/dx*y
  Teuchos::RCP<Epetra_Operator> dgdx_lo = dgdx.getLinearOp();
  Teuchos::RCP<Epetra_MultiVector> dgdx_mv = dgdx.getMultiVector();
  if (dgdx_lo != Teuchos::null) {
    Teuchos::RCP<Epetra_MultiVector> z = 
      Teuchos::rcp(new Epetra_MultiVector(*g_map, m));
    dgdx_lo->Apply(*y, *z);
    Result.Update(-1.0, *z, 1.0);
  }
  else if (dgdx_mv != Teuchos::null) {
    EpetraExt::ModelEvaluator::EDerivativeMultiVectorOrientation dgdx_orient =
      dgdx.getMultiVectorOrientation();
    if (dgdx_orient == EpetraExt::ModelEvaluator::DERIV_MV_BY_COL)
      Result.Multiply('N', 'N', -1.0, *dgdx_mv, *y, 1.0);
    else
      Result.Multiply('T', 'N', -1.0, *dgdx_mv, *y, 1.0);
  }
  else 
    TEST_FOR_EXCEPTION(
      true, std::logic_error,
      std::endl <<"Piro::Epetra::SensitivityOperator::applyForward():  " <<
      "dgdx is NULL" << std::endl);

  return 0;
}

int
Piro::Epetra::SensitivityOperator::
applyAdjoint(const Epetra_MultiVector& Input, 
	     Epetra_MultiVector& Result) const
{
  int m = Input.NumVectors();

  // Note:  Here we are assuming SetUseTranspose(true) has been called on
  // each operator

  // Calculate Result = dg/dp^T*Input
  Teuchos::RCP<Epetra_Operator> dgdp_lo = dgdp.getLinearOp();
  Teuchos::RCP<Epetra_MultiVector> dgdp_mv = dgdp.getMultiVector();
  if (dgdp_lo != Teuchos::null) {
    dgdp_lo->Apply(Input, Result);
  }
  else if (dgdp_mv != Teuchos::null) {
    EpetraExt::ModelEvaluator::EDerivativeMultiVectorOrientation dgdp_orient =
      dgdp.getMultiVectorOrientation();
    if (dgdp_orient == EpetraExt::ModelEvaluator::DERIV_MV_BY_COL)
      Result.Multiply('T', 'N', 1.0, *dgdp_mv, Input, 0.0);
    else
      Result.Multiply('N', 'N', 1.0, *dgdp_mv, Input, 0.0);
  }
  else 
    TEST_FOR_EXCEPTION(
      true, std::logic_error,
      std::endl <<"Piro::Epetra::SensitivityOperator::applyForward():  " <<
      "dgdp is NULL" << std::endl);

  // Compute w = dg/dx^T*Input
  Teuchos::RCP<Epetra_Operator> dgdx_lo = dgdx.getLinearOp();
  Teuchos::RCP<Epetra_MultiVector> dgdx_mv = dgdx.getMultiVector();
  Teuchos::RCP<Epetra_MultiVector> w;
  if (dgdx_lo != Teuchos::null) {
     w = Teuchos::rcp(new Epetra_MultiVector(dgdx_lo->OperatorRangeMap(), m));
    dgdx_lo->Apply(Input, *w);
  }
  else if (dgdx_mv != Teuchos::null) {
    EpetraExt::ModelEvaluator::EDerivativeMultiVectorOrientation dgdx_orient =
      dgdx.getMultiVectorOrientation();
    if (dgdx_orient == EpetraExt::ModelEvaluator::DERIV_MV_BY_COL) {
      Epetra_LocalMap x_map(dgdx_mv->NumVectors(), 0, dgdx_mv->Comm());
      w = Teuchos::rcp(new Epetra_MultiVector(x_map, m));
      w->Multiply('T', 'N', 1.0, *dgdx_mv, Input, 0.0);
    }
    else {
      w = Teuchos::rcp(new Epetra_MultiVector(dgdx_mv->Map(), m));
      w->Multiply('N', 'N', 1.0, *dgdx_mv, Input, 0.0);
    }
  }
  else 
    TEST_FOR_EXCEPTION(
      true, std::logic_error,
      std::endl <<"Piro::Epetra::SensitivityOperator::applyForward():  " <<
      "dgdx is NULL" << std::endl);

  // Calculate y = -J^{-T}*w
  Teuchos::RCP<Epetra_MultiVector> y = 
    Teuchos::rcp(new Epetra_MultiVector(w->Map(), m));
  for (int col=0; col<m; col++) {
    Teuchos::RCP<Epetra_Vector> ww = Teuchos::rcp((*w)(col),false);
    Teuchos::RCP<Epetra_Vector> yy = Teuchos::rcp((*y)(col),false);
    NOX::Epetra::Vector w_nox(ww, NOX::Epetra::Vector::CreateView, 
			      NOX::DeepCopy);
    NOX::Epetra::Vector y_nox(yy, NOX::Epetra::Vector::CreateView, 
			      NOX::DeepCopy);
    tls_strategy->applyJacobianTransposeInverse(*piroParams, w_nox, y_nox);
  }

  // Compute Result = Result - df/dp^T*y
  Teuchos::RCP<Epetra_Operator> dfdp_lo = dfdp.getLinearOp();
  Teuchos::RCP<Epetra_MultiVector> dfdp_mv = dfdp.getMultiVector();
  if (dfdp_lo != Teuchos::null) {
    Teuchos::RCP<Epetra_MultiVector> z = 
      Teuchos::rcp(new Epetra_MultiVector(*p_map, m));
    dfdp_lo->Apply(*y, *z);
    Result.Update(-1.0, *z, 1.0);
  }
  else if (dfdp_mv != Teuchos::null) {
    EpetraExt::ModelEvaluator::EDerivativeMultiVectorOrientation dfdp_orient =
      dfdp.getMultiVectorOrientation();
    if (dfdp_orient == EpetraExt::ModelEvaluator::DERIV_MV_BY_COL)
      Result.Multiply('T', 'N', -1.0, *dfdp_mv, *y, 1.0);
    else 
      Result.Multiply('N', 'N', -1.0, *dfdp_mv, *y, 1.0);
  }
  else 
    TEST_FOR_EXCEPTION(
      true, std::logic_error,
      std::endl <<"Piro::Epetra::SensitivityOperator::applyForward():  " <<
      "dfdp is NULL" << std::endl);

  return 0;
}
