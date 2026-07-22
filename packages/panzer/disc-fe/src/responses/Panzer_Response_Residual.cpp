// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Panzer_Response_Residual.hpp"

#include "Panzer_ThyraObjFactory.hpp"
#include "Panzer_ThyraObjContainer.hpp"

#include "Thyra_VectorSpaceBase.hpp"

namespace panzer {

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

Teuchos::RCP<Thyra::VectorBase<panzer::Traits::RealType> > 
Response_Residual<panzer::Traits::Residual>::
getGhostedResidual() const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  typedef ThyraObjContainer<panzer::Traits::RealType> TOC;

  // if already computed, uses that ghosted vector
  if(ghostedResidual_!=Teuchos::null)
    return ghostedResidual_;

  // otherwise, allocate a new ghosted vector
  RCP<LinearObjContainer> loc = linObjFactory_->buildGhostedLinearObjContainer();
  linObjFactory_->initializeGhostedContainer(LinearObjContainer::F,*loc);
 
  RCP<TOC> th_loc = rcp_dynamic_cast<TOC>(loc);
  return th_loc->get_f_th();
}

Teuchos::RCP<Thyra::VectorBase<panzer::Traits::RealType> > 
Response_Residual<panzer::Traits::Residual>::
getResidual() const
{
  return residual_;
}

void 
Response_Residual<panzer::Traits::Residual>::
setResidual(const Teuchos::RCP<Thyra::VectorBase<panzer::Traits::RealType> > & res)
{
  residual_ = res;
}

Teuchos::RCP<Thyra::VectorBase<panzer::Traits::RealType> > 
Response_Residual<panzer::Traits::Residual>::
allocateResidualVector() const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  typedef ThyraObjFactory<panzer::Traits::RealType> ObjFactory;

  RCP<const ObjFactory> objFactory = rcp_dynamic_cast<const ObjFactory>(linObjFactory_);
  return Thyra::createMember(objFactory->getThyraRangeSpace());
}

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

Teuchos::RCP<Thyra::LinearOpBase<panzer::Traits::RealType> > 
Response_Residual<panzer::Traits::Jacobian>::
getGhostedJacobian() const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  typedef ThyraObjContainer<panzer::Traits::RealType> TOC;

  // if already computed, uses that ghosted vector
  if(ghostedJacobian_!=Teuchos::null)
    return ghostedJacobian_;

  // otherwise, allocate a new ghosted vector
  RCP<LinearObjContainer> loc = linObjFactory_->buildGhostedLinearObjContainer();
  linObjFactory_->initializeGhostedContainer(LinearObjContainer::Mat,*loc);
 
  RCP<TOC> th_loc = rcp_dynamic_cast<TOC>(loc);
  return th_loc->get_A_th();
}

Teuchos::RCP<Thyra::LinearOpBase<panzer::Traits::RealType> > 
Response_Residual<panzer::Traits::Jacobian>::
getJacobian() const
{
  return jacobian_;
}

void 
Response_Residual<panzer::Traits::Jacobian>::
setJacobian(const Teuchos::RCP<Thyra::LinearOpBase<panzer::Traits::RealType> > & jac)
{
  jacobian_ = jac;
}

Teuchos::RCP<Thyra::LinearOpBase<panzer::Traits::RealType> > 
Response_Residual<panzer::Traits::Jacobian>::
allocateJacobian() const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  typedef ThyraObjFactory<panzer::Traits::RealType> ObjFactory;

  RCP<const ObjFactory> objFactory = rcp_dynamic_cast<const ObjFactory>(linObjFactory_);
  return objFactory->getThyraMatrix();
}

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

#ifdef Panzer_BUILD_HESSIAN_SUPPORT

Teuchos::RCP<Thyra::LinearOpBase<panzer::Traits::RealType> > 
Response_Residual<panzer::Traits::Hessian>::
getGhostedHessian() const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  typedef ThyraObjContainer<panzer::Traits::RealType> TOC;

  // if already computed, uses that ghosted vector
  if(ghostedHessian_!=Teuchos::null)
    return ghostedHessian_;

  // otherwise, allocate a new ghosted vector
  RCP<LinearObjContainer> loc = linObjFactory_->buildGhostedLinearObjContainer();
  linObjFactory_->initializeGhostedContainer(LinearObjContainer::Mat,*loc);
 
  RCP<TOC> th_loc = rcp_dynamic_cast<TOC>(loc);
  return th_loc->get_A_th();
}

Teuchos::RCP<Thyra::LinearOpBase<panzer::Traits::RealType> > 
Response_Residual<panzer::Traits::Hessian>::
getHessian() const
{
  return hessian_;
}

void 
Response_Residual<panzer::Traits::Hessian>::
setHessian(const Teuchos::RCP<Thyra::LinearOpBase<panzer::Traits::RealType> > & jac)
{
  hessian_ = jac;
}

Teuchos::RCP<Thyra::LinearOpBase<panzer::Traits::RealType> > 
Response_Residual<panzer::Traits::Hessian>::
allocateHessian() const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;

  typedef ThyraObjFactory<panzer::Traits::RealType> ObjFactory;

  RCP<const ObjFactory> objFactory = rcp_dynamic_cast<const ObjFactory>(linObjFactory_);
  return objFactory->getThyraMatrix();
}
#endif

//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

} // end namespace panzer
