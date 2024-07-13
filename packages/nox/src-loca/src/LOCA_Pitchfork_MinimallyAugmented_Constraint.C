// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_Pitchfork_MinimallyAugmented_Constraint.H"
#include "LOCA_Pitchfork_MinimallyAugmented_AbstractGroup.H"
#include "LOCA_BorderedSolver_AbstractStrategy.H"
#include "LOCA_Parameter_SublistParser.H"
#include "LOCA_GlobalData.H"
#include "LOCA_ErrorCheck.H"
#include "LOCA_Factory.H"
#include "NOX_Utils.H"
#include "Teuchos_ParameterList.hpp"

LOCA::Pitchfork::MinimallyAugmented::Constraint::
Constraint(
    const Teuchos::RCP<LOCA::GlobalData>& global_data,
    const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
    const Teuchos::RCP<Teuchos::ParameterList>& pfParams,
    const Teuchos::RCP<LOCA::Pitchfork::MinimallyAugmented::AbstractGroup>& g,
    const Teuchos::RCP<const NOX::Abstract::Vector>& psi,
    int bif_param) :
  LOCA::TurningPoint::MinimallyAugmented::Constraint(global_data, topParams,
                             pfParams, g, bif_param),
  pf_grp(g),
  psi_vector(psi),
  dgdx(psi->createMultiVector(2, NOX::ShapeCopy)),
  pf_constraints(2, 1)
{
}

LOCA::Pitchfork::MinimallyAugmented::Constraint::
Constraint(const LOCA::Pitchfork::MinimallyAugmented::Constraint& source,
       NOX::CopyType type) :
  LOCA::TurningPoint::MinimallyAugmented::Constraint(source, type),
  pf_grp(Teuchos::null),
  psi_vector(source.psi_vector),
  dgdx(source.dgdx->clone(type)),
  pf_constraints(source.pf_constraints)
{
  // We don't explicitly copy the group because the constrained group
  // will do that
}

LOCA::Pitchfork::MinimallyAugmented::Constraint::
~Constraint()
{
}

void
LOCA::Pitchfork::MinimallyAugmented::Constraint::
setGroup(const Teuchos::RCP<LOCA::TurningPoint::MinimallyAugmented::AbstractGroup>& g)
{
  LOCA::TurningPoint::MinimallyAugmented::Constraint::setGroup(g);
  pf_grp = Teuchos::rcp_dynamic_cast<LOCA::Pitchfork::MinimallyAugmented::AbstractGroup>(g,true);
}

void
LOCA::Pitchfork::MinimallyAugmented::Constraint::
copy(const LOCA::MultiContinuation::ConstraintInterface& src)
{
  const LOCA::Pitchfork::MinimallyAugmented::Constraint& source =
  dynamic_cast<const LOCA::Pitchfork::MinimallyAugmented::Constraint&>(src);

  if (this != &source) {
    LOCA::TurningPoint::MinimallyAugmented::Constraint::copy(src);
    psi_vector = source.psi_vector;
    *dgdx = *source.dgdx;
    pf_constraints.assign(source.pf_constraints);

    // We don't explicitly copy the group because the constrained group
    // will do that
  }
}

Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface>
LOCA::Pitchfork::MinimallyAugmented::Constraint::
clone(NOX::CopyType type) const
{
  return Teuchos::rcp(new Constraint(*this, type));
}

int
LOCA::Pitchfork::MinimallyAugmented::Constraint::
numConstraints() const
{
  return 2;
}

NOX::Abstract::Group::ReturnType
LOCA::Pitchfork::MinimallyAugmented::Constraint::
computeConstraints()
{
  if (isValidConstraints)
    return NOX::Abstract::Group::Ok;

  // Compute sigma
  NOX::Abstract::Group::ReturnType status =
    LOCA::TurningPoint::MinimallyAugmented::Constraint::computeConstraints();
  pf_constraints(0,0) = constraints(0,0);


  // Compute <psi,x>
  pf_constraints(1,0) = pf_grp->innerProduct(*psi_vector, pf_grp->getX());

  return status;
}

NOX::Abstract::Group::ReturnType
LOCA::Pitchfork::MinimallyAugmented::Constraint::
computeDX()
{
  if (isValidDX)
    return NOX::Abstract::Group::Ok;

  // Compute sigma_x
  NOX::Abstract::Group::ReturnType status =
    LOCA::TurningPoint::MinimallyAugmented::Constraint::computeDX();
  (*dgdx)[0] = (*sigma_x)[0];

  // Compute <psi,x>_x = psi
  (*dgdx)[1] = *psi_vector;

  return status;
}

NOX::Abstract::Group::ReturnType
LOCA::Pitchfork::MinimallyAugmented::Constraint::
computeDP(const std::vector<int>& paramIDs,
      NOX::Abstract::MultiVector::DenseMatrix& dgdp,
      bool isValidG)
{
  // Compute sigma_p
  NOX::Abstract::MultiVector::DenseMatrix dgdp_sub(Teuchos::View, dgdp, 1,
                           paramIDs.size()+1, 0, 0);
  NOX::Abstract::Group::ReturnType status =
    LOCA::TurningPoint::MinimallyAugmented::Constraint::computeDP(paramIDs,
                                  dgdp_sub,
                                  isValidG);

  // Compute <psi,x>_p
  if (!isValidG)
    dgdp(1,0) = pf_grp->innerProduct(*psi_vector, pf_grp->getX());
  for (unsigned int i=0; i<paramIDs.size(); i++)
    dgdp(1,i+1) = 0.0;

  return status;
}

const NOX::Abstract::MultiVector::DenseMatrix&
LOCA::Pitchfork::MinimallyAugmented::Constraint::
getConstraints() const
{
  return pf_constraints;
}

const NOX::Abstract::MultiVector*
LOCA::Pitchfork::MinimallyAugmented::Constraint::
getDX() const
{
  return dgdx.get();
}
