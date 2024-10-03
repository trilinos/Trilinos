// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_ParameterList.hpp"
#include "LOCA_MultiContinuation_NaturalGroup.H"
#include "LOCA_MultiContinuation_NaturalConstraint.H"
#include "LOCA_MultiContinuation_ConstrainedGroup.H"

LOCA::MultiContinuation::NaturalGroup::NaturalGroup(
      const Teuchos::RCP<LOCA::GlobalData>& global_data,
      const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
      const Teuchos::RCP<Teuchos::ParameterList>& _continuationParams,
      const Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>& grp,
      const Teuchos::RCP<LOCA::MultiPredictor::AbstractStrategy>& pred,
      const std::vector<int>& paramIDs)
  : LOCA::MultiContinuation::ExtendedGroup(global_data, topParams,
                       _continuationParams,
                       grp, pred, paramIDs)
{
  bool skip_dfdp = continuationParams->get("Skip Parameter Derivative", true);
  Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface> cons
    = Teuchos::rcp(new LOCA::MultiContinuation::NaturalConstraint(
    globalData, Teuchos::rcp(this, false)));
  LOCA::MultiContinuation::ExtendedGroup::setConstraints(cons, skip_dfdp);
}

LOCA::MultiContinuation::NaturalGroup::NaturalGroup(
             const LOCA::MultiContinuation::NaturalGroup& source,
             NOX::CopyType type)
  : LOCA::MultiContinuation::ExtendedGroup(source, type)
{
  Teuchos::rcp_dynamic_cast<LOCA::MultiContinuation::NaturalConstraint>(conGroup->getConstraints())->setNaturalGroup(Teuchos::rcp(this, false));
}


LOCA::MultiContinuation::NaturalGroup::~NaturalGroup()
{
}

NOX::Abstract::Group&
LOCA::MultiContinuation::NaturalGroup::operator=(
                      const NOX::Abstract::Group& source)
{
  copy(source);
  return *this;
}

Teuchos::RCP<NOX::Abstract::Group>
LOCA::MultiContinuation::NaturalGroup::clone(NOX::CopyType type) const
{
  return Teuchos::rcp(new LOCA::MultiContinuation::NaturalGroup(*this, type));
}

void
LOCA::MultiContinuation::NaturalGroup::copy(const NOX::Abstract::Group& src)
{

  // Protect against A = A
  if (this != &src) {
    LOCA::MultiContinuation::ExtendedGroup::copy(src);
  }
}


