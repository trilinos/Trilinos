// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "LOCA_Thyra_GroupWrapper.H"              // class definition

LOCA::Thyra::GroupWrapper::GroupWrapper(
        const Teuchos::RCP<LOCA::GlobalData>& global_data,
        const NOX::Thyra::Vector& initial_guess,
        const Teuchos::RCP< ::Thyra::ModelEvaluator<double> >& model,
        const LOCA::ParameterVector& p,
        int p_index,
        bool impl_dfdp) :
  NOX::Thyra::Group(initial_guess, model),
  LOCA::Abstract::Group(global_data),
  LOCA::Thyra::Group(global_data, initial_guess, model, p, p_index, impl_dfdp)
{
}

LOCA::Thyra::GroupWrapper::GroupWrapper(const LOCA::Thyra::GroupWrapper& source,
               NOX::CopyType type) :
  NOX::Thyra::Group(source, type),
  LOCA::Abstract::Group(source, type),
  LOCA::Thyra::Group(source, type)
{
}

LOCA::Thyra::GroupWrapper::~GroupWrapper()
{
}

LOCA::Thyra::GroupWrapper&
LOCA::Thyra::GroupWrapper::operator=(const LOCA::Thyra::GroupWrapper& source)
{
  if (this != &source) {
    LOCA::Thyra::Group::operator=(source);
  }
  return *this;
}

NOX::Abstract::Group&
LOCA::Thyra::GroupWrapper::operator=(const NOX::Abstract::Group& source)
{
  operator=(dynamic_cast<const GroupWrapper&> (source));
  return *this;
}

NOX::Abstract::Group&
LOCA::Thyra::GroupWrapper::operator=(const NOX::Thyra::Group& source)
{
  operator=(dynamic_cast<const GroupWrapper&> (source));
  return *this;
}

LOCA::Thyra::Group&
LOCA::Thyra::GroupWrapper::operator=(const LOCA::Thyra::Group& source)
{
  operator=(dynamic_cast<const GroupWrapper&> (source));
  return *this;
}

Teuchos::RCP<NOX::Abstract::Group>
LOCA::Thyra::GroupWrapper::clone(NOX::CopyType type) const
{
  return Teuchos::rcp(new LOCA::Thyra::GroupWrapper(*this, type));
}

