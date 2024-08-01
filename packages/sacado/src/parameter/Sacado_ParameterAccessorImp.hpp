// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "Sacado_ParameterRegistration.hpp"

template <typename EvalType, typename EvalTypeTraits>
void Sacado::ParameterAccessor<EvalType,EvalTypeTraits>::
registerSacadoParameter(const std::string& name,
                        ParamLib& paramLib)
{
  pr_.push_back(Teuchos::rcp(
                  new Sacado::ParameterRegistration<EvalType,EvalTypeTraits>(
                    name, this, paramLib)));
}

template <typename EvalType, typename EvalTypeTraits>
void Sacado::ParameterAccessor<EvalType,EvalTypeTraits>::
registerSacadoParameter(const std::string& name,
                        const Teuchos::RCP<ParamLib>& paramLib)
{
  pr_.push_back(Teuchos::rcp(
                  new Sacado::ParameterRegistration<EvalType,EvalTypeTraits>(
                    name, this, paramLib)));
}
