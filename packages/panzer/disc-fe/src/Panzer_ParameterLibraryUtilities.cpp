// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_PARAMETER_LIBRARY_UTILITIES_CPP
#define PANZER_PARAMETER_LIBRARY_UTILITIES_CPP

#include "Panzer_ParameterLibraryUtilities.hpp"

namespace panzer {

struct RegisterScalarParameter {
  std::string name;
  Teuchos::Ptr<panzer::ParamLib> pl;

  template <typename T>
  void apply() const
  { createAndRegisterScalarParameter<T>(name,*pl); }
};

void registerScalarParameter(const std::string name,panzer::ParamLib& pl,double realValue)
{
  RegisterScalarParameter rsp; 
  rsp.name = name;
  rsp.pl = Teuchos::ptrFromRef(pl);

  rsp.apply<panzer::Traits::Residual>();
  rsp.apply<panzer::Traits::Jacobian>();
  rsp.apply<panzer::Traits::Tangent>();
#ifdef    Panzer_BUILD_HESSIAN_SUPPORT
  rsp.apply<panzer::Traits::Hessian>();
#endif // Panzer_BUILD_HESSIAN_SUPPORT

  pl.setRealValueForAllTypes(name,realValue);
}

}

#endif
