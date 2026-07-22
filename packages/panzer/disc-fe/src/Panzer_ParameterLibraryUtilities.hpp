// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_PARAMETER_LIBRARY_UTILITIES_HPP
#define PANZER_PARAMETER_LIBRARY_UTILITIES_HPP

#include "Panzer_ParameterLibrary.hpp"
#include "Panzer_ScalarParameterEntry.hpp"
#include "Teuchos_RCP.hpp"

namespace panzer {

  /** \brief Allocates a parameter entry and registers with parameter library
      
  \relates ParameterLibraryAcceptor
  */
  template<typename EvaluationType>
  Teuchos::RCP<panzer::ScalarParameterEntry<EvaluationType> >
  createAndRegisterScalarParameter(const std::string name,
				   panzer::ParamLib& pl);

  template<typename EvaluationType>
  Teuchos::RCP<panzer::ScalarParameterEntry<EvaluationType> >
  accessScalarParameter(const std::string name, panzer::ParamLib& pl);

  void registerScalarParameter(const std::string name,panzer::ParamLib& pl,double realValue);
  

}

#include "Panzer_ParameterLibraryUtilities_impl.hpp"

#endif
