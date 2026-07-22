// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_PARAMETER_LIBRARY_UTILITIES_IMPL_HPP
#define PANZER_PARAMETER_LIBRARY_UTILITIES_IMPL_HPP

namespace panzer {

  template<typename EvaluationType>
  Teuchos::RCP<panzer::ScalarParameterEntry<EvaluationType> >
  createAndRegisterScalarParameter(const std::string name, 
				   panzer::ParamLib& pl)
  {
    if (!pl.isParameter(name))
      pl.addParameterFamily(name,true,false);
    
    Teuchos::RCP<panzer::ScalarParameterEntry<EvaluationType> > entry;
    
    if (pl.isParameterForType<EvaluationType>(name)) {
      Teuchos::RCP<Sacado::ScalarParameterEntry<EvaluationType,panzer::EvaluationTraits> > sacado_entry =
	pl.getEntry<EvaluationType>(name);
      entry = Teuchos::rcp_dynamic_cast<panzer::ScalarParameterEntry<EvaluationType> >(sacado_entry);
    }
    else {
      entry = Teuchos::rcp(new panzer::ScalarParameterEntry<EvaluationType>);
      entry->setValue(NAN);
      pl.addEntry<EvaluationType>(name,entry);
    }

    return entry;
  }

  template<typename EvaluationType>
  Teuchos::RCP<panzer::ScalarParameterEntry<EvaluationType> >
  accessScalarParameter(const std::string name, panzer::ParamLib& pl)
  {
    Teuchos::RCP<Sacado::ScalarParameterEntry<EvaluationType,panzer::EvaluationTraits> > sacado_entry =
      pl.getEntry<EvaluationType>(name);
    return Teuchos::rcp_dynamic_cast<panzer::ScalarParameterEntry<EvaluationType> >(sacado_entry,true);
  }

}

#endif
