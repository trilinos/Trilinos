// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Step01_ClosureModel_Factory_hpp__
#define __Step01_ClosureModel_Factory_hpp__

#include "Panzer_ClosureModel_Factory.hpp"

namespace user_app {

template<typename EvalT>
class ClosureModelFactory : public panzer::ClosureModelFactory<EvalT> {
public:

  Teuchos::RCP< std::vector< Teuchos::RCP<PHX::Evaluator<panzer::Traits> > > >
  buildClosureModels(const std::string& model_id,
                     const Teuchos::ParameterList& models,
                     const panzer::FieldLayoutLibrary& fl,
                     const Teuchos::RCP<panzer::IntegrationRule>& ir,
                     const Teuchos::ParameterList& default_params,
                     const Teuchos::ParameterList& user_data,
                     const Teuchos::RCP<panzer::GlobalData>& global_data,
                     PHX::FieldManager<panzer::Traits>& fm) const;

};

}

#include "Step01_ClosureModel_Factory_impl.hpp"

#endif
