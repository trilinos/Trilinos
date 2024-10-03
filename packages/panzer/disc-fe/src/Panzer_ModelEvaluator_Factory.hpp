// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_MODEL_EVALUATOR_FACTORY_DECL_HPP
#define PANZER_MODEL_EVALUATOR_FACTORY_DECL_HPP

#include "Teuchos_RCP.hpp"

namespace Thyra {
  template<typename Scalar> class ModelEvaluator;
}

namespace panzer {
  
  template <typename ScalarT, typename LO, typename GO>
  Teuchos::RCP<Thyra::ModelEvaluator<ScalarT> > 
  buildModelEvaluator(const RCP<panzer::FieldManagerBuilder & fmb,
		      const RCP<panzer::BlockedEpetraLinearObjFactory<panzer::Traits,int> >& lof);
    
}

#include "Panzer_ModelEvaluator_Factory_impl.hpp"

#endif
