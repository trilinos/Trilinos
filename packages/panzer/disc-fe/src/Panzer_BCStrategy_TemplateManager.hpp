// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_BCSTRATEGY_TEMPLATE_MANAGER_HPP
#define PANZER_BCSTRATEGY_TEMPLATE_MANAGER_HPP

#include "Phalanx_TemplateManager.hpp"

#include "Panzer_Traits.hpp"
#include "Panzer_BCStrategy_Base.hpp"
#include "Panzer_BCStrategy.hpp"

#include "Sacado_mpl_placeholders.hpp"

namespace panzer {

  template<typename Traits>
    class BCStrategy_TemplateManager : 
    public PHX::TemplateManager<typename Traits::EvalTypes,
				panzer::BCStrategyBase,
                                panzer::BCStrategy<_> > {

  public:

    BCStrategy_TemplateManager() {}

    ~BCStrategy_TemplateManager() {}

  };

} 

#endif 
