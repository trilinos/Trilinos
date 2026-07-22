// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EQUATION_SET_FACTORY_HPP
#define PANZER_EQUATION_SET_FACTORY_HPP

#include "Teuchos_RCP.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_CellData.hpp"
#include "Panzer_EquationSet_TemplateManager.hpp"
#include "Panzer_GlobalData.hpp"

namespace panzer {

  /** \brief Allocates and initializes an equation set template manager

     \param[in] params Input parameters to build the equation set
     \param[in] default_integration_order Default order for the integration rule.  NOTE: individual equation sets can override this based on parameters in the <code>plist</code>
     \param[in] cell_data The cell data
     \param[in] global_data  Global data
     \param[in] build_transient_support If true, the transient evaluators will be built, registered, and required in the Phalanx evaluation graph.

     Returns an RCP to a newly allocated EquationSet_TemplateManager.  
  */
  struct EquationSetFactory {
    virtual ~EquationSetFactory() = 0;

    virtual Teuchos::RCP<panzer::EquationSet_TemplateManager<panzer::Traits> >
    buildEquationSet(const Teuchos::RCP<Teuchos::ParameterList>& params,
		     const int& default_integration_order,
		     const panzer::CellData& cell_data,
		     const Teuchos::RCP<panzer::GlobalData>& global_data,
		     bool build_transient_support) const = 0;

  };

  inline EquationSetFactory::~EquationSetFactory() {}
  
}

#endif
