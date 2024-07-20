// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_BCSTRATEGY_FACTORY_HPP
#define PANZER_BCSTRATEGY_FACTORY_HPP

#include "Teuchos_RCP.hpp"
#include "Panzer_Traits.hpp"

namespace panzer {
  
  class BC;
  template<typename T> class BCStrategy_TemplateManager;
  struct GlobalData;

  /** \brief Interface for constructing a BCStrategy_TemplateManager

     \param bc [in] A description of the boundary condition.
     \param globa_data [in] a nonnull rcp to a GlobalData object.

     Returns a nonnull RCP to the BCStrategy_TemplateManager.  The
     object should throw an exception if the BCStrategy object fails
     to build correctly.
  */
  struct BCStrategyFactory {

    BCStrategyFactory() {}
    virtual ~BCStrategyFactory() {}

    virtual Teuchos::RCP<panzer::BCStrategy_TemplateManager<panzer::Traits> >
    buildBCStrategy(const panzer::BC& bc,
		    const Teuchos::RCP<panzer::GlobalData>& global_data) const = 0;

  };
  
}

#endif
