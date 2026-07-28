// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_BCSTRATEGY_FACTORY_COMPOSITE_DECL_HPP
#define PANZER_BCSTRATEGY_FACTORY_COMPOSITE_DECL_HPP

#include "Panzer_BCStrategy_Factory.hpp"
#include "Teuchos_RCP.hpp"
#include <vector>

namespace panzer {
  
  /** \brief A panzer::BCStrategyFactory that contains 
    * an internal vector of user BC strategy factories.
    *
    * buildBCStrategy() tries each constituent factory in order and
    * returns the first non-null BC strategy produced; an exception is
    * thrown if none of the factories recognize the given boundary
    * condition's strategy identifier.
    */
  struct BCFactoryComposite : public panzer::BCStrategyFactory {

  public:

    /// \brief Construct from the list of factories to try, in order.
    BCFactoryComposite(const std::vector<Teuchos::RCP<panzer::BCStrategyFactory> >& factories);

    /** \brief Builds the BC strategy using the first constituent factory that recognizes the given boundary condition.
      *
      * \param[in] bc the boundary condition to build a strategy for; its strategy identifier selects which constituent factory (if any) handles it.
      * \param[in] global_data global data passed through to the constituent factory that builds the strategy.
      */
    Teuchos::RCP<panzer::BCStrategy_TemplateManager<panzer::Traits> >
    buildBCStrategy(const panzer::BC& bc, const Teuchos::RCP<panzer::GlobalData>& global_data) const;

  private:

    std::vector<Teuchos::RCP<panzer::BCStrategyFactory> > m_bc_strategy_factories;
    
  };
  
}

#endif
