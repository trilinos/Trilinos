// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef   __myBCStrategyFactory_hpp__
#define   __myBCStrategyFactory_hpp__

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// Panzer
#include "Panzer_BCStrategy_Factory.hpp"
#include "Panzer_BCStrategy_Factory_Defines.hpp"
#include "Panzer_BCStrategy_TemplateManager.hpp"
#include "Panzer_GlobalData.hpp"
#include "Panzer_Traits.hpp"

// Teuchos
#include "Teuchos_RCP.hpp"

// Files specifically for this example.
#include "myBCStrategy.hpp"

// Create an object that can build our boundary condition strategy.
PANZER_DECLARE_BCSTRATEGY_TEMPLATE_BUILDER(MyBCStrategy, MyBCStrategy)

/**
 *  \brief Our boundary condition strategy factory.
 *
 *  This class is used to build any boundary conditions we'd like to support.
 */
class MyBCStrategyFactory
  :
  public panzer::BCStrategyFactory
{
  public:

    /**
     * 	\brief Build our boundary condition strategy.
     *
     * 	Create all of the various boundary conditions that we'll support.  For
     * 	our specific example, there's only one (MyBCStrategy), but you could
     * 	have as many as you like.
     *
     * 	\param[in] bc          A boundary condition corresponding to one of the
     * 	                       BCs specified in the input XML file.
     * 	\param[in] global_data The global data object, which stores the
     * 	                       parameter library and the default output stream.
     *
     * 	\returns An object that manages all the different boundary conditions
     * 	         we've constructed.
     */
    Teuchos::RCP<panzer::BCStrategy_TemplateManager<panzer::Traits>>
    buildBCStrategy(
      const panzer::BC&                       bc,
      const Teuchos::RCP<panzer::GlobalData>& global_data) const
    {
      using panzer::BCStrategy_TemplateManager;
      using panzer::Traits;
      using std::logic_error;
      using Teuchos::RCP;
      using Teuchos::rcp;

      // Create the BCStrategy_TemplateManager and add our BCStrategy to it.
      RCP<BCStrategy_TemplateManager<Traits>> bcs_tm =
        rcp(new BCStrategy_TemplateManager<Traits>);
      bool found(false);
      PANZER_BUILD_BCSTRATEGY_OBJECTS("MyBCStrategy", MyBCStrategy)
      TEUCHOS_TEST_FOR_EXCEPTION(not found, logic_error, "Error:  The "       \
        "BCStrategy called \"" << bc.strategy() << "\" is not a valid "       \
        "identifier in the BCStrategyFactory.  Either add a valid "           \
        "implementation to your factory or fix your input file.  The "        \
        "relevant boundary condition is:\n\n" << bc << std::endl)
      return bcs_tm;
    } // end of buildBCStrategy()

}; // end of class MyBCStrategyFactory

#endif // __myBCStrategyFactory_hpp__
