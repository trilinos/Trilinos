// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef   __myEquationSetFactory_hpp__
#define   __myEquationSetFactory_hpp__

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// Panzer
#include "Panzer_EquationSet_Factory.hpp"
#include "Panzer_EquationSet_Factory_Defines.hpp"
#include "Panzer_CellData.hpp"

// Files specifically for this example.
#include "myEquationSet.hpp"

// Create an object that can build our equation set.
PANZER_DECLARE_EQSET_TEMPLATE_BUILDER(MyEquationSet, MyEquationSet)

/**
 * 	\brief Our equation set factory.
 *
 * 	This class is used to build any equation sets we'd like to support.
 */
class
MyEquationSetFactory
  :
  public panzer::EquationSetFactory
{
  public:

    /**
     *  \brief Build our equation sets.
     *
     *  This routine creates all our various equation sets and adds them to an
     *  `EquationSet_TemplateManager`.
     *
     *  \param[in] params                    One of the entries (which is
     *                                       itself an unnamed `ParameterList`)
     *                                       in our physics block
     *                                       `ParameterList` in the input XML
     *                                       file.
     *  \param[in] default_integration_order An integration order to use if one
     *                                       is not specified in the input XML
     *                                       file.
     *  \param[in] cell_data                 This is passed to the
     *                                       `EquationSet_DefaultImpl`
     *                                       constructor.
     *  \param[in] global_data               This is passed to the
     *                                       `EquationSet_DefaultImpl`
     *                                       constructor.
     *  \param[in] build_transient_support   This is passed to the
     *                                       `EquationSet_DefaultImpl`
     *                                       constructor.
     *
     *  \returns An object that manages all the different equation sets we've
     *           constructed.
     */
    Teuchos::RCP<panzer::EquationSet_TemplateManager<panzer::Traits>>
    buildEquationSet(
      const Teuchos::RCP<Teuchos::ParameterList>& params,
      const int&                                  default_integration_order,
      const panzer::CellData&                     cell_data,
      const Teuchos::RCP<panzer::GlobalData>&     global_data,
      const bool                                  build_transient_support)
      const
    {
      using panzer::EquationSet_TemplateManager;
      using panzer::Traits;
      using std::logic_error;
      using std::string;
      using Teuchos::RCP;
      using Teuchos::rcp;
      RCP<EquationSet_TemplateManager<Traits>> eq_set =
        rcp(new EquationSet_TemplateManager<Traits>);
      bool found(false);
      PANZER_BUILD_EQSET_OBJECTS("MyEquationSet", MyEquationSet)
      if (not found)
      {
        string msg("Error:  The \"Equation Set\" with \"Type\" = \"" +
          params->get<string>("Type") + "\" is not a valid equation set "     \
          "identifier.  Please supply the correct factory.\n");
        TEUCHOS_TEST_FOR_EXCEPTION(true, logic_error, msg);
      } // end if (not found)
      return eq_set;
    } // end of buildEquationSet()

}; // end of class MyEquationSetFactory

#endif //  __myEquationSetFactory_hpp__
