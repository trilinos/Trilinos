// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef   __myClosureModelFactory_hpp__
#define   __myClosureModelFactory_hpp__

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// Panzer
#include "Panzer_ClosureModel_Factory.hpp"

/**
 *  \brief Our closure model factory.
 *
 *  This class is used to build our closure models.
 */
template<typename EvalT>
class MyClosureModelFactory
  :
  public panzer::ClosureModelFactory<EvalT>
{
  public:
    
    /**
     *  \brief Build the closure models.
     *
     *  This routine builds all the evaluators for all the various closure
     *  models we'll support.  In our case, there is only one, corresponding to
     *  our source term.
     *
     *  \param[in] modelId       The closure model ID, which is the "name" of a
     *                           `ParameterList` in the "Closure Models"
     *                           `ParameterList` in the input XML file.
     *  \param[in] models        The "Closure Models" `ParameterList` from the
     *                           input XML file.
     *  \param[in] fl            This is unused in this routine, though it is
     *                           part of the `ClosureModelFactory` interface.
     *  \param[in] ir            The integration rule that is used in creating
     *                           our closure model objects.
     *  \param[in] defaultParams This is unused in this routine, though it is
     *                           part of the `ClosureModelFactory` interface.
     *  \param[in] userData      This is unused in this routine, though it is
     *                           part of the `ClosureModelFactory` interface.
     *  \param[in] globalData    This is unused in this routine, though it is
     *                           part of the `ClosureModelFactory` interface.
     *  \param[in] fm            This is unused in this routine, though it is
     *                           part of the `ClosureModelFactory` interface.
     *
     *  \returns A list of evaluators corresponding to all the various closure
     *           models we support.
     */
    Teuchos::RCP<std::vector<Teuchos::RCP<PHX::Evaluator<panzer::Traits>>>>
    buildClosureModels(
      const std::string&                           modelId,
      const Teuchos::ParameterList&                models,
      const panzer::FieldLayoutLibrary&            fl,
      const Teuchos::RCP<panzer::IntegrationRule>& ir,
      const Teuchos::ParameterList&                defaultParams,
      const Teuchos::ParameterList&                userData,
      const Teuchos::RCP<panzer::GlobalData>&      globalData,
      PHX::FieldManager<panzer::Traits>&           fm) const;

}; // end of class MyClosureModelFactory

#include "myClosureModelFactoryImpl.hpp"

#endif // __myClosureModelFactory_hpp__
