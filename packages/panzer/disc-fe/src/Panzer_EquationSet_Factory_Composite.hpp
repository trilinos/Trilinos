// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Panzer_EquationSet_Factory.hpp"
#include "Teuchos_RCP.hpp"

namespace Teuchos {
  class ParameterList;
}

namespace panzer {

  class CellData;
  struct GlobalData;

  /** \brief A panzer::EquationSetFactory that delegates to a list of
    * other equation set factories.
    *
    * buildEquationSet() tries each constituent factory in order and
    * returns the first non-null equation set produced; an exception is
    * thrown if none of the factories can build one for the given input
    * parameter list.
    */
  class EquationSet_FactoryComposite : public panzer::EquationSetFactory {

    std::vector<Teuchos::RCP<panzer::EquationSetFactory> > m_factories;

  public:

    /// \brief Construct from the list of factories to try, in order.
    EquationSet_FactoryComposite(const std::vector<Teuchos::RCP<panzer::EquationSetFactory> >& factories);

    /** \brief Builds the equation set using the first constituent factory that succeeds for the given input parameter list.
      *
      * \param[in] input_plist the equation set's ParameterList, tried against each constituent factory in turn.
      * \param[in] default_integration_rule the integration order to use if input_plist does not specify one.
      * \param[in] cell_data the cell topology/dimension information for the physics block being built.
      * \param[in] global_data global data passed through to the constituent factory that builds the equation set.
      * \param[in] build_transient_support if true, the built equation set should also register time-derivative terms.
      */
    Teuchos::RCP<panzer::EquationSet_TemplateManager<panzer::Traits> >
    buildEquationSet(const Teuchos::RCP<Teuchos::ParameterList>& input_plist,
		     const int& default_integration_rule,
		     const panzer::CellData& cell_data,
		     const Teuchos::RCP<panzer::GlobalData>& global_data,
		     const bool build_transient_support) const;

  };

}

