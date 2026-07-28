// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _MiniEM_EquationSetFactory_hpp_
#define _MiniEM_EquationSetFactory_hpp_

#include "Panzer_EquationSet_Factory.hpp"
#include "Panzer_EquationSet_Factory_Defines.hpp"
#include "Panzer_CellData.hpp"

#include "MiniEM_EquationSet_Maxwell.hpp"
#include "MiniEM_EquationSet_Darcy.hpp"

#include "MiniEM_AuxiliaryEquationSet_MassMatrix.hpp"
#include "MiniEM_AuxiliaryEquationSet_SchurComplement.hpp"
#include "MiniEM_AuxiliaryEquationSet_DarcySchurComplement.hpp"
#include "MiniEM_AuxiliaryEquationSet_ProjectedSchurComplement.hpp"
#include "MiniEM_AuxiliaryEquationSet_ProjectedDarcySchurComplement.hpp"
#include "MiniEM_AuxiliaryEquationSet_WeakGradient.hpp"
#include "MiniEM_AuxiliaryEquationSet_MACROS.hpp"

namespace mini_em {

  PANZER_DECLARE_EQSET_TEMPLATE_BUILDER(EquationSet_Maxwell, EquationSet_Maxwell)

  PANZER_DECLARE_EQSET_TEMPLATE_BUILDER(EquationSet_Darcy, EquationSet_Darcy)

  AUX_DECLARE_EQSET_TEMPLATE_BUILDER(AuxiliaryEquationSet_MassMatrix, AuxiliaryEquationSet_MassMatrix)

  AUX_DECLARE_EQSET_TEMPLATE_BUILDER(AuxiliaryEquationSet_SchurComplement, AuxiliaryEquationSet_SchurComplement)

  AUX_DECLARE_EQSET_TEMPLATE_BUILDER(AuxiliaryEquationSet_DarcySchurComplement, AuxiliaryEquationSet_DarcySchurComplement)

  AUX_DECLARE_EQSET_TEMPLATE_BUILDER(AuxiliaryEquationSet_ProjectedSchurComplement, AuxiliaryEquationSet_ProjectedSchurComplement)

  AUX_DECLARE_EQSET_TEMPLATE_BUILDER(AuxiliaryEquationSet_ProjectedDarcySchurComplement, AuxiliaryEquationSet_ProjectedDarcySchurComplement)

  AUX_DECLARE_EQSET_TEMPLATE_BUILDER(AuxiliaryEquationSet_WeakGradient, AuxiliaryEquationSet_WeakGradient)

  /** \brief mini-em's panzer::EquationSetFactory, recognizing the
    * "Maxwell" and "Darcy" primary equation set identifiers (see
    * mini_em::EquationSet_Maxwell, mini_em::EquationSet_Darcy) and the
    * "Auxiliary Mass Matrix", "Auxiliary SchurComplement", "Auxiliary
    * DarcySchurComplement", "Auxiliary ProjectedSchurComplement",
    * "Auxiliary ProjectedDarcySchurComplement", and "Auxiliary Weak
    * Gradient" auxiliary equation set identifiers used to build
    * preconditioner operators.
    */
  class EquationSetFactory : public panzer::EquationSetFactory {

  public:

    /// \brief Construct with the global evaluation data container that auxiliary equation sets scatter their operators into (unused by the primary Maxwell/Darcy equation sets).
    EquationSetFactory(const Teuchos::RCP<panzer::GlobalEvaluationDataContainer> & gedc=Teuchos::null) :
       m_gedc(gedc) {}

    /** \brief Builds the named equation set ("Maxwell", "Darcy", or one of the "Auxiliary ..." identifiers per the class description); throws std::logic_error for any other "Type".
      *
      * \param[in] params the equation set's ParameterList; its "Type" entry selects which equation set is built.
      * \param[in] default_integration_order the integration order to use if params does not specify one.
      * \param[in] cell_data the cell topology/dimension information for the physics block being built.
      * \param[in] global_data global data passed through to the built equation set.
      * \param[in] build_transient_support if true, the built equation set should also register time-derivative terms.
      */
    Teuchos::RCP<panzer::EquationSet_TemplateManager<panzer::Traits> >
    buildEquationSet(const Teuchos::RCP<Teuchos::ParameterList>& params,
        const int& default_integration_order,
        const panzer::CellData& cell_data,
        const Teuchos::RCP<panzer::GlobalData>& global_data,
        const bool build_transient_support) const
    {
      Teuchos::RCP<panzer::EquationSet_TemplateManager<panzer::Traits> > eq_set= 
          Teuchos::rcp(new panzer::EquationSet_TemplateManager<panzer::Traits>);

      bool found = false;

      PANZER_BUILD_EQSET_OBJECTS("Maxwell",              EquationSet_Maxwell)

      PANZER_BUILD_EQSET_OBJECTS("Darcy",                EquationSet_Darcy)

      AUX_BUILD_EQSET_OBJECTS("Auxiliary Mass Matrix",   AuxiliaryEquationSet_MassMatrix)

      AUX_BUILD_EQSET_OBJECTS("Auxiliary SchurComplement",   AuxiliaryEquationSet_SchurComplement)

      AUX_BUILD_EQSET_OBJECTS("Auxiliary DarcySchurComplement",   AuxiliaryEquationSet_DarcySchurComplement)

      AUX_BUILD_EQSET_OBJECTS("Auxiliary ProjectedSchurComplement",   AuxiliaryEquationSet_ProjectedSchurComplement)

      AUX_BUILD_EQSET_OBJECTS("Auxiliary ProjectedDarcySchurComplement",   AuxiliaryEquationSet_ProjectedDarcySchurComplement)

      AUX_BUILD_EQSET_OBJECTS("Auxiliary Weak Gradient", AuxiliaryEquationSet_WeakGradient)

      if (!found) {
        std::string msg = "Error - the \"Equation Set\" with \"Type\" = \"" + params->get<std::string>("Type") +
            "\" is not a valid equation set identifier. Please supply the correct factory.\n";
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
      }

      return eq_set;
    }

  private:

    Teuchos::RCP<panzer::GlobalEvaluationDataContainer> m_gedc;

  };

}

#endif
