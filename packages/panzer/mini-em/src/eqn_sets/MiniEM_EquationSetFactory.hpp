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

  class EquationSetFactory : public panzer::EquationSetFactory {

  public:

    EquationSetFactory(const Teuchos::RCP<panzer::GlobalEvaluationDataContainer> & gedc=Teuchos::null) :
       m_gedc(gedc) {}

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
