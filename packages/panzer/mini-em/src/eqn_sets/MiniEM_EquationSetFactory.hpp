#ifndef _MiniEM_EquationSetFactory_hpp_
#define _MiniEM_EquationSetFactory_hpp_

#include "Panzer_EquationSet_Factory.hpp"
#include "Panzer_EquationSet_Factory_Defines.hpp"
#include "Panzer_CellData.hpp"

#include "MiniEM_EquationSet_Maxwell.hpp"

#include "MiniEM_AuxiliaryEquationSet_MassMatrix.hpp"
#include "MiniEM_AuxiliaryEquationSet_CurlCurl.hpp"
#include "MiniEM_AuxiliaryEquationSet_WeakGradient.hpp"
#include "MiniEM_AuxiliaryEquationSet_MACROS.hpp"

namespace mini_em {

  PANZER_DECLARE_EQSET_TEMPLATE_BUILDER(EquationSet_Maxwell, EquationSet_Maxwell)

  AUX_DECLARE_EQSET_TEMPLATE_BUILDER(AuxiliaryEquationSet_MassMatrix, AuxiliaryEquationSet_MassMatrix)

  AUX_DECLARE_EQSET_TEMPLATE_BUILDER(AuxiliaryEquationSet_CurlCurl, AuxiliaryEquationSet_CurlCurl)

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

      AUX_BUILD_EQSET_OBJECTS("Auxiliary Mass Matrix",   AuxiliaryEquationSet_MassMatrix)

      AUX_BUILD_EQSET_OBJECTS("Auxiliary Curl Curl",   AuxiliaryEquationSet_CurlCurl)

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
