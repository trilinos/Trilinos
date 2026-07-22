// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*
 * SillyEquationSet.hpp
 *
 *  Created on: Sep 9, 2016
 *      Author: mbetten
 */

#ifndef TRILINOS_PACKAGES_PANZER_ADAPTERS_STK_TEST_FACE_TO_ELEM_SILLYEQUATIONSET_HPP_
#define TRILINOS_PACKAGES_PANZER_ADAPTERS_STK_TEST_FACE_TO_ELEM_SILLYEQUATIONSET_HPP_



#include "Panzer_Normals.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_CellData.hpp"
#include "Panzer_Workset.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_CommonArrayFactories.hpp"
#include "Panzer_EquationSet_Factory.hpp"
#include "Panzer_EquationSet_Factory_Defines.hpp"
#include "Panzer_CellData.hpp"
#include "Phalanx_FieldManager.hpp"
#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_SetupUtilities.hpp"
#include "Panzer_STKConnManager.hpp"
#include "Panzer_EquationSet_DefaultImpl.hpp"



template <typename EvalT>
class EquationSet_MeshCoords : public panzer::EquationSet_DefaultImpl<EvalT> {

public:

  EquationSet_MeshCoords(const Teuchos::RCP<Teuchos::ParameterList>& params,
      const int& default_integration_order,
      const panzer::CellData& cell_data,
      const Teuchos::RCP<panzer::GlobalData>& gd,
      const bool build_transient_support) :
        panzer::EquationSet_DefaultImpl<EvalT>(params,default_integration_order,cell_data,gd,build_transient_support ) {
    std::string basis_type = "HGrad"; // use nodal linears
    int basis_order = 1;
    std::string model_id = params->get<std::string>("Model ID");
    int integration_order = params->get<int>("Integration Order");

    dimension_ = cell_data.baseCellDimension();
    std::string dof_names[3] = {"COORDX","COORDY","COORDZ"};
    std::vector<std::string> coord;
    for(int i=0;i<dimension_;i++) {
      std::string dof_name = dof_names[i];
      coord.push_back(dof_name);

      this->addDOF(dof_name,basis_type,basis_order,integration_order);
    }
    this->setupDOFs();

  }

  void buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& /* fm */,
      const panzer::FieldLibrary& /* field_library */,
      const Teuchos::ParameterList& /* user_data */) const{}
private:
  int dimension_;

};

  PANZER_DECLARE_EQSET_TEMPLATE_BUILDER(EquationSet_MeshCoords,
      EquationSet_MeshCoords)


class MyFactory : public panzer::EquationSetFactory {

  public:

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
      PANZER_BUILD_EQSET_OBJECTS("MeshCoords", EquationSet_MeshCoords)
      if (!found)
      {
        std::string msg = "Error - the \"Equation Set\" called \"" +
          params->get<std::string>("Type") + "\" is not a valid equation " + 
          "set identifier.  Please supply the correct factory.\n";
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, msg);
      }
      return eq_set;
    }

};


#endif /* TRILINOS_PACKAGES_PANZER_ADAPTERS_STK_TEST_FACE_TO_ELEM_SILLYEQUATIONSET_HPP_ */
