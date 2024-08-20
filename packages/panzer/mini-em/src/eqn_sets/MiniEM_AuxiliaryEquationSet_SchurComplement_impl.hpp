// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _MiniEM_AuxiliaryEquationSet_SchurComplement_impl_hpp_
#define _MiniEM_AuxiliaryEquationSet_SchurComplement_impl_hpp_

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"

#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_FieldManager.hpp"

#include "Panzer_String_Utilities.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"

// include evaluators here
#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_Integrator_BasisTimesVector.hpp"
#include "Panzer_Integrator_BasisTimesTensorTimesVector.hpp"
#include "Panzer_Integrator_CurlBasisDotVector.hpp"
#include "Panzer_ScalarToVector.hpp"
#include "Panzer_Sum.hpp"
#include "Panzer_Constant.hpp"
#include "Panzer_BlockedDOFManager.hpp"
#include "Panzer_BlockedTpetraLinearObjFactory.hpp"
#ifdef PANZER_HAVE_EPETRA_STACK
#include "Panzer_BlockedEpetraLinearObjFactory.hpp"
#endif
#include "Panzer_ReorderADValues_Evaluator.hpp"
#include "Panzer_LOCPair_GlobalEvaluationData.hpp"

// ***********************************************************************
template <typename EvalT>
mini_em::AuxiliaryEquationSet_SchurComplement<EvalT>::
AuxiliaryEquationSet_SchurComplement(
                   const Teuchos::RCP<panzer::GlobalEvaluationDataContainer> & gedc,
                   const Teuchos::RCP<Teuchos::ParameterList>& params,
		   const int& default_integration_order,
		   const panzer::CellData& cell_data,
		   const Teuchos::RCP<panzer::GlobalData>& global_data,
		   const bool build_transient_support) :
  panzer::EquationSet_DefaultImpl<EvalT>(params, default_integration_order, cell_data, global_data, build_transient_support)
{
  // Options
  {
    Teuchos::ParameterList valid_parameters;
    this->setDefaultValidParameters(valid_parameters);
    valid_parameters.set("Model ID","","Closure model id associated with this equation set");
    valid_parameters.set("DOF Name","","Name of DOF to construct time derivative for");
    valid_parameters.set("Permittivity","epsilon","Permittivity");
    valid_parameters.set("Conductivity","sigma","Conductivity");
    valid_parameters.set("Inverse Permeability","1/mu","Inverse Permeability");
    valid_parameters.set("Basis Type","HCurl","Type of Basis to use");
    valid_parameters.set("Basis Order",1,"Order of the basis");
    valid_parameters.set("Integration Order",default_integration_order,"Order of the integration rule");

    params->validateParametersAndSetDefaults(valid_parameters);
  }
  std::string model_id = params->get<std::string>("Model ID");
  dof_name = params->get<std::string>("DOF Name");

  std::string basis_type = params->get<std::string>("Basis Type");
  int basis_order = params->get<int>("Basis Order");
  int integration_order = params->get<int>("Integration Order");

  permittivity_ = params->get<std::string>("Permittivity");
  conductivity_ = params->get<std::string>("Conductivity");
  inversePermeability_ = params->get<std::string>("Inverse Permeability");

  // ********************
  // Assemble DOF names and Residual names
  // ********************

  m_gedc = gedc;

  m_dof_names = Teuchos::rcp(new std::vector<std::string>);
  m_dof_names->push_back(dof_name);

  this->addDOF(dof_name,basis_type,basis_order,integration_order,"AUX_MASS_RESIDUAL_"+dof_name);
  // this->addDOF(dof_name,basis_type,basis_order,integration_order);
  this->addDOFCurl(dof_name,"Curl_"+dof_name);

  this->addDOFTimeDerivative(dof_name);

  this->addClosureModel(model_id);

  this->setupDOFs();
}

// ***********************************************************************
template <typename EvalT>
void mini_em::AuxiliaryEquationSet_SchurComplement<EvalT>::
buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm,
				      const panzer::FieldLibrary& /* field_library */,
				      const Teuchos::ParameterList& /* user_data */) const
{
  using panzer::BasisIRLayout;
  using panzer::EvaluatorStyle;
  using panzer::IntegrationRule;
  using panzer::Integrator_CurlBasisDotVector;
  using panzer::Traits;
  using PHX::Evaluator;
  using std::string;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  RCP<IntegrationRule> ir = this->getIntRuleForDOF(dof_name);
  RCP<BasisIRLayout> basis = this->getBasisIRLayoutForDOF(dof_name);

  std::vector<std::string> residual_operator_names;
  {
    {
      std::string resid = "AUX_SCHURCOMPLEMENT_RESIDUAL_TIME_OP_"+dof_name;
      ParameterList p("Time Derivative " + dof_name);
      p.set("Residual Name", resid);
      p.set("Value Name", dof_name);
      p.set("Basis", basis);
      p.set("IR", ir);
      p.set("Multiplier", 1.0);
      const std::vector<std::string> fieldMultiplier = {"1/dt",permittivity_};
      p.set("Field Multipliers",Teuchos::rcpFromRef(fieldMultiplier));
      RCP<Evaluator<Traits> > op = rcp(new
        panzer::Integrator_BasisTimesVector<EvalT, Traits>(p));
      fm.template registerEvaluator<EvalT>(op);
      residual_operator_names.push_back(resid);
    }
  }
  {
      std::string resid="AUX_SCHURCOMPLEMENT_RESIDUAL_CONDUCTIVITY_"+dof_name;
      ParameterList p("Conductivity "+dof_name);
      p.set("Residual Name", resid);
      p.set("Value Name", dof_name);
      p.set("Basis", basis);
      p.set("IR", ir);
      p.set("Tensor Name", conductivity_);
      RCP<PHX::Evaluator<panzer::Traits> > op =
        rcp(new panzer::Integrator_BasisTimesTensorTimesVector<EvalT, panzer::Traits>(p));

      this->template registerEvaluator<EvalT>(fm, op);
      residual_operator_names.push_back(resid);
    }
  // Curl curl Operator
  {
    {
      std::string resid = "AUX_SCHURCOMPLEMENT_RESIDUAL_CURLCURL_"+dof_name;
      ParameterList p("Curl Curl " + dof_name);
      p.set("Residual Name", resid);
      p.set("Value Name", "Curl_"+dof_name);
      p.set("Basis", basis);
      p.set("IR", ir);
      p.set("Multiplier", 1.0);
      const std::vector<std::string> fieldMultiplier = {"dt", inversePermeability_};
      p.set("Field Multipliers", Teuchos::rcpFromRef(fieldMultiplier));
      RCP<Evaluator<Traits> > op = rcp(new
        Integrator_CurlBasisDotVector<EvalT, Traits>(p));
      fm.template registerEvaluator<EvalT>(op);
      residual_operator_names.push_back(resid);
    }
  }

  const std::string residualField = "AUX_SCHURCOMPLEMENT_RESIDUAL_"+dof_name;
  this->buildAndRegisterResidualSummationEvaluator(fm,dof_name,residual_operator_names, residualField);
}

// ***********************************************************************
template <typename EvalT >
void mini_em::AuxiliaryEquationSet_SchurComplement<EvalT>::
buildAndRegisterScatterEvaluators(PHX::FieldManager<panzer::Traits>& /* fm */,
				  const panzer::FieldLibrary& /* field_library */,
                                  const panzer::LinearObjFactory<panzer::Traits> & /* lof */,
                                  const Teuchos::ParameterList& /* user_data */) const
{
}

// ***********************************************************************
template < >
void mini_em::AuxiliaryEquationSet_SchurComplement<panzer::Traits::Jacobian>::
buildAndRegisterScatterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
				  const panzer::FieldLibrary& field_library,
                                  const panzer::LinearObjFactory<panzer::Traits> & lof,
                                  const Teuchos::ParameterList& /* user_data */) const
{
   using Teuchos::RCP;
   using Teuchos::rcp;

   typedef panzer::Traits::Jacobian EvalT;

   typedef double Scalar;
#ifdef PANZER_HAVE_EPETRA_STACK
   typedef int LocalOrdinalEpetra;
#endif
   typedef int LocalOrdinalTpetra;
   typedef panzer::GlobalOrdinal GlobalOrdinalTpetra;

   typedef typename panzer::BlockedTpetraLinearObjFactory<panzer::Traits,Scalar,LocalOrdinalTpetra,GlobalOrdinalTpetra> blockedTpetraLinObjFactory;
   typedef typename panzer::TpetraLinearObjFactory<panzer::Traits,Scalar,LocalOrdinalTpetra,GlobalOrdinalTpetra> tpetraLinObjFactory;
#ifdef PANZER_HAVE_EPETRA_STACK
   typedef typename panzer::BlockedEpetraLinearObjFactory<panzer::Traits,LocalOrdinalEpetra> blockedEpetraLinObjFactory;
   typedef typename panzer::BlockedEpetraLinearObjFactory<panzer::Traits,LocalOrdinalEpetra> epetraLinObjFactory;
#endif

   PANZER_FUNC_TIME_MONITOR_DIFF("mini_em::AuxEqSet_SchurComplement::buildAndRegisterScatterEvaluators()",scatter_eval);

   std::string fieldStr = (*this->m_dof_names)[0];
   const std::string residualField = "AUX_SCHURCOMPLEMENT_RESIDUAL_"+dof_name;
   int pFieldNum;
   int blockIndex;

   std::string outPrefix = "ScatterReordered_";

   Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > nlof;

   // must be able to cast to a block linear object factory
   Teuchos::RCP<const blockedTpetraLinObjFactory> tblof
      = Teuchos::rcp_dynamic_cast<const blockedTpetraLinObjFactory>(Teuchos::rcpFromRef(lof));
#ifdef PANZER_HAVE_EPETRA_STACK
   Teuchos::RCP<const blockedEpetraLinObjFactory> eblof
      = Teuchos::rcp_dynamic_cast<const blockedEpetraLinObjFactory>(Teuchos::rcpFromRef(lof));
#endif

   if(tblof != Teuchos::null) {
     Teuchos::RCP<const panzer::BlockedDOFManager> blockedDOFMngr = tblof->getGlobalIndexer();
     TEUCHOS_ASSERT(blockedDOFMngr!=Teuchos::null);

     pFieldNum = blockedDOFMngr->getFieldNum(fieldStr);
     blockIndex = blockedDOFMngr->getFieldBlock(pFieldNum);

     // get the unique global indexer for just this field
     Teuchos::RCP<panzer::GlobalIndexer> ugi = blockedDOFMngr->getFieldDOFManagers()[blockIndex];

     // build a new linear object factory
     nlof = Teuchos::rcp(new tpetraLinObjFactory(Teuchos::rcp(new Teuchos::MpiComm<int>(tblof->getComm())).getConst(),ugi));

     // first build a reordering evaluator to take it to the new sub global indexer
     {
        std::vector<Teuchos::RCP<PHX::DataLayout> > fieldLayouts;
        fieldLayouts.push_back(field_library.lookupBasis(fieldStr)->functional);

        std::vector<std::string> resNames;
        resNames.push_back(residualField);

        RCP< PHX::Evaluator<panzer::Traits> > op = Teuchos::rcp(
              new panzer::ReorderADValues_Evaluator<EvalT,panzer::Traits>(outPrefix,
                                                    resNames,
                                                    fieldLayouts,
                                                    this->getElementBlockId(),
                                                    *blockedDOFMngr,
                                                    *ugi));

        fm.registerEvaluator<EvalT>(op);
     }
#ifdef PANZER_HAVE_EPETRA_STACK
   } else if(eblof != Teuchos::null) {
     Teuchos::RCP<const panzer::BlockedDOFManager> blockedDOFMngr = eblof->getGlobalIndexer();
     TEUCHOS_ASSERT(blockedDOFMngr!=Teuchos::null);

     pFieldNum = blockedDOFMngr->getFieldNum(fieldStr);
     blockIndex = blockedDOFMngr->getFieldBlock(pFieldNum);

     // get the unique global indexer for just this field
     Teuchos::RCP<panzer::GlobalIndexer> ugi = blockedDOFMngr->getFieldDOFManagers()[blockIndex];

     // build a new linear object factory
     nlof = Teuchos::rcp(new epetraLinObjFactory(Teuchos::rcp(new Teuchos::MpiComm<int>(eblof->getComm())).getConst(),ugi));

     // first build a reordering evaluator to take it to the new sub global indexer
     {
        std::vector<Teuchos::RCP<PHX::DataLayout> > fieldLayouts;
        fieldLayouts.push_back(field_library.lookupBasis(fieldStr)->functional);

        std::vector<std::string> resNames;
        resNames.push_back(residualField);

        RCP< PHX::Evaluator<panzer::Traits> > op = Teuchos::rcp(
              new panzer::ReorderADValues_Evaluator<EvalT,panzer::Traits>(outPrefix,
                                                    resNames,
                                                    fieldLayouts,
                                                    this->getElementBlockId(),
                                                    *blockedDOFMngr,
                                                    *ugi));

        fm.registerEvaluator<EvalT>(op);
     }
#endif
   } else
     TEUCHOS_ASSERT(false);

   {
      RCP<std::map<std::string,std::string> > resToField
         = rcp(new std::map<std::string,std::string>);
      (*resToField)[outPrefix+residualField] = dof_name;

      RCP<std::vector<std::string> > inFieldNames
         = rcp(new std::vector<std::string>);
      inFieldNames->push_back(outPrefix+residualField);

      std::string scatterName = "AUX_"+dof_name+"_SchurComplement";

      Teuchos::ParameterList p("Scatter:" + scatterName);
      p.set("Scatter Name", scatterName);
      p.set("Basis",field_library.lookupBasis(fieldStr));
      p.set("Dependent Names", inFieldNames);
      p.set("Dependent Map", resToField);
      p.set("Global Data Key", "SchurComplement " + dof_name + " Scatter Container");

      RCP< PHX::Evaluator<panzer::Traits> > op = nlof->buildScatter<EvalT>(p);

      fm.registerEvaluator<EvalT>(op);

      PHX::Tag<EvalT::ScalarT> tag(scatterName, Teuchos::rcp(new PHX::MDALayout<panzer::Dummy>(0)));
      fm.requireField<EvalT>(tag);
   }

   if(!m_gedc->containsDataObject("SchurComplement " + dof_name + " Scatter Container")) {
      Teuchos::RCP<panzer::GlobalEvaluationData> dataObject
         = Teuchos::rcp(new panzer::LOCPair_GlobalEvaluationData(nlof,panzer::LinearObjContainer::Mat));
      m_gedc->addDataObject("SchurComplement " + dof_name + " Scatter Container",dataObject);
   }
}

// ***********************************************************************
#endif
