// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _MiniEM_AuxiliaryEquationSet_MassMatrix_impl_hpp_
#define _MiniEM_AuxiliaryEquationSet_MassMatrix_impl_hpp_

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
#include "Panzer_Integrator_BasisTimesScalar.hpp"
#include "Panzer_Integrator_BasisTimesVector.hpp"
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
mini_em::AuxiliaryEquationSet_MassMatrix<EvalT>::
AuxiliaryEquationSet_MassMatrix(
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
    valid_parameters.set("Multiplier",1.0,"Scale the operator");
    valid_parameters.set("Field Multipliers","","Scale the operator");
    valid_parameters.set("Basis Type","HGrad","Type of Basis to use");
    valid_parameters.set("Basis Order",1,"Order of the basis");
    valid_parameters.set("Integration Order",default_integration_order,"Order of the integration rule");
    valid_parameters.set("Operator Label","","Optional label for the operator");

    params->validateParametersAndSetDefaults(valid_parameters);
  }
  std::string model_id = params->get<std::string>("Model ID");
  dof_name = params->get<std::string>("DOF Name");
  multiplier = params->get<double>("Multiplier");
  std::stringstream ss(params->get<std::string>("Field Multipliers"));
  std::string item;
  Teuchos::RCP<std::vector<std::string> > fieldMultipliersNonConst = Teuchos::rcp(new std::vector<std::string>);
  panzer::StringTokenizer(*fieldMultipliersNonConst,params->get<std::string>("Field Multipliers"),",",true);
  // while (std::getline(ss, item, ','))
  //   fieldMultipliersNonConst->push_back(item);
  fieldMultipliers = fieldMultipliersNonConst.getConst();
  std::string basis_type = params->get<std::string>("Basis Type");
  int basis_order = params->get<int>("Basis Order");
  int integration_order = params->get<int>("Integration Order");
  opLabel = params->get<std::string>("Operator Label");

  // ********************
  // Assemble DOF names and Residual names
  // ********************

  m_gedc = gedc;

  m_dof_names = Teuchos::rcp(new std::vector<std::string>);
  m_dof_names->push_back(dof_name);

  this->addDOF(dof_name,basis_type,basis_order,integration_order,"AUX_MASS_RESIDUAL_"+dof_name);
  this->addDOFGrad(dof_name,"GRAD_"+dof_name);

  this->addClosureModel(model_id);

  this->setupDOFs();
}

// ***********************************************************************
template <typename EvalT>
void mini_em::AuxiliaryEquationSet_MassMatrix<EvalT>::
buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm,
				      const panzer::FieldLibrary& /* field_library */,
				      const Teuchos::ParameterList& /* user_data */) const
{
  using panzer::BasisIRLayout;
  using panzer::EvaluatorStyle;
  using panzer::IntegrationRule;
  using panzer::Integrator_BasisTimesScalar;
  using panzer::Integrator_BasisTimesVector;
  using panzer::Traits;
  using PHX::Evaluator;
  using std::string;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  
  RCP<IntegrationRule> ir = this->getIntRuleForDOF(dof_name); 
  RCP<BasisIRLayout> basis = this->getBasisIRLayoutForDOF(dof_name); 

  // Mass Operator
  {
    if (basis->getBasis()->isScalarBasis())
    {
      ParameterList p("Mass Matrix " + dof_name + " Residual");
      p.set("Residual Name", "AUX_MASS_RESIDUAL_"+opLabel+dof_name);
      p.set("Value Name", dof_name);
      p.set("Basis", basis);
      p.set("IR", ir);
      p.set("Multiplier", multiplier);
      if (fieldMultipliers != Teuchos::null)
        p.set("Field Multipliers", fieldMultipliers);
      RCP<Evaluator<Traits>> op = rcp(new
        Integrator_BasisTimesScalar<EvalT, Traits>(p));
      fm.template registerEvaluator<EvalT>(op);
    }
    else if (basis->getBasis()->isVectorBasis())
    { 
      ParameterList p("Mass Matrix " + dof_name + " Residual");
      p.set("Residual Name", "AUX_MASS_RESIDUAL_"+opLabel+dof_name);
      p.set("Value Name", dof_name);
      p.set("Basis", basis);
      p.set("IR", ir);
      p.set("Multiplier", multiplier);
      if (fieldMultipliers != Teuchos::null)
        p.set("Field Multipliers", fieldMultipliers);
      RCP<Evaluator<Traits>> op = rcp(new
        Integrator_BasisTimesVector<EvalT, Traits>(p));
      fm.template registerEvaluator<EvalT>(op);
    }
    else
      TEUCHOS_ASSERT(false);
  }
}

// ***********************************************************************
template <typename EvalT >
void mini_em::AuxiliaryEquationSet_MassMatrix<EvalT>::
buildAndRegisterScatterEvaluators(PHX::FieldManager<panzer::Traits>& /* fm */,
				  const panzer::FieldLibrary& /* field_library */,
                                  const panzer::LinearObjFactory<panzer::Traits> & /* lof */,
                                  const Teuchos::ParameterList& /* user_data */) const
{
}

// ***********************************************************************
template < >
void mini_em::AuxiliaryEquationSet_MassMatrix<panzer::Traits::Jacobian>::
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

   std::string fieldStr = (*this->m_dof_names)[0];
   const std::string residualField = "AUX_MASS_RESIDUAL_"+opLabel+dof_name;
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

   std::string scatterContainerName = "Mass Matrix " + opLabel + dof_name + " Scatter Container";
   {
      RCP<std::map<std::string,std::string> > resToField 
         = rcp(new std::map<std::string,std::string>);
      (*resToField)[outPrefix+residualField] = dof_name;

      RCP<std::vector<std::string> > inFieldNames
         = rcp(new std::vector<std::string>);
      inFieldNames->push_back(outPrefix+residualField);

      std::string scatterName = "AUX_"+ opLabel + dof_name + "_MassMatrix";

      Teuchos::ParameterList p("Scatter:" + scatterName);
      p.set("Scatter Name", scatterName);
      p.set("Basis",field_library.lookupBasis(fieldStr));
      p.set("Dependent Names", inFieldNames);
      p.set("Dependent Map", resToField);
      p.set("Global Data Key", scatterContainerName);

      RCP< PHX::Evaluator<panzer::Traits> > op = nlof->buildScatter<EvalT>(p);

      fm.registerEvaluator<EvalT>(op);

      PHX::Tag<EvalT::ScalarT> tag(scatterName, Teuchos::rcp(new PHX::MDALayout<panzer::Dummy>(0)));
      fm.requireField<EvalT>(tag);
   }
   
   if(!m_gedc->containsDataObject(scatterContainerName)) {
      Teuchos::RCP<panzer::GlobalEvaluationData> dataObject 
         = Teuchos::rcp(new panzer::LOCPair_GlobalEvaluationData(nlof,panzer::LinearObjContainer::Mat));
      m_gedc->addDataObject(scatterContainerName,dataObject);
   }
}

// ***********************************************************************
#endif
