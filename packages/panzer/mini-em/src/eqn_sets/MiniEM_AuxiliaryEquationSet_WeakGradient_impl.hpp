// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _MiniEM_AuxiliaryEquationSet_WeakGradient_impl_hpp_
#define _MiniEM_AuxiliaryEquationSet_WeakGradient_impl_hpp_

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"

#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_FieldManager.hpp"

#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"

// include evaluators here
#include "Panzer_Integrator_BasisTimesVector.hpp"
#include "Panzer_ScalarToVector.hpp"
#include "Panzer_Sum.hpp"
#include "Panzer_Constant.hpp"
#include "Panzer_BlockedDOFManager.hpp"
#ifdef PANZER_HAVE_EPETRA_STACK
#include "Panzer_BlockedEpetraLinearObjFactory.hpp"
#endif
#include "Panzer_BlockedTpetraLinearObjFactory.hpp"
#include "Panzer_TpetraLinearObjFactory.hpp"
#include "Panzer_ReorderADValues_Evaluator.hpp"
#include "Panzer_LOCPair_GlobalEvaluationData.hpp"

// ***********************************************************************
template <typename EvalT>
mini_em::AuxiliaryEquationSet_WeakGradient<EvalT>::
AuxiliaryEquationSet_WeakGradient(
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
    valid_parameters.set("Vector Name","","Name of vector field for domain space");
    valid_parameters.set("Scalar Name","","Name of scalar field for range space");
    valid_parameters.set("Basis Type","HCurl","Type of Basis to use");
    valid_parameters.set("Basis Order",1,"Order of the basis");
    valid_parameters.set("Integration Order",default_integration_order,"Order of the integration rule");

    params->validateParametersAndSetDefaults(valid_parameters);
  }
  vector_name = params->get<std::string>("Vector Name");
  scalar_name = params->get<std::string>("Scalar Name");
  std::string basis_type = params->get<std::string>("Basis Type");
  int basis_order = params->get<int>("Basis Order");
  int integration_order = params->get<int>("Integration Order");

  // ********************
  // Assemble DOF names and Residual names
  // ********************
  m_gedc = gedc;
  m_dof_names = Teuchos::rcp(new std::vector<std::string>);
  m_dof_names->push_back(vector_name);

  this->addDOF(vector_name,basis_type,basis_order,integration_order,"AUX_RESIDUAL_"+vector_name);
  
  this->setupDOFs();
}

// ***********************************************************************
template <typename EvalT>
void mini_em::AuxiliaryEquationSet_WeakGradient<EvalT>::
buildAndRegisterEquationSetEvaluators(PHX::FieldManager<panzer::Traits>& fm,
				      const panzer::FieldLibrary& /* field_library */,
				      const Teuchos::ParameterList& /* user_data */) const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  
  Teuchos::RCP<panzer::IntegrationRule> ir = this->getIntRuleForDOF(vector_name); 
  Teuchos::RCP<panzer::BasisIRLayout> basis = this->getBasisIRLayoutForDOF(vector_name); 

  // ********************
  // Gradient Equation
  // ********************

  {
      ParameterList p("Weak Gradient Residual");
      p.set("Residual Name", "AUX_RESIDUAL_"+vector_name);
      p.set("Value Name", "GRAD_"+scalar_name);
      p.set("Basis", basis);
      p.set("IR", ir);
      p.set("Multiplier", 1.0);
      
      RCP< PHX::Evaluator<panzer::Traits> > op = 
	rcp(new panzer::Integrator_BasisTimesVector<EvalT,panzer::Traits>(p));
      
      fm.template registerEvaluator<EvalT>(op);
  
  }

}

// ***********************************************************************
template <typename EvalT >
void mini_em::AuxiliaryEquationSet_WeakGradient<EvalT>::
buildAndRegisterScatterEvaluators(PHX::FieldManager<panzer::Traits>& /* fm */,
				  const panzer::FieldLibrary& /* field_library */,
                                  const panzer::LinearObjFactory<panzer::Traits> & /* lof */,
                                  const Teuchos::ParameterList& /* user_data */) const
{
}

// ***********************************************************************
template < >
void mini_em::AuxiliaryEquationSet_WeakGradient<panzer::Traits::Jacobian>::
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

   PANZER_FUNC_TIME_MONITOR_DIFF("mini_em::AuxEqSet_WeakGradient::buildAndRegisterScatterEvaluators()",scatter_eval);

   std::string fieldStr = (*this->m_dof_names)[0];
   int pFieldNum;
   int uFieldNum;
   int rowBlockIndex;
   int colBlockIndex;

   const std::string outPrefix = "ScatterReordered_";
   const std::string residualField = "AUX_RESIDUAL_"+vector_name; 
   const std::string scatterName = "AUX_"+vector_name+"_WeakGradient";

   Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > nlof;

   // must be able to cast to a block linear object factory
   Teuchos::RCP<const blockedTpetraLinObjFactory> tblof
      = Teuchos::rcp_dynamic_cast<const blockedTpetraLinObjFactory>(Teuchos::rcpFromRef(lof));
#ifdef PANZER_HAVE_EPETRA_STACK
   Teuchos::RCP<const blockedEpetraLinObjFactory> eblof
      = Teuchos::rcp_dynamic_cast<const blockedEpetraLinObjFactory>(Teuchos::rcpFromRef(lof));
#endif

   if(tblof != Teuchos::null) {

     Teuchos::RCP<const panzer::BlockedDOFManager> blockedDOFMngr;
     Teuchos::RCP<panzer::GlobalIndexer> rowUgi;
     Teuchos::RCP<panzer::GlobalIndexer> colUgi;

     blockedDOFMngr = tblof->getGlobalIndexer();
     TEUCHOS_ASSERT(blockedDOFMngr!=Teuchos::null); 

     uFieldNum = blockedDOFMngr->getFieldNum(fieldStr);
     pFieldNum = blockedDOFMngr->getFieldNum(scalar_name);
     rowBlockIndex = blockedDOFMngr->getFieldBlock(uFieldNum);
     colBlockIndex = blockedDOFMngr->getFieldBlock(pFieldNum);

     // get the unique global indexer for just this field
     rowUgi = blockedDOFMngr->getFieldDOFManagers()[rowBlockIndex];
     colUgi = blockedDOFMngr->getFieldDOFManagers()[colBlockIndex];
   
     // build a new tpetra linear object factory
     nlof = Teuchos::rcp(new tpetraLinObjFactory(Teuchos::rcp(new Teuchos::MpiComm<int>(tblof->getComm())).getConst(),rowUgi,colUgi));

     // first build a reordering evaluator to take it to the new sub global indexer
     {
       std::vector<Teuchos::RCP<PHX::DataLayout> > fieldLayouts;
       fieldLayouts.push_back(field_library.lookupBasis(fieldStr)->functional);

       std::vector<std::string> resNames;
       resNames.push_back(this->m_provided_dofs_desc.begin()->second.residualName.second);

       RCP< PHX::Evaluator<panzer::Traits> > op = Teuchos::rcp(
                                                               new panzer::ReorderADValues_Evaluator<EvalT,panzer::Traits>(outPrefix,
                                                                                                                           resNames,
                                                                                                                           fieldLayouts,
                                                                                                                           this->getElementBlockId(),
                                                                                                                           *blockedDOFMngr,
                                                                                                                           *colUgi));

       fm.registerEvaluator<EvalT>(op);
     }
#ifdef PANZER_HAVE_EPETRA_STACK
   } else if(eblof != Teuchos::null) {
     Teuchos::RCP<const panzer::BlockedDOFManager> blockedDOFMngr;
     Teuchos::RCP<panzer::GlobalIndexer> rowUgi;
     Teuchos::RCP<panzer::GlobalIndexer> colUgi;

     blockedDOFMngr = eblof->getGlobalIndexer();
     TEUCHOS_ASSERT(blockedDOFMngr!=Teuchos::null);

     uFieldNum = blockedDOFMngr->getFieldNum(fieldStr);
     pFieldNum = blockedDOFMngr->getFieldNum(scalar_name);
     rowBlockIndex = blockedDOFMngr->getFieldBlock(uFieldNum);
     colBlockIndex = blockedDOFMngr->getFieldBlock(pFieldNum);

     // get the unique global indexer for just this field
     rowUgi = blockedDOFMngr->getFieldDOFManagers()[rowBlockIndex];
     colUgi = blockedDOFMngr->getFieldDOFManagers()[colBlockIndex];

     // build a new epetra linear object factory
     nlof = Teuchos::rcp(new epetraLinObjFactory(Teuchos::rcp(new Teuchos::MpiComm<int>(eblof->getComm())).getConst(),rowUgi,colUgi));

     // first build a reordering evaluator to take it to the new sub global indexer
     {
       std::vector<Teuchos::RCP<PHX::DataLayout> > fieldLayouts;
       fieldLayouts.push_back(field_library.lookupBasis(fieldStr)->functional);

       std::vector<std::string> resNames;
       resNames.push_back(this->m_provided_dofs_desc.begin()->second.residualName.second);

       RCP< PHX::Evaluator<panzer::Traits> > op = Teuchos::rcp(
                                                               new panzer::ReorderADValues_Evaluator<EvalT,panzer::Traits>(outPrefix,
                                                                                                                           resNames,
                                                                                                                           fieldLayouts,
                                                                                                                           this->getElementBlockId(),
                                                                                                                           *blockedDOFMngr,
                                                                                                                           *colUgi));

       fm.registerEvaluator<EvalT>(op);
     }
#endif
   } else
     TEUCHOS_ASSERT(false);

   {
     RCP<std::map<std::string,std::string> > resToField 
       = rcp(new std::map<std::string,std::string>);
     (*resToField)[outPrefix+residualField] = vector_name;
   
     RCP<std::vector<std::string> > inFieldNames
       = rcp(new std::vector<std::string>);
     inFieldNames->push_back(outPrefix+residualField);
   
     Teuchos::ParameterList p("Scatter");
     p.set("Scatter Name", scatterName);
     p.set("Basis",field_library.lookupBasis(fieldStr));
     p.set("Dependent Names", inFieldNames);
     p.set("Dependent Map", resToField);
     p.set("Global Data Key", "Weak Gradient Scatter Container");
   
     RCP< PHX::Evaluator<panzer::Traits> > op = nlof->buildScatter<EvalT>(p);
   
     fm.registerEvaluator<EvalT>(op);
   
     PHX::Tag<EvalT::ScalarT> tag(scatterName, Teuchos::rcp(new PHX::MDALayout<panzer::Dummy>(0)));
     fm.requireField<EvalT>(tag);
   }
   
   if(!m_gedc->containsDataObject("Weak Gradient Scatter Container")) {
      Teuchos::RCP<panzer::GlobalEvaluationData> dataObject 
         = Teuchos::rcp(new panzer::LOCPair_GlobalEvaluationData(nlof,panzer::LinearObjContainer::Mat));
      m_gedc->addDataObject("Weak Gradient Scatter Container",dataObject);
   }
}


// ***********************************************************************
#endif
