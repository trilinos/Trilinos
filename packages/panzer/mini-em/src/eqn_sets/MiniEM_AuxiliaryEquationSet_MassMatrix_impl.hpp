#ifndef _MiniEM_AuxiliaryEquationSet_MassMatrix_impl_hpp_
#define _MiniEM_AuxiliaryEquationSet_MassMatrix_impl_hpp_

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"

#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_FieldManager.hpp"

#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"

// include evaluators here
#include "Panzer_Integrator_BasisTimesScalar.hpp"
#include "Panzer_Integrator_BasisTimesVector.hpp"
#include "Panzer_Integrator_GradBasisDotVector.hpp"
#include "Panzer_ScalarToVector.hpp"
#include "Panzer_Sum.hpp"
#include "Panzer_Constant.hpp"
#include "Panzer_BlockedDOFManager.hpp"
#include "Panzer_BlockedTpetraLinearObjFactory.hpp"
#include "Panzer_TpetraLinearObjFactory.hpp"
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
    valid_parameters.set("DOF Name","","Name of DOF to construct time derivative for");
    valid_parameters.set("Multiplier",1.0,"Scale the operator");
    valid_parameters.set("Basis Type","HGrad","Type of Basis to use");
    valid_parameters.set("Basis Order",1,"Order of the basis");
    valid_parameters.set("Integration Order",default_integration_order,"Order of the integration rule");

    params->validateParametersAndSetDefaults(valid_parameters);
  }
  dof_name = params->get<std::string>("DOF Name");
  multiplier = params->get<double>("Multiplier");
  std::string basis_type = params->get<std::string>("Basis Type");
  int basis_order = params->get<int>("Basis Order");
  int integration_order = params->get<int>("Integration Order");

  // ********************
  // Assemble DOF names and Residual names
  // ********************

  m_gedc = gedc;

  m_dof_names = Teuchos::rcp(new std::vector<std::string>);
  m_dof_names->push_back(dof_name);

  this->addDOF(dof_name,basis_type,basis_order,integration_order,"AUX_MASS_RESIDUAL_"+dof_name);
  this->addDOFGrad(dof_name,"GRAD_"+dof_name);

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
      string resName("AUX_MASS_RESIDUAL_" + dof_name), valName(dof_name);
      RCP<Evaluator<Traits>> op = rcp(new
        Integrator_BasisTimesScalar<EvalT, Traits>(EvaluatorStyle::EVALUATES,
        resName, valName, *basis, *ir, multiplier));
      fm.template registerEvaluator<EvalT>(op);
    }
    else if (basis->getBasis()->isVectorBasis())
    { 
      ParameterList p("Mass Matrix " + dof_name + " Residual");
      p.set("Residual Name", "AUX_MASS_RESIDUAL_"+dof_name);
      p.set("Value Name", dof_name);
      p.set("Basis", basis);
      p.set("IR", ir);
      p.set("Multiplier", multiplier);
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

   std::string fieldStr = (*this->m_dof_names)[0];
   const std::string residualField = "AUX_MASS_RESIDUAL_"+dof_name; 
   int pFieldNum;
   int blockIndex;

   // must be able to cast to a block linear object factory 
   Teuchos::RCP<const panzer::BlockedTpetraLinearObjFactory<panzer::Traits,double,int,panzer::Ordinal64> > tblof
      = Teuchos::rcp_dynamic_cast<const panzer::BlockedTpetraLinearObjFactory<panzer::Traits,double,int,panzer::Ordinal64> >(Teuchos::rcpFromRef(lof)); 

   std::string outPrefix = "ScatterReordered_";


   Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > nlof;

   if(tblof != Teuchos::null){
     Teuchos::RCP<const panzer::BlockedDOFManager<int,panzer::Ordinal64> > blockedDOFMngr = tblof->getGlobalIndexer();
     TEUCHOS_ASSERT(blockedDOFMngr!=Teuchos::null); 

     pFieldNum = blockedDOFMngr->getFieldNum(fieldStr);
     blockIndex = blockedDOFMngr->getFieldBlock(pFieldNum);

     // get the unique global indexer for just this field
     Teuchos::RCP<panzer::UniqueGlobalIndexer<int,panzer::Ordinal64> > ugi = blockedDOFMngr->getFieldDOFManagers()[blockIndex];
 
     // build a new tpetra linear object factory 
     nlof = Teuchos::rcp(new panzer::TpetraLinearObjFactory<panzer::Traits,double,int,panzer::Ordinal64>(Teuchos::rcp(new Teuchos::MpiComm<int>(tblof->getComm())).getConst(),ugi));

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

   } else
     TEUCHOS_ASSERT(false); 


   {
      RCP<std::map<std::string,std::string> > resToField 
         = rcp(new std::map<std::string,std::string>);
      (*resToField)[outPrefix+residualField] = dof_name;

      RCP<std::vector<std::string> > inFieldNames
         = rcp(new std::vector<std::string>);
      inFieldNames->push_back(outPrefix+residualField);

      std::string scatterName = "AUX_"+dof_name+"_MassMatrix";

      Teuchos::ParameterList p("Scatter:" + scatterName);
      p.set("Scatter Name", scatterName);
      p.set("Basis",field_library.lookupBasis(fieldStr));
      p.set("Dependent Names", inFieldNames);
      p.set("Dependent Map", resToField);
      p.set("Global Data Key", "Mass Matrix " + dof_name + " Scatter Container");

      RCP< PHX::Evaluator<panzer::Traits> > op = nlof->buildScatter<EvalT>(p);

      fm.registerEvaluator<EvalT>(op);

      PHX::Tag<EvalT::ScalarT> tag(scatterName, Teuchos::rcp(new PHX::MDALayout<panzer::Dummy>(0)));
      fm.requireField<EvalT>(tag);
   }
   
   if(!m_gedc->containsDataObject("Mass Matrix " + dof_name + " Scatter Container")) {
      Teuchos::RCP<panzer::GlobalEvaluationData> dataObject 
         = Teuchos::rcp(new panzer::LOCPair_GlobalEvaluationData(nlof,panzer::LinearObjContainer::Mat));
      m_gedc->addDataObject("Mass Matrix " + dof_name + " Scatter Container",dataObject);
   }
}

// ***********************************************************************
#endif
