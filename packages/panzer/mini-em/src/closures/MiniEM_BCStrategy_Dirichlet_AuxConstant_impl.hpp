// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MiniEM_BC_STRATEGY_DIRICHLET_AuxCONSTANT_IMPL_HPP
#define MiniEM_BC_STRATEGY_DIRICHLET_AuxCONSTANT_IMPL_HPP

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_FieldManager.hpp"
#include "Panzer_PhysicsBlock.hpp"

#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Dirichlet_Residual.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_BlockedDOFManager.hpp"
#include "Panzer_BlockedTpetraLinearObjFactory.hpp"
#include "Panzer_TpetraLinearObjFactory.hpp"
#include "Panzer_ReorderADValues_Evaluator.hpp"

// Evaluators
#include "Panzer_Constant.hpp"

// ***********************************************************************
template <typename EvalT>
mini_em::BCStrategy_Dirichlet_AuxConstant<EvalT>::
BCStrategy_Dirichlet_AuxConstant(const panzer::BC& bc,
    const Teuchos::RCP<panzer::GlobalData>& /* global_data */) 
  : panzer::BCStrategy<EvalT>(bc)
{
  TEUCHOS_ASSERT(this->m_bc.strategy() == "AuxConstant");
}

// ***********************************************************************
template <typename EvalT>
void mini_em::BCStrategy_Dirichlet_AuxConstant<EvalT>::
setup(const panzer::PhysicsBlock& side_pb,
      const Teuchos::ParameterList& /* user_data */)
{
  using Teuchos::RCP;
  using std::vector;
  using std::string;
  using std::pair;

  operatorName_ = this->m_bc.equationSetName();
  value_ = this->m_bc.params()->template get<double>("Value");
  fieldName_ = this->m_bc.params()->template get<std::string>("Field Name");

  // find the basis for this dof 
  const std::vector<std::pair<std::string,RCP<panzer::PureBasis> > >& dofs = side_pb.getProvidedDOFs();

  for (std::vector<std::pair<std::string,RCP<panzer::PureBasis> > >::const_iterator dof_it = 
	 dofs.begin(); dof_it != dofs.end(); ++dof_it) {
    if (dof_it->first == fieldName_)
      this->basis_ = dof_it->second;
  }

  TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::is_null(this->basis_), std::runtime_error,
		     "Error the name \"" << this->m_bc.equationSetName()
		     << "\" is not a valid DOF for the boundary condition:\n"
		     << this->m_bc << "\n");
}

// ***********************************************************************
template <typename EvalT>
void mini_em::BCStrategy_Dirichlet_AuxConstant<EvalT>::
buildAndRegisterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
			   const panzer::PhysicsBlock& /* pb */,
			   const panzer::ClosureModelFactory_TemplateManager<panzer::Traits>& /* factory */,
			   const Teuchos::ParameterList& /* models */,
			   const Teuchos::ParameterList& /* user_data */) const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // Evaluator for Constant dirichlet BCs
  // provide a constant target value to map into residual
  {
    ParameterList p("BC Constant Dirichlet");
    p.set("Name", "AuxConstant_" + this->m_bc.equationSetName()+"_"+fieldName_);
    p.set("Data Layout", basis_->functional);
    p.set("Value", this->m_bc.params()->template get<double>("Value"));
    
    RCP< PHX::Evaluator<panzer::Traits> > op = 
      rcp(new panzer::Constant<EvalT,panzer::Traits>(p));
    
    fm.template registerEvaluator<EvalT>(op);
  }

  {
    Teuchos::ParameterList p("Dirichlet Residual: "+operatorName_ + " " +fieldName_);
    p.set("Residual Name", "Residual_AuxConstant_" + this->m_bc.equationSetName()+"_"+fieldName_);
    p.set("DOF Name", fieldName_);
    p.set("Value Name", "AuxConstant_" + this->m_bc.equationSetName()+"_"+fieldName_);
    p.set("Data Layout", basis_->functional);

    RCP< PHX::Evaluator<panzer::Traits> > op = 
      rcp(new panzer::DirichletResidual<EvalT,panzer::Traits>(p));
    
    fm.template registerEvaluator<EvalT>(op);
  }
}

template <typename EvalT>
void mini_em::BCStrategy_Dirichlet_AuxConstant<EvalT>::
buildAndRegisterGatherScatterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                        const panzer::PhysicsBlock& side_pb,
                                        const panzer::LinearObjFactory<panzer::Traits> & lof,
                                        const Teuchos::ParameterList& user_data) const {}

template <typename EvalT>
void mini_em::BCStrategy_Dirichlet_AuxConstant<EvalT>::
buildAndRegisterScatterEvaluators(PHX::FieldManager<panzer::Traits>& /* fm */,
                                  const panzer::PhysicsBlock& /* side_pb */,
				  const panzer::LinearObjFactory<panzer::Traits> & /* lof */,
				  const Teuchos::ParameterList& /* user_data */) const {}

template <typename EvalT>
void mini_em::BCStrategy_Dirichlet_AuxConstant<EvalT>::
buildAndRegisterGatherAndOrientationEvaluators(PHX::FieldManager<panzer::Traits>& /* fm */,
                                                 const panzer::PhysicsBlock& /* side_pb */,
                                                 const panzer::LinearObjFactory<panzer::Traits> & /* lof */,
                                                 const Teuchos::ParameterList& /* user_data */) const {}

template < >
void mini_em::BCStrategy_Dirichlet_AuxConstant<panzer::Traits::Jacobian>::
buildAndRegisterGatherScatterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                        const panzer::PhysicsBlock& side_pb,
                                        const panzer::LinearObjFactory<panzer::Traits> & lof,
                                        const Teuchos::ParameterList& user_data) const
{
  buildAndRegisterGatherAndOrientationEvaluators(fm,side_pb,lof,user_data);
  buildAndRegisterScatterEvaluators(fm,side_pb,lof,user_data);
}

template < >
void mini_em::BCStrategy_Dirichlet_AuxConstant<panzer::Traits::Jacobian>::
buildAndRegisterGatherAndOrientationEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                                 const panzer::PhysicsBlock& /* side_pb */,
                                                 const panzer::LinearObjFactory<panzer::Traits> & lof,
                                                 const Teuchos::ParameterList& /* user_data */) const 
{
  typedef panzer::Traits::Jacobian EvalT;

  using Teuchos::RCP;
  using Teuchos::rcp;
 
  // Gather
  {
    
    Teuchos::ParameterList p("BC Gather");
    
    RCP<std::vector<std::string> > gather_names_vec = rcp(new std::vector<std::string>);
    gather_names_vec->push_back(fieldName_);
    
    p.set("DOF Names", gather_names_vec);
    p.set("Indexer Names", gather_names_vec);
    p.set("Basis", basis_);
    
    RCP< PHX::Evaluator<panzer::Traits> > op = lof.buildGather<EvalT>(p);
    
    fm.registerEvaluator<EvalT>(op);
  }

}

template < >
void mini_em::BCStrategy_Dirichlet_AuxConstant<panzer::Traits::Jacobian>::
buildAndRegisterScatterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                  const panzer::PhysicsBlock& side_pb,
				  const panzer::LinearObjFactory<panzer::Traits> & lof,
				  const Teuchos::ParameterList& /* user_data */) const 
{
  typedef panzer::Traits::Jacobian EvalT;

  using Teuchos::RCP;
  using Teuchos::rcp;

  // must be able to cast to a block linear object factory
  const panzer::BlockedTpetraLinearObjFactory<panzer::Traits,double,int,panzer::GlobalOrdinal> & blof
     = Teuchos::dyn_cast<const panzer::BlockedTpetraLinearObjFactory<panzer::Traits,double,int,panzer::GlobalOrdinal> >(lof); 
  Teuchos::RCP<const panzer::BlockedDOFManager> blockedDOFMngr = blof.getGlobalIndexer();
  TEUCHOS_ASSERT(blockedDOFMngr!=Teuchos::null); 

  int fieldNum = blockedDOFMngr->getFieldNum(fieldName_);
  int blockIndex = blockedDOFMngr->getFieldBlock(fieldNum);

  // get the unique global indexer for just this field
  Teuchos::RCP<panzer::GlobalIndexer> ugi = blockedDOFMngr->getFieldDOFManagers()[blockIndex];
   
  // build a new epetra linear object factory 
  Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > elof
     = Teuchos::rcp(new panzer::TpetraLinearObjFactory<panzer::Traits,double,int,panzer::GlobalOrdinal>(Teuchos::rcp(new Teuchos::MpiComm<int>(blof.getComm())).getConst(),ugi));

  const std::string residualField = "Residual_AuxConstant_" + this->m_bc.equationSetName()+"_"+fieldName_;
  const std::string outPrefix = "ScatterReordered_";
  const std::string scatterName = "Scatter_" + residualField;

  // first build a reordering evaluator to take it to the new sub global indexer
  {
     std::vector<Teuchos::RCP<PHX::DataLayout> > fieldLayouts; 
     fieldLayouts.push_back(basis_->functional);

     std::vector<std::string> resNames; 
     resNames.push_back(residualField);

     RCP< PHX::Evaluator<panzer::Traits> > op = Teuchos::rcp(
           new panzer::ReorderADValues_Evaluator<EvalT,panzer::Traits>(outPrefix,
                                                 resNames,
                                                 fieldLayouts,
                                                 side_pb.elementBlockID(),
                                                 *blockedDOFMngr,
                                                 *ugi));

     fm.registerEvaluator<EvalT>(op);
  }

  {
    RCP<std::map<std::string,std::string> > resToField 
       = rcp(new std::map<std::string,std::string>);
    (*resToField)[outPrefix+residualField] = fieldName_;

    RCP<std::vector<std::string> > inFieldNames
       = rcp(new std::vector<std::string>);
    inFieldNames->push_back(outPrefix+residualField);

    Teuchos::ParameterList p("Scatter");
    p.set("Scatter Name", scatterName);
    p.set("Basis",basis_);
    p.set("Dependent Names", inFieldNames);
    p.set("Dependent Map", resToField);
    p.set("Preserve Diagonal", true);
    p.set<int>("Side Subcell Dimension", side_pb.getBaseCellTopology().getDimension() - 1);
    p.set<int>("Local Side ID", side_pb.cellData().side());
    p.set("Global Data Key", operatorName_+" Scatter Container");
    p.set("Check Apply BC",false);

    RCP< PHX::Evaluator<panzer::Traits> > op = elof->buildScatterDirichlet<EvalT>(p);

    fm.registerEvaluator<EvalT>(op);

    PHX::Tag<EvalT::ScalarT> tag(scatterName, Teuchos::rcp(new PHX::MDALayout<panzer::Dummy>(0)));
    fm.requireField<EvalT>(tag);
  }
}

#endif
