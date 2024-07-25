// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_BCSTRATEGY_DIRICHLET_DEFAULT_IMPL_IMPL_HPP
#define PANZER_BCSTRATEGY_DIRICHLET_DEFAULT_IMPL_IMPL_HPP

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"

#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_FieldManager.hpp"

#include "Phalanx_MDField.hpp"
#include "Phalanx_DataLayout.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"

// Evaluators
#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_Dirichlet_Residual.hpp"
#include "Panzer_Dirichlet_Residual_EdgeBasis.hpp"
#include "Panzer_Dirichlet_Residual_FaceBasis.hpp"

#ifdef PANZER_HAVE_EPETRA_STACK
#include "Panzer_GatherSolution_Epetra.hpp"
#include "Panzer_ScatterDirichletResidual_Epetra.hpp"
#endif

#include "Panzer_GatherBasisCoordinates.hpp"
#include "Panzer_BasisValues_Evaluator.hpp"
#include "Panzer_PointValues_Evaluator.hpp"
#include "Panzer_DOF.hpp"
#include "Panzer_DOF_PointValues.hpp"

// ***********************************************************************
template <typename EvalT>
panzer::BCStrategy_Dirichlet_DefaultImpl<EvalT>::
BCStrategy_Dirichlet_DefaultImpl(const panzer::BC& bc,
				 const Teuchos::RCP<panzer::GlobalData>& global_data,
				 const bool in_check_apply_bc) :
  panzer::BCStrategy<EvalT>(bc),
  panzer::GlobalDataAcceptorDefaultImpl(global_data),
  check_apply_bc(in_check_apply_bc),
  descriptor_map_built(false)
{

}

// ***********************************************************************
template <typename EvalT>
panzer::BCStrategy_Dirichlet_DefaultImpl<EvalT>::
~BCStrategy_Dirichlet_DefaultImpl()
{

}

// ***********************************************************************
template <typename EvalT>
void panzer::BCStrategy_Dirichlet_DefaultImpl<EvalT>::
buildAndRegisterGatherScatterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
			                const panzer::PhysicsBlock& pb,
				        const panzer::LinearObjFactory<panzer::Traits> & lof,
					const Teuchos::ParameterList& user_data) const
{
  // for deprecated interface support
  buildDescriptorMapFromVectors();

  buildAndRegisterGatherAndOrientationEvaluators(fm,pb,lof,user_data);
  buildAndRegisterScatterEvaluators(fm,pb,lof,user_data);
}

// ***********************************************************************

template <typename EvalT>
void panzer::BCStrategy_Dirichlet_DefaultImpl<EvalT>::
buildAndRegisterScatterEvaluators(PHX::FieldManager<panzer::Traits>& fm,
                                  const panzer::PhysicsBlock& pb,
			          const LinearObjFactory<panzer::Traits> & lof,
			          const Teuchos::ParameterList& /* user_data */) const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::vector;
  using std::map;
  using std::string;
  using std::pair;

  // for deprecated interface support
  buildDescriptorMapFromVectors();

  // Scatter
  // for (map<string,string>::const_iterator res_to_dof = residual_to_dof_names_map.begin();
  //      res_to_dof != residual_to_dof_names_map.end(); ++res_to_dof) {
  for(DescriptorIterator itr=m_provided_dofs_desc.begin();
      itr!=m_provided_dofs_desc.end(); ++itr) {

    const DOFDescriptor & desc = itr->second;

    // there is no residual to scatter
    if(!desc.residualName.first)
      continue;

    // std::string residualName = res_to_dof->first;
    // std::string dofName = res_to_dof->second;
    std::string dofName      = desc.dofName;
    std::string residualName = desc.residualName.second;

    ParameterList p("Scatter: "+residualName + " to " + dofName);

    // Set name
    string scatter_field_name = "Dummy Scatter: " + this->m_bc.identifier() + residualName;
    p.set("Scatter Name", scatter_field_name);

    // Set basis
    const vector<pair<string,RCP<panzer::PureBasis> > >& dofBasisPair = pb.getProvidedDOFs();
    RCP<panzer::PureBasis> basis;
    for (vector<pair<string,RCP<panzer::PureBasis> > >::const_iterator it =
	   dofBasisPair.begin(); it != dofBasisPair.end(); ++it) {
      if (it->first == dofName)
	basis = it->second;
    }

    TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::is_null(basis), std::runtime_error,
		       "Error the name \"" << dofName
		       << "\" is not a valid DOF for the boundary condition:\n"
		       << this->m_bc << "\n");

    p.set("Basis", basis);

    RCP<vector<string> > residual_names = rcp(new vector<string>);
    residual_names->push_back(residualName);
    p.set("Dependent Names", residual_names);

    RCP<map<string,string> > names_map = rcp(new map<string,string>);
    names_map->insert(std::make_pair(residualName,dofName));
    p.set("Dependent Map", names_map);

    TEUCHOS_TEST_FOR_EXCEPTION(!pb.cellData().isSide(), std::logic_error,
		       "Error - physics block is not a side set!");

    p.set<int>("Side Subcell Dimension",
	       pb.getBaseCellTopology().getDimension() - 1);
    p.set<int>("Local Side ID", pb.cellData().side());

    p.set("Check Apply BC",check_apply_bc);

    RCP< PHX::Evaluator<panzer::Traits> > op = lof.buildScatterDirichlet<EvalT>(p);

    this->template registerEvaluator<EvalT>(fm, op);

    // Require variables
    {
      using panzer::Dummy;
      PHX::Tag<typename EvalT::ScalarT> tag(scatter_field_name,
					    rcp(new PHX::MDALayout<Dummy>(0)));
      fm.template requireField<EvalT>(tag);
    }

  }
}

// ***********************************************************************

template <typename EvalT>
void panzer::BCStrategy_Dirichlet_DefaultImpl<EvalT>::
buildAndRegisterGatherAndOrientationEvaluators(PHX::FieldManager<panzer::Traits>& fm,
	                                       const panzer::PhysicsBlock& pb,
					       const LinearObjFactory<panzer::Traits> & lof,
					       const Teuchos::ParameterList& /* user_data */) const
{
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using std::vector;
  using std::map;
  using std::string;
  using std::pair;

  // for deprecated interface support
  buildDescriptorMapFromVectors();

  // volume cell data object (used for handling vector valued fields)
  panzer::CellData cellData(pb.cellData().numCells(),pb.cellData().getCellTopology());

  // **************************
  // Coordinates for basis functions (no integration points needed)
  // **************************
  {
    const std::map<std::string,Teuchos::RCP<panzer::PureBasis> > & bases = pb.getBases();
    for (std::map<std::string,Teuchos::RCP<panzer::PureBasis> >::const_iterator it=bases.begin();
         it!=bases.end();it++) {

       Teuchos::RCP<panzer::PureBasis> basis = it->second;

       // add basis coordinates no matter what
       {
         RCP< PHX::Evaluator<panzer::Traits> > basis_op
            = rcp(new panzer::GatherBasisCoordinates<EvalT,panzer::Traits>(*it->second));
         this->template registerEvaluator<EvalT>(fm, basis_op);
       }

       // add point values and basis values
       if(basis->isVectorBasis()) {
         RCP<const panzer::PointRule> pointRule = rcp(new panzer::PointRule(basis->name()+":BasisPoints",basis->cardinality(),cellData));

         {
           RCP< PHX::Evaluator<panzer::Traits> > eval
             = rcp(new panzer::PointValues_Evaluator<EvalT,panzer::Traits>(pointRule,basis));

           this->template registerEvaluator<EvalT>(fm, eval);
         }
         {
           // note basis values are not constructed!
           RCP< PHX::Evaluator<panzer::Traits> > eval
             = rcp(new panzer::BasisValues_Evaluator<EvalT,panzer::Traits>(pointRule,basis,false));

           this->template registerEvaluator<EvalT>(fm, eval);
         }
       }
    }
  }

  // Gather
  for(DescriptorIterator itr=m_provided_dofs_desc.begin();
      itr!=m_provided_dofs_desc.end(); ++itr) {

    // get the dofName from the descriptor
    std::string dofName = itr->second.dofName;
    std::string fieldDof = !itr->second.timeDerivative.first
                              ? itr->second.dofName : itr->second.timeDerivative.second;

    const vector<pair<string,RCP<panzer::PureBasis> > >& dofBasisPair = pb.getProvidedDOFs();
    RCP<panzer::PureBasis> basis;
    for (vector<pair<string,RCP<panzer::PureBasis> > >::const_iterator it =
	   dofBasisPair.begin(); it != dofBasisPair.end(); ++it) {
      if (it->first == dofName)
	basis = it->second;
    }

    TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::is_null(basis), std::runtime_error,
		       "Error the name \"" << dofName
		       << "\" is not a valid DOF for the boundary condition:\n"
		       << this->m_bc << "\n");

    {
      ParameterList p("BC Gather");
      RCP<vector<string> > gather_field_names_vec = rcp(new vector<string>);
      RCP<vector<string> > gather_names_vec = rcp(new vector<string>);
      gather_field_names_vec->push_back(fieldDof);
      gather_names_vec->push_back(dofName);

      p.set("DOF Names", gather_field_names_vec);
      p.set("Indexer Names", gather_names_vec);
      p.set("Basis", basis);
      p.set("Use Time Derivative Solution Vector",itr->second.timeDerivative.first);

      RCP< PHX::Evaluator<panzer::Traits> > op = lof.buildGather<EvalT>(p);
      this->template registerEvaluator<EvalT>(fm, op);
    }

    if(basis->requiresOrientations())  {
      ParameterList p("Gather Orientation");
      RCP<vector<string> > gather_field_names_vec = rcp(new vector<string>);
      RCP<vector<string> > gather_names_vec = rcp(new vector<string>);
      gather_field_names_vec->push_back(fieldDof);
      gather_names_vec->push_back(dofName);

      p.set("DOF Names", gather_field_names_vec);
      p.set("Indexer Names", gather_names_vec);
      p.set("Basis", basis);

      RCP< PHX::Evaluator<panzer::Traits> > op = lof.buildGatherOrientation<EvalT>(p);

      this->template registerEvaluator<EvalT>(fm, op);
    }

    // evaluator a vector basis at the basis points
    if(basis->isVectorBasis()) {
      RCP<const panzer::PointRule> pointRule = rcp(new panzer::PointRule(basis->name()+":BasisPoints",basis->cardinality(),cellData));

      ParameterList p;
      p.set("Name",fieldDof);
      p.set("Basis",basis.getConst());
      p.set("Point Rule",pointRule);

      RCP< PHX::Evaluator<panzer::Traits> > eval
             = rcp(new panzer::DOF_PointValues<EvalT,panzer::Traits>(p));
      this->template registerEvaluator<EvalT>(fm, eval);
    }

  }

  // Dirichlet Residual: residual = dof_value - target_value
  for(DescriptorIterator itr=m_provided_dofs_desc.begin();
      itr!=m_provided_dofs_desc.end(); ++itr) {

    const DOFDescriptor & desc = itr->second;

    // there is no residual to scatter
    if(!desc.residualName.first)
      continue;

    // std::string dofName = desc.dofName;
    std::string dofName = !itr->second.timeDerivative.first
                              ? itr->second.dofName : itr->second.timeDerivative.second;
    std::string residualName = desc.residualName.second;
    std::string targetName = desc.targetName.second;

    const vector<pair<string,RCP<panzer::PureBasis> > >& dofBasisPair = pb.getProvidedDOFs();
    RCP<panzer::PureBasis> basis;
    for (vector<pair<string,RCP<panzer::PureBasis> > >::const_iterator it =
	   dofBasisPair.begin(); it != dofBasisPair.end(); ++it) {
      if (it->first == itr->second.dofName)
	basis = it->second;
    }

    TEUCHOS_TEST_FOR_EXCEPTION(Teuchos::is_null(basis), std::runtime_error,
		       "Error the name \"" << itr->second.dofName
		       << "\" is not a valid DOF for the boundary condition:\n"
		       << this->m_bc << "\n");

    if(basis->isScalarBasis() || desc.coefficientResidual) {
      ParameterList p("Dirichlet Residual: "+residualName + " to " + dofName);
      p.set("Residual Name", residualName);
      p.set("DOF Name", dofName);
      p.set("Value Name", targetName);
      p.set("Data Layout", basis->functional);

      RCP< PHX::Evaluator<panzer::Traits> > op =
        rcp(new panzer::DirichletResidual<EvalT,panzer::Traits>(p));

      this->template registerEvaluator<EvalT>(fm, op);
    }
    // This assumes that dofs on faces are all named "<dof>_face"
    else if(basis->isVectorBasis()&&basis->supportsDiv()) {
      RCP<const panzer::PointRule> pointRule = rcp(new panzer::PointRule(basis->name()+":BasisPoints",basis->cardinality(),cellData));

      ParameterList p;
      p.set("Residual Name", residualName);
      p.set("DOF Name", dofName);
      p.set("Value Name", targetName);
      p.set("Basis", basis.getConst());
      p.set("Point Rule", pointRule);

      RCP< PHX::Evaluator<panzer::Traits> > op =
        rcp(new panzer::DirichletResidual_FaceBasis<EvalT,panzer::Traits>(p));

      this->template registerEvaluator<EvalT>(fm, op);
    }
    else if(basis->isVectorBasis()) {
      RCP<const panzer::PointRule> pointRule = rcp(new panzer::PointRule(basis->name()+":BasisPoints",basis->cardinality(),cellData));

      ParameterList p;
      p.set("Residual Name", residualName);
      p.set("DOF Name", dofName);
      p.set("Value Name", targetName);
      p.set("Basis", basis.getConst());
      p.set("Point Rule", pointRule);

      RCP< PHX::Evaluator<panzer::Traits> > op =
        rcp(new panzer::DirichletResidual_EdgeBasis<EvalT,panzer::Traits>(p));

      this->template registerEvaluator<EvalT>(fm, op);
    }

  }
}

// ***********************************************************************

template <typename EvalT>
void panzer::BCStrategy_Dirichlet_DefaultImpl<EvalT>::
buildDescriptorMapFromVectors() const
{
  if(descriptor_map_built)
    return;

  // build the reverse lookup (this assumes that the map is invertible, though it should be!)
  std::map<std::string,std::string> dof_names_to_residual_map;
  for(std::map<std::string,std::string>::const_iterator itr=residual_to_dof_names_map.begin();
      itr!=residual_to_dof_names_map.end();++itr) {
    dof_names_to_residual_map[itr->second] = itr->first;
  }

  for(std::size_t i=0;i<required_dof_names.size();i++) {
    std::string dof_name      = required_dof_names[i];

    // add the DOF right away
    const_cast<BCStrategy_Dirichlet_DefaultImpl<EvalT> *>(this)->addDOF(dof_name);

    // add target information if its required
    if(dof_names_to_residual_map.find(dof_name)!=dof_names_to_residual_map.end()) {
      std::string residual_name = dof_names_to_residual_map[dof_name];
      std::string target_name   = residual_to_target_field_map.find(residual_name)->second;

      // add the descriptor from the data structures: Note the const_cast is neccessary
      // only because this is protecting deprecated code
      const_cast<BCStrategy_Dirichlet_DefaultImpl<EvalT> *>(this)->
      addTarget(target_name,
                dof_name,
                residual_name);
    }
  }

  descriptor_map_built = true;
}

// ***********************************************************************

template <typename EvalT>
void panzer::BCStrategy_Dirichlet_DefaultImpl<EvalT>::
addDOF(const std::string & dofName)
{
  DescriptorIterator itr = m_provided_dofs_desc.find(dofName);
  TEUCHOS_ASSERT(itr==m_provided_dofs_desc.end());

  DOFDescriptor & desc = m_provided_dofs_desc[dofName];
  desc.dofName = dofName;
}

// ***********************************************************************

template <typename EvalT>
void panzer::BCStrategy_Dirichlet_DefaultImpl<EvalT>::
addTarget(const std::string & targetName,
          const std::string & dofName,
          const std::string & residualName)
{
  typename std::map<std::string,DOFDescriptor>::iterator itr = m_provided_dofs_desc.find(dofName);
  TEUCHOS_ASSERT(itr!=m_provided_dofs_desc.end());

  DOFDescriptor & desc = itr->second;
  desc.dofName = dofName;
  desc.coefficientResidual = false;
  desc.targetName = std::make_pair(true,targetName);
  desc.residualName = (residualName=="") ? std::make_pair(true,"RESIDUAL_"+dofName)
                                         : std::make_pair(true,residualName);
}

template <typename EvalT>
void panzer::BCStrategy_Dirichlet_DefaultImpl<EvalT>::
addCoefficientTarget(const std::string & targetName,
                     const std::string & dofName,
                     const std::string & residualName)
{
  typename std::map<std::string,DOFDescriptor>::iterator itr = m_provided_dofs_desc.find(dofName);
  TEUCHOS_ASSERT(itr!=m_provided_dofs_desc.end());

  DOFDescriptor & desc = itr->second;
  desc.dofName = dofName;
  desc.coefficientResidual = true;
  desc.targetName = std::make_pair(true,targetName);
  desc.residualName = (residualName=="") ? std::make_pair(true,"RESIDUAL_"+dofName)
                                         : std::make_pair(true,residualName);
}

// ***********************************************************************

template <typename EvalT>
void panzer::BCStrategy_Dirichlet_DefaultImpl<EvalT>::
addDotTarget(const std::string & targetName,
             const std::string & dofName,
             const std::string & dotName,
             const std::string & residualName)
{
  typename std::map<std::string,DOFDescriptor>::iterator itr = m_provided_dofs_desc.find(dofName);
  TEUCHOS_ASSERT(itr!=m_provided_dofs_desc.end());

  DOFDescriptor & desc = itr->second;
  desc.dofName = dofName;
  desc.targetName = std::make_pair(true,targetName);
  desc.timeDerivative = (dotName=="") ? std::make_pair(true,"DXDT_"+dofName)
                                      : std::make_pair(true,dotName);
  desc.residualName = (residualName=="") ? std::make_pair(true,"RESIDUAL_"+dofName)
                                         : std::make_pair(true,residualName);
}

// ***********************************************************************

#endif
