// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_EVALUATOR_TRANSIENT_BASISTIMESSCALAR_IMPL_HPP
#define PANZER_EVALUATOR_TRANSIENT_BASISTIMESSCALAR_IMPL_HPP

#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Kokkos_ViewFactory.hpp"

namespace panzer {

//**********************************************************************
template<typename EvalT, typename Traits>
Integrator_TransientBasisTimesScalar<EvalT, Traits>::
Integrator_TransientBasisTimesScalar(
  const Teuchos::ParameterList& p) :
  residual( p.get<std::string>("Residual Name"), 
	    p.get< Teuchos::RCP<panzer::BasisIRLayout> >("Basis")->functional),
  scalar( p.get<std::string>("Value Name"), 
	  p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR")->dl_scalar),
  basis_name(p.get< Teuchos::RCP<panzer::BasisIRLayout> >("Basis")->name())
{
  Teuchos::RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  Teuchos::RCP<const PureBasis> basis 
     = p.get< Teuchos::RCP<BasisIRLayout> >("Basis")->getBasis();

  // Verify that this basis supports the gradient operation
  TEUCHOS_TEST_FOR_EXCEPTION(!basis->isScalarBasis(),std::logic_error,
                             "Integrator_TransientBasisTimesScalar: Basis of type \"" << basis->name() << "\" is not a "
                             "scalar basis");

  this->addEvaluatedField(residual);
  this->addDependentField(scalar);
    
  multiplier = p.get<double>("Multiplier");


  if (p.isType<Teuchos::RCP<const std::vector<std::string> > >("Field Multipliers")) {
    const std::vector<std::string>& field_multiplier_names = 
      *(p.get<Teuchos::RCP<const std::vector<std::string> > >("Field Multipliers"));

    for (std::vector<std::string>::const_iterator name = 
	   field_multiplier_names.begin(); 
	 name != field_multiplier_names.end(); ++name) {
      PHX::MDField<const ScalarT,Cell,IP> tmp_field(*name, p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR")->dl_scalar);
      field_multipliers.push_back(tmp_field);
    }
  }

  for (typename std::vector<PHX::MDField<const ScalarT,Cell,IP> >::iterator field = field_multipliers.begin();
       field != field_multipliers.end(); ++field)
    this->addDependentField(*field);

  std::string n = "Integrator_TransientBasisTimesScalar: " + residual.fieldTag().name();
  this->setName(n);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
Integrator_TransientBasisTimesScalar<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& fm)
{
  this->utils.setFieldData(residual,fm);
  this->utils.setFieldData(scalar,fm);
  
  num_nodes = residual.extent(1);
  num_qp = scalar.extent(1);

  basis_index = panzer::getBasisIndex(basis_name, (*sd.worksets_)[0], this->wda);

  tmp = Kokkos::createDynRankView(residual.get_static_view(),"tmp",scalar.extent(0), num_qp); 
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
Integrator_TransientBasisTimesScalar<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{ 
  if (workset.evaluate_transient_terms) {
    
   // for (int i=0; i < residual.size(); ++i)
   //   residual[i] = 0.0;
    
   Kokkos::deep_copy (residual.get_static_view(), ScalarT(0.0));

    for (index_t cell = 0; cell < workset.num_cells; ++cell) {
      for (std::size_t qp = 0; qp < num_qp; ++qp) {
	tmp(cell,qp) = multiplier * scalar(cell,qp);
	for (typename std::vector<PHX::MDField<const ScalarT,Cell,IP> >::iterator field = field_multipliers.begin();
	     field != field_multipliers.end(); ++field)
	  tmp(cell,qp) = tmp(cell,qp) * (*field)(cell,qp);  
      }
    }

    if(workset.num_cells>0)
      Intrepid2::FunctionSpaceTools<PHX::exec_space>::
        integrate<ScalarT>(residual.get_view(),
                           tmp, 
			   (this->wda(workset).bases[basis_index])->weighted_basis_scalar.get_view());
  }
}

//**********************************************************************
template<typename EvalT, typename TRAITS>
Teuchos::RCP<Teuchos::ParameterList> 
Integrator_TransientBasisTimesScalar<EvalT, TRAITS>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set<std::string>("Residual Name", "?");
  p->set<std::string>("Value Name", "?");
  Teuchos::RCP<panzer::BasisIRLayout> basis;
  p->set("Basis", basis);
  Teuchos::RCP<panzer::IntegrationRule> ir;
  p->set("IR", ir);
  p->set<double>("Multiplier", 1.0);
  Teuchos::RCP<const std::vector<std::string> > fms;
  p->set("Field Multipliers", fms);
  return p;
}

//**********************************************************************

}

#endif

