// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef __Panzer_Integrator_DivBasisTimesScalar_impl_hpp__ 
#define __Panzer_Integrator_DivBasisTimesScalar_impl_hpp__ 

#include "Intrepid2_FunctionSpaceTools.hpp"

#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Kokkos_ViewFactory.hpp"

namespace panzer {

//**********************************************************************
template<typename EvalT, typename Traits>
Integrator_DivBasisTimesScalar<EvalT, Traits>::
Integrator_DivBasisTimesScalar(
  const Teuchos::ParameterList& p) :
  residual( p.get<std::string>("Residual Name"), 
	    p.get< Teuchos::RCP<panzer::BasisIRLayout> >("Basis")->functional),
  basis_name(p.get< Teuchos::RCP<panzer::BasisIRLayout> >("Basis")->name())
{
  Teuchos::RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  Teuchos::RCP<const PureBasis> basis 
     = p.get< Teuchos::RCP<BasisIRLayout> >("Basis")->getBasis();

  // Verify that this basis supports the curl operation
  TEUCHOS_TEST_FOR_EXCEPTION(!basis->supportsDiv(),std::logic_error,
                             "Integrator_DivBasisTimesScalar: Basis of type \"" << basis->name() << "\" does not support DIV.");
  TEUCHOS_TEST_FOR_EXCEPTION(!basis->requiresOrientations(),std::logic_error,
                             "Integration_DivBasisTimesScalar: Basis of type \"" << basis->name() << "\" should require orientations. So we are throwing.");

  scalar = PHX::MDField<const ScalarT,Cell,IP>( p.get<std::string>("Value Name"),
                                          p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR")->dl_scalar );

  this->addEvaluatedField(residual);
  this->addDependentField(scalar);
  
  multiplier = p.get<double>("Multiplier");
  if (p.isType<Teuchos::RCP<const std::vector<std::string> > >("Field Multipliers")) {

    const std::vector<std::string>& field_multiplier_names = 
      *(p.get<Teuchos::RCP<const std::vector<std::string> > >("Field Multipliers"));

    for (const std::string & name : field_multiplier_names) {
      PHX::MDField<const ScalarT,Cell,IP> tmp_field(name, p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR")->dl_scalar);
      field_multipliers.push_back(tmp_field);
    }
  }

  for (const auto & field : field_multipliers)
    this->addDependentField(field);

  std::string n = 
    "Integrator_DivBasisTimesScalar: " + residual.fieldTag().name();

  this->setName(n);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
Integrator_DivBasisTimesScalar<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& fm)
{
  this->utils.setFieldData(residual,fm);
  this->utils.setFieldData(scalar,fm);

  for (auto & field : field_multipliers)
    this->utils.setFieldData(field,fm);

  num_nodes = residual.extent(1);
  num_qp = scalar.extent(1);

  basis_index = panzer::getBasisIndex(basis_name, (*sd.worksets_)[0], this->wda);

  tmp = Kokkos::createDynRankView(residual.get_static_view(),"tmp",scalar.extent(0), num_qp); 
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
Integrator_DivBasisTimesScalar<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{ 
  // zero the reisdual
  residual.deep_copy(ScalarT(0.0));
  
  for (index_t cell = 0; cell < workset.num_cells; ++cell) {
    for (std::size_t qp = 0; qp < num_qp; ++qp) {
      ScalarT tmpVar = 1.0;
      for (const auto & field : field_multipliers)
        tmpVar = tmpVar * field(cell,qp);

      // no dimension to loop over for scalar fields
      tmp(cell,qp) = multiplier * tmpVar * scalar(cell,qp);
    }
  }
  
  {
    // const Kokkos::DynRankView<double,PHX::Device> & weighted_div_basis = (this->wda(workset).bases[basis_index])->weighted_div_basis;
    const BasisValues2<double> & bv = *this->wda(workset).bases[basis_index];

    for (index_t cell = 0; cell < workset.num_cells; ++cell)
      for (std::size_t basis = 0; basis < num_nodes; ++basis) {
        for (std::size_t qp = 0; qp < num_qp; ++qp)
          residual(cell,basis) += tmp(cell,qp)*bv.weighted_div_basis(cell,basis,qp);
      }
  }
/*
  if(workset.num_cells>0) {
     Intrepid2::FunctionSpaceTools::
       integrate<ScalarT>(residual, tmp, 
                       this->wda(workset).bases[basis_index]->weighted_div_basis, 
		       Intrepid2::COMP_CPP);
  }
*/
}

//**********************************************************************

template<typename EvalT, typename TRAITS>
Teuchos::RCP<Teuchos::ParameterList> 
Integrator_DivBasisTimesScalar<EvalT, TRAITS>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set<std::string>("Residual Name", "?");
  p->set<std::string>("Value Name", "?");
  p->set<std::string>("Test Field Name", "?");
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
