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

#ifndef PANZER_EVALUATOR_BASISTIMESSCALAR_IMPL_HPP
#define PANZER_EVALUATOR_BASISTIMESSCALAR_IMPL_HPP

#include "Intrepid_FunctionSpaceTools.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Workset_Utilities.hpp"

#define PANZER_USE_FAST_QUAD 1
// #define PANZER_USE_FAST_QUAD 0

namespace panzer {

//**********************************************************************
PHX_EVALUATOR_CTOR(Integrator_BasisTimesScalar,p) :
  residual( p.get<std::string>("Residual Name"), 
	    p.get< Teuchos::RCP<panzer::BasisIRLayout> >("Basis")->functional),
  scalar( p.get<std::string>("Value Name"), 
	  p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR")->dl_scalar),
  basis_name(p.get< Teuchos::RCP<panzer::BasisIRLayout> >("Basis")->name())
{
  using Teuchos::RCP;

  Teuchos::RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  Teuchos::RCP<const PureBasis> basis 
     = p.get< Teuchos::RCP<BasisIRLayout> >("Basis")->getBasis();

  // Verify that this basis supports the gradient operation
  TEUCHOS_TEST_FOR_EXCEPTION(!basis->isScalarBasis(),std::logic_error,
                             "Integrator_BasisTimesScalar: Basis of type \"" << basis->name() << "\" is not "
                             "a scalar basis");

  this->addEvaluatedField(residual);
  this->addDependentField(scalar);
    
  multiplier = p.get<double>("Multiplier");


  // build field multpliers if vector is nonnull (defaults to null)
  if(p.isType<Teuchos::RCP<const std::vector<std::string> > >("Field Multipliers")) {
    RCP<const std::vector<std::string> > field_multiplier_names = 
      p.get<Teuchos::RCP<const std::vector<std::string> > >("Field Multipliers");
    if(field_multiplier_names!=Teuchos::null) {
      for (std::vector<std::string>::const_iterator name = 
  	     field_multiplier_names->begin(); 
  	   name != field_multiplier_names->end(); ++name) {
        PHX::MDField<ScalarT,Cell,IP> tmp_field(*name, p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR")->dl_scalar);
        field_multipliers.push_back(tmp_field);
      }
    }
  }

  // add dependent field multiplers
  for (typename std::vector<PHX::MDField<ScalarT,Cell,IP> >::iterator field = field_multipliers.begin();
       field != field_multipliers.end(); ++field)
    this->addDependentField(*field);

  std::string n = "Integrator_BasisTimesScalar: " + residual.fieldTag().name();
  this->setName(n);
}

//**********************************************************************
PHX_POST_REGISTRATION_SETUP(Integrator_BasisTimesScalar,sd,fm)
{
  this->utils.setFieldData(residual,fm);
  this->utils.setFieldData(scalar,fm);
  
  for (typename std::vector<PHX::MDField<ScalarT,Cell,IP> >::iterator field = field_multipliers.begin();
       field != field_multipliers.end(); ++field)
    this->utils.setFieldData(*field,fm);

  num_nodes = residual.dimension(1);
  num_qp = scalar.dimension(1);

  basis_index = panzer::getBasisIndex(basis_name, (*sd.worksets_)[0]);

  tmp = Intrepid::FieldContainer<ScalarT>(scalar.dimension(0), num_qp); 
}

//**********************************************************************
PHX_EVALUATE_FIELDS(Integrator_BasisTimesScalar,workset)
{ 
  for (int i=0; i < residual.size(); ++i)
    residual[i] = 0.0;

#if PANZER_USE_FAST_QUAD
  // do a scaled copy
  for (int i=0; i < scalar.size(); ++i)
    tmp[i] = multiplier * scalar[i];

  for (typename std::vector<PHX::MDField<ScalarT,Cell,IP> >::iterator field = field_multipliers.begin();
	   field != field_multipliers.end(); ++field) {
    PHX::MDField<ScalarT,Cell,IP> field_data = *field;

    for (int i=0; i < field_data.size(); ++i)
      tmp[i] *= field_data[i];
  }

  const Intrepid::FieldContainer<double> & weighted_basis = workset.bases[basis_index]->weighted_basis;

  for (std::size_t cell = 0; cell < workset.num_cells; ++cell) {
    for (std::size_t basis = 0; basis < num_nodes; ++basis) {
      for (std::size_t qp = 0; qp < num_qp; ++qp) {
        residual(cell,basis) += tmp(cell,qp)*weighted_basis(cell,basis,qp);
      }
    }
  }

#else
  for (std::size_t cell = 0; cell < workset.num_cells; ++cell) {
    for (std::size_t qp = 0; qp < num_qp; ++qp) {
      tmp(cell,qp) = multiplier * scalar(cell,qp);
      for (typename std::vector<PHX::MDField<ScalarT,Cell,IP> >::iterator field = field_multipliers.begin();
	   field != field_multipliers.end(); ++field)
	tmp(cell,qp) = tmp(cell,qp) * (*field)(cell,qp);  
    }
  }

  if(workset.num_cells>0)
     Intrepid::FunctionSpaceTools::
       integrate<ScalarT>(residual, tmp, 
   		          (workset.bases[basis_index])->weighted_basis, 
		          Intrepid::COMP_BLAS);
#endif
}

//**********************************************************************
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList> 
Integrator_BasisTimesScalar<EvalT, Traits>::getValidParameters() const
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

