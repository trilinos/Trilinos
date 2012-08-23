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

#ifndef PANZER_EVALUATOR_SCALAR_IMPL_HPP
#define PANZER_EVALUATOR_SCALAR_IMPL_HPP

#include "Intrepid_FunctionSpaceTools.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"

namespace panzer {

//**********************************************************************
PHX_EVALUATOR_CTOR(Integrator_Scalar,p) : quad_index(-1)
{
  Teuchos::RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  Teuchos::RCP<panzer::IntegrationRule> ir = p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR");
  quad_order = ir->cubature_degree;

  Teuchos::RCP<PHX::DataLayout> dl_cell = Teuchos::rcp(new PHX::MDALayout<Cell>(ir->dl_scalar->dimension(0)));
  integral = PHX::MDField<ScalarT,Cell>( p.get<std::string>("Integral Name"), dl_cell);
  scalar = PHX::MDField<ScalarT,Cell,IP>( p.get<std::string>("Integrand Name"), ir->dl_scalar);

  this->addEvaluatedField(integral);
  this->addDependentField(scalar);
    
  multiplier = 1.0;
  if(p.isType<double>("Multiplier"))
     multiplier = p.get<double>("Multiplier");

  if (p.isType<Teuchos::RCP<const std::vector<std::string> > >("Field Multipliers")) {
    const std::vector<std::string>& field_multiplier_names = 
      *(p.get<Teuchos::RCP<const std::vector<std::string> > >("Field Multipliers"));

    for (std::vector<std::string>::const_iterator name = 
	   field_multiplier_names.begin(); 
	 name != field_multiplier_names.end(); ++name) {
      PHX::MDField<ScalarT,Cell,IP> tmp_field(*name, p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR")->dl_scalar);
      field_multipliers.push_back(tmp_field);
    }
  }

  for (typename std::vector<PHX::MDField<ScalarT,Cell,IP> >::iterator field = field_multipliers.begin();
       field != field_multipliers.end(); ++field)
    this->addDependentField(*field);

  std::string n = "Integrator_Scalar: " + integral.fieldTag().name();
  this->setName(n);
}

//**********************************************************************
PHX_POST_REGISTRATION_SETUP(Integrator_Scalar,sd,fm)
{
  this->utils.setFieldData(integral,fm);
  this->utils.setFieldData(scalar,fm);
  
  for (typename std::vector<PHX::MDField<ScalarT,Cell,IP> >::iterator field = field_multipliers.begin();
       field != field_multipliers.end(); ++field)
    this->utils.setFieldData(*field,fm);

  num_qp = scalar.dimension(1);

  tmp = Intrepid::FieldContainer<ScalarT>(scalar.dimension(0), num_qp); 

  quad_index =  panzer::getIntegrationRuleIndex(quad_order,(*sd.worksets_)[0]);
}

//**********************************************************************
PHX_EVALUATE_FIELDS(Integrator_Scalar,workset)
{ 
  for (int i=0; i < integral.size(); ++i)
    integral[i] = 0.0;

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
       integrate<ScalarT>(integral, tmp, 
   		          (workset.int_rules[quad_index])->weighted_measure, 
		          Intrepid::COMP_BLAS);
}

//**********************************************************************
template<typename EvalT, typename Traits>
Teuchos::RCP<Teuchos::ParameterList> 
Integrator_Scalar<EvalT, Traits>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set<std::string>("Integral Name", "?");
  p->set<std::string>("Integrand Name", "?");

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
