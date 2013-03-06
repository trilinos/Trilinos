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

#ifndef PANZER_EVALUATOR_GRADBASISDOTVECTOR_IMPL_HPP
#define PANZER_EVALUATOR_GRADBASISDOTVECTOR_IMPL_HPP

#include "Intrepid_FunctionSpaceTools.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Workset_Utilities.hpp"

#define PANZER_USE_FAST_QUAD 1
// #define PANZER_USE_FAST_QUAD 0

namespace panzer {

//**********************************************************************
PHX_EVALUATOR_CTOR(Integrator_GradBasisDotVector,p) :
  residual( p.get<std::string>("Residual Name"), 
	    p.get< Teuchos::RCP<panzer::BasisIRLayout> >("Basis")->functional),
  flux( p.get<std::string>("Flux Name"), 
	p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR")->dl_vector ),
  basis_name(p.get< Teuchos::RCP<panzer::BasisIRLayout> >("Basis")->name())
{
  Teuchos::RCP<const PureBasis> basis 
     = p.get< Teuchos::RCP<BasisIRLayout> >("Basis")->getBasis();

  // Verify that this basis supports the gradient operation
  TEUCHOS_TEST_FOR_EXCEPTION(!basis->supportsGrad(),std::logic_error,
                             "Integrator_GradBasisDotVector: Basis of type \"" << basis->name() << "\" does not support GRAD");

  this->addEvaluatedField(residual);
  this->addDependentField(flux);
  
  multiplier = p.get<double>("Multiplier");

  if (p.isType<Teuchos::RCP<const std::vector<std::string> > >("Field Multipliers")) 
  {
    const std::vector<std::string>& field_multiplier_names = 
      *(p.get<Teuchos::RCP<const std::vector<std::string> > >("Field Multipliers"));

    for (std::vector<std::string>::const_iterator name = field_multiplier_names.begin(); 
      name != field_multiplier_names.end(); ++name) 
    {
      PHX::MDField<ScalarT,Cell,IP> tmp_field(*name, p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR")->dl_scalar);
      field_multipliers.push_back(tmp_field);
    }
  }

  for (typename std::vector<PHX::MDField<ScalarT,Cell,IP> >::iterator field = field_multipliers.begin();
       field != field_multipliers.end(); ++field)
    this->addDependentField(*field);

  std::string n = 
    "Integrator_GradBasisDotVector: " + residual.fieldTag().name();

  this->setName(n);
}

//**********************************************************************
PHX_POST_REGISTRATION_SETUP(Integrator_GradBasisDotVector,sd,fm)
{
  this->utils.setFieldData(residual,fm);
  this->utils.setFieldData(flux,fm);

  for (typename std::vector<PHX::MDField<ScalarT,Cell,IP> >::iterator field = field_multipliers.begin();
       field != field_multipliers.end(); ++field)
    this->utils.setFieldData(*field,fm);

  num_nodes = residual.dimension(1);
  num_qp = flux.dimension(1);
  num_dim = flux.dimension(2);

  basis_index = panzer::getBasisIndex(basis_name, (*sd.worksets_)[0]);

  tmp = Intrepid::FieldContainer<ScalarT>(flux.dimension(0), num_qp, num_dim); 
}

//**********************************************************************
PHX_EVALUATE_FIELDS(Integrator_GradBasisDotVector,workset)
{ 
  for (int i=0; i < residual.size(); ++i)
    residual[i] = 0.0;

#if PANZER_USE_FAST_QUAD
  // do a scaled copy
  for (int i=0; i < flux.size(); ++i)
    tmp[i] = multiplier * flux[i];

  for (typename std::vector<PHX::MDField<ScalarT,Cell,IP> >::iterator field = field_multipliers.begin();
       field != field_multipliers.end(); ++field) {
    PHX::MDField<ScalarT,Cell,IP> field_data = *field;

    for (std::size_t cell = 0; cell < workset.num_cells; ++cell) {
      for (std::size_t qp = 0; qp < num_qp; ++qp) {
        ScalarT tmpVar = field_data(cell,qp);  

        for (std::size_t dim = 0; dim < num_dim; ++dim)
          tmp(cell,qp,dim) *= tmpVar;
      }
    }
  } 

  const Intrepid::FieldContainer<double> & weighted_grad_basis = workset.bases[basis_index]->weighted_grad_basis;

  // perform integration and vector dot product (at the same time! whoah!)
  for (std::size_t cell = 0; cell < workset.num_cells; ++cell) {
    for (std::size_t basis = 0; basis < num_nodes; ++basis) {
      for (std::size_t qp = 0; qp < num_qp; ++qp) {
        for (std::size_t dim = 0; dim < num_dim; ++dim)
          residual(cell,basis) += tmp(cell,qp,dim)*weighted_grad_basis(cell,basis,qp,dim);
      }
    }
  }
#else
  for (std::size_t cell = 0; cell < workset.num_cells; ++cell)
  {
    for (std::size_t qp = 0; qp < num_qp; ++qp)
    {
      ScalarT tmpVar = 1.0;
      for (typename std::vector<PHX::MDField<ScalarT,Cell,IP> >::iterator field = field_multipliers.begin();
           field != field_multipliers.end(); ++field)
        tmpVar = tmpVar * (*field)(cell,qp);  

      for (std::size_t dim = 0; dim < num_dim; ++dim)
        tmp(cell,qp,dim) = multiplier * tmpVar * flux(cell,qp,dim);
    }
  }
  
  if(workset.num_cells>0)
     Intrepid::FunctionSpaceTools::
       integrate<ScalarT>(residual, tmp, 
   		       (workset.bases[basis_index])->weighted_grad_basis, 
		       Intrepid::COMP_BLAS);
#endif
}

//**********************************************************************

}

#endif
