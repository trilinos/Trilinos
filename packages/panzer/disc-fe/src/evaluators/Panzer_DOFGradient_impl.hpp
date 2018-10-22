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

#ifndef PANZER_DOF_GRADIENT_IMPL_HPP
#define PANZER_DOF_GRADIENT_IMPL_HPP

#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"

namespace panzer {

namespace {

template <typename ScalarT,typename ArrayT>
struct evaluateGrad_withSens {
  evaluateGrad_withSens (PHX::MDField<ScalarT> dof_grad,
      PHX::MDField<const ScalarT,Cell,Point> dof_value,
      const ArrayT  grad_basis) :
        dof_grad_(dof_grad),dof_value_(dof_value),grad_basis_(grad_basis)
      {}
  KOKKOS_INLINE_FUNCTION
  void operator() (const size_t &cell) const
  {
    // evaluate at quadrature points
    int numFields = grad_basis_.extent(1);
    int numPoints = grad_basis_.extent(2);
    int spaceDim  = grad_basis_.extent(3);
    for (int pt=0; pt<numPoints; pt++) {
      for (int d=0; d<spaceDim; d++) {
        // first initialize to the right thing (prevents over writing with 0)
        // then loop over one less basis function
        dof_grad_(cell,pt,d) = dof_value_(cell, 0) * grad_basis_(cell, 0, pt, d);
        for (int bf=1; bf<numFields; bf++)
          dof_grad_(cell,pt,d) += dof_value_(cell, bf) * grad_basis_(cell, bf, pt, d);
      }
    }
  }
  PHX::MDField<ScalarT>  dof_grad_;
  PHX::MDField<const ScalarT,Cell,Point>  dof_value_;
  const ArrayT grad_basis_;

};

}

//**********************************************************************
template<typename EvalT, typename Traits>
DOFGradient<EvalT, Traits>::
DOFGradient(
  const Teuchos::ParameterList& p) :
  use_descriptors_(false),
  dof_value( p.get<std::string>("Name"), 
	     p.get< Teuchos::RCP<panzer::BasisIRLayout> >("Basis")->functional),
  dof_gradient( p.get<std::string>("Gradient Name"), 
		p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR")->dl_vector ),
  basis_name(p.get< Teuchos::RCP<panzer::BasisIRLayout> >("Basis")->name())
{
  Teuchos::RCP<const PureBasis> basis 
     = p.get< Teuchos::RCP<BasisIRLayout> >("Basis")->getBasis();

  // Verify that this basis supports the gradient operation
  TEUCHOS_TEST_FOR_EXCEPTION(!basis->supportsGrad(),std::logic_error,
                             "DOFGradient: Basis of type \"" << basis->name() << "\" does not support GRAD");

  this->addEvaluatedField(dof_gradient);
  this->addDependentField(dof_value);
  
  std::string n = "DOFGradient: " + dof_gradient.fieldTag().name() + " ("+PHX::typeAsString<EvalT>()+")";
  this->setName(n);
}

//**********************************************************************
template<typename EvalT, typename TRAITS>                   
DOFGradient<EvalT, TRAITS>::
DOFGradient(const PHX::FieldTag & input,
            const PHX::FieldTag & output,
            const panzer::BasisDescriptor & bd,
            const panzer::IntegrationDescriptor & id)
  : use_descriptors_(true)
  , bd_(bd) 
  , id_(id) 
  , dof_value(input)
  , dof_gradient(output)
{
  // Verify that this basis supports the gradient operation
  TEUCHOS_TEST_FOR_EXCEPTION(bd_.getType()=="HGrad",std::logic_error,
                             "DOFGradient: Basis of type \"" << bd_.getType() << "\" does not support GRAD");

  this->addEvaluatedField(dof_gradient);
  this->addDependentField(dof_value);
  
  std::string n = "DOFGradient: " + dof_gradient.fieldTag().name() + " ("+PHX::typeAsString<EvalT>()+")";
  this->setName(n);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
DOFGradient<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& fm)
{
  this->utils.setFieldData(dof_value,fm);
  this->utils.setFieldData(dof_gradient,fm);

  if(not use_descriptors_)
    basis_index = panzer::getBasisIndex(basis_name, (*sd.worksets_)[0], this->wda);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
DOFGradient<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{ 
  if (workset.num_cells == 0 )
    return;

  const panzer::BasisValues2<double> & basisValues = use_descriptors_ ?  this->wda(workset).getBasisValues(bd_,id_)
                                                                      : *this->wda(workset).bases[basis_index];

  typedef decltype(basisValues.grad_basis) ArrayT;

  evaluateGrad_withSens<ScalarT, ArrayT> eval(dof_gradient,dof_value,basisValues.grad_basis);
  Kokkos::parallel_for(workset.num_cells, eval);
}

//**********************************************************************

}

#endif
