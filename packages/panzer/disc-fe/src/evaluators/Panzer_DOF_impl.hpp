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

#ifndef PANZER_DOF_IMPL_HPP
#define PANZER_DOF_IMPL_HPP

#include <algorithm>
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_CommonArrayFactories.hpp"
#include "Panzer_DOF_Functors.hpp"

#include "Intrepid2_FunctionSpaceTools.hpp"

namespace panzer {

//**********************************************************************
//* DOF evaluator
//**********************************************************************

//**********************************************************************
// MOST EVALUATION TYPES
//**********************************************************************

//**********************************************************************
template<typename EvalT, typename TRAITS>                   
DOF<EvalT, TRAITS>::
DOF(const Teuchos::ParameterList & p) :
  use_descriptors_(false),
  dof_basis( p.get<std::string>("Name"), 
	     p.get< Teuchos::RCP<BasisIRLayout> >("Basis")->functional),
  basis_name(p.get< Teuchos::RCP<BasisIRLayout> >("Basis")->name())
{
  Teuchos::RCP<const PureBasis> basis 
     = p.get< Teuchos::RCP<BasisIRLayout> >("Basis")->getBasis();
  is_vector_basis = basis->isVectorBasis();

  // swap between scalar basis value, or vector basis value
  if(basis->isScalarBasis()) {
     dof_ip_scalar = PHX::MDField<ScalarT,Cell,Point>(
                p.get<std::string>("Name"), 
     	        p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR")->dl_scalar);
     this->addEvaluatedField(dof_ip_scalar);
  }
  else if(basis->isVectorBasis()) {
     dof_ip_vector = PHX::MDField<ScalarT,Cell,Point,Dim>(
                p.get<std::string>("Name"), 
     	        p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR")->dl_vector);
     this->addEvaluatedField(dof_ip_vector);
  }
  else
  { TEUCHOS_ASSERT(false); }

  this->addDependentField(dof_basis);

  std::string n = "DOF: " + dof_basis.fieldTag().name() + " ("+PHX::typeAsString<EvalT>()+")";
  this->setName(n);
}

//**********************************************************************
template<typename EvalT, typename TRAITS>                   
DOF<EvalT, TRAITS>::
DOF(const PHX::FieldTag & input,
    const PHX::FieldTag & output,
    const panzer::BasisDescriptor & bd,
    const panzer::IntegrationDescriptor & id)
  : use_descriptors_(true)
  , bd_(bd) 
  , id_(id) 
  , dof_basis(input)
{
  TEUCHOS_ASSERT(bd.getType()=="HGrad" || bd.getType()=="HCurl" || 
                 bd.getType()=="HDiv" || bd.getType()=="Const")

  is_vector_basis = (bd.getType()=="HCurl" || bd.getType()=="HDiv");
 
  // swap between scalar basis value, or vector basis value
  if(not is_vector_basis) {
     dof_ip_scalar = output;
     this->addEvaluatedField(dof_ip_scalar);
  }
  else {
     dof_ip_vector = output;
     this->addEvaluatedField(dof_ip_vector);
  }

  this->addDependentField(dof_basis);

  std::string n = "DOF: " + dof_basis.fieldTag().name() + " ("+PHX::typeAsString<EvalT>()+")";
  this->setName(n);
}

//**********************************************************************
template<typename EvalT, typename TRAITS>                   
void DOF<EvalT, TRAITS>::
postRegistrationSetup(typename TRAITS::SetupData sd,
                      PHX::FieldManager<TRAITS>& fm)
{
  this->utils.setFieldData(dof_basis,fm);
  if(is_vector_basis)
    this->utils.setFieldData(dof_ip_vector,fm);
  else
    this->utils.setFieldData(dof_ip_scalar,fm);

  // descriptors don't access the basis values in the same way
  if(not use_descriptors_)
    basis_index = panzer::getBasisIndex(basis_name, (*sd.worksets_)[0], this->wda);
}

//**********************************************************************
template<typename EvalT, typename TRAITS>                   
void DOF<EvalT, TRAITS>::
evaluateFields(typename TRAITS::EvalData workset)
{ 
  const panzer::BasisValues2<double> & basisValues = use_descriptors_ ?  this->wda(workset).getBasisValues(bd_,id_)
                                                                      : *this->wda(workset).bases[basis_index];

  if(is_vector_basis) {
    int spaceDim  = basisValues.basis_vector.extent(3);
    if(spaceDim==3) {
      dof_functors::EvaluateDOFWithSens_Vector<ScalarT,typename BasisValues2<double>::Array_CellBasisIPDim,3> functor(dof_basis,dof_ip_vector,basisValues.basis_vector);
      Kokkos::parallel_for(workset.num_cells,functor);
    }
    else {
      dof_functors::EvaluateDOFWithSens_Vector<ScalarT,typename BasisValues2<double>::Array_CellBasisIPDim,2> functor(dof_basis,dof_ip_vector,basisValues.basis_vector);
      Kokkos::parallel_for(workset.num_cells,functor);
    }
  }
  else {
    dof_functors::EvaluateDOFWithSens_Scalar<ScalarT,typename BasisValues2<double>::Array_CellBasisIP> functor(dof_basis,dof_ip_scalar,basisValues.basis_scalar);
    Kokkos::parallel_for(workset.num_cells,functor);
  }
}

//**********************************************************************

//**********************************************************************
// JACOBIAN EVALUATION TYPES
//**********************************************************************

//**********************************************************************
template<typename TRAITS>                   
DOF<typename TRAITS::Jacobian, TRAITS>::
DOF(const Teuchos::ParameterList & p) :
  use_descriptors_(false),
  dof_basis( p.get<std::string>("Name"), 
	     p.get< Teuchos::RCP<BasisIRLayout> >("Basis")->functional),
  basis_name(p.get< Teuchos::RCP<BasisIRLayout> >("Basis")->name())
{
  Teuchos::RCP<const PureBasis> basis 
     = p.get< Teuchos::RCP<BasisIRLayout> >("Basis")->getBasis();
  is_vector_basis = basis->isVectorBasis();

  if(p.isType<Teuchos::RCP<const std::vector<int> > >("Jacobian Offsets Vector")) {
    const std::vector<int> & offsets = *p.get<Teuchos::RCP<const std::vector<int> > >("Jacobian Offsets Vector");

    // allocate and copy offsets vector to Kokkos array
    offsets_array = Kokkos::View<int*,PHX::Device>("offsets",offsets.size());
    for(std::size_t i=0;i<offsets.size();i++)
      offsets_array(i) = offsets[i];

    accelerate_jacobian_enabled = true;  // short cut for identity matrix

    // get the sensitivities name that is valid for accelerated jacobians
    sensitivities_name = true;
    if (p.isType<std::string>("Sensitivities Name"))
      sensitivities_name = p.get<std::string>("Sensitivities Name");
  }
  else
    accelerate_jacobian_enabled = false; // don't short cut for identity matrix

  // swap between scalar basis value, or vector basis value
  if(basis->isScalarBasis()) {
     dof_ip_scalar = PHX::MDField<ScalarT,Cell,Point>(
                p.get<std::string>("Name"), 
     	        p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR")->dl_scalar);
     this->addEvaluatedField(dof_ip_scalar);
  }
  else if(basis->isVectorBasis()) {
     dof_ip_vector = PHX::MDField<ScalarT,Cell,Point,Dim>(
                p.get<std::string>("Name"), 
     	        p.get< Teuchos::RCP<panzer::IntegrationRule> >("IR")->dl_vector);
     this->addEvaluatedField(dof_ip_vector);
  }
  else
  { TEUCHOS_ASSERT(false); }

  this->addDependentField(dof_basis);

  std::string n = "DOF: " + dof_basis.fieldTag().name() 
                          + ( accelerate_jacobian_enabled ? " accel_jac " : "slow_jac" ) 
                          + " ("+PHX::typeAsString<panzer::Traits::Jacobian>()+")";
  this->setName(n);
}

//**********************************************************************
template<typename TRAITS>                   
DOF<typename TRAITS::Jacobian, TRAITS>::
DOF(const PHX::FieldTag & input,
    const PHX::FieldTag & output,
    const panzer::BasisDescriptor & bd,
    const panzer::IntegrationDescriptor & id)
  : use_descriptors_(true)
  , bd_(bd) 
  , id_(id) 
  , dof_basis(input)
{
  TEUCHOS_ASSERT(bd.getType()=="HGrad" || bd.getType()=="HCurl" || 
                 bd.getType()=="HDiv" || bd.getType()=="Const")

  accelerate_jacobian_enabled = false; // don't short cut for identity matrix

  is_vector_basis = (bd.getType()=="HCurl" || bd.getType()=="HDiv");
 
  // swap between scalar basis value, or vector basis value
  if(not is_vector_basis) {
     dof_ip_scalar = output;
     this->addEvaluatedField(dof_ip_scalar);
  }
  else {
     dof_ip_vector = output;
     this->addEvaluatedField(dof_ip_vector);
  }

  this->addDependentField(dof_basis);

  std::string n = "DOF: " + dof_basis.fieldTag().name() + " slow_jac(descriptor) ("+PHX::typeAsString<typename TRAITS::Jacobian>()+")";
  this->setName(n);
}

//**********************************************************************
template<typename TRAITS>                   
void DOF<typename TRAITS::Jacobian, TRAITS>::
postRegistrationSetup(typename TRAITS::SetupData sd,
                      PHX::FieldManager<TRAITS>& fm)
{
  this->utils.setFieldData(dof_basis,fm);
  if(is_vector_basis)
    this->utils.setFieldData(dof_ip_vector,fm);
  else
    this->utils.setFieldData(dof_ip_scalar,fm);

  // descriptors don't access the basis values in the same way
  if(not use_descriptors_)
    basis_index = panzer::getBasisIndex(basis_name, (*sd.worksets_)[0], this->wda);
}

// **********************************************************************
template<typename TRAITS>
void DOF<typename TRAITS::Jacobian, TRAITS>::
preEvaluate(typename TRAITS::PreEvalData d)
{
  // if sensitivities were requrested for this field enable accelerated
  // jacobian calculations
  accelerate_jacobian = false;
  if(accelerate_jacobian_enabled && d.first_sensitivities_name==sensitivities_name) {
    accelerate_jacobian = true;
  }
}

//**********************************************************************
template<typename TRAITS>                   
void DOF<typename TRAITS::Jacobian, TRAITS>::
evaluateFields(typename TRAITS::EvalData workset)
{ 
  const panzer::BasisValues2<double> & basisValues = use_descriptors_ ?  this->wda(workset).getBasisValues(bd_,id_)
                                                                      : *this->wda(workset).bases[basis_index];

  if(is_vector_basis) {
    if(accelerate_jacobian) {
      int spaceDim  = basisValues.basis_vector.extent(3);
      if(spaceDim==3) {
        dof_functors::EvaluateDOFFastSens_Vector<ScalarT,typename BasisValues2<double>::Array_CellBasisIPDim,3> functor(dof_basis,dof_ip_vector,offsets_array,basisValues.basis_vector);
        Kokkos::parallel_for(workset.num_cells,functor);
      }
      else {
        dof_functors::EvaluateDOFFastSens_Vector<ScalarT,typename BasisValues2<double>::Array_CellBasisIPDim,2> functor(dof_basis,dof_ip_vector,offsets_array,basisValues.basis_vector);
        Kokkos::parallel_for(workset.num_cells,functor);
      }
    }
    else {
      int spaceDim  = basisValues.basis_vector.extent(3);
      if(spaceDim==3) {
        dof_functors::EvaluateDOFWithSens_Vector<ScalarT,typename BasisValues2<double>::Array_CellBasisIPDim,3> functor(dof_basis,dof_ip_vector,basisValues.basis_vector);
        Kokkos::parallel_for(workset.num_cells,functor);
      }
      else {
        dof_functors::EvaluateDOFWithSens_Vector<ScalarT,typename BasisValues2<double>::Array_CellBasisIPDim,2> functor(dof_basis,dof_ip_vector,basisValues.basis_vector);
        Kokkos::parallel_for(workset.num_cells,functor);
      }
    }
  }
  else {
    if(accelerate_jacobian) {
      dof_functors::EvaluateDOFFastSens_Scalar<ScalarT,typename BasisValues2<double>::Array_CellBasisIP> functor(dof_basis,dof_ip_scalar,offsets_array,basisValues.basis_scalar);
      Kokkos::parallel_for(workset.num_cells,functor);
    }
    else {
      dof_functors::EvaluateDOFWithSens_Scalar<ScalarT,typename BasisValues2<double>::Array_CellBasisIP> functor(dof_basis,dof_ip_scalar,basisValues.basis_scalar);
      Kokkos::parallel_for(workset.num_cells,functor);
    }
  }
}

}

#endif
