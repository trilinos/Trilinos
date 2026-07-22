// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_DOF_IMPL_HPP
#define PANZER_DOF_IMPL_HPP

#include <algorithm>
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_BasisIRLayout.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_CommonArrayFactories.hpp"
#include "Panzer_DOF_Functors.hpp"
#include "Panzer_HierarchicParallelism.hpp"

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

  std::string n = "DOF: " + dof_basis.fieldTag().name() + " ("+PHX::print<EvalT>()+")";
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

  std::string n = "DOF: " + dof_basis.fieldTag().name() + " ("+PHX::print<EvalT>()+")";
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

  const auto policy = panzer::HP::inst().teamPolicy<ScalarT,PHX::exec_space>(workset.num_cells);
  const bool use_shared_memory = panzer::HP::inst().useSharedMemory<ScalarT>();

  if(is_vector_basis) {
    using Array=typename BasisValues2<double>::ConstArray_CellBasisIPDim;
    Array array = use_descriptors_ ? basisValues.getVectorBasisValues(false) : Array(basisValues.basis_vector);
    const int spaceDim  = array.extent(3);
    if(spaceDim==3) {
      dof_functors::EvaluateDOFWithSens_Vector<ScalarT,Array,3> functor(dof_basis.get_static_view(),dof_ip_vector.get_static_view(),array,use_shared_memory);
      Kokkos::parallel_for(this->getName(),policy,functor);
    }
    else {
      dof_functors::EvaluateDOFWithSens_Vector<ScalarT,Array,2> functor(dof_basis.get_static_view(),dof_ip_vector.get_static_view(),array,use_shared_memory);
      Kokkos::parallel_for(this->getName(),policy,functor);
    }

  }
  else {
    using Array=typename BasisValues2<double>::ConstArray_CellBasisIP;
    Array interpolation_array = use_descriptors_ ? basisValues.getBasisValues(false) : Array(basisValues.basis_scalar);
    dof_functors::EvaluateDOFWithSens_Scalar<ScalarT,Array> functor(dof_basis,dof_ip_scalar,interpolation_array);
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
    offsets_array = PHX::View<int*>("offsets",offsets.size());
    auto offsets_array_h = Kokkos::create_mirror_view(offsets_array);
    for(std::size_t i=0;i<offsets.size();i++)
      offsets_array_h(i) = offsets[i];
    Kokkos::deep_copy(offsets_array, offsets_array_h);

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
                          + " ("+PHX::print<panzer::Traits::Jacobian>()+")";
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

  std::string n = "DOF: " + dof_basis.fieldTag().name() + " slow_jac(descriptor) ("+PHX::print<typename TRAITS::Jacobian>()+")";
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
      using Array=typename BasisValues2<double>::ConstArray_CellBasisIPDim;
      Array array = use_descriptors_ ? basisValues.getVectorBasisValues(false) : Array(basisValues.basis_vector);
      const int spaceDim  = array.extent(3);
      if(spaceDim==3) {
        dof_functors::EvaluateDOFFastSens_Vector<ScalarT,Array,3> functor(dof_basis,dof_ip_vector,offsets_array,array);
        Kokkos::parallel_for(workset.num_cells,functor);
      }
      else {
        dof_functors::EvaluateDOFFastSens_Vector<ScalarT,Array,2> functor(dof_basis,dof_ip_vector,offsets_array,array);
        Kokkos::parallel_for(workset.num_cells,functor);
      }
    }
    else {
      const bool use_shared_memory = panzer::HP::inst().useSharedMemory<ScalarT>();
      const auto policy = panzer::HP::inst().teamPolicy<ScalarT,PHX::exec_space>(workset.num_cells);
      using Array=typename BasisValues2<double>::ConstArray_CellBasisIPDim;
      Array array = use_descriptors_ ? basisValues.getVectorBasisValues(false) : Array(basisValues.basis_vector);
      const int spaceDim  = array.extent(3);
      if(spaceDim==3) {
        dof_functors::EvaluateDOFWithSens_Vector<ScalarT,Array,3> functor(dof_basis.get_static_view(),dof_ip_vector.get_static_view(),array,use_shared_memory);
	Kokkos::parallel_for(this->getName(),policy,functor);
      }
      else {
        dof_functors::EvaluateDOFWithSens_Vector<ScalarT,Array,2> functor(dof_basis.get_static_view(),dof_ip_vector.get_static_view(),array,use_shared_memory);
	Kokkos::parallel_for(this->getName(),policy,functor);
      }
    }
  }
  else {
    using Array=typename BasisValues2<double>::ConstArray_CellBasisIP;
    Array array = use_descriptors_ ? basisValues.getBasisValues(false) : Array(basisValues.basis_scalar);
    if(accelerate_jacobian) {
      dof_functors::EvaluateDOFFastSens_Scalar<ScalarT,Array> functor(dof_basis,dof_ip_scalar,offsets_array,array);
      Kokkos::parallel_for(workset.num_cells,functor);
    }
    else {
      dof_functors::EvaluateDOFWithSens_Scalar<ScalarT,Array> functor(dof_basis,dof_ip_scalar,array);
      Kokkos::parallel_for(workset.num_cells,functor);
    }
  }
}

}

#endif
