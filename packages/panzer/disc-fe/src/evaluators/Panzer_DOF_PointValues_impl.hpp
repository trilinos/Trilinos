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
//* DOF_PointValues evaluator
//**********************************************************************

//**********************************************************************
// MOST EVALUATION TYPES
//**********************************************************************

//**********************************************************************
template<typename EvalT, typename TRAITS>                   
DOF_PointValues<EvalT, TRAITS>::
DOF_PointValues(const Teuchos::ParameterList & p)
{
  const std::string fieldName = p.get<std::string>("Name");
  basis = p.get< Teuchos::RCP<const PureBasis> >("Basis");
  Teuchos::RCP<const PointRule> pointRule = p.get< Teuchos::RCP<const PointRule> >("Point Rule");
  is_vector_basis = basis->isVectorBasis();

  std::string evalName = fieldName+"_"+pointRule->getName();
  if(p.isType<bool>("Use DOF Name")) {
    if(p.get<bool>("Use DOF Name"))
      evalName = fieldName;
  }

  dof_basis = PHX::MDField<const ScalarT,Cell,Point>(fieldName, basis->functional);

  this->addDependentField(dof_basis);

  // setup all basis fields that are required
  Teuchos::RCP<BasisIRLayout> layout = Teuchos::rcp(new BasisIRLayout(basis,*pointRule));
  basisValues = Teuchos::rcp(new BasisValues2<double>(basis->name()+"_"+pointRule->getName()+"_"));
  basisValues->setupArrays(layout,false);

  // the field manager will allocate all of these field
  // swap between scalar basis value, or vector basis value
  if(basis->isScalarBasis()) {
     dof_ip_scalar = PHX::MDField<ScalarT,Cell,Point>(
                evalName,
     	        pointRule->dl_scalar);
     this->addEvaluatedField(dof_ip_scalar);
     this->addNonConstDependentField(basisValues->basis_ref_scalar);
     this->addNonConstDependentField(basisValues->basis_scalar);
  }
  else if(basis->isVectorBasis()) {
     dof_ip_vector = PHX::MDField<ScalarT,Cell,Point,Dim>(
                evalName,
     	        pointRule->dl_vector);
     this->addEvaluatedField(dof_ip_vector);
     this->addNonConstDependentField(basisValues->basis_ref_vector);
     this->addNonConstDependentField(basisValues->basis_vector);
  }
  else
  { TEUCHOS_ASSERT(false); }

  std::string n = "DOF_PointValues: " + dof_basis.fieldTag().name();
  this->setName(n);
}

//**********************************************************************
template<typename EvalT, typename TRAITS>                   
void DOF_PointValues<EvalT, TRAITS>::
postRegistrationSetup(typename TRAITS::SetupData /* sd */,
                      PHX::FieldManager<TRAITS>& fm)
{
  if(!is_vector_basis) {
    this->utils.setFieldData(basisValues->basis_ref_scalar,fm);
    this->utils.setFieldData(basisValues->basis_scalar,fm);
  }
  else {
    this->utils.setFieldData(basisValues->basis_ref_vector,fm);      
    this->utils.setFieldData(basisValues->basis_vector,fm);           
  }
}

//**********************************************************************
template<typename EvalT, typename TRAITS>                   
void DOF_PointValues<EvalT, TRAITS>::
evaluateFields(typename TRAITS::EvalData workset)
{ 
  const int vector_size = panzer::HP::inst().vectorSize<ScalarT>();

  if(is_vector_basis) {
    int spaceDim  = basisValues->basis_vector.extent(3);
    if(spaceDim==3) {
      dof_functors::EvaluateDOFWithSens_Vector<ScalarT,typename BasisValues2<double>::Array_CellBasisIPDim,3> functor(dof_basis.get_static_view(),dof_ip_vector.get_static_view(),basisValues->basis_vector);
      Kokkos::parallel_for(Kokkos::TeamPolicy<PHX::Device>(workset.num_cells,Kokkos::AUTO(),vector_size),functor);
    }
    else {
      dof_functors::EvaluateDOFWithSens_Vector<ScalarT,typename BasisValues2<double>::Array_CellBasisIPDim,2> functor(dof_basis.get_static_view(),dof_ip_vector.get_static_view(),basisValues->basis_vector);
      Kokkos::parallel_for(Kokkos::TeamPolicy<PHX::Device>(workset.num_cells,Kokkos::AUTO(),vector_size),functor);
    }
  }
  else {
    dof_functors::EvaluateDOFWithSens_Scalar<ScalarT,typename BasisValues2<double>::Array_CellBasisIP> functor(dof_basis,dof_ip_scalar,basisValues->basis_scalar);
    Kokkos::parallel_for(workset.num_cells,functor);
  }
}

//**********************************************************************

//**********************************************************************
// JACOBIAN EVALUATION TYPES
//**********************************************************************

//**********************************************************************
template<typename TRAITS>                   
DOF_PointValues<typename TRAITS::Jacobian, TRAITS>::
DOF_PointValues(const Teuchos::ParameterList & p)
{
  const std::string fieldName = p.get<std::string>("Name");
  basis = p.get< Teuchos::RCP<const PureBasis> >("Basis");
  Teuchos::RCP<const PointRule> pointRule = p.get< Teuchos::RCP<const PointRule> >("Point Rule");
  is_vector_basis = basis->isVectorBasis();

  if(p.isType<Teuchos::RCP<const std::vector<int> > >("Jacobian Offsets Vector")) {
    const std::vector<int> & offsets = *p.get<Teuchos::RCP<const std::vector<int> > >("Jacobian Offsets Vector");

    // allocate and copy offsets vector to Kokkos array
    offsets_array = PHX::View<int*>("offsets",offsets.size());
    for(std::size_t i=0;i<offsets.size();i++)
      offsets_array(i) = offsets[i];

    accelerate_jacobian = true;  // short cut for identity matrix
  }
  else
    accelerate_jacobian = false; // don't short cut for identity matrix

  std::string evalName = fieldName+"_"+pointRule->getName();
  if(p.isType<bool>("Use DOF Name")) {
    if(p.get<bool>("Use DOF Name"))
      evalName = fieldName;
  }

  dof_basis = PHX::MDField<const ScalarT,Cell,Point>(fieldName, basis->functional);

  this->addDependentField(dof_basis);

  // setup all basis fields that are required
  Teuchos::RCP<BasisIRLayout> layout = Teuchos::rcp(new BasisIRLayout(basis,*pointRule));
  basisValues = Teuchos::rcp(new BasisValues2<double>(basis->name()+"_"+pointRule->getName()+"_"));
  basisValues->setupArrays(layout,false);

  // the field manager will allocate all of these field
  // swap between scalar basis value, or vector basis value
  if(basis->isScalarBasis()) {
     dof_ip_scalar = PHX::MDField<ScalarT,Cell,Point>(
                evalName,
     	        pointRule->dl_scalar);
     this->addEvaluatedField(dof_ip_scalar);
     this->addNonConstDependentField(basisValues->basis_ref_scalar);
     this->addNonConstDependentField(basisValues->basis_scalar);
  }
  else if(basis->isVectorBasis()) {
     dof_ip_vector = PHX::MDField<ScalarT,Cell,Point,Dim>(
                evalName,
     	        pointRule->dl_vector);
     this->addEvaluatedField(dof_ip_vector);
     this->addNonConstDependentField(basisValues->basis_ref_vector);
     this->addNonConstDependentField(basisValues->basis_vector);
  }
  else
  { TEUCHOS_ASSERT(false); }

  std::string n = "DOF_PointValues: " + dof_basis.fieldTag().name() + " Jacobian";
  this->setName(n);
}

//**********************************************************************
template<typename TRAITS>                   
void DOF_PointValues<typename TRAITS::Jacobian, TRAITS>::
postRegistrationSetup(typename TRAITS::SetupData /* sd */,
                      PHX::FieldManager<TRAITS>& fm)
{
  if(!is_vector_basis) {
    this->utils.setFieldData(basisValues->basis_ref_scalar,fm);      
    this->utils.setFieldData(basisValues->basis_scalar,fm);           
  }
  else {
    this->utils.setFieldData(basisValues->basis_ref_vector,fm);      
    this->utils.setFieldData(basisValues->basis_vector,fm);           
  }
}

//**********************************************************************
template<typename TRAITS>                   
void DOF_PointValues<typename TRAITS::Jacobian, TRAITS>::
evaluateFields(typename TRAITS::EvalData workset)
{ 
  const int vector_size = panzer::HP::inst().vectorSize<ScalarT>();

  if(is_vector_basis) {
    if(accelerate_jacobian) {
      int spaceDim  = basisValues->basis_vector.extent(3);
      if(spaceDim==3) {
        dof_functors::EvaluateDOFFastSens_Vector<ScalarT,typename BasisValues2<double>::Array_CellBasisIPDim,3> functor(dof_basis,dof_ip_vector,offsets_array,basisValues->basis_vector);
        Kokkos::parallel_for(workset.num_cells,functor);
      }
      else {
        dof_functors::EvaluateDOFFastSens_Vector<ScalarT,typename BasisValues2<double>::Array_CellBasisIPDim,2> functor(dof_basis,dof_ip_vector,offsets_array,basisValues->basis_vector);
        Kokkos::parallel_for(workset.num_cells,functor);
      }
    }
    else {
      int spaceDim  = basisValues->basis_vector.extent(3);
      if(spaceDim==3) {
        dof_functors::EvaluateDOFWithSens_Vector<ScalarT,typename BasisValues2<double>::Array_CellBasisIPDim,3> functor(dof_basis.get_static_view(),dof_ip_vector.get_static_view(),basisValues->basis_vector);
	Kokkos::parallel_for(Kokkos::TeamPolicy<PHX::Device>(workset.num_cells,Kokkos::AUTO(),vector_size),functor);
      }
      else {
        dof_functors::EvaluateDOFWithSens_Vector<ScalarT,typename BasisValues2<double>::Array_CellBasisIPDim,2> functor(dof_basis.get_static_view(),dof_ip_vector.get_static_view(),basisValues->basis_vector);
	Kokkos::parallel_for(Kokkos::TeamPolicy<PHX::Device>(workset.num_cells,Kokkos::AUTO(),vector_size),functor);
      }
    }
  }
  else {
    if(accelerate_jacobian) {
      dof_functors::EvaluateDOFFastSens_Scalar<ScalarT,typename BasisValues2<double>::Array_CellBasisIP> functor(dof_basis,dof_ip_scalar,offsets_array,basisValues->basis_scalar);
      Kokkos::parallel_for(workset.num_cells,functor);
    }
    else {
      dof_functors::EvaluateDOFWithSens_Scalar<ScalarT,typename BasisValues2<double>::Array_CellBasisIP> functor(dof_basis,dof_ip_scalar,basisValues->basis_scalar);
      Kokkos::parallel_for(workset.num_cells,functor);
    }
  }
}

}

#endif
