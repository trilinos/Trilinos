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

#include "Intrepid_FunctionSpaceTools.hpp"

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

  dof_basis = PHX::MDField<ScalarT,Cell,Point>(fieldName, basis->functional);

  this->addDependentField(dof_basis);

  // setup all basis fields that are required
  Teuchos::RCP<BasisIRLayout> layout = Teuchos::rcp(new BasisIRLayout(basis,*pointRule));
  basisValues = Teuchos::rcp(new BasisValues2<ScalarT>(basis->name()+"_"+pointRule->getName()+"_"));
  basisValues->setupArrays(layout,false);

  // the field manager will allocate all of these field
  // swap between scalar basis value, or vector basis value
  if(basis->isScalarBasis()) {
     dof_ip_scalar = PHX::MDField<ScalarT,Cell,Point>(
                evalName,
     	        pointRule->dl_scalar);
     this->addEvaluatedField(dof_ip_scalar);

     this->addDependentField(basisValues->basis_ref_scalar); 
     this->addDependentField(basisValues->basis_scalar); 
  }
  else if(basis->isVectorBasis()) {
     dof_ip_vector = PHX::MDField<ScalarT,Cell,Point,Dim>(
                evalName,
     	        pointRule->dl_vector);
     this->addEvaluatedField(dof_ip_vector);

     this->addDependentField(basisValues->basis_ref_vector); 
     this->addDependentField(basisValues->basis_vector); 
  }
  else
  { TEUCHOS_ASSERT(false); }

  std::string n = "DOF_PointValues: " + dof_basis.fieldTag().name();
  this->setName(n);
}

//**********************************************************************
template<typename EvalT, typename TRAITS>                   
void DOF_PointValues<EvalT, TRAITS>::
postRegistrationSetup(typename TRAITS::SetupData sd,
                      PHX::FieldManager<TRAITS>& fm)
{
  this->utils.setFieldData(dof_basis,fm);

  if(!is_vector_basis) {
    this->utils.setFieldData(dof_ip_scalar,fm);

    // setup the pointers for the basis values data structure
    this->utils.setFieldData(basisValues->basis_ref_scalar,fm);      
    this->utils.setFieldData(basisValues->basis_scalar,fm);           
  }
  else {
    this->utils.setFieldData(dof_ip_vector,fm);

    // setup the pointers for the basis values data structure
    this->utils.setFieldData(basisValues->basis_ref_vector,fm);      
    this->utils.setFieldData(basisValues->basis_vector,fm);           
  }
}

//**********************************************************************
template<typename EvalT, typename TRAITS>                   
void DOF_PointValues<EvalT, TRAITS>::
evaluateFields(typename TRAITS::EvalData workset)
{ 
  // evaluateDOF_withSens(dof_basis,dof_ip,dof_orientation,is_vector_basis,workset.num_cells,basisValues.basis);

  if(is_vector_basis) {
    int spaceDim  = basisValues->basis_vector.dimension(3);
    if(spaceDim==3) {
      dof_functors::EvaluateDOFWithSens_Vector<ScalarT,typename BasisValues2<ScalarT>::Array_CellBasisIPDim,3> functor(dof_basis,dof_ip_vector,basisValues->basis_vector);
      Kokkos::parallel_for(workset.num_cells,functor);
    }
    else {
      dof_functors::EvaluateDOFWithSens_Vector<ScalarT,typename BasisValues2<ScalarT>::Array_CellBasisIPDim,2> functor(dof_basis,dof_ip_vector,basisValues->basis_vector);
      Kokkos::parallel_for(workset.num_cells,functor);
    }
  }
  else {
    dof_functors::EvaluateDOFWithSens_Scalar<ScalarT,typename BasisValues2<ScalarT>::Array_CellBasisIP> functor(dof_basis,dof_ip_scalar,basisValues->basis_scalar);
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
    offsets_array = Kokkos::View<int*,PHX::Device>("offsets",offsets.size());
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

  dof_basis = PHX::MDField<ScalarT,Cell,Point>(fieldName, basis->functional);

  this->addDependentField(dof_basis);

  // setup all basis fields that are required
  Teuchos::RCP<BasisIRLayout> layout = Teuchos::rcp(new BasisIRLayout(basis,*pointRule));
  basisValues = Teuchos::rcp(new BasisValues2<ScalarT>(basis->name()+"_"+pointRule->getName()+"_"));
  basisValues->setupArrays(layout,false);

  // the field manager will allocate all of these field
  // swap between scalar basis value, or vector basis value
  if(basis->isScalarBasis()) {
     dof_ip_scalar = PHX::MDField<ScalarT,Cell,Point>(
                evalName,
     	        pointRule->dl_scalar);
     this->addEvaluatedField(dof_ip_scalar);

     this->addDependentField(basisValues->basis_ref_scalar); 
     this->addDependentField(basisValues->basis_scalar); 
  }
  else if(basis->isVectorBasis()) {
     dof_ip_vector = PHX::MDField<ScalarT,Cell,Point,Dim>(
                evalName,
     	        pointRule->dl_vector);
     this->addEvaluatedField(dof_ip_vector);

     this->addDependentField(basisValues->basis_ref_vector); 
     this->addDependentField(basisValues->basis_vector); 
  }
  else
  { TEUCHOS_ASSERT(false); }

  std::string n = "DOF_PointValues: " + dof_basis.fieldTag().name() + " Jacobian";
  this->setName(n);
}

//**********************************************************************
template<typename TRAITS>                   
void DOF_PointValues<typename TRAITS::Jacobian, TRAITS>::
postRegistrationSetup(typename TRAITS::SetupData sd,
                      PHX::FieldManager<TRAITS>& fm)
{
  this->utils.setFieldData(dof_basis,fm);

  if(!is_vector_basis) {
    this->utils.setFieldData(dof_ip_scalar,fm);

    // setup the pointers for the basis values data structure
    this->utils.setFieldData(basisValues->basis_ref_scalar,fm);      
    this->utils.setFieldData(basisValues->basis_scalar,fm);           
  }
  else {
    this->utils.setFieldData(dof_ip_vector,fm);

    // setup the pointers for the basis values data structure
    this->utils.setFieldData(basisValues->basis_ref_vector,fm);      
    this->utils.setFieldData(basisValues->basis_vector,fm);           
  }
}

//**********************************************************************
template<typename TRAITS>                   
void DOF_PointValues<typename TRAITS::Jacobian, TRAITS>::
evaluateFields(typename TRAITS::EvalData workset)
{ 
  if(is_vector_basis) {
    if(accelerate_jacobian) {
      int spaceDim  = basisValues->basis_vector.dimension(3);
      if(spaceDim==3) {
        dof_functors::EvaluateDOFFastSens_Vector<ScalarT,typename BasisValues2<ScalarT>::Array_CellBasisIPDim,3> functor(dof_basis,dof_ip_vector,offsets_array,basisValues->basis_vector);
        Kokkos::parallel_for(workset.num_cells,functor);
      }
      else {
        dof_functors::EvaluateDOFFastSens_Vector<ScalarT,typename BasisValues2<ScalarT>::Array_CellBasisIPDim,2> functor(dof_basis,dof_ip_vector,offsets_array,basisValues->basis_vector);
        Kokkos::parallel_for(workset.num_cells,functor);
      }
    }
    else {
      int spaceDim  = basisValues->basis_vector.dimension(3);
      if(spaceDim==3) {
        dof_functors::EvaluateDOFWithSens_Vector<ScalarT,typename BasisValues2<ScalarT>::Array_CellBasisIPDim,3> functor(dof_basis,dof_ip_vector,basisValues->basis_vector);
        Kokkos::parallel_for(workset.num_cells,functor);
      }
      else {
        dof_functors::EvaluateDOFWithSens_Vector<ScalarT,typename BasisValues2<ScalarT>::Array_CellBasisIPDim,2> functor(dof_basis,dof_ip_vector,basisValues->basis_vector);
        Kokkos::parallel_for(workset.num_cells,functor);
      }
    }
  }
  else {
    if(accelerate_jacobian) {
      dof_functors::EvaluateDOFFastSens_Scalar<ScalarT,typename BasisValues2<ScalarT>::Array_CellBasisIP> functor(dof_basis,dof_ip_scalar,offsets_array,basisValues->basis_scalar);
      Kokkos::parallel_for(workset.num_cells,functor);
    }
    else {
      dof_functors::EvaluateDOFWithSens_Scalar<ScalarT,typename BasisValues2<ScalarT>::Array_CellBasisIP> functor(dof_basis,dof_ip_scalar,basisValues->basis_scalar);
      Kokkos::parallel_for(workset.num_cells,functor);
    }
  }
}

}

#endif
