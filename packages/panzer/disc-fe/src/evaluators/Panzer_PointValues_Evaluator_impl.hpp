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

#ifndef PANZER_PointValues_Evaluator_IMPL_HPP
#define PANZER_PointValues_Evaluator_IMPL_HPP

#include <algorithm>
#include "Panzer_PointRule.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_CommonArrayFactories.hpp"

namespace panzer {

//**********************************************************************
template<typename EvalT, typename Traits>
PointValues_Evaluator<EvalT, Traits>::
PointValues_Evaluator(
  const Teuchos::ParameterList& p)
{
  basis_index = 0;

  Teuchos::RCP<const panzer::PointRule> pointRule 
     = p.get< Teuchos::RCP<const panzer::PointRule> >("Point Rule");
  Teuchos::RCP<const Kokkos::DynRankView<double,PHX::Device> > userArray
     = p.get<Teuchos::RCP<const Kokkos::DynRankView<double,PHX::Device> > >("Point Array");

  initialize(pointRule,userArray.ptr(),Teuchos::null);
}

//**********************************************************************
template <typename EvalT, typename TRAITST>
PointValues_Evaluator<EvalT,TRAITST>::PointValues_Evaluator(const Teuchos::RCP<const panzer::PointRule> & pointRule,
                                                            const Kokkos::DynRankView<double,PHX::Device> & userArray)
{
  basis_index = 0;

  initialize(pointRule,Teuchos::ptrFromRef(userArray),Teuchos::null);
}

//**********************************************************************
template <typename EvalT, typename TRAITST>
PointValues_Evaluator<EvalT,TRAITST>::PointValues_Evaluator(const Teuchos::RCP<const panzer::PointRule> & pointRule,
                                                            const PHX::MDField<double, panzer::IP, panzer::Dim> & userArray)
{
  basis_index = 0;

  initialize(pointRule,Teuchos::ptrFromRef(userArray),Teuchos::null);
}

//**********************************************************************
template <typename EvalT, typename TRAITST>
PointValues_Evaluator<EvalT,TRAITST>::PointValues_Evaluator(const Teuchos::RCP<const panzer::PointRule> & pointRule,
                                                            const Teuchos::RCP<const panzer::PureBasis> & pureBasis)
{
  basis_index = 0;

  Teuchos::Ptr<const PHX::MDField<double, panzer::IP, panzer::Dim> > userArray;
  initialize(pointRule,userArray,pureBasis);
}

//**********************************************************************
template <typename EvalT, typename TRAITST>
template <typename ArrayT>
void PointValues_Evaluator<EvalT,TRAITST>::initialize(const Teuchos::RCP<const panzer::PointRule> & pointRule,
                                                      const Teuchos::Ptr<const ArrayT> & userArray,
                                                      const Teuchos::RCP<const panzer::PureBasis> & pureBasis)
{
  basis = pureBasis;

  if(userArray!=Teuchos::null && basis==Teuchos::null) 
    useBasisValuesRefArray = false;
  else if(userArray==Teuchos::null && basis!=Teuchos::null) 
    useBasisValuesRefArray = true;
  else {
    // this is a conflicting request, throw an exception
    TEUCHOS_ASSERT(false);
  }

  // copy user array data
  if(userArray!=Teuchos::null) {
    TEUCHOS_ASSERT(userArray->rank()==2);
    MDFieldArrayFactory md_af("refPointArray",true);

    refPointArray = md_af.buildStaticArray<double,NODE,Dim>("refPointArray",userArray->dimension(0),userArray->dimension(1));
    // TEUCHOS_ASSERT(refPointArray.size()==userArray->size());
    for(int i=0;i<userArray->extent_int(0);i++)
      for(int j=0;j<userArray->extent_int(1);j++)
        refPointArray(i,j) = (*userArray)(i,j); 
  }

  // setup all fields to be evaluated and constructed
  pointValues = PointValues2<ScalarT>(pointRule->getName()+"_",false);
  pointValues.setupArrays(pointRule);

  // the field manager will allocate all of these field
  this->addEvaluatedField(pointValues.coords_ref);
  this->addEvaluatedField(pointValues.node_coordinates);
  this->addEvaluatedField(pointValues.jac);
  this->addEvaluatedField(pointValues.jac_inv);
  this->addEvaluatedField(pointValues.jac_det);
  this->addEvaluatedField(pointValues.point_coords);

  std::string n = "PointValues_Evaluator: " + pointRule->getName();
  this->setName(n);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
PointValues_Evaluator<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& fm)
{
  // setup the pointers for evaluation
  this->utils.setFieldData(pointValues.coords_ref,fm);
  this->utils.setFieldData(pointValues.node_coordinates,fm);
  this->utils.setFieldData(pointValues.jac,fm);
  this->utils.setFieldData(pointValues.jac_inv,fm);
  this->utils.setFieldData(pointValues.jac_det,fm);
  this->utils.setFieldData(pointValues.point_coords,fm);

  if(useBasisValuesRefArray) {
    basis_index = panzer::getPureBasisIndex(basis->name(), (*sd.worksets_)[0], this->wda);

    // basis better have coordinates if you want to use them! Assertion to protect
    // a silent failure.
    TEUCHOS_ASSERT(basis->supportsBasisCoordinates());
  }
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
PointValues_Evaluator<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{ 
  if(useBasisValuesRefArray) {
    panzer::BasisValues2<double> & basisValues = *this->wda(workset).bases[basis_index];

    // evaluate the point values (construct jacobians etc...)
    pointValues.evaluateValues(this->wda(workset).cell_vertex_coordinates,
                               basisValues.basis_coordinates_ref);
  }
  else {
    // evaluate the point values (construct jacobians etc...)
    pointValues.evaluateValues(this->wda(workset).cell_vertex_coordinates,refPointArray);
  }
}

//**********************************************************************

}

#endif
