// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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

    refPointArray = md_af.buildStaticArray<double,NODE,Dim>("refPointArray",userArray->extent(0),userArray->extent(1));
    Kokkos::deep_copy(PHX::as_view(refPointArray), PHX::as_view(*userArray));

  }

  // setup all fields to be evaluated and constructed
  pointValues = PointValues2<double>(pointRule->getName()+"_",false);
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
    pointValues.evaluateValues(this->wda(workset).cell_node_coordinates,
                               basisValues.basis_coordinates_ref,
                               workset.num_cells);
  }
  else {
    // evaluate the point values (construct jacobians etc...)
    pointValues.evaluateValues(this->wda(workset).cell_node_coordinates,refPointArray,
                               workset.num_cells);
  }
}

//**********************************************************************

}

#endif
