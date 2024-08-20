// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_SubcellSum_impl_hpp__
#define __Panzer_SubcellSum_impl_hpp__

#include "Panzer_PureBasis.hpp"
#include "Panzer_Workset_Utilities.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"

#include "Phalanx_DataLayout_MDALayout.hpp"

namespace panzer {

//**********************************************************************
template<typename EvalT, typename Traits>
SubcellSum<EvalT, Traits>::
SubcellSum(
  const Teuchos::ParameterList& p) 
  : evaluateOnClosure_(false)
{
  Teuchos::RCP<Teuchos::ParameterList> valid_params = this->getValidParameters();
  p.validateParameters(*valid_params);

  const std::string inName = p.get<std::string>("Field Name");
  const std::string outName = p.get<std::string>("Sum Name");
  Teuchos::RCP<const PureBasis> basis = p.get< Teuchos::RCP<const PureBasis> >("Basis");
  multiplier = p.get<double>("Multiplier");
  if(p.isType<bool>("Evaluate On Closure"))
    evaluateOnClosure_ = p.get<bool>("Evaluate On Closure");

  inField = PHX::MDField<const ScalarT,Cell,BASIS>( inName, basis->functional);
  outField = PHX::MDField<ScalarT,Cell>( outName, basis->cell_data);

  this->addDependentField(inField);
  this->addEvaluatedField(outField);

  // build a field pattern object so that looking up closure indices is easy
  fieldPattern_ = Teuchos::rcp(new Intrepid2FieldPattern(basis->getIntrepid2Basis<PHX::exec_space,double,double>()));
    
  std::string n = "SubcellSum: " + outField.fieldTag().name();
  this->setName(n);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
SubcellSum<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{ 
  std::vector<int> indices;
 
  // figure out which indices to sum (this can be made more efficient by 
  // simply saving the indices and only updating if the subcell dimension
  // and index changes)
  if(evaluateOnClosure_)
    fieldPattern_->getSubcellClosureIndices(workset.subcell_dim,this->wda(workset).subcell_index,indices);
  else
    indices = fieldPattern_->getSubcellIndices(workset.subcell_dim,this->wda(workset).subcell_index);

  auto outField_h = Kokkos::create_mirror_view(outField.get_static_view());
  auto inField_h = Kokkos::create_mirror_view(inField.get_static_view());
  Kokkos::deep_copy(inField_h, inField.get_static_view());
  for(index_t c=0;c<workset.num_cells;c++) {
    outField_h(c) = 0.0; // initialize field 

    // sum over all relevant indices for this subcell
    for(std::size_t i=0;i<indices.size();i++)
      outField_h(c) += inField_h(c,indices[i]);
 
    // scale by what ever the user wants
    outField_h(c) *= multiplier;
  }
  Kokkos::deep_copy(outField.get_static_view(), outField_h);
}

//**********************************************************************
template<typename EvalT, typename TRAITS>
Teuchos::RCP<Teuchos::ParameterList> 
SubcellSum<EvalT, TRAITS>::getValidParameters() const
{
  Teuchos::RCP<Teuchos::ParameterList> p = Teuchos::rcp(new Teuchos::ParameterList);
  p->set<std::string>("Sum Name", "?");
  p->set<std::string>("Field Name", "?");
  p->set<double>("Multiplier",1.0);
  p->set<bool>("Evaluate On Closure",false);

  Teuchos::RCP<const panzer::PureBasis> basis;
  p->set("Basis", basis);

  return p;
}

//**********************************************************************

}

#endif
