// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __PANZER_STK_GatherRefCoords_impl_hpp__
#define __PANZER_STK_GatherRefCoords_impl_hpp__

#include "Teuchos_Assert.hpp"
#include "Phalanx_DataLayout.hpp"

#include "Panzer_BasisIRLayout.hpp"

#include "Teuchos_FancyOStream.hpp"

// **********************************************************************
// Specialization: Residual
// **********************************************************************

template<typename EvalT, typename Traits> 
panzer_stk::GatherRefCoords<EvalT, Traits>::
GatherRefCoords(const Teuchos::RCP<const STK_Interface> & mesh,
                const panzer::BasisIRLayout & basis,
                const std::string & fieldName)
{ 
  using panzer::Cell;
  using panzer::NODE;
 
  mesh_ = mesh;

  coordField_ = PHX::MDField<ScalarT,panzer::Cell,panzer::NODE,panzer::Dim>(fieldName,basis.functional_grad);
  this->addEvaluatedField(coordField_);

  this->setName("Gather STK Fields");
}

// **********************************************************************
template<typename EvalT, typename Traits> 
void panzer_stk::GatherRefCoords<EvalT, Traits>::
evaluateFields(typename Traits::EvalData workset)
{ 
   const std::vector<stk::mesh::Entity> & localElements = *mesh_->getElementsOrderedByLID();
 
   // for convenience pull out some objects from workset
   const std::vector<std::size_t> & localCellIds = this->wda(workset).cell_local_ids;

   // convert to a vector of entity objects
   std::vector<stk::mesh::Entity> selected_elements;
   for(std::size_t cell=0;cell<localCellIds.size();cell++)
     selected_elements.push_back(localElements[localCellIds[cell]]);
   
   auto coordField = coordField_.get_view();
   mesh_->getElementVertices_FromCoordsNoResize(selected_elements,coordField);
}

#endif
