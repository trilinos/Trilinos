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

#ifndef PANZER_STK_GATHER_FIELDS_IMPL_HPP
#define PANZER_STK_GATHER_FIELDS_IMPL_HPP

#include "Teuchos_Assert.hpp"
#include "Phalanx_DataLayout.hpp"

#include "Panzer_BasisIRLayout.hpp"

#include "Teuchos_FancyOStream.hpp"

// **********************************************************************
// Specialization: Residual
// **********************************************************************

template<typename Traits>
panzer_stk::GatherFields<panzer::Traits::Residual, Traits>::
  GatherFields(const Teuchos::RCP<const STK_Interface> & mesh,const Teuchos::ParameterList& p)
{ 
  using panzer::Cell;
  using panzer::NODE;
 
  mesh_ = mesh;

  const std::vector<std::string>& names = 
    *(p.get< Teuchos::RCP< std::vector<std::string> > >("Field Names"));

  Teuchos::RCP<panzer::BasisIRLayout> basis = 
    p.get< Teuchos::RCP<panzer::BasisIRLayout> >("Basis");

  gatherFields_.resize(names.size());
  stkFields_.resize(names.size());
  for (std::size_t fd = 0; fd < names.size(); ++fd) {
    gatherFields_[fd] = 
      PHX::MDField<ScalarT,Cell,NODE>(names[fd],basis->functional);
    this->addEvaluatedField(gatherFields_[fd]);
  }

  this->setName("Gather STK Fields");
}

// **********************************************************************
template<typename Traits>
void panzer_stk::GatherFields<panzer::Traits::Residual, Traits>::
postRegistrationSetup(typename Traits::SetupData d, 
		      PHX::FieldManager<Traits>& fm)
{
  for (std::size_t fd = 0; fd < gatherFields_.size(); ++fd) {
    std::string fieldName = gatherFields_[fd].fieldTag().name();

    stkFields_[fd] = mesh_->getMetaData()->get_field<VariableField>(fieldName);

    // setup the field data object
    this->utils.setFieldData(gatherFields_[fd],fm);
  }
}

// **********************************************************************
template<typename Traits>
void panzer_stk::GatherFields<panzer::Traits::Residual, Traits>::
evaluateFields(typename Traits::EvalData workset)
{ 
   const std::vector<stk::mesh::Entity*> & localElements = *mesh_->getElementsOrderedByLID();
 
   // for convenience pull out some objects from workset
   const std::vector<std::size_t> & localCellIds = workset.cell_local_ids;
 
   // gather operation for each cell in workset
   for(std::size_t worksetCellIndex=0;worksetCellIndex<localCellIds.size();++worksetCellIndex) {
      std::size_t cellLocalId = localCellIds[worksetCellIndex];
      stk::mesh::PairIterRelation relations = localElements[cellLocalId]->relations(mesh_->getNodeRank());
 
      // loop over the fields to be gathered
      for (std::size_t fieldIndex=0; fieldIndex<gatherFields_.size();fieldIndex++) {
         VariableField * field = stkFields_[fieldIndex];

         std::size_t basisCnt = gatherFields_[fieldIndex].dimension(1);

         // loop over basis functions and fill the fields
         for(std::size_t basis=0;basis<basisCnt;basis++) {
            stk::mesh::Entity * node = relations[basis].entity();
            stk::mesh::EntityArray<VariableField> fieldData(*field,*node);
            (gatherFields_[fieldIndex])(worksetCellIndex,basis) = fieldData(); // from STK
         }
      }
   }
}

#endif
