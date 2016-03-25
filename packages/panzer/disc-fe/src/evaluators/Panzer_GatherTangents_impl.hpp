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

#ifndef PANZER_GATHER_TANGENTS_IMPL_HPP
#define PANZER_GATHER_TANGENTS_IMPL_HPP

#include "Teuchos_Assert.hpp"
#include "Phalanx_DataLayout.hpp"

#include "Panzer_PureBasis.hpp"
#include "Panzer_CommonArrayFactories.hpp"

#include "Teuchos_FancyOStream.hpp"

template<typename EvalT,typename Traits>
panzer::GatherTangents<EvalT, Traits>::
GatherTangents(
  const Teuchos::ParameterList& p)
{ 
  dof_name = (p.get< std::string >("DOF Name"));

  if(p.isType< Teuchos::RCP<PureBasis> >("Basis"))
    basis = p.get< Teuchos::RCP<PureBasis> >("Basis");
  else
    basis = p.get< Teuchos::RCP<const PureBasis> >("Basis");

  pointRule = p.get<Teuchos::RCP<const PointRule> >("Point Rule");

  Teuchos::RCP<PHX::DataLayout> basis_layout         = basis->functional;
  Teuchos::RCP<PHX::DataLayout> vector_layout_vector = basis->functional_grad;

  // some sanity checks
  TEUCHOS_ASSERT(basis->isVectorBasis());

  // setup the orientation field
  std::string orientationFieldName = basis->name() + " Orientation";
  dof_orientation = PHX::MDField<ScalarT,Cell,NODE>(orientationFieldName, basis_layout);

  // setup all basis fields that are required
  MDFieldArrayFactory af_pv(pointRule->getName()+"_");

  // setup all fields to be evaluated and constructed
  pointValues.setupArrays(pointRule,af_pv);

  // the field manager will allocate all of these field
  this->addDependentField(dof_orientation);
  this->addDependentField(pointValues.jac);

  gatherFieldTangents = PHX::MDField<ScalarT,Cell,NODE,Dim>(dof_name+"_Tangents",vector_layout_vector);
  this->addEvaluatedField(gatherFieldTangents);

  this->setName("Gather Tangents");
}

// **********************************************************************
template<typename EvalT,typename Traits>
void panzer::GatherTangents<EvalT, Traits>::
postRegistrationSetup(typename Traits::SetupData d, 
		      PHX::FieldManager<Traits>& fm)
{
  // setup the field data object
  this->utils.setFieldData(gatherFieldTangents,fm);
  this->utils.setFieldData(dof_orientation,fm);
  this->utils.setFieldData(pointValues.jac,fm);

  edgeTan = Intrepid2::FieldContainer<ScalarT>(gatherFieldTangents.dimension(0),gatherFieldTangents.dimension(1),gatherFieldTangents.dimension(2));
}

// **********************************************************************
template<typename EvalT,typename Traits>
void panzer::GatherTangents<EvalT, Traits>::
evaluateFields(typename Traits::EvalData workset)
{ 

  if(workset.num_cells<=0)
    return;
  else {
    const shards::CellTopology & parentCell = *basis->getCellTopology();
    int cellDim = parentCell.getDimension();
    int numEdges = gatherFieldTangents.dimension(1);

    refEdgeTan = Intrepid2::FieldContainer<ScalarT>(numEdges,cellDim);

    for(int i=0;i<numEdges;i++) {
      Intrepid2::FieldContainer<double> refEdgeTan_local(cellDim);
      Intrepid2::CellTools<double>::getReferenceEdgeTangent(refEdgeTan_local, i, parentCell);

      for(int d=0;d<cellDim;d++)
        refEdgeTan(i,d) = refEdgeTan_local(d);
    }

    // Loop over workset faces and edge points
    for(std::size_t c=0;c<workset.num_cells;c++) {
      for(int pt = 0; pt < numEdges; pt++) {

        // Apply parent cell Jacobian to ref. edge tangent
        for(int i = 0; i < cellDim; i++) {
          edgeTan(c, pt, i) = 0.0;
          for(int j = 0; j < cellDim; j++){
            edgeTan(c, pt, i) +=  pointValues.jac(c, pt, i, j)*refEdgeTan(pt,j);
          }// for j
        }// for i
      }// for pt
    }// for pCell

    // Multiply tangent by orientation
    for(std::size_t c=0;c<workset.num_cells;c++) {
      for(int b=0;b<gatherFieldTangents.dimension(1);b++) {
        for(int d=0;d<gatherFieldTangents.dimension(2);d++) {
          gatherFieldTangents(c,b,d) = edgeTan(c,b,d)*dof_orientation(c,b); 
        }
      }
    }
  }

}

#endif
