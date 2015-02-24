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

#ifndef PANZER_DIRICHLET_RESIDUAL_EDGEBASIS_IMPL_HPP
#define PANZER_DIRICHLET_RESIDUAL_EDGEBASIS_IMPL_HPP

#include <cstddef>
#include <string>
#include <vector>

#include "Intrepid_CellTools.hpp"

#include "Panzer_CommonArrayFactories.hpp"

namespace panzer {

//**********************************************************************
PHX_EVALUATOR_CTOR(DirichletResidual_EdgeBasis,p)
{
  std::string residual_name = p.get<std::string>("Residual Name");

  basis = p.get<Teuchos::RCP<const panzer::PureBasis> >("Basis");
  pointRule = p.get<Teuchos::RCP<const panzer::PointRule> >("Point Rule");

  std::string field_name = p.get<std::string>("DOF Name");
  std::string dof_name = field_name+"_"+pointRule->getName(); 
  std::string value_name = p.get<std::string>("Value Name");

  Teuchos::RCP<PHX::DataLayout> basis_layout = basis->functional;
  Teuchos::RCP<PHX::DataLayout> vector_layout_dof    = pointRule->dl_vector;
  Teuchos::RCP<PHX::DataLayout> vector_layout_vector = basis->functional_grad;

  // some sanity checks
  TEUCHOS_ASSERT(basis->isVectorBasis());
  TEUCHOS_ASSERT(basis_layout->dimension(0)==vector_layout_dof->dimension(0));
  TEUCHOS_ASSERT(basis_layout->dimension(1)==vector_layout_dof->dimension(1));
  TEUCHOS_ASSERT(basis->dimension()==vector_layout_dof->dimension(2));
  TEUCHOS_ASSERT(vector_layout_vector->dimension(0)==vector_layout_dof->dimension(0));
  TEUCHOS_ASSERT(vector_layout_vector->dimension(1)==vector_layout_dof->dimension(1));
  TEUCHOS_ASSERT(vector_layout_vector->dimension(2)==vector_layout_dof->dimension(2));

  residual = PHX::MDField<ScalarT,Cell,BASIS>(residual_name, basis_layout);
  dof      = PHX::MDField<ScalarT,Cell,Point,Dim>(dof_name, vector_layout_dof);
  value    = PHX::MDField<ScalarT,Cell,Point,Dim>(value_name, vector_layout_vector);

  // setup the orientation field
  std::string orientationFieldName = basis->name() + " Orientation";
  // std::string orientationFieldName = field_name+" Orientation";
  // if(p.isType<std::string>("Orientation Field Name"))
  //   orientationFieldName = p.get<std::string>("Orientation Field Name");
  dof_orientation = PHX::MDField<ScalarT,Cell,BASIS>(orientationFieldName,
	                                                basis_layout);

  // setup all basis fields that are required
  panzer::MDFieldArrayFactory af_pv(pointRule->getName()+"_");

  // setup all fields to be evaluated and constructed
  pointValues.setupArrays(pointRule,af_pv);

  // the field manager will allocate all of these field
  this->addDependentField(pointValues.jac);

  
  this->addEvaluatedField(residual);
  this->addDependentField(dof);
  this->addDependentField(dof_orientation);
  this->addDependentField(value);
  this->addDependentField(pointValues.jac);
 
  std::string n = "Dirichlet Residual Edge Basis Evaluator";
  this->setName(n);
}

//**********************************************************************
PHX_POST_REGISTRATION_SETUP(DirichletResidual_EdgeBasis,worksets,fm)
{
  this->utils.setFieldData(residual,fm);
  this->utils.setFieldData(dof,fm);
  this->utils.setFieldData(dof_orientation,fm);
  this->utils.setFieldData(value,fm);
  this->utils.setFieldData(pointValues.jac,fm);

  edgeTan = Intrepid::FieldContainer<ScalarT>(dof.dimension(0),dof.dimension(1),dof.dimension(2));
}

//**********************************************************************
PHX_EVALUATE_FIELDS(DirichletResidual_EdgeBasis,workset)
{ 
  if(workset.num_cells<=0)
    return;

  residual.deep_copy(ScalarT(0.0));

  if(workset.subcell_dim==1) {
    Intrepid::CellTools<ScalarT>::getPhysicalEdgeTangents(edgeTan,
                                            pointValues.jac,
                                            workset.subcell_index, 
                                           *basis->getCellTopology());

    for(std::size_t c=0;c<workset.num_cells;c++) {
      for(int b=0;b<dof.dimension(1);b++) {
        for(int d=0;d<dof.dimension(2);d++)
          residual(c,b) += (dof(c,b,d)-value(c,b,d))*edgeTan(c,b,d);
      } 
    }
  }
  else if(workset.subcell_dim==2) {
    // we need to compute the tangents on each edge for each cell.
    // how do we do this????
    const shards::CellTopology & parentCell = *basis->getCellTopology();
    int cellDim = parentCell.getDimension();
    int numEdges = dof.dimension(1);

    refEdgeTan = Intrepid::FieldContainer<ScalarT>(numEdges,cellDim);

    for(int i=0;i<numEdges;i++) {
      Intrepid::FieldContainer<double> refEdgeTan_local(cellDim);
      Intrepid::CellTools<double>::getReferenceEdgeTangent(refEdgeTan_local, i, parentCell);

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

    for(std::size_t c=0;c<workset.num_cells;c++) {
      for(int b=0;b<dof.dimension(1);b++) {
        for(int d=0;d<dof.dimension(2);d++)
          residual(c,b) += (dof(c,b,d)-value(c,b,d))*edgeTan(c,b,d);
      } 
    }

  }
  else {
    // don't know what to do 
    TEUCHOS_ASSERT(false);
  }

  // loop over residuals scaling by orientation. This gurantees
  // everything is oriented in the "positive" direction, this allows
  // sums acrossed processor to be oriented in the same way (right?)
  for(std::size_t c=0;c<workset.num_cells;c++) {
    for(int b=0;b<dof.dimension(1);b++) {
      residual(c,b) *= dof_orientation(c,b);
    }
  }
}

//**********************************************************************

}

#endif
