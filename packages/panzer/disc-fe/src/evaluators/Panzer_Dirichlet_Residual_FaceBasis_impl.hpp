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

#ifndef PANZER_DIRICHLET_RESIDUAL_FACEBASIS_IMPL_HPP
#define PANZER_DIRICHLET_RESIDUAL_FACEBASIS_IMPL_HPP

#include <cstddef>
#include <string>
#include <vector>

#include "Intrepid2_CellTools.hpp"

#include "Kokkos_ViewFactory.hpp"

namespace panzer {

//**********************************************************************
PHX_EVALUATOR_CTOR(DirichletResidual_FaceBasis,p)
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
  TEUCHOS_ASSERT(basis->dimension()==vector_layout_dof->extent_int(2));
  TEUCHOS_ASSERT(vector_layout_vector->dimension(0)==vector_layout_dof->dimension(0));
  TEUCHOS_ASSERT(vector_layout_vector->dimension(1)==vector_layout_dof->dimension(1));
  TEUCHOS_ASSERT(vector_layout_vector->dimension(2)==vector_layout_dof->dimension(2));

  residual = PHX::MDField<ScalarT,Cell,BASIS>(residual_name, basis_layout);
  dof      = PHX::MDField<const ScalarT,Cell,Point,Dim>(dof_name, vector_layout_dof);
  value    = PHX::MDField<const ScalarT,Cell,Point,Dim>(value_name, vector_layout_vector);

  // setup all fields to be evaluated and constructed
  pointValues = PointValues2<ScalarT> (pointRule->getName()+"_",false);
  pointValues.setupArrays(pointRule);

  // the field manager will allocate all of these field
  constJac_ = pointValues.jac;
  this->addDependentField(constJac_);

  
  this->addEvaluatedField(residual);
  this->addDependentField(dof);
  // this->addDependentField(dof_orientation);
  this->addDependentField(value);
 
  std::string n = "Dirichlet Residual Face Basis Evaluator";
  this->setName(n);
}

//**********************************************************************
PHX_POST_REGISTRATION_SETUP(DirichletResidual_FaceBasis, /* worksets */, fm)
{
  this->utils.setFieldData(residual,fm);
  this->utils.setFieldData(dof,fm);
  // this->utils.setFieldData(dof_orientation,fm);
  this->utils.setFieldData(value,fm);
  this->utils.setFieldData(pointValues.jac,fm);

  faceNormal = Kokkos::createDynRankView(residual.get_static_view(),"faceNormal",dof.dimension(0),dof.dimension(1),dof.dimension(2));
}

//**********************************************************************
PHX_EVALUATE_FIELDS(DirichletResidual_FaceBasis,workset)
{ 
  if(workset.num_cells<=0)
    return;
  else {
    Intrepid2::CellTools<PHX::exec_space>::getPhysicalFaceNormals(faceNormal,
                                                                  pointValues.jac.get_view(),
                                                                  this->wda(workset).subcell_index, 
                                                                  *basis->getCellTopology());
  
    for(index_t c=0;c<workset.num_cells;c++) {
      for(int b=0;b<dof.extent_int(1);b++) {
        residual(c,b) = ScalarT(0.0);
        for(int d=0;d<dof.extent_int(2);d++)
          residual(c,b) += (dof(c,b,d)-value(c,b,d))*faceNormal(c,b,d);
      } 
    }
  }
/*
  if(workset.subcell_dim==1) {
    Intrepid2::CellTools<ScalarT>::getPhysicalFaceNormals(faceNormal,
                                            pointValues.jac,
                                            this->wda(workset).subcell_index, 
                                           *basis->getCellTopology());
  
    for(index_t c=0;c<workset.num_cells;c++) {
      for(int b=0;b<dof.extent_int(1);b++) {
        residual(c,b) = ScalarT(0.0);
        for(int d=0;d<dof.extent_int(2);d++)
          residual(c,b) += (dof(c,b,d)-value(c,b,d))*faceNormal(c,b,d);
      } 
    }
  }
  else if(workset.subcell_dim==2) {
    // we need to compute the tangents on each edge for each cell.
    // how do we do this????
    const shards::CellTopology & parentCell = *basis->getCellTopology();
    int cellDim = parentCell.getDimension();
    int numFaces = dof.extent_int(1);

    refFaceNormal = Kokkos::createDynRankView(residual.get_static_view(),"refFaceNormal",numFaces,cellDim);

    for(int i=0;i<numFaces;i++) {
      Kokkos::DynRankView<double,PHX::Device> refFaceNormal_local(cellDim);
      Intrepid2::CellTools<double>::getReferenceFaceNormal(refFaceNormal_local, i, parentCell);

      for(int d=0;d<cellDim;d++) 
        refFaceNormal(i,d) = refFaceNormal_local(d);
    }

    // Loop over workset faces and edge points
    for(index_t c=0;c<workset.num_cells;c++) {
      for(int pt = 0; pt < numFaces; pt++) {

        // Apply parent cell Jacobian to ref. edge tangent
        for(int i = 0; i < cellDim; i++) {
          faceNormal(c, pt, i) = 0.0;
          for(int j = 0; j < cellDim; j++){
            faceNormal(c, pt, i) +=  pointValues.jac(c, pt, i, j)*refFaceNormal(pt,j);
          }// for j
        }// for i
      }// for pt
    }// for pCell

    for(index_t c=0;c<workset.num_cells;c++) {
      for(int b=0;b<dof.extent_int(1);b++) {
        residual(c,b) = ScalarT(0.0);
        for(int d=0;d<dof.extent_int(2);d++)
          residual(c,b) += (dof(c,b,d)-value(c,b,d))*faceNormal(c,b,d);
      } 
    }

  }
  else {
    // don't know what to do 
    TEUCHOS_ASSERT(false);
  }
*/
/*
  // loop over residuals scaling by orientation. This gurantees
  // everything is oriented in the "positive" direction, this allows
  // sums acrossed processor to be oriented in the same way (right?)
  for(index_t c=0;c<workset.num_cells;c++) {
    for(int b=0;b<dof.extent_int(1);b++) {
      residual(c,b) *= dof_orientation(c,b);
    }
  }
*/
}

//**********************************************************************

}

#endif
