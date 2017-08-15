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

#include "Intrepid2_CellTools.hpp"

#include "Panzer_CommonArrayFactories.hpp"
#include "Kokkos_ViewFactory.hpp"

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
  TEUCHOS_ASSERT(Teuchos::as<unsigned>(basis->dimension())==vector_layout_dof->dimension(2));
  TEUCHOS_ASSERT(vector_layout_vector->dimension(0)==vector_layout_dof->dimension(0));
  TEUCHOS_ASSERT(vector_layout_vector->dimension(1)==vector_layout_dof->dimension(1));
  TEUCHOS_ASSERT(vector_layout_vector->dimension(2)==vector_layout_dof->dimension(2));

  residual = PHX::MDField<ScalarT,Cell,BASIS>(residual_name, basis_layout);
  dof      = PHX::MDField<const ScalarT,Cell,Point,Dim>(dof_name, vector_layout_dof);
  value    = PHX::MDField<const ScalarT,Cell,Point,Dim>(value_name, vector_layout_vector);

  // setup all basis fields that are required

  // setup all fields to be evaluated and constructed
  pointValues = PointValues2<ScalarT>(pointRule->getName()+"_",false);
  pointValues.setupArrays(pointRule);

  // the field manager will allocate all of these field
  constJac_ = pointValues.jac;
  this->addDependentField(constJac_);

  
  this->addEvaluatedField(residual);
  this->addDependentField(dof);
  this->addDependentField(value);
 
  std::string n = "Dirichlet Residual Edge Basis Evaluator";
  this->setName(n);
}

//**********************************************************************
PHX_POST_REGISTRATION_SETUP(DirichletResidual_EdgeBasis,sd,fm)
{
  orientations = sd.orientations_;

  this->utils.setFieldData(residual,fm);
  this->utils.setFieldData(dof,fm);
  this->utils.setFieldData(value,fm);
  this->utils.setFieldData(pointValues.jac,fm);

  // C, P, D
  edgeTan  = Kokkos::createDynRankView(residual.get_static_view(),"edgeTan" ,dof.dimension(0),dof.dimension(1),dof.dimension(2));
  faceTanU = Kokkos::createDynRankView(residual.get_static_view(),"faceTanU",dof.dimension(0),dof.dimension(1),dof.dimension(2));
  faceTanV = Kokkos::createDynRankView(residual.get_static_view(),"faceTanV",dof.dimension(0),dof.dimension(1),dof.dimension(2));
}

//**********************************************************************
PHX_EVALUATE_FIELDS(DirichletResidual_EdgeBasis,workset)
{ 
  if(workset.num_cells<=0)
    return;
  else {
    // dofs shoiuld be already account orientations, we don't need to apply orientation separately
    // we need to distinguish if this ordinal is just within edge or face or edge + face.
    const int subcellOrd = this->wda(workset).subcell_index;
    const int edgeOrd = subcellOrd, faceOrd = -1;

    const auto cellTopo = *basis->getCellTopology();
    const auto worksetJacobians = pointValues.jac.get_view();

    const int numEdges = cellTopo.getEdgeCount(), numFaces = cellTopo.getFaceCount();

    auto intrepid_basis = basis->getIntrepid2Basis();
    WorksetDetails & details = workset;
    
    if (edgeOrd > 0 && intrepid_basis->getDofCount(1, edgeOrd)) {
      Intrepid2::CellTools<PHX::exec_space>::getPhysicalEdgeTangents(edgeTan,
                                                                     worksetJacobians,
                                                                     edgeOrd,
                                                                     cellTopo);
      int edgeOrts[12] = {};
      for(index_t c=0;c<workset.num_cells;c++) {
        orientations->at(details.cell_local_ids[c]).getEdgeOrientation(edgeOrts, numEdges);
        const int signVal = edgeOrts[edgeOrd] ? 1.0 : -1.0;
        for(int b=0;b<dof.extent_int(1);b++) {
          residual(c,b) = ScalarT(0.0);
          for(int d=0;d<dof.extent_int(2);d++)
            residual(c,b) += (dof(c,b,d)-value(c,b,d))*edgeTan(c,b,d);
          residual(c,b) *= signVal;
        }
      }
      
    } 
    else if (faceOrd > 0 && intrepid_basis->getDofCount(2, faceOrd)) {
      // Intrepid2::CellTools<PHX::exec_space>::getPhysicalFaceTangents(faceTanU,
      //                                                                faceTanV,
      //                                                                worksetJacobians,
      //                                                                faceOrd,
      //                                                                cellTopo);
      // for(index_t c=0;c<workset.num_cells;c++) {
      //   for(int b=0;b<dof.extent_int(1);b++) {
      //     residual(c,b) = ScalarT(0.0);
      //     for(int d=0;d<dof.extent_int(2);d++) {
      //       residual(c,b) += (dof(c,b,d)-value(c,b,d))*faceTanU(c,b,d);
      //       residual(c,b) += (dof(c,b,d)-value(c,b,d))*faceTanV(c,b,d);
      //     }
      //   }
      // }      
    }
  }
  
}
  
//**********************************************************************

}

#endif
