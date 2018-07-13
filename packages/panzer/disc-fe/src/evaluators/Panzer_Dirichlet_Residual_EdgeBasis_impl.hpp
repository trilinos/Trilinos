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

#include "Intrepid2_Kernels.hpp"
#include "Intrepid2_CellTools.hpp"
#include "Intrepid2_OrientationTools.hpp"

#include "Phalanx_TypeStrings.hpp"

#include "Panzer_CommonArrayFactories.hpp"
#include "Kokkos_ViewFactory.hpp"

namespace panzer {

//**********************************************************************
template<typename EvalT, typename Traits>
DirichletResidual_EdgeBasis<EvalT, Traits>::
DirichletResidual_EdgeBasis(
  const Teuchos::ParameterList& p)
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
  TEUCHOS_ASSERT(basis_layout->extent(0)==vector_layout_dof->extent(0));
  TEUCHOS_ASSERT(basis_layout->extent(1)==vector_layout_dof->extent(1));
  TEUCHOS_ASSERT(Teuchos::as<unsigned>(basis->dimension())==vector_layout_dof->extent(2));
  TEUCHOS_ASSERT(vector_layout_vector->extent(0)==vector_layout_dof->extent(0));
  TEUCHOS_ASSERT(vector_layout_vector->extent(1)==vector_layout_dof->extent(1));
  TEUCHOS_ASSERT(vector_layout_vector->extent(2)==vector_layout_dof->extent(2));

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
template<typename EvalT, typename Traits>
void
DirichletResidual_EdgeBasis<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& fm)
{
  orientations = sd.orientations_;

  this->utils.setFieldData(residual,fm);
  this->utils.setFieldData(dof,fm);
  this->utils.setFieldData(value,fm);
  this->utils.setFieldData(pointValues.jac,fm);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
DirichletResidual_EdgeBasis<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{ 
  const int numCells = workset.num_cells;
  if(numCells <= 0)
    return;
  else {
    residual.deep_copy(ScalarT(0.0));

    // dofs are already oriented but tangent directions are not oriented

    const int subcellDim = workset.subcell_dim;
    const int subcellOrd = this->wda(workset).subcell_index;
    
    const auto cellTopo = *basis->getCellTopology();
    const auto worksetJacobians = pointValues.jac.get_view();

    const int cellDim = cellTopo.getDimension();

    auto intrepid_basis = basis->getIntrepid2Basis();
    const WorksetDetails & details = workset;

    const bool is_normalize = true;
    auto work = Kokkos::createDynRankView(residual.get_static_view(),"work", 4, cellDim);

    // compute residual
    switch (subcellDim) {
    case 1: {  // 2D element Tri and Quad
      if (intrepid_basis->getDofCount(1, subcellOrd)) {
        auto phyEdgeTan = Kokkos::subview(work, 0, Kokkos::ALL());
        auto ortEdgeTan = Kokkos::subview(work, 1, Kokkos::ALL());
        
        const int ndofsEdge = intrepid_basis->getDofCount(1, subcellOrd);
        const int numEdges = cellTopo.getEdgeCount();
        /* */ int edgeOrts[4] = {};
        for(index_t c=0;c<workset.num_cells;c++) {
          orientations->at(details.cell_local_ids[c]).getEdgeOrientation(edgeOrts, numEdges);
          
          Intrepid2::Orientation::getReferenceEdgeTangent(ortEdgeTan,
                                                          subcellOrd,
                                                          cellTopo,
                                                          edgeOrts[subcellOrd],
                                                          is_normalize);

          for (int i=0;i<ndofsEdge;++i) {
            const int b = intrepid_basis->getDofOrdinal(1, subcellOrd, i);
            auto J = Kokkos::subview(worksetJacobians, c, b, Kokkos::ALL(), Kokkos::ALL());
            Intrepid2::Kernels::Serial::matvec_product(phyEdgeTan, J, ortEdgeTan);            
            
            for(int d=0;d<cellDim;d++) {
              residual(c,b) += (dof(c,b,d)-value(c,b,d))*phyEdgeTan(d);
            }
          }
        }
      }
      break;
    } 
    case 2: { // 3D element Tet and Hex
      const int numEdges = cellTopo.getEdgeCount();
      const int numFaces = cellTopo.getFaceCount();
      
      {
        auto phyEdgeTan = Kokkos::subview(work, 0, Kokkos::ALL());
        auto ortEdgeTan = Kokkos::subview(work, 1, Kokkos::ALL());

        const int numEdgesOfFace= cellTopo.getEdgeCount(2, subcellOrd);

        int edgeOrts[12] = {};
        for(index_t c=0;c<workset.num_cells;c++) {
          for (int i=0;i<numEdgesOfFace;++i) {

            const int edgeOrd = Intrepid2::Orientation::getEdgeOrdinalOfFace(i, subcellOrd, cellTopo);
            const int b = edgeOrd;
            orientations->at(details.cell_local_ids[c]).getEdgeOrientation(edgeOrts, numEdges);
            
            Intrepid2::Orientation::getReferenceEdgeTangent(ortEdgeTan,
                                                            edgeOrd,
                                                            cellTopo,
                                                            edgeOrts[edgeOrd],
                                                            is_normalize);
            
            // for(int b=0;b<dof.extent_int(1);b++) 
            {
              auto J = Kokkos::subview(worksetJacobians, c, b, Kokkos::ALL(), Kokkos::ALL());
              Intrepid2::Kernels::Serial::matvec_product(phyEdgeTan, J, ortEdgeTan);
              
              for(int d=0;d<dof.extent_int(2);d++) {
                residual(c,b) += (dof(c,b,d)-value(c,b,d))*phyEdgeTan(d);
              }
            }
          }
        }
      }

      if (intrepid_basis->getDofCount(2, subcellOrd)) {
        auto phyFaceTanU = Kokkos::subview(work, 0, Kokkos::ALL());
        auto ortFaceTanU = Kokkos::subview(work, 1, Kokkos::ALL());
        auto phyFaceTanV = Kokkos::subview(work, 2, Kokkos::ALL());
        auto ortFaceTanV = Kokkos::subview(work, 3, Kokkos::ALL());
        
        int faceOrts[6] = {};
        for(index_t c=0;c<workset.num_cells;c++) {
          orientations->at(details.cell_local_ids[c]).getFaceOrientation(faceOrts, numFaces);
          Intrepid2::Orientation::getReferenceFaceTangents(ortFaceTanU,
                                                           ortFaceTanV,
                                                           subcellOrd,
                                                           cellTopo,
                                                           faceOrts[subcellOrd],
                                                           is_normalize);
          
          for(int b=0;b<dof.extent_int(1);b++) {
            auto J = Kokkos::subview(worksetJacobians, c, b, Kokkos::ALL(), Kokkos::ALL());
            Intrepid2::Kernels::Serial::matvec_product(phyFaceTanU, J, ortFaceTanU);
            Intrepid2::Kernels::Serial::matvec_product(phyFaceTanV, J, ortFaceTanV);
            
            for(int d=0;d<dof.extent_int(2);d++) {
              residual(c,b) += (dof(c,b,d)-value(c,b,d))*phyFaceTanU(d);
              residual(c,b) += (dof(c,b,d)-value(c,b,d))*phyFaceTanV(d);
            }
          }
        }
      }
      
      break;
    }
    }
  }
  
}
  
//**********************************************************************

}

#endif
