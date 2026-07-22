// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_DIRICHLET_RESIDUAL_EDGEBASIS_IMPL_HPP
#define PANZER_DIRICHLET_RESIDUAL_EDGEBASIS_IMPL_HPP

#include <cstddef>
#include <string>
#include <vector>

#include "Intrepid2_Kernels.hpp"
#include "Intrepid2_CellTools.hpp"
#include "Intrepid2_OrientationTools.hpp"

#include "Phalanx_Print.hpp"

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
  pointValues = PointValues2<double>(pointRule->getName()+"_",false);
  pointValues.setupArrays(pointRule);

  // the field manager will allocate all of these field
  this->addNonConstDependentField(pointValues.jac);
  
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
  this->utils.setFieldData(pointValues.jac,fm);

  const auto cellTopo = *basis->getCellTopology();
  const int edgeDim = 1;
  const int faceDim = 2;
  if(cellTopo.getDimension() > edgeDim)
    edgeParam = Intrepid2::RefSubcellParametrization<Kokkos::HostSpace>::get(edgeDim, cellTopo.getKey());

  if(cellTopo.getDimension() > faceDim)
    faceParam = Intrepid2::RefSubcellParametrization<Kokkos::HostSpace>::get(faceDim, cellTopo.getKey());
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
    const int edgeDim = 1;
    const int faceDim = 2;

    auto intrepid_basis = basis->getIntrepid2Basis();
    const WorksetDetails & details = workset;

    //const bool is_normalize = true;
    auto work = Kokkos::create_mirror_view(Kokkos::createDynRankView(residual.get_static_view(),"work", 4, cellDim));

    // compute residual
    auto residual_h = Kokkos::create_mirror_view(residual.get_static_view());
    auto dof_h = Kokkos::create_mirror_view(dof.get_static_view());
    auto value_h = Kokkos::create_mirror_view(value.get_static_view());
    Kokkos::deep_copy(dof_h, dof.get_static_view());
    Kokkos::deep_copy(value_h, value.get_static_view());
    auto worksetJacobians_h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),worksetJacobians);
    switch (subcellDim) {
    case 1: {  // 2D element Tri and Quad
      if (intrepid_basis->getDofCount(1, subcellOrd)) {
        auto ortEdgeTan = Kokkos::subview(work, 0, Kokkos::ALL());
        auto phyEdgeTan = Kokkos::subview(work, 1, Kokkos::ALL());
        
        const int ndofsEdge = intrepid_basis->getDofCount(1, subcellOrd);
        const int numEdges = cellTopo.getEdgeCount();
        /* */ int edgeOrts[4] = {};
	
        for(index_t c=0;c<workset.num_cells;c++) {
          orientations->at(details.cell_local_ids[c]).getEdgeOrientation(edgeOrts, numEdges);
          
          Intrepid2::Impl::OrientationTools::getRefSubcellTangents(work,
              edgeParam,
              cellTopo.getKey(edgeDim,subcellOrd),
              subcellOrd,
              edgeOrts[subcellOrd]);

          for (int i=0;i<ndofsEdge;++i) {
            const int b = intrepid_basis->getDofOrdinal(1, subcellOrd, i);
            auto J = Kokkos::subview(worksetJacobians_h, c, b, Kokkos::ALL(), Kokkos::ALL());
            Intrepid2::Kernels::Serial::matvec_product(phyEdgeTan, J, ortEdgeTan);            
            
            for(int d=0;d<cellDim;d++) {
              residual_h(c,b) += (dof_h(c,b,d)-value_h(c,b,d))*phyEdgeTan(d);
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

        auto ortEdgeTan = Kokkos::subview(work, 0, Kokkos::ALL());
        auto phyEdgeTan = Kokkos::subview(work, 1, Kokkos::ALL());

        const int numEdgesOfFace= cellTopo.getEdgeCount(2, subcellOrd);

        int edgeOrts[12] = {};
        for(index_t c=0;c<workset.num_cells;c++) {
          for (int i=0;i<numEdgesOfFace;++i) {

            const int edgeOrd = Intrepid2::Orientation::getEdgeOrdinalOfFace(i, subcellOrd, cellTopo);
            const int b = edgeOrd;
            orientations->at(details.cell_local_ids[c]).getEdgeOrientation(edgeOrts, numEdges);
            
            Intrepid2::Impl::OrientationTools::getRefSubcellTangents(work,
                edgeParam,
                cellTopo.getKey(edgeDim,edgeOrd),
                edgeOrd,
                edgeOrts[edgeOrd]);

            {
              auto J = Kokkos::subview(worksetJacobians_h, c, b, Kokkos::ALL(), Kokkos::ALL());
              Intrepid2::Kernels::Serial::matvec_product(phyEdgeTan, J, ortEdgeTan);
              
              for(int d=0;d<dof.extent_int(2);d++) {
                residual_h(c,b) += (dof_h(c,b,d)-value_h(c,b,d))*phyEdgeTan(d);
              }
            }
          }
        }
      }

      if (intrepid_basis->getDofCount(2, subcellOrd)) {
        auto ortFaceTanU = Kokkos::subview(work, 0, Kokkos::ALL());
        auto ortFaceTanV = Kokkos::subview(work, 1, Kokkos::ALL());
        auto phyFaceTanU = Kokkos::subview(work, 2, Kokkos::ALL());
        auto phyFaceTanV = Kokkos::subview(work, 3, Kokkos::ALL());
        
        int faceOrts[6] = {};
        for(index_t c=0;c<workset.num_cells;c++) {
          orientations->at(details.cell_local_ids[c]).getFaceOrientation(faceOrts, numFaces);

          Intrepid2::Impl::OrientationTools::getRefSubcellTangents(work,
              faceParam,
              cellTopo.getKey(faceDim,subcellOrd),
              subcellOrd,
              faceOrts[subcellOrd]);

          for(int b=0;b<dof.extent_int(1);b++) {
            auto J = Kokkos::subview(worksetJacobians, c, b, Kokkos::ALL(), Kokkos::ALL());
            Intrepid2::Kernels::Serial::matvec_product(phyFaceTanU, J, ortFaceTanU);
            Intrepid2::Kernels::Serial::matvec_product(phyFaceTanV, J, ortFaceTanV);
            
            for(int d=0;d<dof.extent_int(2);d++) {
              residual_h(c,b) += (dof_h(c,b,d)-value_h(c,b,d))*phyFaceTanU(d);
              residual_h(c,b) += (dof_h(c,b,d)-value_h(c,b,d))*phyFaceTanV(d);
            }
          }
        }
      }
      
      break;
    }
    }
    Kokkos::deep_copy(residual.get_static_view(), residual_h);
  }
  
}
  
//**********************************************************************

}

#endif
