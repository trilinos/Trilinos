// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_GATHER_NORMALS_IMPL_HPP
#define PANZER_GATHER_NORMALS_IMPL_HPP

#include "Teuchos_Assert.hpp"
#include "Phalanx_DataLayout.hpp"

#include "Panzer_PureBasis.hpp"
#include "Kokkos_ViewFactory.hpp"

#include "Intrepid2_Kernels.hpp"
#include "Intrepid2_OrientationTools.hpp"

#include "Teuchos_FancyOStream.hpp"

template<typename EvalT,typename Traits>
panzer::GatherNormals<EvalT, Traits>::
GatherNormals(
  const Teuchos::ParameterList& p)
{
  dof_name_ = (p.get< std::string >("DOF Name"));

  if(p.isType< Teuchos::RCP<PureBasis> >("Basis"))
    basis_ = p.get< Teuchos::RCP<PureBasis> >("Basis");
  else
    basis_ = p.get< Teuchos::RCP<const PureBasis> >("Basis");

  pointRule_ = p.get<Teuchos::RCP<const PointRule> >("Point Rule");

  Teuchos::RCP<PHX::DataLayout> basis_layout         = basis_->functional;
  Teuchos::RCP<PHX::DataLayout> vector_layout_vector = basis_->functional_grad;

  // some sanity checks
  TEUCHOS_ASSERT(basis_->isVectorBasis());

  // setup the orientation field
  std::string orientationFieldName = basis_->name() + " Orientation";
  // setup all fields to be evaluated and constructed
  pointValues_ = panzer::PointValues2<double> (pointRule_->getName()+"_",false);
  pointValues_.setupArrays(pointRule_);

  // the field manager will allocate all of these field
  constJac_ = pointValues_.jac;
  this->addDependentField(constJac_);

  gatherFieldNormals_ = PHX::MDField<ScalarT,Cell,NODE,Dim>(dof_name_+"_Normals",vector_layout_vector);
  this->addEvaluatedField(gatherFieldNormals_);

  this->setName("Gather Normals");
}

// **********************************************************************
template<typename EvalT,typename Traits>
void panzer::GatherNormals<EvalT, Traits>::
postRegistrationSetup(typename Traits::SetupData d,
		      PHX::FieldManager<Traits>& fm)
{
  auto orientations = d.orientations_;
  orientations_ = Kokkos::View<Intrepid2::Orientation*>("orientations_",orientations->size());
  auto orientations_host = Kokkos::create_mirror_view(orientations_);
  for (size_t i=0; i < orientations->size(); ++i)
    orientations_host(i) = (*orientations)[i];
  Kokkos::deep_copy(orientations_,orientations_host);

  this->utils.setFieldData(pointValues_.jac,fm);

  const shards::CellTopology & parentCell = *basis_->getCellTopology();
  int sideDim = parentCell.getDimension()-1;
  sideParam_ = Intrepid2::RefSubcellParametrization<PHX::Device>::get(sideDim, parentCell.getKey());

  int numFaces = gatherFieldNormals_.extent(1);
  keys_ = Kokkos::View<unsigned int*>("parentCell.keys",numFaces);
  auto keys_host = Kokkos::create_mirror_view(keys_);
  for (int i=0; i < numFaces; ++i)
    keys_host(i) = parentCell.getKey(sideDim,i);
  Kokkos::deep_copy(keys_,keys_host);

  // allocate space that is sized correctly for AD
  int cellDim = parentCell.getDimension();
  refEdges_ = Kokkos::createDynRankViewWithType<Kokkos::DynRankView<ScalarT,PHX::Device>>(gatherFieldNormals_.get_static_view(),"ref_edges", (*d.worksets_)[0].num_cells, sideDim, cellDim);
  phyEdges_ = Kokkos::createDynRankViewWithType<Kokkos::DynRankView<ScalarT,PHX::Device>>(gatherFieldNormals_.get_static_view(),"phy_edges", (*d.worksets_)[0].num_cells, sideDim, cellDim);
}

// **********************************************************************
template<typename EvalT,typename Traits>
void panzer::GatherNormals<EvalT, Traits>::
evaluateFields(typename Traits::EvalData workset)
{

  if(workset.num_cells<=0)
    return;

  int numFaces = gatherFieldNormals_.extent(1);
  const auto worksetJacobians = pointValues_.jac.get_view();
  const auto cell_local_ids = workset.getLocalCellIDs();
  auto gatherFieldNormals = gatherFieldNormals_;
  auto sideParam = sideParam_;
  auto keys = keys_;
  auto orientations = orientations_;
  auto refEdges = refEdges_;
  auto phyEdges = phyEdges_;

  // Loop over workset faces and edge points
  Kokkos::parallel_for("panzer::GatherNormals",workset.num_cells,KOKKOS_LAMBDA(const int c){
    int faceOrts[6] = {};
    orientations(cell_local_ids(c)).getFaceOrientation(faceOrts, numFaces);

    for(int pt = 0; pt < numFaces; pt++) {
      auto ortEdgeTan_U = Kokkos::subview(refEdges, c, 0, Kokkos::ALL());
      auto ortEdgeTan_V = Kokkos::subview(refEdges, c, 1, Kokkos::ALL());

      auto tmpRefEdges = Kokkos::subview(refEdges, c, Kokkos::ALL(), Kokkos::ALL());
      Intrepid2::Impl::OrientationTools::getRefSubcellTangents(tmpRefEdges,
                                                               sideParam,
                                                               keys(pt),
                                                               pt,
                                                               faceOrts[pt]);

      auto phyEdgeTan_U = Kokkos::subview(phyEdges, c, 0, Kokkos::ALL());
      auto phyEdgeTan_V = Kokkos::subview(phyEdges, c, 1, Kokkos::ALL());
      auto J = Kokkos::subview(worksetJacobians, c, pt, Kokkos::ALL(), Kokkos::ALL());

      Intrepid2::Kernels::Serial::matvec_product(phyEdgeTan_U, J, ortEdgeTan_U);
      Intrepid2::Kernels::Serial::matvec_product(phyEdgeTan_V, J, ortEdgeTan_V);

      // take the cross product of the two vectors
      gatherFieldNormals(c,pt,0) = (phyEdgeTan_U(1)*phyEdgeTan_V(2) - phyEdgeTan_U(2)*phyEdgeTan_V(1));
      gatherFieldNormals(c,pt,1) = (phyEdgeTan_U(2)*phyEdgeTan_V(0) - phyEdgeTan_U(0)*phyEdgeTan_V(2));
      gatherFieldNormals(c,pt,2) = (phyEdgeTan_U(0)*phyEdgeTan_V(1) - phyEdgeTan_U(1)*phyEdgeTan_V(0));
    }
  });

}

#endif
