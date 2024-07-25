// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_GATHER_TANGENTS_IMPL_HPP
#define PANZER_GATHER_TANGENTS_IMPL_HPP

#include "Teuchos_Assert.hpp"
#include "Phalanx_DataLayout.hpp"

#include "Intrepid2_Kernels.hpp"
#include "Intrepid2_OrientationTools.hpp"

#include "Panzer_PureBasis.hpp"
#include "Kokkos_ViewFactory.hpp"

#include "Teuchos_FancyOStream.hpp"

template<typename EvalT,typename Traits>
panzer::GatherTangents<EvalT, Traits>::
GatherTangents(
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

  gatherFieldTangents_ = PHX::MDField<ScalarT,Cell,NODE,Dim>(dof_name_+"_Tangents",vector_layout_vector);
  this->addEvaluatedField(gatherFieldTangents_);

  this->setName("Gather Tangents");
}

// **********************************************************************
template<typename EvalT,typename Traits>
void panzer::GatherTangents<EvalT, Traits>::
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
  const int edgeDim = 1;
  edgeParam_ = Intrepid2::RefSubcellParametrization<PHX::Device>::get(edgeDim, parentCell.getKey());

  const int numEdges = gatherFieldTangents_.extent(1);
  keys_ = Kokkos::View<unsigned int*>("parentCell.keys",numEdges);
  auto keys_host = Kokkos::create_mirror_view(keys_);
  for (int i=0; i < numEdges; ++i)
    keys_host(i) = parentCell.getKey(edgeDim,i);
  Kokkos::deep_copy(keys_,keys_host);
}

// **********************************************************************
template<typename EvalT,typename Traits>
void panzer::GatherTangents<EvalT, Traits>::
evaluateFields(typename Traits::EvalData workset)
{ 

  if(workset.num_cells<=0)
    return;
  else {
    const shards::CellTopology & parentCell = *basis_->getCellTopology();
    const int cellDim = parentCell.getDimension();
    const int numEdges = gatherFieldTangents_.extent(1);

    auto workspace_tmp = Kokkos::createDynRankView(gatherFieldTangents_.get_static_view(),"workspace", workset.num_cells,4, cellDim);
    const auto worksetJacobians = pointValues_.jac.get_view();
    const auto cell_local_ids = workset.getLocalCellIDs();
    auto gatherFieldTangents = gatherFieldTangents_.get_static_view();
    auto keys = keys_;
    auto orientations = orientations_;
    auto edgeParam = edgeParam_;

    // Loop over workset faces and edge points
    Kokkos::parallel_for("panzer::GatherTangets",workset.num_cells,KOKKOS_LAMBDA(const int c){
      int edgeOrts[12] = {};
      orientations(cell_local_ids(c)).getEdgeOrientation(edgeOrts, numEdges);

      auto workspace = Kokkos::subview(workspace_tmp,c,Kokkos::ALL(),Kokkos::ALL());
      for(int pt = 0; pt < numEdges; pt++) {
        auto phyEdgeTan = Kokkos::subview(gatherFieldTangents, c, pt, Kokkos::ALL());
        auto ortEdgeTan = Kokkos::subview(workspace_tmp,c,0,Kokkos::ALL());

        Intrepid2::Impl::OrientationTools::getRefSubcellTangents(
            workspace,
            edgeParam,
            keys(pt),
            pt,
            edgeOrts[pt]);

        auto J = Kokkos::subview(worksetJacobians, c, pt, Kokkos::ALL(), Kokkos::ALL());
        Intrepid2::Kernels::Serial::matvec_product(phyEdgeTan, J, ortEdgeTan);            

      }// for pt
    });// for pCell

  }

}

#endif
