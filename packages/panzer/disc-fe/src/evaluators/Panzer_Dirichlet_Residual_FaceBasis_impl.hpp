// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
template<typename EvalT, typename Traits>
DirichletResidual_FaceBasis<EvalT, Traits>::
DirichletResidual_FaceBasis(
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
  TEUCHOS_ASSERT(basis->dimension()==vector_layout_dof->extent_int(2));
  TEUCHOS_ASSERT(vector_layout_vector->extent(0)==vector_layout_dof->extent(0));
  TEUCHOS_ASSERT(vector_layout_vector->extent(1)==vector_layout_dof->extent(1));
  TEUCHOS_ASSERT(vector_layout_vector->extent(2)==vector_layout_dof->extent(2));

  residual = PHX::MDField<ScalarT,Cell,BASIS>(residual_name, basis_layout);
  dof      = PHX::MDField<const ScalarT,Cell,Point,Dim>(dof_name, vector_layout_dof);
  value    = PHX::MDField<const ScalarT,Cell,Point,Dim>(value_name, vector_layout_vector);

  // setup all fields to be evaluated and constructed
  pointValues = PointValues2<double> (pointRule->getName()+"_",false);
  pointValues.setupArrays(pointRule);

  // the field manager will allocate all of these field
  constJac_ = pointValues.jac;
  this->addDependentField(constJac_);

  
  this->addEvaluatedField(residual);
  this->addDependentField(dof);
  this->addDependentField(value);
 
  std::string n = "Dirichlet Residual Face Basis Evaluator";
  this->setName(n);
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
DirichletResidual_FaceBasis<EvalT, Traits>::
postRegistrationSetup(
  typename Traits::SetupData sd,
  PHX::FieldManager<Traits>& fm)
{
  orientations = sd.orientations_;
  this->utils.setFieldData(pointValues.jac,fm);
  faceNormal = Kokkos::createDynRankView(residual.get_static_view(),"faceNormal",dof.extent(0),dof.extent(1),dof.extent(2));
}

//**********************************************************************
template<typename EvalT, typename Traits>
void
DirichletResidual_FaceBasis<EvalT, Traits>::
evaluateFields(
  typename Traits::EvalData workset)
{ 
  // basic cell topology data
  const shards::CellTopology & parentCell = *basis->getCellTopology();
  const int cellDim = parentCell.getDimension();

  // face only, subcellOrd does not include edge count
  const int subcellDim = 2;
  const int subcellOrd = this->wda(workset).subcell_index;

  const int numFaces = parentCell.getSubcellCount(subcellDim);  
  const int numFaceDofs = dof.extent_int(1);

  TEUCHOS_ASSERT(cellDim == dof.extent_int(2));

  if(workset.num_cells<=0)
    return;
  else {
    Intrepid2::CellTools<PHX::exec_space>::getPhysicalFaceNormals(faceNormal,
                                                                  pointValues.jac.get_view(),
                                                                  subcellOrd,
                                                                  parentCell);
    PHX::Device().fence();

    const auto subcellBaseTopo = shards::CellTopology(parentCell.getBaseCellTopologyData(subcellDim, subcellOrd));
    TEUCHOS_ASSERT(subcellBaseTopo.getBaseKey() == shards::Triangle<>::key ||
                   subcellBaseTopo.getBaseKey() == shards::Quadrilateral<>::key);

    const WorksetDetails & details = workset;
    const auto subcellVertexCount = static_cast<Intrepid2::ordinal_type>(subcellBaseTopo.getVertexCount());

    // Copy orientations to device.
    Kokkos::View<Intrepid2::Orientation*,PHX::Device> orientations_device(Kokkos::view_alloc("orientations_device",Kokkos::WithoutInitializing),orientations->size());
    auto orientations_host = Kokkos::create_mirror_view(Kokkos::WithoutInitializing,orientations_device);
    for (size_t i=0; i < orientations_host.extent(0); ++i)
      orientations_host(i) = orientations->at(i);
    Kokkos::deep_copy(orientations_device,orientations_host);

    // Local temporaries for device lambda capture
    const auto cell_local_ids_k = details.cell_local_ids_k;
    auto residual_local = residual;
    auto dof_local = dof;
    auto value_local = value;
    auto faceNormal_local = faceNormal;

    Kokkos::parallel_for("panzer::DirichletRsidual_FaceBasis::evalauteFields",
                         workset.num_cells,
                         KOKKOS_LAMBDA(const index_t c) {
      const auto ort = orientations_device(cell_local_ids_k[c]);
      Intrepid2::ordinal_type faceOrts[6] = {};
      ort.getFaceOrientation(faceOrts, numFaces); 

      // vertex count represent rotation count before it flips
      const double ortVal = faceOrts[subcellOrd] < subcellVertexCount ? 1.0 : -1.0;
      for(int b=0;b<numFaceDofs;b++) {
        residual_local(c,b) = 0.0;
        for(int d=0;d<cellDim;d++)
          residual_local(c,b) += (dof_local(c,b,d)-value_local(c,b,d))*faceNormal_local(c,b,d);
        residual_local(c,b) *= ortVal;
      }
    });
  }
}

//**********************************************************************

}

#endif
