// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_RESPONSE_SCATTER_EVALUATOR_EXTREMEVALUE_IMPL_HPP
#define PANZER_RESPONSE_SCATTER_EVALUATOR_EXTREMEVALUE_IMPL_HPP

#include <iostream>
#include <string>

#include "PanzerDiscFE_config.hpp"

#include "Phalanx_Evaluator_Macros.hpp"
#include "Phalanx_MDField.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"

#include "Panzer_CellData.hpp"
#include "Panzer_PointRule.hpp"
#include "Panzer_CommonArrayFactories.hpp"
#include "Panzer_ResponseBase.hpp"
#include "Panzer_Dimension.hpp"
#include "Panzer_GlobalEvaluationDataContainer.hpp"

#include "Intrepid2_FunctionSpaceTools.hpp"

#include "Thyra_SpmdVectorBase.hpp"
#include "Teuchos_ArrayRCP.hpp"

#include "Kokkos_ViewFactory.hpp"

namespace panzer {

template<typename EvalT, typename Traits, typename LO, typename GO>
ResponseScatterEvaluator_ProbeBase<EvalT,Traits,LO,GO>::
ResponseScatterEvaluator_ProbeBase(
  const std::string & responseName,
  const std::string & fieldName,
  const int fieldComponent,
  const Teuchos::Array<double>& point,
  const IntegrationRule & ir,
  const Teuchos::RCP<const PureBasis>& basis,
  const Teuchos::RCP<const panzer::GlobalIndexer> & indexer,
  const Teuchos::RCP<ProbeScatterBase> & probeScatter)
  : responseName_(responseName)
  , fieldName_(fieldName)
  , fieldComponent_(fieldComponent)
  , point_(point)
  , basis_(basis)
  , topology_(ir.topology)
  , globalIndexer_(indexer)
  , scatterObj_(probeScatter)
  , haveProbe_(false)
  , cellIndex_(-1)
  , workset_id_(0)
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  // the field manager will allocate all of these fields
  field_ = PHX::MDField<const ScalarT,Cell,BASIS>(fieldName,basis_->functional);
  this->addDependentField(field_);

  num_basis = basis->cardinality();
  num_dim = basis->dimension();
  TEUCHOS_ASSERT(num_dim == static_cast<size_t>(point_.size()));

  basis_values_ = Kokkos::DynRankView<double,PHX::Device>(
    "basis_values", 1, num_basis, 1); // Cell, Basis, Point

  // build dummy target tag
  std::string dummyName =
    ResponseBase::buildLookupName(responseName) + " dummy target";
  RCP<PHX::DataLayout> dl_dummy = rcp(new PHX::MDALayout<panzer::Dummy>(0));
  scatterHolder_ = rcp(new PHX::Tag<ScalarT>(dummyName,dl_dummy));
  this->addEvaluatedField(*scatterHolder_);

  std::string n = "Probe Response Scatter: " + responseName;
  this->setName(n);
}

template<typename EvalT, typename Traits, typename LO, typename GO>
void ResponseScatterEvaluator_ProbeBase<EvalT,Traits,LO,GO>::
postRegistrationSetup(typename Traits::SetupData sd,
                      PHX::FieldManager<Traits>& )
{
  for (const auto& workset : *sd.worksets_) {
    this->findCellAndComputeBasisValues(workset);
    if (haveProbe_)
      break;
  }
}

template<typename EvalT, typename Traits, typename LO, typename GO>
void ResponseScatterEvaluator_ProbeBase<EvalT,Traits,LO,GO>::
preEvaluate(typename Traits::PreEvalData d)
{
  // extract linear object container
  responseObj_ =
    Teuchos::rcp_dynamic_cast<Response_Probe<EvalT> >(
      d.gedc->getDataObject(ResponseBase::buildLookupName(responseName_)),
      true);
}

template<typename EvalT, typename Traits, typename LO, typename GO>
bool ResponseScatterEvaluator_ProbeBase<EvalT,Traits,LO,GO>::
findCellAndComputeBasisValues(typename Traits::EvalData d)
{
  // This evaluator needs to run on host until checkPointwiseInclusion
  // is moved to device.
  using HostSpace = Kokkos::DefaultHostExecutionSpace;
  using CTD = Intrepid2::CellTools<HostSpace>;
  using FST = Intrepid2::FunctionSpaceTools<HostSpace>;

  // Find which cell contains our point
  const int num_points = 1;
  Kokkos::DynRankView<int,HostSpace> inCell("inCell", this->wda(d).cell_node_coordinates.extent_int(0), num_points);
  Kokkos::DynRankView<double,HostSpace> physical_points_cell("physical_points_cell", this->wda(d).cell_node_coordinates.extent_int(0), num_points, num_dim);
  auto tmp_point = point_;
  {
    Kokkos::MDRangePolicy<HostSpace,Kokkos::Rank<2>> policy({0,0},{d.num_cells,static_cast<decltype(d.num_cells)>(num_dim)});
    Kokkos::parallel_for("copy node coords",policy,[&](const int cell, const int dim){
      physical_points_cell(cell,0,dim) = tmp_point[dim];
    });
    HostSpace().fence();

    auto cell_coords = this->wda(d).cell_node_coordinates.get_view();
    auto cell_coords_host = Kokkos::create_mirror_view(cell_coords);
    Kokkos::deep_copy(cell_coords_host, cell_coords);

    const double tol = 1.0e-12;
    CTD::checkPointwiseInclusion(inCell,
                                 physical_points_cell,
                                 cell_coords_host,
                                 *topology_,
                                 tol);
  }

  for (index_t cell=0; cell<static_cast<int>(d.num_cells); ++cell) {
    if (inCell(cell,0) == 1) {
      cellIndex_ = cell;
      workset_id_ = d.getIdentifier();
      haveProbe_ = true;
      break;
    }
  }

  // If no cell does, we're done
  if (!haveProbe_) {
    return false;
  }

  // Map point to reference frame
  const size_t num_nodes = this->wda(d).cell_node_coordinates.extent(1);
  Kokkos::DynRankView<double,HostSpace> cell_coords("cell_coords", 1, int(num_nodes), int(num_dim)); // <C,B,D>
  auto cnc_host = Kokkos::create_mirror_view(this->wda(d).cell_node_coordinates.get_view());
  Kokkos::deep_copy(cnc_host,this->wda(d).cell_node_coordinates.get_view());
  for (size_t i=0; i<num_nodes; ++i) {
    for (size_t j=0; j<num_dim; ++j) {
      cell_coords(0,i,j) = cnc_host(cellIndex_,i,j);
    }
  }
  Kokkos::DynRankView<double,HostSpace> physical_points("physical_points", 1, 1, num_dim); // <C,P,D>
  for (size_t i=0; i<num_dim; ++i)
    physical_points(0,0,i) = physical_points_cell(0,0,i);

  Kokkos::DynRankView<double,HostSpace> reference_points("reference_points", 1, 1, num_dim); // <C,P,D>
  CTD::mapToReferenceFrame(reference_points, physical_points, cell_coords, *topology_);

  Kokkos::DynRankView<double,HostSpace> reference_points_cell("reference_points_cell", 1, num_dim); // <P,D>
  for (size_t i=0; i<num_dim; ++i)
    reference_points_cell(0,i) = reference_points(0,0,i);

  // Compute basis functions at point
  if (basis_->getElementSpace() == PureBasis::CONST ||
      basis_->getElementSpace() == PureBasis::HGRAD) {

    // Evaluate basis at reference values
    Kokkos::DynRankView<double,HostSpace> ref_basis_values("ref_basis_values", num_basis, 1); // <B,P>
    basis_->getIntrepid2Basis<HostSpace,double,double>()->getValues(ref_basis_values,
                                                                    reference_points_cell,
                                                                    Intrepid2::OPERATOR_VALUE);

    // Apply transformation to physical frame
    auto basis_values_host = Kokkos::create_mirror_view(basis_values_);
    FST::HGRADtransformVALUE<double>(basis_values_host, ref_basis_values);
    Kokkos::deep_copy(basis_values_,basis_values_host);
  }
  else if (basis_->getElementSpace() == PureBasis::HCURL ||
           basis_->getElementSpace() == PureBasis::HDIV) {

    // Evaluate basis at reference values
    Kokkos::DynRankView<double,HostSpace> ref_basis_values("ref_basis_values", num_basis, 1, num_dim); // <B,P,D>
    basis_->getIntrepid2Basis<HostSpace,double,double>()->getValues(ref_basis_values,
                                                                    reference_points_cell,
                                                                    Intrepid2::OPERATOR_VALUE);

    // Apply transformation to physical frame
    Kokkos::DynRankView<double,HostSpace> jac("jac", 1, 1, num_dim, num_dim); // <C,P,D,D>
    CTD::setJacobian(jac, reference_points, cell_coords, *topology_);
    Kokkos::DynRankView<double,HostSpace> basis_values_vec("basis_values_vec", 1, num_basis, 1, num_dim); // <C,B,P,D>
    if (basis_->getElementSpace() == PureBasis::HCURL) {
      Kokkos::DynRankView<double,HostSpace> jac_inv("jac_inv", 1, 1, num_dim, num_dim); // <C,P,D,D>
      CTD::setJacobianInv(jac_inv, jac);
      FST::HCURLtransformVALUE<double>(basis_values_vec, jac_inv,
                                       ref_basis_values);
    }
    else {
      Kokkos::DynRankView<double,HostSpace> jac_det("jac_det", 1, 1); // <C,P>
      CTD::setJacobianDet(jac_det, jac);
      FST::HDIVtransformVALUE<double>(basis_values_vec, jac, jac_det,
                                      ref_basis_values);
    }

    // Compute element orientations
    std::vector<double> orientation;
    globalIndexer_->getElementOrientation(cellIndex_, orientation);
    std::string blockId = this->wda(d).block_id;
    int fieldNum = globalIndexer_->getFieldNum(fieldName_);
    const std::vector<int> & elmtOffset = globalIndexer_->getGIDFieldOffsets(blockId,fieldNum);

    // Extract component of basis
    for (size_t i=0; i<num_basis; ++i) {
      int offset = elmtOffset[i];
      basis_values_(0,i,0) = orientation[offset] * basis_values_vec(0,i,0,fieldComponent_);
    }

  }

  return true;
}

template<typename EvalT, typename Traits, typename LO, typename GO>
void ResponseScatterEvaluator_ProbeBase<EvalT,Traits,LO,GO>::
evaluateFields(typename Traits::EvalData d)
{
  using HostSpace = Kokkos::DefaultHostExecutionSpace;

  if ( !haveProbe_ ||
       (haveProbe_ && d.getIdentifier() != workset_id_) )
    return;

  auto field_coeffs_host = Kokkos::create_mirror_view(field_.get_view());
  Kokkos::deep_copy(field_coeffs_host,field_.get_view());

  auto field_coeffs_host_subview = Kokkos::subview(field_coeffs_host,std::pair<int,int>(cellIndex_,cellIndex_+1),Kokkos::ALL);

  auto field_val = Kokkos::createDynRankViewWithType<Kokkos::DynRankView<ScalarT,HostSpace>>(field_coeffs_host, "field_val_at_point", 1, 1); // <C,P>

  auto basis_values_host = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),basis_values_);

  Intrepid2::FunctionSpaceTools<HostSpace>::evaluate(field_val, field_coeffs_host_subview, basis_values_host);
  responseObj_->value = field_val(0,0);
  responseObj_->have_probe = true;
}

template <typename LO, typename GO>
void ResponseScatterEvaluator_Probe<panzer::Traits::Jacobian,panzer::Traits,LO,GO>::
evaluateFields(panzer::Traits::EvalData d)
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Thyra::SpmdVectorBase;

  TEUCHOS_ASSERT(this->scatterObj_!=Teuchos::null);

  Base::evaluateFields(d);

  // grab local data for inputing
  Teuchos::ArrayRCP<double> local_dgdx;
  RCP<SpmdVectorBase<double> > dgdx =
    rcp_dynamic_cast<SpmdVectorBase<double> >(this->responseObj_->getGhostedVector());
  dgdx->getNonconstLocalData(ptrFromRef(local_dgdx));
  TEUCHOS_ASSERT(!local_dgdx.is_null());

  this->scatterObj_->scatterDerivative(this->responseObj_->value,
                                       this->cellIndex_,
                                       this->responseObj_->have_probe,
                                       d,this->wda,
                                       local_dgdx);
}

}

#endif
