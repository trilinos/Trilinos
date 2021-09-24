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

#include "PanzerDiscFE_config.hpp"
#include "Panzer_Traits.hpp"

#include "Panzer_CommonArrayFactories.hpp"
#include "Kokkos_ViewFactory.hpp"
#include "Panzer_OrientationsInterface.hpp"

#include "Intrepid2_Utils.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_Orientation.hpp"
#include "Intrepid2_OrientationTools.hpp"

// FIXME: There are some calls in Intrepid2 that require non-const arrays when they should be const - search for PHX::getNonConstDynRankViewFromConstMDField
#include "Phalanx_GetNonConstDynRankViewFromConstMDField.hpp"

namespace panzer {
namespace {

template<typename Scalar>
void
applyOrientationsImpl(const int num_cells,
                      Kokkos::DynRankView<Scalar, PHX::Device> view,
                      Kokkos::DynRankView<Intrepid2::Orientation, PHX::Device> orientations,
                      const typename BasisValues2<Scalar>::IntrepidBasis & basis)
{
  using ots=Intrepid2::OrientationTools<PHX::Device>;

  auto sub_orientations = Kokkos::subview(orientations, std::make_pair(0,num_cells));

  // There are two options based on rank
  if(view.rank() == 3){
    // Grab subview of object to re-orient and create a copy of it
    auto sub_view = Kokkos::subview(view, std::make_pair(0,num_cells), Kokkos::ALL(), Kokkos::ALL());
    auto sub_view_clone = Kokkos::createDynRankView(view, "sub_view_clone", sub_view.extent(0), sub_view.extent(1), sub_view.extent(2));
    Kokkos::deep_copy(sub_view_clone, sub_view);

    // Apply the orientations to the subview
    ots::modifyBasisByOrientation(sub_view, sub_view_clone, sub_orientations, &basis);
  } else if (view.rank() == 4){
    // Grab subview of object to re-orient and create a copy of it
    auto sub_view = Kokkos::subview(view, std::make_pair(0,num_cells), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
    auto sub_view_clone = Kokkos::createDynRankView(view, "sub_view_clone", sub_view.extent(0), sub_view.extent(1), sub_view.extent(2), sub_view.extent(3));
    Kokkos::deep_copy(sub_view_clone, sub_view);

    // Apply the orientations to the subview
    ots::modifyBasisByOrientation(sub_view, sub_view_clone, sub_orientations, &basis);
  } else
    throw std::logic_error("applyOrientationsImpl : Unknown view of rank " + std::to_string(view.rank()));
}

template<typename Scalar>
void
applyOrientationsImpl(const int num_cells,
                      Kokkos::DynRankView<Scalar, PHX::Device> view,
                      const std::vector<Intrepid2::Orientation> & orientations,
                      const typename BasisValues2<Scalar>::IntrepidBasis & basis)
{

  // Move orientations vector to device
  Kokkos::DynRankView<Intrepid2::Orientation,PHX::Device> device_orientations("drv_orts", num_cells);
  auto host_orientations = Kokkos::create_mirror_view(device_orientations);
  for(int i=0; i < num_cells; ++i)
    host_orientations(i) = orientations[i];
  Kokkos::deep_copy(device_orientations,host_orientations);

  // Call the device orientations applicator
  applyOrientationsImpl(num_cells, view, device_orientations, basis);
}

}


template <typename Scalar>
panzer::BasisValues2<Scalar>::
BasisValues2(const std::string & pre,
             const bool allocArrays,
             const bool buildWeighted)
    : compute_derivatives(true)
    , build_weighted(buildWeighted)
    , alloc_arrays(allocArrays)
    , prefix(pre)
    , ddims_(1,0)
    , references_evaluated(false)
    , orientations_applied_(false)
    , num_cells_(0)
    , num_evaluate_cells_(0)
    , is_uniform_(false)

{
  // Default all lazy evaluated components to not-evaluated
  basis_ref_scalar_evaluated_ = false;
  basis_scalar_evaluated_ = false;
  basis_ref_vector_evaluated_ = false;
  basis_vector_evaluated_ = false;
  grad_basis_ref_evaluated_ = false;
  grad_basis_evaluated_ = false;
  curl_basis_ref_scalar_evaluated_ = false;
  curl_basis_scalar_evaluated_ = false;
  curl_basis_ref_vector_evaluated_ = false;
  curl_basis_vector_evaluated_ = false;
  div_basis_ref_evaluated_ = false;
  div_basis_evaluated_ = false;
  weighted_basis_scalar_evaluated_ = false;
  weighted_basis_vector_evaluated_ = false;
  weighted_grad_basis_evaluated_ = false;
  weighted_curl_basis_scalar_evaluated_ = false;
  weighted_curl_basis_vector_evaluated_ = false;
  weighted_div_basis_evaluated_ = false;
  basis_coordinates_ref_evaluated_ = false;
  basis_coordinates_evaluated_ = false;
}

template <typename Scalar>
void panzer::BasisValues2<Scalar>::
evaluateValues(const PHX::MDField<Scalar,IP,Dim> & cub_points,
               const PHX::MDField<Scalar,Cell,IP,Dim,Dim> & jac,
               const PHX::MDField<Scalar,Cell,IP> & jac_det,
               const PHX::MDField<Scalar,Cell,IP,Dim,Dim> & jac_inv,
               const int in_num_cells)
{

  build_weighted = false;

  setupUniform(basis_layout, cub_points, jac, jac_det, jac_inv, in_num_cells);

  const auto elmtspace = getElementSpace();
  const int num_dims = jac.extent(2);

  // Evaluate Reference values
  if(elmtspace == PureBasis::HGRAD or elmtspace == PureBasis::CONST or elmtspace == PureBasis::HVOL)
    getBasisValuesRef(true,true);

  if(elmtspace == PureBasis::HDIV or elmtspace == PureBasis::HCURL)
    getVectorBasisValuesRef(true,true);

  if(elmtspace == PureBasis::HGRAD and compute_derivatives)
    getGradBasisValuesRef(true,true);

  if(elmtspace == PureBasis::HCURL and compute_derivatives){
    if(num_dims == 2)
      getCurl2DVectorBasisRef(true,true);
    else if(num_dims == 3)
      getCurlVectorBasisRef(true,true);
  }

  if(elmtspace == PureBasis::HDIV and compute_derivatives)
    getDivVectorBasisRef(true, true);

  references_evaluated = true;

  // Evaluate real space values
  if(elmtspace == PureBasis::HGRAD or elmtspace == PureBasis::CONST or elmtspace == PureBasis::HVOL)
    getBasisValues(false,true,true);

  if(elmtspace == PureBasis::HDIV or elmtspace == PureBasis::HCURL)
    getVectorBasisValues(false,true,true);

  if(elmtspace == PureBasis::HGRAD and compute_derivatives)
    getGradBasisValues(false,true,true);

  if(elmtspace == PureBasis::HCURL and compute_derivatives){
    if(num_dims == 2)
      getCurl2DVectorBasis(false,true,true);
    else if(num_dims == 3)
      getCurlVectorBasis(false,true,true);
  }

  if(elmtspace == PureBasis::HDIV and compute_derivatives)
    getDivVectorBasis(false,true,true);

  // Orientations have been applied only if this pointer is valid
  orientations_applied_ = orientations_ != Teuchos::null;
}

template <typename Scalar>
void panzer::BasisValues2<Scalar>::
evaluateValues(const PHX::MDField<Scalar,IP,Dim> & cub_points,
               const PHX::MDField<Scalar,Cell,IP,Dim,Dim> & jac,
               const PHX::MDField<Scalar,Cell,IP> & jac_det,
               const PHX::MDField<Scalar,Cell,IP,Dim,Dim> & jac_inv,
               const PHX::MDField<Scalar,Cell,IP> & weighted_measure,
               const PHX::MDField<Scalar,Cell,NODE,Dim> & vertex_coordinates,
               bool use_vertex_coordinates,
               const int in_num_cells)
{

  // This is the base call that will add all the non-weighted versions
  evaluateValues(cub_points, jac, jac_det, jac_inv, in_num_cells);
  if(weighted_measure.size() > 0)
    setWeightedMeasure(weighted_measure);

  cell_vertex_coordinates_ = vertex_coordinates;

  // Add the weighted versions of all basis components - this will add to the non-weighted versions generated by evaluateValues above

  const auto elmtspace = getElementSpace();
  const int num_dims = jac.extent(2);

  if(elmtspace == PureBasis::HGRAD or elmtspace == PureBasis::CONST or elmtspace == PureBasis::HVOL)
    getBasisValues(true,true,true);

  if(elmtspace == PureBasis::HDIV or elmtspace == PureBasis::HCURL)
    getVectorBasisValues(true,true,true);

  if(elmtspace == PureBasis::HGRAD and compute_derivatives)
    getGradBasisValues(true,true,true);

  if(elmtspace == PureBasis::HCURL and compute_derivatives){
    if(num_dims == 2)
      getCurl2DVectorBasis(true,true,true);
    else if(num_dims == 3)
      getCurlVectorBasis(true,true,true);
  }

  if(elmtspace == PureBasis::HDIV and compute_derivatives)
    getDivVectorBasis(true,true,true);

  // Add the vertex components
  if(use_vertex_coordinates){
    getBasisCoordinatesRef(true,true);
    getBasisCoordinates(true,true);
  }

  // Orientations have been applied only if this pointer is valid
  orientations_applied_ = orientations_ != Teuchos::null;
}

template <typename Scalar>
void panzer::BasisValues2<Scalar>::
evaluateValues(const PHX::MDField<Scalar,Cell,IP,Dim> & cub_points,
               const PHX::MDField<Scalar,Cell,IP,Dim,Dim> & jac,
               const PHX::MDField<Scalar,Cell,IP> & jac_det,
               const PHX::MDField<Scalar,Cell,IP,Dim,Dim> & jac_inv,
               const PHX::MDField<Scalar,Cell,IP> & weighted_measure,
               const PHX::MDField<Scalar,Cell,NODE,Dim> & vertex_coordinates,
               bool use_vertex_coordinates,
               const int in_num_cells)
{

  cell_vertex_coordinates_ = vertex_coordinates;

  setup(basis_layout, cub_points, jac, jac_det, jac_inv, in_num_cells);
  if(weighted_measure.size() > 0)
    setWeightedMeasure(weighted_measure);

  const auto elmtspace = getElementSpace();
  const int num_dims = jac.extent(2);

  if(elmtspace == PureBasis::HGRAD or elmtspace == PureBasis::CONST or elmtspace == PureBasis::HVOL){
    getBasisValues(false,true,true);
    if(build_weighted) getBasisValues(true,true,true);
  }

  if(elmtspace == PureBasis::HDIV or elmtspace == PureBasis::HCURL){
    getVectorBasisValues(false,true,true);
    if(build_weighted) getVectorBasisValues(true,true,true);
  }

  if(elmtspace == PureBasis::HGRAD and compute_derivatives){
    getGradBasisValues(false,true,true);
    if(build_weighted) getGradBasisValues(true,true,true);
  }

  if(elmtspace == PureBasis::HCURL and compute_derivatives){
    if(num_dims == 2){
      getCurl2DVectorBasis(false,true,true);
      if(build_weighted) getCurl2DVectorBasis(true,true,true);
    } else if(num_dims == 3) {
      getCurlVectorBasis(false,true,true);
      if(build_weighted) getCurlVectorBasis(true,true,true);
    }
  }

  if(elmtspace == PureBasis::HDIV and compute_derivatives){
    getDivVectorBasis(false,true,true);
    if(build_weighted) getDivVectorBasis(true,true,true);
  }

  // Add the vertex components
  if(use_vertex_coordinates){
    getBasisCoordinatesRef(true,true);
    getBasisCoordinates(true,true);
  }

  // Orientations have been applied only if this pointer is valid
  orientations_applied_ = orientations_ != Teuchos::null;

}

template <typename Scalar>
void panzer::BasisValues2<Scalar>::
evaluateValuesCV(const PHX::MDField<Scalar,Cell,IP,Dim> & cub_points,
                 const PHX::MDField<Scalar,Cell,IP,Dim,Dim> & jac,
                 const PHX::MDField<Scalar,Cell,IP> & jac_det,
                 const PHX::MDField<Scalar,Cell,IP,Dim,Dim> & jac_inv)
{

  PHX::MDField<Scalar,Cell,IP> weighted_measure;
  PHX::MDField<Scalar,Cell,NODE,Dim> vertex_coordinates;
  evaluateValues(cub_points,jac,jac_det,jac_inv,weighted_measure,vertex_coordinates,false,jac.extent(0));

}

template <typename Scalar>
void panzer::BasisValues2<Scalar>::
evaluateValuesCV(const PHX::MDField<Scalar,Cell,IP,Dim> & cell_cub_points,
                 const PHX::MDField<Scalar,Cell,IP,Dim,Dim> & jac,
                 const PHX::MDField<Scalar,Cell,IP> & jac_det,
                 const PHX::MDField<Scalar,Cell,IP,Dim,Dim> & jac_inv,
                 const PHX::MDField<Scalar,Cell,NODE,Dim> & vertex_coordinates,
                 bool use_vertex_coordinates,
                 const int in_num_cells)
{
  PHX::MDField<Scalar,Cell,IP> weighted_measure;
  evaluateValues(cell_cub_points,jac,jac_det,jac_inv,weighted_measure,vertex_coordinates,use_vertex_coordinates,in_num_cells);
}

template <typename Scalar>
void panzer::BasisValues2<Scalar>::
evaluateBasisCoordinates(const PHX::MDField<Scalar,Cell,NODE,Dim> & vertex_coordinates,
                         const int in_num_cells)
{
  num_evaluate_cells_ = in_num_cells < 0 ? vertex_coordinates.extent(0) : in_num_cells;
  cell_vertex_coordinates_ = vertex_coordinates;

  getBasisCoordinates(true,true);
}

// method for applying orientations
template <typename Scalar>
void BasisValues2<Scalar>::
applyOrientations(const std::vector<Intrepid2::Orientation> & orientations,
                  const int in_num_cells)
{
  if (!intrepid_basis->requireOrientation())
    return;

  // We only allow the orientations to be applied once
  TEUCHOS_ASSERT(not orientations_applied_);

  const PureBasis::EElementSpace elmtspace = getElementSpace();

  const int num_cell_basis_layout = in_num_cells < 0 ? basis_layout->numCells() : in_num_cells;
  const int num_cell_orientation = orientations.size();
  const int num_cells  = num_cell_basis_layout < num_cell_orientation ? num_cell_basis_layout : num_cell_orientation;
  const int num_dim   = basis_layout->dimension();

  // Copy orientations to device
  Kokkos::DynRankView<Intrepid2::Orientation,PHX::Device> device_orientations("device_orientations", num_cells);
  auto host_orientations = Kokkos::create_mirror_view(device_orientations);
  for(int i=0; i < num_cells; ++i)
    host_orientations(i) = orientations[i];
  Kokkos::deep_copy(device_orientations,host_orientations);

  if (elmtspace == PureBasis::HGRAD){
    applyOrientationsImpl<Scalar>(num_cells, basis_scalar.get_view(), device_orientations, *intrepid_basis);
    if(build_weighted) applyOrientationsImpl<Scalar>(num_cells, weighted_basis_scalar.get_view(), device_orientations, *intrepid_basis);
    if(compute_derivatives) applyOrientationsImpl<Scalar>(num_cells, grad_basis.get_view(), device_orientations, *intrepid_basis);
    if(compute_derivatives and build_weighted) applyOrientationsImpl<Scalar>(num_cells, weighted_grad_basis.get_view(), device_orientations, *intrepid_basis);
  }

  if(elmtspace == PureBasis::HCURL){
    applyOrientationsImpl<Scalar>(num_cells, basis_vector.get_view(), device_orientations, *intrepid_basis);
    if(build_weighted) applyOrientationsImpl<Scalar>(num_cells, weighted_basis_vector.get_view(), device_orientations, *intrepid_basis);
    if(num_dim == 2){
      if(compute_derivatives) applyOrientationsImpl<Scalar>(num_cells, curl_basis_scalar.get_view(), device_orientations, *intrepid_basis);
      if(compute_derivatives and build_weighted) applyOrientationsImpl<Scalar>(num_cells, weighted_curl_basis_scalar.get_view(), device_orientations, *intrepid_basis);
    }
    if(num_dim == 3){
      if(compute_derivatives) applyOrientationsImpl<Scalar>(num_cells, curl_basis_vector.get_view(), device_orientations, *intrepid_basis);
      if(compute_derivatives and build_weighted) applyOrientationsImpl<Scalar>(num_cells, weighted_curl_basis_vector.get_view(), device_orientations, *intrepid_basis);
    }
  }

  if(elmtspace == PureBasis::HDIV){
    applyOrientationsImpl<Scalar>(num_cells, basis_vector.get_view(), device_orientations, *intrepid_basis);
    if(build_weighted) applyOrientationsImpl<Scalar>(num_cells, weighted_basis_vector.get_view(), device_orientations, *intrepid_basis);
    if(compute_derivatives) applyOrientationsImpl<Scalar>(num_cells, div_basis.get_view(), device_orientations, *intrepid_basis);
    if(compute_derivatives and build_weighted) applyOrientationsImpl<Scalar>(num_cells, weighted_div_basis.get_view(), device_orientations, *intrepid_basis);
  }

  orientations_applied_ = true;
}

// method for applying orientations
template <typename Scalar>
void BasisValues2<Scalar>::
applyOrientations(const PHX::MDField<const Scalar,Cell,BASIS> & orientations)
{
  TEUCHOS_TEST_FOR_EXCEPT_MSG(true,"panzer::BasisValues2::applyOrientations : this should not be called.");
}

template <typename Scalar>
PureBasis::EElementSpace BasisValues2<Scalar>::getElementSpace() const
{ return basis_layout->getBasis()->getElementSpace(); }

template <typename Scalar>
void panzer::BasisValues2<Scalar>::
setupArrays(const Teuchos::RCP<const panzer::BasisIRLayout>& layout,
            bool computeDerivatives)
{
  MDFieldArrayFactory af(prefix,alloc_arrays);

  compute_derivatives = computeDerivatives;
  basis_layout = layout;
  num_cells_ = basis_layout->numCells();
  Teuchos::RCP<const panzer::PureBasis> basisDesc = layout->getBasis();

  // for convience pull out basis and quadrature information
  int num_quad = layout->numPoints();
  int dim      = basisDesc->dimension();
  int card     = basisDesc->cardinality();
  int numcells = basisDesc->numCells();
  panzer::PureBasis::EElementSpace elmtspace = basisDesc->getElementSpace();
  Teuchos::RCP<const shards::CellTopology> cellTopo = basisDesc->getCellTopology();

  intrepid_basis = basisDesc->getIntrepid2Basis<PHX::Device::execution_space,Scalar,Scalar>();

  // allocate field containers
  // field sizes defined by http://trilinos.sandia.gov/packages/docs/dev/packages/intrepid/doc/html/basis_page.html#basis_md_array_sec

  // compute basis fields
  if(elmtspace==panzer::PureBasis::HGRAD) {
     // HGRAD is a nodal field

     // build values
     ///////////////////////////////////////////////
     basis_ref_scalar = af.buildStaticArray<Scalar,BASIS,IP>("basis_ref",card,num_quad); // F, P
     basis_scalar = af.buildStaticArray<Scalar,Cell,BASIS,IP>("basis",numcells,card,num_quad);

     if(build_weighted)
       weighted_basis_scalar = af.buildStaticArray<Scalar,Cell,BASIS,IP>("weighted_basis",numcells,card,num_quad);

     // build gradients
     ///////////////////////////////////////////////

     if(compute_derivatives) {
       grad_basis_ref = af.buildStaticArray<Scalar,BASIS,IP,Dim>("grad_basis_ref",card,num_quad,dim); // F, P, D
       grad_basis = af.buildStaticArray<Scalar,Cell,BASIS,IP,Dim>("grad_basis",numcells,card,num_quad,dim);

       if(build_weighted)
         weighted_grad_basis = af.buildStaticArray<Scalar,Cell,BASIS,IP,Dim>("weighted_grad_basis",numcells,card,num_quad,dim);
     }

     // build curl
     ///////////////////////////////////////////////

     // nothing - HGRAD does not support CURL operation
  }
  else if(elmtspace==panzer::PureBasis::HCURL) {
     // HCURL is a vector field

     // build values
     ///////////////////////////////////////////////

     basis_ref_vector = af.buildStaticArray<Scalar,BASIS,IP,Dim>("basis_ref",card,num_quad,dim); // F, P, D
     basis_vector = af.buildStaticArray<Scalar,Cell,BASIS,IP,Dim>("basis",numcells,card,num_quad,dim);

     if(build_weighted)
       weighted_basis_vector = af.buildStaticArray<Scalar,Cell,BASIS,IP,Dim>("weighted_basis",numcells,card,num_quad,dim);

     // build gradients
     ///////////////////////////////////////////////

     // nothing - HCURL does not support GRAD operation

     // build curl
     ///////////////////////////////////////////////

     if(compute_derivatives) {
       if(dim==2) {
          // curl of HCURL basis is not dimension dependent
          curl_basis_ref_scalar = af.buildStaticArray<Scalar,BASIS,IP>("curl_basis_ref",card,num_quad); // F, P
          curl_basis_scalar = af.buildStaticArray<Scalar,Cell,BASIS,IP>("curl_basis",numcells,card,num_quad);

          if(build_weighted)
            weighted_curl_basis_scalar = af.buildStaticArray<Scalar,Cell,BASIS,IP>("weighted_curl_basis",numcells,card,num_quad);
       }
       else if(dim==3){
          curl_basis_ref_vector = af.buildStaticArray<Scalar,BASIS,IP,Dim>("curl_basis_ref",card,num_quad,dim); // F, P, D
          curl_basis_vector = af.buildStaticArray<Scalar,Cell,BASIS,IP,Dim>("curl_basis",numcells,card,num_quad,dim);

          if(build_weighted)
            weighted_curl_basis_vector = af.buildStaticArray<Scalar,Cell,BASIS,IP,Dim>("weighted_curl_basis",numcells,card,num_quad,dim);
       }
       else { TEUCHOS_ASSERT(false); } // what do I do with 1D?
     }
  }
  else if(elmtspace==panzer::PureBasis::HDIV) {
     // HDIV is a vector field

     // build values
     ///////////////////////////////////////////////

     basis_ref_vector = af.buildStaticArray<Scalar,BASIS,IP,Dim>("basis_ref",card,num_quad,dim); // F, P, D
     basis_vector = af.buildStaticArray<Scalar,Cell,BASIS,IP,Dim>("basis",numcells,card,num_quad,dim);

     if(build_weighted)
       weighted_basis_vector = af.buildStaticArray<Scalar,Cell,BASIS,IP,Dim>("weighted_basis",numcells,card,num_quad,dim);

     // build gradients
     ///////////////////////////////////////////////

     // nothing - HCURL does not support GRAD operation

     // build curl
     ///////////////////////////////////////////////

     // nothing - HDIV does not support CURL operation

     // build div
     ///////////////////////////////////////////////

     if(compute_derivatives) {
       div_basis_ref = af.buildStaticArray<Scalar,BASIS,IP>("div_basis_ref",card,num_quad); // F, P
       div_basis = af.buildStaticArray<Scalar,Cell,BASIS,IP>("div_basis",numcells,card,num_quad);

       if(build_weighted)
         weighted_div_basis = af.buildStaticArray<Scalar,Cell,BASIS,IP>("weighted_div_basis",numcells,card,num_quad);
     }
  }
  else if(elmtspace==panzer::PureBasis::CONST || elmtspace==panzer::PureBasis::HVOL) {
     // CONST is a nodal field

     // build values
     ///////////////////////////////////////////////
     basis_ref_scalar = af.buildStaticArray<Scalar,BASIS,IP>("basis_ref",card,num_quad); // F, P
     basis_scalar = af.buildStaticArray<Scalar,Cell,BASIS,IP>("basis",numcells,card,num_quad);

     if(build_weighted)
       weighted_basis_scalar = af.buildStaticArray<Scalar,Cell,BASIS,IP>("weighted_basis",numcells,card,num_quad);

     // build gradients
     ///////////////////////////////////////////////

     // nothing - CONST does not support GRAD operation

     // build curl
     ///////////////////////////////////////////////

     // nothing - CONST does not support CURL operation

     // build div
     ///////////////////////////////////////////////

     // nothing - CONST does not support DIV operation
  }
  else { TEUCHOS_ASSERT(false); }

  basis_coordinates_ref = af.buildStaticArray<Scalar,BASIS,Dim>("basis_coordinates_ref",card,dim);
  basis_coordinates = af.buildStaticArray<Scalar,Cell,BASIS,Dim>("basis_coordinates",numcells,card,dim);
}


//=======================================================================================================
// New Interface


template <typename Scalar>
void
BasisValues2<Scalar>::
setup(const Teuchos::RCP<const panzer::BasisIRLayout> & basis,
      PHX::MDField<const Scalar, Cell, IP, Dim>         reference_points,
      PHX::MDField<const Scalar, Cell, IP, Dim, Dim>    point_jacobian,
      PHX::MDField<const Scalar, Cell, IP>              point_jacobian_determinant,
      PHX::MDField<const Scalar, Cell, IP, Dim, Dim>    point_jacobian_inverse,
      const int                                         num_evaluated_cells)
{
  basis_layout = basis;
  intrepid_basis = basis->getBasis()->getIntrepid2Basis<PHX::Device::execution_space,Scalar,Scalar>();
  num_cells_ = basis_layout->numCells();
  num_evaluate_cells_ = num_evaluated_cells >= 0 ? num_evaluated_cells : num_cells_;
  build_weighted = false;
  is_uniform_ = false;

  cubature_points_ref_ = reference_points;
  cubature_jacobian_ = point_jacobian;
  cubature_jacobian_determinant_ = point_jacobian_determinant;
  cubature_jacobian_inverse_ = point_jacobian_inverse;

  // Reset internal data
  resetArrays();

}

template <typename Scalar>
void
BasisValues2<Scalar>::
setupUniform(const Teuchos::RCP<const panzer::BasisIRLayout> &  basis,
             PHX::MDField<const Scalar, IP, Dim>                reference_points,
             PHX::MDField<const Scalar, Cell, IP, Dim, Dim>     point_jacobian,
             PHX::MDField<const Scalar, Cell, IP>               point_jacobian_determinant,
             PHX::MDField<const Scalar, Cell, IP, Dim, Dim>     point_jacobian_inverse,
             const int                                          num_evaluated_cells)
{
  basis_layout = basis;
  intrepid_basis = basis->getBasis()->getIntrepid2Basis<PHX::Device::execution_space,Scalar,Scalar>();
  num_cells_ = basis_layout->numCells();
  num_evaluate_cells_ = num_evaluated_cells >= 0 ? num_evaluated_cells : num_cells_;
  cubature_points_uniform_ref_ = reference_points;
  build_weighted = false;
  is_uniform_ = true;

  cubature_jacobian_ = point_jacobian;
  cubature_jacobian_determinant_ = point_jacobian_determinant;
  cubature_jacobian_inverse_ = point_jacobian_inverse;

  // Reset internal data
  resetArrays();

}

template <typename Scalar>
void
BasisValues2<Scalar>::
setOrientations(const Teuchos::RCP<const OrientationsInterface> & orientations,
                const int num_orientations_cells)
{
  if(num_orientations_cells < 0)
    num_orientations_cells_ = num_evaluate_cells_;
  else
    num_orientations_cells_ = num_orientations_cells;
  if(orientations == Teuchos::null){
    orientations_applied_ = false;
    orientations_ = Teuchos::null;
    // Normally we would reset arrays here, but it seems to causes a lot of problems
  } else {
    orientations_ = orientations;
    orientations_applied_ = true;
    // Setting orientations means that we need to reset our arrays
    resetArrays();
  }
}

template <typename Scalar>
void
BasisValues2<Scalar>::
setCellVertexCoordinates(PHX::MDField<Scalar,Cell,NODE,Dim> vertex_coordinates)
{
  cell_vertex_coordinates_ = vertex_coordinates;
}

template <typename Scalar>
void
BasisValues2<Scalar>::
resetArrays()
{
  // Turn off all evaluated fields (forces re-evaluation)
  basis_ref_scalar_evaluated_ = false;
  basis_scalar_evaluated_ = false;
  basis_ref_vector_evaluated_ = false;
  basis_vector_evaluated_ = false;
  grad_basis_ref_evaluated_ = false;
  grad_basis_evaluated_ = false;
  curl_basis_ref_scalar_evaluated_ = false;
  curl_basis_scalar_evaluated_ = false;
  curl_basis_ref_vector_evaluated_ = false;
  curl_basis_vector_evaluated_ = false;
  div_basis_ref_evaluated_ = false;
  div_basis_evaluated_ = false;
  weighted_basis_scalar_evaluated_ = false;
  weighted_basis_vector_evaluated_ = false;
  weighted_grad_basis_evaluated_ = false;
  weighted_curl_basis_scalar_evaluated_ = false;
  weighted_curl_basis_vector_evaluated_ = false;
  weighted_div_basis_evaluated_ = false;
  basis_coordinates_ref_evaluated_ = false;
  basis_coordinates_evaluated_ = false;

  // TODO: Enable this feature - requires the old interface to go away
  // De-allocate arrays if necessary
//  if(not alloc_arrays){
//    basis_ref_scalar = Array_BasisIP();
//    basis_scalar = Array_CellBasisIP();
//    basis_ref_vector = Array_BasisIPDim();
//    basis_vector = Array_CellBasisIPDim();
//    grad_basis_ref = Array_BasisIPDim();
//    grad_basis = Array_CellBasisIPDim();
//    curl_basis_ref_scalar = Array_BasisIP();
//    curl_basis_scalar = Array_CellBasisIP();
//    curl_basis_ref_vector = Array_BasisIPDim();
//    curl_basis_vector = Array_CellBasisIPDim();
//    div_basis_ref = Array_BasisIP();
//    div_basis = Array_CellBasisIP();
//    weighted_basis_scalar = Array_CellBasisIP();
//    weighted_basis_vector = Array_CellBasisIPDim();
//    weighted_grad_basis = Array_CellBasisIPDim();
//    weighted_curl_basis_scalar = Array_CellBasisIP();
//    weighted_curl_basis_vector = Array_CellBasisIPDim();
//    weighted_div_basis = Array_CellBasisIP();
//    basis_coordinates_ref = Array_BasisDim();
//    basis_coordinates = Array_CellBasisDim();
//  }
}

template <typename Scalar>
void
BasisValues2<Scalar>::
setWeightedMeasure(PHX::MDField<const Scalar, Cell, IP> weighted_measure)
{
  TEUCHOS_TEST_FOR_EXCEPT_MSG(build_weighted,
                              "BasisValues2::setWeightedMeasure : Weighted measure already set. Can only set weighted measure once after setup or setupUniform have beens called.");
  cubature_weights_ = weighted_measure;
  build_weighted = true;
}

// If the array is allocated we can not reassign it - this means we
// have to deep copy into it. The use of deep_copy is need when an
// applicaiton or the library has cached a BasisValues2 member view
// outside of the BasisValues object. This could happen when the basis
// values object is embedded in an evalautor for mesh motion.
#define PANZER_CACHE_DATA(name) \
  if(cache) { \
    if(name.size()==tmp_##name.size()){ \
      Kokkos::deep_copy(name.get_view(), tmp_##name.get_view()); \
    } else { \
      name = tmp_##name; \
    } \
    name##_evaluated_ = true; \
  }

template <typename Scalar>
typename BasisValues2<Scalar>::ConstArray_BasisDim
BasisValues2<Scalar>::
getBasisCoordinatesRef(const bool cache,
                       const bool force) const
{
  // Check if array already exists
  if(basis_coordinates_ref_evaluated_ and not force)
    return basis_coordinates_ref;

  MDFieldArrayFactory af(prefix,getExtendedDimensions(),true);

  const int num_card  = basis_layout->cardinality();
  const int num_dim   = basis_layout->dimension();

  using coordsScalarType = typename Intrepid2::Basis<PHX::Device::execution_space,Scalar,Scalar>::scalarType;
  auto tmp_basis_coordinates_ref = af.buildStaticArray<coordsScalarType,BASIS,Dim>("basis_coordinates_ref", num_card, num_dim);
  intrepid_basis->getDofCoords(tmp_basis_coordinates_ref.get_view());
  PHX::Device().fence();

  // Store for later if cache is enabled
  PANZER_CACHE_DATA(basis_coordinates_ref)

  return tmp_basis_coordinates_ref;
}

template <typename Scalar>
typename BasisValues2<Scalar>::ConstArray_BasisIP
BasisValues2<Scalar>::
getBasisValuesRef(const bool cache,
                  const bool force) const
{
  // Check if array already exists
  if(basis_ref_scalar_evaluated_ and not force)
    return basis_ref_scalar;

  // Reference quantities only exist for a uniform reference space
  TEUCHOS_ASSERT(hasUniformReferenceSpace());

  // Make sure basis is valid
  PureBasis::EElementSpace elmtspace = getElementSpace();
  TEUCHOS_ASSERT(elmtspace==PureBasis::HGRAD || elmtspace==PureBasis::CONST || elmtspace==PureBasis::HVOL);

  MDFieldArrayFactory af(prefix,getExtendedDimensions(),true);

  const int num_quad  = basis_layout->numPoints();
  const int num_card  = basis_layout->cardinality();

  auto tmp_basis_ref_scalar = af.buildStaticArray<Scalar,BASIS,IP>("dyn_basis_ref_scalar",num_card,num_quad);
  auto cubature_points_uniform_ref = PHX::getNonConstDynRankViewFromConstMDField(cubature_points_uniform_ref_);

  intrepid_basis->getValues(tmp_basis_ref_scalar.get_view(), cubature_points_uniform_ref, Intrepid2::OPERATOR_VALUE);
  PHX::Device().fence();

  // Store for later if cache is enabled
  PANZER_CACHE_DATA(basis_ref_scalar);

  return tmp_basis_ref_scalar;
}

template <typename Scalar>
typename BasisValues2<Scalar>::ConstArray_BasisIPDim
BasisValues2<Scalar>::
getVectorBasisValuesRef(const bool cache,
                        const bool force) const
{
  // Check if array already exists
  if(basis_ref_vector_evaluated_ and not force)
    return basis_ref_vector;

  // Reference quantities only exist for a uniform reference space
  TEUCHOS_ASSERT(hasUniformReferenceSpace());

  // Make sure basis is valid
  PureBasis::EElementSpace elmtspace = getElementSpace();
  TEUCHOS_ASSERT(elmtspace==PureBasis::HDIV || elmtspace==PureBasis::HCURL);

  MDFieldArrayFactory af(prefix,getExtendedDimensions(),true);

  const int num_quad  = basis_layout->numPoints();
  const int num_card  = basis_layout->cardinality();
  const int num_dim   = basis_layout->dimension();

  auto tmp_basis_ref_vector = af.buildStaticArray<Scalar,BASIS,IP,Dim>("dyn_basis_ref_vector",num_card,num_quad,num_dim);
  auto cubature_points_uniform_ref = PHX::getNonConstDynRankViewFromConstMDField(cubature_points_uniform_ref_);

  intrepid_basis->getValues(tmp_basis_ref_vector.get_view(),cubature_points_uniform_ref,Intrepid2::OPERATOR_VALUE);
  PHX::Device().fence();

  // Store for later if cache is enabled
  PANZER_CACHE_DATA(basis_ref_vector);

  return tmp_basis_ref_vector;
}

template <typename Scalar>
typename BasisValues2<Scalar>::ConstArray_BasisIPDim
BasisValues2<Scalar>::
getGradBasisValuesRef(const bool cache,
                      const bool force) const
{
  // Check if array already exists
  if(grad_basis_ref_evaluated_ and not force)
    return grad_basis_ref;

  // Reference quantities only exist for a uniform reference space
  TEUCHOS_ASSERT(hasUniformReferenceSpace());

  // Make sure basis is valid
  PureBasis::EElementSpace elmtspace = getElementSpace();
  TEUCHOS_ASSERT(elmtspace==PureBasis::HGRAD);

  MDFieldArrayFactory af(prefix,getExtendedDimensions(),true);

  const int num_quad  = basis_layout->numPoints();
  const int num_card  = basis_layout->cardinality();
  const int num_dim   = basis_layout->dimension();

  auto tmp_grad_basis_ref = af.buildStaticArray<Scalar,BASIS,IP,Dim>("dyn_basis_ref_vector",num_card,num_quad,num_dim);
  auto cubature_points_uniform_ref = PHX::getNonConstDynRankViewFromConstMDField(cubature_points_uniform_ref_);

  intrepid_basis->getValues(tmp_grad_basis_ref.get_view(), cubature_points_uniform_ref, Intrepid2::OPERATOR_GRAD);
  PHX::Device().fence();

  // Store for later if cache is enabled
  PANZER_CACHE_DATA(grad_basis_ref);

  return tmp_grad_basis_ref;
}

template <typename Scalar>
typename BasisValues2<Scalar>::ConstArray_BasisIP
BasisValues2<Scalar>::
getCurl2DVectorBasisRef(const bool cache,
                        const bool force) const
{
  // Check if array already exists
  if(curl_basis_ref_scalar_evaluated_ and not force)
    return curl_basis_ref_scalar;

  // Reference quantities only exist for a uniform reference space
  TEUCHOS_ASSERT(hasUniformReferenceSpace());
  TEUCHOS_ASSERT(basis_layout->dimension() == 2);

  // Make sure basis is valid
  PureBasis::EElementSpace elmtspace = getElementSpace();
  TEUCHOS_ASSERT(elmtspace==PureBasis::HCURL);

  MDFieldArrayFactory af(prefix,getExtendedDimensions(),true);

  const int num_quad  = basis_layout->numPoints();
  const int num_card  = basis_layout->cardinality();

  auto tmp_curl_basis_ref_scalar = af.buildStaticArray<Scalar,BASIS,IP>("dyn_curl_basis_ref_scalar",num_card,num_quad);
  auto cubature_points_uniform_ref = PHX::getNonConstDynRankViewFromConstMDField(cubature_points_uniform_ref_);

  intrepid_basis->getValues(tmp_curl_basis_ref_scalar.get_view(), cubature_points_uniform_ref, Intrepid2::OPERATOR_CURL);
  PHX::Device().fence();

  // Store for later if cache is enabled
  PANZER_CACHE_DATA(curl_basis_ref_scalar);

  return tmp_curl_basis_ref_scalar;
}

template <typename Scalar>
typename BasisValues2<Scalar>::ConstArray_BasisIPDim
BasisValues2<Scalar>::
getCurlVectorBasisRef(const bool cache,
                      const bool force) const
{
  // Check if array already exists
  if(curl_basis_ref_vector_evaluated_ and not force)
    return curl_basis_ref_vector;

  // Reference quantities only exist for a uniform reference space
  TEUCHOS_ASSERT(hasUniformReferenceSpace());
  TEUCHOS_ASSERT(basis_layout->dimension() == 3);

  // Make sure basis is valid
  PureBasis::EElementSpace elmtspace = getElementSpace();
  TEUCHOS_ASSERT(elmtspace==PureBasis::HCURL);

  MDFieldArrayFactory af(prefix,getExtendedDimensions(),true);

  const int num_quad  = basis_layout->numPoints();
  const int num_card  = basis_layout->cardinality();
  const int num_dim   = basis_layout->dimension();

  auto tmp_curl_basis_ref_vector = af.buildStaticArray<Scalar,BASIS,IP,Dim>("dyn_curl_basis_ref_vector",num_card,num_quad,num_dim);
  auto cubature_points_uniform_ref = PHX::getNonConstDynRankViewFromConstMDField(cubature_points_uniform_ref_);

  intrepid_basis->getValues(tmp_curl_basis_ref_vector.get_view(), cubature_points_uniform_ref, Intrepid2::OPERATOR_CURL);
  PHX::Device().fence();

  // Store for later if cache is enabled
  PANZER_CACHE_DATA(curl_basis_ref_vector);

  return tmp_curl_basis_ref_vector;
}

template <typename Scalar>
typename BasisValues2<Scalar>::ConstArray_BasisIP
BasisValues2<Scalar>::
getDivVectorBasisRef(const bool cache,
                     const bool force) const
{
  // Check if array already exists
  if(div_basis_ref_evaluated_ and not force)
    return div_basis_ref;

  // Reference quantities only exist for a uniform reference space
  TEUCHOS_ASSERT(hasUniformReferenceSpace());

  // Make sure basis is valid
  PureBasis::EElementSpace elmtspace = getElementSpace();
  TEUCHOS_ASSERT(elmtspace==PureBasis::HDIV);

  MDFieldArrayFactory af(prefix,getExtendedDimensions(),true);

  const int num_quad  = basis_layout->numPoints();
  const int num_card  = basis_layout->cardinality();

  auto tmp_div_basis_ref = af.buildStaticArray<Scalar,BASIS,IP>("dyn_div_basis_ref_scalar",num_card,num_quad);
  auto cubature_points_uniform_ref = PHX::getNonConstDynRankViewFromConstMDField(cubature_points_uniform_ref_);

  intrepid_basis->getValues(tmp_div_basis_ref.get_view(), cubature_points_uniform_ref, Intrepid2::OPERATOR_DIV);
  PHX::Device().fence();

  // Store for later if cache is enabled
  PANZER_CACHE_DATA(div_basis_ref);

  return tmp_div_basis_ref;
}

template <typename Scalar>
typename BasisValues2<Scalar>::ConstArray_CellBasisDim
BasisValues2<Scalar>::
getBasisCoordinates(const bool cache,
                    const bool force) const
{
  // Check if array already exists
  if(basis_coordinates_evaluated_ and not force)
    return basis_coordinates;

  TEUCHOS_ASSERT(cell_vertex_coordinates_.size() > 0);

  const std::pair<int,int> cell_range(0,num_evaluate_cells_);
  const auto s_vertex_coordinates = Kokkos::subview(cell_vertex_coordinates_.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL());

  MDFieldArrayFactory af(prefix,getExtendedDimensions(),true);

  const int num_card  = basis_layout->cardinality();
  const int num_dim   = basis_layout->dimension();

  auto tmp_basis_coordinates = af.buildStaticArray<Scalar, Cell, BASIS, IP>("basis_coordinates",  num_cells_, num_card, num_dim);
  auto s_aux = Kokkos::subview(tmp_basis_coordinates.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL());

  // Don't forget that since we are not caching this, we have to make sure the managed view remains alive while we use the non-const wrapper
  auto const_bcr = getBasisCoordinatesRef(false);
  auto bcr = PHX::getNonConstDynRankViewFromConstMDField(const_bcr);

  Intrepid2::CellTools<PHX::Device::execution_space> cell_tools;
  cell_tools.mapToPhysicalFrame(s_aux, bcr, s_vertex_coordinates, intrepid_basis->getBaseCellTopology());
  PHX::Device().fence();

  // Store for later if cache is enabled
  PANZER_CACHE_DATA(basis_coordinates);

  return tmp_basis_coordinates;
}

template <typename Scalar>
typename BasisValues2<Scalar>::ConstArray_CellBasisIP
BasisValues2<Scalar>::
getBasisValues(const bool weighted,
               const bool cache,
               const bool force) const
{
  if(weighted){
    if(weighted_basis_scalar_evaluated_ and not force)
      return weighted_basis_scalar;
  } else
    if(basis_scalar_evaluated_ and not force)
      return basis_scalar;

  MDFieldArrayFactory af(prefix,getExtendedDimensions(),true);

  const int num_cells  = num_cells_;
  const int num_points = basis_layout->numPoints();
  const int num_card   = basis_layout->cardinality();
  const int num_dim    = basis_layout->dimension();

  if(weighted){

    TEUCHOS_ASSERT(cubature_weights_.size() > 0);

    // Get the basis_scalar values - do not cache
    const auto bv = getBasisValues(false, force);

    // Apply the weighted measure (cubature weights)
    auto tmp_weighted_basis_scalar = af.buildStaticArray<Scalar, Cell, BASIS, IP>("weighted_basis_scalar",  num_cells, num_card, num_points);

    const std::pair<int,int> cell_range(0,num_evaluate_cells_);
    auto s_aux = Kokkos::subview(tmp_weighted_basis_scalar.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL());
    auto s_cw = Kokkos::subview(cubature_weights_.get_view(),cell_range,Kokkos::ALL());
    auto s_bv = Kokkos::subview(bv.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL());

    using fst=Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>;
    fst::multiplyMeasure(s_aux,s_cw,s_bv);

    // NOTE: Weighted path has orientations already applied so doesn't
    // need the applyOrientations call at the bottom of this function.

    // Store for later if cache is enabled
    PANZER_CACHE_DATA(weighted_basis_scalar);

    return tmp_weighted_basis_scalar;

  } else {

    const auto element_space = getElementSpace();
    TEUCHOS_ASSERT(element_space == PureBasis::HVOL || element_space == PureBasis::HGRAD || element_space == PureBasis::CONST);

    // HVol requires the jacobian determinant
    if(element_space == PureBasis::HVOL){
      TEUCHOS_ASSERT(cubature_jacobian_determinant_.size() > 0);
    }

    auto cell_basis_ref_scalar = af.buildStaticArray<Scalar,BASIS,IP>("cell_basis_ref_scalar",num_card,num_points);
    auto tmp_basis_scalar = af.buildStaticArray<Scalar,Cell,BASIS,IP>("basis_scalar",num_cells,num_card,num_points);

    if(hasUniformReferenceSpace()){

      auto cubature_points_uniform_ref = PHX::getNonConstDynRankViewFromConstMDField(cubature_points_uniform_ref_);

      // Apply a single reference representation to all cells
      intrepid_basis->getValues(cell_basis_ref_scalar.get_view(),cubature_points_uniform_ref,Intrepid2::OPERATOR_VALUE);

      const std::pair<int,int> cell_range(0,num_evaluate_cells_);
      auto s_aux = Kokkos::subview(tmp_basis_scalar.get_view(), cell_range, Kokkos::ALL(), Kokkos::ALL());

      // Apply transformation (HGRAD version is just a copy operation)
      using fst=Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>;
      if(element_space == PureBasis::HVOL){
        auto s_cjd = Kokkos::subview(cubature_jacobian_determinant_.get_view(), cell_range, Kokkos::ALL());
        fst::HVOLtransformVALUE(s_aux,s_cjd,cell_basis_ref_scalar.get_view());
      }else if(element_space == PureBasis::HGRAD || element_space == PureBasis::CONST)
        fst::HGRADtransformVALUE(s_aux,cell_basis_ref_scalar.get_view());

      PHX::Device().fence();

    } else {

      // This is ugly. The algorithm is restricted to host/serial due
      // to a call to intrepid tools that require uniform reference
      // representation. For DG, CVFEM and sidesets this reference is
      // nonuniform.

      // Local allocation used for each cell
      auto cell_basis_scalar = af.buildStaticArray<Scalar,Cell,BASIS,IP>("cell_basis_scalar",1,num_card,num_points);
      auto cell_cub_points = af.buildStaticArray<Scalar,IP,Dim>("cell_cub_points",num_points,num_dim);
      auto cell_jac_det = af.buildStaticArray<Scalar,Cell,IP>("cell_jac_det",1,num_points);

      // The array factory is difficult to extend to host space
      // without extra template magic and changing a ton of code in a
      // non-backwards compatible way, so we use some of the arrays
      // above only to get derivative array sized correctly and then
      // create the mirror on host.
      auto cell_basis_scalar_host = Kokkos::create_mirror_view(cell_basis_scalar.get_view());
      auto cell_cub_points_host = Kokkos::create_mirror_view(cell_cub_points.get_view());
      auto cell_jac_det_host = Kokkos::create_mirror_view(cell_jac_det.get_view());
      auto cell_basis_ref_scalar_host = Kokkos::create_mirror_view(cell_basis_ref_scalar.get_view());
      auto cubature_points_ref_host = Kokkos::create_mirror_view(cubature_points_ref_.get_view());
      Kokkos::deep_copy(cubature_points_ref_host,cubature_points_ref_.get_view());
      auto cubature_jacobian_determinant_host = Kokkos::create_mirror_view(cubature_jacobian_determinant_.get_view());
      Kokkos::deep_copy(cubature_jacobian_determinant_host,cubature_jacobian_determinant_.get_view());
      auto tmp_basis_scalar_host = Kokkos::create_mirror_view(tmp_basis_scalar.get_view());

      // We have to iterate through cells and apply a separate reference representation for each
      for(int cell=0; cell<num_evaluate_cells_; ++cell){

        // Load external into cell-local arrays
        for(int p=0;p<num_points;++p)
          for(int d=0;d<num_dim;++d)
            cell_cub_points_host(p,d)=cubature_points_ref_host(cell,p,d);
        Kokkos::deep_copy(cell_cub_points.get_view(),cell_cub_points_host);

        // Convert to mesh space. The intrepid_basis is virtual and
        // built in another class, short of adding a "clone on host"
        // method, we will have to make this call on device. This is a
        // major performance hit.
        intrepid_basis->getValues(cell_basis_ref_scalar.get_view(),cell_cub_points.get_view(),Intrepid2::OPERATOR_VALUE);
        Kokkos::deep_copy(cell_basis_ref_scalar_host,cell_basis_ref_scalar.get_view());

        using fst=Intrepid2::FunctionSpaceTools<Kokkos::HostSpace>;

        if(element_space == PureBasis::HVOL){
          // Need the jacobian determinant for HVOL
          for(int p=0;p<num_points;++p)
            cell_jac_det_host(0,p)=cubature_jacobian_determinant_host(cell,p);

          fst::HVOLtransformVALUE(cell_basis_scalar_host,cell_jac_det_host,cell_basis_ref_scalar_host);
        } else if(element_space == PureBasis::HGRAD || element_space == PureBasis::CONST)
          fst::HGRADtransformVALUE(cell_basis_scalar_host,cell_basis_ref_scalar_host);

        // Copy cell quantity back into main array
        for(int b=0; b<num_card; ++b)
          for(int p=0; p<num_points; ++p)
            tmp_basis_scalar_host(cell,b,p) = cell_basis_scalar_host(0,b,p);

        Kokkos::deep_copy(tmp_basis_scalar.get_view(),tmp_basis_scalar_host);
      }
    }

    // NOTE: weighted already has orientations applied, so this code
    // should only be reached if non-weighted. This is a by-product of
    // the two construction paths in panzer. Unification should help
    // fix the logic.

    if(orientations_ != Teuchos::null)
      applyOrientationsImpl<Scalar>(num_orientations_cells_, tmp_basis_scalar.get_view(), *orientations_->getOrientations(), *intrepid_basis);

    // Store for later if cache is enabled
    PANZER_CACHE_DATA(basis_scalar);

    return tmp_basis_scalar;

  }

}

template <typename Scalar>
typename BasisValues2<Scalar>::ConstArray_CellBasisIPDim
BasisValues2<Scalar>::
getVectorBasisValues(const bool weighted,
                     const bool cache,
                     const bool force) const
{
  if(weighted){
    if(weighted_basis_vector_evaluated_ and not force)
      return weighted_basis_vector;
  } else
    if(basis_vector_evaluated_ and not force)
      return basis_vector;

  MDFieldArrayFactory af(prefix,getExtendedDimensions(),true);

  const int num_cells  = num_cells_;
  const int num_points = basis_layout->numPoints();
  const int num_card   = basis_layout->cardinality();
  const int num_dim    = basis_layout->dimension();

  if(weighted){

    TEUCHOS_ASSERT(cubature_weights_.size() > 0);

    // Get the basis_scalar values - do not cache
    const auto bv = getVectorBasisValues(false, force);

    // Apply the weighted measure (cubature weights)
    auto tmp_weighted_basis_vector = af.buildStaticArray<Scalar, Cell, BASIS, IP, Dim>("weighted_basis_vector",  num_cells, num_card, num_points, num_dim);

    const std::pair<int,int> cell_range(0,num_evaluate_cells_);
    auto s_aux = Kokkos::subview(tmp_weighted_basis_vector.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
    auto s_cw = Kokkos::subview(cubature_weights_.get_view(),cell_range,Kokkos::ALL());
    auto s_bv = Kokkos::subview(bv.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());

    using fst=Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>;
    fst::multiplyMeasure(s_aux, s_cw, s_bv);

    // Store for later if cache is enabled
    PANZER_CACHE_DATA(weighted_basis_vector);

    return tmp_weighted_basis_vector;

  } else { // non-weighted

    const auto element_space = getElementSpace();
    TEUCHOS_ASSERT(element_space == PureBasis::HCURL || element_space == PureBasis::HDIV);
    TEUCHOS_ASSERT(num_dim != 1);

    // HDIV and HCURL have unique jacobian requirements
    if(element_space == PureBasis::HCURL){
      TEUCHOS_ASSERT(cubature_jacobian_inverse_.size() > 0);
    } else if(element_space == PureBasis::HDIV){
      TEUCHOS_ASSERT(cubature_jacobian_.size() > 0 && cubature_jacobian_determinant_.size() > 0);
    }

    auto cell_basis_ref_vector = af.buildStaticArray<Scalar,BASIS,IP,Dim>("cell_basis_ref_scalar",num_card,num_points,num_dim);
    auto tmp_basis_vector = af.buildStaticArray<Scalar,Cell,BASIS,IP,Dim>("basis_vector",num_cells,num_card,num_points,num_dim);

    if(hasUniformReferenceSpace()){

      auto cubature_points_uniform_ref = PHX::getNonConstDynRankViewFromConstMDField(cubature_points_uniform_ref_);

      // Apply a single reference representation to all cells
      intrepid_basis->getValues(cell_basis_ref_vector.get_view(),cubature_points_uniform_ref,Intrepid2::OPERATOR_VALUE);

      const std::pair<int,int> cell_range(0,num_evaluate_cells_);
      auto s_aux = Kokkos::subview(tmp_basis_vector.get_view(), cell_range, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());

      // Apply transformation (HGRAD version is just a copy operation)
      using fst=Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>;
      if(element_space == PureBasis::HCURL){

        auto s_jac_inv = Kokkos::subview(cubature_jacobian_inverse_.get_view(), cell_range, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());

        fst::HCURLtransformVALUE(s_aux,s_jac_inv,cell_basis_ref_vector.get_view());
      } else if(element_space == PureBasis::HDIV){

        auto s_jac = Kokkos::subview(cubature_jacobian_.get_view(), cell_range, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
        auto s_jac_det = Kokkos::subview(cubature_jacobian_determinant_.get_view(), cell_range, Kokkos::ALL());

        fst::HDIVtransformVALUE(s_aux,s_jac, s_jac_det, cell_basis_ref_vector.get_view());
      }

      PHX::Device().fence();

    } else {

      // This is ugly. The algorithm is restricted to host/serial due
      // to intrepid tools that requiring uniform reference
      // representation. For DG, CVFEM and sidesets this reference is
      // nonuniform.

      // Local allocation used for each cell
      auto cell_basis_vector = af.buildStaticArray<Scalar,Cell,BASIS,IP,Dim>("cell_basis_vector",1,num_card,num_points,num_dim);
      auto cell_cub_points = af.buildStaticArray<Scalar,IP,Dim>("cell_cub_points",num_points,num_dim);
      auto cell_jac_det = af.buildStaticArray<Scalar,Cell,IP>("cell_jac_det",1,num_points);
      auto cell_jac = af.buildStaticArray<Scalar,Cell,IP,Dim,Dim>("cell_jac",1,num_points,num_dim,num_dim);
      auto cell_jac_inv = af.buildStaticArray<Scalar,Cell,IP,Dim,Dim>("cell_jac_inv",1,num_points,num_dim,num_dim);

      // The array factory is difficult to extend to host space
      // without extra template magic and changing a ton of code in a
      // non-backwards compatible way, so we use some of the arrays
      // above only to get derivative array sized correctly and then
      // create the mirror on host.
      auto cell_basis_vector_host = Kokkos::create_mirror_view(cell_basis_vector.get_view());
      auto cell_cub_points_host = Kokkos::create_mirror_view(cell_cub_points.get_view());
      auto cell_jac_det_host = Kokkos::create_mirror_view(cell_jac_det.get_view());
      auto cell_jac_host = Kokkos::create_mirror_view(cell_jac.get_view());
      auto cell_jac_inv_host = Kokkos::create_mirror_view(cell_jac_inv.get_view());
      auto cell_basis_ref_vector_host = Kokkos::create_mirror_view(cell_basis_ref_vector.get_view());
      auto cubature_points_ref_host = Kokkos::create_mirror_view(cubature_points_ref_.get_view());
      Kokkos::deep_copy(cubature_points_ref_host,cubature_points_ref_.get_view());
      auto cubature_jacobian_host = Kokkos::create_mirror_view(cubature_jacobian_.get_view());
      Kokkos::deep_copy(cubature_jacobian_host,cubature_jacobian_.get_view());
      auto cubature_jacobian_inverse_host = Kokkos::create_mirror_view(cubature_jacobian_inverse_.get_view());
      Kokkos::deep_copy(cubature_jacobian_inverse_host,cubature_jacobian_inverse_.get_view());
      auto cubature_jacobian_determinant_host = Kokkos::create_mirror_view(cubature_jacobian_determinant_.get_view());
      Kokkos::deep_copy(cubature_jacobian_determinant_host,cubature_jacobian_determinant_.get_view());
      auto tmp_basis_vector_host = Kokkos::create_mirror_view(tmp_basis_vector.get_view());

      // We have to iterate through cells and apply a separate reference representation for each
      for(int cell=0; cell<num_evaluate_cells_; ++cell){

        // Load external into cell-local arrays
        for(int p=0;p<num_points;++p)
          for(int d=0;d<num_dim;++d)
            cell_cub_points_host(p,d)=cubature_points_ref_host(cell,p,d);
        Kokkos::deep_copy(cell_cub_points.get_view(),cell_cub_points_host);

        // Convert to mesh space. The intrepid_basis is virtual and
        // built in another class, short of adding a "clone on host"
        // method, we will have to make this call on device. This is a
        // major performance hit.
        intrepid_basis->getValues(cell_basis_ref_vector.get_view(),cell_cub_points.get_view(),Intrepid2::OPERATOR_VALUE);
        Kokkos::deep_copy(cell_basis_ref_vector_host,cell_basis_ref_vector.get_view());

        using fst=Intrepid2::FunctionSpaceTools<Kokkos::HostSpace>;

        if(element_space == PureBasis::HCURL){
          // Need the jacobian inverse for HCurl
          for(int p=0;p<num_points;++p)
            for(int i=0; i<num_dim; ++i)
              for(int j=0; j<num_dim; ++j)
                cell_jac_inv_host(0,p,i,j)=cubature_jacobian_inverse_host(cell,p,i,j);

          fst::HCURLtransformVALUE(cell_basis_vector_host,cell_jac_inv_host,cell_basis_ref_vector_host);
        } else {
          // Need the jacobian for HDiv
          for(int p=0;p<num_points;++p)
            for(int i=0; i<num_dim; ++i)
              for(int j=0; j<num_dim; ++j)
                cell_jac_host(0,p,i,j)=cubature_jacobian_host(cell,p,i,j);

          for(int p=0;p<num_points;++p)
            cell_jac_det_host(0,p)=cubature_jacobian_determinant_host(cell,p);

          fst::HDIVtransformVALUE(cell_basis_vector_host,cell_jac_host,cell_jac_det_host,cell_basis_ref_vector_host);
        }

        // Copy cell quantity back into main array
        for(int b=0; b<num_card; ++b)
          for(int p=0; p<num_points; ++p)
            for(int d=0; d<num_dim; ++d)
            tmp_basis_vector_host(cell,b,p,d) = cell_basis_vector_host(0,b,p,d);

        Kokkos::deep_copy(tmp_basis_vector.get_view(),tmp_basis_vector_host);
      }
    }

    if(orientations_ != Teuchos::null)
      applyOrientationsImpl<Scalar>(num_orientations_cells_, tmp_basis_vector.get_view(), *orientations_->getOrientations(), *intrepid_basis);

    // Store for later if cache is enabled
    PANZER_CACHE_DATA(basis_vector);

    return tmp_basis_vector;

  }

}

template <typename Scalar>
typename BasisValues2<Scalar>::ConstArray_CellBasisIPDim
BasisValues2<Scalar>::
getGradBasisValues(const bool weighted,
                   const bool cache,
                   const bool force) const
{
  if(weighted){
    if(weighted_grad_basis_evaluated_ and not force)
      return weighted_grad_basis;
  } else
    if(grad_basis_evaluated_ and not force)
      return grad_basis;

  MDFieldArrayFactory af(prefix,getExtendedDimensions(),true);

  const int num_cells  = num_cells_;
  const int num_points = basis_layout->numPoints();
  const int num_card   = basis_layout->cardinality();
  const int num_dim    = basis_layout->dimension();

  if(weighted){

    TEUCHOS_ASSERT(cubature_weights_.size() > 0);

    // Get the basis_scalar values - do not cache
    const auto bv = getGradBasisValues(false, force);

    // Apply the weighted measure (cubature weights)
    auto tmp_weighted_grad_basis = af.buildStaticArray<Scalar, Cell, BASIS, IP, Dim>("weighted_grad_basis",  num_cells, num_card, num_points, num_dim);

    const std::pair<int,int> cell_range(0,num_evaluate_cells_);
    auto s_aux = Kokkos::subview(tmp_weighted_grad_basis.get_view(), cell_range, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
    auto s_cw = Kokkos::subview(cubature_weights_.get_view(), cell_range, Kokkos::ALL());
    auto s_bv = Kokkos::subview(bv.get_view(), cell_range, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());

    using fst=Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>;
    fst::multiplyMeasure(s_aux,s_cw,s_bv);

    // Store for later if cache is enabled
    PANZER_CACHE_DATA(weighted_grad_basis);

    return tmp_weighted_grad_basis;

  } else {

    TEUCHOS_ASSERT(cubature_jacobian_inverse_.size() > 0);

    const auto element_space = getElementSpace();
    TEUCHOS_ASSERT(element_space == PureBasis::CONST || element_space == PureBasis::HGRAD);

    auto cell_grad_basis_ref = af.buildStaticArray<Scalar,BASIS,IP,Dim>("cell_grad_basis_ref",num_card,num_points,num_dim);
    auto tmp_grad_basis = af.buildStaticArray<Scalar,Cell,BASIS,IP,Dim>("basis_scalar",num_cells,num_card,num_points,num_dim);

    if(hasUniformReferenceSpace()){

      auto cubature_points_uniform_ref = PHX::getNonConstDynRankViewFromConstMDField(cubature_points_uniform_ref_);

      // Apply a single reference representation to all cells
      intrepid_basis->getValues(cell_grad_basis_ref.get_view(),cubature_points_uniform_ref,Intrepid2::OPERATOR_GRAD);

      const std::pair<int,int> cell_range(0,num_evaluate_cells_);
      auto s_aux = Kokkos::subview(tmp_grad_basis.get_view(), cell_range, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
      auto s_jac_inv = Kokkos::subview(cubature_jacobian_inverse_.get_view(), cell_range, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());

      // Apply transformation
      using fst=Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>;
      fst::HGRADtransformGRAD(s_aux, s_jac_inv,cell_grad_basis_ref.get_view());

      PHX::Device().fence();

    } else {

      // This is ugly. The algorithm is restricted to host/serial due
      // to intrepid tools that requiring uniform reference
      // representation. For DG, CVFEM and sidesets this reference is
      // nonuniform.

      // Local allocation used for each cell
      auto cell_grad_basis = af.buildStaticArray<Scalar,Cell,BASIS,IP,Dim>("cell_grad_basis",1,num_card,num_points,num_dim);
      auto cell_cub_points = af.buildStaticArray<Scalar,IP,Dim>("cell_cub_points",num_points,num_dim);
      auto cell_jac_inv = af.buildStaticArray<Scalar,Cell,IP,Dim,Dim>("cell_jac_inv",1,num_points,num_dim,num_dim);

      // The array factory is difficult to extend to host space
      // without extra template magic and changing a ton of code in a
      // non-backwards compatible way, so we use some of the arrays
      // above only to get derivative array sized correctly and then
      // create the mirror on host.

      // auto cell_basis_ref_vector = af.buildStaticArray<Scalar,BASIS,IP,Dim>("cell_basis_ref_scalar",num_card,num_points,num_dim);

      auto cell_cub_points_host = Kokkos::create_mirror_view(cell_cub_points.get_view());
      auto cubature_points_ref_host = Kokkos::create_mirror_view(cubature_points_ref_.get_view());
      Kokkos::deep_copy(cubature_points_ref_host,cubature_points_ref_.get_view());
      auto cell_jac_inv_host = Kokkos::create_mirror_view(cell_jac_inv.get_view());
      auto cubature_jacobian_inverse_host = Kokkos::create_mirror_view(cubature_jacobian_inverse_.get_view());
      Kokkos::deep_copy(cubature_jacobian_inverse_host,cubature_jacobian_inverse_.get_view());
      auto cell_grad_basis_ref_host = Kokkos::create_mirror_view(cell_grad_basis_ref.get_view());
      auto cell_grad_basis_host = Kokkos::create_mirror_view(cell_grad_basis.get_view());
      auto tmp_grad_basis_host = Kokkos::create_mirror_view(tmp_grad_basis.get_view());

      // We have to iterate through cells and apply a separate reference representation for each
      for(int cell=0; cell<num_evaluate_cells_; ++cell){

        // Load external into cell-local arrays
        for(int p=0;p<num_points;++p)
          for(int d=0;d<num_dim;++d)
            cell_cub_points_host(p,d)=cubature_points_ref_host(cell,p,d);

        // Convert to mesh space. The intrepid_basis is virtual and
        // built in another class, short of adding a "clone on host"
        // method, we will have to make this call on device. This is a
        // major performance hit.
        Kokkos::deep_copy(cell_cub_points.get_view(),cell_cub_points_host);
        intrepid_basis->getValues(cell_grad_basis_ref.get_view(),cell_cub_points.get_view(),Intrepid2::OPERATOR_GRAD);
        Kokkos::deep_copy(cell_grad_basis_ref_host,cell_grad_basis_ref.get_view());

        for(int p=0;p<num_points;++p)
          for(int d=0;d<num_dim;++d)
            for(int d2=0;d2<num_dim;++d2)
              cell_jac_inv_host(0,p,d,d2)=cubature_jacobian_inverse_host(cell,p,d,d2);

        using fst=Intrepid2::FunctionSpaceTools<Kokkos::HostSpace>;
        fst::HGRADtransformGRAD(cell_grad_basis_host,cell_jac_inv_host,cell_grad_basis_ref_host);
        // PHX::Device().fence();

        // Copy cell quantity back into main array
        for(int b=0; b<num_card; ++b)
          for(int p=0; p<num_points; ++p)
            for(int d=0; d<num_dim; ++d)
              tmp_grad_basis_host(cell,b,p,d) = cell_grad_basis_host(0,b,p,d);
        Kokkos::deep_copy(tmp_grad_basis.get_view(),tmp_grad_basis_host);
      }
    }

    if(orientations_ != Teuchos::null)
      applyOrientationsImpl<Scalar>(num_orientations_cells_, tmp_grad_basis.get_view(), *orientations_->getOrientations(), *intrepid_basis);

    // Store for later if cache is enabled
    PANZER_CACHE_DATA(grad_basis);

    return tmp_grad_basis;

  }

}

template <typename Scalar>
typename BasisValues2<Scalar>::ConstArray_CellBasisIP
BasisValues2<Scalar>::
getCurl2DVectorBasis(const bool weighted,
                     const bool cache,
                     const bool force) const
{
  if(weighted){
    if(weighted_curl_basis_scalar_evaluated_ and not force)
      return weighted_curl_basis_scalar;
  } else
    if(curl_basis_scalar_evaluated_ and not force)
      return curl_basis_scalar;

  MDFieldArrayFactory af(prefix,getExtendedDimensions(),true);

  const int num_cells  = num_cells_;
  const int num_points = basis_layout->numPoints();
  const int num_card   = basis_layout->cardinality();
  const int num_dim    = basis_layout->dimension();

  if(weighted){

    TEUCHOS_ASSERT(cubature_weights_.size() > 0);

    // Get the basis_scalar values - do not cache
    const auto bv = getCurl2DVectorBasis(false, force);

    // Apply the weighted measure (cubature weights)
    auto tmp_weighted_curl_basis_scalar = af.buildStaticArray<Scalar, Cell, BASIS, IP>("weighted_curl_basis_scalar",  num_cells, num_card, num_points);

    const std::pair<int,int> cell_range(0,num_evaluate_cells_);
    auto s_aux = Kokkos::subview(tmp_weighted_curl_basis_scalar.get_view(), cell_range, Kokkos::ALL(), Kokkos::ALL());
    auto s_cw = Kokkos::subview(cubature_weights_.get_view(), cell_range, Kokkos::ALL());
    auto s_bv = Kokkos::subview(bv.get_view(), cell_range, Kokkos::ALL(), Kokkos::ALL());

    using fst=Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>;
    fst::multiplyMeasure(s_aux,s_cw,s_bv);

    // Store for later if cache is enabled
    PANZER_CACHE_DATA(weighted_curl_basis_scalar);

    return tmp_weighted_curl_basis_scalar;

  } else {

    TEUCHOS_ASSERT(cubature_jacobian_determinant_.size() > 0);
    TEUCHOS_ASSERT(num_dim == 2);

    const auto element_space = getElementSpace();
    TEUCHOS_ASSERT(element_space == PureBasis::HCURL);

    auto cell_curl_basis_ref_scalar =  af.buildStaticArray<Scalar,BASIS,IP>("cell_curl_basis_ref_scalar",num_card,num_points);
    auto tmp_curl_basis_scalar = af.buildStaticArray<Scalar,Cell,BASIS,IP>("curl_basis_scalar",num_cells,num_card,num_points);

    if(hasUniformReferenceSpace()){

      auto cubature_points_uniform_ref = PHX::getNonConstDynRankViewFromConstMDField(cubature_points_uniform_ref_);

      intrepid_basis->getValues(cell_curl_basis_ref_scalar.get_view(),cubature_points_uniform_ref,Intrepid2::OPERATOR_CURL);

      const std::pair<int,int> cell_range(0,num_evaluate_cells_);
      auto s_aux = Kokkos::subview(tmp_curl_basis_scalar.get_view(), cell_range, Kokkos::ALL(), Kokkos::ALL());
      auto s_jac_det = Kokkos::subview(cubature_jacobian_determinant_.get_view(), cell_range, Kokkos::ALL());

      // note only volume deformation is needed!
      // this relates directly to this being in
      // the divergence space in 2D!
      using fst=Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>;
      fst::HDIVtransformDIV(s_aux,s_jac_det,cell_curl_basis_ref_scalar.get_view());
      PHX::Device().fence();

    } else {

      // This is ugly. The algorithm is restricted to host/serial due
      // to intrepid tools that requiring uniform reference
      // representation. For DG, CVFEM and sidesets this reference is
      // nonuniform.

      // Local allocation used for each cell
      auto cell_curl = af.buildStaticArray<Scalar,Cell,BASIS,IP>("cell_curl",1,num_card,num_points);
      auto cell_cub_points = af.buildStaticArray<Scalar,IP,Dim>("cell_cub_points",num_points,num_dim);
      auto cell_jac_det = af.buildStaticArray<Scalar,Cell,IP>("cell_jac_det",1,num_points);

      // The array factory is difficult to extend to host space
      // without extra template magic and changing a ton of code in a
      // non-backwards compatible way, so we use some of the arrays
      // above only to get derivative array sized correctly and then
      // create the mirror on host.
      auto cell_cub_points_host = Kokkos::create_mirror_view(cell_cub_points.get_view());
      auto cubature_points_ref_host = Kokkos::create_mirror_view(cubature_points_ref_.get_view());
      Kokkos::deep_copy(cubature_points_ref_host,cubature_points_ref_.get_view());
      auto cell_jac_det_host = Kokkos::create_mirror_view(cell_jac_det.get_view());
      auto cubature_jacobian_determinant_host = Kokkos::create_mirror_view(cubature_jacobian_determinant_.get_view());
      Kokkos::deep_copy(cubature_jacobian_determinant_host,cubature_jacobian_determinant_.get_view());
      auto cell_curl_basis_ref_scalar_host = Kokkos::create_mirror_view(cell_curl_basis_ref_scalar.get_view());
      auto cell_curl_host = Kokkos::create_mirror_view(cell_curl.get_view());
      auto tmp_curl_basis_scalar_host = Kokkos::create_mirror_view(tmp_curl_basis_scalar.get_view());

      // We have to iterate through cells and apply a separate reference representation for each
      for(int cell=0; cell<num_evaluate_cells_; ++cell){

        // Load external into cell-local arrays
        for(int p=0;p<num_points;++p)
          for(int d=0;d<num_dim;++d)
            cell_cub_points_host(p,d)=cubature_points_ref_host(cell,p,d);
        for(int p=0;p<num_points;++p)
          cell_jac_det_host(0,p)=cubature_jacobian_determinant_host(cell,p);

        
        // Convert to mesh space. The intrepid_basis is virtual and
        // built in another class, short of adding a "clone on host"
        // method, we will have to make this call on device. This is a
        // major performance hit.
        Kokkos::deep_copy(cell_cub_points.get_view(),cell_cub_points_host);
        intrepid_basis->getValues(cell_curl_basis_ref_scalar.get_view(),cell_cub_points.get_view(),Intrepid2::OPERATOR_CURL);
        Kokkos::deep_copy(cell_curl_basis_ref_scalar_host,cell_curl_basis_ref_scalar.get_view());

        // note only volume deformation is needed!
        // this relates directly to this being in
        // the divergence space in 2D!
        using fst=Intrepid2::FunctionSpaceTools<Kokkos::HostSpace>;
        fst::HDIVtransformDIV(cell_curl_host,cell_jac_det_host,cell_curl_basis_ref_scalar_host);
        PHX::Device().fence();

        // Copy cell quantity back into main array
        for(int b=0; b<num_card; ++b)
          for(int p=0; p<num_points; ++p)
            tmp_curl_basis_scalar_host(cell,b,p) = cell_curl_host(0,b,p);
        Kokkos::deep_copy(tmp_curl_basis_scalar.get_view(),tmp_curl_basis_scalar_host);
      }
    }

    if(orientations_ != Teuchos::null)
      applyOrientationsImpl<Scalar>(num_orientations_cells_, tmp_curl_basis_scalar.get_view(), *orientations_->getOrientations(), *intrepid_basis);

    // Store for later if cache is enabled
    PANZER_CACHE_DATA(curl_basis_scalar);

    return tmp_curl_basis_scalar;

  }

}

template <typename Scalar>
typename BasisValues2<Scalar>::ConstArray_CellBasisIPDim
BasisValues2<Scalar>::
getCurlVectorBasis(const bool weighted,
                   const bool cache,
                   const bool force) const
{
  if(weighted){
    if(weighted_curl_basis_vector_evaluated_ and not force)
      return weighted_curl_basis_vector;
  } else
    if(curl_basis_vector_evaluated_ and not force)
      return curl_basis_vector;

  MDFieldArrayFactory af(prefix,getExtendedDimensions(),true);

  const int num_cells  = num_cells_;
  const int num_points = basis_layout->numPoints();
  const int num_card   = basis_layout->cardinality();
  const int num_dim    = basis_layout->dimension();

  if(weighted){

    TEUCHOS_ASSERT(cubature_weights_.size() > 0);

    // Get the basis_scalar values - do not cache
    const auto bv = getCurlVectorBasis(false, force);

    // Apply the weighted measure (cubature weights)
    auto tmp_weighted_curl_basis_vector = af.buildStaticArray<Scalar, Cell, BASIS, IP, Dim>("weighted_curl_basis_vector",  num_cells, num_card, num_points, num_dim);

    const std::pair<int,int> cell_range(0,num_evaluate_cells_);
    auto s_aux = Kokkos::subview(tmp_weighted_curl_basis_vector.get_view(), cell_range, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
    auto s_cw = Kokkos::subview(cubature_weights_.get_view(), cell_range, Kokkos::ALL());
    auto s_bv = Kokkos::subview(bv.get_view(), cell_range, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());

    using fst=Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>;
    fst::multiplyMeasure(s_aux, s_cw, s_bv);

    // Store for later if cache is enabled
    PANZER_CACHE_DATA(weighted_curl_basis_vector);

    return tmp_weighted_curl_basis_vector;

  } else {

    TEUCHOS_ASSERT(cubature_jacobian_determinant_.size() > 0);
    TEUCHOS_ASSERT(num_dim == 3);

    const auto element_space = getElementSpace();
    TEUCHOS_ASSERT(element_space == PureBasis::HCURL);

    auto cell_curl_basis_ref_vector =  af.buildStaticArray<Scalar,BASIS,IP,Dim>("cell_curl_basis_ref_vector",num_card,num_points,num_dim);
    auto tmp_curl_basis_vector = af.buildStaticArray<Scalar,Cell,BASIS,IP,Dim>("curl_basis_vector",num_cells,num_card,num_points,num_dim);

    if(hasUniformReferenceSpace()){

      auto cubature_points_uniform_ref = PHX::getNonConstDynRankViewFromConstMDField(cubature_points_uniform_ref_);

      intrepid_basis->getValues(cell_curl_basis_ref_vector.get_view(),cubature_points_uniform_ref,Intrepid2::OPERATOR_CURL);

      const std::pair<int,int> cell_range(0,num_evaluate_cells_);
      auto s_aux = Kokkos::subview(tmp_curl_basis_vector.get_view(), cell_range, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
      auto s_jac = Kokkos::subview(cubature_jacobian_.get_view(), cell_range, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
      auto s_jac_det = Kokkos::subview(cubature_jacobian_determinant_.get_view(), cell_range, Kokkos::ALL());

      using fst=Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>;
      fst::HCURLtransformCURL(s_aux, s_jac, s_jac_det,cell_curl_basis_ref_vector.get_view());
      PHX::Device().fence();

    } else {

      // This is ugly. The algorithm is restricted to host/serial due
      // to intrepid tools that requiring uniform reference
      // representation. For DG, CVFEM and sidesets this reference is
      // nonuniform.

      // Local allocation used for each cell
      auto cell_curl = af.buildStaticArray<Scalar,Cell,BASIS,IP,Dim>("cell_curl",1,num_card,num_points,num_dim);
      auto cell_cub_points = af.buildStaticArray<Scalar,IP,Dim>("cell_cub_points",num_points,num_dim);
      auto cell_jac = af.buildStaticArray<Scalar,Cell,IP,Dim,Dim>("cell_jac",1,num_points,num_dim,num_dim);
      auto cell_jac_det = af.buildStaticArray<Scalar,Cell,IP>("cell_jac_det",1,num_points);

      // The array factory is difficult to extend to host space
      // without extra template magic and changing a ton of code in a
      // non-backwards compatible way, so we use some of the arrays
      // above only to get derivative array sized correctly and then
      // create the mirror on host.
      auto cell_cub_points_host = Kokkos::create_mirror_view(cell_cub_points.get_view());
      auto cubature_points_ref_host = Kokkos::create_mirror_view(cubature_points_ref_.get_view());
      Kokkos::deep_copy(cubature_points_ref_host,cubature_points_ref_.get_view());
      auto cell_jac_host = Kokkos::create_mirror_view(cell_jac.get_view());
      auto cubature_jacobian_host = Kokkos::create_mirror_view(cubature_jacobian_.get_view());
      Kokkos::deep_copy(cubature_jacobian_host,cubature_jacobian_.get_view());
      auto cell_jac_det_host = Kokkos::create_mirror_view(cell_jac_det.get_view());
      auto cubature_jacobian_determinant_host = Kokkos::create_mirror_view(cubature_jacobian_determinant_.get_view());
      Kokkos::deep_copy(cubature_jacobian_determinant_host,cubature_jacobian_determinant_.get_view());
      auto cell_curl_basis_ref_vector_host = Kokkos::create_mirror_view(cell_curl_basis_ref_vector.get_view());
      auto cell_curl_host = Kokkos::create_mirror_view(cell_curl.get_view());
      auto tmp_curl_basis_vector_host = Kokkos::create_mirror_view(tmp_curl_basis_vector.get_view());

      // We have to iterate through cells and apply a separate reference representation for each
      for(int cell=0; cell<num_evaluate_cells_; ++cell){

        // Load external into cell-local arrays
        for(int p=0;p<num_points;++p)
          for(int d=0;d<num_dim;++d)
            cell_cub_points_host(p,d)=cubature_points_ref_host(cell,p,d);
        for(int p=0;p<num_points;++p)
          for(int d=0;d<num_dim;++d)
            for(int d2=0;d2<num_dim;++d2)
              cell_jac_host(0,p,d,d2)=cubature_jacobian_host(cell,p,d,d2);
        for(int p=0;p<num_points;++p)
          cell_jac_det_host(0,p)=cubature_jacobian_determinant_host(cell,p);

        // Convert to mesh space. The intrepid_basis is virtual and
        // built in another class, short of adding a "clone on host"
        // method, we will have to make this call on device. This is a
        // major performance hit.
        Kokkos::deep_copy(cell_cub_points.get_view(),cell_cub_points_host);
        intrepid_basis->getValues(cell_curl_basis_ref_vector.get_view(),cell_cub_points.get_view(),Intrepid2::OPERATOR_CURL);
        Kokkos::deep_copy(cell_curl_basis_ref_vector_host,cell_curl_basis_ref_vector.get_view());

        using fst=Intrepid2::FunctionSpaceTools<Kokkos::HostSpace>;
        fst::HCURLtransformCURL(cell_curl_host,cell_jac_host,cell_jac_det_host,cell_curl_basis_ref_vector_host);
        // PHX::Device().fence();

        // Copy cell quantity back into main array
        for(int b=0; b<num_card; ++b)
          for(int p=0; p<num_points; ++p)
            for(int d=0;d<num_dim;++d)
              tmp_curl_basis_vector_host(cell,b,p,d) = cell_curl_host(0,b,p,d);
        Kokkos::deep_copy(tmp_curl_basis_vector.get_view(),tmp_curl_basis_vector_host);
      }
    }

    if(orientations_ != Teuchos::null)
      applyOrientationsImpl<Scalar>(num_orientations_cells_, tmp_curl_basis_vector.get_view(), *orientations_->getOrientations(), *intrepid_basis);

    // Store for later if cache is enabled
    PANZER_CACHE_DATA(curl_basis_vector);

    return tmp_curl_basis_vector;

  }

}

template <typename Scalar>
typename BasisValues2<Scalar>::ConstArray_CellBasisIP
BasisValues2<Scalar>::
getDivVectorBasis(const bool weighted,
                  const bool cache,
                  const bool force) const
{
  if(weighted){
    if(weighted_div_basis_evaluated_ and not force)
      return weighted_div_basis;
  } else
    if(div_basis_evaluated_ and not force)
      return div_basis;

  MDFieldArrayFactory af(prefix,getExtendedDimensions(),true);

  const int num_cells  = num_cells_;
  const int num_points = basis_layout->numPoints();
  const int num_card   = basis_layout->cardinality();
  const int num_dim    = basis_layout->dimension();

  if(weighted){

    TEUCHOS_ASSERT(cubature_weights_.size() > 0);

    // Get the basis_scalar values - do not cache
    const auto bv = getDivVectorBasis(false, force);

    // Apply the weighted measure (cubature weights)
    auto tmp_weighted_div_basis = af.buildStaticArray<Scalar, Cell, BASIS, IP>("weighted_div_basis",  num_cells, num_card, num_points);

    const std::pair<int,int> cell_range(0,num_evaluate_cells_);
    auto s_aux = Kokkos::subview(tmp_weighted_div_basis.get_view(), cell_range, Kokkos::ALL(), Kokkos::ALL());
    auto s_cw = Kokkos::subview(cubature_weights_.get_view(), cell_range, Kokkos::ALL());
    auto s_bv = Kokkos::subview(bv.get_view(), cell_range, Kokkos::ALL(), Kokkos::ALL());

    using fst=Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>;
    fst::multiplyMeasure(s_aux, s_cw, s_bv);

    // Store for later if cache is enabled
    PANZER_CACHE_DATA(weighted_div_basis);

    return tmp_weighted_div_basis;

  } else {

    TEUCHOS_ASSERT(cubature_jacobian_determinant_.size() > 0);

    const auto element_space = getElementSpace();
    TEUCHOS_ASSERT(element_space == PureBasis::HDIV);

    auto cell_div_basis_ref = af.buildStaticArray<Scalar,BASIS,IP>("cell_div_basis_ref",num_card,num_points);
    auto tmp_div_basis = af.buildStaticArray<Scalar,Cell,BASIS,IP>("div_basis",num_cells,num_card,num_points);

    if(hasUniformReferenceSpace()){

      auto cubature_points_uniform_ref = PHX::getNonConstDynRankViewFromConstMDField(cubature_points_uniform_ref_);

      intrepid_basis->getValues(cell_div_basis_ref.get_view(),cubature_points_uniform_ref,Intrepid2::OPERATOR_DIV);

      const std::pair<int,int> cell_range(0,num_evaluate_cells_);
      auto s_aux = Kokkos::subview(tmp_div_basis.get_view(), cell_range, Kokkos::ALL(), Kokkos::ALL());
      auto s_jac_det = Kokkos::subview(cubature_jacobian_determinant_.get_view(), cell_range, Kokkos::ALL());

      using fst=Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>;
      fst::HDIVtransformDIV(s_aux,s_jac_det,cell_div_basis_ref.get_view());
      PHX::Device().fence();

    } else {

      // This is ugly. The algorithm is restricted to host/serial due
      // to intrepid tools that requiring uniform reference
      // representation. For DG, CVFEM and sidesets this reference is
      // nonuniform.

      // Local allocation used for each cell
      auto cell_div_basis = af.buildStaticArray<Scalar,Cell,BASIS,IP>("cell_div_basis",1,num_card,num_points);
      auto cell_cub_points = af.buildStaticArray<Scalar,IP,Dim>("cell_cub_points",num_points,num_dim);
      auto cell_jac_det = af.buildStaticArray<Scalar,Cell,IP>("cell_jac_det",1,num_points);

      // The array factory is difficult to extend to host space
      // without extra template magic and changing a ton of code in a
      // non-backwards compatible way, so we use some of the arrays
      // above only to get derivative array sized correctly and then
      // create the mirror on host.
      auto cell_cub_points_host = Kokkos::create_mirror_view(cell_cub_points.get_view());
      auto cubature_points_ref_host = Kokkos::create_mirror_view(cubature_points_ref_.get_view());
      Kokkos::deep_copy(cubature_points_ref_host,cubature_points_ref_.get_view());
      auto cell_jac_det_host = Kokkos::create_mirror_view(cell_jac_det.get_view());
      auto cubature_jacobian_determinant_host = Kokkos::create_mirror_view(cubature_jacobian_determinant_.get_view());
      Kokkos::deep_copy(cubature_jacobian_determinant_host,cubature_jacobian_determinant_.get_view());
      auto cell_div_basis_ref_host = Kokkos::create_mirror_view(cell_div_basis_ref.get_view());
      auto cell_div_basis_host = Kokkos::create_mirror_view(cell_div_basis.get_view());
      auto tmp_div_basis_host = Kokkos::create_mirror_view(tmp_div_basis.get_view());

      // We have to iterate through cells and apply a separate reference representation for each
      for(int cell=0; cell<num_evaluate_cells_; ++cell){

        // Load external into cell-local arrays
        for(int p=0;p<num_points;++p)
          for(int d=0;d<num_dim;++d)
            cell_cub_points_host(p,d)=cubature_points_ref_host(cell,p,d);
        for(int p=0;p<num_points;++p)
          cell_jac_det_host(0,p)=cubature_jacobian_determinant_host(cell,p);

        // Convert to mesh space. The intrepid_basis is virtual and
        // built in another class, short of adding a "clone on host"
        // method, we will have to make this call on device. This is a
        // major performance hit.
        Kokkos::deep_copy(cell_cub_points.get_view(),cell_cub_points_host);
        intrepid_basis->getValues(cell_div_basis_ref.get_view(),cell_cub_points.get_view(),Intrepid2::OPERATOR_DIV);
        Kokkos::deep_copy(cell_div_basis_ref_host,cell_div_basis_ref.get_view());

        using fst=Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>;
        fst::HDIVtransformDIV(cell_div_basis.get_view(),cell_jac_det.get_view(),cell_div_basis_ref.get_view());
	Kokkos::deep_copy(cell_div_basis_host, cell_div_basis.get_static_view());


        // Copy cell quantity back into main array
        for(int b=0; b<num_card; ++b)
          for(int p=0; p<num_points; ++p)
            tmp_div_basis_host(cell,b,p) = cell_div_basis_host(0,b,p);
        Kokkos::deep_copy(tmp_div_basis.get_view(),tmp_div_basis_host);
      }
    }

    if(orientations_ != Teuchos::null)
      applyOrientationsImpl<Scalar>(num_orientations_cells_, tmp_div_basis.get_view(), *orientations_->getOrientations(), *intrepid_basis);

    // Store for later if cache is enabled
    PANZER_CACHE_DATA(div_basis);

    return tmp_div_basis;

  }

}

//=======================================================================================================

} // namespace panzer
