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

#include "Intrepid2_Utils.hpp"
#include "Intrepid2_FunctionSpaceTools.hpp"
#include "Intrepid2_Orientation.hpp"
#include "Intrepid2_OrientationTools.hpp"


namespace panzer {

template <typename Scalar>
void panzer::BasisValues2<Scalar>::
evaluateValues(const PHX::MDField<Scalar,IP,Dim,void,void,void,void,void,void> & cub_points,
               const PHX::MDField<Scalar,Cell,IP,Dim,Dim,void,void,void,void> & jac,
               const PHX::MDField<Scalar,Cell,IP,void,void,void,void,void,void> & jac_det,
               const PHX::MDField<Scalar,Cell,IP,Dim,Dim,void,void,void,void> & jac_inv,
               const int in_num_cells)
{
  PHX::MDField<Scalar,Cell,IP> weighted_measure;
  PHX::MDField<Scalar,Cell,NODE,Dim> vertex_coordinates;
  build_weighted = false;
  evaluateValues(cub_points,jac,jac_det,jac_inv,weighted_measure,vertex_coordinates,false,in_num_cells);
}

template <typename Scalar>
void panzer::BasisValues2<Scalar>::
evaluateValues(const PHX::MDField<Scalar,IP,Dim,void,void,void,void,void,void> & cub_points,
               const PHX::MDField<Scalar,Cell,IP,Dim,Dim,void,void,void,void> & jac,
               const PHX::MDField<Scalar,Cell,IP,void,void,void,void,void,void> & jac_det,
               const PHX::MDField<Scalar,Cell,IP,Dim,Dim,void,void,void,void> & jac_inv,
               const PHX::MDField<Scalar,Cell,IP> & weighted_measure,
               const PHX::MDField<Scalar,Cell,NODE,Dim> & vertex_coordinates,
               bool use_vertex_coordinates,
               const int in_num_cells)
{
  MDFieldArrayFactory af("",ddims_,true);

  int num_dim   = basis_layout->dimension();

  // currently this just copies on the basis objects are converted
  // in intrepid
  evaluateReferenceValues(cub_points,compute_derivatives,use_vertex_coordinates);

  const int num_cells = in_num_cells < 0 ? jac.extent(0) : in_num_cells;
  const std::pair<int,int> cell_range(0,num_cells);

  PureBasis::EElementSpace elmtspace = getElementSpace();
  if(elmtspace==PureBasis::CONST ||
     elmtspace==PureBasis::HGRAD) {
    auto s_basis_scalar = Kokkos::subview(basis_scalar.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL());
    Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>::
      HGRADtransformVALUE(s_basis_scalar,
                          basis_ref_scalar.get_view());

    if(build_weighted) {
      auto s_weighted_basis_scalar = Kokkos::subview(weighted_basis_scalar.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL());
      auto s_weighted_measure = Kokkos::subview(weighted_measure.get_view(),cell_range,Kokkos::ALL());
      Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>::
        multiplyMeasure(s_weighted_basis_scalar,
                        s_weighted_measure,
                        s_basis_scalar);
    }
  }
  else if(elmtspace==PureBasis::HCURL) {
    auto s_basis_vector = Kokkos::subview(basis_vector.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
    auto s_jac_inv = Kokkos::subview(jac_inv.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
    Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>::
      HCURLtransformVALUE(s_basis_vector,
                          s_jac_inv,
                          basis_ref_vector.get_view());

    if(build_weighted) {
      auto s_weighted_basis_vector = Kokkos::subview(weighted_basis_vector.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
      auto s_weighted_measure = Kokkos::subview(weighted_measure.get_view(),cell_range,Kokkos::ALL());
      Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>::
        multiplyMeasure(s_weighted_basis_vector,
                        s_weighted_measure,
                        s_basis_vector);
    }
  }
  else if(elmtspace==PureBasis::HDIV)
  {
    auto s_basis_vector = Kokkos::subview(basis_vector.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
    auto s_jac = Kokkos::subview(jac.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
    auto s_jac_det = Kokkos::subview(jac_det.get_view(),cell_range,Kokkos::ALL());
    Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>::
      HDIVtransformVALUE(s_basis_vector,
                         s_jac,
                         s_jac_det,
                         basis_ref_vector.get_view());

    if(build_weighted) {
      auto s_weighted_basis_vector = Kokkos::subview(weighted_basis_vector.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
      auto s_weighted_measure = Kokkos::subview(weighted_measure.get_view(),cell_range,Kokkos::ALL());
      Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>::
        multiplyMeasure(s_weighted_basis_vector,
                        s_weighted_measure,
                        s_basis_vector);
    }
  }
  else if(elmtspace==PureBasis::HVOL)
  {
    auto s_basis_scalar = Kokkos::subview(basis_scalar.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL());
    auto s_jac_det = Kokkos::subview(jac_det.get_view(),cell_range,Kokkos::ALL());
    Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>::
      HVOLtransformVALUE(s_basis_scalar,
                         s_jac_det,
                         basis_ref_scalar.get_view());

    if(build_weighted) {
      auto s_weighted_basis_scalar = Kokkos::subview(weighted_basis_scalar.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL());
      auto s_weighted_measure = Kokkos::subview(weighted_measure.get_view(),cell_range,Kokkos::ALL());
      Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>::
        multiplyMeasure(s_weighted_basis_scalar,
                        s_weighted_measure,
                        s_basis_scalar);
    }
  }
  else { TEUCHOS_ASSERT(false); }

  if(elmtspace==PureBasis::HGRAD && compute_derivatives) {
    auto s_grad_basis = Kokkos::subview(grad_basis.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
    auto s_jac_inv = Kokkos::subview(jac_inv.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
    Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>::
      HGRADtransformGRAD(s_grad_basis,
                         s_jac_inv,
                         grad_basis_ref.get_view());

    if(build_weighted) {
      auto s_weighted_grad_basis = Kokkos::subview(weighted_grad_basis.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
      auto s_weighted_measure = Kokkos::subview(weighted_measure.get_view(),cell_range,Kokkos::ALL());
      Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>::
        multiplyMeasure(s_weighted_grad_basis,
                        s_weighted_measure,
                        s_grad_basis);
    }
  }
  else if(elmtspace==PureBasis::HCURL && num_dim==2 && compute_derivatives) {
    auto s_curl_basis_scalar = Kokkos::subview(curl_basis_scalar.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL());
    auto s_jac_det = Kokkos::subview(jac_det.get_view(),cell_range,Kokkos::ALL());
    Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>::
      HDIVtransformDIV(s_curl_basis_scalar,
                       s_jac_det,   // note only volume deformation is needed!
                                    // this relates directly to this being in
                                    // the divergence space in 2D!
                       curl_basis_ref_scalar.get_view());

    if(build_weighted) {
      auto s_weighted_curl_basis_scalar = Kokkos::subview(weighted_curl_basis_scalar.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL());
      auto s_weighted_measure = Kokkos::subview(weighted_measure.get_view(),cell_range,Kokkos::ALL());
      Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>::
        multiplyMeasure(s_weighted_curl_basis_scalar,
                        s_weighted_measure,
                        s_curl_basis_scalar);
    }
  }
  else if(elmtspace==PureBasis::HCURL && num_dim==3 && compute_derivatives) {
    auto s_curl_basis_vector = Kokkos::subview(curl_basis_vector.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
    auto s_jac = Kokkos::subview(jac.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
    auto s_jac_det = Kokkos::subview(jac_det.get_view(),cell_range,Kokkos::ALL());
    Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>::
      HCURLtransformCURL(s_curl_basis_vector,
                         s_jac,
                         s_jac_det,
                         curl_basis_ref_vector.get_view());

    if(build_weighted) {
      auto s_weighted_curl_basis_vector = Kokkos::subview(weighted_curl_basis_vector.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
      auto s_weighted_measure = Kokkos::subview(weighted_measure.get_view(),cell_range,Kokkos::ALL());
      Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>::
        multiplyMeasure(s_weighted_curl_basis_vector,
                        s_weighted_measure,
                        s_curl_basis_vector);
    }
  }
  else if(elmtspace==PureBasis::HDIV && compute_derivatives) {
    auto s_div_basis = Kokkos::subview(div_basis.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL());
    auto s_jac_det = Kokkos::subview(jac_det.get_view(),cell_range,Kokkos::ALL());
    Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>::
      HDIVtransformDIV(s_div_basis,
                       s_jac_det,
                       div_basis_ref.get_view());

    if(build_weighted) {
      auto s_weighted_div_basis = Kokkos::subview(weighted_div_basis.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL());
      auto s_weighted_measure = Kokkos::subview(weighted_measure.get_view(),cell_range,Kokkos::ALL());
      Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>::
        multiplyMeasure(s_weighted_div_basis,
                        s_weighted_measure,
                        s_div_basis);
    }
  }

  // If basis supports coordinate values at basis points, then
  // compute these values
  if(use_vertex_coordinates) {
    // Teuchos::RCP<Intrepid2::DofCoordsInterface<ArrayDynamic> > coords
    //     = Teuchos::rcp_dynamic_cast<Intrepid2::DofCoordsInterface<ArrayDynamic> >(intrepid_basis);
    // if (!Teuchos::is_null(coords)) {
/*
      ArrayDynamic dyn_basis_coordinates_ref = af.buildArray<Scalar,BASIS,Dim>("basis_coordinates_ref",basis_coordinates_ref.extent(0),basis_coordinates_ref.extent(1));
      coords->getDofCoords(dyn_basis_coordinates_ref);

      // fill in basis coordinates
      for (size_type i = 0; i < basis_coordinates_ref.extent(0); ++i)
        for (size_type j = 0; j < basis_coordinates_ref.extent(1); ++j)
           basis_coordinates_ref(i,j) = dyn_basis_coordinates_ref(i,j);
*/

    auto s_basis_coordinates = Kokkos::subview(basis_coordinates.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL());
    auto s_vertex_coordinates = Kokkos::subview(vertex_coordinates.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL());
    Intrepid2::CellTools<PHX::Device::execution_space> cell_tools;
    cell_tools.mapToPhysicalFrame(s_basis_coordinates,
                                  basis_coordinates_ref.get_view(),
                                  s_vertex_coordinates,
                                  intrepid_basis->getBaseCellTopology());
  }
}







template <typename Scalar>
void panzer::BasisValues2<Scalar>::
evaluateBasisCoordinates(const PHX::MDField<Scalar,Cell,NODE,Dim> & vertex_coordinates,
                         const int in_num_cells)
{
  MDFieldArrayFactory af("",ddims_,true);

  // intrepid_basis->getDofCoords requires DynRankView, but basis_coordinates_ref is more of a static View
  // We use an auxiliary 'dyn' array to get around this
  using coordsScalarType = typename Intrepid2::Basis<PHX::Device::execution_space,Scalar,Scalar>::scalarType;
  auto dyn_basis_coordinates_ref = af.buildArray<coordsScalarType,BASIS,Dim>("basis_coordinates_ref",
                                                                             basis_coordinates_ref.extent(0),
                                                                             basis_coordinates_ref.extent(1));
  intrepid_basis->getDofCoords(dyn_basis_coordinates_ref.get_view());

  // fill in basis coordinates
  for (int i = 0; i < basis_coordinates_ref.extent_int(0); ++i)
    for (int j = 0; j < basis_coordinates_ref.extent_int(1); ++j)
      basis_coordinates_ref(i,j) = dyn_basis_coordinates_ref(i,j);

  const int num_cells = in_num_cells < 0 ? vertex_coordinates.extent(0) : in_num_cells;
  const std::pair<int,int> cell_range(0,num_cells);
  const auto s_basis_coordinates = Kokkos::subview(basis_coordinates.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL());
  const auto s_vertex_coordinates = Kokkos::subview(vertex_coordinates.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL());

  Intrepid2::CellTools<PHX::Device::execution_space> cell_tools;
  cell_tools.mapToPhysicalFrame(s_basis_coordinates,
                                basis_coordinates_ref.get_view(),
                                s_vertex_coordinates,
                                intrepid_basis->getBaseCellTopology());
}





template <typename Scalar>
void panzer::BasisValues2<Scalar>::
evaluateValues(const PHX::MDField<Scalar,Cell,IP,Dim,void,void,void,void,void> & cub_points,
               const PHX::MDField<Scalar,Cell,IP,Dim,Dim,void,void,void,void> & jac,
               const PHX::MDField<Scalar,Cell,IP,void,void,void,void,void,void> & jac_det,
               const PHX::MDField<Scalar,Cell,IP,Dim,Dim,void,void,void,void> & jac_inv,
               const PHX::MDField<Scalar,Cell,IP> & weighted_measure,
               const PHX::MDField<Scalar,Cell,NODE,Dim> & vertex_coordinates,
               bool use_vertex_coordinates,
               const int in_num_cells)
{

  PureBasis::EElementSpace elmtspace = getElementSpace();

  if(elmtspace == PureBasis::CONST){
    evaluateValues_Const(cub_points,jac_inv,weighted_measure,in_num_cells);
  } else if(elmtspace == PureBasis::HVOL){
    evaluateValues_HVol(cub_points,jac_det,jac_inv,weighted_measure,in_num_cells);
  } else if(elmtspace == PureBasis::HGRAD){
    evaluateValues_HGrad(cub_points,jac_inv,weighted_measure,in_num_cells);
  } else if(elmtspace == PureBasis::HCURL){
    evaluateValues_HCurl(cub_points,jac,jac_det,jac_inv,weighted_measure,in_num_cells);
  } else if(elmtspace == PureBasis::HDIV){
    evaluateValues_HDiv(cub_points,jac,jac_det,weighted_measure,in_num_cells);
  } else {
    TEUCHOS_TEST_FOR_EXCEPT_MSG(true,"panzer::BasisValues2::evaluateValues : Element space not recognized.");
  }

  if(use_vertex_coordinates) {
    TEUCHOS_TEST_FOR_EXCEPT_MSG(elmtspace == PureBasis::CONST,"panzer::BasisValues2::evaluateValues : Const basis cannot have basis coordinates.");
    evaluateBasisCoordinates(vertex_coordinates);
  }

}

template <typename Scalar>
void panzer::BasisValues2<Scalar>::
evaluateValues_Const(const PHX::MDField<Scalar,Cell,IP,Dim,void,void,void,void,void> & cub_points,
                     const PHX::MDField<Scalar,Cell,IP,Dim,Dim,void,void,void,void> & jac_inv,
                     const PHX::MDField<Scalar,Cell,IP> & weighted_measure,
                     const int in_num_cells)
{

  TEUCHOS_ASSERT(getElementSpace() == PureBasis::CONST);

  typedef Intrepid2::FunctionSpaceTools<PHX::Device::execution_space> fst;
  MDFieldArrayFactory af("",ddims_,true);

  const panzer::PureBasis & basis = *(basis_layout->getBasis());

  const int num_points = basis_layout->numPoints();
  const int num_basis  = basis.cardinality();
  const int num_dim    = basis_layout->dimension();
  const int num_cells = in_num_cells < 0 ? basis_layout->numCells() : in_num_cells;

  auto cell_basis_scalar = af.buildStaticArray<Scalar,Cell,BASIS,IP>("cell_basis_scalar",1,num_basis,num_points);
  auto cell_cub_points = af.buildStaticArray<Scalar,IP,Dim>("cell_cub_points",num_points,num_dim);
  auto cell_grad_basis = af.buildStaticArray<Scalar,Cell,BASIS,IP,Dim>("cell_grad_basis",1,num_basis,num_points,num_dim);
  auto cell_jac_inv = af.buildStaticArray<Scalar,Cell,IP,Dim,Dim>("cell_jac_inv",1,num_points,num_dim,num_dim);

  auto cell_basis_ref_scalar = af.buildStaticArray<Scalar,BASIS,IP>("cell_basis_ref_scalar",num_basis,num_points);
  auto cell_grad_basis_ref = af.buildStaticArray<Scalar,BASIS,IP,Dim>("cell_grad_basis_ref",num_basis,num_points,num_dim);

  for(int cell=0;cell<num_cells;++cell){

    // =============================================
    // Load external into cell-local arrays

    for(int p=0;p<num_points;++p)
      for(int d=0;d<num_dim;++d)
        for(int d2=0;d2<num_dim;++d2)
          cell_jac_inv(0,p,d,d2)=jac_inv(cell,p,d,d2);
    for(int p=0;p<num_points;++p)
      for(int d=0;d<num_dim;++d)
        cell_cub_points(p,d)=cub_points(cell,p,d);

    // =============================================
    // Load Reference Values

    intrepid_basis->getValues(cell_basis_ref_scalar.get_view(),cell_cub_points.get_view(),Intrepid2::OPERATOR_VALUE);

    if(compute_derivatives){
      Kokkos::deep_copy(cell_grad_basis_ref.get_view(),0.0);
    }

    // =============================================
    // Transform reference values to physical values

    fst::HGRADtransformVALUE(cell_basis_scalar.get_view(),cell_basis_ref_scalar.get_view());
    for(int b=0;b<num_basis;++b)
      for(int p=0;p<num_points;++p)
        basis_scalar(cell,b,p)=cell_basis_scalar(0,b,p);

    if(compute_derivatives){
        fst::HGRADtransformGRAD(cell_grad_basis.get_view(),cell_jac_inv.get_view(),cell_grad_basis_ref.get_view());
        for(int b=0;b<num_basis;++b)
          for(int p=0;p<num_points;++p)
            for(int d=0;d<num_dim;++d)
              grad_basis(cell,b,p,d)=cell_grad_basis(0,b,p,d);
    }
    // =============================================
  }


  if(build_weighted){
    const std::pair<int,int> cell_range(0,num_cells);
    auto s_weighted_basis_scalar = Kokkos::subview(weighted_basis_scalar.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL());
    auto s_weighted_measure = Kokkos::subview(weighted_measure.get_view(),cell_range,Kokkos::ALL());
    auto s_basis_scalar = Kokkos::subview(basis_scalar.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL());
    fst::multiplyMeasure(s_weighted_basis_scalar,s_weighted_measure,s_basis_scalar);
    if(compute_derivatives){
      auto s_weighted_grad_basis = Kokkos::subview(weighted_grad_basis.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
      auto s_grad_basis = Kokkos::subview(grad_basis.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
      fst::multiplyMeasure(s_weighted_grad_basis,s_weighted_measure,s_grad_basis);
    }
  }


}

template <typename Scalar>
void panzer::BasisValues2<Scalar>::
evaluateValues_HVol(const PHX::MDField<Scalar,Cell,IP,Dim,void,void,void,void,void> & cub_points,
                    const PHX::MDField<Scalar,Cell,IP,void,void,void,void,void,void> & jac_det,
                    const PHX::MDField<Scalar,Cell,IP,Dim,Dim,void,void,void,void> & jac_inv,
                    const PHX::MDField<Scalar,Cell,IP> & weighted_measure,
                    const int in_num_cells)
{

  TEUCHOS_ASSERT(getElementSpace() == PureBasis::HVOL);

  typedef Intrepid2::FunctionSpaceTools<PHX::Device::execution_space> fst;
  MDFieldArrayFactory af("",ddims_,true);

  const panzer::PureBasis & basis = *(basis_layout->getBasis());

  const int num_points = basis_layout->numPoints();
  const int num_basis  = basis.cardinality();
  const int num_dim    = basis_layout->dimension();
  const int num_cells = in_num_cells < 0 ? basis_layout->numCells() : in_num_cells;

  auto cell_basis_scalar = af.buildStaticArray<Scalar,Cell,BASIS,IP>("cell_basis_scalar",1,num_basis,num_points);
  auto cell_cub_points = af.buildStaticArray<Scalar,IP,Dim>("cell_cub_points",num_points,num_dim);
  auto cell_jac_det = af.buildStaticArray<Scalar,Cell,IP>("cell_jac_det",1,num_points);

  auto cell_basis_ref_scalar = af.buildStaticArray<Scalar,BASIS,IP>("cell_basis_ref_scalar",num_basis,num_points);

  for(int cell=0;cell<num_cells;++cell){

    // =============================================
    // Load external into cell-local arrays

    for(int p=0;p<num_points;++p)
      cell_jac_det(0,p)=jac_det(cell,p);
    for(int p=0;p<num_points;++p)
      for(int d=0;d<num_dim;++d)
        cell_cub_points(p,d)=cub_points(cell,p,d);

    // =============================================
    // Load Reference Values

    intrepid_basis->getValues(cell_basis_ref_scalar.get_view(),cell_cub_points.get_view(),Intrepid2::OPERATOR_VALUE);

    // =============================================
    // Transform reference values to physical values

    fst::HVOLtransformVALUE(cell_basis_scalar.get_view(),cell_jac_det.get_view(),cell_basis_ref_scalar.get_view());
    for(int b=0;b<num_basis;++b)
      for(int p=0;p<num_points;++p)
        basis_scalar(cell,b,p)=cell_basis_scalar(0,b,p);

    // =============================================
  }


  if(build_weighted) {
    const std::pair<int,int> cell_range(0,num_cells);
    auto s_weighted_basis_scalar = Kokkos::subview(weighted_basis_scalar.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL());
    auto s_weighted_measure = Kokkos::subview(weighted_measure.get_view(),cell_range,Kokkos::ALL());
    auto s_basis_scalar = Kokkos::subview(basis_scalar.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL());
    fst::multiplyMeasure(s_weighted_basis_scalar,s_weighted_measure,s_basis_scalar);
  }

}

template <typename Scalar>
void panzer::BasisValues2<Scalar>::
evaluateValues_HGrad(const PHX::MDField<Scalar,Cell,IP,Dim,void,void,void,void,void> & cub_points,
                     const PHX::MDField<Scalar,Cell,IP,Dim,Dim,void,void,void,void> & jac_inv,
                     const PHX::MDField<Scalar,Cell,IP> & weighted_measure,
                     const int in_num_cells)
{

  TEUCHOS_ASSERT(getElementSpace() == PureBasis::HGRAD);

  typedef Intrepid2::FunctionSpaceTools<PHX::Device::execution_space> fst;
  MDFieldArrayFactory af("",ddims_,true);

  const panzer::PureBasis & basis = *(basis_layout->getBasis());

  const int num_points = basis_layout->numPoints();
  const int num_basis  = basis.cardinality();
  const int num_dim    = basis_layout->dimension();
  const int num_cells = in_num_cells < 0 ? cub_points.extent(0) : in_num_cells;

  auto cell_basis_scalar = af.buildStaticArray<Scalar,Cell,BASIS,IP>("cell_basis_scalar",1,num_basis,num_points);
  auto cell_cub_points = af.buildStaticArray<Scalar,IP,Dim>("cell_cub_points",num_points,num_dim);
  auto cell_grad_basis = af.buildStaticArray<Scalar,Cell,BASIS,IP,Dim>("cell_grad_basis",1,num_basis,num_points,num_dim);
  auto cell_jac_inv = af.buildStaticArray<Scalar,Cell,IP,Dim,Dim>("cell_jac_inv",1,num_points,num_dim,num_dim);

  auto cell_basis_ref_scalar = af.buildStaticArray<Scalar,BASIS,IP>("cell_basis_ref_scalar",num_basis,num_points);
  auto cell_grad_basis_ref = af.buildStaticArray<Scalar,BASIS,IP,Dim>("cell_grad_basis_ref",num_basis,num_points,num_dim);

  for(int cell=0;cell<num_cells;++cell){

    // =============================================
    // Load external into cell-local arrays

    for(int p=0;p<num_points;++p)
      for(int d=0;d<num_dim;++d)
        for(int d2=0;d2<num_dim;++d2)
          cell_jac_inv(0,p,d,d2)=jac_inv(cell,p,d,d2);
    for(int p=0;p<num_points;++p)
      for(int d=0;d<num_dim;++d)
        cell_cub_points(p,d)=cub_points(cell,p,d);

    // =============================================
    // Load Reference Values

    intrepid_basis->getValues(cell_basis_ref_scalar.get_view(),cell_cub_points.get_view(),Intrepid2::OPERATOR_VALUE);

    if(compute_derivatives){
      intrepid_basis->getValues(cell_grad_basis_ref.get_view(),cell_cub_points.get_view(),Intrepid2::OPERATOR_GRAD);
    }

    // =============================================
    // Transform reference values to physical values

    fst::HGRADtransformVALUE(cell_basis_scalar.get_view(),cell_basis_ref_scalar.get_view());
    for(int b=0;b<num_basis;++b)
      for(int p=0;p<num_points;++p)
        basis_scalar(cell,b,p)=cell_basis_scalar(0,b,p);

    if(compute_derivatives){
        fst::HGRADtransformGRAD(cell_grad_basis.get_view(),cell_jac_inv.get_view(),cell_grad_basis_ref.get_view());
        for(int b=0;b<num_basis;++b)
          for(int p=0;p<num_points;++p)
            for(int d=0;d<num_dim;++d)
              grad_basis(cell,b,p,d)=cell_grad_basis(0,b,p,d);
    }
    // =============================================
  }

  if(build_weighted){
    const std::pair<int,int> cell_range(0,num_cells);
    auto s_weighted_basis_scalar = Kokkos::subview(weighted_basis_scalar.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL());
    auto s_weighted_measure = Kokkos::subview(weighted_measure.get_view(),cell_range,Kokkos::ALL());
    auto s_basis_scalar = Kokkos::subview(basis_scalar.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL());
    fst::multiplyMeasure(s_weighted_basis_scalar,s_weighted_measure,s_basis_scalar);
    if(compute_derivatives){
      auto s_weighted_grad_basis = Kokkos::subview(weighted_grad_basis.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
      auto s_grad_basis = Kokkos::subview(grad_basis.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
      fst::multiplyMeasure(s_weighted_grad_basis,s_weighted_measure,s_grad_basis);
    }
  }

}


template <typename Scalar>
void panzer::BasisValues2<Scalar>::
evaluateValues_HCurl(const PHX::MDField<Scalar,Cell,IP,Dim,void,void,void,void,void> & cub_points,
               const PHX::MDField<Scalar,Cell,IP,Dim,Dim,void,void,void,void> & jac,
               const PHX::MDField<Scalar,Cell,IP,void,void,void,void,void,void> & jac_det,
               const PHX::MDField<Scalar,Cell,IP,Dim,Dim,void,void,void,void> & jac_inv,
                     const PHX::MDField<Scalar,Cell,IP> & weighted_measure,
                     const int in_num_cells)
{

  TEUCHOS_ASSERT(getElementSpace() == PureBasis::HCURL);


  typedef Intrepid2::FunctionSpaceTools<PHX::Device::execution_space> fst;
  MDFieldArrayFactory af("",ddims_,true);

  const panzer::PureBasis & basis = *(basis_layout->getBasis());

  const int num_points = basis_layout->numPoints();
  const int num_basis  = basis.cardinality();
  const int num_dim    = basis_layout->dimension();
  const int num_cells = in_num_cells < 0 ? basis_layout->numCells() : in_num_cells;

  auto cell_cub_points = af.buildStaticArray<Scalar,IP,Dim>("cell_cub_points",num_points,num_dim);
  auto cell_jac = af.buildStaticArray<Scalar,Cell,IP,Dim,Dim>("cell_jac",1,num_points,num_dim,num_dim);
  auto cell_jac_inv = af.buildStaticArray<Scalar,Cell,IP,Dim,Dim>("cell_jac_inv",1,num_points,num_dim,num_dim);
  auto cell_jac_det = af.buildStaticArray<Scalar,Cell,IP>("cell_jac_det",1,num_points);

  auto cell_basis_vector = af.buildStaticArray<Scalar,Cell,BASIS,IP,Dim>("cell_basis_vector",1,num_basis,num_points,num_dim);
  auto cell_curl_basis_scalar = af.buildStaticArray<Scalar,Cell,BASIS,IP>("cell_curl_basis_scalar",1,num_basis,num_points);
  auto cell_curl_basis_vector = af.buildStaticArray<Scalar,Cell,BASIS,IP,Dim>("cell_curl_basis_vector",1,num_basis,num_points,num_dim);

  auto cell_curl_basis_ref = af.buildArray<Scalar,BASIS,IP,Dim>("cell_curl_basis_ref",num_basis,num_points,num_dim);
  auto cell_curl_basis_ref_scalar =  af.buildStaticArray<Scalar,BASIS,IP>("cell_curl_basis_ref_scalar",num_basis,num_points);
  auto cell_basis_ref_vector = af.buildArray<Scalar,BASIS,IP,Dim>("cell_basis_ref_vector",num_basis,num_points,num_dim);

  for(int cell=0;cell<num_cells;++cell){

    // =============================================
    // Load external into cell-local arrays

    for(int p=0;p<num_points;++p)
      for(int d=0;d<num_dim;++d)
        for(int d2=0;d2<num_dim;++d2)
          cell_jac(0,p,d,d2)=jac(cell,p,d,d2);
    for(int p=0;p<num_points;++p)
      for(int d=0;d<num_dim;++d)
        for(int d2=0;d2<num_dim;++d2)
          cell_jac_inv(0,p,d,d2)=jac_inv(cell,p,d,d2);
    for(int p=0;p<num_points;++p)
      cell_jac_det(0,p)=jac_det(cell,p);
    for(int p=0;p<num_points;++p)
      for(int d=0;d<num_dim;++d)
        cell_cub_points(p,d)=cub_points(cell,p,d);

    // =============================================
    // Load Reference Values

    intrepid_basis->getValues(cell_basis_ref_vector.get_view(),cell_cub_points.get_view(),Intrepid2::OPERATOR_VALUE);

    if(compute_derivatives){
      if(num_dim==2){
        intrepid_basis->getValues(cell_curl_basis_ref_scalar.get_view(),cell_cub_points.get_view(),Intrepid2::OPERATOR_CURL);
      } else if(num_dim==3){
        intrepid_basis->getValues(cell_curl_basis_ref.get_view(),cell_cub_points.get_view(),Intrepid2::OPERATOR_CURL);
      }
    }

    // =============================================
    // Transform reference values to physical values

    fst::HCURLtransformVALUE(cell_basis_vector.get_view(),cell_jac_inv.get_view(),cell_basis_ref_vector.get_view());
    for(int b=0;b<num_basis;++b)
      for(int p=0;p<num_points;++p)
        for(int d=0;d<num_dim;++d)
          basis_vector(cell,b,p,d)=cell_basis_vector(0,b,p,d);

    if(compute_derivatives){
      if(num_dim==2){
        // note only volume deformation is needed!
        // this relates directly to this being in
        // the divergence space in 2D!
        fst::HDIVtransformDIV(cell_curl_basis_scalar.get_view(),cell_jac_det.get_view(),cell_curl_basis_ref_scalar.get_view());
        for(int b=0;b<num_basis;++b)
          for(int p=0;p<num_points;++p)
            curl_basis_scalar(cell,b,p)=cell_curl_basis_scalar(0,b,p);
      } else if(num_dim==3) {
        fst::HCURLtransformCURL(cell_curl_basis_vector.get_view(),cell_jac.get_view(),cell_jac_det.get_view(),cell_curl_basis_ref.get_view());
        for(int b=0;b<num_basis;++b)
          for(int p=0;p<num_points;++p)
            for(int d=0;d<num_dim;++d)
              curl_basis_vector(cell,b,p,d)=cell_curl_basis_vector(0,b,p,d);
      } else {
        TEUCHOS_TEST_FOR_EXCEPT_MSG(true,"panzer::BasisValues2::evaluateValues_HCurl : HCurl only setup for 2D and 3D.");
      }
    }
  }

  if(build_weighted){
    const std::pair<int,int> cell_range(0,num_cells);
    auto s_weighted_basis_vector = Kokkos::subview(weighted_basis_vector.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
    auto s_weighted_measure = Kokkos::subview(weighted_measure.get_view(),cell_range,Kokkos::ALL());
    auto s_basis_vector = Kokkos::subview(basis_vector.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
    fst::multiplyMeasure(s_weighted_basis_vector,s_weighted_measure,s_basis_vector);
    if(compute_derivatives){
      if(num_dim==2){
        auto s_weighted_curl_basis_scalar = Kokkos::subview(weighted_curl_basis_scalar.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL());
        auto s_curl_basis_scalar = Kokkos::subview(curl_basis_scalar.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL());
        fst::multiplyMeasure(s_weighted_curl_basis_scalar,s_weighted_measure,s_curl_basis_scalar);
      } else if(num_dim==3){
        auto s_weighted_curl_basis_vector = Kokkos::subview(weighted_curl_basis_vector.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
        auto s_curl_basis_vector = Kokkos::subview(curl_basis_vector.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
        fst::multiplyMeasure(s_weighted_curl_basis_vector,s_weighted_measure,s_curl_basis_vector);
      }
    }
  }

}

template <typename Scalar>
void panzer::BasisValues2<Scalar>::
evaluateValues_HDiv(const PHX::MDField<Scalar,Cell,IP,Dim,void,void,void,void,void> & cub_points,
                    const PHX::MDField<Scalar,Cell,IP,Dim,Dim,void,void,void,void> & jac,
                    const PHX::MDField<Scalar,Cell,IP,void,void,void,void,void,void> & jac_det,
                    const PHX::MDField<Scalar,Cell,IP> & weighted_measure,
                    const int in_num_cells)
{

  TEUCHOS_ASSERT(getElementSpace() == PureBasis::HDIV);

  typedef Intrepid2::FunctionSpaceTools<PHX::Device::execution_space> fst;
  MDFieldArrayFactory af("",ddims_,true);

  const panzer::PureBasis & basis = *(basis_layout->getBasis());

  const int num_points = basis_layout->numPoints();
  const int num_basis  = basis.cardinality();
  const int num_dim    = basis_layout->dimension();
  const int num_cells = in_num_cells < 0 ? basis_layout->numCells() : in_num_cells;

  auto cell_cub_points = af.buildStaticArray<Scalar,IP,Dim>("cell_cub_points",num_points,num_dim);
  auto cell_jac = af.buildStaticArray<Scalar,Cell,IP,Dim,Dim>("cell_jac",1,num_points,num_dim,num_dim);
  auto cell_jac_det = af.buildStaticArray<Scalar,Cell,IP>("cell_jac_det",1,num_points);

  auto cell_basis_vector = af.buildStaticArray<Scalar,Cell,BASIS,IP,Dim>("cell_basis_vector",1,num_basis,num_points,num_dim);
  auto cell_div_basis = af.buildStaticArray<Scalar,Cell,BASIS,IP>("cell_div_basis",1,num_basis,num_points);

  auto cell_basis_ref_vector = af.buildArray<Scalar,BASIS,IP,Dim>("cell_basis_ref_vector",num_basis,num_points,num_dim);
  auto cell_div_basis_ref =  af.buildStaticArray<Scalar,BASIS,IP>("cell_div_basis_ref",num_basis,num_points);

  for(int cell=0;cell<num_cells;++cell){

    // =============================================
    // Load external into cell-local arrays

    for(int p=0;p<num_points;++p)
      for(int d=0;d<num_dim;++d)
        for(int d2=0;d2<num_dim;++d2)
          cell_jac(0,p,d,d2)=jac(cell,p,d,d2);
    for(int p=0;p<num_points;++p)
      cell_jac_det(0,p)=jac_det(cell,p);
    for(int p=0;p<num_points;++p)
      for(int d=0;d<num_dim;++d)
        cell_cub_points(p,d)=cub_points(cell,p,d);
    // =============================================
    // Load Reference Values

    intrepid_basis->getValues(cell_basis_ref_vector.get_view(),cell_cub_points.get_view(),Intrepid2::OPERATOR_VALUE);

    if(compute_derivatives){
      intrepid_basis->getValues(cell_div_basis_ref.get_view(),cell_cub_points.get_view(),Intrepid2::OPERATOR_DIV);
    }

    // =============================================
    // Transform reference values to physical values

    fst::HDIVtransformVALUE(cell_basis_vector.get_view(),cell_jac.get_view(),cell_jac_det.get_view(),cell_basis_ref_vector.get_view());
    for(int b=0;b<num_basis;++b)
      for(int p=0;p<num_points;++p)
        for(int d=0;d<num_dim;++d)
          basis_vector(cell,b,p,d)=cell_basis_vector(0,b,p,d);

    if(compute_derivatives){
      fst::HDIVtransformDIV(cell_div_basis.get_view(),cell_jac_det.get_view(),cell_div_basis_ref.get_view());
      for(int b=0;b<num_basis;++b)
        for(int p=0;p<num_points;++p)
          div_basis(cell,b,p)=cell_div_basis(0,b,p);
    }
  }

  if(build_weighted){
    const std::pair<int,int> cell_range(0,num_cells);
    auto s_weighted_basis_vector = Kokkos::subview(weighted_basis_vector.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
    auto s_weighted_measure = Kokkos::subview(weighted_measure.get_view(),cell_range,Kokkos::ALL());
    auto s_basis_vector = Kokkos::subview(basis_vector.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL(),Kokkos::ALL());
    fst::multiplyMeasure(s_weighted_basis_vector,s_weighted_measure,s_basis_vector);
    if(compute_derivatives){
      auto s_weighted_div_basis = Kokkos::subview(weighted_div_basis.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL());
      auto s_div_basis = Kokkos::subview(div_basis.get_view(),cell_range,Kokkos::ALL(),Kokkos::ALL());
      fst::multiplyMeasure(s_weighted_div_basis,s_weighted_measure,s_div_basis);
    }
  }

}











































template <typename Scalar>
void panzer::BasisValues2<Scalar>::
evaluateValuesCV(const PHX::MDField<Scalar,Cell,IP,Dim,void,void,void,void,void> & cell_cub_points,
                 const PHX::MDField<Scalar,Cell,IP,Dim,Dim,void,void,void,void> & jac,
                 const PHX::MDField<Scalar,Cell,IP,void,void,void,void,void,void> & jac_det,
                 const PHX::MDField<Scalar,Cell,IP,Dim,Dim,void,void,void,void> & jac_inv)
{
  MDFieldArrayFactory af("",ddims_,true);

  int num_ip    = basis_layout->numPoints();
  int num_card  = basis_layout->cardinality();
  int num_dim   = basis_layout->dimension();

  size_type num_cells = jac.extent(0);

  PureBasis::EElementSpace elmtspace = getElementSpace();
  ArrayDynamic dyn_cub_points = af.buildArray<Scalar,IP,Dim>("dyn_cub_points", num_ip, num_dim);

  // Integration points are located on physical cells rather than reference cells,
  // so we evaluate the basis in a loop over cells.
  for (size_type icell = 0; icell < num_cells; ++icell)
  {
    for (int ip = 0; ip < num_ip; ++ip)
      for (int d = 0; d < num_dim; ++d)
         dyn_cub_points(ip,d) = cell_cub_points(icell,ip,d);

    if(elmtspace==PureBasis::CONST) {
       ArrayDynamic dyn_basis_ref_scalar = af.buildArray<Scalar,BASIS,IP>("dyn_basis_ref_scalar",num_card,num_ip);

       intrepid_basis->getValues(dyn_basis_ref_scalar.get_view(),
                                 dyn_cub_points.get_view(),
                                 Intrepid2::OPERATOR_VALUE);

       // transform values method just transfers values to array with cell index - no need to call
       for (int b = 0; b < num_card; ++b)
         for (int ip = 0; ip < num_ip; ++ip)
           basis_scalar(icell,b,ip) = dyn_basis_ref_scalar(b,ip);

    }
    if(elmtspace==PureBasis::HVOL) {
      ArrayDynamic dyn_basis_ref_scalar = af.buildArray<Scalar,BASIS,IP>("dyn_basis_ref_scalar",num_card,num_ip);

      intrepid_basis->getValues(dyn_basis_ref_scalar.get_view(),
                                dyn_cub_points.get_view(),
                                Intrepid2::OPERATOR_VALUE);

      int one_cell= 1;
      ArrayDynamic dyn_basis_scalar = af.buildArray<Scalar,Cell,BASIS,IP>("dyn_basis_vector",one_cell,num_card,num_ip);
      ArrayDynamic dyn_jac = af.buildArray<Scalar,Cell,IP,Dim,Dim>("dyn_jac",one_cell,num_ip,num_dim,num_dim);
      ArrayDynamic dyn_jac_det = af.buildArray<Scalar,Cell,IP>("dyn_jac_det",one_cell,num_ip);

      int cellInd = 0;
      for (int ip = 0; ip < num_ip; ++ip)
        dyn_jac_det(cellInd,ip) = jac_det(icell,ip);

      Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>::HVOLtransformVALUE(dyn_basis_scalar.get_view(),
                                                                                      dyn_jac_det.get_view(),
                                                                                      dyn_basis_ref_scalar.get_view());

       for (int b = 0; b < num_card; ++b)
         for (int ip = 0; ip < num_ip; ++ip)
           basis_scalar(icell,b,ip) = dyn_basis_scalar(0,b,ip);

    }
    if(elmtspace==PureBasis::HGRAD) {
       ArrayDynamic dyn_basis_ref_scalar = af.buildArray<Scalar,BASIS,IP>("dyn_basis_ref_scalar",num_card,num_ip);

       intrepid_basis->getValues(dyn_basis_ref_scalar.get_view(),
                                 dyn_cub_points.get_view(),
                                 Intrepid2::OPERATOR_VALUE);

       // transform values method just transfers values to array with cell index - no need to call
       for (int b = 0; b < num_card; ++b)
         for (int ip = 0; ip < num_ip; ++ip)
           basis_scalar(icell,b,ip) = dyn_basis_ref_scalar(b,ip);

       if(compute_derivatives) {

          int one_cell = 1;
          ArrayDynamic dyn_grad_basis_ref = af.buildArray<Scalar,BASIS,IP,Dim>("dyn_grad_basis_ref",num_card,num_ip,num_dim);
          ArrayDynamic dyn_grad_basis = af.buildArray<Scalar,Cell,BASIS,IP,Dim>("dyn_grad_basis",one_cell,num_card,num_ip,num_dim);
          ArrayDynamic dyn_jac_inv = af.buildArray<Scalar,Cell,IP,Dim,Dim>("dyn_jac_inv",one_cell,num_ip,num_dim,num_dim);

          intrepid_basis->getValues(dyn_grad_basis_ref.get_view(),
                                    dyn_cub_points.get_view(),
                                    Intrepid2::OPERATOR_GRAD);

          int cellInd = 0;
          for (int ip = 0; ip < num_ip; ++ip)
             for (int d1 = 0; d1 < num_dim; ++d1)
               for (int d2 = 0; d2 < num_dim; ++d2)
                  dyn_jac_inv(cellInd,ip,d1,d2) = jac_inv(icell,ip,d1,d2);

          Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>::HGRADtransformGRAD<Scalar>(dyn_grad_basis.get_view(),
                                                                                                  dyn_jac_inv.get_view(),
                                                                                                  dyn_grad_basis_ref.get_view());

          for (int b = 0; b < num_card; ++b)
            for (int ip = 0; ip < num_ip; ++ip)
              for (int d = 0; d < num_dim; ++d)
                 grad_basis(icell,b,ip,d) = dyn_grad_basis(0,b,ip,d);

        }
    }
    else if(elmtspace==PureBasis::HCURL) {
      ArrayDynamic dyn_basis_ref_vector = af.buildArray<Scalar,BASIS,IP,Dim>("dyn_basis_ref_vector",num_card,num_ip,num_dim);

      intrepid_basis->getValues(dyn_basis_ref_vector.get_view(),
                                dyn_cub_points.get_view(),
                                Intrepid2::OPERATOR_VALUE);

      const int one_cell = 1;
      ArrayDynamic dyn_basis_vector = af.buildArray<Scalar,Cell,BASIS,IP,Dim>("dyn_basis_vector",one_cell,num_card,num_ip,num_dim);
      ArrayDynamic dyn_jac_inv = af.buildArray<Scalar,Cell,IP,Dim,Dim>("dyn_jac_inv",one_cell,num_ip,num_dim,num_dim);

      const int cellInd = 0;
      for (int ip = 0; ip < num_ip; ++ip)
        for (int d1 = 0; d1 < num_dim; ++d1)
          for (int d2 = 0; d2 < num_dim; ++d2)
              dyn_jac_inv(cellInd,ip,d1,d2) = jac_inv(icell,ip,d1,d2);

      Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>::HCURLtransformVALUE(dyn_basis_vector.get_view(),
                                                                                       dyn_jac_inv.get_view(),
                                                                                       dyn_basis_ref_vector.get_view());

      for (int b = 0; b < num_card; ++b)
        for (int ip = 0; ip < num_ip; ++ip)
          for (int d = 0; d < num_dim; ++d)
             basis_vector(icell,b,ip,d) = dyn_basis_vector(0,b,ip,d);

      if(compute_derivatives && num_dim ==2) {

          ArrayDynamic dyn_curl_basis_ref_scalar = af.buildArray<Scalar,BASIS,IP>("dyn_curl_basis_ref_scalar",num_card,num_ip);
          ArrayDynamic dyn_curl_basis_scalar = af.buildArray<Scalar,Cell,BASIS,IP>("dyn_curl_basis_scalar",one_cell,num_card,num_ip);
          ArrayDynamic dyn_jac_det = af.buildArray<Scalar,Cell,IP>("dyn_jac_det",one_cell,num_ip);

          intrepid_basis->getValues(dyn_curl_basis_ref_scalar.get_view(),
                                    dyn_cub_points.get_view(),
                                    Intrepid2::OPERATOR_CURL);

          for (int ip = 0; ip < num_ip; ++ip)
              dyn_jac_det(cellInd,ip) = jac_det(icell,ip);

          Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>::HDIVtransformDIV(dyn_curl_basis_scalar.get_view(),
                                                                                        dyn_jac_det.get_view(),
                                                                                        dyn_curl_basis_ref_scalar.get_view());

          for (int b = 0; b < num_card; ++b)
            for (int ip = 0; ip < num_ip; ++ip)
                curl_basis_scalar(icell,b,ip) = dyn_curl_basis_scalar(0,b,ip);

      }
      if(compute_derivatives && num_dim ==3) {

          ArrayDynamic dyn_curl_basis_ref = af.buildArray<Scalar,BASIS,IP,Dim>("dyn_curl_basis_ref_vector",num_card,num_ip,num_dim);
          ArrayDynamic dyn_curl_basis = af.buildArray<Scalar,Cell,BASIS,IP,Dim>("dyn_curl_basis_vector",one_cell,num_card,num_ip,num_dim);
          ArrayDynamic dyn_jac_det = af.buildArray<Scalar,Cell,IP>("dyn_jac_det",one_cell,num_ip);
          ArrayDynamic dyn_jac = af.buildArray<Scalar,Cell,IP,Dim,Dim>("dyn_jac",one_cell,num_ip,num_dim,num_dim);

          intrepid_basis->getValues(dyn_curl_basis_ref.get_view(),
                                    dyn_cub_points.get_view(),
                                    Intrepid2::OPERATOR_CURL);

          for (int ip = 0; ip < num_ip; ++ip)
          {
             dyn_jac_det(cellInd,ip) = jac_det(icell,ip);
             for (int d1 = 0; d1 < num_dim; ++d1)
                for (int d2 = 0; d2 < num_dim; ++d2)
                  dyn_jac(cellInd,ip,d1,d2) = jac(icell,ip,d1,d2);
          }

          Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>::HCURLtransformCURL(dyn_curl_basis.get_view(),
                                                                                          dyn_jac.get_view(),
                                                                                          dyn_jac_det.get_view(),
                                                                                          dyn_curl_basis_ref.get_view());

          for (int b = 0; b < num_card; ++b)
            for (int ip = 0; ip < num_ip; ++ip)
               for (int d = 0; d < num_dim; ++d)
                  curl_basis_vector(icell,b,ip,d) = dyn_curl_basis(0,b,ip,d);

      }

    }
    else if(elmtspace==PureBasis::HDIV) {

      ArrayDynamic dyn_basis_ref_vector = af.buildArray<Scalar,BASIS,IP,Dim>("dyn_basis_ref_vector",num_card,num_ip,num_dim);

      intrepid_basis->getValues(dyn_basis_ref_vector.get_view(),
                                dyn_cub_points.get_view(),
                                Intrepid2::OPERATOR_VALUE);

      int one_cell= 1;
      ArrayDynamic dyn_basis_vector = af.buildArray<Scalar,Cell,BASIS,IP,Dim>("dyn_basis_vector",one_cell,num_card,num_ip,num_dim);
      ArrayDynamic dyn_jac = af.buildArray<Scalar,Cell,IP,Dim,Dim>("dyn_jac",one_cell,num_ip,num_dim,num_dim);
      ArrayDynamic dyn_jac_det = af.buildArray<Scalar,Cell,IP>("dyn_jac_det",one_cell,num_ip);

      int cellInd = 0;
      for (int ip = 0; ip < num_ip; ++ip)
      {
        dyn_jac_det(cellInd,ip) = jac_det(icell,ip);
        for (int d1 = 0; d1 < num_dim; ++d1)
          for (int d2 = 0; d2 < num_dim; ++d2)
              dyn_jac(cellInd,ip,d1,d2) = jac(icell,ip,d1,d2);
      }

      Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>::HDIVtransformVALUE(dyn_basis_vector.get_view(),
                                                                                      dyn_jac.get_view(),
                                                                                      dyn_jac_det.get_view(),
                                                                                      dyn_basis_ref_vector.get_view());

       for (int b = 0; b < num_card; ++b)
         for (int ip = 0; ip < num_ip; ++ip)
           for (int d = 0; d < num_dim; ++d)
              basis_vector(icell,b,ip,d) = dyn_basis_vector(0,b,ip,d);

       if(compute_derivatives) {

           ArrayDynamic dyn_div_basis_ref = af.buildArray<Scalar,BASIS,IP>("dyn_div_basis_ref_scalar",num_card,num_ip);
           ArrayDynamic dyn_div_basis = af.buildArray<Scalar,Cell,BASIS,IP>("dyn_div_basis_scalar",one_cell,num_card,num_ip);

           intrepid_basis->getValues(dyn_div_basis_ref.get_view(),
                                     dyn_cub_points.get_view(),
                                     Intrepid2::OPERATOR_DIV);

           Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>::HDIVtransformDIV<Scalar>(dyn_div_basis.get_view(),
                                                                                                 dyn_jac_det.get_view(),
                                                                                                 dyn_div_basis_ref.get_view());

           for (int b = 0; b < num_card; ++b)
             for (int ip = 0; ip < num_ip; ++ip)
                 div_basis(icell,b,ip) = dyn_div_basis(0,b,ip);

        }

    }
    else { TEUCHOS_ASSERT(false); }

  } // cell loop

}

template <typename Scalar>
void panzer::BasisValues2<Scalar>::
evaluateReferenceValues(const PHX::MDField<Scalar,IP,Dim> & cub_points,bool in_compute_derivatives,bool use_vertex_coordinates)
{
  MDFieldArrayFactory af("",ddims_,true);

  int num_quad    = basis_layout->numPoints();
  int num_dim   = basis_layout->dimension();
  int num_card  = basis_layout->cardinality();

  ArrayDynamic dyn_cub_points = af.buildArray<Scalar,IP,Dim>("dyn_cub_points",  num_quad,num_dim);

  for (int ip = 0; ip < num_quad; ++ip)
    for (int d = 0; d < num_dim; ++d)
      dyn_cub_points(ip,d) = cub_points(ip,d);

  PureBasis::EElementSpace elmtspace = getElementSpace();
  if(elmtspace==PureBasis::HGRAD || elmtspace==PureBasis::CONST || elmtspace==PureBasis::HVOL) {
    ArrayDynamic dyn_basis_ref_scalar = af.buildArray<Scalar,BASIS,IP>("dyn_basis_ref_scalar",num_card,num_quad);

    intrepid_basis->getValues(dyn_basis_ref_scalar.get_view(),
                              dyn_cub_points.get_view(),
                              Intrepid2::OPERATOR_VALUE);

    for (int b = 0; b < num_card; ++b)
      for (int ip = 0; ip < num_quad; ++ip)
        basis_ref_scalar(b,ip) = dyn_basis_ref_scalar(b,ip);
  }
  else if(elmtspace==PureBasis::HDIV || elmtspace==PureBasis::HCURL) {
    ArrayDynamic dyn_basis_ref_vector = af.buildArray<Scalar,BASIS,IP,Dim>("dyn_basis_ref_vector",num_card,num_quad,num_dim);

    intrepid_basis->getValues(dyn_basis_ref_vector.get_view(),
                              dyn_cub_points.get_view(),
                              Intrepid2::OPERATOR_VALUE);

    for (int b = 0; b < num_card; ++b)
      for (int ip = 0; ip < num_quad; ++ip)
        for (int d = 0; d < num_dim; ++d)
           basis_ref_vector(b,ip,d) = dyn_basis_ref_vector(b,ip,d);
  }
  else { TEUCHOS_ASSERT(false); }

  if(elmtspace==PureBasis::HGRAD && in_compute_derivatives) {
    ArrayDynamic dyn_grad_basis_ref = af.buildArray<Scalar,BASIS,IP,Dim>("dyn_basis_ref_vector",num_card,num_quad,num_dim);

    intrepid_basis->getValues(dyn_grad_basis_ref.get_view(),
                              dyn_cub_points.get_view(),
                              Intrepid2::OPERATOR_GRAD);

    for (int b = 0; b < num_card; ++b)
      for (int ip = 0; ip < num_quad; ++ip)
        for (int d = 0; d < num_dim; ++d)
           grad_basis_ref(b,ip,d) = dyn_grad_basis_ref(b,ip,d);
  }
  else if(elmtspace==PureBasis::HCURL && in_compute_derivatives && num_dim==2) {
    ArrayDynamic dyn_curl_basis_ref = af.buildArray<Scalar,BASIS,IP>("dyn_curl_basis_ref_scalar",num_card,num_quad);

    intrepid_basis->getValues(dyn_curl_basis_ref.get_view(),
                              dyn_cub_points.get_view(),
                              Intrepid2::OPERATOR_CURL);

    for (int b = 0; b < num_card; ++b)
      for (int ip = 0; ip < num_quad; ++ip)
        curl_basis_ref_scalar(b,ip) = dyn_curl_basis_ref(b,ip);
  }
  else if(elmtspace==PureBasis::HCURL && in_compute_derivatives && num_dim==3) {
    ArrayDynamic dyn_curl_basis_ref = af.buildArray<Scalar,BASIS,IP,Dim>("dyn_curl_basis_ref_vector",num_card,num_quad,num_dim);

    intrepid_basis->getValues(dyn_curl_basis_ref.get_view(),
                              dyn_cub_points.get_view(),
                              Intrepid2::OPERATOR_CURL);

    for (int b = 0; b < num_card; ++b)
      for (int ip = 0; ip < num_quad; ++ip)
        for (int d = 0; d < num_dim; ++d)
           curl_basis_ref_vector(b,ip,d) = dyn_curl_basis_ref(b,ip,d);
  }
  else if(elmtspace==PureBasis::HDIV && in_compute_derivatives) {
    ArrayDynamic dyn_div_basis_ref = af.buildArray<Scalar,BASIS,IP>("dyn_div_basis_ref_scalar",num_card,num_quad);

    intrepid_basis->getValues(dyn_div_basis_ref.get_view(),
                              dyn_cub_points.get_view(),
                              Intrepid2::OPERATOR_DIV);

    for (int b = 0; b < num_card; ++b)
      for (int ip = 0; ip < num_quad; ++ip)
        div_basis_ref(b,ip) = dyn_div_basis_ref(b,ip);
  }


  if(use_vertex_coordinates) {
    // Intrepid removes fad types from the coordinate scalar type. We
    // pull the actual field scalar type from the basis object to be
    // consistent.
    if (elmtspace != PureBasis::CONST) {
      using coordsScalarType = typename Intrepid2::Basis<PHX::Device::execution_space,Scalar,Scalar>::scalarType;
      auto dyn_basis_coordinates_ref = af.buildArray<coordsScalarType,BASIS,Dim>("basis_coordinates_ref",
                                                                                 basis_coordinates_ref.extent(0),
                                                                                 basis_coordinates_ref.extent(1));
      intrepid_basis->getDofCoords(dyn_basis_coordinates_ref.get_view());

      // fill in basis coordinates
      for (int i = 0; i < basis_coordinates_ref.extent_int(0); ++i)
        for (int j = 0; j < basis_coordinates_ref.extent_int(1); ++j)
          basis_coordinates_ref(i,j) = dyn_basis_coordinates_ref(i,j);
    }
  }

  references_evaluated = true;
}

// method for applying orientations
template <typename Scalar>
void BasisValues2<Scalar>::
applyOrientations(const std::vector<Intrepid2::Orientation> & orientations,
                  const int in_num_cells)
{
  if (!intrepid_basis->requireOrientation())
    return;

  typedef Intrepid2::OrientationTools<PHX::Device> ots;
  const PureBasis::EElementSpace elmtspace = getElementSpace();

  // orientation (right now std vector) is created using push_back method.
  // thus, its size is the actual number of elements to be applied.
  // on the other hand, basis_layout num cell indicates workset size.
  // to get right size of cells, use minimum of them.
  const int num_cell_basis_layout = in_num_cells < 0 ? basis_layout->numCells() : in_num_cells;
  const int num_cell_orientation = orientations.size();
  const int num_cell  = num_cell_basis_layout < num_cell_orientation ? num_cell_basis_layout : num_cell_orientation;
  const int num_dim   = basis_layout->dimension();
  const Kokkos::pair<int,int> range_cell(0, num_cell);

  // Kokkos::DynRankView<Intrepid2::Orientation,PHX::Device>
  //   drv_orts((Intrepid2::Orientation*)orientations.data(), num_cell);
  Kokkos::DynRankView<Intrepid2::Orientation,PHX::Device> drv_orts("drv_orts", num_cell);
  auto host_drv_orts = Kokkos::create_mirror_view(drv_orts);
  for (size_t i=0; i < drv_orts.size(); ++i)
    host_drv_orts(i) = orientations[i];
  Kokkos::deep_copy(drv_orts,host_drv_orts);
  PHX::Device::fence();

  ///
  /// HGRAD elements
  ///
  if (elmtspace==PureBasis::HGRAD) {
    {
      {
        auto drv_basis_scalar = Kokkos::subview(basis_scalar.get_view(), range_cell, Kokkos::ALL(), Kokkos::ALL());
        auto drv_basis_scalar_tmp = Kokkos::createDynRankView(basis_scalar.get_view(),
                                                              "drv_basis_scalar_tmp",
                                                              drv_basis_scalar.extent(0),  // C
                                                              drv_basis_scalar.extent(1),  // F
                                                              drv_basis_scalar.extent(2)); // P
        Kokkos::deep_copy(drv_basis_scalar_tmp, drv_basis_scalar);
        ots::modifyBasisByOrientation(drv_basis_scalar,
                                      drv_basis_scalar_tmp,
                                      drv_orts,
                                      intrepid_basis);
      }
      if(build_weighted) {
        auto drv_basis_scalar = Kokkos::subview(weighted_basis_scalar.get_view(), range_cell, Kokkos::ALL(), Kokkos::ALL());
        auto drv_basis_scalar_tmp = Kokkos::createDynRankView(weighted_basis_scalar.get_view(),
                                                              "drv_basis_scalar_tmp",
                                                              drv_basis_scalar.extent(0),  // C
                                                              drv_basis_scalar.extent(1),  // F
                                                              drv_basis_scalar.extent(2)); // P
        Kokkos::deep_copy(drv_basis_scalar_tmp, drv_basis_scalar);
        ots::modifyBasisByOrientation(drv_basis_scalar,
                                      drv_basis_scalar_tmp,
                                      drv_orts,
                                      intrepid_basis);
      }

    }

    if (compute_derivatives) {
      {
        auto drv_grad_basis = Kokkos::subview(grad_basis.get_view(), range_cell, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
        auto drv_grad_basis_tmp = Kokkos::createDynRankView(grad_basis.get_view(),
                                                            "drv_grad_basis_tmp",
                                                            drv_grad_basis.extent(0),  // C
                                                            drv_grad_basis.extent(1),  // F
                                                            drv_grad_basis.extent(2),  // P
                                                            drv_grad_basis.extent(3)); // D
        Kokkos::deep_copy(drv_grad_basis_tmp, drv_grad_basis);
        ots::modifyBasisByOrientation(drv_grad_basis,
                                      drv_grad_basis_tmp,
                                      drv_orts,
                                      intrepid_basis);
      }
      if(build_weighted) {
        auto drv_grad_basis = Kokkos::subview(weighted_grad_basis.get_view(), range_cell, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
        auto drv_grad_basis_tmp = Kokkos::createDynRankView(weighted_grad_basis.get_view(),
                                                            "drv_grad_basis_tmp",
                                                            drv_grad_basis.extent(0),  // C
                                                            drv_grad_basis.extent(1),  // F
                                                            drv_grad_basis.extent(2),  // P
                                                            drv_grad_basis.extent(3)); // D
        Kokkos::deep_copy(drv_grad_basis_tmp, drv_grad_basis);
        ots::modifyBasisByOrientation(drv_grad_basis,
                                      drv_grad_basis_tmp,
                                      drv_orts,
                                      intrepid_basis);
      }
    }
  }

  ///
  /// hcurl 2d elements
  ///
  else if (elmtspace==PureBasis::HCURL && num_dim==2) {
    {
      {
        auto drv_basis_vector = Kokkos::subview(basis_vector.get_view(), range_cell, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
        auto drv_basis_vector_tmp = Kokkos::createDynRankView(basis_vector.get_view(),
                                                              "drv_basis_vector_tmp",
                                                              drv_basis_vector.extent(0),  // C
                                                              drv_basis_vector.extent(1),  // F
                                                              drv_basis_vector.extent(2),  // P
                                                              drv_basis_vector.extent(3)); // D
        Kokkos::deep_copy(drv_basis_vector_tmp, drv_basis_vector);
        ots::modifyBasisByOrientation(drv_basis_vector,
                                      drv_basis_vector_tmp,
                                      drv_orts,
                                      intrepid_basis);
      }
      if(build_weighted) {
        auto drv_basis_vector = Kokkos::subview(weighted_basis_vector.get_view(), range_cell, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
        auto drv_basis_vector_tmp = Kokkos::createDynRankView(weighted_basis_vector.get_view(),
                                                              "drv_basis_vector_tmp",
                                                              drv_basis_vector.extent(0),  // C
                                                              drv_basis_vector.extent(1),  // F
                                                              drv_basis_vector.extent(2),  // P
                                                              drv_basis_vector.extent(3)); // D
        Kokkos::deep_copy(drv_basis_vector_tmp, drv_basis_vector);
        ots::modifyBasisByOrientation(drv_basis_vector,
                                      drv_basis_vector_tmp,
                                      drv_orts,
                                      intrepid_basis);
      }
    }

    if (compute_derivatives) {
      {
        auto drv_curl_basis_scalar = Kokkos::subview(curl_basis_scalar.get_view(), range_cell, Kokkos::ALL(), Kokkos::ALL());
        auto drv_curl_basis_scalar_tmp = Kokkos::createDynRankView(curl_basis_scalar.get_view(),
                                                                   "drv_curl_basis_scalar_tmp",
                                                                   drv_curl_basis_scalar.extent(0),  // C
                                                                   drv_curl_basis_scalar.extent(1),  // F
                                                                   drv_curl_basis_scalar.extent(2));  // P
        Kokkos::deep_copy(drv_curl_basis_scalar_tmp, drv_curl_basis_scalar);
        ots::modifyBasisByOrientation(drv_curl_basis_scalar,
                                      drv_curl_basis_scalar_tmp,
                                      drv_orts,
                                      intrepid_basis);
      }

      if(build_weighted) {
        auto drv_curl_basis_scalar = Kokkos::subview(weighted_curl_basis_scalar.get_view(), range_cell, Kokkos::ALL(), Kokkos::ALL());
        auto drv_curl_basis_scalar_tmp = Kokkos::createDynRankView(weighted_curl_basis_scalar.get_view(),
                                                                   "drv_curl_basis_scalar_tmp",
                                                                   drv_curl_basis_scalar.extent(0),  // C
                                                                   drv_curl_basis_scalar.extent(1),  // F
                                                                   drv_curl_basis_scalar.extent(2));  // P
        Kokkos::deep_copy(drv_curl_basis_scalar_tmp, drv_curl_basis_scalar);
        ots::modifyBasisByOrientation(drv_curl_basis_scalar,
                                      drv_curl_basis_scalar_tmp,
                                      drv_orts,
                                      intrepid_basis);
      }
    }
  }

  ///
  /// hcurl 3d elements
  ///
  else if (elmtspace==PureBasis::HCURL && num_dim==3) {
    {
      {
        auto drv_basis_vector = Kokkos::subview(basis_vector.get_view(), range_cell, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
        auto drv_basis_vector_tmp = Kokkos::createDynRankView(basis_vector.get_view(),
                                                              "drv_basis_vector_tmp",
                                                              drv_basis_vector.extent(0),  // C
                                                              drv_basis_vector.extent(1),  // F
                                                              drv_basis_vector.extent(2),  // P
                                                              drv_basis_vector.extent(3)); // D
        Kokkos::deep_copy(drv_basis_vector_tmp, drv_basis_vector);
        ots::modifyBasisByOrientation(drv_basis_vector,
                                      drv_basis_vector_tmp,
                                      drv_orts,
                                      intrepid_basis);
      }
      if(build_weighted) {
        auto drv_basis_vector = Kokkos::subview(weighted_basis_vector.get_view(), range_cell, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
        auto drv_basis_vector_tmp = Kokkos::createDynRankView(weighted_basis_vector.get_view(),
                                                              "drv_basis_vector_tmp",
                                                              drv_basis_vector.extent(0),  // C
                                                              drv_basis_vector.extent(1),  // F
                                                              drv_basis_vector.extent(2),  // P
                                                              drv_basis_vector.extent(3)); // D
        Kokkos::deep_copy(drv_basis_vector_tmp, drv_basis_vector);
        ots::modifyBasisByOrientation(drv_basis_vector,
                                      drv_basis_vector_tmp,
                                      drv_orts,
                                      intrepid_basis);
      }
    }

    if (compute_derivatives) {
      {
        auto drv_curl_basis_vector = Kokkos::subview(curl_basis_vector.get_view(), range_cell, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
        auto drv_curl_basis_vector_tmp = Kokkos::createDynRankView(curl_basis_vector.get_view(),
                                                                   "drv_curl_basis_vector_tmp",
                                                                   drv_curl_basis_vector.extent(0),  // C
                                                                   drv_curl_basis_vector.extent(1),  // F
                                                                   drv_curl_basis_vector.extent(2),  // P
                                                                   drv_curl_basis_vector.extent(3));  // D
        Kokkos::deep_copy(drv_curl_basis_vector_tmp, drv_curl_basis_vector);
        ots::modifyBasisByOrientation(drv_curl_basis_vector,
                                      drv_curl_basis_vector_tmp,
                                      drv_orts,
                                      intrepid_basis);
      }
      if(build_weighted) {
        auto drv_curl_basis_vector = Kokkos::subview(weighted_curl_basis_vector.get_view(), range_cell, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
        auto drv_curl_basis_vector_tmp = Kokkos::createDynRankView(weighted_curl_basis_vector.get_view(),
                                                                   "drv_curl_basis_vector_tmp",
                                                                   drv_curl_basis_vector.extent(0),  // C
                                                                   drv_curl_basis_vector.extent(1),  // F
                                                                   drv_curl_basis_vector.extent(2),  // P
                                                                   drv_curl_basis_vector.extent(3));  // D
        Kokkos::deep_copy(drv_curl_basis_vector_tmp, drv_curl_basis_vector);
        ots::modifyBasisByOrientation(drv_curl_basis_vector,
                                      drv_curl_basis_vector_tmp,
                                      drv_orts,
                                      intrepid_basis);
      }
    }
  }
  ///
  /// hdiv elements (2d and 3d)
  ///
  else if (elmtspace==PureBasis::HDIV) {
    {
      {
        auto drv_basis_vector = Kokkos::subview(basis_vector.get_view(), range_cell, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
        auto drv_basis_vector_tmp = Kokkos::createDynRankView(basis_vector.get_view(),
                                                              "drv_basis_vector_tmp",
                                                              drv_basis_vector.extent(0),  // C
                                                              drv_basis_vector.extent(1),  // F
                                                              drv_basis_vector.extent(2),  // P
                                                              drv_basis_vector.extent(3)); // D
        Kokkos::deep_copy(drv_basis_vector_tmp, drv_basis_vector);
        ots::modifyBasisByOrientation(drv_basis_vector,
                                      drv_basis_vector_tmp,
                                      drv_orts,
                                      intrepid_basis);
      }
      if(build_weighted) {
        auto drv_basis_vector = Kokkos::subview(weighted_basis_vector.get_view(), range_cell, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
        auto drv_basis_vector_tmp = Kokkos::createDynRankView(weighted_basis_vector.get_view(),
                                                              "drv_basis_vector_tmp",
                                                              drv_basis_vector.extent(0),  // C
                                                              drv_basis_vector.extent(1),  // F
                                                              drv_basis_vector.extent(2),  // P
                                                              drv_basis_vector.extent(3)); // D
        Kokkos::deep_copy(drv_basis_vector_tmp, drv_basis_vector);
        ots::modifyBasisByOrientation(drv_basis_vector,
                                      drv_basis_vector_tmp,
                                      drv_orts,
                                      intrepid_basis);
      }
    }
    if (compute_derivatives) {
      {
        auto drv_div_basis = Kokkos::subview(div_basis.get_view(), range_cell, Kokkos::ALL(), Kokkos::ALL());
        auto drv_div_basis_tmp = Kokkos::createDynRankView(div_basis.get_view(),
                                                           "drv_div_basis_tmp",
                                                           drv_div_basis.extent(0),  // C
                                                           drv_div_basis.extent(1),  // F
                                                           drv_div_basis.extent(2));  // P
        Kokkos::deep_copy(drv_div_basis_tmp, drv_div_basis);
        ots::modifyBasisByOrientation(drv_div_basis,
                                      drv_div_basis_tmp,
                                      drv_orts,
                                      intrepid_basis);
      }
      if(build_weighted) {
        auto drv_div_basis = Kokkos::subview(weighted_div_basis.get_view(), range_cell, Kokkos::ALL(), Kokkos::ALL());
        auto drv_div_basis_tmp = Kokkos::createDynRankView(weighted_div_basis.get_view(),
                                                           "drv_div_basis_tmp",
                                                           drv_div_basis.extent(0),  // C
                                                           drv_div_basis.extent(1),  // F
                                                           drv_div_basis.extent(2));  // P
        Kokkos::deep_copy(drv_div_basis_tmp, drv_div_basis);
        ots::modifyBasisByOrientation(drv_div_basis,
                                      drv_div_basis_tmp,
                                      drv_orts,
                                      intrepid_basis);
      }
    }
  }
}

// method for applying orientations
template <typename Scalar>
void BasisValues2<Scalar>::
applyOrientations(const PHX::MDField<const Scalar,Cell,BASIS> & orientations)
{

  TEUCHOS_TEST_FOR_EXCEPT_MSG(true,"panzer::BasisValues2::applyOrientations : this should not be called.");

  int num_cell  = orientations.extent(0);
  int num_basis = orientations.extent(1);
  int num_dim   = basis_layout->dimension();
  int num_ip    = basis_layout->numPoints();
  PureBasis::EElementSpace elmtspace = getElementSpace();

  if(elmtspace==PureBasis::HCURL && num_dim==2) {

    // setup the orientations for the trial space
    // Intrepid2::FunctionSpaceTools::applyFieldSigns<Scalar>(basis_vector,orientations);

    for (int c=0; c<num_cell; c++)
      for (int b=0; b<num_basis; b++)
        for (int p=0; p<num_ip; p++)
          for (int d=0; d<num_dim; d++)
           basis_vector(c, b, p, d) *= orientations(c, b);

    if(compute_derivatives) {
      // Intrepid2::FunctionSpaceTools::applyFieldSigns<Scalar>(curl_basis_scalar,orientations);
      for (int c=0; c<num_cell; c++)
        for (int b=0; b<num_basis; b++)
          for (int p=0; p<num_ip; p++)
            curl_basis_scalar(c, b, p) *= orientations(c, b);
    }

    // setup the orientations for the test space
    if(build_weighted) {
      Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>::applyFieldSigns(weighted_basis_vector.get_view(),orientations.get_view());

      if(compute_derivatives)
        Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>::applyFieldSigns(weighted_curl_basis_scalar.get_view(),orientations.get_view());
    }
  }
  else if(elmtspace==PureBasis::HCURL && num_dim==3) {

    // setup the orientations for the trial space
    // Intrepid2::FunctionSpaceTools::applyFieldSigns<Scalar>(basis_vector,orientations);

    for (int c=0; c<num_cell; c++)
      for (int b=0; b<num_basis; b++)
        for (int p=0; p<num_ip; p++)
          for (int d=0; d<num_dim; d++)
           basis_vector(c, b, p, d) *= orientations(c, b);

    if(compute_derivatives) {
      // Intrepid2::FunctionSpaceTools::applyFieldSigns<Scalar>(curl_basis_vector,orientations);
      for (int c=0; c<num_cell; c++)
        for (int b=0; b<num_basis; b++)
          for (int p=0; p<num_ip; p++)
            for (int d=0; d<num_dim; d++)
              curl_basis_vector(c, b, p,d) *= orientations(c, b);
    }

    // setup the orientations for the test space
    if(build_weighted) {
      Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>::applyFieldSigns(weighted_basis_vector.get_view(),orientations.get_view());

      if(compute_derivatives)
        Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>::applyFieldSigns(weighted_curl_basis_vector.get_view(),orientations.get_view());
    }
  }
  else if(elmtspace==PureBasis::HDIV) {
    // setup the orientations for the trial space
    // Intrepid2::FunctionSpaceTools::applyFieldSigns<Scalar>(basis_vector,orientations);

    for (int c=0; c<num_cell; c++)
      for (int b=0; b<num_basis; b++)
        for (int p=0; p<num_ip; p++)
          for (int d=0; d<num_dim; d++)
           basis_vector(c, b, p, d) *= orientations(c, b);

    if(compute_derivatives) {
      // Intrepid2::FunctionSpaceTools::applyFieldSigns<Scalar>(div_basis,orientations);

      for (int c=0; c<num_cell; c++)
        for (int b=0; b<num_basis; b++)
          for (int p=0; p<num_ip; p++)
            div_basis(c, b, p) *= orientations(c, b);
    }

    // setup the orientations for the test space
    if(build_weighted) {
      Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>::applyFieldSigns(weighted_basis_vector.get_view(),orientations.get_view());

      if(compute_derivatives)
        Intrepid2::FunctionSpaceTools<PHX::Device::execution_space>::applyFieldSigns(weighted_div_basis.get_view(),orientations.get_view());
    }
  }
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

} // namespace panzer
