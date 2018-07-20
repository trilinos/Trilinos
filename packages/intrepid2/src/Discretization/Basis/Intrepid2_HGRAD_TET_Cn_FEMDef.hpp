// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov), or
//                    Mauro Perego  (mperego@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid2_HGRAD_TET_Cn_FEMDef.hpp
    \brief  Definition file for FEM basis functions of degree n for H(grad) functions on TET cells.
    \author Created by R. Kirby and P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
 */

#ifndef __INTREPID2_HGRAD_TET_CN_FEM_DEF_HPP__
#define __INTREPID2_HGRAD_TET_CN_FEM_DEF_HPP__

#include "Intrepid2_HGRAD_TET_Cn_FEM_ORTH.hpp"
#include "Intrepid2_HGRAD_TET_Cn_FEM.hpp"

namespace Intrepid2 {

// -------------------------------------------------------------------------------------
namespace Impl {

template<EOperator opType>
template<typename outputViewType,
typename inputViewType,
typename workViewType,
typename vinvViewType>
KOKKOS_INLINE_FUNCTION
void
Basis_HGRAD_TET_Cn_FEM::Serial<opType>::
getValues(       outputViewType output,
    const inputViewType  input,
    workViewType   work,
    const vinvViewType   vinv ) {

  constexpr ordinal_type spaceDim = 3;
  const ordinal_type
  card = vinv.extent(0),
  npts = input.extent(0);

  // compute order
  ordinal_type order = 0;
  for (ordinal_type p=0;p<=Parameters::MaxOrder;++p) {
    if (card == Intrepid2::getPnCardinality<spaceDim>(p)) {
      order = p;
      break;
    }
  }

  typedef typename Kokkos::DynRankView<typename workViewType::value_type, typename workViewType::memory_space> viewType;
  auto vcprop = Kokkos::common_view_alloc_prop(work);
  auto ptr = work.data();

  switch (opType) {
  case OPERATOR_VALUE: {
    const viewType phis(Kokkos::view_wrap(ptr, vcprop), card, npts);
    viewType dummyView;

    Impl::Basis_HGRAD_TET_Cn_FEM_ORTH::
    Serial<opType>::getValues(phis, input, dummyView, order);

    for (ordinal_type i=0;i<card;++i)
      for (ordinal_type j=0;j<npts;++j) {
        output.access(i,j) = 0.0;
        for (ordinal_type k=0;k<card;++k)
          output.access(i,j) += vinv(k,i)*phis.access(k,j);
      }
    break;
  }
  case OPERATOR_GRAD:
  case OPERATOR_D1: {
    const viewType phis(Kokkos::view_wrap(ptr, vcprop), card, npts, spaceDim);
    ptr += card*npts*spaceDim*get_dimension_scalar(work);
    const viewType workView(Kokkos::view_wrap(ptr, vcprop), card, npts, spaceDim+1);
    Impl::Basis_HGRAD_TET_Cn_FEM_ORTH::
    Serial<opType>::getValues(phis, input, workView, order);

    for (ordinal_type i=0;i<card;++i)
      for (ordinal_type j=0;j<npts;++j)
        for (ordinal_type k=0;k<spaceDim;++k) {
          output.access(i,j,k) = 0.0;
          for (ordinal_type l=0;l<card;++l)
            output.access(i,j,k) += vinv(l,i)*phis.access(l,j,k);
        }
    break;
  }
  case OPERATOR_D2:
  case OPERATOR_D3:
  case OPERATOR_D4:
  case OPERATOR_D5:
  case OPERATOR_D6:
  case OPERATOR_D7:
  case OPERATOR_D8:
  case OPERATOR_D9:
  case OPERATOR_D10: {
    const ordinal_type dkcard = getDkCardinality<opType,spaceDim>(); //(orDn + 1);
    const viewType phis(Kokkos::view_wrap(ptr, vcprop), card, npts, dkcard);
    viewType dummyView;

    Impl::Basis_HGRAD_TET_Cn_FEM_ORTH::
    Serial<opType>::getValues(phis, input, dummyView, order);

    for (ordinal_type i=0;i<card;++i)
      for (ordinal_type j=0;j<npts;++j)
        for (ordinal_type k=0;k<dkcard;++k) {
          output.access(i,j,k) = 0.0;
          for (ordinal_type l=0;l<card;++l)
            output.access(i,j,k) += vinv(l,i)*phis.access(l,j,k);
        }
    break;
  }
  default: {
    INTREPID2_TEST_FOR_ABORT( true,
        ">>> ERROR (Basis_HGRAD_TET_Cn_FEM): Operator type not implemented");
  }
  }
}

template<typename SpT, ordinal_type numPtsPerEval,
typename outputValueValueType, class ...outputValueProperties,
typename inputPointValueType,  class ...inputPointProperties,
typename vinvValueType,        class ...vinvProperties>
void
Basis_HGRAD_TET_Cn_FEM::
getValues(       Kokkos::DynRankView<outputValueValueType,outputValueProperties...> outputValues,
    const Kokkos::DynRankView<inputPointValueType, inputPointProperties...>  inputPoints,
    const Kokkos::DynRankView<vinvValueType,       vinvProperties...>        vinv,
    const EOperator operatorType) {
  typedef          Kokkos::DynRankView<outputValueValueType,outputValueProperties...>         outputValueViewType;
  typedef          Kokkos::DynRankView<inputPointValueType, inputPointProperties...>          inputPointViewType;
  typedef          Kokkos::DynRankView<vinvValueType,       vinvProperties...>                vinvViewType;
  typedef typename ExecSpace<typename inputPointViewType::execution_space,SpT>::ExecSpaceType ExecSpaceType;

  // loopSize corresponds to cardinality
  const auto loopSizeTmp1 = (inputPoints.extent(0)/numPtsPerEval);
  const auto loopSizeTmp2 = (inputPoints.extent(0)%numPtsPerEval != 0);
  const auto loopSize = loopSizeTmp1 + loopSizeTmp2;
  Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);

  typedef typename inputPointViewType::value_type inputPointType;

  const ordinal_type cardinality = outputValues.extent(0);
  const ordinal_type spaceDim = 3;

  auto vcprop = Kokkos::common_view_alloc_prop(inputPoints);
  typedef typename Kokkos::DynRankView< inputPointType, typename inputPointViewType::memory_space> workViewType;

  switch (operatorType) {
  case OPERATOR_VALUE: {
    workViewType  work(Kokkos::view_alloc("Basis_HGRAD_TET_Cn_FEM::getValues::work", vcprop), cardinality, inputPoints.extent(0));
    typedef Functor<outputValueViewType,inputPointViewType,vinvViewType, workViewType,
        OPERATOR_VALUE,numPtsPerEval> FunctorType;
    Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints, vinv, work) );
    break;
  }
  case OPERATOR_GRAD:
  case OPERATOR_D1: {
    workViewType  work(Kokkos::view_alloc("Basis_HGRAD_TET_Cn_FEM::getValues::work", vcprop), cardinality*(2*spaceDim+1), inputPoints.extent(0));
    typedef Functor<outputValueViewType,inputPointViewType,vinvViewType, workViewType,
        OPERATOR_D1,numPtsPerEval> FunctorType;
    Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints, vinv, work) );
    break;
  }
  case OPERATOR_D2: {
    typedef Functor<outputValueViewType,inputPointViewType,vinvViewType, workViewType,
        OPERATOR_D2,numPtsPerEval> FunctorType;
    workViewType  work(Kokkos::view_alloc("Basis_HGRAD_TET_Cn_FEM::getValues::work", vcprop), cardinality*outputValues.extent(2), inputPoints.extent(0));
    Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints, vinv, work) );
    break;
  }
  /*  case OPERATOR_D3: {
        typedef Functor<outputValueViewType,inputPointViewType,vinvViewType, workViewType
            OPERATOR_D3,numPtsPerEval> FunctorType;
        workViewType  work(Kokkos::view_alloc("Basis_HGRAD_TET_Cn_FEM::getValues::work", vcprop), cardinality, inputPoints.extent(0), outputValues.extent(2));
        Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints, vinv, work) );
        break;
      }*/
  default: {
    INTREPID2_TEST_FOR_EXCEPTION( true , std::invalid_argument,
        ">>> ERROR (Basis_HGRAD_TET_Cn_FEM): Operator type not implemented" );
  }
  }
}
}

// -------------------------------------------------------------------------------------
template<typename SpT, typename OT, typename PT>
Basis_HGRAD_TET_Cn_FEM<SpT,OT,PT>::
Basis_HGRAD_TET_Cn_FEM( const ordinal_type order,
    const EPointType   pointType ) {
  constexpr ordinal_type spaceDim = 3;

  this->basisCardinality_  = Intrepid2::getPnCardinality<spaceDim>(order); // bigN
  this->basisDegree_       = order; // small n
  this->basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Tetrahedron<4> >() );
  this->basisType_         = BASIS_FEM_FIAT;
  this->basisCoordinates_  = COORDINATES_CARTESIAN;

  const ordinal_type card = this->basisCardinality_;

  // points are computed in the host and will be copied
  Kokkos::DynRankView<scalarType,typename SpT::array_layout,Kokkos::HostSpace>
  dofCoords("Hgrad::Tet::Cn::dofCoords", card, spaceDim);

  // Basis-dependent initializations
  constexpr ordinal_type tagSize  = 4;        // size of DoF tag, i.e., number of fields in the tag
  constexpr ordinal_type maxCard = Intrepid2::getPnCardinality<spaceDim, Parameters::MaxOrder>();
  ordinal_type tags[maxCard][tagSize];

  // construct lattice

  const ordinal_type numEdges = this->basisCellTopology_.getEdgeCount();
  const ordinal_type numFaces = this->basisCellTopology_.getFaceCount();

  shards::CellTopology edgeTop(shards::getCellTopologyData<shards::Line<2> >() );
  shards::CellTopology faceTop(shards::getCellTopologyData<shards::Triangle<3> >() );

  const int numVertexes = PointTools::getLatticeSize( this->basisCellTopology_ ,
      1 ,
      0 );

  const int numPtsPerEdge = PointTools::getLatticeSize( edgeTop ,
      order ,
      1 );

  const int numPtsPerFace = PointTools::getLatticeSize( faceTop ,
      order ,
      1 );

  const int numPtsPerCell = PointTools::getLatticeSize( this->basisCellTopology_ ,
      order ,
      1 );

  Kokkos::DynRankView<scalarType,typename SpT::array_layout,Kokkos::HostSpace> vertexes("Hcurl::Tet::In::vertexes", numVertexes , spaceDim );
  Kokkos::DynRankView<scalarType,typename SpT::array_layout,Kokkos::HostSpace> linePts("Hcurl::Tet::In::linePts", numPtsPerEdge , 1 );
  Kokkos::DynRankView<scalarType,typename SpT::array_layout,Kokkos::HostSpace> triPts("Hcurl::Tet::In::triPts", numPtsPerFace , 2 );

  // construct lattice
  const ordinal_type offset = 1;


  PointTools::getLattice( vertexes,
      this->basisCellTopology_ ,
      1, 0,
      pointType );

  PointTools::getLattice( linePts,
      edgeTop,
      order, offset,
      pointType );

  PointTools::getLattice( triPts,
      faceTop,
      order, offset,
      pointType );

  // holds the image of the line points
  Kokkos::DynRankView<scalarType,typename SpT::array_layout,Kokkos::HostSpace> edgePts("Hcurl::Tet::In::edgePts", numPtsPerEdge , spaceDim );
  Kokkos::DynRankView<scalarType,typename SpT::array_layout,Kokkos::HostSpace> facePts("Hcurl::Tet::In::facePts", numPtsPerFace , spaceDim );

  for (ordinal_type i=0;i<numVertexes;i++) {
    auto i_card=i;
    for(ordinal_type k=0; k<spaceDim; ++k)
      dofCoords(i_card,k) = vertexes(i,k);
    tags[i_card][0] = 0; // vertex dof
    tags[i_card][1] = i; // vertex id
    tags[i_card][2] = 0; // local dof id
    tags[i_card][3] = 1; // total vert dof
  }


  // these are tangents scaled by the appropriate edge lengths.
  for (ordinal_type i=0;i<numEdges;i++) {  // loop over edges
    CellTools<Kokkos::HostSpace::execution_space>::mapToReferenceSubcell( edgePts ,
        linePts ,
        1 ,
        i ,
        this->basisCellTopology_ );


    // loop over points (rows of V2)
    for (ordinal_type j=0;j<numPtsPerEdge;j++) {

      const ordinal_type i_card = numVertexes + numPtsPerEdge*i+j;

      //save dof coordinates and coefficients
      for(ordinal_type k=0; k<spaceDim; ++k)
        dofCoords(i_card,k) = edgePts(j,k);

      tags[i_card][0] = 1; // edge dof
      tags[i_card][1] = i; // edge id
      tags[i_card][2] = j; // local dof id
      tags[i_card][3] = numPtsPerEdge; // total edge dof

    }
  }

  if(numPtsPerFace >0) {//handle faces if needed  (order >1)

    for (ordinal_type i=0;i<numFaces;i++) {  // loop over faces

      CellTools<Kokkos::HostSpace::execution_space>::mapToReferenceSubcell( facePts ,
          triPts ,
          2 ,
          i ,
          this->basisCellTopology_ );
      for (ordinal_type j=0;j<numPtsPerFace;j++) {

        const ordinal_type i_card = numVertexes+numEdges*numPtsPerEdge+numPtsPerFace*i+j;

        //save dof coordinates
        for(ordinal_type k=0; k<spaceDim; ++k)
          dofCoords(i_card,k) = facePts(j,k);

        tags[i_card][0] = 2; // face dof
        tags[i_card][1] = i; // face id
        tags[i_card][2] = j; // local face id
        tags[i_card][3] = numPtsPerFace; // total face dof
      }
    }
  }


  // internal dof, if needed
  if (numPtsPerCell > 0) {
    Kokkos::DynRankView<scalarType,typename SpT::array_layout,Kokkos::HostSpace>
    cellPoints( "Hcurl::Tet::In::cellPoints", numPtsPerCell , spaceDim );
    PointTools::getLattice( cellPoints ,
        this->basisCellTopology_ ,
        order,
        1 ,
        pointType );

    // copy values into right positions of V2
    for (ordinal_type j=0;j<numPtsPerCell;j++) {

      const ordinal_type i_card = numVertexes+numEdges*numPtsPerEdge+numFaces*numPtsPerFace+j;

      //save dof coordinates
      for(ordinal_type dim=0; dim<spaceDim; ++dim)
        dofCoords(i_card,dim) = cellPoints(j,dim);

      tags[i_card][0] = spaceDim; // elem dof
      tags[i_card][1] = 0; // elem id
      tags[i_card][2] = j; // local dof id
      tags[i_card][3] = numPtsPerCell; // total vert dof
    }
  }

  this->dofCoords_ = Kokkos::create_mirror_view(typename SpT::memory_space(), dofCoords);
  Kokkos::deep_copy(this->dofCoords_, dofCoords);

  // form Vandermonde matrix.  Actually, this is the transpose of the VDM,
  // so we transpose on copy below.
  const ordinal_type lwork = card*card;
  Kokkos::DynRankView<scalarType,Kokkos::LayoutLeft,Kokkos::HostSpace>
  vmat("Hgrad::Tet::Cn::vmat", card, card),
  work("Hgrad::Tet::Cn::work", lwork),
  ipiv("Hgrad::Tet::Cn::ipiv", card);

  Impl::Basis_HGRAD_TET_Cn_FEM_ORTH::getValues<Kokkos::HostSpace::execution_space,Parameters::MaxNumPtsPerBasisEval>(vmat, dofCoords, order, OPERATOR_VALUE);

  ordinal_type info = 0;
  Teuchos::LAPACK<ordinal_type,scalarType> lapack;

  lapack.GETRF(card, card,
      vmat.data(), vmat.stride_1(),
      (ordinal_type*)ipiv.data(),
      &info);

  INTREPID2_TEST_FOR_EXCEPTION( info != 0,
      std::runtime_error ,
      ">>> ERROR: (Intrepid2::Basis_HGRAD_TET_Cn_FEM) lapack.GETRF returns nonzero info." );

  lapack.GETRI(card,
      vmat.data(), vmat.stride_1(),
      (ordinal_type*)ipiv.data(),
      work.data(), lwork,
      &info);

  INTREPID2_TEST_FOR_EXCEPTION( info != 0,
      std::runtime_error ,
      ">>> ERROR: (Intrepid2::Basis_HGRAD_TET_Cn_FEM) lapack.GETRI returns nonzero info." );

  // create host mirror
  Kokkos::DynRankView<scalarType,typename SpT::array_layout,Kokkos::HostSpace>
  vinv("Hgrad::Line::Cn::vinv", card, card);

  for (ordinal_type i=0;i<card;++i)
    for (ordinal_type j=0;j<card;++j)
      vinv(i,j) = vmat(j,i);

  this->vinv_ = Kokkos::create_mirror_view(typename SpT::memory_space(), vinv);
  Kokkos::deep_copy(this->vinv_ , vinv);

  // initialize tags
  {
    // Basis-dependent initializations
    const ordinal_type posScDim = 0;        // position in the tag, counting from 0, of the subcell dim
    const ordinal_type posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
    const ordinal_type posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell

    ordinal_type_array_1d_host tagView(&tags[0][0], card*tagSize);

    // Basis-independent function sets tag and enum data in tagToOrdinal_ and ordinalToTag_ arrays:
    // tags are constructed on host
    this->setOrdinalTagData(this->tagToOrdinal_,
        this->ordinalToTag_,
        tagView,
        this->basisCardinality_,
        tagSize,
        posScDim,
        posScOrd,
        posDfOrd);
  }
}
} // namespace Intrepid2
#endif
