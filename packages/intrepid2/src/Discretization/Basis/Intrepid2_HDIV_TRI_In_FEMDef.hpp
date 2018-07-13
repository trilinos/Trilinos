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

/** \file   Intrepid2_HDIV_TRI_In_FEMDef.hpp
    \brief  Definition file for FEM basis functions of degree n for H(div) functions on TRI cells.
    \author Created by R. Kirby and P. Bochev and D. Ridzal.
    Kokkorized by Kyungjoo Kim
 */

#ifndef __INTREPID2_HDIV_TRI_IN_FEM_DEF_HPP__
#define __INTREPID2_HDIV_TRI_IN_FEM_DEF_HPP__

#include "Intrepid2_HGRAD_TRI_Cn_FEM_ORTH.hpp"
#include "Intrepid2_CubatureDirectTriDefault.hpp"

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
Basis_HDIV_TRI_In_FEM::Serial<opType>::
getValues( /* */ outputViewType output,
    const inputViewType  input,
    /* */ workViewType   work,
    const vinvViewType   coeffs ) {

  constexpr ordinal_type spaceDim = 2;
  const ordinal_type
  cardPn = coeffs.extent(0)/spaceDim,
  card = coeffs.extent(1),
  npts = input.extent(0);

  // compute order
  ordinal_type order = 0;
  for (ordinal_type p=0;p<=Parameters::MaxOrder;++p) {
    if (card == CardinalityHDivTri(p)) {
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
    workViewType dummyView;

    Impl::Basis_HGRAD_TRI_Cn_FEM_ORTH::
    Serial<opType>::getValues(phis, input, dummyView, order);

    for (ordinal_type i=0;i<card;++i)
      for (ordinal_type j=0;j<npts;++j)
        for (ordinal_type d=0;d<spaceDim;++d) {
          output.access(i,j,d) = 0.0;
          for (ordinal_type k=0;k<cardPn;++k)
            output.access(i,j,d) += coeffs(k+d*cardPn,i) * phis.access(k,j);
        }
    break;
  }
  case OPERATOR_DIV: {
    const viewType phis(Kokkos::view_wrap(ptr, vcprop), card, npts, spaceDim);
    ptr += card*npts*spaceDim*get_dimension_scalar(work);
    const viewType workView(Kokkos::view_wrap(ptr, vcprop), card, npts, spaceDim+1);

    Impl::Basis_HGRAD_TRI_Cn_FEM_ORTH::
    Serial<OPERATOR_GRAD>::getValues(phis, input, workView, order);

    for (ordinal_type i=0;i<card;++i)
      for (ordinal_type j=0;j<npts;++j) {
        output.access(i,j) = 0.0;
        for (ordinal_type k=0; k<cardPn; ++k)
          for (ordinal_type d=0; d<spaceDim; ++d)
            output.access(i,j) += coeffs(k+d*cardPn,i)*phis.access(k,j,d);
      }
    break;
  }
  default: {
    INTREPID2_TEST_FOR_ABORT( true,
        ">>> ERROR (Basis_HDIV_TRI_In_FEM): Operator type not implemented");
  }
  }
}

template<typename SpT, ordinal_type numPtsPerEval,
typename outputValueValueType, class ...outputValueProperties,
typename inputPointValueType,  class ...inputPointProperties,
typename vinvValueType,        class ...vinvProperties>
void
Basis_HDIV_TRI_In_FEM::
getValues( /* */ Kokkos::DynRankView<outputValueValueType,outputValueProperties...> outputValues,
    const Kokkos::DynRankView<inputPointValueType, inputPointProperties...>  inputPoints,
    const Kokkos::DynRankView<vinvValueType,       vinvProperties...>        coeffs,
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
  const ordinal_type spaceDim = 2;

  auto vcprop = Kokkos::common_view_alloc_prop(inputPoints);
  typedef typename Kokkos::DynRankView< inputPointType, typename inputPointViewType::memory_space> workViewType;

  switch (operatorType) {
  case OPERATOR_VALUE: {
    workViewType  work(Kokkos::view_alloc("Basis_HDIV_TRI_In_FEM::getValues::work", vcprop), cardinality, inputPoints.extent(0));
    typedef Functor<outputValueViewType,inputPointViewType,vinvViewType, workViewType,
        OPERATOR_VALUE,numPtsPerEval> FunctorType;
    Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints, coeffs, work) );
    break;
  }
  case OPERATOR_DIV: {
    workViewType  work(Kokkos::view_alloc("Basis_HDIV_TRI_In_FEM::getValues::work", vcprop), cardinality*(2*spaceDim+1), inputPoints.extent(0));
    typedef Functor<outputValueViewType,inputPointViewType,vinvViewType, workViewType,
        OPERATOR_DIV,numPtsPerEval> FunctorType;
    Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints, coeffs, work) );
    break;
  }
  default: {
    INTREPID2_TEST_FOR_EXCEPTION( true , std::invalid_argument,
        ">>> ERROR (Basis_HDIV_TRI_In_FEM): Operator type not implemented" );
  }
  }
}
}

// -------------------------------------------------------------------------------------
template<typename SpT, typename OT, typename PT>
Basis_HDIV_TRI_In_FEM<SpT,OT,PT>::
Basis_HDIV_TRI_In_FEM( const ordinal_type order,
    const EPointType   pointType ) {

  constexpr ordinal_type spaceDim = 2;
  this->basisCardinality_  = CardinalityHDivTri(order);
  this->basisDegree_       = order; // small n
  this->basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Triangle<3> >() );
  this->basisType_         = BASIS_FEM_FIAT;
  this->basisCoordinates_  = COORDINATES_CARTESIAN;

  const ordinal_type card = this->basisCardinality_;

  const ordinal_type  cardPn = Intrepid2::getPnCardinality<spaceDim>(order); // dim of (P_{n}) -- smaller space
  const ordinal_type  cardPnm1 = Intrepid2::getPnCardinality<spaceDim>(order-1);  // dim of (P_{n-1}) -- smaller space
  const ordinal_type  cardPnm2 = Intrepid2::getPnCardinality<spaceDim>(order-2); // dim of (P_{n-2}) -- smaller space
  const ordinal_type  cardVecPn = spaceDim*cardPn;  // dim of (P_{n})^2 -- larger space
  const ordinal_type  cardVecPnm1 = spaceDim*cardPnm1;   // dim of (P_{n-1})^2 -- smaller space


  // Basis-dependent initializations
  constexpr ordinal_type tagSize  = 4;        // size of DoF tag, i.e., number of fields in the tag
  constexpr ordinal_type maxCard = CardinalityHDivTri(Parameters::MaxOrder);
  ordinal_type tags[maxCard][tagSize];

  // points are computed in the host and will be copied
  Kokkos::DynRankView<scalarType,typename SpT::array_layout,Kokkos::HostSpace>
  dofCoords("Hdiv::Tri::In::dofCoords", card, spaceDim);

  Kokkos::DynRankView<scalarType,typename SpT::array_layout,Kokkos::HostSpace>
  dofCoeffs("Hdiv::Tri::In::dofCoeffs", card, spaceDim);

  Kokkos::DynRankView<scalarType,typename SpT::array_layout,Kokkos::HostSpace>
  coeffs("Hdiv::Tri::In::coeffs", cardVecPn, card);

  // first, need to project the basis for RT space onto the
  // orthogonal basis of degree n
  // get coefficients of PkHx

  const ordinal_type lwork = card*card;
  Kokkos::DynRankView<scalarType,typename SpT::array_layout,Kokkos::HostSpace>
  V1("Hdiv::Tri::In::V1", cardVecPn, card);

  // basis for the space is
  // { (phi_i,0) }_{i=0}^{cardPnm1-1} ,
  // { (0,phi_i) }_{i=0}^{cardPnm1-1} ,
  // { (x,y) . phi_i}_{i=cardPnm2}^{cardPnm1-1}
  // columns of V1 are expansion of this basis in terms of the basis
  // for P_{n}^2

  // these two loops get the first two sets of basis functions
  for (ordinal_type i=0;i<cardPnm1;i++) {
    V1(i,i) = 1.0;
    V1(cardPn+i,cardPnm1+i) = 1.0;
  }

  // now I need to integrate { (x,y) phi } against the big basis
  // first, get a cubature rule.
  CubatureDirectTriDefault<Kokkos::HostSpace::execution_space,scalarType,scalarType> myCub( 2 * order );
  Kokkos::DynRankView<scalarType,typename SpT::array_layout,Kokkos::HostSpace> cubPoints("Hdiv::Tri::In::cubPoints", myCub.getNumPoints() , spaceDim );
  Kokkos::DynRankView<scalarType,typename SpT::array_layout,Kokkos::HostSpace> cubWeights("Hdiv::Tri::In::cubWeights", myCub.getNumPoints() );
  myCub.getCubature( cubPoints , cubWeights );

  // tabulate the scalar orthonormal basis at cubature points
  Kokkos::DynRankView<scalarType,typename SpT::array_layout,Kokkos::HostSpace> phisAtCubPoints("Hdiv::Tri::In::phisAtCubPoints", cardPn , myCub.getNumPoints() );
  Impl::Basis_HGRAD_TRI_Cn_FEM_ORTH::getValues<Kokkos::HostSpace::execution_space,Parameters::MaxNumPtsPerBasisEval>(phisAtCubPoints, cubPoints, order, OPERATOR_VALUE);

  // now do the integration
  for (ordinal_type i=0;i<order;i++) {
    for (ordinal_type j=0;j<cardPn;j++) { // int (x,y) phi_i \cdot (phi_j,phi_{j+cardPn})
      V1(j,cardVecPnm1+i) = 0.0;
      for (ordinal_type d=0; d< spaceDim; ++d)
        for (ordinal_type k=0;k<myCub.getNumPoints();k++) {
          V1(j+d*cardPn,cardVecPnm1+i) +=
              cubWeights(k) * cubPoints(k,d)
              * phisAtCubPoints(cardPnm2+i,k)
              * phisAtCubPoints(j,k);
        }
    }
  }

  // next, apply the RT nodes (rows) to the basis for (P_n)^2 (columns)
  Kokkos::DynRankView<scalarType,typename SpT::array_layout,Kokkos::HostSpace>
  V2("Hdiv::Tri::In::V2", card ,cardVecPn);

  const ordinal_type numEdges = this->basisCellTopology_.getEdgeCount();

  shards::CellTopology edgeTop(shards::getCellTopologyData<shards::Line<2> >() );

  const int numPtsPerEdge = PointTools::getLatticeSize( edgeTop ,
      order+1 ,
      1 );

  // first numEdges * degree nodes are normals at each edge
  // get the points on the line
  Kokkos::DynRankView<scalarType,typename SpT::array_layout,Kokkos::HostSpace> linePts("Hdiv::Tri::In::linePts", numPtsPerEdge , 1 );

  // construct lattice
  const ordinal_type offset = 1;
  PointTools::getLattice( linePts,
      edgeTop,
      order+1, offset,
      pointType );

  // holds the image of the line points
  Kokkos::DynRankView<scalarType,typename SpT::array_layout,Kokkos::HostSpace> edgePts("Hdiv::Tri::In::edgePts", numPtsPerEdge , spaceDim );
  Kokkos::DynRankView<scalarType,typename SpT::array_layout,Kokkos::HostSpace> phisAtEdgePoints("Hdiv::Tri::In::phisAtEdgePoints", cardPn , numPtsPerEdge );
  Kokkos::DynRankView<scalarType,typename SpT::array_layout,Kokkos::HostSpace> edgeNormal("Hcurl::Tri::In::edgeNormal", spaceDim );

  // these are normal scaled by the appropriate edge lengths.
  for (ordinal_type edge=0;edge<numEdges;edge++) {  // loop over edges
    CellTools<Kokkos::HostSpace::execution_space>::getReferenceSideNormal( edgeNormal ,
        edge ,
        this->basisCellTopology_ );

    /* multiply by measure of reference edge so that magnitude of the edgeTan is equal to the edge measure */
    const scalarType refEdgeMeasure = 2.0;
    for (ordinal_type j=0;j<spaceDim;j++)
      edgeNormal(j) *= refEdgeMeasure;


    CellTools<Kokkos::HostSpace::execution_space>::mapToReferenceSubcell( edgePts ,
        linePts ,
        1 ,
        edge ,
        this->basisCellTopology_ );

    Impl::Basis_HGRAD_TRI_Cn_FEM_ORTH::getValues<Kokkos::HostSpace::execution_space,Parameters::MaxNumPtsPerBasisEval>(phisAtEdgePoints , edgePts, order, OPERATOR_VALUE);

    // loop over points (rows of V2)
    for (ordinal_type j=0;j<numPtsPerEdge;j++) {

      const ordinal_type i_card = numPtsPerEdge*edge+j;

      // loop over orthonormal basis functions (columns of V2)
      for (ordinal_type k=0;k<cardPn;k++) {
        // loop over space dimension
        for (ordinal_type l=0; l<spaceDim; l++)
          V2(i_card,k+l*cardPn) = edgeNormal(l) * phisAtEdgePoints(k,j);
      }


      //save dof coordinates and coefficients
      for(ordinal_type l=0; l<spaceDim; ++l) {
        dofCoords(i_card,l) = edgePts(j,l);
        dofCoeffs(i_card,l) = edgeNormal(l);
      }

      tags[i_card][0] = 1; // edge dof
      tags[i_card][1] = edge; // edge id
      tags[i_card][2] = j; // local dof id
      tags[i_card][3] = numPtsPerEdge; // total vert dof

    }


  }

  // remaining nodes are divided into two pieces:  point value of x
  // components and point values of y components.  These are
  // evaluated at the interior of a lattice of degree + 1, For then
  // the degree == 1 space corresponds classicaly to RT0 and so gets
  // no internal nodes, and degree == 2 corresponds to RT1 and needs
  // one internal node per vector component.
  const ordinal_type numPtsPerCell = PointTools::getLatticeSize( this->basisCellTopology_ ,
      order + 1 ,
      1 );

  if (numPtsPerCell > 0) {
    Kokkos::DynRankView<scalarType,typename SpT::array_layout,Kokkos::HostSpace>
    internalPoints( "Hdiv::Tri::In::internalPoints", numPtsPerCell , spaceDim );
    PointTools::getLattice( internalPoints ,
        this->basisCellTopology_ ,
        order + 1 ,
        1 ,
        pointType );

    Kokkos::DynRankView<scalarType,typename SpT::array_layout,Kokkos::HostSpace>
    phisAtInternalPoints("Hdiv::Tri::In::phisAtInternalPoints", cardPn , numPtsPerCell );
    Impl::Basis_HGRAD_TRI_Cn_FEM_ORTH::getValues<Kokkos::HostSpace::execution_space,Parameters::MaxNumPtsPerBasisEval>( phisAtInternalPoints , internalPoints , order, OPERATOR_VALUE );

    // copy values into right positions of V2
    for (ordinal_type j=0;j<numPtsPerCell;j++) {

      const ordinal_type i_card = numEdges*order+j;

      for (ordinal_type k=0;k<cardPn;k++) {
        for (ordinal_type l=0;l<spaceDim;l++) {
          V2(i_card+l*numPtsPerCell,l*cardPn+k) = phisAtInternalPoints(k,j);
        }
      }

      //save dof coordinates and coefficients
      for(ordinal_type d=0; d<spaceDim; ++d) {
        for(ordinal_type l=0; l<spaceDim; ++l) {
          dofCoords(i_card+d*numPtsPerCell,l) = internalPoints(j,l);
          dofCoeffs(i_card+d*numPtsPerCell,l) = (l==d);
        }

        tags[i_card+d*numPtsPerCell][0] = spaceDim; // elem dof
        tags[i_card+d*numPtsPerCell][1] = 0; // elem id
        tags[i_card+d*numPtsPerCell][2] = spaceDim*j+d; // local dof id
        tags[i_card+d*numPtsPerCell][3] = spaceDim*numPtsPerCell; // total vert dof
      }
    }
  }

  // form Vandermonde matrix.  Actually, this is the transpose of the VDM,
  // so we transpose on copy below.
  Kokkos::DynRankView<scalarType,Kokkos::LayoutLeft,Kokkos::HostSpace>
  vmat("Hdiv::Tri::In::vmat", card, card),
  work("Hdiv::Tri::In::work", lwork),
  ipiv("Hdiv::Tri::In::ipiv", card);

  //vmat' = V2*V1;
  for(ordinal_type i=0; i< card; ++i) {
    for(ordinal_type j=0; j< card; ++j) {
      scalarType s=0;
      for(ordinal_type k=0; k< cardVecPn; ++k)
        s += V2(i,k)*V1(k,j);
      vmat(i,j) = s;
    }
  }

  ordinal_type info = 0;
  Teuchos::LAPACK<ordinal_type,scalarType> lapack;

  lapack.GETRF(card, card,
      vmat.data(), vmat.stride_1(),
      (ordinal_type*)ipiv.data(),
      &info);

  INTREPID2_TEST_FOR_EXCEPTION( info != 0,
      std::runtime_error ,
      ">>> ERROR: (Intrepid2::Basis_HDIV_TRI_In_FEM) lapack.GETRF returns nonzero info." );

  lapack.GETRI(card,
      vmat.data(), vmat.stride_1(),
      (ordinal_type*)ipiv.data(),
      work.data(), lwork,
      &info);

  INTREPID2_TEST_FOR_EXCEPTION( info != 0,
      std::runtime_error ,
      ">>> ERROR: (Intrepid2::Basis_HDIV_TRI_In_FEM) lapack.GETRI returns nonzero info." );

  for (ordinal_type i=0;i<cardVecPn;++i)
    for (ordinal_type j=0;j<card;++j){
      scalarType s=0;
      for(ordinal_type k=0; k< card; ++k)
        s += V1(i,k)*vmat(k,j);
      coeffs(i,j) = s;
    }

  this->coeffs_ = Kokkos::create_mirror_view(typename SpT::memory_space(), coeffs);
  Kokkos::deep_copy(this->coeffs_ , coeffs);

  this->dofCoords_ = Kokkos::create_mirror_view(typename SpT::memory_space(), dofCoords);
  Kokkos::deep_copy(this->dofCoords_, dofCoords);

  this->dofCoeffs_ = Kokkos::create_mirror_view(typename SpT::memory_space(), dofCoeffs);
  Kokkos::deep_copy(this->dofCoeffs_, dofCoeffs);


  // set tags
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
