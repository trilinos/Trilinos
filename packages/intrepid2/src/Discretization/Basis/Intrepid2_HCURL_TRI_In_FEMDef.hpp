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

/** \file   Intrepid2_HCURL_TRI_In_FEMDef.hpp
    \brief  Definition file for FEM basis functions of degree n for H(curl) functions on TRI.
    \author Created by R. Kirby and P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_HCURL_TRI_IN_FEM_DEF_HPP__
#define __INTREPID2_HCURL_TRI_IN_FEM_DEF_HPP__

#include "Intrepid2_HGRAD_TRI_Cn_FEM_ORTH.hpp"
#include "Intrepid2_CubatureDirectTriDefault.hpp"

namespace Intrepid2 {

  // -------------------------------------------------------------------------------------

  namespace Impl {

    template<EOperator opType>
    template<typename OutputViewType,
             typename inputViewType,
             typename workViewType,
             typename vinvViewType>
    KOKKOS_INLINE_FUNCTION
    void
    Basis_HCURL_TRI_In_FEM::Serial<opType>::
    getValues(       OutputViewType output,
               const inputViewType  input,
                     workViewType   work,
               const vinvViewType   coeffs ) {

      constexpr ordinal_type spaceDim = 2;
      const ordinal_type
        cardPn = coeffs.extent(0)/spaceDim,
        card = coeffs.extent(1),
        npts = input.extent(0);

      // compute order
      ordinal_type order = 0;
      for (ordinal_type p=0;p<=Parameters::MaxOrder;++p) {
        if (card == CardinalityHCurlTri(p)) {
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
                output.access(i,j,d) += coeffs(k+d*cardPn,i) * phis(k,j);
            }
        break;
      }
      case OPERATOR_CURL: {
        const viewType phis(Kokkos::view_wrap(ptr, vcprop), card, npts, spaceDim);
        ptr += card*npts*spaceDim*get_dimension_scalar(work);
        const viewType workView(Kokkos::view_wrap(ptr, vcprop), card, npts, spaceDim+1);

        Impl::Basis_HGRAD_TRI_Cn_FEM_ORTH::
          Serial<OPERATOR_GRAD>::getValues(phis, input, workView, order);

        for (ordinal_type i=0;i<card;++i)
          for (ordinal_type j=0;j<npts;++j) {
            output.access(i,j) = 0.0;
            for (ordinal_type k=0; k<cardPn; ++k)
              output.access(i,j) += - coeffs(k,i)*phis(k,j,1)              // - dy of x component
                + coeffs(k+cardPn,i)*phis(k,j,0);      // dx of y component
          }
        break;
      }
      default: {
        INTREPID2_TEST_FOR_ABORT( true,
                                  ">>> ERROR (Basis_HCURL_TRI_In_FEM): Operator type not implemented");
      }
      }
    }

    template<typename DT, ordinal_type numPtsPerEval,
             typename outputValueValueType, class ...outputValueProperties,
             typename inputPointValueType,  class ...inputPointProperties,
             typename vinvValueType,        class ...vinvProperties>
    void
    Basis_HCURL_TRI_In_FEM::
    getValues(       Kokkos::DynRankView<outputValueValueType,outputValueProperties...> outputValues,
               const Kokkos::DynRankView<inputPointValueType, inputPointProperties...>  inputPoints,
               const Kokkos::DynRankView<vinvValueType,       vinvProperties...>        coeffs,
               const EOperator operatorType) {
      typedef          Kokkos::DynRankView<outputValueValueType,outputValueProperties...>         outputValueViewType;
      typedef          Kokkos::DynRankView<inputPointValueType, inputPointProperties...>          inputPointViewType;
      typedef          Kokkos::DynRankView<vinvValueType,       vinvProperties...>                vinvViewType;
      typedef typename ExecSpace<typename inputPointViewType::execution_space,typename DT::execution_space>::ExecSpaceType ExecSpaceType;

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
        workViewType  work(Kokkos::view_alloc("Basis_HCURL_TRI_In_FEM::getValues::work", vcprop), cardinality, inputPoints.extent(0));
        typedef Functor<outputValueViewType,inputPointViewType,vinvViewType, workViewType,
          OPERATOR_VALUE,numPtsPerEval> FunctorType;
        Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints, coeffs, work) );
        break;
      }
      case OPERATOR_CURL: {
        workViewType  work(Kokkos::view_alloc("Basis_HCURL_TRI_In_FEM::getValues::work", vcprop), cardinality*(2*spaceDim+1), inputPoints.extent(0));
        typedef Functor<outputValueViewType,inputPointViewType,vinvViewType, workViewType,
          OPERATOR_CURL,numPtsPerEval> FunctorType;
        Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints, coeffs, work) );
        break;
      }
      default: {
        INTREPID2_TEST_FOR_EXCEPTION( true , std::invalid_argument,
                                      ">>> ERROR (Basis_HCURL_TRI_In_FEM): Operator type not implemented" );
      }
      }
    }
  }

  // -------------------------------------------------------------------------------------
  template<typename DT, typename OT, typename PT>
  Basis_HCURL_TRI_In_FEM<DT,OT,PT>::
  Basis_HCURL_TRI_In_FEM( const ordinal_type order,
                          const EPointType   pointType ) {

    constexpr ordinal_type spaceDim = 2;
    this->basisCardinality_  = CardinalityHCurlTri(order);
    this->basisDegree_       = order; // small n
    this->basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Triangle<3> >() );
    this->basisType_         = BASIS_FEM_LAGRANGIAN;
    this->basisCoordinates_  = COORDINATES_CARTESIAN;
    this->functionSpace_     = FUNCTION_SPACE_HCURL;
    pointType_ = pointType;

    const ordinal_type card = this->basisCardinality_;

    const ordinal_type  cardPn = Intrepid2::getPnCardinality<spaceDim>(order); // dim of (P_{n}) -- smaller space
    const ordinal_type  cardPnm1 = Intrepid2::getPnCardinality<spaceDim>(order-1);  // dim of (P_{n-1}) -- smaller space
    const ordinal_type  cardPnm2 = Intrepid2::getPnCardinality<spaceDim>(order-2); // dim of (P_{n-2}) -- smaller space
    const ordinal_type  cardVecPn = spaceDim*cardPn;  // dim of (P_{n})^2 -- larger space
    const ordinal_type  cardVecPnm1 = spaceDim*cardPnm1;   // dim of (P_{n-1})^2 -- smaller space


    // Basis-dependent initializations
    constexpr ordinal_type tagSize  = 4;        // size of DoF tag, i.e., number of fields in the tag
    constexpr ordinal_type maxCard = CardinalityHCurlTri(Parameters::MaxOrder);
    ordinal_type tags[maxCard][tagSize];

    // points are computed in the host and will be copied
    Kokkos::DynRankView<scalarType,typename DT::execution_space::array_layout,Kokkos::HostSpace>
      dofCoords("Hcurl::Tri::In::dofCoords", card, spaceDim);

    Kokkos::DynRankView<scalarType,typename DT::execution_space::array_layout,Kokkos::HostSpace>
      coeffs("Hcurl::Tri::In::coeffs", cardVecPn, card);

    Kokkos::DynRankView<scalarType,typename DT::execution_space::array_layout,Kokkos::HostSpace>
      dofCoeffs("Hcurl::Tri::In::dofCoeffs", card, spaceDim);

    // first, need to project the basis for RT space onto the
    // orthogonal basis of degree n
    // get coefficients of PkHx

    const ordinal_type lwork = card*card;
    Kokkos::DynRankView<scalarType,typename DT::execution_space::array_layout,Kokkos::HostSpace>
      V1("Hcurl::Tri::In::V1", cardVecPn, card);

    // basis for the space is
    // { (phi_i,0) }_{i=0}^{cardPnm1-1} ,
    // { (0,phi_i) }_{i=0}^{cardPnm1-1} ,
    // { (x,y) \times phi_i}_{i=cardPnm2}^{cardPnm1-1}
    // { (x,y) \times phi = (y phi , -x \phi)
    // columns of V1 are expansion of this basis in terms of the basis
    // for P_{n}^2

    // these two loops get the first two sets of basis functions
    for (ordinal_type i=0;i<cardPnm1;i++)
      for (ordinal_type d=0;d<spaceDim;d++)
        V1(d*cardPn+i,d*cardPnm1+i) = 1.0;


    // now I need to integrate { (x,y) \times phi } against the big basis
    // first, get a cubature rule.
    CubatureDirectTriDefault<Kokkos::HostSpace::execution_space,scalarType,scalarType> myCub( 2 * order );
    Kokkos::DynRankView<scalarType,typename DT::execution_space::array_layout,Kokkos::HostSpace> cubPoints("Hcurl::Tri::In::cubPoints", myCub.getNumPoints() , spaceDim );
    Kokkos::DynRankView<scalarType,typename DT::execution_space::array_layout,Kokkos::HostSpace> cubWeights("Hcurl::Tri::In::cubWeights", myCub.getNumPoints() );
    myCub.getCubature( cubPoints , cubWeights );

    // tabulate the scalar orthonormal basis at cubature points
    Kokkos::DynRankView<scalarType,typename DT::execution_space::array_layout,Kokkos::HostSpace> phisAtCubPoints("Hcurl::Tri::In::phisAtCubPoints", cardPn , myCub.getNumPoints() );
    Impl::Basis_HGRAD_TRI_Cn_FEM_ORTH::getValues<Kokkos::HostSpace::execution_space,Parameters::MaxNumPtsPerBasisEval>(phisAtCubPoints, cubPoints, order, OPERATOR_VALUE);

    // now do the integration
    for (ordinal_type i=0;i<order;i++) {
      for (ordinal_type j=0;j<cardPn;j++) { // int (x,y) phi_i \cdot (phi_j,phi_{j+cardPn})
        for (ordinal_type k=0;k<myCub.getNumPoints();k++) {
          V1(j,cardVecPnm1+i) -=
            cubWeights(k) * cubPoints(k,1)
            * phisAtCubPoints(cardPnm2+i,k)
            * phisAtCubPoints(j,k);
          V1(j+cardPn,cardVecPnm1+i) +=
            cubWeights(k) * cubPoints(k,0)
            * phisAtCubPoints(cardPnm2+i,k)
            * phisAtCubPoints(j,k);
        }
      }
    }

    // next, apply the RT nodes (rows) to the basis for (P_n)^2 (columns)
    Kokkos::DynRankView<scalarType,typename DT::execution_space::array_layout,Kokkos::HostSpace>
      V2("Hcurl::Tri::In::V2", card ,cardVecPn);

    const ordinal_type numEdges = this->basisCellTopology_.getEdgeCount();

    shards::CellTopology edgeTop(shards::getCellTopologyData<shards::Line<2> >() );

    const int numPtsPerEdge = PointTools::getLatticeSize( edgeTop ,
                                                          order+1 ,
                                                          1 );

    // first numEdges * degree nodes are tangents at each edge
    // get the points on the line
    Kokkos::DynRankView<scalarType,typename DT::execution_space::array_layout,Kokkos::HostSpace> linePts("Hcurl::Tri::In::linePts", numPtsPerEdge , 1 );

    // construct lattice
    const ordinal_type offset = 1;
    PointTools::getLattice( linePts,
                            edgeTop,
                            order+1, offset,
                            pointType );

    // holds the image of the line points
    Kokkos::DynRankView<scalarType,typename DT::execution_space::array_layout,Kokkos::HostSpace> edgePts("Hcurl::Tri::In::edgePts", numPtsPerEdge , spaceDim );
    Kokkos::DynRankView<scalarType,typename DT::execution_space::array_layout,Kokkos::HostSpace> phisAtEdgePoints("Hcurl::Tri::In::phisAtEdgePoints", cardPn , numPtsPerEdge );
    Kokkos::DynRankView<scalarType,typename DT::execution_space::array_layout,Kokkos::HostSpace> edgeTan("Hcurl::Tri::In::edgeTan", spaceDim );

    // these are tangents scaled by the appropriate edge lengths.
    for (ordinal_type edge=0;edge<numEdges;edge++) {  // loop over edges
      CellTools<Kokkos::HostSpace>::getReferenceEdgeTangent( edgeTan ,
                                                                              edge ,
                                                                              this->basisCellTopology_ );

      CellTools<Kokkos::HostSpace>::mapToReferenceSubcell( edgePts ,
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
          V2(i_card,k) = edgeTan(0) * phisAtEdgePoints(k,j);
          V2(i_card,k+cardPn) = edgeTan(1) * phisAtEdgePoints(k,j);
        }


        //save dof coordinates
        for(ordinal_type k=0; k<spaceDim; ++k) {
          dofCoords(i_card,k) = edgePts(j,k);
          dofCoeffs(i_card,k) = edgeTan(k);
        }

        tags[i_card][0] = 1; // edge dof
        tags[i_card][1] = edge; // edge id
        tags[i_card][2] = j; // local dof id
        tags[i_card][3] = numPtsPerEdge; // total edge dof

      }


    }

    // remaining nodes are x- and y- components at internal points (this code is same as HDIV).
    //These are evaluated at the interior of a lattice of degree + 1, For then
    // the degree == 1 space corresponds classicaly to RT0 and so gets
    // no internal nodes, and degree == 2 corresponds to RT1 and needs
    // one internal node per vector component.
    const ordinal_type numPtsPerCell = PointTools::getLatticeSize( this->basisCellTopology_ ,
                                                                   order + 1 ,
                                                                   1 );

    if (numPtsPerCell > 0) {
      Kokkos::DynRankView<scalarType,typename DT::execution_space::array_layout,Kokkos::HostSpace>
        internalPoints( "Hcurl::Tri::In::internalPoints", numPtsPerCell , spaceDim );
      PointTools::getLattice( internalPoints ,
                              this->basisCellTopology_ ,
                              order + 1 ,
                              1 ,
                              pointType );

      Kokkos::DynRankView<scalarType,typename DT::execution_space::array_layout,Kokkos::HostSpace>
        phisAtInternalPoints("Hcurl::Tri::In::phisAtInternalPoints", cardPn , numPtsPerCell );
      Impl::Basis_HGRAD_TRI_Cn_FEM_ORTH::getValues<Kokkos::HostSpace::execution_space,Parameters::MaxNumPtsPerBasisEval>( phisAtInternalPoints , internalPoints , order, OPERATOR_VALUE );

      // copy values into right positions of V2
      for (ordinal_type j=0;j<numPtsPerCell;j++) {

        const ordinal_type i_card = numEdges*order+spaceDim*j;

        for (ordinal_type k=0;k<cardPn;k++) {
          // x component
          V2(i_card,k) = phisAtInternalPoints(k,j);
          // y component
          V2(i_card+1,cardPn+k) = phisAtInternalPoints(k,j);
        }

        //save dof coordinates
        for(ordinal_type d=0; d<spaceDim; ++d) {
          for(ordinal_type dim=0; dim<spaceDim; ++dim) {
            dofCoords(i_card+d,dim) = internalPoints(j,dim);
            dofCoeffs(i_card+d,dim) = (d==dim);
          }

          tags[i_card+d][0] = spaceDim; // elem dof
          tags[i_card+d][1] = 0; // elem id
          tags[i_card+d][2] = spaceDim*j+d; // local dof id
          tags[i_card+d][3] = spaceDim*numPtsPerCell; // total vert dof
        }
      }
    }

    // form Vandermonde matrix.  Actually, this is the transpose of the VDM,
    // so we transpose on copy below.
    Kokkos::DynRankView<scalarType,Kokkos::LayoutLeft,Kokkos::HostSpace>
      vmat("Hcurl::Tri::In::vmat", card, card),
      work("Hcurl::Tri::In::work", lwork),
      ipiv("Hcurl::Tri::In::ipiv", card);

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
                                  ">>> ERROR: (Intrepid2::Basis_HCURL_TRI_In_FEM) lapack.GETRF returns nonzero info." );

    lapack.GETRI(card,
                 vmat.data(), vmat.stride_1(),
                 (ordinal_type*)ipiv.data(),
                 work.data(), lwork,
                 &info);

    INTREPID2_TEST_FOR_EXCEPTION( info != 0,
                                  std::runtime_error ,
                                  ">>> ERROR: (Intrepid2::Basis_HCURL_TRI_In_FEM) lapack.GETRI returns nonzero info." );

    for (ordinal_type i=0;i<cardVecPn;++i)
      for (ordinal_type j=0;j<card;++j){
        scalarType s=0;
        for(ordinal_type k=0; k< card; ++k)
          s += V1(i,k)*vmat(k,j);
        coeffs(i,j) = s;
      }

    this->coeffs_ = Kokkos::create_mirror_view(typename DT::memory_space(), coeffs);
    Kokkos::deep_copy(this->coeffs_ , coeffs);

    this->dofCoords_ = Kokkos::create_mirror_view(typename DT::memory_space(), dofCoords);
    Kokkos::deep_copy(this->dofCoords_, dofCoords);

    this->dofCoeffs_ = Kokkos::create_mirror_view(typename DT::memory_space(), dofCoeffs);
    Kokkos::deep_copy(this->dofCoeffs_, dofCoeffs);


    // set tags
    {
      // Basis-dependent initializations
      const ordinal_type posScDim = 0;        // position in the tag, counting from 0, of the subcell dim
      const ordinal_type posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
      const ordinal_type posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell

      OrdinalTypeArray1DHost tagView(&tags[0][0], card*tagSize);

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
