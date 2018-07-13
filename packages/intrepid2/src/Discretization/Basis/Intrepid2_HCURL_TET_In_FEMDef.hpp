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

/** \file   Intrepid2_HCURL_TET_In_FEMDef.hpp
    \brief  Definition file for FEM basis functions of degree n for H(curl) functions on TET.
    \author Created by R. Kirby and P. Bochev and D. Ridzal.
    Kokkorized by Kyungjoo Kim
 */

#ifndef __INTREPID2_HCURL_TET_IN_FEM_DEF_HPP__
#define __INTREPID2_HCURL_TET_IN_FEM_DEF_HPP__

#include "Intrepid2_HGRAD_TET_Cn_FEM_ORTH.hpp"
#include "Intrepid2_CubatureDirectTetDefault.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"

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
Basis_HCURL_TET_In_FEM::Serial<opType>::
getValues(       outputViewType output,
    const inputViewType  input,
          workViewType   work,
    const vinvViewType   coeffs ) {

  constexpr ordinal_type spaceDim = 3;
  const ordinal_type
  cardPn = coeffs.extent(0)/spaceDim,
  card = coeffs.extent(1),
  npts = input.extent(0);

  // compute order
  ordinal_type order = 0;
  for (ordinal_type p=0;p<=Parameters::MaxOrder;++p) {
    if (card == CardinalityHCurlTet(p)) {
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

    Impl::Basis_HGRAD_TET_Cn_FEM_ORTH::
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

    Impl::Basis_HGRAD_TET_Cn_FEM_ORTH::
    Serial<OPERATOR_GRAD>::getValues(phis, input, workView, order);

    for (ordinal_type i=0;i<card;++i) {
      for (ordinal_type j=0;j<npts;++j) {
        for (ordinal_type d=0; d< spaceDim; ++d) {
          output.access(i,j,d) = 0.0;
          ordinal_type d1 = (d+1) % spaceDim, d2 = (d+2) % spaceDim;
          for (ordinal_type k=0; k<cardPn; ++k)   //\sum_k (coeffs_k, coeffs_{k+cardPn}, coeffs_{k+2 cardPn}) \times phis_kj  (cross product)
            output.access(i,j,d) += coeffs(k+d2*cardPn,i)*phis(k,j,d1)
            -coeffs(k+d1*cardPn,i)*phis(k,j,d2);
        }
      }
    }
    break;
  }
  default: {
    INTREPID2_TEST_FOR_ABORT( true,
        ">>> ERROR (Basis_HCURL_TET_In_FEM): Operator type not implemented");
  }
  }
}

template<typename SpT, ordinal_type numPtsPerEval,
typename outputValueValueType, class ...outputValueProperties,
typename inputPointValueType,  class ...inputPointProperties,
typename vinvValueType,        class ...vinvProperties>
void
Basis_HCURL_TET_In_FEM::
getValues(       Kokkos::DynRankView<outputValueValueType,outputValueProperties...> outputValues,
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
  const ordinal_type spaceDim = 3;

  auto vcprop = Kokkos::common_view_alloc_prop(inputPoints);
  typedef typename Kokkos::DynRankView< inputPointType, typename inputPointViewType::memory_space> workViewType;

  switch (operatorType) {
  case OPERATOR_VALUE: {
    workViewType  work(Kokkos::view_alloc("Basis_HCURL_TET_In_FEM::getValues::work", vcprop), cardinality, inputPoints.extent(0));
    typedef Functor<outputValueViewType,inputPointViewType,vinvViewType, workViewType,
        OPERATOR_VALUE,numPtsPerEval> FunctorType;
    Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints, coeffs, work) );
    break;
  }
  case OPERATOR_CURL: {
    workViewType  work(Kokkos::view_alloc("Basis_HCURL_TET_In_FEM::getValues::work", vcprop), cardinality*(2*spaceDim+1), inputPoints.extent(0));
    typedef Functor<outputValueViewType,inputPointViewType,vinvViewType, workViewType,
        OPERATOR_CURL,numPtsPerEval> FunctorType;
    Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints, coeffs, work) );
    break;
  }
  default: {
    INTREPID2_TEST_FOR_EXCEPTION( true , std::invalid_argument,
        ">>> ERROR (Basis_HCURL_TET_In_FEM): Operator type not implemented" );
  }
  }
}
}

// -------------------------------------------------------------------------------------
template<typename SpT, typename OT, typename PT>
Basis_HCURL_TET_In_FEM<SpT,OT,PT>::
Basis_HCURL_TET_In_FEM( const ordinal_type order,
    const EPointType   pointType ) {

  constexpr ordinal_type spaceDim = 3;
  this->basisCardinality_  = CardinalityHCurlTet(order);
  this->basisDegree_       = order; // small n
  this->basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Tetrahedron<4> >() );
  this->basisType_         = BASIS_FEM_FIAT;
  this->basisCoordinates_  = COORDINATES_CARTESIAN;

  const ordinal_type card = this->basisCardinality_;

  const ordinal_type  cardPn = Intrepid2::getPnCardinality<spaceDim>(order); // dim of (P_{n}) -- smaller space
  const ordinal_type  cardPnm1 = Intrepid2::getPnCardinality<spaceDim>(order-1);  // dim of (P_{n-1}) -- smaller space
  const ordinal_type  cardPnm2 = Intrepid2::getPnCardinality<spaceDim>(order-2); // dim of (P_{n-2}) -- smaller space
  const ordinal_type  cardVecPn = spaceDim*cardPn;  // dim of (P_{n})^2 -- larger space
  const ordinal_type  cardVecPnm1 = spaceDim*cardPnm1;   // dim of (P_{n-1})^2 -- smaller space
  const ordinal_type  cardPnm1H = cardPnm1-cardPnm2; //Homogeneous polynomial of order (n-1)



  // Basis-dependent initializations
  constexpr ordinal_type tagSize  = 4;        // size of DoF tag, i.e., number of fields in the tag
  constexpr ordinal_type maxCard = CardinalityHCurlTet(Parameters::MaxOrder);
  ordinal_type tags[maxCard][tagSize];

  // points are computed in the host and will be copied
  Kokkos::DynRankView<scalarType,typename SpT::array_layout,Kokkos::HostSpace>
  dofCoords("Hcurl::Tet::In::dofCoords", card, spaceDim);

  Kokkos::DynRankView<scalarType,typename SpT::array_layout,Kokkos::HostSpace>
  coeffs("Hcurl::Tet::In::coeffs", cardVecPn, card);

  Kokkos::DynRankView<scalarType,typename SpT::array_layout,Kokkos::HostSpace>
  dofCoeffs("Hcurl::Tet::In::dofCoeffs", card, spaceDim);

  // first, need to project the basis for RT space onto the
  // orthogonal basis of degree n
  // get coefficients of PkHx

  Kokkos::DynRankView<scalarType,Kokkos::LayoutLeft,Kokkos::HostSpace> //use LayoutLeft for Lapack
  V1("Hcurl::Tet::In::V1", cardVecPn, cardVecPnm1 + spaceDim*cardPnm1H);


  // these two loops get the first three sets of basis functions
  for (ordinal_type i=0;i<cardPnm1;i++)
    for (ordinal_type d=0;d<spaceDim;d++)
      V1(i+d*cardPn,i+d*cardPnm1) = 1.0;


  // now I need to integrate { (x,y) \times phi } against the big basis
  // first, get a cubature rule.
  CubatureDirectTetDefault<Kokkos::HostSpace::execution_space,scalarType,scalarType> myCub( 2 * order );
  Kokkos::DynRankView<scalarType,typename SpT::array_layout,Kokkos::HostSpace> cubPoints("Hcurl::Tet::In::cubPoints", myCub.getNumPoints() , spaceDim );
  Kokkos::DynRankView<scalarType,typename SpT::array_layout,Kokkos::HostSpace> cubWeights("Hcurl::Tet::In::cubWeights", myCub.getNumPoints() );
  myCub.getCubature( cubPoints , cubWeights );

  // tabulate the scalar orthonormal basis at cubature points
  Kokkos::DynRankView<scalarType,typename SpT::array_layout,Kokkos::HostSpace> phisAtCubPoints("Hcurl::Tet::In::phisAtCubPoints", cardPn , myCub.getNumPoints() );
  Impl::Basis_HGRAD_TET_Cn_FEM_ORTH::getValues<Kokkos::HostSpace::execution_space,Parameters::MaxNumPtsPerBasisEval>(phisAtCubPoints, cubPoints, order, OPERATOR_VALUE);

  // Integrate (x psi_j, y psi_j, z psi_j) \times (phi_i, phi_{i+cardPn}, phi_{i+2 cardPn})    cross product. psi are homogeneous polynomials of order (n-1)
  for (ordinal_type i=0;i<cardPn;i++) {
    for (ordinal_type j=0;j<cardPnm1H;j++) { // loop over homogeneous polynomials
      for (ordinal_type d=0; d< spaceDim; ++d) {
        scalarType integral(0);
        for (ordinal_type k=0;k<myCub.getNumPoints();k++)
          integral += cubWeights(k) * cubPoints(k,d)
          * phisAtCubPoints(cardPnm2+j,k)
          * phisAtCubPoints(i,k);
        ordinal_type d1 = (d+1) % spaceDim, d2 = (d+2) % spaceDim;
        V1(i+d2*cardPn,cardVecPnm1+d1*cardPnm1H + j) = -integral;
        V1(i+d1*cardPn,cardVecPnm1+d2*cardPnm1H + j) = integral;
      }
    }
  }





  // now I need to set up an SVD to get a basis for the space
  Kokkos::DynRankView<scalarType,Kokkos::LayoutLeft,Kokkos::HostSpace>
  S("Hcurl::Tet::In::S", cardVecPn,1),
  U("Hcurl::Tet::In::U", cardVecPn, cardVecPn),
  Vt("Hcurl::Tet::In::Vt", cardVecPn, cardVecPn),
  work("Hcurl::Tet::In::work", 5*cardVecPn,1),
  rWork("Hcurl::Tet::In::rW", 1,1);



  ordinal_type info = 0;
  Teuchos::LAPACK<ordinal_type,scalarType> lapack;


  lapack.GESVD( 'A',
      'N',
      V1.extent(0) ,
      V1.extent(1) ,
      V1.data() ,
      V1.stride_1() ,
      S.data() ,
      U.data() ,
      U.stride_1() ,
      Vt.data() ,
      Vt.stride_1() ,
      work.data() ,
      5*cardVecPn ,
      rWork.data() ,
      &info );


#ifdef HAVE_INTREPID2_DEBUG
  ordinal_type num_nonzero_sv = 0;
  for (int i=0;i<cardVecPn;i++)
    num_nonzero_sv += (S(i,0) > tolerence());

  INTREPID2_TEST_FOR_EXCEPTION( num_nonzero_sv != card, std::invalid_argument,
      ">>> ERROR: (Intrepid2::Basis_HCURL_TET_In_FEM( order, pointType), Matrix V1 should have rank equal to the cardinality of HCURL space");
#endif

  // next, apply the RT nodes (rows) to the basis for (P_n)^2 (columns)
  Kokkos::DynRankView<scalarType,typename SpT::array_layout,Kokkos::HostSpace>
  V2("Hcurl::Tet::In::V2", card ,cardVecPn);

  const ordinal_type numEdges = this->basisCellTopology_.getEdgeCount();
  const ordinal_type numFaces = this->basisCellTopology_.getFaceCount();

  // first numEdges * degree nodes are normals at each edge
  // get the points on the line

  shards::CellTopology edgeTop(shards::getCellTopologyData<shards::Line<2> >() );
  shards::CellTopology faceTop(shards::getCellTopologyData<shards::Triangle<3> >() );

  const int numPtsPerEdge = PointTools::getLatticeSize( edgeTop ,
      order+1 ,
      1 );

  const int numPtsPerFace = PointTools::getLatticeSize( faceTop ,
      order+1 ,
      1 );

  const int numPtsPerCell = PointTools::getLatticeSize( this->basisCellTopology_ ,
      order+1 ,
      1 );

  Kokkos::DynRankView<scalarType,typename SpT::array_layout,Kokkos::HostSpace> linePts("Hcurl::Tet::In::linePts", numPtsPerEdge , 1 );
  Kokkos::DynRankView<scalarType,typename SpT::array_layout,Kokkos::HostSpace> triPts("Hcurl::Tet::In::triPts", numPtsPerFace , 2 );

  // construct lattice
  const ordinal_type offset = 1;



  PointTools::getLattice( linePts,
      edgeTop,
      order+1, offset,
      pointType );

  PointTools::getLattice( triPts,
      faceTop,
      order+1, offset,
      pointType );

  // holds the image of the line points
  Kokkos::DynRankView<scalarType,typename SpT::array_layout,Kokkos::HostSpace> edgePts("Hcurl::Tet::In::edgePts", numPtsPerEdge , spaceDim );
  Kokkos::DynRankView<scalarType,typename SpT::array_layout,Kokkos::HostSpace> facePts("Hcurl::Tet::In::facePts", numPtsPerFace , spaceDim );
  Kokkos::DynRankView<scalarType,typename SpT::array_layout,Kokkos::HostSpace> phisAtEdgePoints("Hcurl::Tet::In::phisAtEdgePoints", cardPn , numPtsPerEdge );
  Kokkos::DynRankView<scalarType,typename SpT::array_layout,Kokkos::HostSpace> phisAtFacePoints("Hcurl::Tet::In::phisAtFacePoints", cardPn , numPtsPerFace);

  Kokkos::DynRankView<scalarType,typename SpT::array_layout,Kokkos::HostSpace> edgeTan("Hcurl::Tet::In::edgeTan", spaceDim );

  // these are tangents scaled by the appropriate edge lengths.
  for (ordinal_type i=0;i<numEdges;i++) {  // loop over edges
    CellTools<Kokkos::HostSpace::execution_space>::getReferenceEdgeTangent( edgeTan ,
        i ,
        this->basisCellTopology_ );
    /* multiply by measure of reference edge so that magnitude of the edgeTan is equal to the edge measure */
    const scalarType refEdgeMeasure = 2.0;
    for (ordinal_type j=0;j<spaceDim;j++)
      edgeTan(j) *= refEdgeMeasure;

    CellTools<Kokkos::HostSpace::execution_space>::mapToReferenceSubcell( edgePts ,
        linePts ,
        1 ,
        i ,
        this->basisCellTopology_ );

    Impl::Basis_HGRAD_TET_Cn_FEM_ORTH::getValues<Kokkos::HostSpace::execution_space,Parameters::MaxNumPtsPerBasisEval>(phisAtEdgePoints , edgePts, order, OPERATOR_VALUE);

    // loop over points (rows of V2)
    for (ordinal_type j=0;j<numPtsPerEdge;j++) {

      const ordinal_type i_card = numPtsPerEdge*i+j;

      // loop over orthonormal basis functions (columns of V2)
      for (ordinal_type k=0;k<cardPn;k++)
        for (ordinal_type d=0;d<spaceDim;d++)
          V2(i_card,k+d*cardPn) = edgeTan(d) * phisAtEdgePoints(k,j);

      //save dof coordinates and coefficients
      for(ordinal_type k=0; k<spaceDim; ++k) {
        dofCoords(i_card,k) = edgePts(j,k);
        dofCoeffs(i_card,k) = edgeTan(k);
      }

      tags[i_card][0] = 1; // edge dof
      tags[i_card][1] = i; // edge id
      tags[i_card][2] = j; // local dof id
      tags[i_card][3] = numPtsPerEdge; // total vert dof

    }
  }

  if(numPtsPerFace >0) {//handle faces if needed  (order >1)
    Kokkos::DynRankView<scalarType,typename SpT::array_layout,Kokkos::HostSpace> faceTan1("Hcurl::Tet::In::edgeTan", spaceDim );
    Kokkos::DynRankView<scalarType,typename SpT::array_layout,Kokkos::HostSpace> faceTan2("Hcurl::Tet::In::edgeTan", spaceDim );

    for (ordinal_type i=0;i<numFaces;i++) {  // loop over faces
      CellTools<Kokkos::HostSpace::execution_space>::getReferenceFaceTangents( faceTan1 ,
          faceTan2,
          i ,
          this->basisCellTopology_ );

      CellTools<Kokkos::HostSpace::execution_space>::mapToReferenceSubcell( facePts ,
          triPts ,
          2 ,
          i ,
          this->basisCellTopology_ );

      Impl::Basis_HGRAD_TET_Cn_FEM_ORTH::getValues<Kokkos::HostSpace::execution_space,Parameters::MaxNumPtsPerBasisEval>(phisAtFacePoints , facePts, order, OPERATOR_VALUE);

      // loop over points (rows of V2)
      for (ordinal_type j=0;j<numPtsPerFace;j++) {

        const ordinal_type i_card = numEdges*numPtsPerEdge+2*numPtsPerFace*i+2*j;

        // loop over orthonormal basis functions (columns of V2)
        for (ordinal_type k=0;k<cardPn;k++)
          for (ordinal_type d=0;d<spaceDim;d++)  {
            V2(i_card,k+d*cardPn) = faceTan1(d) * phisAtFacePoints(k,j);
            V2(i_card+1,k+d*cardPn) = faceTan2(d) * phisAtFacePoints(k,j);
          }

        //save dof coordinates
        for(ordinal_type k=0; k<spaceDim; ++k) {
          dofCoords(i_card,k) = facePts(j,k);
          dofCoords(i_card+1,k) = facePts(j,k);
          dofCoeffs(i_card,k) = faceTan1(k);
          dofCoeffs(i_card+1,k) = faceTan2(k);
        }

        tags[i_card][0] = 2; // face dof
        tags[i_card][1] = i; // face id
        tags[i_card][2] = 2*j; // local face id
        tags[i_card][3] = 2*numPtsPerFace; // total face dof

        tags[i_card+1][0] = 2; // face dof
        tags[i_card+1][1] = i; // face id
        tags[i_card+1][2] = 2*j+1; // local face id
        tags[i_card+1][3] = 2*numPtsPerFace; // total face dof

      }
    }
  }


  // internal dof, if needed
  if (numPtsPerCell > 0) {
    Kokkos::DynRankView<scalarType,typename SpT::array_layout,Kokkos::HostSpace>
    cellPoints( "Hcurl::Tet::In::cellPoints", numPtsPerCell , spaceDim );
    PointTools::getLattice( cellPoints ,
        this->basisCellTopology_ ,
        order + 1 ,
        1 ,
        pointType );

    Kokkos::DynRankView<scalarType,typename SpT::array_layout,Kokkos::HostSpace>
    phisAtCellPoints("Hcurl::Tet::In::phisAtCellPoints", cardPn , numPtsPerCell );
    Impl::Basis_HGRAD_TET_Cn_FEM_ORTH::getValues<Kokkos::HostSpace::execution_space,Parameters::MaxNumPtsPerBasisEval>( phisAtCellPoints , cellPoints , order, OPERATOR_VALUE );

    // copy values into right positions of V2
    for (ordinal_type j=0;j<numPtsPerCell;j++) {

      const ordinal_type i_card = numEdges*numPtsPerEdge+2*numFaces*numPtsPerFace+j;

      for (ordinal_type k=0;k<cardPn;k++)
        for (ordinal_type d=0;d<spaceDim;d++)
          V2(i_card+d*numPtsPerCell,d*cardPn+k) = phisAtCellPoints(k,j);


      //save dof coordinates
      for(ordinal_type d=0; d<spaceDim; ++d) {
        for(ordinal_type dim=0; dim<spaceDim; ++dim) {
          dofCoords(i_card+d*numPtsPerCell,dim) = cellPoints(j,dim);
          dofCoeffs(i_card+d*numPtsPerCell,dim) = (d==dim);
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
  const ordinal_type lwork = card*card;
  Kokkos::DynRankView<scalarType,Kokkos::LayoutLeft,Kokkos::HostSpace>
  vmat("Hcurl::Tet::In::vmat", card, card),
  work1("Hcurl::Tet::In::work", lwork),
  ipiv("Hcurl::Tet::In::ipiv", card);

  //vmat = V2*U;
  for(ordinal_type i=0; i< card; ++i) {
    for(ordinal_type j=0; j< card; ++j) {
      scalarType s=0;
      for(ordinal_type k=0; k< cardVecPn; ++k)
        s += V2(i,k)*U(k,j);
      vmat(i,j) = s;
    }
  }

  info = 0;

  lapack.GETRF(card, card,
      vmat.data(), vmat.stride_1(),
      (ordinal_type*)ipiv.data(),
      &info);

  INTREPID2_TEST_FOR_EXCEPTION( info != 0,
      std::runtime_error ,
      ">>> ERROR: (Intrepid2::Basis_HCURL_TET_In_FEM) lapack.GETRF returns nonzero info." );

  lapack.GETRI(card,
      vmat.data(), vmat.stride_1(),
      (ordinal_type*)ipiv.data(),
      work1.data(), lwork,
      &info);

  INTREPID2_TEST_FOR_EXCEPTION( info != 0,
      std::runtime_error ,
      ">>> ERROR: (Intrepid2::Basis_HCURL_TET_In_FEM) lapack.GETRI returns nonzero info." );

  for (ordinal_type i=0;i<cardVecPn;++i) {
    for (ordinal_type j=0;j<card;++j){
      scalarType s=0;
      for(ordinal_type k=0; k< card; ++k)
        s += U(i,k)*vmat(k,j);
      coeffs(i,j) = s;
    }
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
