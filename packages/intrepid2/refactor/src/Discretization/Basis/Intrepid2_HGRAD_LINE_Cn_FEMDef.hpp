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

/** \file   Intrepid_HGRAD_LINE_Cn_FEM_Def.hpp
    \brief  Definition file for FEM basis functions of degree n for H(grad) functions on LINE.
    \author Created by R. Kirby and P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_HGRAD_LINE_CN_FEM_DEF_HPP__
#define __INTREPID2_HGRAD_LINE_CN_FEM_DEF_HPP__


namespace Intrepid2 {

  // -------------------------------------------------------------------------------------
  template<typename SpT, typename OT, typename PT>
  template<EOperator opType>
  template<typename outputValueType, class ...outputProperties,
           typename phisValueType,   class ...phisProperties,
           typename vinvValueType,   class ...vinvProperties>
  KOKKOS_INLINE_FUNCTION
  void
  Basis_HGRAD_LINE_Cn_FEM<SpT,OT,PT>::Serial<opType>::
  getValues( /**/  Kokkos::DynRankView<outputValueType,outputProperties...>  output,
             const Kokkos::DynRankView<phisValueType,phisProperties...>      phis,
             const Kokkos::DynRankView<vinvValueType,vinvPointProperties...> vinv ){    
    switch (opType) {
    case OPERATOR_VALUE: {
      const auto card = phis.dimension(0);      
      const auto np   = phis.dimension(1);      
      for (auto i=0;i<np;++i) {
        output(i) = 0.0;
        for (auto j=0;j<card;++j) 
          output(i) += vinv(j)*phis(j,i);
      }
      break;
    }
    case OPERATOR_GRAD:
    case OPERATOR_D1:
    case OPERATOR_D2:
    case OPERATOR_D3:
    case OPERATOR_D4:
    case OPERATOR_D5:
    case OPERATOR_D6:
    case OPERATOR_D7:
    case OPERATOR_D8:
    case OPERATOR_D9:
    case OPERATOR_D10: {
      const auto card   = phis.dimension(0);      
      const auto np     = phis.dimension(1);      
      const auto dkcard = phis.dimension(2);      
      for (auto i=0;i<np;++i) 
        for (auto j=0;j<dkcard;++j) {
          output(i,j) = 0.0;
          for (auto k=0;k<card;++k) 
            output(i,j) += vinv(k)*phis(k,i,j);
        }
      break;
    }
    default: {
      INTREPID2_TEST_FOR_ABORT( true,
                                ">>> ERROR: (Intrepid2::Basis_HGRAD_LINE_Cn_FEM::Serial::getValues) operator is not supported." );
    }
    }
  }

  // -------------------------------------------------------------------------------------

  template<typename SpT, typename OT, typename PT>
  Basis_HGRAD_LINE_Cn_FEM<SpT,OT,PT>::
  Basis_HGRAD_LINE_Cn_FEM( const ordinal_type order,
                           const EPointType   pointType ) 
    : phis_(order) {
    this->basisCardinality_  = order+1;
    this->basisDegree_       = order;
    this->basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Line<2> >() );
    this->basisType_         = BASIS_FEM_FIAT;
    this->basisCoordinates_  = COORDINATES_CARTESIAN;

    const auto card = this->basisCardinality_;
    
    // points are computed in the host and will be copied 
    Kokkos::DynRankView<PT,typename SpT::array_layout,Kokkos::HostSpace>
      latticePts("Hgrad::Line::Cn::latticePts", card, 1);

    typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;
    switch (pointType) {
    case POINTTYPE_EQUISPACED:
    case POINTTYPE_WARPBLEND: {
      // two vertices
      latticePts(0,0) = -1.0;
      latticePts(1,0) =  1.0;

      // internal points
      auto pts = Kokkos::subdynrankview(latticePts, range_type(2, card), Kokkos::ALL());
      PointTools::getLattice( pts,
                              this->basisCellTopology_, 
                              order, 1, 
                              pointType );
      break;
    }
    case POINTTYPE_GUASS: {
      // internal points only
      PointTools::getGaussPoints( latticePts, 
                                  order );
      break;
    }
    default: {
      INTREPID2_TEST_FOR_EXCEPTION( !isValidPointType(pointType),
                                    std::invalid_argument , 
                                    ">>> ERROR: (Intrepid2::Basis_HGRAD_LINE_Cn_FEM) invalid pointType." );
    }
    }

    this->latticePts_ = Kokkos::create_mirror_view(typename SpT::memory_space(), latticePts);
    Kokkos::deep_copy(this->latticePts_, latticePts);
    
    // form Vandermonde matrix; actually, this is the transpose of the VDM,
    // this matrix is used in LAPACK so it should be column major and left layout
    const ordinal_type lwork = card*card;
    Kokkos::DynRankView<PT,Kokkkos::LeftLayout,Kokkos::HostSpace>
      Vmat("Hgrad::Line::Cn::Vmat", card, card), Work("Hgrad::Line::Cn::Work", lwork);

    Basis_HGRAD_LINE_Cn_FEM_JACOBI<Kokkos::HostSpace,OT,PT> phis(order);
    phis.impl_.getValues( Vmat, latticePts, OPERATOR_VALUE );

    ordinal_type ipiv[MaxOrder+1] = {}, info = 0;
    Teuchos::LAPACK<ordinal_type,PT> lapack;
    lapack.GETRI(card, 
                 Vmat.data(), Vmat.stride(0),
                 ipiv,
                 Work.data(), lsize,
                 &info);

    INTREPID2_TEST_FOR_EXCEPTION( info != 0,
                                  std::runtime_error , 
                                  ">>> ERROR: (Intrepid2::Basis_HGRAD_LINE_Cn_FEM) lapack.GETRI returns nonzero info." );
    
    // create host mirror 
    Kokkos::DynRankView<PT,typename SpT::array_layout,Kokkos::HostSpace>
      Vinv("Hgrad::Line::Cn::Vinv", card, card);

    for (auto j=0;j<card;++j)
      for (auto i=0;i<card;++i)
        Vinv(i,j) = Vmat(j,i);

    this->Vinv_ = Kokkos::create_mirror_view(typename SpT::memory_space(), Vinv);

    // initialize tags
    {
      const bool is_vertex_included = (pointType != POINTTYPE_GUASS)

      // Basis-dependent initializations
      const ordinal_type tagSize  = 4;        // size of DoF tag, i.e., number of fields in the tag
      const ordinal_type posScDim = 0;        // position in the tag, counting from 0, of the subcell dim 
      const ordinal_type posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
      const ordinal_type posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell
      
      ordinal_type tags[MaxOrder+1][4];

      // now we check the points for association 
      ordinal_type ibeg = 0;
      if (is_vertex_included) {
        tags[0][0] = 0; // vertex dof
        tags[0][1] = 0; // vertex id
        tags[0][2] = 0; // local dof id
        tags[0][3] = 1; // total number of dofs in this vertex

        tags[1][0] = 0; // vertex dof
        tags[1][1] = 1; // vertex id
        tags[1][2] = 0; // local dof id
        tags[1][3] = 1; // total number of dofs in this vertex

        const auto iend = card - 2;
        for (auto i=0;i<iend;++i) {
          tags[i][0] = 1;    // edge dof
          tags[i][1] = 0;    // edge id
          tags[i][2] = i;    // local dof id
          tags[i][3] = iend; // total number of dofs in this edge
        }
      } else {
        for (auto i=0;i<card;++i) {
          tags[i][0] = 1;    // edge dof
          tags[i][1] = 0;    // edge id
          tags[i][2] = i;    // local dof id
          tags[i][3] = card; // total number of dofs in this edge
        }
      }

      ordinal_type_array_1d_host tagView(&tag[0][0], card*4);

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
  
  
  template<typename SpT, typename OT, typename PT>
  template<typename outputValueValueType, class ...outputValueProperties,
           typename inputPointValueType,  class ...inputPointProperties>
  void
  Basis_HGRAD_LINE_Cn_FEM<SpT,OT,PT>::Internal::
  getValues( /**/  Kokkos::DynRankView<outputValueValueType,outputValueProperties...> outputValues,
             const Kokkos::DynRankView<inputPointValueType, inputPointProperties...>  inputPoints,
             const EOperator operatorType ) const {
#ifdef HAVE_INTREPID2_DEBUG
    Intrepid2::getValues_HGRAD_Args(outputValues,
                                    inputPoints,
                                    operatorType,
                                    obj_->getBaseCellTopology(),
                                    obj_->getCardinality() );
#endif

    typedef          Kokkos::DynRankView<outputValueValueType,outputValueProperties...>         outputValueViewType;
    typedef          Kokkos::DynRankView<inputPointValueType, inputPointProperties...>          inputPointViewType;
    typedef typename ExecSpace<typename inputPointViewType::execution_space,SpT>::ExecSpaceType ExecSpaceType;

    // loopSize corresponds to cardinality
    const auto loopSize = outputValues.dimension(0);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);

    const auto nbf = obj_->getCardinality();
    const auto npts = inputPoints.dimension(0);

    switch (operatorType) {
    case OPERATOR_VALUE: {
      Kokkos::DynRankView<OT,SpT> phisCur("Hgrad::Line::getValues::phisCur", nbf, npts);
      phis_.impl_.getValues( phisCur, inputPoints, operatorType );
      
      typedef Functor<outputValueViewType,inputPointViewType,OPERATOR_VALUE> FunctorType;
      Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints, phisCur) );
      break;
    }
    case OPERATOR_GRAD:
    case OPERATOR_D1:
    case OPERATOR_D2:
    case OPERATOR_D3:
    case OPERATOR_D4:
    case OPERATOR_D5:
    case OPERATOR_D6:
    case OPERATOR_D7:
    case OPERATOR_D8:
    case OPERATOR_D9:
    case OPERATOR_D10: {
      const auto dkcard = ( operatorType == OPERATOR_GRAD ? getDkCardinality(OPERATOR_D1, 1) :
                            /**/                            getDkCardinality(operatorType,1) );

      Kokkos::DynRankView<OT,SpT> phisCur("Hgrad::Line::getValues::phisCur", nbf, npts, dkcard);          
      phis_.impl_.getValues( phisCur , inputPoints , operatorType );

      typedef Functor<outputValueViewType,inputPointViewType,OPERATOR_MAX> FunctorType;
      Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints, phisCur) );
      break;
    }
    default: {
      INTREPID2_TEST_FOR_EXCEPTION( true , std::invalid_argument,
                                    ">>> ERROR (Basis_HGRAD_LINE_Cn_FEM): Operator type not implemented" );
      break;
    }
    }
  }
  

  template<typename SpT, typename OT, typename PT>
  template<typename dofCoordValueType, class ...dofCoordProperties>
  void
  Basis_HGRAD_LINE_Cn_FEM<SpT,OT,PT>::Internal::
  getDofCoords( Kokkos::DynRankView<dofCoordValueType,dofCoordProperties...> dofCoords ) const {
    
#ifdef HAVE_INTREPID2_DEBUG
    // Verify rank of output array.
    INTREPID2_TEST_FOR_EXCEPTION( dofCoords.rank() != 2, std::invalid_argument,
                                  ">>> ERROR: (Intrepid2::Basis_HGRAD_LINE_Cn_FEM::getDofCoords) rank = 2 required for dofCoords array");
    // Verify 0th dimension of output array.
    INTREPID2_TEST_FOR_EXCEPTION( dofCoords.dimension(0) != obj_->basisCardinality_, std::invalid_argument,
                                  ">>> ERROR: (Intrepid2::Basis_HGRAD_LINE_Cn_FEM::getDofCoords) mismatch in number of dof and 0th dimension of dofCoords array");
    // Verify 1st dimension of output array.
    INTREPID2_TEST_FOR_EXCEPTION( dofCoords.dimension(1) != obj_->basisCellTopology_.getDimension(), std::invalid_argument,
                                  ">>> ERROR: (Intrepid2::Basis_HGRAD_LINE_Cn_FEM::getDofCoords) incorrect reference cell (1st) dimension in dofCoords array");
#endif
    Kokkos::deep_copy(dofCoords, obj_->latticePts_);
  }
  
  
}// namespace Intrepid2
#endif















