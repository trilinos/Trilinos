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

/** \file   Intrepid_HDIV_QUAD_Bn_FEMDef.hpp
    \brief  Definition file for FEM basis functions of degree n for H(div) functions on QUAD cells.
    \author Created by R. Kirby, P. Bochev, D. Ridzal and K. Peterson.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_HDIV_QUAD_BN_FEM_DEF_HPP__
#define __INTREPID2_HDIV_QUAD_BN_FEM_DEF_HPP__

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
    Basis_HDIV_QUAD_Bn_FEM::Serial<opType>::
    getValues( /**/  outputViewType output,
               const inputViewType  input,
               /**/  workViewType   work,
               const vinvViewType   vinvBubble,
               const ordinal_type   ort ) {
      const ordinal_type 
        cardBubble = vinvBubble.dimension(0),
        npts = input.dimension(0);

      typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;
      const auto input_x = Kokkos::subview(input, Kokkos::ALL(), range_type(0,1));
      const auto input_y = Kokkos::subview(input, Kokkos::ALL(), range_type(1,2));

      switch (opType) {
      case OPERATOR_VALUE: {
        auto ptr = work.data();
        
        Kokkos::DynRankView<typename workViewType::value_type,
            typename workViewType::memory_space> workBubble(ptr, cardBubble, npts);
        ptr += (cardBubble*npts);
        
        Kokkos::DynRankView<typename workViewType::value_type,
            typename workViewType::memory_space> output_x(ptr, cardBubble, npts);
        ptr += (cardBubble*npts);

        Kokkos::DynRankView<typename workViewType::value_type,
            typename workViewType::memory_space> output_y(ptr, cardBubble, npts);
        ptr += (cardBubble*npts);

        // tensor product
        ordinal_type idx = 0;
        Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
          getValues(output_x, input_x, workBubble, vinvBubble);
        
        Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
          getValues(output_y, input_y, workBubble, vinvBubble);

        const double signVal = (ort > 3 ? -1.0 : 1.0);
        
        // normal component (bubbleBasis(y) bubbleBasis(x))
        for (ordinal_type j=0;j<cardBubble;++j) // y
          for (ordinal_type i=0;i<cardBubble;++i,++idx) // x
            for (ordinal_type k=0;k<npts;++k) 
              output(idx,k) = signVal*output_x(i,k)*output_y(j,k);
        break;
      }
      default: {
        INTREPID2_TEST_FOR_ABORT( true,
                                  ">>> ERROR: (Intrepid2::Basis_HDIV_QUAD_Bn_FEM::Serial::getValues) operator is not supported" );
      }
      }
    }

    template<typename SpT, ordinal_type numPtsPerEval,
             typename outputValueValueType, class ...outputValueProperties,
             typename inputPointValueType,  class ...inputPointProperties,
             typename vinvValueType,        class ...vinvProperties>
    void
    Basis_HDIV_QUAD_Bn_FEM::
    getValues( /**/  Kokkos::DynRankView<outputValueValueType,outputValueProperties...> outputValues,
               const Kokkos::DynRankView<inputPointValueType, inputPointProperties...>  inputPoints,
               const Kokkos::DynRankView<vinvValueType,       vinvProperties...>        vinvBubble,
               const EOperator operatorType,
               const ordinal_type ort ) {
      typedef          Kokkos::DynRankView<outputValueValueType,outputValueProperties...>         outputValueViewType;
      typedef          Kokkos::DynRankView<inputPointValueType, inputPointProperties...>          inputPointViewType;
      typedef          Kokkos::DynRankView<vinvValueType,       vinvProperties...>                vinvViewType;
      typedef typename ExecSpace<typename inputPointViewType::execution_space,SpT>::ExecSpaceType ExecSpaceType;

      // loopSize corresponds to cardinality
      const auto loopSizeTmp1 = (inputPoints.dimension(0)/numPtsPerEval);
      const auto loopSizeTmp2 = (inputPoints.dimension(0)%numPtsPerEval != 0);
      const auto loopSize = loopSizeTmp1 + loopSizeTmp2;
      Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
      
      switch (operatorType) {
      case OPERATOR_VALUE: {
        typedef Functor<outputValueViewType,inputPointViewType,vinvViewType,
            OPERATOR_VALUE,numPtsPerEval> FunctorType;
        Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints, vinvBubble, ort ) );
        break;
      }
      default: {
        INTREPID2_TEST_FOR_EXCEPTION( true , std::invalid_argument,
                                      ">>> ERROR (Basis_HDIV_QUAD_Bn_FEM): Operator type not implemented" );
        break;
      }
      }
    }
  }

  // -------------------------------------------------------------------------------------
  template<typename SpT, typename OT, typename PT>
  Basis_HDIV_QUAD_Bn_FEM<SpT,OT,PT>::
  Basis_HDIV_QUAD_Bn_FEM( const ordinal_type order,
                          const EPointType   pointType ) {

    INTREPID2_TEST_FOR_EXCEPTION( !(pointType == POINTTYPE_EQUISPACED ||
                                    pointType == POINTTYPE_WARPBLEND), std::invalid_argument,
                                  ">>> ERROR (Basis_HDIV_QUAD_In_FEM): pointType must be either equispaced or warpblend.");

    // this should be in host
    Basis_HGRAD_LINE_Cn_FEM<SpT,OT,PT> bubbleBasis( order - 1, POINTTYPE_GAUSS );

    const ordinal_type cardBubble = bubbleBasis.getCardinality();
    this->vinvBubble_ = Kokkos::DynRankView<OT,SpT>("Hdiv::Quad::Bn::vinvBubble", cardBubble, cardBubble);

    bubbleBasis.getVandermondeInverse(this->vinvBubble_);

    this->basisCardinality_  = cardBubble*cardBubble;
    this->basisDegree_       = order;
    this->basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<4> >() );
    this->basisType_         = BASIS_FEM_FIAT;
    this->basisCoordinates_  = COORDINATES_CARTESIAN;

    // initialize tags
    {
      // Basis-dependent initializations
      const ordinal_type tagSize  = 4;        // size of DoF tag, i.e., number of fields in the tag
      const ordinal_type posScDim = 0;        // position in the tag, counting from 0, of the subcell dim
      const ordinal_type posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
      const ordinal_type posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell

      // An array with local DoF tags assigned to the basis functions, in the order of their local enumeration
      constexpr ordinal_type maxCardLine = Parameters::MaxOrder + 1;
      ordinal_type tags[maxCardLine*maxCardLine][4];

      {
        ordinal_type idx = 0;
        ///
        /// Product rule: y -> x, x-normal first
        ///

        // since there are x/y components in the interior
        // dof sum should be computed before the information
        // is assigned to tags
        const ordinal_type intr_ndofs = cardBubble*cardBubble;

        // normal component (bubbleBasis(y) bubbleBasis(x))
        for (ordinal_type j=0;j<cardBubble;++j) { // y
          const auto tag_y = bubbleBasis.getDofTag(j);
          for (ordinal_type i=0;i<cardBubble;++i,++idx) { // x
            const auto tag_x = bubbleBasis.getDofTag(i);
            
            // interior
            tags[idx][0] = 2; // interior dof
            tags[idx][1] = 0;
            tags[idx][2] = tag_x(2) + tag_x(3)*tag_y(2); // local dof id
            tags[idx][3] = intr_ndofs; // total number of dofs interior nodes
          }
        }
        
        INTREPID2_TEST_FOR_EXCEPTION( idx != this->basisCardinality_ , std::runtime_error,
                                      ">>> ERROR (Basis_HDIV_QUAD_Bn_FEM): " \
                                      "counted tag index is not same as cardinality." );
      }
      
      ordinal_type_array_1d_host tagView(&tags[0][0], this->basisCardinality_*4);

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
    
    // dofCoords on host and create its mirror view to device
    Kokkos::DynRankView<typename scalarViewType::value_type,typename SpT::array_layout,Kokkos::HostSpace>
      dofCoordsHost("dofCoordsHost", this->basisCardinality_, this->basisCellTopology_.getDimension());

    Kokkos::DynRankView<typename scalarViewType::value_type,SpT>
      dofCoordsBubble("dofCoordsBubble", cardBubble, 1);

    bubbleBasis.getDofCoords(dofCoordsBubble);
    auto dofCoordsBubbleHost = Kokkos::create_mirror_view(Kokkos::HostSpace(), dofCoordsBubble);
    Kokkos::deep_copy(dofCoordsBubbleHost, dofCoordsBubble);
    
    {
      ordinal_type idx = 0;

      // normal component (bubbleBasis(y) bubbleBasis(x))
      for (ordinal_type j=0;j<cardBubble;++j) { // y
        for (ordinal_type i=0;i<cardBubble;++i,++idx) { // x
          dofCoordsHost(idx,0) = dofCoordsBubbleHost(i,0);
          dofCoordsHost(idx,1) = dofCoordsBubbleHost(j,0);
        }
      }
    }
    
    this->dofCoords_ = Kokkos::create_mirror_view(typename SpT::memory_space(), dofCoordsHost);
    Kokkos::deep_copy(this->dofCoords_, dofCoordsHost);
  }
  
}

#endif
