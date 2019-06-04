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

/** \file   Intrepid2_HDIV_QUAD_In_FEMDef.hpp
    \brief  Definition file for FEM basis functions of degree n for H(div) functions on QUAD cells.
    \author Created by R. Kirby, P. Bochev, D. Ridzal and K. Peterson.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_HDIV_QUAD_IN_FEM_DEF_HPP__
#define __INTREPID2_HDIV_QUAD_IN_FEM_DEF_HPP__

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
    Basis_HDIV_QUAD_In_FEM::Serial<opType>::
    getValues(       outputViewType output,
               const inputViewType  input,
                     workViewType   work,
               const vinvViewType   vinvLine,
               const vinvViewType   vinvBubble) {
      const ordinal_type cardLine = vinvLine.extent(0);
      const ordinal_type cardBubble = vinvBubble.extent(0);

      const ordinal_type npts = input.extent(0);

      typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;
      const auto input_x = Kokkos::subview(input, Kokkos::ALL(), range_type(0,1));
      const auto input_y = Kokkos::subview(input, Kokkos::ALL(), range_type(1,2));

      const int dim_s = get_dimension_scalar(work);
      auto ptr0 = work.data();
      auto ptr1 = work.data()+cardLine*npts*dim_s;
      auto ptr2 = work.data()+2*cardLine*npts*dim_s;
      
      
      typedef typename Kokkos::DynRankView<typename workViewType::value_type, typename workViewType::memory_space> viewType;
      auto vcprop = Kokkos::common_view_alloc_prop(work);

      switch (opType) {
      case OPERATOR_VALUE: {
        viewType workLine(Kokkos::view_wrap(ptr0, vcprop), cardLine, npts);
        viewType outputLine(Kokkos::view_wrap(ptr1, vcprop), cardLine, npts);
        viewType outputBubble(Kokkos::view_wrap(ptr2, vcprop), cardBubble, npts);

        // tensor product
        ordinal_type idx = 0;
        {
          Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
            getValues(outputBubble, input_x, workLine, vinvBubble);

          Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
            getValues(outputLine, input_y, workLine, vinvLine);

          // x component (lineBasis(y) bubbleBasis(x))
          const auto output_x = outputBubble;
          const auto output_y = outputLine;

          for (ordinal_type j=0;j<cardLine;++j) // y
            for (ordinal_type i=0;i<cardBubble;++i,++idx) // x
              for (ordinal_type k=0;k<npts;++k) {
                output.access(idx,k,0) = 0.0;
                output.access(idx,k,1) = output_x.access(i,k)*output_y.access(j,k);
              }
        }
        {
          Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
            getValues(outputBubble, input_y, workLine, vinvBubble);

          Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
            getValues(outputLine, input_x, workLine, vinvLine);

          // y component (bubbleBasis(y) lineBasis(x))
          const auto output_x = outputLine;
          const auto output_y = outputBubble;
          for (ordinal_type j=0;j<cardBubble;++j) // y
            for (ordinal_type i=0;i<cardLine;++i,++idx) // x
              for (ordinal_type k=0;k<npts;++k) {
                output.access(idx,k,0) = output_x.access(i,k)*output_y.access(j,k);
                output.access(idx,k,1) = 0.0;
              }
        }
        break;
      }
      case OPERATOR_DIV: {
        ordinal_type idx = 0;
        { // x - component
          viewType workLine(Kokkos::view_wrap(ptr0, vcprop), cardLine, npts);
          // x bubble value
          viewType output_x(Kokkos::view_wrap(ptr2, vcprop), cardBubble, npts);
          // y line grad
          viewType output_y(Kokkos::view_wrap(ptr1, vcprop), cardLine, npts,1);

          Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
            getValues(output_x, input_x, workLine, vinvBubble);

          Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_Dn>::
            getValues(output_y, input_y, workLine, vinvLine, 1);

          // tensor product (extra dimension of ouput x and y are ignored)
          for (ordinal_type j=0;j<cardLine;++j) // y
            for (ordinal_type i=0;i<cardBubble;++i,++idx) // x
              for (ordinal_type k=0;k<npts;++k)
                output.access(idx,k) = output_x.access(i,k)*output_y.access(j,k,0);
        }
        { // y - component
          viewType workLine(Kokkos::view_wrap(ptr0, vcprop), cardLine, npts);
          // x line grad
          viewType output_x(Kokkos::view_wrap(ptr1, vcprop), cardLine, npts,1);
          // y bubble value
          viewType output_y(Kokkos::view_wrap(ptr2, vcprop), cardBubble, npts);

          Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
            getValues(output_y, input_y, workLine, vinvBubble);

          Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_Dn>::
            getValues(output_x, input_x, workLine, vinvLine, 1);

          // tensor product (extra dimension of ouput x and y are ignored)
          for (ordinal_type j=0;j<cardBubble;++j) // y
            for (ordinal_type i=0;i<cardLine;++i,++idx) // x
              for (ordinal_type k=0;k<npts;++k)
                output.access(idx,k) = output_x.access(i,k,0)*output_y.access(j,k);
        }
        break;
      }
      default: {
        INTREPID2_TEST_FOR_ABORT( true,
                                  ">>> ERROR: (Intrepid2::Basis_HDIV_QUAD_In_FEM::Serial::getValues) operator is not supported" );
      }
      }
    }

    template<typename SpT, ordinal_type numPtsPerEval,
             typename outputValueValueType, class ...outputValueProperties,
             typename inputPointValueType,  class ...inputPointProperties,
             typename vinvValueType,        class ...vinvProperties>
    void
    Basis_HDIV_QUAD_In_FEM::
    getValues(       Kokkos::DynRankView<outputValueValueType,outputValueProperties...> outputValues,
               const Kokkos::DynRankView<inputPointValueType, inputPointProperties...>  inputPoints,
               const Kokkos::DynRankView<vinvValueType,       vinvProperties...>        vinvLine,
               const Kokkos::DynRankView<vinvValueType,       vinvProperties...>        vinvBubble,
               const EOperator operatorType ) {
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
      //get basis order based on basis cardinality.
      ordinal_type order = 0;   // = std::sqrt(cardinality/2);
      ordinal_type cardBubble;  // = std::sqrt(cardinality/2);
      ordinal_type cardLine;  // = cardBubble+1;      
      do {
        cardBubble = Intrepid2::getPnCardinality<1>(order);
        cardLine = Intrepid2::getPnCardinality<1>(++order);
      } while((2*cardBubble*cardLine !=  cardinality) && (order != Parameters::MaxOrder));

      auto vcprop = Kokkos::common_view_alloc_prop(inputPoints);
      typedef typename Kokkos::DynRankView< inputPointType, typename inputPointViewType::memory_space> workViewType;

      switch (operatorType) {
      case OPERATOR_VALUE: {
        auto workSize = Serial<OPERATOR_VALUE>::getWorkSizePerPoint(order);
        workViewType  work(Kokkos::view_alloc("Basis_HDIV_QUAD_In_FEM::getValues::work", vcprop), workSize, inputPoints.extent(0));
        typedef Functor<outputValueViewType,inputPointViewType,vinvViewType, workViewType,
            OPERATOR_VALUE,numPtsPerEval> FunctorType;
        Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints, vinvLine, vinvBubble, work) );
        break;
      }
      case OPERATOR_DIV: {
        auto workSize = Serial<OPERATOR_DIV>::getWorkSizePerPoint(order);
        workViewType  work(Kokkos::view_alloc("Basis_HDIV_QUAD_In_FEM::getValues::work", vcprop), workSize, inputPoints.extent(0));
        typedef Functor<outputValueViewType,inputPointViewType,vinvViewType, workViewType,
            OPERATOR_DIV,numPtsPerEval> FunctorType;
        Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints, vinvLine, vinvBubble, work) );
        break;
      }
      default: {
        INTREPID2_TEST_FOR_EXCEPTION( true , std::invalid_argument,
                                      ">>> ERROR (Basis_HDIV_QUAD_In_FEM): Operator type not implemented" );
      }
      }
    }
  }

  // -------------------------------------------------------------------------------------
  template<typename SpT, typename OT, typename PT>
  Basis_HDIV_QUAD_In_FEM<SpT,OT,PT>::
  Basis_HDIV_QUAD_In_FEM( const ordinal_type order,
                           const EPointType   pointType ) {

    INTREPID2_TEST_FOR_EXCEPTION( !(pointType == POINTTYPE_EQUISPACED ||
                                    pointType == POINTTYPE_WARPBLEND), std::invalid_argument,
                                  ">>> ERROR (Basis_HDIV_QUAD_In_FEM): pointType must be either equispaced or warpblend.");

    // this should be in host
    Basis_HGRAD_LINE_Cn_FEM<SpT,OT,PT> lineBasis( order, pointType );
    Basis_HGRAD_LINE_Cn_FEM<SpT,OT,PT> bubbleBasis( order - 1, POINTTYPE_GAUSS );

    const ordinal_type
      cardLine = lineBasis.getCardinality(),
      cardBubble = bubbleBasis.getCardinality();

    this->vinvLine_   = Kokkos::DynRankView<typename scalarViewType::value_type,SpT>("Hdiv::Quad::In::vinvLine", cardLine, cardLine);
    this->vinvBubble_ = Kokkos::DynRankView<typename scalarViewType::value_type,SpT>("Hdiv::Quad::In::vinvBubble", cardBubble, cardBubble);

    lineBasis.getVandermondeInverse(this->vinvLine_);
    bubbleBasis.getVandermondeInverse(this->vinvBubble_);

    this->basisCardinality_  = 2*cardLine*cardBubble;
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
      constexpr ordinal_type maxCardBubble = Parameters::MaxOrder;
      ordinal_type tags[2*maxCardLine*maxCardBubble][4];

      const ordinal_type edge_x[2] = {0,2};
      const ordinal_type edge_y[2] = {3,1};
      {
        ordinal_type idx = 0;

        /// 
        /// Product rule: y -> x, x-normal first
        ///

        // since there are x/y components in the interior
        // dof sum should be computed before the information
        // is assigned to tags
        const ordinal_type
          intr_ndofs_per_direction = (cardLine-2)*cardBubble,
          intr_ndofs = 2*intr_ndofs_per_direction;

        // x component (lineBasis(y) bubbleBasis(x))
        for (ordinal_type j=0;j<cardLine;++j) { // y
          const auto tag_y = lineBasis.getDofTag(j);
          for (ordinal_type i=0;i<cardBubble;++i,++idx) { // x
            const auto tag_x = bubbleBasis.getDofTag(i);

            if (tag_x(0) == 1 && tag_y(0) == 0) {
              // edge: x edge, y vert
              tags[idx][0] = 1; // edge dof
              tags[idx][1] = edge_x[tag_y(1)];
              tags[idx][2] = tag_x(2); // local dof id
              tags[idx][3] = tag_x(3); // total number of dofs in this vertex
            } else {
              // interior
              tags[idx][0] = 2; // interior dof
              tags[idx][1] = 0;
              tags[idx][2] = tag_x(2) + tag_x(3)*tag_y(2); // local dof id
              tags[idx][3] = intr_ndofs; // total number of dofs in this vertex
            }
          }
        }

        // y component (bubbleBasis(y) lineBasis(x))
        for (ordinal_type j=0;j<cardBubble;++j) { // y
          const auto tag_y = bubbleBasis.getDofTag(j);
          for (ordinal_type i=0;i<cardLine;++i,++idx) { // x
            const auto tag_x = lineBasis.getDofTag(i);

            if (tag_x(0) == 0 && tag_y(0) == 1) {
              // edge: x vert, y edge
              tags[idx][0] = 1; // edge dof
              tags[idx][1] = edge_y[tag_x(1)];
              tags[idx][2] = tag_y(2); // local dof id
              tags[idx][3] = tag_y(3); // total number of dofs in this vertex
            } else {
              // interior
              tags[idx][0] = 2; // interior dof
              tags[idx][1] = 0;
              tags[idx][2] = intr_ndofs_per_direction + tag_x(2) + tag_x(3)*tag_y(2); // local dof id
              tags[idx][3] = intr_ndofs; // total number of dofs in this vertex
            }
          }
        }
        INTREPID2_TEST_FOR_EXCEPTION( idx != this->basisCardinality_ , std::runtime_error,
                                      ">>> ERROR (Basis_HDIV_QUAD_In_FEM): " \
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

    // dofCoeffs on host and create its mirror view to device
    Kokkos::DynRankView<typename scalarViewType::value_type,typename SpT::array_layout,Kokkos::HostSpace>
      dofCoeffsHost("dofCoeffsHost", this->basisCardinality_, this->basisCellTopology_.getDimension());

    Kokkos::DynRankView<typename scalarViewType::value_type,SpT>
      dofCoordsLine("dofCoordsLine", cardLine, 1),
      dofCoordsBubble("dofCoordsBubble", cardBubble, 1);

    lineBasis.getDofCoords(dofCoordsLine);
    auto dofCoordsLineHost = Kokkos::create_mirror_view(Kokkos::HostSpace(), dofCoordsLine);
    Kokkos::deep_copy(dofCoordsLineHost, dofCoordsLine);

    bubbleBasis.getDofCoords(dofCoordsBubble);
    auto dofCoordsBubbleHost = Kokkos::create_mirror_view(Kokkos::HostSpace(), dofCoordsBubble);
    Kokkos::deep_copy(dofCoordsBubbleHost, dofCoordsBubble);

    {
      ordinal_type idx = 0;

      // x component (lineBasis(y) bubbleBasis(x))
      for (ordinal_type j=0;j<cardLine;++j) { // y
        for (ordinal_type i=0;i<cardBubble;++i,++idx) { // x
          dofCoordsHost(idx,0) = dofCoordsBubbleHost(i,0);
          dofCoordsHost(idx,1) = dofCoordsLineHost(j,0);
          dofCoeffsHost(idx,1) = 1.0;
        }
      }

      // y component (bubbleBasis(y) lineBasis(x))
      for (ordinal_type j=0;j<cardBubble;++j) { // y
        for (ordinal_type i=0;i<cardLine;++i,++idx) { // x
          dofCoordsHost(idx,0) = dofCoordsLineHost(i,0);
          dofCoordsHost(idx,1) = dofCoordsBubbleHost(j,0);
          dofCoeffsHost(idx,0) = 1.0;
        }
      }
    }

    this->dofCoords_ = Kokkos::create_mirror_view(typename SpT::memory_space(), dofCoordsHost);
    Kokkos::deep_copy(this->dofCoords_, dofCoordsHost);

    this->dofCoeffs_ = Kokkos::create_mirror_view(typename SpT::memory_space(), dofCoeffsHost);
    Kokkos::deep_copy(this->dofCoeffs_, dofCoeffsHost);
  }

}

#endif
