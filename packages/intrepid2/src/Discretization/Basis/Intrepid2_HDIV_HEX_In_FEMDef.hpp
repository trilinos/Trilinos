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

/** \file   Intrepid2_HDIV_HEX_In_FEMDef.hpp
    \brief  Definition file for FEM basis functions of degree n for H(div) functions on HEX cells.
    \author Created by R. Kirby, P. Bochev, D. Ridzal and K. Peterson.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_HDIV_HEX_IN_FEM_DEF_HPP__
#define __INTREPID2_HDIV_HEX_IN_FEM_DEF_HPP__

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
    Basis_HDIV_HEX_In_FEM::Serial<opType>::
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
      const auto input_z = Kokkos::subview(input, Kokkos::ALL(), range_type(2,3));

      const ordinal_type dim_s = get_dimension_scalar(work);
      auto ptr0 = work.data();
      auto ptr1 = work.data()+cardLine*npts*dim_s;
      auto ptr2 = work.data()+2*cardLine*npts*dim_s;
      auto ptr3 = work.data()+(2*cardLine+cardBubble)*npts*dim_s;      
      
      typedef typename Kokkos::DynRankView<typename workViewType::value_type, typename workViewType::memory_space> viewType;
      auto vcprop = Kokkos::common_view_alloc_prop(work);

      switch (opType) {
      case OPERATOR_VALUE: {
        viewType workLine(Kokkos::view_wrap(ptr0, vcprop), cardLine, npts);
        viewType outputLine(Kokkos::view_wrap(ptr1, vcprop), cardLine, npts);
        viewType outputBubble_A(Kokkos::view_wrap(ptr2, vcprop), cardBubble, npts);
        viewType outputBubble_B(Kokkos::view_wrap(ptr3, vcprop), cardBubble, npts);
        
        // tensor product
        ordinal_type idx = 0;
        {
          Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
            getValues(outputLine, input_x, workLine, vinvLine);

          Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
            getValues(outputBubble_A, input_y, workLine, vinvBubble);

          Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
            getValues(outputBubble_B, input_z, workLine, vinvBubble);


          // x component (bubbleBasis(z) bubbleBasis(y) lineBasis(x))
          const auto output_x = outputLine;
          const auto output_y = outputBubble_A;
          const auto output_z = outputBubble_B;

          for (ordinal_type k=0;k<cardBubble;++k) // z
            for (ordinal_type j=0;j<cardBubble;++j) // y
              for (ordinal_type i=0;i<cardLine;++i,++idx) // x
                for (ordinal_type l=0;l<npts;++l) {
                  output.access(idx,l,0) = output_x.access(i,l)*output_y.access(j,l)*output_z.access(k,l);
                  output.access(idx,l,1) = 0.0;
                  output.access(idx,l,2) = 0.0;
                }
        }
        {
          Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
            getValues(outputBubble_A, input_x, workLine, vinvBubble);

          Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
            getValues(outputLine, input_y, workLine, vinvLine);

          //Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
          //  getValues(outputBubble_B, input_z, workLine, vinvBubble);

          // y component (lineBasis(z) bubbleBasis(y) lineBasis(x))
          const auto output_x = outputBubble_A;
          const auto output_y = outputLine;
          const auto output_z = outputBubble_B;

          for (ordinal_type k=0;k<cardBubble;++k) // z
            for (ordinal_type j=0;j<cardLine;++j) // y
              for (ordinal_type i=0;i<cardBubble;++i,++idx) // x
                for (ordinal_type l=0;l<npts;++l) {
                  output.access(idx,l,0) = 0.0;
                  output.access(idx,l,1) = output_x.access(i,l)*output_y.access(j,l)*output_z.access(k,l);
                  output.access(idx,l,2) = 0.0;
                }
        }
        {
          //Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
          //  getValues(outputBubble_A, input_x, workLine, vinvBubble);

          Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
            getValues(outputBubble_B, input_y, workLine, vinvBubble);

          Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
            getValues(outputLine, input_z, workLine, vinvLine);

          // z component (bubbleBasis(z) bubbleBasis(y) lineBasis(x))
          const auto output_x = outputBubble_A;
          const auto output_y = outputBubble_B;
          const auto output_z = outputLine;

          for (ordinal_type k=0;k<cardLine;++k) // z
            for (ordinal_type j=0;j<cardBubble;++j) // y
              for (ordinal_type i=0;i<cardBubble;++i,++idx) // x
                for (ordinal_type l=0;l<npts;++l) {
                  output.access(idx,l,0) = 0.0;
                  output.access(idx,l,1) = 0.0;
                  output.access(idx,l,2) = output_x.access(i,l)*output_y.access(j,l)*output_z.access(k,l);
                }
        }
        break;
      }
      case OPERATOR_DIV: {      
        viewType workLine(Kokkos::view_wrap(ptr0, vcprop), cardLine, npts);
        // A line value
        viewType outputBubble_A(Kokkos::view_wrap(ptr2, vcprop), cardBubble, npts);
        // B line value
        viewType outputBubble_B(Kokkos::view_wrap(ptr3, vcprop), cardBubble, npts);
        // Line grad
        viewType outputLine(Kokkos::view_wrap(ptr1, vcprop), cardLine, npts, 1);
        
        // tensor product
        ordinal_type idx = 0;

        { // x - component
          Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_Dn>::
            getValues(outputLine, input_x, workLine, vinvLine, 1);

          Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
            getValues(outputBubble_A, input_y, workLine, vinvBubble);

          Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
            getValues(outputBubble_B, input_z, workLine, vinvBubble);

          // x component (bubbleBasis(z) bubbleBasis(y) lineBasis(x))
          const auto output_dx = outputLine;
          const auto output_y  = outputBubble_A;
          const auto output_z  = outputBubble_B;

          for (ordinal_type k=0;k<cardBubble;++k) // z
            for (ordinal_type j=0;j<cardBubble;++j) // y
              for (ordinal_type i=0;i<cardLine;++i,++idx) // x
                for (ordinal_type l=0;l<npts;++l)
                  output.access(idx,l) =  output_dx.access(i,l,0)*output_y.access (j,l)  *output_z.access(k,l);
        }
        { // y - component
          Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
            getValues(outputBubble_A, input_x, workLine, vinvBubble);

          Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_Dn>::
            getValues(outputLine, input_y, workLine, vinvLine, 1);

          // Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
          //   getValues(outputBubble_B, input_z, workLine, vinvBubble);

          //(bubbleBasis(z) lineBasis(y) bubbleBasis(x))
          const auto output_x  = outputBubble_A;
          const auto output_dy = outputLine;
          const auto output_z  = outputBubble_B;

          for (ordinal_type k=0;k<cardBubble;++k) // z
            for (ordinal_type j=0;j<cardLine;++j) // y
              for (ordinal_type i=0;i<cardBubble;++i,++idx) // x
                for (ordinal_type l=0;l<npts;++l)
                  output.access(idx,l) = output_x.access(i,l)*output_dy.access(j,l,0)*output_z.access(k,l);
        }
        { // z - component
          // Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
          //   getValues(outputBubble_A, input_x, workLine, vinvBubble);

          Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
            getValues(outputBubble_B, input_y, workLine, vinvBubble);

          Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_Dn>::
            getValues(outputLine, input_z, workLine, vinvLine, 1);

          // (bubbleBasis(z) bubbleBasis(y) lineBasis(x))
          const auto output_x  = outputBubble_A;
          const auto output_y  = outputBubble_B;
          const auto output_dz = outputLine;

          for (ordinal_type k=0;k<cardLine;++k) // z
            for (ordinal_type j=0;j<cardBubble;++j) // y
              for (ordinal_type i=0;i<cardBubble;++i,++idx) // x
                for (ordinal_type l=0;l<npts;++l)
                  output.access(idx,l) = output_x.access(i,l)*output_y.access(j,l)*output_dz.access(k,l,0);
        }
        break;
      }
      default: {
        INTREPID2_TEST_FOR_ABORT( true,
                                  ">>> ERROR: (Intrepid2::Basis_HDIV_HEX_In_FEM::Serial::getValues) operator is not supported" );
      }
      }
    }

    template<typename SpT, ordinal_type numPtsPerEval,
             typename outputValueValueType, class ...outputValueProperties,
             typename inputPointValueType,  class ...inputPointProperties,
             typename vinvValueType,        class ...vinvProperties>
    void
    Basis_HDIV_HEX_In_FEM::
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
      ordinal_type order = 0;
      ordinal_type cardBubble;  // = std::cbrt(cardinality/3);
      ordinal_type cardLine;  // = cardBubble+1;      
      do {
        cardBubble = Intrepid2::getPnCardinality<1>(order);
        cardLine = Intrepid2::getPnCardinality<1>(++order);
      } while((3*cardBubble*cardBubble*cardLine !=  cardinality) && (order != Parameters::MaxOrder));

      auto vcprop = Kokkos::common_view_alloc_prop(inputPoints);
      typedef typename Kokkos::DynRankView< inputPointType, typename inputPointViewType::memory_space> workViewType;

      switch (operatorType) {
      case OPERATOR_VALUE: {
        auto workSize = Serial<OPERATOR_VALUE>::getWorkSizePerPoint(order);
        workViewType  work(Kokkos::view_alloc("Basis_HDIV_HEX_In_FEM::getValues::work", vcprop), workSize, inputPoints.extent(0));
        typedef Functor<outputValueViewType,inputPointViewType,vinvViewType, workViewType,
            OPERATOR_VALUE,numPtsPerEval> FunctorType;
        Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints, vinvLine, vinvBubble, work) );
        break;
      }
      case OPERATOR_DIV: {
        auto workSize = Serial<OPERATOR_DIV>::getWorkSizePerPoint(order);
        workViewType  work(Kokkos::view_alloc("Basis_HDIV_HEX_In_FEM::getValues::work", vcprop), workSize, inputPoints.extent(0));
        typedef Functor<outputValueViewType,inputPointViewType,vinvViewType, workViewType,
            OPERATOR_DIV,numPtsPerEval> FunctorType;
        Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints, vinvLine, vinvBubble, work) );
        break;
      }
      default: {
        INTREPID2_TEST_FOR_EXCEPTION( true , std::invalid_argument,
                                      ">>> ERROR (Basis_HDIV_HEX_In_FEM): Operator type not implemented" );
        // break;  commented out since exception is thrown
      }
      }
    }
  }

  // -------------------------------------------------------------------------------------

  template<typename SpT, typename OT, typename PT>
  Basis_HDIV_HEX_In_FEM<SpT,OT,PT>::
  Basis_HDIV_HEX_In_FEM( const ordinal_type order,
                         const EPointType   pointType ) {

    INTREPID2_TEST_FOR_EXCEPTION( !(pointType == POINTTYPE_EQUISPACED ||
                                    pointType == POINTTYPE_WARPBLEND), std::invalid_argument,
                                  ">>> ERROR (Basis_HDIV_HEX_In_FEM): pointType must be either equispaced or warpblend.");

    // this should be created in host and vinv should be deep copied into device space
    Basis_HGRAD_LINE_Cn_FEM<SpT,OT,PT> lineBasis( order, pointType );
    Basis_HGRAD_LINE_Cn_FEM<SpT,OT,PT> bubbleBasis( order - 1, POINTTYPE_GAUSS );

    const ordinal_type
      cardLine = lineBasis.getCardinality(),
      cardBubble = bubbleBasis.getCardinality();

    this->vinvLine_   = Kokkos::DynRankView<typename scalarViewType::value_type,SpT>("Hcurl::Hex::In::vinvLine", cardLine, cardLine);
    this->vinvBubble_ = Kokkos::DynRankView<typename scalarViewType::value_type,SpT>("Hcurl::Hex::In::vinvBubble", cardBubble, cardBubble);

    lineBasis.getVandermondeInverse(this->vinvLine_);
    bubbleBasis.getVandermondeInverse(this->vinvBubble_);

    this->basisCardinality_  = 3*cardLine*cardBubble*cardBubble;
    this->basisDegree_       = order;
    this->basisCellTopology_ = shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<8> >() );
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
      ordinal_type tags[3*maxCardLine*maxCardLine*maxCardLine][4];

      const ordinal_type face_yz[2] = {3, 1};
      const ordinal_type face_xz[2] = {0, 2};
      const ordinal_type face_xy[2] = {4, 5};

      {
        ordinal_type idx = 0;

        ///
        /// Product rule: z -> y -> x, x-normal first
        ///

        // since there are x/y components in the interior
        // dof sum should be computed before the information
        // is assigned to tags
        const ordinal_type
          intr_ndofs_per_direction = (cardLine-2)*cardBubble*cardBubble,
          intr_ndofs = 3*intr_ndofs_per_direction;

        // x component (bubbleBasis(z) bubbleBasis(y) lineBasis(x))
        for (ordinal_type k=0;k<cardBubble;++k) { // z
          const auto tag_z = bubbleBasis.getDofTag(k);
          for (ordinal_type j=0;j<cardBubble;++j) { // y
            const auto tag_y = bubbleBasis.getDofTag(j);
            for (ordinal_type i=0;i<cardLine;++i,++idx) { // x
              const auto tag_x = lineBasis.getDofTag(i);

              if (tag_x(0) == 0 && tag_y(0) == 1 && tag_z(0) == 1) {
                // face, x vert, y edge, z edge
                tags[idx][0] = 2; // face dof
                tags[idx][1] = face_yz[tag_x(1)]; // face id
                tags[idx][2] = tag_y(2) + tag_y(3)*tag_z(2); // local dof id
                tags[idx][3] = tag_y(3)*tag_z(3); // total number of dofs in this vertex
              } else {
                // interior
                tags[idx][0] = 3; // interior dof
                tags[idx][1] = 0;
                tags[idx][2] = tag_x(2) + tag_x(3)*tag_y(2) + tag_x(3)*tag_y(3)*tag_z(2);  // local dof id
                tags[idx][3] = intr_ndofs; // total number of dofs in this vertex
              }
            }
          }
        }

        // y component (bubbleBasis(z) lineBasis(y) bubbleBasis(x))
        for (ordinal_type k=0;k<cardBubble;++k) { // z
          const auto tag_z = bubbleBasis.getDofTag(k);
          for (ordinal_type j=0;j<cardLine;++j) { // y
            const auto tag_y = lineBasis.getDofTag(j);
            for (ordinal_type i=0;i<cardBubble;++i,++idx) { // x
              const auto tag_x = bubbleBasis.getDofTag(i);

              if (tag_x(0) == 1 && tag_y(0) == 0 && tag_z(0) == 1) {
                // face, x edge, y vert, z edge
                tags[idx][0] = 2; // face dof
                tags[idx][1] = face_xz[tag_y(1)]; // face id
                tags[idx][2] = tag_x(2) + tag_x(3)*tag_z(2); // local dof id
                tags[idx][3] = tag_x(3)*tag_z(3); // total number of dofs in this vertex
              } else {
                // interior
                tags[idx][0] = 3; // interior dof
                tags[idx][1] = 0;
                tags[idx][2] = intr_ndofs_per_direction + tag_x(2) + tag_x(3)*tag_y(2) + tag_x(3)*tag_y(3)*tag_z(2); // local dof id
                tags[idx][3] = intr_ndofs; // total number of dofs in this vertex
              }
            }
          }
        }

        // z component (bubbleBasis(z) lineBasis(y) lineBasis(x))
        for (ordinal_type k=0;k<cardLine;++k) { // y
          const auto tag_z = lineBasis.getDofTag(k);
          for (ordinal_type j=0;j<cardBubble;++j) { // z
            const auto tag_y = bubbleBasis.getDofTag(j);
            for (ordinal_type i=0;i<cardBubble;++i,++idx) { // x
              const auto tag_x = bubbleBasis.getDofTag(i);

              if (tag_x(0) == 1 && tag_y(0) == 1 && tag_z(0) == 0) {
                // face, x edge, y edge, z vert
                tags[idx][0] = 2; // face dof
                tags[idx][1] = face_xy[tag_z(1)]; // face id
                tags[idx][2] = tag_x(2) + tag_x(3)*tag_y(2); // local dof id
                tags[idx][3] = tag_x(3)*tag_y(3); // total number of dofs in this vertex
              } else {
                // interior
                tags[idx][0] = 3; // interior dof
                tags[idx][1] = 0;
                tags[idx][2] = 2*intr_ndofs_per_direction + tag_x(2) + tag_x(3)*tag_y(2) + tag_x(3)*tag_y(3)*tag_z(2); // local dof id
                tags[idx][3] = intr_ndofs; // total number of dofs in this vertex
              }
            }
          }
        }

        INTREPID2_TEST_FOR_EXCEPTION( idx != this->basisCardinality_ , std::runtime_error,
                                      ">>> ERROR (Basis_HDIV_HEX_In_FEM): " \
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

      // x component (lineBasis(z) lineBasis(y) bubbleBasis(x))
      for (ordinal_type k=0;k<cardBubble;++k) { // z
        for (ordinal_type j=0;j<cardBubble;++j) { // y
          for (ordinal_type i=0;i<cardLine;++i,++idx) { // x
            dofCoordsHost(idx,0) = dofCoordsLineHost(i,0);
            dofCoordsHost(idx,1) = dofCoordsBubbleHost(j,0);
            dofCoordsHost(idx,2) = dofCoordsBubbleHost(k,0);
            dofCoeffsHost(idx,0) = 1.0;
          }
        }
      }

      // y component (lineBasis(z) bubbleBasis(y) lineBasis(x))
      for (ordinal_type k=0;k<cardBubble;++k) { // z
        for (ordinal_type j=0;j<cardLine;++j) { // y
          for (ordinal_type i=0;i<cardBubble;++i,++idx) { // x
            dofCoordsHost(idx,0) = dofCoordsBubbleHost(i,0);
            dofCoordsHost(idx,1) = dofCoordsLineHost(j,0);
            dofCoordsHost(idx,2) = dofCoordsBubbleHost(k,0);
            dofCoeffsHost(idx,1) = 1.0;
          }
        }
      }

      // z component (bubbleBasis(z) lineBasis(y) lineBasis(x))
      for (ordinal_type k=0;k<cardLine;++k) { // z
        for (ordinal_type j=0;j<cardBubble;++j) { // y
          for (ordinal_type i=0;i<cardBubble;++i,++idx) { // x
            dofCoordsHost(idx,0) = dofCoordsBubbleHost(i,0);
            dofCoordsHost(idx,1) = dofCoordsBubbleHost(j,0);
            dofCoordsHost(idx,2) = dofCoordsLineHost(k,0);
            dofCoeffsHost(idx,2) = 1.0;
          }
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
