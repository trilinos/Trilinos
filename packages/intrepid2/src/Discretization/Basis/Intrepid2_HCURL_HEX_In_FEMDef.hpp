// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_HCURL_HEX_In_FEMDef.hpp
    \brief  Definition file for FEM basis functions of degree n for H(curl) functions on HEX cells.
    \author Created by R. Kirby, P. Bochev, D. Ridzal and K. Peterson.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_HCURL_HEX_IN_FEM_DEF_HPP__
#define __INTREPID2_HCURL_HEX_IN_FEM_DEF_HPP__

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
    Basis_HCURL_HEX_In_FEM::Serial<opType>::
    getValues(       OutputViewType output,
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
      auto ptr1 = work.data() + cardLine*npts*dim_s;
      auto ptr2 = work.data() + 2*cardLine*npts*dim_s;
      auto ptr3 = work.data() + 3*cardLine*npts*dim_s;   
      
      typedef typename Kokkos::DynRankView<typename workViewType::value_type, typename workViewType::memory_space> viewType;
      auto vcprop = Kokkos::common_view_alloc_prop(work);

      switch (opType) {
      case OPERATOR_VALUE: {

        viewType workLine(Kokkos::view_wrap(ptr0, vcprop), cardLine, npts);
        viewType outputLine_A(Kokkos::view_wrap(ptr1, vcprop), cardLine, npts);
        viewType outputLine_B(Kokkos::view_wrap(ptr2, vcprop), cardLine, npts);
        viewType outputBubble(Kokkos::view_wrap(ptr3, vcprop), cardBubble, npts);
        
        // tensor product
        ordinal_type idx = 0;
        {
          Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
            getValues(outputBubble, input_x, workLine, vinvBubble);

          Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
            getValues(outputLine_A, input_y, workLine, vinvLine);

          Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
            getValues(outputLine_B, input_z, workLine, vinvLine);

          // x component (lineBasis(z) lineBasis(y) bubbleBasis(x))
          const auto output_x = outputBubble;
          const auto output_y = outputLine_A;
          const auto output_z = outputLine_B;

          for (ordinal_type k=0;k<cardLine;++k) // z
            for (ordinal_type j=0;j<cardLine;++j) // y
              for (ordinal_type i=0;i<cardBubble;++i,++idx) // x
                for (ordinal_type l=0;l<npts;++l) {
                  output.access(idx,l,0) = output_x.access(i,l)*output_y.access(j,l)*output_z.access(k,l);
                  output.access(idx,l,1) = 0.0;
                  output.access(idx,l,2) = 0.0;
                }
        }
        {
          Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
            getValues(outputLine_A, input_x, workLine, vinvLine);

          Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
            getValues(outputBubble, input_y, workLine, vinvBubble);

          //Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
          //  getValues(outputLine_B, input_z, workLine, vinvLine);

          // y component (lineBasis(z) bubbleBasis(y) lineBasis(x))
          const auto output_x = outputLine_A;
          const auto output_y = outputBubble;
          const auto output_z = outputLine_B;

          for (ordinal_type k=0;k<cardLine;++k) // z
            for (ordinal_type j=0;j<cardBubble;++j) // y
              for (ordinal_type i=0;i<cardLine;++i,++idx) // x
                for (ordinal_type l=0;l<npts;++l) {
                  output.access(idx,l,0) = 0.0;
                  output.access(idx,l,1) = output_x.access(i,l)*output_y.access(j,l)*output_z.access(k,l);
                  output.access(idx,l,2) = 0.0;
                }
        }
        {
          //Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
          //  getValues(outputLine_A, input_x, workLine, vinvLine);

          Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
            getValues(outputLine_B, input_y, workLine, vinvLine);

          Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
            getValues(outputBubble, input_z, workLine, vinvBubble);

          // z component (bubbleBasis(z) bubbleBasis(y) lineBasis(x))
          const auto output_x = outputLine_A;
          const auto output_y = outputLine_B;
          const auto output_z = outputBubble;

          for (ordinal_type k=0;k<cardBubble;++k) // z
            for (ordinal_type j=0;j<cardLine;++j) // y
              for (ordinal_type i=0;i<cardLine;++i,++idx) // x
                for (ordinal_type l=0;l<npts;++l) {
                  output.access(idx,l,0) = 0.0;
                  output.access(idx,l,1) = 0.0;
                  output.access(idx,l,2) = output_x.access(i,l)*output_y.access(j,l)*output_z.access(k,l);
                }
        }
        break;
      }
      case OPERATOR_CURL: {
          
         auto ptr4 = work.data() + 4*cardLine*npts*dim_s;
         auto ptr5 = work.data() + 5*cardLine*npts*dim_s;
         
        viewType workLine(Kokkos::view_wrap(ptr0, vcprop), cardLine, npts);
        viewType outputLine_A(Kokkos::view_wrap(ptr1, vcprop), cardLine, npts);
        viewType outputLine_B(Kokkos::view_wrap(ptr2, vcprop), cardLine, npts);
        viewType outputLine_DA(Kokkos::view_wrap(ptr3, vcprop), cardLine, npts, 1);
        viewType outputLine_DB(Kokkos::view_wrap(ptr4, vcprop), cardLine, npts, 1);
        viewType outputBubble(Kokkos::view_wrap(ptr5, vcprop), cardBubble, npts);

        // tensor product
        ordinal_type idx = 0;

        { // x - component
          Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
            getValues(outputBubble, input_x, workLine, vinvBubble);

          Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
            getValues(outputLine_A, input_y, workLine, vinvLine);

          Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_Dn>::
            getValues(outputLine_DA, input_y, workLine, vinvLine, 1);

          Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
            getValues(outputLine_B, input_z, workLine, vinvLine);

          Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_Dn>::
            getValues(outputLine_DB, input_z, workLine, vinvLine, 1);

          // x component (lineBasis(z) lineBasis(y) bubbleBasis(x))
          const auto output_x  = outputBubble;
          const auto output_y  = outputLine_A;
          const auto output_dy = outputLine_DA;
          const auto output_z  = outputLine_B;
          const auto output_dz = outputLine_DB;

          for (ordinal_type k=0;k<cardLine;++k) // z
            for (ordinal_type j=0;j<cardLine;++j) // y
              for (ordinal_type i=0;i<cardBubble;++i,++idx) // x
                for (ordinal_type l=0;l<npts;++l) {
                  output.access(idx,l,0) =  0.0;
                  output.access(idx,l,1) =  output_x.access(i,l)*output_y.access (j,l)  *output_dz.access(k,l,0);
                  output.access(idx,l,2) = -output_x.access(i,l)*output_dy.access(j,l,0)*output_z.access (k,l);
                }
        }
        { // y - component
          Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
            getValues(outputLine_A, input_x, workLine, vinvLine);

          Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_Dn>::
            getValues(outputLine_DA, input_x, workLine, vinvLine, 1);

          Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
            getValues(outputBubble, input_y, workLine, vinvBubble);

          //Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
          //  getValues(outputLine_B, input_z, workLine, vinvLine);

          //Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_Dn>::
          //  getValues(outputLine_DB, input_z, workLine, vinvLine, 1);

          // x component (lineBasis(z) lineBasis(y) bubbleBasis(x))
          const auto output_x  = outputLine_A;
          const auto output_dx = outputLine_DA;
          const auto output_y  = outputBubble;
          const auto output_z  = outputLine_B;
          const auto output_dz = outputLine_DB;

          for (ordinal_type k=0;k<cardLine;++k) // z
            for (ordinal_type j=0;j<cardBubble;++j) // y
              for (ordinal_type i=0;i<cardLine;++i,++idx) // x
                for (ordinal_type l=0;l<npts;++l) {
                  output.access(idx,l,0) = -output_x.access (i,l)  *output_y.access(j,l)*output_dz.access(k,l,0);
                  output.access(idx,l,1) =  0.0;
                  output.access(idx,l,2) =  output_dx(i,l,0)*output_y.access(j,l)*output_z.access (k,l);
                }
        }
        { // z - component
          // Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
          //   getValues(outputLine_A, input_x, workLine, vinvLine);

          // Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_Dn>::
          //   getValues(outputLine_DA, input_x, workLine, vinvLine, 1);

          Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
            getValues(outputLine_B, input_y, workLine, vinvLine);

          Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_Dn>::
            getValues(outputLine_DB, input_y, workLine, vinvLine, 1);

          Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
            getValues(outputBubble, input_z, workLine, vinvBubble);

          // x component (lineBasis(z) lineBasis(y) bubbleBasis(x))
          const auto output_x  = outputLine_A;
          const auto output_dx = outputLine_DA;
          const auto output_y  = outputLine_B;
          const auto output_dy = outputLine_DB;
          const auto output_z  = outputBubble;

          for (ordinal_type k=0;k<cardBubble;++k) // z
            for (ordinal_type j=0;j<cardLine;++j) // y
              for (ordinal_type i=0;i<cardLine;++i,++idx) // x
                for (ordinal_type l=0;l<npts;++l) {
                  output.access(idx,l,0) =  output_x.access (i,l)  *output_dy.access(j,l,0)*output_z.access(k,l);
                  output.access(idx,l,1) = -output_dx(i,l,0)*output_y.access (j,l)  *output_z.access(k,l);
                  output.access(idx,l,2) =  0.0;
                }
        }
        break;
      }
      default: {
        INTREPID2_TEST_FOR_ABORT( true,
                                  ">>> ERROR: (Intrepid2::Basis_HCURL_HEX_In_FEM::Serial::getValues) operator is not supported" );
      }
      }
    }
    template<typename DT, ordinal_type numPtsPerEval,
             typename outputValueValueType, class ...outputValueProperties,
             typename inputPointValueType,  class ...inputPointProperties,
             typename vinvValueType,        class ...vinvProperties>
    void
    Basis_HCURL_HEX_In_FEM::
    getValues(       Kokkos::DynRankView<outputValueValueType,outputValueProperties...> outputValues,
               const Kokkos::DynRankView<inputPointValueType, inputPointProperties...>  inputPoints,
               const Kokkos::DynRankView<vinvValueType,       vinvProperties...>        vinvLine,
               const Kokkos::DynRankView<vinvValueType,       vinvProperties...>        vinvBubble,
               const EOperator operatorType ) {
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
      //get basis order based on basis cardinality.
      ordinal_type order = 0;
      ordinal_type cardBubble;  // = std::cbrt(cardinality/3);
      ordinal_type cardLine;  // = cardBubble+1;      
      do {
        cardBubble = Intrepid2::getPnCardinality<1>(order);
        cardLine = Intrepid2::getPnCardinality<1>(++order);
      } while((3*cardBubble*cardLine*cardLine !=  cardinality) && (order != Parameters::MaxOrder));
      auto vcprop = Kokkos::common_view_alloc_prop(inputPoints);
      typedef typename Kokkos::DynRankView< inputPointType, typename inputPointViewType::memory_space> workViewType;

      switch (operatorType) {
      case OPERATOR_VALUE: {
        auto workSize = Serial<OPERATOR_VALUE>::getWorkSizePerPoint(order);
        workViewType  work(Kokkos::view_alloc("Basis_CURL_HEX_In_FEM::getValues::work", vcprop), workSize, inputPoints.extent(0));
        typedef Functor<outputValueViewType,inputPointViewType,vinvViewType, workViewType,
            OPERATOR_VALUE,numPtsPerEval> FunctorType;
        Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints, vinvLine, vinvBubble, work) );
        break;
      }
      case OPERATOR_CURL: {
        auto workSize = Serial<OPERATOR_CURL>::getWorkSizePerPoint(order);
        workViewType  work(Kokkos::view_alloc("Basis_CURL_HEX_In_FEM::getValues::work", vcprop), workSize, inputPoints.extent(0));
        typedef Functor<outputValueViewType,inputPointViewType,vinvViewType, workViewType,
            OPERATOR_CURL,numPtsPerEval> FunctorType;
        Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints, vinvLine, vinvBubble, work) );
        break;
      }
      default: {
        INTREPID2_TEST_FOR_EXCEPTION( true , std::invalid_argument,
                                      ">>> ERROR (Basis_HCURL_HEX_In_FEM): Operator type not implemented" );
        // break;commented out since exception is thrown
      }
      }
    }
  }

  // -------------------------------------------------------------------------------------

  template<typename DT, typename OT, typename PT>
  Basis_HCURL_HEX_In_FEM<DT,OT,PT>::
  Basis_HCURL_HEX_In_FEM( const ordinal_type order,
                          const EPointType   pointType ) {

    INTREPID2_TEST_FOR_EXCEPTION( !(pointType == POINTTYPE_EQUISPACED ||
                                    pointType == POINTTYPE_WARPBLEND), std::invalid_argument,
                                  ">>> ERROR (Basis_HCURL_HEX_In_FEM): pointType must be either equispaced or warpblend.");

    // this should be created in host and vinv should be deep copied into device space
    Basis_HGRAD_LINE_Cn_FEM<DT,OT,PT> lineBasis( order, pointType );
    Basis_HVOL_LINE_Cn_FEM<DT,OT,PT> bubbleBasis( order - 1, POINTTYPE_GAUSS );

    const ordinal_type
      cardLine = lineBasis.getCardinality(),
      cardBubble = bubbleBasis.getCardinality();

    this->vinvLine_   = Kokkos::DynRankView<typename ScalarViewType::value_type,DT>("Hcurl::Hex::In::vinvLine", cardLine, cardLine);
    this->vinvBubble_ = Kokkos::DynRankView<typename ScalarViewType::value_type,DT>("Hcurl::Hex::In::vinvBubble", cardBubble, cardBubble);

    lineBasis.getVandermondeInverse(this->vinvLine_);
    bubbleBasis.getVandermondeInverse(this->vinvBubble_);

    const ordinal_type spaceDim = 3;
    this->basisCardinality_     = 3*cardLine*cardLine*cardBubble;
    this->basisDegree_          = order;
    this->basisCellTopologyKey_ = shards::Hexahedron<8>::key;
    this->basisType_            = BASIS_FEM_LAGRANGIAN;
    this->basisCoordinates_     = COORDINATES_CARTESIAN;
    this->functionSpace_        = FUNCTION_SPACE_HCURL;
    pointType_ = pointType;

    // initialize tags
    {
      // Basis-dependent initializations
      const ordinal_type tagSize  = 4;        // size of DoF tag, i.e., number of fields in the tag
      const ordinal_type posScDim = 0;        // position in the tag, counting from 0, of the subcell dim
      const ordinal_type posScOrd = 1;        // position in the tag, counting from 0, of the subcell ordinal
      const ordinal_type posDfOrd = 2;        // position in the tag, counting from 0, of DoF ordinal relative to the subcell

      // Note: the only reason why equispaced can't support higher order than Parameters::MaxOrder appears to be the fact that the tags below get stored into a fixed-length array.
      // TODO: relax the maximum order requirement by setting up tags in a different container, perhaps directly into an OrdinalTypeArray1DHost (tagView, below).  (As of this writing (1/25/22), looks like other nodal bases do this in a similar way -- those should be fixed at the same time; maybe search for Parameters::MaxOrder.)
      INTREPID2_TEST_FOR_EXCEPTION( order > Parameters::MaxOrder, std::invalid_argument, "polynomial order exceeds the max supported by this class");
      
      // An array with local DoF tags assigned to the basis functions, in the order of their local enumeration
      constexpr ordinal_type maxCardLine = Parameters::MaxOrder + 1;
      ordinal_type tags[3*maxCardLine*maxCardLine*maxCardLine][4];

      const ordinal_type edge_x[2][2] = { {0, 4}, {2, 6} };
      const ordinal_type edge_y[2][2] = { {3, 7}, {1, 5} };
      const ordinal_type edge_z[2][2] = { {8,11}, {9,10} };

      const ordinal_type face_yz[2] = {3, 1};
      const ordinal_type face_xz[2] = {0, 2};
      const ordinal_type face_xy[2] = {4, 5};

      {
        ordinal_type idx = 0;

        ///
        /// Product rule: z -> y -> x, x-tangent first
        ///

        // since there are x/y components in the interior
        // dof sum should be computed before the information
        // is assigned to tags
        const ordinal_type
          face_ndofs_per_direction = (cardLine-2)*cardBubble,
          face_ndofs = 2*face_ndofs_per_direction,
          intr_ndofs_per_direction = (cardLine-2)*(cardLine-2)*cardBubble,
          intr_ndofs = 3*intr_ndofs_per_direction;

        // x component (lineBasis(z) lineBasis(y) bubbleBasis(x))
        for (ordinal_type k=0;k<cardLine;++k) { // z
          const auto tag_z = lineBasis.getDofTag(k);
          for (ordinal_type j=0;j<cardLine;++j) { // y
            const auto tag_y = lineBasis.getDofTag(j);
            for (ordinal_type i=0;i<cardBubble;++i,++idx) { // x
              const auto tag_x = bubbleBasis.getDofTag(i);

              if (tag_x(0) == 1 && tag_y(0) == 0 && tag_z(0) == 0) {
                // edge: x edge, y vert, z vert
                tags[idx][0] = 1; // edge dof
                tags[idx][1] = edge_x[tag_y(1)][tag_z(1)]; // edge id
                tags[idx][2] = tag_x(2); // local dof id
                tags[idx][3] = tag_x(3); // total number of dofs in this edge
              } else if (tag_x(0) == 1 && tag_y(0) == 0 && tag_z(0) == 1) {
                // face, x edge, y vert, z edge
                tags[idx][0] = 2; // face dof
                tags[idx][1] = face_xz[tag_y(1)]; // face id
                tags[idx][2] = tag_x(2) + tag_x(3)*tag_z(2); // local dof id
                tags[idx][3] = face_ndofs; // total number of dofs in this face
              } else if (tag_x(0) == 1 && tag_y(0) == 1 && tag_z(0) == 0) {
                // face, x edge, y edge, z vert
                tags[idx][0] = 2; // face dof
                tags[idx][1] = face_xy[tag_z(1)]; // face id
                tags[idx][2] = tag_x(2) + tag_x(3)*tag_y(2); // local dof id
                tags[idx][3] = face_ndofs; // total number of dofs in this face
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

        // y component (lineBasis(z) bubbleBasis(y) lineBasis(x))
        for (ordinal_type k=0;k<cardLine;++k) { // z
          const auto tag_z = lineBasis.getDofTag(k);
          for (ordinal_type j=0;j<cardBubble;++j) { // y
            const auto tag_y = bubbleBasis.getDofTag(j);
            for (ordinal_type i=0;i<cardLine;++i,++idx) { // x
              const auto tag_x = lineBasis.getDofTag(i);

              if (tag_x(0) == 0 && tag_y(0) == 1 && tag_z(0) == 0) {
                // edge: x vert, y edge, z vert
                tags[idx][0] = 1; // edge dof
                tags[idx][1] = edge_y[tag_x(1)][tag_z(1)]; // edge id
                tags[idx][2] = tag_y(2); // local dof id
                tags[idx][3] = tag_y(3); // total number of dofs in this vertex
              } else if (tag_x(0) == 0 && tag_y(0) == 1 && tag_z(0) == 1) {
                // face, x vert, y edge, z edge
                tags[idx][0] = 2; // face dof
                tags[idx][1] = face_yz[tag_x(1)]; // face id
                tags[idx][2] = tag_y(2) + tag_y(3)*tag_z(2); // local dof id
                tags[idx][3] = face_ndofs; // total number of dofs in this vertex
              } else if (tag_x(0) == 1 && tag_y(0) == 1 && tag_z(0) == 0) {
                // face, x edge, y edge, z vert
                tags[idx][0] = 2; // face dof
                tags[idx][1] = face_xy[tag_z(1)]; // face id
                tags[idx][2] = face_ndofs_per_direction + tag_x(2) + tag_x(3)*tag_y(2); // local dof id
                tags[idx][3] = face_ndofs; // total number of dofs in this vertex
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
        for (ordinal_type k=0;k<cardBubble;++k) { // y
          const auto tag_z = bubbleBasis.getDofTag(k);
          for (ordinal_type j=0;j<cardLine;++j) { // z
            const auto tag_y = lineBasis.getDofTag(j);
            for (ordinal_type i=0;i<cardLine;++i,++idx) { // x
              const auto tag_x = lineBasis.getDofTag(i);

              if (tag_x(0) == 0 && tag_y(0) == 0 && tag_z(0) == 1) {
                // edge: x vert, y vert, z edge
                tags[idx][0] = 1; // edge dof
                tags[idx][1] = edge_z[tag_x(1)][tag_y(1)]; // edge id
                tags[idx][2] = tag_z(2); // local dof id
                tags[idx][3] = tag_z(3); // total number of dofs in this vertex
              } else if (tag_x(0) == 0 && tag_y(0) == 1 && tag_z(0) == 1) {
                // face, x vert, y edge, z edge
                tags[idx][0] = 2; // face dof
                tags[idx][1] = face_yz[tag_x(1)]; // face id
                tags[idx][2] = face_ndofs_per_direction + tag_y(2) + tag_y(3)*tag_z(2); // local dof id
                tags[idx][3] = face_ndofs; // total number of dofs in this vertex
              } else if (tag_x(0) == 1 && tag_y(0) == 0 && tag_z(0) == 1) {
                // face, x edge, y vert, z edge
                tags[idx][0] = 2; // face dof
                tags[idx][1] = face_xz[tag_y(1)]; // face id
                tags[idx][2] = face_ndofs_per_direction + tag_x(2) + tag_x(3)*tag_z(2); // local dof id
                tags[idx][3] = face_ndofs; // total number of dofs in this vertex
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
                                      ">>> ERROR (Basis_HCURL_HEX_In_FEM): " \
                                      "counted tag index is not same as cardinality." );
      }

      OrdinalTypeArray1DHost tagView(&tags[0][0], this->basisCardinality_*4);

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
    Kokkos::DynRankView<typename ScalarViewType::value_type,typename DT::execution_space::array_layout,Kokkos::HostSpace>
      dofCoordsHost("dofCoordsHost", this->basisCardinality_, spaceDim);

    // dofCoords on host and create its mirror view to device
    Kokkos::DynRankView<typename ScalarViewType::value_type,typename DT::execution_space::array_layout,Kokkos::HostSpace>
      dofCoeffsHost("dofCoeffsHost", this->basisCardinality_, spaceDim);

    Kokkos::DynRankView<typename ScalarViewType::value_type,DT>
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
      for (ordinal_type k=0;k<cardLine;++k) { // z
        for (ordinal_type j=0;j<cardLine;++j) { // y
          for (ordinal_type i=0;i<cardBubble;++i,++idx) { // x
            dofCoordsHost(idx,0) = dofCoordsBubbleHost(i,0);
            dofCoordsHost(idx,1) = dofCoordsLineHost(j,0);
            dofCoordsHost(idx,2) = dofCoordsLineHost(k,0);
            dofCoeffsHost(idx,0) = 1.0;
          }
        }
      }

      // y component (lineBasis(z) bubbleBasis(y) lineBasis(x))
      for (ordinal_type k=0;k<cardLine;++k) { // z
        for (ordinal_type j=0;j<cardBubble;++j) { // y
          for (ordinal_type i=0;i<cardLine;++i,++idx) { // x
            dofCoordsHost(idx,0) = dofCoordsLineHost(i,0);
            dofCoordsHost(idx,1) = dofCoordsBubbleHost(j,0);
            dofCoordsHost(idx,2) = dofCoordsLineHost(k,0);
            dofCoeffsHost(idx,1) = 1.0;
          }
        }
      }

      // z component (bubbleBasis(z) lineBasis(y) lineBasis(x))
      for (ordinal_type k=0;k<cardBubble;++k) { // z
        for (ordinal_type j=0;j<cardLine;++j) { // y
          for (ordinal_type i=0;i<cardLine;++i,++idx) { // x
            dofCoordsHost(idx,0) = dofCoordsLineHost(i,0);
            dofCoordsHost(idx,1) = dofCoordsLineHost(j,0);
            dofCoordsHost(idx,2) = dofCoordsBubbleHost(k,0);
            dofCoeffsHost(idx,2) = 1.0;
          }
        }
      }
    }

    this->dofCoords_ = Kokkos::create_mirror_view(typename DT::memory_space(), dofCoordsHost);
    Kokkos::deep_copy(this->dofCoords_, dofCoordsHost);

    this->dofCoeffs_ = Kokkos::create_mirror_view(typename DT::memory_space(), dofCoeffsHost);
    Kokkos::deep_copy(this->dofCoeffs_, dofCoeffsHost);
  }
}

#endif
