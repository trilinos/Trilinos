// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_HVOL_QUAD_Cn_FEMDef.hpp
    \brief  Definition file for FEM basis functions of degree n for H(vol) functions on QUAD
    \author Created by M. Perego, based on the Intrepid2::HGRAD_QUAD_Cn_FEM class
*/

#ifndef __INTREPID2_HVOL_QUAD_CN_FEM_DEF_HPP__
#define __INTREPID2_HVOL_QUAD_CN_FEM_DEF_HPP__

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
    Basis_HVOL_QUAD_Cn_FEM::Serial<opType>::
    getValues(       OutputViewType output,
               const inputViewType  input,
                     workViewType   work,
               const vinvViewType   vinv,
               const ordinal_type   operatorDn ) {
      ordinal_type opDn = operatorDn;
      
      const ordinal_type cardLine = vinv.extent(0);
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
        viewType work_line(Kokkos::view_wrap(ptr0, vcprop), cardLine, npts);
        viewType output_x(Kokkos::view_wrap(ptr1, vcprop), cardLine, npts);
        viewType output_y(Kokkos::view_wrap(ptr2, vcprop), cardLine, npts);
        
        Impl::Basis_HVOL_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
          getValues(output_x, input_x, work_line, vinv);

        Impl::Basis_HVOL_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
          getValues(output_y, input_y, work_line, vinv);

        // tensor product
        ordinal_type idx = 0;
        for (ordinal_type j=0;j<cardLine;++j) // y
          for (ordinal_type i=0;i<cardLine;++i,++idx)  // x
            for (ordinal_type k=0;k<npts;++k)
              output.access(idx,k) = output_x.access(i,k)*output_y.access(j,k);
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
      case OPERATOR_D10:
        opDn = getOperatorOrder(opType);
      case OPERATOR_Dn: {
        const auto dkcard = opDn + 1;
        for (auto l=0;l<dkcard;++l) {
          viewType work_line(Kokkos::view_wrap(ptr0, vcprop), cardLine, npts);
          
          viewType output_x, output_y;
          
          const auto mult_x = opDn - l;
          const auto mult_y = l;
          
          if (mult_x) {
            output_x = viewType(Kokkos::view_wrap(ptr1, vcprop), cardLine, npts, 1);
            Impl::Basis_HVOL_LINE_Cn_FEM::Serial<OPERATOR_Dn>::
              getValues(output_x, input_x, work_line, vinv, mult_x);                           
          } else {
            output_x = viewType(Kokkos::view_wrap(ptr1, vcprop), cardLine, npts);
            Impl::Basis_HVOL_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
              getValues(output_x, input_x, work_line, vinv);                           
          }

          if (mult_y) {
            output_y = viewType(Kokkos::view_wrap(ptr2, vcprop), cardLine, npts, 1);
            Impl::Basis_HVOL_LINE_Cn_FEM::Serial<OPERATOR_Dn>::
              getValues(output_y, input_y, work_line, vinv, mult_y);                           
          } else {
            output_y = viewType(Kokkos::view_wrap(ptr2, vcprop), cardLine, npts);
            Impl::Basis_HVOL_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
              getValues(output_y, input_y, work_line, vinv);                           
          }

          // tensor product (extra dimension of ouput x and y are ignored)
          ordinal_type idx = 0;
          for (ordinal_type j=0;j<cardLine;++j) // y
            for (ordinal_type i=0;i<cardLine;++i,++idx)  // x
              for (ordinal_type k=0;k<npts;++k)
                output.access(idx,k,l) = output_x.access(i,k,0)*output_y.access(j,k,0);
        }
        break;
      }
      default: {
        INTREPID2_TEST_FOR_ABORT( true,
                                  ">>> ERROR: (Intrepid2::Basis_HVOL_QUAD_Cn_FEM::Serial::getValues) operator is not supported" );
      }
      }
    }
    
    template<typename DT, ordinal_type numPtsPerEval,
             typename outputValueValueType, class ...outputValueProperties,
             typename inputPointValueType,  class ...inputPointProperties,
             typename vinvValueType,        class ...vinvProperties>
    void
    Basis_HVOL_QUAD_Cn_FEM::
    getValues(       Kokkos::DynRankView<outputValueValueType,outputValueProperties...> outputValues,
               const Kokkos::DynRankView<inputPointValueType, inputPointProperties...>  inputPoints,
               const Kokkos::DynRankView<vinvValueType,       vinvProperties...>        vinv,
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
      const ordinal_type cardLine = std::sqrt(cardinality);
      const ordinal_type workSize = 3*cardLine;

      auto vcprop = Kokkos::common_view_alloc_prop(inputPoints);
      typedef typename Kokkos::DynRankView< inputPointType, typename inputPointViewType::memory_space> workViewType;
      workViewType  work(Kokkos::view_alloc("Basis_HVOL_QUAD_Cn_FEM::getValues::work", vcprop), workSize, inputPoints.extent(0));

      switch (operatorType) {
      case OPERATOR_VALUE: {
        typedef Functor<outputValueViewType,inputPointViewType,vinvViewType,workViewType,
            OPERATOR_VALUE,numPtsPerEval> FunctorType;
        Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints, vinv, work) );
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
        typedef Functor<outputValueViewType,inputPointViewType,vinvViewType,workViewType,
            OPERATOR_Dn,numPtsPerEval> FunctorType;
        Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints, vinv, work,
                                                  getOperatorOrder(operatorType)) );
        break;
      }
      default: {
        INTREPID2_TEST_FOR_EXCEPTION( true , std::invalid_argument,
                                      ">>> ERROR (Basis_HVOL_QUAD_Cn_FEM): Operator type not implemented" );
        // break;commented out because exception 
      }
      }
    }
  }

  // -------------------------------------------------------------------------------------
  template<typename DT, typename OT, typename PT>
  Basis_HVOL_QUAD_Cn_FEM<DT,OT,PT>::
  Basis_HVOL_QUAD_Cn_FEM( const ordinal_type order,
                           const EPointType   pointType ) {
    // INTREPID2_TEST_FOR_EXCEPTION( !(pointType == POINTTYPE_EQUISPACED ||
    //                                 pointType == POINTTYPE_WARPBLEND), std::invalid_argument,
    //                               ">>> ERROR (Basis_HVOL_QUAD_Cn_FEM): pointType must be either equispaced or warpblend." );

    // this should be in host
    Basis_HVOL_LINE_Cn_FEM<DT,OT,PT> lineBasis( order, pointType );
    const auto cardLine = lineBasis.getCardinality();
    
    this->pointType_ = pointType;
    this->vinv_ = Kokkos::DynRankView<typename ScalarViewType::value_type,DT>("HVOL::Quad::Cn::vinv", cardLine, cardLine);
    lineBasis.getVandermondeInverse(this->vinv_);

    const ordinal_type spaceDim = 2;
    this->basisCardinality_     = cardLine*cardLine;
    this->basisDegree_          = order;
    this->basisCellTopologyKey_ = shards::Quadrilateral<4>::key;
    this->basisType_            = BASIS_FEM_LAGRANGIAN;
    this->basisCoordinates_     = COORDINATES_CARTESIAN;
    this->functionSpace_        = FUNCTION_SPACE_HVOL;

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
        for (ordinal_type j=0;j<cardLine;++j) { // y      
          const auto tag_y = lineBasis.getDofTag(j);
          for (ordinal_type i=0;i<cardLine;++i,++idx) { // x
            const auto tag_x = lineBasis.getDofTag(i);          
            
            // interior
            tags[idx][0] = 2; // interior dof
            tags[idx][1] = 0;
            tags[idx][2] = tag_x(2) + tag_x(3)*tag_y(2); // local dof id
            tags[idx][3] = tag_x(3)*tag_y(3); // total number of dofs in this vertex
          }
        }
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

    Kokkos::DynRankView<typename ScalarViewType::value_type,DT>
      dofCoordsLine("dofCoordsLine", cardLine, 1);

    lineBasis.getDofCoords(dofCoordsLine);
    auto dofCoordsLineHost = Kokkos::create_mirror_view(dofCoordsLine);
    Kokkos::deep_copy(dofCoordsLineHost, dofCoordsLine);
    {
      ordinal_type idx = 0;
      for (ordinal_type j=0;j<cardLine;++j) { // y      
        for (ordinal_type i=0;i<cardLine;++i,++idx) { // x
          dofCoordsHost(idx,0) = dofCoordsLineHost(i,0);
          dofCoordsHost(idx,1) = dofCoordsLineHost(j,0);
        }
      }
    }

    this->dofCoords_ = Kokkos::create_mirror_view(typename DT::memory_space(), dofCoordsHost);
    Kokkos::deep_copy(this->dofCoords_, dofCoordsHost);
  }
  
}

#endif
