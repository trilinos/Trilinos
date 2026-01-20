// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_HGRAD_QUAD_Cn_FEMDef.hpp
    \brief  Definition file for FEM basis functions of degree n for H(grad) functions on QUAD cells.
    \author Created by R. Kirby.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_HGRAD_QUAD_CN_FEM_DEF_HPP__
#define __INTREPID2_HGRAD_QUAD_CN_FEM_DEF_HPP__

namespace Intrepid2 {

  // -------------------------------------------------------------------------------------
  namespace Impl {

    template<EOperator OpType>
    template<typename OutputViewType,
             typename InputViewType,
             typename WorkViewType,
             typename VinvViewType>
    KOKKOS_INLINE_FUNCTION
    void
    Basis_HGRAD_QUAD_Cn_FEM::Serial<OpType>::
    getValues(       OutputViewType output,
               const InputViewType  input,
                     WorkViewType   work,
               const VinvViewType   vinv,
               const ordinal_type   operatorDn ) {
      ordinal_type opDn = operatorDn;
      
      const ordinal_type cardLine = vinv.extent(0);
      const ordinal_type npts = input.extent(0);

      typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;
      const auto input_x = Kokkos::subview(input, Kokkos::ALL(), range_type(0,1));
      const auto input_y = Kokkos::subview(input, Kokkos::ALL(), range_type(1,2));

      const int dim_s = get_dimension_scalar(input);
      auto ptr0 = work.data();
      auto ptr1 = work.data()+cardLine*npts*dim_s;
      auto ptr2 = work.data()+2*cardLine*npts*dim_s;
      
      typedef typename Kokkos::DynRankView<typename InputViewType::value_type, typename WorkViewType::memory_space> ViewType;
      
      switch (OpType) {
      case OPERATOR_VALUE: {
        ViewType work_line = createMatchingUnmanagedView<ViewType>(input, ptr0, cardLine, npts);
        ViewType output_x = createMatchingUnmanagedView<ViewType>(input, ptr1, cardLine, npts);
        ViewType output_y = createMatchingUnmanagedView<ViewType>(input, ptr2, cardLine, npts);
        
        Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
          getValues(output_x, input_x, work_line, vinv);

        Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
          getValues(output_y, input_y, work_line, vinv);

        // tensor product
        ordinal_type idx = 0;
        for (ordinal_type j=0;j<cardLine;++j) // y
          for (ordinal_type i=0;i<cardLine;++i,++idx)  // x
            for (ordinal_type k=0;k<npts;++k)
              output.access(idx,k) = output_x.access(i,k)*output_y.access(j,k);
        break;
      }
      case OPERATOR_CURL: {
        for (auto l=0;l<2;++l) {
          ViewType work_line = createMatchingUnmanagedView<ViewType>(input, ptr0, cardLine, npts);
          
          ViewType output_x, output_y;
          
          typename WorkViewType::value_type s = 0.0;
          if (l) {
            // l = 1
            output_x = createMatchingUnmanagedView<ViewType>(input, ptr1, cardLine, npts, 1);
            Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_Dn>::
              getValues(output_x, input_x, work_line, vinv, 1);                           

            output_y = createMatchingUnmanagedView<ViewType>(input, ptr2, cardLine, npts);
            Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
              getValues(output_y, input_y, work_line, vinv);                           

            s = -1.0;
          } else {
            // l = 0
            output_x = createMatchingUnmanagedView<ViewType>(input, ptr1, cardLine, npts);
            Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
              getValues(output_x, input_x, work_line, vinv);                           

            output_y = createMatchingUnmanagedView<ViewType>(input, ptr2, cardLine, npts, 1);
            Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_Dn>::
              getValues(output_y, input_y, work_line, vinv, 1);                           

            s = 1.0;
          }

          // tensor product (extra dimension of ouput x and y are ignored)
          ordinal_type idx = 0;
          for (ordinal_type j=0;j<cardLine;++j) // y
            for (ordinal_type i=0;i<cardLine;++i,++idx)  // x
              for (ordinal_type k=0;k<npts;++k)
                output.access(idx,k,l) = s*output_x.access(i,k,0)*output_y.access(j,k,0);
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
      case OPERATOR_D10:
        opDn = getOperatorOrder(OpType);
      case OPERATOR_Dn: {
        const auto dkcard = opDn + 1;
        for (auto l=0;l<dkcard;++l) {
          ViewType work_line = createMatchingUnmanagedView<ViewType>(input, ptr0, cardLine, npts);
          
          ViewType output_x, output_y;
          
          const auto mult_x = opDn - l;
          const auto mult_y = l;
          
          if (mult_x) {
            output_x = createMatchingUnmanagedView<ViewType>(input, ptr1, cardLine, npts, 1);
            Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_Dn>::
              getValues(output_x, input_x, work_line, vinv, mult_x);                           
          } else {
            output_x = createMatchingUnmanagedView<ViewType>(input, ptr1, cardLine, npts);
            Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
              getValues(output_x, input_x, work_line, vinv);                           
          }

          if (mult_y) {
            output_y = createMatchingUnmanagedView<ViewType>(input, ptr2, cardLine, npts, 1);
            Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_Dn>::
              getValues(output_y, input_y, work_line, vinv, mult_y);                           
          } else {
            output_y = createMatchingUnmanagedView<ViewType>(input, ptr2, cardLine, npts);
            Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
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
                                  ">>> ERROR: (Intrepid2::Basis_HGRAD_QUAD_Cn_FEM::Serial::getValues) operator is not supported" );
      }
      }
    }
    
    template<typename DT, ordinal_type numPtsPerEval,
             typename outputValueValueType, class ...outputValueProperties,
             typename inputPointValueType,  class ...inputPointProperties,
             typename vinvValueType,        class ...vinvProperties>
    void
    Basis_HGRAD_QUAD_Cn_FEM::
    getValues( const typename DT::execution_space& space,
                     Kokkos::DynRankView<outputValueValueType,outputValueProperties...> outputValues,
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
      Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(space, 0, loopSize);

      //typedef typename inputPointViewType::value_type inputPointType;

      const ordinal_type cardinality = outputValues.extent(0);
      const ordinal_type cardLine = std::sqrt(cardinality);
      const ordinal_type workSize = 3*cardLine;
      
      auto work = createMatchingDynRankView(inputPoints, "Basis_HGRAD_QUAD_Cn_FEM::getValues::work", workSize, inputPoints.extent(0));

      switch (operatorType) {
      case OPERATOR_VALUE: {
        typedef Functor<outputValueViewType,inputPointViewType,vinvViewType,decltype(work),
            OPERATOR_VALUE,numPtsPerEval> FunctorType;
        Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints, vinv, work) );
        break;
      }
      case OPERATOR_CURL: {
        typedef Functor<outputValueViewType,inputPointViewType,vinvViewType,decltype(work),
            OPERATOR_CURL,numPtsPerEval> FunctorType;
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
        typedef Functor<outputValueViewType,inputPointViewType,vinvViewType,decltype(work),
            OPERATOR_Dn,numPtsPerEval> FunctorType;
        Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints, vinv, work,
                                                  getOperatorOrder(operatorType)) );
        break;
      }
      default: {
        INTREPID2_TEST_FOR_EXCEPTION( true , std::invalid_argument,
                                      ">>> ERROR (Basis_HGRAD_QUAD_Cn_FEM): Operator type not implemented" );
        // break;commented out because exception 
      }
      }
    }
  }

  // -------------------------------------------------------------------------------------
  template<typename DT, typename OT, typename PT>
  Basis_HGRAD_QUAD_Cn_FEM<DT,OT,PT>::
  Basis_HGRAD_QUAD_Cn_FEM( const ordinal_type order,
                           const EPointType   pointType ) {
    // INTREPID2_TEST_FOR_EXCEPTION( !(pointType == POINTTYPE_EQUISPACED ||
    //                                 pointType == POINTTYPE_WARPBLEND), std::invalid_argument,
    //                               ">>> ERROR (Basis_HGRAD_QUAD_Cn_FEM): pointType must be either equispaced or warpblend." );

    // this should be in host
    Basis_HGRAD_LINE_Cn_FEM<DT,OT,PT> lineBasis( order, pointType );
    const auto cardLine = lineBasis.getCardinality();
    
    this->vinv_ = Kokkos::DynRankView<typename ScalarViewType::value_type,DT>("Hgrad::Quad::Cn::vinv", cardLine, cardLine);         
    lineBasis.getVandermondeInverse(this->vinv_);
    
    const ordinal_type spaceDim = 2;
    this->basisCardinality_     = cardLine*cardLine;
    this->basisDegree_          = order;
    this->basisCellTopologyKey_ = shards::Quadrilateral<4>::key;
    this->basisType_            = BASIS_FEM_LAGRANGIAN;
    this->basisCoordinates_     = COORDINATES_CARTESIAN;
    this->functionSpace_        = FUNCTION_SPACE_HGRAD;
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
      ordinal_type tags[maxCardLine*maxCardLine][4];

      const ordinal_type vert[2][2] = { {0,1}, {3,2} }; //[y][x]

      const ordinal_type edge_x[2] = {0,2}; 
      const ordinal_type edge_y[2] = {3,1}; 
      {
        ordinal_type idx = 0;
        for (ordinal_type j=0;j<cardLine;++j) { // y      
          const auto tag_y = lineBasis.getDofTag(j);
          for (ordinal_type i=0;i<cardLine;++i,++idx) { // x
            const auto tag_x = lineBasis.getDofTag(i);          
            
            if (tag_x(0) == 0 && tag_y(0) == 0) {
              // vertices
              tags[idx][0] = 0; // vertex dof
              tags[idx][1] = vert[tag_y(1)][tag_x(1)]; // vertex id
              tags[idx][2] = 0; // local dof id
              tags[idx][3] = 1; // total number of dofs in this vertex
            } else if (tag_x(0) == 1 && tag_y(0) == 0) {
              // edge: x edge, y vert
              tags[idx][0] = 1; // edge dof
              tags[idx][1] = edge_x[tag_y(1)];
              tags[idx][2] = tag_x(2); // local dof id
              tags[idx][3] = tag_x(3); // total number of dofs in this vertex
            } else if (tag_x(0) == 0 && tag_y(0) == 1) {
              // edge: x vert, y edge
              tags[idx][0] = 1; // edge dof
              tags[idx][1] = edge_y[tag_x(1)];
              tags[idx][2] = tag_y(2); // local dof id
              tags[idx][3] = tag_y(3); // total number of dofs in this vertex
            } else {
              // interior
              tags[idx][0] = 2; // interior dof
              tags[idx][1] = 0;
              tags[idx][2] = tag_x(2) + tag_x(3)*tag_y(2); // local dof id
              tags[idx][3] = tag_x(3)*tag_y(3); // total number of dofs in this vertex
            }
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

  template<typename DT, typename OT, typename PT>
  void 
  Basis_HGRAD_QUAD_Cn_FEM<DT,OT,PT>::getScratchSpaceSize(       
                                    ordinal_type& perTeamSpaceSize,
                                    ordinal_type& perThreadSpaceSize,
                              const PointViewType inputPoints,
                              const EOperator operatorType) const {
    perTeamSpaceSize = 0;
    perThreadSpaceSize = 3*this->vinv_.extent(0)*get_dimension_scalar(inputPoints)*sizeof(typename BasisBase::scalarType);
  }

  template<typename DT, typename OT, typename PT>
  KOKKOS_INLINE_FUNCTION
  void 
  Basis_HGRAD_QUAD_Cn_FEM<DT,OT,PT>::getValues(       
          OutputViewType outputValues,
      const PointViewType  inputPoints,
      const EOperator operatorType,
      const typename Kokkos::TeamPolicy<typename DT::execution_space>::member_type& team_member,
      const typename DT::execution_space::scratch_memory_space & scratchStorage, 
      const ordinal_type subcellDim,
      const ordinal_type subcellOrdinal) const {

      INTREPID2_TEST_FOR_ABORT( !((subcellDim == -1) && (subcellOrdinal == -1)),
        ">>> ERROR: (Intrepid2::Basis_HGRAD_QUAD_Cn_FEM::getValues), The capability of selecting subsets of basis functions has not been implemented yet.");

      const int numPoints = inputPoints.extent(0);
      using ScalarType = typename ScalarTraits<typename PointViewType::value_type>::scalar_type;
      using WorkViewType = Kokkos::DynRankView< ScalarType,typename DT::execution_space::scratch_memory_space,Kokkos::MemoryTraits<Kokkos::Unmanaged> >;
      ordinal_type sizePerPoint = 3*this->vinv_.extent(0)*get_dimension_scalar(inputPoints);
      WorkViewType workView(scratchStorage, sizePerPoint*team_member.team_size());
      using range_type = Kokkos::pair<ordinal_type,ordinal_type>;

      switch(operatorType) {
        case OPERATOR_VALUE:
          Kokkos::parallel_for (Kokkos::TeamThreadRange (team_member, numPoints), [=, &vinv_ = this->vinv_] (ordinal_type& pt) {
            auto       output = Kokkos::subview( outputValues, Kokkos::ALL(), range_type  (pt,pt+1), Kokkos::ALL() );
            const auto input  = Kokkos::subview( inputPoints,                 range_type(pt, pt+1), Kokkos::ALL() );
            WorkViewType  work(workView.data() + sizePerPoint*team_member.team_rank(), sizePerPoint);
            Impl::Basis_HGRAD_QUAD_Cn_FEM::Serial<OPERATOR_VALUE>::getValues( output, input, work, vinv_ );
          });
          break;
        case OPERATOR_GRAD:
          Kokkos::parallel_for (Kokkos::TeamThreadRange (team_member, numPoints), [=, &vinv_ = this->vinv_] (ordinal_type& pt) {
            auto       output = Kokkos::subview( outputValues, Kokkos::ALL(), range_type(pt,pt+1), Kokkos::ALL() );
            const auto input  = Kokkos::subview( inputPoints,                 range_type(pt,pt+1), Kokkos::ALL() );
            WorkViewType  work(workView.data() + sizePerPoint*team_member.team_rank(), sizePerPoint);
            Impl::Basis_HGRAD_QUAD_Cn_FEM::Serial<OPERATOR_GRAD>::getValues( output, input, work, vinv_ );
          });
          break;
          case OPERATOR_CURL:
          Kokkos::parallel_for (Kokkos::TeamThreadRange (team_member, numPoints), [=, &vinv_ = this->vinv_] (ordinal_type& pt) {
            auto       output = Kokkos::subview( outputValues, Kokkos::ALL(), range_type(pt,pt+1), Kokkos::ALL() );
            const auto input  = Kokkos::subview( inputPoints,                 range_type(pt,pt+1), Kokkos::ALL() );
            WorkViewType  work(workView.data() + sizePerPoint*team_member.team_rank(), sizePerPoint);
            Impl::Basis_HGRAD_QUAD_Cn_FEM::Serial<OPERATOR_CURL>::getValues( output, input, work, vinv_ );
          });
          break;
        default: {          
          INTREPID2_TEST_FOR_ABORT( true,
            ">>> ERROR (Basis_HGRAD_QUAD_Cn_FEM): getValues not implemented for this operator");
          }
    }
  }

} // namespace Intrepid2

#endif
