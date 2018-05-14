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

/** \file   Intrepid2_HVOL_HEX_Cn_FEMDef.hpp
    \brief  Definition file for FEM basis functions of degree n for H(vol) functions on HEX cells
    \author Created by M. Perego, based on the Intrepid2::HGRAD_HEX_Cn_FEM class
*/

#ifndef __INTREPID2_HVOL_HEX_CN_FEMDEF_HPP__
#define __INTREPID2_HVOL_HEX_CN_FEMDEF_HPP__

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
    Basis_HVOL_HEX_Cn_FEM::Serial<opType>::
    getValues(       outputViewType output,
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
      const auto input_z = Kokkos::subview(input, Kokkos::ALL(), range_type(2,3));

      const ordinal_type dim_s = get_dimension_scalar(work);
      auto ptr0 = work.data();
      auto ptr1 = work.data()+cardLine*npts*dim_s;
      auto ptr2 = work.data()+2*cardLine*npts*dim_s;
      auto ptr3 = work.data()+3*cardLine*npts*dim_s;      
      
      typedef typename Kokkos::DynRankView<typename workViewType::value_type, typename workViewType::memory_space> viewType;
      auto vcprop = Kokkos::common_view_alloc_prop(work);

      switch (opType) {
      case OPERATOR_VALUE: {
        viewType work_line(Kokkos::view_wrap(ptr0, vcprop), cardLine, npts);
        viewType output_x(Kokkos::view_wrap(ptr1, vcprop), cardLine, npts);
        viewType output_y(Kokkos::view_wrap(ptr2, vcprop), cardLine, npts);
        viewType output_z(Kokkos::view_wrap(ptr3, vcprop), cardLine, npts);
        
        Impl::Basis_HVOL_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
          getValues(output_x, input_x, work_line, vinv);

        Impl::Basis_HVOL_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
          getValues(output_y, input_y, work_line, vinv);

        Impl::Basis_HVOL_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
          getValues(output_z, input_z, work_line, vinv);

        // tensor product
        ordinal_type idx = 0;
        for (ordinal_type k=0;k<cardLine;++k) // z
          for (ordinal_type j=0;j<cardLine;++j) // y
            for (ordinal_type i=0;i<cardLine;++i,++idx)  // x
              for (ordinal_type l=0;l<npts;++l)
                output.access(idx,l) = output_x.access(i,l)*output_y.access(j,l)*output_z.access(k,l);
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
        const ordinal_type dkcard = opDn + 1;

        ordinal_type d = 0;
        for (ordinal_type l1=0;l1<dkcard;++l1)
          for (ordinal_type l0=0;l0<(l1+1);++l0) {
            const ordinal_type mult_x = (opDn - l1);
            const ordinal_type mult_y = l1 - l0;
            const ordinal_type mult_z = l0;

            //std::cout << " l0, l1 = " << l0 << " " << l1 << std::endl;
            //std::cout << " x , y , z = " << mult_x << " " << mult_y << " " << mult_z << std::endl;

            if (mult_x < 0) {
              // pass
            } else {              
              viewType work_line(Kokkos::view_wrap(ptr0, vcprop), cardLine, npts);
              decltype(work_line)  output_x, output_y, output_z;
                  
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

              if (mult_z) {
                output_z = viewType(Kokkos::view_wrap(ptr3, vcprop), cardLine, npts, 1);
                Impl::Basis_HVOL_LINE_Cn_FEM::Serial<OPERATOR_Dn>::
                  getValues(output_z, input_z, work_line, vinv, mult_z);
              } else {
                output_z = viewType(Kokkos::view_wrap(ptr3, vcprop), cardLine, npts);
                Impl::Basis_HVOL_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
                  getValues(output_z, input_z, work_line, vinv);
              }

              // tensor product (extra dimension of ouput x,y and z are ignored)              
              ordinal_type idx = 0;
              for (ordinal_type k=0;k<cardLine;++k) // z
                for (ordinal_type j=0;j<cardLine;++j) // y
                  for (ordinal_type i=0;i<cardLine;++i,++idx)  // x
                    for (ordinal_type l=0;l<npts;++l)
                      output.access(idx,l,d) = output_x.access(i,l,0)*output_y.access(j,l,0)*output_z.access(k,l,0);
              ++d;
            }
          }
        break;
      }
      default: {
        INTREPID2_TEST_FOR_ABORT( true , 
                                  ">>> ERROR (Basis_HVOL_HEX_Cn_FEM): Operator type not implemented");
        break;
      }
      }
    }
  
    template<typename SpT, ordinal_type numPtsPerEval,
             typename outputValueValueType, class ...outputValueProperties,
             typename inputPointValueType,  class ...inputPointProperties,
             typename vinvValueType,        class ...vinvProperties>
    void
    Basis_HVOL_HEX_Cn_FEM::
    getValues(       Kokkos::DynRankView<outputValueValueType,outputValueProperties...> outputValues,
               const Kokkos::DynRankView<inputPointValueType, inputPointProperties...>  inputPoints,
               const Kokkos::DynRankView<vinvValueType,       vinvProperties...>        vinv,
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
      const ordinal_type cardLine = std::cbrt(cardinality);
      const ordinal_type workSize = 4*cardLine;

      auto vcprop = Kokkos::common_view_alloc_prop(inputPoints);
      typedef typename Kokkos::DynRankView< inputPointType, typename inputPointViewType::memory_space> workViewType;
      workViewType  work(Kokkos::view_alloc("Basis_HVOL_HEX_Cn_FEM::getValues::work", vcprop), workSize, inputPoints.extent(0));

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
                                      ">>> ERROR (Basis_HVOL_HEX_Cn_FEM): Operator type not implemented" );
        // break; commented out since exception is thrown
      }
      }
    }
  }

  // -------------------------------------------------------------------------------------
  template<typename SpT, typename OT, typename PT>
  Basis_HVOL_HEX_Cn_FEM<SpT,OT,PT>::
  Basis_HVOL_HEX_Cn_FEM( const ordinal_type order,
                          const EPointType   pointType ) {

    // this should be in host
    Basis_HVOL_LINE_Cn_FEM<SpT,OT,PT> lineBasis( order, pointType );
    const auto cardLine = lineBasis.getCardinality();

    this->vinv_ = Kokkos::DynRankView<typename scalarViewType::value_type,SpT>("HVOL::HEX::Cn::vinv", cardLine, cardLine);
    lineBasis.getVandermondeInverse(this->vinv_);
    
    this->basisCardinality_  = cardLine*cardLine*cardLine;
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
      ordinal_type tags[maxCardLine*maxCardLine*maxCardLine][4];

      {
        ordinal_type idx = 0;
        for (auto k=0;k<cardLine;++k) { // z
          const auto tag_z = lineBasis.getDofTag(k);
          for (ordinal_type j=0;j<cardLine;++j) { // y
            const auto tag_y = lineBasis.getDofTag(j);
            for (ordinal_type i=0;i<cardLine;++i,++idx) { // x
              const auto tag_x = lineBasis.getDofTag(i);
              
              // interior
              tags[idx][0] = 3; // interior dof
              tags[idx][1] = 0;
              tags[idx][2] = tag_x(2) + tag_x(3)*tag_y(2) + tag_x(3)*tag_y(3)*tag_z(2); // local dof id
              tags[idx][3] = tag_x(3)*tag_y(3)*tag_z(3); // total number of dofs in this vertex
            }
          }
        }
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
      dofCoordsLine("dofCoordsLine", cardLine, 1);

    lineBasis.getDofCoords(dofCoordsLine);
    auto dofCoordsLineHost = Kokkos::create_mirror_view(Kokkos::HostSpace(), dofCoordsLine);
    Kokkos::deep_copy(dofCoordsLineHost, dofCoordsLine);
    {
      ordinal_type idx = 0;
      for (auto k=0;k<cardLine;++k) { // z
        for (ordinal_type j=0;j<cardLine;++j) { // y
          for (ordinal_type i=0;i<cardLine;++i,++idx) { // x
            dofCoordsHost(idx,0) = dofCoordsLineHost(i,0);
            dofCoordsHost(idx,1) = dofCoordsLineHost(j,0);
            dofCoordsHost(idx,2) = dofCoordsLineHost(k,0);
          }
        }
      }

    }

    this->dofCoords_ = Kokkos::create_mirror_view(typename SpT::memory_space(), dofCoordsHost);
    Kokkos::deep_copy(this->dofCoords_, dofCoordsHost);
  }
  
}// namespace Intrepid2

#endif
