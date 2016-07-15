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

/** \file   Intrepid_HGRAD_QUAD_Cn_FEMDef.hpp
    \brief  Definition file for the Intrepid2::HGRAD_QUAD_Cn_FEM class.
    \author Created by R. Kirby.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_HGRAD_QUAD_CN_FEM_DEF_HPP__
#define __INTREPID2_HGRAD_QUAD_CN_FEM_DEF_HPP__

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
    Basis_HGRAD_QUAD_Cn_FEM::Serial<opType>::
    getValues( /**/  outputViewType output,
               const inputViewType  input,
               /**/  workViewType   work,
               const vinvViewType   vinv,
               const ordinal_type   operatorDn ) {
      ordinal_type opDn = operatorDn;
      
      const auto card = vinv.dimension(0);
      const auto npts = input.dimension(0);

      typedef Kokkos::pair<ordinal_type,ordinal_type> range_type;
      const auto input_x = Kokkos::subdynrankview(input, Kokkos::ALL(), range_type(0,1));
      const auto input_y = Kokkos::subdynrankview(input, Kokkos::ALL(), range_type(1,2));
      
      switch (opType) {
      case OPERATOR_VALUE: {
        auto ptr = work.data();
        
        Kokkos::DynRankView<typename workViewType::value_type,
            typename workViewType::memory_space> work_line(ptr, card, npts);
        ptr += (card*npts);

        Kokkos::DynRankView<typename workViewType::value_type,
            typename workViewType::memory_space> output_x(ptr, card, npts);
        ptr += (card*npts);

        Kokkos::DynRankView<typename workViewType::value_type,
            typename workViewType::memory_space> output_y(ptr, card, npts);
        ptr += (card*npts);
        
        Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
          getValues(output_x, input_x, work_line, vinv);

        Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
          getValues(output_y, input_y, work_line, vinv);

        // tensor product
        ordinal_type idx = 0;
        for (auto j=0;j<card;++j) // y
          for (auto i=0;i<card;++i,++idx)  // x
            for (auto k=0;k<npts;++k) 
              output(idx,k) = output_x(i,k)*output_y(j,k);
        break;
      }
      case OPERATOR_CURL: {
        for (auto l=0;l<2;++l) {
          auto ptr = work.data();
          
          Kokkos::DynRankView<typename workViewType::value_type,
            typename workViewType::memory_space> work_line(ptr, card, npts);
          ptr += (card*npts);
          
          Kokkos::DynRankView<typename workViewType::value_type,
            typename workViewType::memory_space,Kokkos::MemoryUnmanaged> output_x, output_y;
          
          typename workViewType::value_type s = 0.0;
          if (l) {
            // l = 1
            output_x = Kokkos::DynRankView<typename workViewType::value_type,
              typename workViewType::memory_space>(ptr, card, npts, 1);
            ptr += (card*npts);
            Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_Dn>::
              getValues(output_x, input_x, work_line, vinv, 1);                           

            output_y = Kokkos::DynRankView<typename workViewType::value_type,
              typename workViewType::memory_space>(ptr, card, npts);
            ptr += (card*npts);
            Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
              getValues(output_y, input_y, work_line, vinv);                           

            s = -1.0;
          } else {
            // l = 0
            output_x = Kokkos::DynRankView<typename workViewType::value_type,
              typename workViewType::memory_space>(ptr, card, npts);
            ptr += (card*npts);
            Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
              getValues(output_x, input_x, work_line, vinv);                           

            output_y = Kokkos::DynRankView<typename workViewType::value_type,
              typename workViewType::memory_space>(ptr, card, npts, 1);
            ptr += (card*npts);
            Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_Dn>::
              getValues(output_y, input_y, work_line, vinv, 1);                           

            s = 1.0;
          }

          // tensor product (extra dimension of ouput x and y are ignored)
          ordinal_type idx = 0;
          for (auto j=0;j<card;++j) // y
            for (auto i=0;i<card;++i,++idx)  // x
              for (auto k=0;k<npts;++k) 
                output(idx,k,l) = s*output_x(i,k,0)*output_y(j,k,0);
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
        opDn = getOperatorOrder(opType);
      case OPERATOR_Dn: {
        const auto dkcard = opDn + 1;
        for (auto l=0;l<dkcard;++l) {
          auto ptr = work.data();
          
          Kokkos::DynRankView<typename workViewType::value_type,
            typename workViewType::memory_space> work_line(ptr, card, npts);
          ptr += (card*npts);
          
          Kokkos::DynRankView<typename workViewType::value_type,
            typename workViewType::memory_space,Kokkos::MemoryUnmanaged> output_x, output_y;
          
          const auto mult_x = opDn - l;
          const auto mult_y = l;
          
          if (mult_x) {
            output_x = Kokkos::DynRankView<typename workViewType::value_type,
              typename workViewType::memory_space>(ptr, card, npts, 1);
            ptr += (card*npts);
            Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_Dn>::
              getValues(output_x, input_x, work_line, vinv, mult_x);                           
          } else {
            output_x = Kokkos::DynRankView<typename workViewType::value_type,
              typename workViewType::memory_space>(ptr, card, npts);
            ptr += (card*npts);
            Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
              getValues(output_x, input_x, work_line, vinv);                           
          }

          if (mult_y) {
            output_y = Kokkos::DynRankView<typename workViewType::value_type,
              typename workViewType::memory_space>(ptr, card, npts, 1);
            ptr += (card*npts);
            Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_Dn>::
              getValues(output_y, input_y, work_line, vinv, mult_y);                           
          } else {
            output_y = Kokkos::DynRankView<typename workViewType::value_type,
              typename workViewType::memory_space>(ptr, card, npts);
            ptr += (card*npts);
            Impl::Basis_HGRAD_LINE_Cn_FEM::Serial<OPERATOR_VALUE>::
              getValues(output_y, input_y, work_line, vinv);                           
          }

          // tensor product (extra dimension of ouput x and y are ignored)
          ordinal_type idx = 0;
          for (auto j=0;j<card;++j) // y
            for (auto i=0;i<card;++i,++idx)  // x
              for (auto k=0;k<npts;++k) 
                output(idx,k,l) = output_x(i,k,0)*output_y(j,k,0);
        }
        break;
      }
      default: {
        INTREPID2_TEST_FOR_ABORT( true,
                                  ">>> ERROR: (Intrepid2::Basis_HGRAD_QUAD_Cn_FEM::Serial::getValues) operator is not supported" );
      }
      }
    }
    
    template<typename SpT, ordinal_type numPtsPerEval,
             typename outputValueValueType, class ...outputValueProperties,
             typename inputPointValueType,  class ...inputPointProperties,
             typename vinvValueType,        class ...vinvProperties>
    void
    Basis_HGRAD_QUAD_Cn_FEM::
    getValues( /**/  Kokkos::DynRankView<outputValueValueType,outputValueProperties...> outputValues,
               const Kokkos::DynRankView<inputPointValueType, inputPointProperties...>  inputPoints,
               const Kokkos::DynRankView<vinvValueType,       vinvProperties...>        vinv,
               const EOperator operatorType ) {
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
        Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints, vinv) );
        break;
      }
      case OPERATOR_CURL: {
        typedef Functor<outputValueViewType,inputPointViewType,vinvViewType,
            OPERATOR_CURL,numPtsPerEval> FunctorType;
        Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints, vinv) );
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
        typedef Functor<outputValueViewType,inputPointViewType,vinvViewType,
            OPERATOR_Dn,numPtsPerEval> FunctorType;
        Kokkos::parallel_for( policy, FunctorType(outputValues, inputPoints, vinv,
                                                  getOperatorOrder(operatorType)) );
        break;
      }
      default: {
        INTREPID2_TEST_FOR_EXCEPTION( true , std::invalid_argument,
                                      ">>> ERROR (Basis_HGRAD_QUAD_Cn_FEM): Operator type not implemented" );
        break;
      }
      }
    }
  }

  // -------------------------------------------------------------------------------------
  template<typename SpT, typename OT, typename PT>
  Basis_HGRAD_QUAD_Cn_FEM<SpT,OT,PT>::
  Basis_HGRAD_QUAD_Cn_FEM( const ordinal_type order,
                           const EPointType   pointType ) {
    // this should be in host
    Basis_HGRAD_LINE_Cn_FEM<SpT,OT,PT> lineBasis( order, pointType );
    const auto cardLine = lineBasis.getCardinality();
    
    this->vinv_ = Kokkos::DynRankView<OT,SpT>("Hgrad::Quad::Cn::vinv", cardLine, cardLine);         
    lineBasis.getVandermondeInverse(this->vinv_);

    this->basisCardinality_  = cardLine*cardLine;
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

      const ordinal_type vert[2][2] = { {0,1}, {3,2} }; //[y][x]
      {
        ordinal_type idx = 0;
        for (auto j=0;j<cardLine;++j) { // y      
          const auto tag_y = lineBasis.getDofTag(j);
          for (auto i=0;i<cardLine;++i,++idx) { // x
            const auto tag_x = lineBasis.getDofTag(i);          
            
            if (tag_x(0) == 0 && tag_y(0) == 0) {
              // vertices
              tags[idx][0] = 0; // vertex dof
              tags[idx][1] = vert[tag_y(1)][tag_x(1)]; // vertex id
              tags[idx][2] = 0; // local dof id
              tags[idx][3] = 1; // total number of dofs in this vertex
            } else if (tag_x(0) == 1 && tag_y(0) == 0) {
              // horizontal edge 
              tags[idx][0] = 1; // edge dof
              tags[idx][1] = (tag_y(1) == 0 ? 1 : 3);
              tags[idx][2] = tag_x(2); // local dof id
              tags[idx][3] = tag_x(3); // total number of dofs in this vertex
            } else if (tag_x(0) == 0 && tag_y(0) == 1) {
              // vertical edge 
              tags[idx][0] = 1; // edge dof
              tags[idx][1] = (tag_x(1) == 0 ? 4 : 2);
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
    Kokkos::DynRankView<PT,typename SpT::array_layout,Kokkos::HostSpace>
      dofCoordsHost("dofCoordsHost", this->basisCardinality_, this->basisCellTopology_.getDimension());

    Kokkos::DynRankView<PT,SpT>
      dofCoordsLine("dofCoordsLine", cardLine, 1);

    lineBasis.getDofCoords(dofCoordsLine);
    auto dofCoordsLineHost = Kokkos::create_mirror_view(Kokkos::HostSpace(), dofCoordsLine);
    Kokkos::deep_copy(dofCoordsLineHost, dofCoordsLine);
    {
      ordinal_type idx = 0;
      for (auto j=0;j<cardLine;++j) { // y      
        for (auto i=0;i<cardLine;++i,++idx) { // x
          dofCoordsHost(idx,0) = dofCoordsLineHost(i,0);
          dofCoordsHost(idx,1) = dofCoordsLineHost(j,0);
        }
      }
    }

    this->dofCoords_ = Kokkos::create_mirror_view(typename SpT::memory_space(), dofCoordsHost);
    Kokkos::deep_copy(this->dofCoords_, dofCoordsHost);
  }
  
}

#endif
