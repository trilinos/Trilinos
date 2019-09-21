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

/** \file   Intrepid2_HGRAD_HEX_C2_FEM.hpp
    \brief  Header file for the Intrepid2::Basis_HGRAD_HEX_C2_FEM class.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
 */

#ifndef __INTREPID2_HGRAD_HEX_C2_FEM_HPP__
#define __INTREPID2_HGRAD_HEX_C2_FEM_HPP__

#include "Intrepid2_Basis.hpp"

namespace Intrepid2 {

  /** \class  Intrepid2::Basis_HGRAD_HEX_C2_FEM
      \brief  Implementation of the default H(grad)-compatible FEM basis of degree 2 on Hexahedron cell

      Implements Lagrangian basis of degree 2 on the reference Hexahedron cell. The basis has
      cardinality 27 and spans a COMPLETE tri-quadratic polynomial space. Basis functions are dual
      to a unisolvent set of degrees-of-freedom (DoF) defined and enumerated as follows:

      \verbatim
      =================================================================================================
      |         |           degree-of-freedom-tag table                    |                           |
      |   DoF   |----------------------------------------------------------|      DoF definition       |
      | ordinal |  subc dim    | subc ordinal | subc DoF ord |subc num DoF |                           |
      |=========|==============|==============|==============|=============|===========================|
      |    0    |       0      |       0      |       0      |      1      |   L_0(u) = u(-1,-1,-1)    |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |    1    |       0      |       1      |       0      |      1      |   L_1(u) = u( 1,-1,-1)    |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |    2    |       0      |       2      |       0      |      1      |   L_2(u) = u( 1, 1,-1)    |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |    3    |       0      |       3      |       0      |      1      |   L_3(u) = u(-1, 1,-1)    |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |    4    |       0      |       4      |       0      |      1      |   L_4(u) = u(-1,-1, 1)    |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |    5    |       0      |       5      |       0      |      1      |   L_5(u) = u( 1,-1, 1)    |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |    6    |       0      |       6      |       0      |      1      |   L_6(u) = u( 1, 1, 1)    |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |    7    |       0      |       7      |       0      |      1      |   L_7(u) = u(-1, 1, 1)    |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |    8    |       1      |       0      |       0      |      1      |   L_8(u) = u( 0,-1,-1)    |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |    9    |       1      |       1      |       0      |      1      |   L_9(u) = u( 1, 0,-1)    |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |   10    |       1      |       2      |       0      |      1      |   L_10(u) = u( 0, 1,-1)   |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |   11    |       1      |       3      |       0      |      1      |   L_11(u) = u(-1, 0,-1)   |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |   12    |       1      |       8      |       0      |      1      |   L_12(u) = u(-1,-1, 0)   |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |   13    |       1      |       9      |       0      |      1      |   L_13(u) = u( 1,-1, 0)   |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |   14    |       1      |      10      |       0      |      1      |   L_14(u) = u( 1, 1, 0)   |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |   15    |       1      |      11      |       0      |      1      |   L_15(u) = u(-1, 1, 0)   |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |   16    |       1      |       4      |       0      |      1      |   L_16(u) = u( 0,-1, 1)   |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |   17    |       1      |       5      |       0      |      1      |   L_17(u) = u( 1, 0, 1)   |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |   18    |       1      |       6      |       0      |      1      |   L_18(u) = u( 0, 1, 1)   |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |   19    |       1      |       7      |       0      |      1      |   L_19(u) = u(-1, 0, 1)   |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |   20    |       3      |       0      |       0      |      1      |   L_20(u) = u( 0, 0, 0)   |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |   21    |       2      |       4      |       0      |      1      |   L_21(u) = u( 0, 0,-1)   |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |   22    |       2      |       5      |       0      |      1      |   L_22(u) = u( 0, 0, 1)   |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |   23    |       2      |       3      |       0      |      1      |   L_23(u) = u(-1, 0, 0)   |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |   24    |       2      |       1      |       0      |      1      |   L_24(u) = u( 1, 0, 0)   |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |   25    |       2      |       0      |       0      |      1      |   L_25(u) = u( 0,-1, 0)   |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |   26    |       2      |       2      |       0      |      1      |   L_26(u) = u( 0, 1, 0)   |
      |=========|==============|==============|==============|=============|===========================|
      |   MAX   |  maxScDim=2  |  maxScOrd=12 |  maxDfOrd=0  |      -      |                           |
      |=========|==============|==============|==============|=============|===========================|
      \endverbatim

      \remark   Ordering of DoFs follows the node order in Hexahedron<27> topology. Note that node
      order in this topology does not follow the natural oder of k-subcells where the nodes
      are located, except for nodes 0 to 7 which coincide with the vertices of the base
      Hexahedrn <8> topology. As a result, L_0 to L_7 are associated with nodes 0 to 7, but
      L_8 to L_19 are not associated with edges 0 to 12 in that order.

  */

  namespace Impl {

    /**
      \brief See Intrepid2::Basis_HGRAD_HEX_C2_FEM
    */
    class Basis_HGRAD_HEX_C2_FEM {
    public:
      typedef struct Hexahedron<27> cell_topology_type;
      /**
        \brief See Intrepid2::Basis_HGRAD_HEX_C2_FEM
      */
      template<EOperator opType>
      struct Serial {
        template<typename OutputViewType,
                 typename inputViewType>
        KOKKOS_INLINE_FUNCTION
        static void
        getValues(       OutputViewType output,
                   const inputViewType input );
        
      };
      
      template<typename ExecSpaceType, 
               typename outputValueValueType, class ...outputValueProperties,
               typename inputPointValueType,  class ...inputPointProperties>
      static void
      getValues(       Kokkos::DynRankView<outputValueValueType,outputValueProperties...> outputValues,
                 const Kokkos::DynRankView<inputPointValueType, inputPointProperties...>  inputPoints,
                 const EOperator operatorType);
      
      /**
        \brief See Intrepid2::Basis_HGRAD_HEX_C2_FEM
      */
      template<typename outputValueViewType,
               typename inputPointViewType,
               EOperator opType>
      struct Functor {
              outputValueViewType _outputValues;
        const inputPointViewType  _inputPoints;
        
        KOKKOS_INLINE_FUNCTION
        Functor(       outputValueViewType outputValues_,
                       inputPointViewType  inputPoints_ )
          : _outputValues(outputValues_), _inputPoints(inputPoints_) {}
        
        KOKKOS_INLINE_FUNCTION
        void operator()(const ordinal_type pt) const {
          switch (opType) {
          case OPERATOR_VALUE : {
            auto       output = Kokkos::subview( _outputValues, Kokkos::ALL(), pt );
            const auto input  = Kokkos::subview( _inputPoints,                 pt, Kokkos::ALL() );
            Serial<opType>::getValues( output, input );
            break;
          }
          case OPERATOR_GRAD :
          case OPERATOR_D1 :
          case OPERATOR_D2 :
          case OPERATOR_D3 :
          case OPERATOR_D4 :
          case OPERATOR_MAX : {
            auto       output = Kokkos::subview( _outputValues, Kokkos::ALL(), pt, Kokkos::ALL() );
            const auto input  = Kokkos::subview( _inputPoints,                 pt, Kokkos::ALL() );
            Serial<opType>::getValues( output, input );
            break;
          }
          default: {
            INTREPID2_TEST_FOR_ABORT( opType != OPERATOR_VALUE &&
                                      opType != OPERATOR_GRAD &&
                                      opType != OPERATOR_D1 &&
                                      opType != OPERATOR_D2 &&
                                      opType != OPERATOR_D3 &&
                                      opType != OPERATOR_D4 &&
                                      opType != OPERATOR_MAX,
                                      ">>> ERROR: (Intrepid2::Basis_HGRAD_HEX_C2_FEM::Serial::getValues) operator is not supported");
          }
          }
        }
      };
    };
  }
      
  template<typename ExecSpaceType = void,
           typename outputValueType = double,
           typename pointValueType = double>
  class Basis_HGRAD_HEX_C2_FEM : public Basis<ExecSpaceType,outputValueType,pointValueType> {
  public:
    using OrdinalTypeArray1DHost = typename Basis<ExecSpaceType,outputValueType,pointValueType>::OrdinalTypeArray1DHost;
    using OrdinalTypeArray2DHost = typename Basis<ExecSpaceType,outputValueType,pointValueType>::OrdinalTypeArray2DHost;
    using OrdinalTypeArray3DHost = typename Basis<ExecSpaceType,outputValueType,pointValueType>::OrdinalTypeArray3DHost;

    using ordinal_type_array_1d_host INTREPID2_DEPRECATED_TYPENAME_REPLACEMENT("use OrdinalTypeArray1DHost instead","OrdinalTypeArray1DHost") = OrdinalTypeArray1DHost INTREPID2_DEPRECATED_TYPENAME_TRAILING_ATTRIBUTE("use OrdinalTypeArray1DHost instead");
    using ordinal_type_array_2d_host INTREPID2_DEPRECATED_TYPENAME_REPLACEMENT("use OrdinalTypeArray2DHost instead","OrdinalTypeArray2DHost") = OrdinalTypeArray2DHost INTREPID2_DEPRECATED_TYPENAME_TRAILING_ATTRIBUTE("use OrdinalTypeArray2DHost instead");
    using ordinal_type_array_3d_host INTREPID2_DEPRECATED_TYPENAME_REPLACEMENT("use OrdinalTypeArray3DHost instead","OrdinalTypeArray3DHost") = OrdinalTypeArray3DHost INTREPID2_DEPRECATED_TYPENAME_TRAILING_ATTRIBUTE("use OrdinalTypeArray3DHost instead");

    /** \brief  Constructor.
     */
    Basis_HGRAD_HEX_C2_FEM();

    using OutputViewType = typename Basis<ExecSpaceType,outputValueType,pointValueType>::OutputViewType;
    using PointViewType  = typename Basis<ExecSpaceType,outputValueType,pointValueType>::PointViewType;
    using ScalarViewType = typename Basis<ExecSpaceType,outputValueType,pointValueType>::ScalarViewType;
    
    using outputViewType INTREPID2_DEPRECATED_TYPENAME_REPLACEMENT("use OutputViewType instead","OutputViewType") = OutputViewType INTREPID2_DEPRECATED_TYPENAME_TRAILING_ATTRIBUTE("use OutputViewType instead");
    using pointViewType INTREPID2_DEPRECATED_TYPENAME_REPLACEMENT("use PointViewType instead","PointViewType") = PointViewType INTREPID2_DEPRECATED_TYPENAME_TRAILING_ATTRIBUTE("use PointViewType instead");
    using scalarViewType INTREPID2_DEPRECATED_TYPENAME_REPLACEMENT("use ScalarViewType instead","ScalarViewType") = ScalarViewType INTREPID2_DEPRECATED_TYPENAME_TRAILING_ATTRIBUTE("use ScalarViewType instead");

    using Basis<ExecSpaceType,outputValueType,pointValueType>::getValues;

    virtual
    void
    getValues(       OutputViewType outputValues,
               const PointViewType  inputPoints,
               const EOperator operatorType = OPERATOR_VALUE ) const {
#ifdef HAVE_INTREPID2_DEBUG
      // Verify arguments
      Intrepid2::getValues_HGRAD_Args(outputValues,
                                      inputPoints,
                                      operatorType,
                                      this->getBaseCellTopology(),
                                      this->getCardinality() );
#endif
      Impl::Basis_HGRAD_HEX_C2_FEM::
        getValues<ExecSpaceType>( outputValues,
                                  inputPoints,
                                  operatorType );
    }

    virtual
    void
    getDofCoords( ScalarViewType dofCoords ) const {
#ifdef HAVE_INTREPID2_DEBUG
      // Verify rank of output array.
      INTREPID2_TEST_FOR_EXCEPTION( dofCoords.rank() != 2, std::invalid_argument,
                                    ">>> ERROR: (Intrepid2::Basis_HGRAD_HEX_C2_FEM::getDofCoords) rank = 2 required for dofCoords array");
      // Verify 0th dimension of output array.
      INTREPID2_TEST_FOR_EXCEPTION( static_cast<ordinal_type>(dofCoords.extent(0)) != this->getCardinality(), std::invalid_argument,
                                    ">>> ERROR: (Intrepid2::Basis_HGRAD_HEX_C2_FEM::getDofCoords) mismatch in number of dof and 0th dimension of dofCoords array");
      // Verify 1st dimension of output array.
      INTREPID2_TEST_FOR_EXCEPTION( dofCoords.extent(1) != this->getBaseCellTopology().getDimension(), std::invalid_argument,
                                    ">>> ERROR: (Intrepid2::Basis_HGRAD_HEX_C2_FEM::getDofCoords) incorrect reference cell (1st) dimension in dofCoords array");
#endif
      Kokkos::deep_copy(dofCoords, this->dofCoords_);
    }

    virtual
    void
    getDofCoeffs( ScalarViewType dofCoeffs ) const {
#ifdef HAVE_INTREPID2_DEBUG
      // Verify rank of output array.
      INTREPID2_TEST_FOR_EXCEPTION( dofCoeffs.rank() != 1, std::invalid_argument,
                                    ">>> ERROR: (Intrepid2::Basis_HGRAD_HEX_C2_FEM::getdofCoeffs) rank = 1 required for dofCoeffs array");
      // Verify 0th dimension of output array.
      INTREPID2_TEST_FOR_EXCEPTION( static_cast<ordinal_type>(dofCoeffs.extent(0)) != this->getCardinality(), std::invalid_argument,
                                    ">>> ERROR: (Intrepid2::Basis_HGRAD_HEX_C2_FEM::getdofCoeffs) mismatch in number of dof and 0th dimension of dofCoeffs array");
#endif
      Kokkos::deep_copy(dofCoeffs, 1.0);
    }

    virtual
    const char*
    getName() const {
      return "Intrepid2_HGRAD_HEX_C2_FEM";
    }

  };
}// namespace Intrepid2

#include "Intrepid2_HGRAD_HEX_C2_FEMDef.hpp"

#endif
