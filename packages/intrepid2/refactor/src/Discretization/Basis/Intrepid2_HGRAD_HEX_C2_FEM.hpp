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

/** \file   Intrepid_HGRAD_HEX_C2_FEM.hpp
    \brief  Header file for the Intrepid2::HGRAD_HEX_C2_FEM class.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
 */

#ifndef INTREPID2_HGRAD_HEX_C2_FEM_HPP
#define INTREPID2_HGRAD_HEX_C2_FEM_HPP
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
  
  template<typename ExecSpaceType = void,
           typename outputValueType = double,
           typename pointValueType = double>
  class Basis_HGRAD_HEX_C2_FEM : public Basis<ExecSpaceType,outputValueType,pointValueType> {
  public:
    typedef typename Basis<ExecSpaceType,outputValueType,pointValueType>::ordinal_type_array_1d_host ordinal_type_array_1d_host;
    typedef typename Basis<ExecSpaceType,outputValueType,pointValueType>::ordinal_type_array_2d_host ordinal_type_array_2d_host;
    typedef typename Basis<ExecSpaceType,outputValueType,pointValueType>::ordinal_type_array_3d_host ordinal_type_array_3d_host;

    template<EOperator opType>
    struct Serial {
      template<typename outputValueValueType, class ...outputValueProperties,
               typename inputPointValueType,  class ...inputPointProperties>
      KOKKOS_INLINE_FUNCTION
      static void
      getValues( /**/  Kokkos::DynRankView<outputValueValueType,outputValueProperties...> outputValues,
                 const Kokkos::DynRankView<inputPointValueType, inputPointProperties...>  inputPoints );

    };

    template<typename outputValueViewType,
             typename inputPointViewType,
             EOperator opType>
    struct Functor {
      /**/  outputValueViewType _outputValues;
      const inputPointViewType  _inputPoints;

      KOKKOS_INLINE_FUNCTION
      Functor( /**/  outputValueViewType outputValues_,
               /**/  inputPointViewType  inputPoints_ )
        : _outputValues(outputValues_), _inputPoints(inputPoints_) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const ordinal_type pt) const {
        switch (opType) {
        case OPERATOR_VALUE : {
          auto       output = Kokkos::subdynrankview( _outputValues, Kokkos::ALL(), pt );
          const auto input  = Kokkos::subdynrankview( _inputPoints,                 pt, Kokkos::ALL() );
          Serial<opType>::getValues( output, input );
          break;
        }
        case OPERATOR_GRAD :
        case OPERATOR_D1 :
        case OPERATOR_D2 :
        case OPERATOR_D3 :
        case OPERATOR_D4 :
        case OPERATOR_MAX : {
          auto       output = Kokkos::subdynrankview( _outputValues, Kokkos::ALL(), pt, Kokkos::ALL() );
          const auto input  = Kokkos::subdynrankview( _inputPoints,                 pt, Kokkos::ALL() );
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

    class Internal {
    private:
      Basis_HGRAD_HEX_C2_FEM *obj_;

    public:
      Internal(Basis_HGRAD_HEX_C2_FEM *obj)
        : obj_(obj) {}

      /** \brief  Evaluation of a FEM basis on a <strong>reference Hexahedron</strong> cell.

          Returns values of <var>operatorType</var> acting on FEM basis functions for a set of
          points in the <strong>reference Hexahedron</strong> cell. For rank and dimensions of
          I/O array arguments see Section \ref basis_md_array_sec.

          \param  outputValues      [out] - rank-2 or 3 array with the computed basis values
          \param  inputPoints       [in]  - rank-2 array with dimensions (P,D) containing reference points
          \param  operatorType      [in]  - operator applied to basis functions
      */
      template<typename outputValueValueType, class ...outputValueProperties,
               typename inputPointValueType,  class ...inputPointProperties>
      void
      getValues( /**/  Kokkos::DynRankView<outputValueValueType,outputValueProperties...> outputValues,
                 const Kokkos::DynRankView<inputPointValueType, inputPointProperties...>  inputPoints,
                 const EOperator operatorType  = OPERATOR_VALUE ) const;


      /** \brief  Returns spatial locations (coordinates) of degrees of freedom on a
          <strong>reference Quadrilateral</strong>.

          \param  DofCoords      [out] - array with the coordinates of degrees of freedom,
          dimensioned (F,D)
      */
      template<typename dofCoordValueType, class ...dofCoordProperties>
      void
      getDofCoords( Kokkos::DynRankView<dofCoordValueType,dofCoordProperties...> dofCoords ) const;

    };
    Internal impl_;


  /** \brief  Constructor.
    */
    Basis_HGRAD_HEX_C2_FEM();
    Basis_HGRAD_HEX_C2_FEM(const Basis_HGRAD_HEX_C2_FEM &b)
    : Basis<ExecSpaceType,outputValueType,pointValueType>(b),
      impl_(this) {}

    Basis_HGRAD_HEX_C2_FEM& operator=(const Basis_HGRAD_HEX_C2_FEM &b) {
      if (this != &b) {
        Basis<ExecSpaceType,outputValueType,pointValueType>::operator= (b);
        // do not copy impl
      }
      return *this;
    }
  
    typedef typename Basis<ExecSpaceType,outputValueType,pointValueType>::outputViewType outputViewType;
    typedef typename Basis<ExecSpaceType,outputValueType,pointValueType>::pointViewType  pointViewType;

    virtual
    void
    getValues( /**/  outputViewType outputValues,
               const pointViewType  inputPoints,
               const EOperator operatorType = OPERATOR_VALUE ) const {
      impl_.getValues( outputValues, inputPoints, operatorType );
    }

    virtual
    void
    getDofCoords( pointViewType dofCoords ) const {
      impl_.getDofCoords( dofCoords );
    }
};
}// namespace Intrepid2

#include "Intrepid2_HGRAD_HEX_C2_FEMDef.hpp"

#endif
