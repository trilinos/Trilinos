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

/** \file   Intrepid_HGRAD_QUAD_Cn_FEM.hpp
    \brief  Header file for the Intrepid2::HGRAD_QUAD_Cn_FEM class.
    \author Created by R. Kirby.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_HGRAD_QUAD_CN_FEM_HPP__
#define __INTREPID2_HGRAD_QUAD_CN_FEM_HPP__

#include "Intrepid2_Basis.hpp"
#include "Intrepid2_HGRAD_LINE_Cn_FEM.hpp"

namespace Intrepid2 {
  
  /** \class  Intrepid2::Basis_HGRAD_QUAD_C1_FEM
      \brief  Implementation of the default H(grad)-compatible FEM basis of degree 1 on Quadrilateral cell
      Implements Lagrangian basis of degree n on the reference Quadrilateral cell using
      a tensor product of points
  */
  template<typename ExecSpaceType = void,
           typename outputValueType = double,
           typename pointValueType = double>
  class Basis_HGRAD_QUAD_Cn_FEM
    : public Basis<ExecSpaceType,outputValueType,pointValueType> {
  public:
    typedef typename Basis<ExecSpaceType,outputValueType,pointValueType>::ordinal_type_array_1d_host ordinal_type_array_1d_host;
    typedef typename Basis<ExecSpaceType,outputValueType,pointValueType>::ordinal_type_array_2d_host ordinal_type_array_2d_host;
    typedef typename Basis<ExecSpaceType,outputValueType,pointValueType>::ordinal_type_array_3d_host ordinal_type_array_3d_host;

    // // this is specialized for each operator type
    // template<EOperator opType>
    // struct Serial {
    //   template<typename outputValueViewType,
    //            typename inputPointViewType,
    //            typename phisValueViewType>
    //   KOKKOS_INLINE_FUNCTION
    //   static void
    //   getValues( /**/  outputValueViewType outputValues,
    //              const inputPointViewType  inputPoints,
    //              const phisValueViewType   phisValues );
    // };

    // template<typename outputValueViewType,
    //          typename inputPointViewType,
    //          typename phisValueViewType,
    //          EOperator opType>
    // struct Functor {
    //   /**/  outputValueViewType _outputValues;
    //   const inputPointViewType  _inputPoints;
    //   const phisValueViewType   _phisValues;

    //   KOKKOS_INLINE_FUNCTION
    //   Functor( /**/  outputValueViewType outputValues_,
    //            /**/  inputPointViewType  inputPoints_,
    //            /**/  phisValueViewType   phisValues_ )
    //     : _outputValues(outputValues_), _inputPoints(inputPoints_), _phisValues(phisValues_) {}

    //   KOKKOS_INLINE_FUNCTION
    //   void operator()(const ordinal_type ord) const {
    //     switch (opType) {
    //     case OPERATOR_VALUE : {
    //       auto output = Kokkos::subdynrankview( _outputValues, ord, Kokkos::ALL() );
    //       Serial<opType>::getValues( output, _inputPoints, _phisValues );
    //       break;
    //     }
    //     case OPERATOR_MAX : {
    //       auto output = Kokkos::subdynrankview( _outputValues, ord, Kokkos::ALL(), Kokkos::ALL() );
    //       Serial<opType>::getValues( output, _inputPoints, _phisValues );
    //       break;
    //     }
    //     default: {
    //       INTREPID2_TEST_FOR_ABORT( true,
    //                                 ">>> ERROR: (Intrepid2::Basis_HGRAD_QUAD_Cn_FEM::Functor) operator is not supported");

    //     }
    //     }
    //   }
    // };

    class Internal {
    private:
      Basis_HGRAD_QUAD_Cn_FEM *obj_;

    public:
      Internal(Basis_HGRAD_QUAD_Cn_FEM *obj)
        : obj_(obj) {}

      /** \brief  Evaluation of a FEM basis on a <strong>reference Line</strong> cell.

          Returns values of <var>operatorType</var> acting on FEM basis functions for a set of
          points in the <strong>reference Line</strong> cell. For rank and dimensions of
          I/O array arguments see Section \ref basis_md_array_sec .

          \param  outputValues      [out] - variable rank array with the basis values
          \param  inputPoints       [in]  - rank-2 array (P,D) with the evaluation points
          \param  operatorType      [in]  - the operator acting on the basis functions
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
    Basis_HGRAD_QUAD_Cn_FEM(const ordinal_type order,
                            const EPointType   pointType);
    Basis_HGRAD_QUAD_Cn_FEM(const Basis_HGRAD_QUAD_Cn_FEM &b)
      : Basis<ExecSpaceType,outputValueType,pointValueType>(b),
        impl_(this) {}

    Basis_HGRAD_QUAD_Cn_FEM& operator=(const Basis_HGRAD_QUAD_Cn_FEM &b) {
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

  private:

    /** \brief orthogonal basis */
    Basis_HGRAD_QUAD_Cn_FEM<ExecSpaceType,outputValueType,pointValueType> lineX_, lineY;
  };

}// namespace Intrepid2

#include "Intrepid2_HGRAD_QUAD_Cn_FEMDef.hpp"

#endif
