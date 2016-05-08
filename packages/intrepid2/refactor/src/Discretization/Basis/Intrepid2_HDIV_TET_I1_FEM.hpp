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

/** \file   Intrepid_HDIV_TET_I1_FEM.hpp
    \brief  Header file for the Intrepid2::HDIV_TET_I1_FEM class.
    \author Created by P. Bochev and D. Ridzal and K. Peterson.
 */

#ifndef INTREPID2_HDIV_TET_I1_FEM_HPP
#define INTREPID2_HDIV_TET_I1_FEM_HPP
#include "Intrepid2_Basis.hpp"

namespace Intrepid2 {
  
/** \class  Intrepid2::Basis_HDIV_TET_I1_FEM
    \brief  Implementation of the default H(div)-compatible FEM basis of degree 1 on a Tetrahedron cell.
  
            Implements Raviart-Thomas basis of degree 1 on the reference Tetrahedron cell. The basis has
            cardinality 4 and spans an INCOMPLETE linear polynomial space. Basis functions are dual
            to a unisolvent set of degrees-of-freedom (DoF) defined and enumerated as follows:
  
  \verbatim
  =========================================================================================================
  |         |           degree-of-freedom-tag table                    |                                  |
  |   DoF   |----------------------------------------------------------|       DoF definition             |
  | ordinal |  subc dim    | subc ordinal | subc DoF ord |subc num DoF |                                  |
  |=========|==============|==============|==============|=============|==================================|
  |    0    |       2      |       0      |       0      |      1      | L_0(u) = (u.n)(1/3,0,1/3)        |
  |---------|--------------|--------------|--------------|-------------|----------------------------------|
  |    1    |       2      |       1      |       0      |      1      | L_1(u) = (u.n)(1/3,1/3,1/3)      |
  |---------|--------------|--------------|--------------|-------------|----------------------------------|
  |    2    |       2      |       2      |       0      |      1      | L_2(u) = (u.n)(0,1/3,1/3)        |
  |---------|--------------|--------------|--------------|-------------|----------------------------------|
  |    3    |       2      |       3      |       0      |      1      | L_3(u) = (u.n)(1/3,1/3,0)        |
  |=========|==============|==============|==============|=============|==================================|
  |   MAX   |  maxScDim=2  |  maxScOrd=3  |  maxDfOrd=0  |      -      |                                  |
  |=========|==============|==============|==============|=============|==================================|
  \endverbatim
  
    \remarks
    \li     In the DoF functional \f${\bf n}\f$ is a face normal. Direction of face normals 
            is determined by the right-hand rule applied to faces oriented by their vertex order
            in the cell topology, from face vertex 0 to last face vertex, whereas their length is
            set equal to face area (see http://mathworld.wolfram.com/Right-HandRule.html for definition 
            of right-hand rule). For example, face 1 of all Tetrahedron cells has vertex order {1,2,3} 
            and its right-hand rule normal can be computed, e.g., by the vector product of edge 
            tangents to edges {1,2} and {2,3}. On the reference Tetrahedron the coordinates of 
            face 1 vertices are (1,0,0), (0,1,0), and (0,0,1), the edge tangents are (-1,1,0) and 
            (0,-1,1) and the face normal direction is (-1,1,0) X (0,-1,1) = (1,1,1). Length of this 
            raw face normal is twice the face area of face 1 and so the final face normal to face 1 is
            obtained by scaling the raw normal by 1/2: (1/2,1/2,1/2).
  
    \li     The length of the face normal equals the face area. As a result, the DoF functional 
            is the value of the normal component of a vector field at the face center times the 
            face area. The resulting basis is equivalent to a basis defined by using the face 
            flux as a DoF functional. Note that faces 0, 2, and 3 of reference Tetrahedron<> 
            cells have area 1/2 and face 1 has area Sqrt(3)/2.
  
 */
  
  template<typename ExecSpaceType = void,
           typename outputValueType = double,
           typename pointValueType = double>
  class Basis_HDIV_TET_I1_FEM: public Basis<ExecSpaceType,outputValueType,pointValueType> {
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
          auto       output = Kokkos::subdynrankview( _outputValues, Kokkos::ALL(), pt, Kokkos::ALL() );
          const auto input  = Kokkos::subdynrankview( _inputPoints,                 pt, Kokkos::ALL() );
          Serial<opType>::getValues( output, input );
          break;
        }
        case OPERATOR_DIV : {
          auto       output = Kokkos::subdynrankview( _outputValues, Kokkos::ALL(), pt );
          const auto input  = Kokkos::subdynrankview( _inputPoints,                 pt, Kokkos::ALL() );
          Serial<opType>::getValues( output, input );
          break;
        }
        default: {
          INTREPID2_TEST_FOR_ABORT( opType != OPERATOR_VALUE &&
                                    opType != OPERATOR_DIV,
                                    ">>> ERROR: (Intrepid2::Basis_HDIV_TET_I1_FEM::Serial::getValues) operator is not supported");
        }
        }
      }
    };
  
    class Internal {
    private:
      Basis_HDIV_TET_I1_FEM *obj_;

    public:
      Internal(Basis_HDIV_TET_I1_FEM *obj)
        : obj_(obj) {}

      /** \brief  Evaluation of a FEM basis on a <strong>reference Tetrahedron</strong> cell.

          Returns values of <var>operatorType</var> acting on FEM basis functions for a set of
          points in the <strong>reference Tetrahedron</strong> cell. For rank and dimensions of
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
          <strong>reference Triangule</strong>.

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
    Basis_HDIV_TET_I1_FEM();
    Basis_HDIV_TET_I1_FEM(const Basis_HDIV_TET_I1_FEM &b)
      : Basis<ExecSpaceType,outputValueType,pointValueType>(b),
        impl_(this) {}

    Basis_HDIV_TET_I1_FEM& operator=(const Basis_HDIV_TET_I1_FEM &b) {
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
    
  /** \brief  Returns spatial locations (coordinates) of degrees of freedom on a
              <strong>reference Tetrahedron</strong>.

      \param  DofCoords      [out] - array with the coordinates of degrees of freedom,
                                     dimensioned (F,D)
  */

    virtual
    void
    getDofCoords( pointViewType dofCoords ) const {
      impl_.getDofCoords( dofCoords );
    }

  };

}// namespace Intrepid2

#include "Intrepid2_HDIV_TET_I1_FEMDef.hpp"

#endif
