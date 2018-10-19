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

/** \file   Intrepid2_HGRAD_TRI_Cn_FEM.hpp
    \brief  Header file for the Intrepid2::Basis_HGRAD_TRI_Cn_FEM class.
    \author Created by R. Kirby and P. Bochev and D. Ridzal.
    Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_HGRAD_TRI_CN_FEM_HPP__
#define __INTREPID2_HGRAD_TRI_CN_FEM_HPP__

#include "Intrepid2_Basis.hpp"
#include "Intrepid2_HGRAD_TRI_Cn_FEM_ORTH.hpp"

#include "Intrepid2_PointTools.hpp"
#include "Teuchos_LAPACK.hpp"

namespace Intrepid2 {

  /** \class  Intrepid2::Basis_HGRAD_TRI_Cn_FEM
      \brief  Implementation of the default H(grad)-compatible Lagrange basis of arbitrary degree  on Triangle cell

      Implements Lagrangian basis of degree n on the reference Triangle cell. The basis has
      cardinality (n+1)(n+2)/2 and spans a COMPLETE polynomial space of degree n.
      Basis functions are dual
      to a unisolvent set of degrees-of-freedom (DoF) defined on a lattice of order n (see PointTools).
      In particular, the degrees of freedom are point evaluation at
      \li the vertices
      \li (n-1) points on each edge of the triangle
      \li max((n-1)(n-2)/2,0) points on the inside of the triangle.

      The distribution of these points is specified by the pointType argument to the class constructor.
      Currently, either equispaced lattice points or Warburton's warp-blend points are available.

      The dof are enumerated according to the ordering on the lattice (see PointTools).  In particular,
      dof number 0 is at the bottom left vertex (0,0).  The dof increase
      along the lattice with points along the lines of constant
      x adjacent in the enumeration.
  */

  namespace Impl {

    /**
       \brief See Intrepid2::Basis_HGRAD_TRI_Cn_FEM
    */
    class Basis_HGRAD_TRI_Cn_FEM {
    public:
      typedef struct Triangle<3> cell_topology_type;
      /**
         \brief See Intrepid2::Basis_HGRAD_TRI_Cn_FEM
         work is a rank 1 view having the same value_type of inputPoints
         and having size equal to getWorkSizePerPoint()*inputPoints.extent(0);
      */
      template<EOperator opType>
      struct Serial {
        template<typename outputValueViewType,
                 typename inputPointViewType,
                 typename workViewType,
                 typename vinvViewType>
        KOKKOS_INLINE_FUNCTION
        static void
        getValues(       outputValueViewType outputValues,
                         const inputPointViewType  inputPoints,
                         workViewType        work,
                         const vinvViewType        vinv );
      };

      template<typename ExecSpaceType, ordinal_type numPtsPerEval,
               typename outputValueValueType, class ...outputValueProperties,
               typename inputPointValueType,  class ...inputPointProperties,
               typename vinvValueType,        class ...vinvProperties>
      static void
      getValues(        Kokkos::DynRankView<outputValueValueType,outputValueProperties...> outputValues,
                        const Kokkos::DynRankView<inputPointValueType, inputPointProperties...>  inputPoints,
                        const Kokkos::DynRankView<vinvValueType,       vinvProperties...>        vinv,
                        const EOperator operatorType);

      /**
         \brief See Intrepid2::Basis_HGRAD_TRI_Cn_FEM
      */
      template<typename outputValueViewType,
               typename inputPointViewType,
               typename vinvViewType,
               typename workViewType,
               EOperator opType,
               ordinal_type numPtsEval>
      struct Functor {
        outputValueViewType _outputValues;
        const inputPointViewType  _inputPoints;
        const vinvViewType        _vinv;
        workViewType        _work;

        KOKKOS_INLINE_FUNCTION
        Functor(       outputValueViewType outputValues_,
                       inputPointViewType  inputPoints_,
                       vinvViewType        vinv_,
                       workViewType        work_)
          : _outputValues(outputValues_), _inputPoints(inputPoints_),
            _vinv(vinv_), _work(work_) {}

        KOKKOS_INLINE_FUNCTION
        void operator()(const size_type iter) const {
          const auto ptBegin = Util<ordinal_type>::min(iter*numPtsEval,    _inputPoints.extent(0));
          const auto ptEnd   = Util<ordinal_type>::min(ptBegin+numPtsEval, _inputPoints.extent(0));

          const auto ptRange = Kokkos::pair<ordinal_type,ordinal_type>(ptBegin, ptEnd);
          const auto input   = Kokkos::subview( _inputPoints, ptRange, Kokkos::ALL() );

          typename workViewType::pointer_type ptr = _work.data() + _work.extent(0)*ptBegin*get_dimension_scalar(_work);

          auto vcprop = Kokkos::common_view_alloc_prop(_work);
          workViewType  work(Kokkos::view_wrap(ptr,vcprop), (ptEnd-ptBegin)*_work.extent(0));

          switch (opType) {
          case OPERATOR_VALUE : {
            auto output = Kokkos::subview( _outputValues, Kokkos::ALL(), ptRange );
            Serial<opType>::getValues( output, input, work, _vinv );
            break;
          }
          case OPERATOR_CURL:
          case OPERATOR_D1:
          case OPERATOR_D2:  {
            auto output = Kokkos::subview( _outputValues, Kokkos::ALL(), ptRange, Kokkos::ALL() );
            Serial<opType>::getValues( output, input, work, _vinv );
            break;
          }
          default: {
            INTREPID2_TEST_FOR_ABORT( true,
                                      ">>> ERROR: (Intrepid2::Basis_HGRAD_TRI_Cn_FEM::Functor) operator is not supported");

          }
          }
        }
      };
    };
  }

  template<typename ExecSpaceType = void,
           typename outputValueType = double,
           typename pointValueType = double>
  class Basis_HGRAD_TRI_Cn_FEM
    : public Basis<ExecSpaceType,outputValueType,pointValueType> {
  public:
    typedef typename Basis<ExecSpaceType,outputValueType,pointValueType>::ordinal_type_array_1d_host ordinal_type_array_1d_host;
    typedef typename Basis<ExecSpaceType,outputValueType,pointValueType>::ordinal_type_array_2d_host ordinal_type_array_2d_host;
    typedef typename Basis<ExecSpaceType,outputValueType,pointValueType>::ordinal_type_array_3d_host ordinal_type_array_3d_host;

    typedef typename Basis<ExecSpaceType,outputValueType,pointValueType>::outputViewType outputViewType;
    typedef typename Basis<ExecSpaceType,outputValueType,pointValueType>::pointViewType  pointViewType;
    typedef typename Basis<ExecSpaceType,outputValueType,pointValueType>::scalarViewType  scalarViewType;
    typedef typename Basis<ExecSpaceType,outputValueType,pointValueType>::scalarType  scalarType;

  private:

    /** \brief inverse of Generalized Vandermonde matrix, whose columns store the expansion
        coefficients of the nodal basis in terms of phis_ */
    Kokkos::DynRankView<scalarType,ExecSpaceType> vinv_;

  public:
    /** \brief  Constructor.
     */
    Basis_HGRAD_TRI_Cn_FEM(const ordinal_type order,
                           const EPointType   pointType = POINTTYPE_EQUISPACED);

    using Basis<ExecSpaceType,outputValueType,pointValueType>::getValues;

    virtual
    void
    getValues(       outputViewType outputValues,
                     const pointViewType  inputPoints,
                     const EOperator operatorType = OPERATOR_VALUE) const {
#ifdef HAVE_INTREPID2_DEBUG
      Intrepid2::getValues_HGRAD_Args(outputValues,
                                      inputPoints,
                                      operatorType,
                                      this->getBaseCellTopology(),
                                      this->getCardinality() );
#endif
      constexpr ordinal_type numPtsPerEval = Parameters::MaxNumPtsPerBasisEval;
      Impl::Basis_HGRAD_TRI_Cn_FEM::
        getValues<ExecSpaceType,numPtsPerEval>( outputValues,
                                                inputPoints,
                                                this->vinv_,
                                                operatorType);
    }

    virtual
    void
    getDofCoords( scalarViewType dofCoords ) const {
#ifdef HAVE_INTREPID2_DEBUG
      // Verify rank of output array.
      INTREPID2_TEST_FOR_EXCEPTION( dofCoords.rank() != 2, std::invalid_argument,
                                    ">>> ERROR: (Intrepid2::Basis_HGRAD_TRI_Cn_FEM::getDofCoords) rank = 2 required for dofCoords array");
      // Verify 0th dimension of output array.
      INTREPID2_TEST_FOR_EXCEPTION( static_cast<ordinal_type>(dofCoords.extent(0)) != this->getCardinality(), std::invalid_argument,
                                    ">>> ERROR: (Intrepid2::Basis_HGRAD_TRI_Cn_FEM::getDofCoords) mismatch in number of dof and 0th dimension of dofCoords array");
      // Verify 1st dimension of output array.
      INTREPID2_TEST_FOR_EXCEPTION( dofCoords.extent(1) != this->getBaseCellTopology().getDimension(), std::invalid_argument,
                                    ">>> ERROR: (Intrepid2::Basis_HGRAD_TRI_Cn_FEM::getDofCoords) incorrect reference cell (1st) dimension in dofCoords array");
#endif
      Kokkos::deep_copy(dofCoords, this->dofCoords_);
    }

    virtual
    void
    getDofCoeffs( scalarViewType dofCoeffs ) const {
#ifdef HAVE_INTREPID2_DEBUG
      // Verify rank of output array.
      INTREPID2_TEST_FOR_EXCEPTION( dofCoeffs.rank() != 1, std::invalid_argument,
                                    ">>> ERROR: (Intrepid2::Basis_HGRAD_TRI_Cn_FEM::getdofCoeffs) rank = 1 required for dofCoeffs array");
      // Verify 0th dimension of output array.
      INTREPID2_TEST_FOR_EXCEPTION( static_cast<ordinal_type>(dofCoeffs.extent(0)) != this->getCardinality(), std::invalid_argument,
                                    ">>> ERROR: (Intrepid2::Basis_HGRAD_TRI_Cn_FEM::getdofCoeffs) mismatch in number of dof and 0th dimension of dofCoeffs array");
#endif
      Kokkos::deep_copy(dofCoeffs, 1.0);
    }


    virtual
    const char*
    getName() const {
      return "Intrepid2_HGRAD_TRI_Cn_FEM";
    }

    virtual
    bool
    requireOrientation() const {
      return (this->basisDegree_ > 2);
    }

    void
    getVandermondeInverse( scalarViewType vinv ) const {
      // has to be same rank and dimensions
      Kokkos::deep_copy(vinv, this->vinv_);
    }

    Kokkos::DynRankView<typename scalarViewType::const_value_type,ExecSpaceType>
    getVandermondeInverse() const {
      return vinv_;
    }

    ordinal_type
    getWorkSizePerPoint(const EOperator operatorType) const {
      auto cardinality = getPnCardinality<2>(this->basisDegree_);
      switch (operatorType) {
      case OPERATOR_GRAD:
      case OPERATOR_CURL:
      case OPERATOR_D1:
        return 5*cardinality;
      default:
        return getDkCardinality(operatorType, 2)*cardinality;
      }
    }
    
  };

}// namespace Intrepid2

#include "Intrepid2_HGRAD_TRI_Cn_FEMDef.hpp"

#endif
