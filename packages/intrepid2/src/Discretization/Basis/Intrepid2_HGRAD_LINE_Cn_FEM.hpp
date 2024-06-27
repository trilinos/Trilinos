// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_HGRAD_LINE_Cn_FEM.hpp
    \brief  Header file for the Intrepid2::Basis_HGRAD_LINE_Cn_FEM class.
    \author Created by R. Kirby and P. Bochev and D. Ridzal.
    Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_HGRAD_LINE_CN_FEM_HPP__
#define __INTREPID2_HGRAD_LINE_CN_FEM_HPP__

#include "Intrepid2_Basis.hpp"
#include "Intrepid2_HGRAD_LINE_Cn_FEM_JACOBI.hpp"

#include "Intrepid2_PointTools.hpp"
#include "Teuchos_LAPACK.hpp"

namespace Intrepid2 {

  /** \class  Intrepid2::Basis_HGRAD_LINE_Cn_FEM
      \brief  Implementation of the locally H(grad)-compatible FEM basis of variable order
      on the [-1,1] reference line cell, using Lagrange polynomials. 

      Implements Lagrange basis of variable order \f$n\f$ on
      the reference [-1,1] line cell.  The distribution of the points
      may be equispaced points with our without the endpoints, the Gauss-Legendre
      or Gauss-Lobatto points.  These points are {x_i}_{i=0}^n. with x_i < x_{i+1}

      The basis has cardinality \f$n+1\f$ and spans a COMPLETE linear polynomial space.
      Basis functions are dual to a unisolvent set of degrees of freedom (DoF)
      n_i( psi ) = psi(x_i).  The DoF are ordered by i.  The DoF at points 
      -1 and 1 (if included in {x_i} are attached to the vertices, and the rest of
      the DoF are attached to the edge itself.
  */

  namespace Impl {

    /**
       \brief See Intrepid2::Basis_HGRAD_LINE_Cn_FEM
    */
    class Basis_HGRAD_LINE_Cn_FEM {
    public:
      typedef struct Line<2> cell_topology_type;
      /**
         \brief See Intrepid2::Basis_HGRAD_LINE_Cn_FEM
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
                         const vinvViewType        vinv,
                         const ordinal_type        operatorDn = 0 );
      };

      template<typename DeviceType, ordinal_type numPtsPerEval,
               typename outputValueValueType, class ...outputValueProperties,
               typename inputPointValueType,  class ...inputPointProperties,
               typename vinvValueType,        class ...vinvProperties>
      static void
      getValues( const typename DeviceType::execution_space& space,
                       Kokkos::DynRankView<outputValueValueType,outputValueProperties...> outputValues,
                 const Kokkos::DynRankView<inputPointValueType, inputPointProperties...>  inputPoints,
                 const Kokkos::DynRankView<vinvValueType,       vinvProperties...>        vinv,
                 const EOperator operatorType );

      /**
         \brief See Intrepid2::Basis_HGRAD_LINE_Cn_FEM
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
        const ordinal_type        _opDn;

        KOKKOS_INLINE_FUNCTION
        Functor(       outputValueViewType outputValues_,
                       inputPointViewType  inputPoints_,
                       vinvViewType        vinv_,
                       workViewType        work_,
                       const ordinal_type        opDn_ = 0 )
          : _outputValues(outputValues_), _inputPoints(inputPoints_), 
            _vinv(vinv_), _work(work_), _opDn(opDn_) {}

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
          case OPERATOR_Dn : {
            auto output = Kokkos::subview( _outputValues, Kokkos::ALL(), ptRange, Kokkos::ALL() );
            Serial<opType>::getValues( output, input, work, _vinv, _opDn );
            break;
          }
          default: {
            INTREPID2_TEST_FOR_ABORT( true,
                                      ">>> ERROR: (Intrepid2::Basis_HGRAD_LINE_Cn_FEM::Functor) operator is not supported");

          }
          }
        }
      };
    };
  }

  template<typename DeviceType = void,
           typename outputValueType = double,
           typename pointValueType = double>
  class Basis_HGRAD_LINE_Cn_FEM
    : public Basis<DeviceType,outputValueType,pointValueType> {
  public:
    using BasisBase = Basis<DeviceType,outputValueType,pointValueType>;
    using typename BasisBase::ExecutionSpace;

    using HostBasis = Basis_HGRAD_LINE_Cn_FEM<typename Kokkos::HostSpace::device_type,outputValueType,pointValueType>;

    using typename BasisBase::OrdinalTypeArray1DHost;
    using typename BasisBase::OrdinalTypeArray2DHost;
    using typename BasisBase::OrdinalTypeArray3DHost;

    using typename BasisBase::OutputViewType;
    using typename BasisBase::PointViewType ;
    using typename BasisBase::ScalarViewType;

  private:

    /** \brief inverse of Generalized Vandermonde matrix, whose columns store the expansion
        coefficients of the nodal basis in terms of phis_ */
    Kokkos::DynRankView<typename ScalarViewType::value_type,DeviceType> vinv_;
    EPointType   pointType_;
  public:
    /** \brief  Constructor.
     */
    Basis_HGRAD_LINE_Cn_FEM(const ordinal_type order,
                            const EPointType   pointType = POINTTYPE_EQUISPACED);  

    using BasisBase::getValues;

    virtual
    void
    getValues( const ExecutionSpace& space,
                     OutputViewType outputValues,
               const PointViewType  inputPoints,
               const EOperator operatorType = OPERATOR_VALUE ) const override {
#ifdef HAVE_INTREPID2_DEBUG
      Intrepid2::getValues_HGRAD_Args(outputValues,
                                      inputPoints,
                                      operatorType,
                                      this->getBaseCellTopology(),
                                      this->getCardinality() );
#endif
      constexpr ordinal_type numPtsPerEval = 1;
      Impl::Basis_HGRAD_LINE_Cn_FEM::
        getValues<DeviceType,numPtsPerEval>(space,
                                            outputValues,
                                            inputPoints,
                                            this->vinv_,
                                            operatorType);
    }

    virtual
    void
    getDofCoords( ScalarViewType dofCoords ) const override {
#ifdef HAVE_INTREPID2_DEBUG
      // Verify rank of output array.
      INTREPID2_TEST_FOR_EXCEPTION( rank(dofCoords) != 2, std::invalid_argument,
                                    ">>> ERROR: (Intrepid2::Basis_HGRAD_LINE_Cn_FEM::getDofCoords) rank = 2 required for dofCoords array");
      // Verify 0th dimension of output array.
      INTREPID2_TEST_FOR_EXCEPTION( static_cast<ordinal_type>(dofCoords.extent(0)) != this->getCardinality(), std::invalid_argument,
                                    ">>> ERROR: (Intrepid2::Basis_HGRAD_LINE_Cn_FEM::getDofCoords) mismatch in number of dof and 0th dimension of dofCoords array");
      // Verify 1st dimension of output array.
      INTREPID2_TEST_FOR_EXCEPTION( dofCoords.extent(1) != this->getBaseCellTopology().getDimension(), std::invalid_argument,
                                    ">>> ERROR: (Intrepid2::Basis_HGRAD_LINE_Cn_FEM::getDofCoords) incorrect reference cell (1st) dimension in dofCoords array");
#endif
      Kokkos::deep_copy(dofCoords, this->dofCoords_);
    }

    virtual
    void
    getDofCoeffs( ScalarViewType dofCoeffs ) const override {
#ifdef HAVE_INTREPID2_DEBUG
      // Verify rank of output array.
      INTREPID2_TEST_FOR_EXCEPTION( rank(dofCoeffs) != 1, std::invalid_argument,
                                    ">>> ERROR: (Intrepid2::Basis_HGRAD_LINE_Cn_FEM::getdofCoeffs) rank = 1 required for dofCoeffs array");
      // Verify 0th dimension of output array.
      INTREPID2_TEST_FOR_EXCEPTION( static_cast<ordinal_type>(dofCoeffs.extent(0)) != this->getCardinality(), std::invalid_argument,
                                    ">>> ERROR: (Intrepid2::Basis_HGRAD_LINE_Cn_FEM::getdofCoeffs) mismatch in number of dof and 0th dimension of dofCoeffs array");
#endif
      Kokkos::deep_copy(dofCoeffs, 1.0);
    }

    virtual
    const char*
    getName() const override {
      return "Intrepid2_HGRAD_LINE_Cn_FEM";
    }

    void
    getVandermondeInverse( ScalarViewType vinv ) const {
      // has to be same rank and dimensions
      Kokkos::deep_copy(vinv, this->vinv_);      
    }

    Kokkos::DynRankView<typename ScalarViewType::const_value_type,DeviceType>
    getVandermondeInverse() const {
      return vinv_;
    }

    ordinal_type
    getWorkSizePerPoint(const EOperator operatorType) const {
      return getPnCardinality<1>(this->basisDegree_);
    }

  /** \brief Creates and returns a Basis object whose DeviceType template argument is Kokkos::HostSpace::device_type, but is otherwise identical to this.

         \return Pointer to the new Basis object.
      */
     virtual HostBasisPtr<outputValueType,pointValueType>
     getHostBasis() const override {
       auto hostBasis = Teuchos::rcp(new HostBasis(this->basisDegree_, pointType_));

       return hostBasis;
     }
  };

}// namespace Intrepid2

#include "Intrepid2_HGRAD_LINE_Cn_FEMDef.hpp"

#endif
