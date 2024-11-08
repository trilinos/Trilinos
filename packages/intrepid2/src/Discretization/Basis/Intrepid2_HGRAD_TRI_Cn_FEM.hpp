// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
#include "Intrepid2_HGRAD_LINE_Cn_FEM.hpp"

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
      template<EOperator OpType>
      struct Serial {
        template<typename OutputValueViewType,
                 typename InputPointViewType,
                 typename WorkViewType,
                 typename VinvViewType>
        KOKKOS_INLINE_FUNCTION
        static void
        getValues(       OutputValueViewType outputValues,
                         const InputPointViewType  inputPoints,
                         WorkViewType              work,
                         const VinvViewType        vinv,
                         const ordinal_type        order);
      };

      template<typename DeviceType, ordinal_type numPtsPerEval,
               typename OutputValueValueType, class ...OutputValueProperties,
               typename InputPointValueType,  class ...InputPointProperties,
               typename VinvValueType,        class ...VinvProperties>
      static void
      getValues( const typename DeviceType::execution_space& space,
                       Kokkos::DynRankView<OutputValueValueType,OutputValueProperties...> outputValues,
                 const Kokkos::DynRankView<InputPointValueType, InputPointProperties...>  inputPoints,
                 const Kokkos::DynRankView<VinvValueType,       VinvProperties...>        vinv,
                 const ordinal_type order,
                 const EOperator operatorType);

      /**
         \brief See Intrepid2::Basis_HGRAD_TRI_Cn_FEM
      */
      template<typename OutputValueViewType,
               typename InputPointViewType,
               typename VinvViewType,
               typename WorkViewType,
               EOperator OpType,
               ordinal_type numPtsEval>
      struct Functor {
        OutputValueViewType _outputValues;
        const InputPointViewType  _inputPoints;
        const VinvViewType        _vinv;
        WorkViewType              _work;
        const ordinal_type        _order;

        KOKKOS_INLINE_FUNCTION
        Functor(       OutputValueViewType outputValues_,
                       InputPointViewType  inputPoints_,
                       VinvViewType        vinv_,
                       WorkViewType        work_,
                       ordinal_type        order_)
          : _outputValues(outputValues_), _inputPoints(inputPoints_),
            _vinv(vinv_), _work(work_), _order(order_) {}

        KOKKOS_INLINE_FUNCTION
        void operator()(const size_type iter) const {
          const auto ptBegin = Util<ordinal_type>::min(iter*numPtsEval,    _inputPoints.extent(0));
          const auto ptEnd   = Util<ordinal_type>::min(ptBegin+numPtsEval, _inputPoints.extent(0));

          const auto ptRange = Kokkos::pair<ordinal_type,ordinal_type>(ptBegin, ptEnd);
          const auto input   = Kokkos::subview( _inputPoints, ptRange, Kokkos::ALL() );

          typename WorkViewType::pointer_type ptr = _work.data() + _work.extent(0)*ptBegin*get_dimension_scalar(_work);

          auto vcprop = Kokkos::common_view_alloc_prop(_work);
          WorkViewType  work(Kokkos::view_wrap(ptr,vcprop), (ptEnd-ptBegin)*_work.extent(0));

          switch (OpType) {
          case OPERATOR_VALUE : {
            auto output = Kokkos::subview( _outputValues, Kokkos::ALL(), ptRange );
            Serial<OpType>::getValues( output, input, work, _vinv, _order );
            break;
          }
          case OPERATOR_CURL:
          case OPERATOR_D1:
          case OPERATOR_D2:  {
            auto output = Kokkos::subview( _outputValues, Kokkos::ALL(), ptRange, Kokkos::ALL() );
            Serial<OpType>::getValues( output, input, work, _vinv, _order );
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

  template<typename DeviceType = void,
           typename outputValueType = double,
           typename pointValueType = double>
  class Basis_HGRAD_TRI_Cn_FEM
    : public Basis<DeviceType,outputValueType,pointValueType> {
  public:
    using BasisBase = Basis<DeviceType, outputValueType, pointValueType>;
    using HostBasis = Basis_HGRAD_TRI_Cn_FEM<typename Kokkos::HostSpace::device_type,outputValueType,pointValueType>;
    using typename BasisBase::ExecutionSpace;
    using typename BasisBase::OrdinalTypeArray1DHost;
    using typename BasisBase::OrdinalTypeArray2DHost;
    using typename BasisBase::OrdinalTypeArray3DHost;

    using typename BasisBase::OutputViewType;
    using typename BasisBase::PointViewType;
    using typename BasisBase::ScalarViewType;

    using typename BasisBase::scalarType;

  private:

    /** \brief inverse of Generalized Vandermonde matrix, whose columns store the expansion
        coefficients of the nodal basis in terms of phis_ */
    Kokkos::DynRankView<scalarType,DeviceType> vinv_;

    /** \brief type of lattice used for creating the DoF coordinates  */
    EPointType pointType_;

  public:
    /** \brief  Constructor.
     */
    Basis_HGRAD_TRI_Cn_FEM(const ordinal_type order,
                           const EPointType   pointType = POINTTYPE_EQUISPACED);

    using BasisBase::getValues;

    virtual
    void
    getValues( const ExecutionSpace& space,
                     OutputViewType outputValues,
               const PointViewType  inputPoints,
               const EOperator operatorType = OPERATOR_VALUE) const override {
#ifdef HAVE_INTREPID2_DEBUG
      Intrepid2::getValues_HGRAD_Args(outputValues,
                                      inputPoints,
                                      operatorType,
                                      this->getBaseCellTopology(),
                                      this->getCardinality() );
#endif
      constexpr ordinal_type numPtsPerEval = Parameters::MaxNumPtsPerBasisEval;
      Impl::Basis_HGRAD_TRI_Cn_FEM::
        getValues<DeviceType,numPtsPerEval>(space,
                                            outputValues,
                                            inputPoints,
                                            this->vinv_,
                                            this->basisDegree_,
                                            operatorType);
    }

    virtual void 
    getScratchSpaceSize(      ordinal_type& perTeamSpaceSize,
                              ordinal_type& perThreadSpaceSize,
                        const PointViewType inputPoints,
                        const EOperator operatorType = OPERATOR_VALUE) const override;

    KOKKOS_INLINE_FUNCTION
    virtual void 
    getValues(       
          OutputViewType outputValues,
      const PointViewType  inputPoints,
      const EOperator operatorType,
      const typename Kokkos::TeamPolicy<typename DeviceType::execution_space>::member_type& team_member,
      const typename DeviceType::execution_space::scratch_memory_space & scratchStorage, 
      const ordinal_type subcellDim = -1,
      const ordinal_type subcellOrdinal = -1) const override;



    virtual
    void
    getDofCoords( ScalarViewType dofCoords ) const override {
#ifdef HAVE_INTREPID2_DEBUG
      // Verify rank of output array.
      INTREPID2_TEST_FOR_EXCEPTION( rank(dofCoords) != 2, std::invalid_argument,
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
    getDofCoeffs( ScalarViewType dofCoeffs ) const override {
#ifdef HAVE_INTREPID2_DEBUG
      // Verify rank of output array.
      INTREPID2_TEST_FOR_EXCEPTION( rank(dofCoeffs) != 1, std::invalid_argument,
                                    ">>> ERROR: (Intrepid2::Basis_HGRAD_TRI_Cn_FEM::getdofCoeffs) rank = 1 required for dofCoeffs array");
      // Verify 0th dimension of output array.
      INTREPID2_TEST_FOR_EXCEPTION( static_cast<ordinal_type>(dofCoeffs.extent(0)) != this->getCardinality(), std::invalid_argument,
                                    ">>> ERROR: (Intrepid2::Basis_HGRAD_TRI_Cn_FEM::getdofCoeffs) mismatch in number of dof and 0th dimension of dofCoeffs array");
#endif
      Kokkos::deep_copy(dofCoeffs, 1.0);
    }


    virtual
    const char*
    getName() const override {
      return "Intrepid2_HGRAD_TRI_Cn_FEM";
    }

    virtual
    bool
    requireOrientation() const override {
      return (this->basisDegree_ > 2);
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
    
    /** \brief returns the basis associated to a subCell.

        The bases of the subCell are the restriction to the subCell
        of the bases of the parent cell.
        \param [in] subCellDim - dimension of subCell
        \param [in] subCellOrd - position of the subCell among of the subCells having the same dimension
        \return pointer to the subCell basis of dimension subCellDim and position subCellOrd
     */
    BasisPtr<DeviceType,outputValueType,pointValueType>
      getSubCellRefBasis(const ordinal_type subCellDim, const ordinal_type subCellOrd) const override{
      if(subCellDim == 1) {
        return Teuchos::rcp(new
            Basis_HGRAD_LINE_Cn_FEM<DeviceType,outputValueType,pointValueType>
            (this->basisDegree_, pointType_));
      }
      INTREPID2_TEST_FOR_EXCEPTION(true,std::invalid_argument,"Input parameters out of bounds");
    }

    BasisPtr<typename Kokkos::HostSpace::device_type,outputValueType,pointValueType>
    getHostBasis() const override{
      return Teuchos::rcp(new Basis_HGRAD_TRI_Cn_FEM<typename Kokkos::HostSpace::device_type,outputValueType,pointValueType>(this->basisDegree_, pointType_));
    }
  };

}// namespace Intrepid2

#include "Intrepid2_HGRAD_TRI_Cn_FEMDef.hpp"

#endif
