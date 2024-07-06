// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_HGRAD_QUAD_C2_FEM.hpp
    \brief  Header file for the Intrepid2::Basis_HGRAD_QUAD_C2_FEM class.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
 */

#ifndef __INTREPID2_HGRAD_QUAD_C2_FEM_HPP__
#define __INTREPID2_HGRAD_QUAD_C2_FEM_HPP__

#include "Intrepid2_Basis.hpp"
#include "Intrepid2_HGRAD_LINE_C2_FEM.hpp"

namespace Intrepid2 {

  /** \class  Intrepid2::Basis_HGRAD_QUAD_DEG2_FEM
      \brief  Implementation of the default H(grad)-compatible FEM basis of degree 2 on Quadrilateral cell

      Implements Lagrangian basis of degree 2 on the reference Quadrilateral cell. 
      When the serendipity template argument is false, the basis has
      cardinality 9 and spans a COMPLETE bi-quadratic polynomial space. 
      Note, Basis_HGRAD_QUAD_C2_FEM = Basis_HGRAD_QUAD_DEG2_FEM<false>

      When the serendipity template argument is true, the basis has
      cardinality 8 and spans an INCOMPLETE bi-quadratic polynomial space (Dofs are associated only to vertices and edges). 
      Note, Basis_HGRAD_QUAD_I2_FEM = Basis_HGRAD_QUAD_C2_FEM<true>      
      
      Basis functions are dual
      to a unisolvent set of degrees-of-freedom (DoF) defined and enumerated as follows:

      \verbatim
      =================================================================================================
      |         |           degree-of-freedom-tag table                    |                           |
      |   DoF   |----------------------------------------------------------|      DoF definition       |
      | ordinal |  subc dim    | subc ordinal | subc DoF ord |subc num DoF |                           |
      |=========|==============|==============|==============|=============|===========================|
      |    0    |       0      |       0      |       0      |      1      |   L_0(u) = u(-1,-1)       |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |    1    |       0      |       1      |       0      |      1      |   L_1(u) = u( 1,-1)       |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |    2    |       0      |       2      |       0      |      1      |   L_2(u) = u( 1, 1)       |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |    3    |       0      |       3      |       0      |      1      |   L_3(u) = u(-1, 1)       |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |    4    |       1      |       0      |       0      |      1      |   L_4(u) = u( 0,-1)       |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |    5    |       1      |       1      |       0      |      1      |   L_5(u) = u( 1, 0)       |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |    6    |       1      |       2      |       0      |      1      |   L_6(u) = u( 0, 1)       |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |    7    |       1      |       3      |       0      |      1      |   L_7(u) = u(-1, 0)       |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |    8    |       2      |       0      |       0      |      1      |   L_8(u) = u( 0, 0)       |
      |=========|==============|==============|==============|=============|===========================|
      |   MAX   |  maxScDim=2  |  maxScOrd=3  |  maxDfOrd=0  |     -       |                           |
      |=========|==============|==============|==============|=============|===========================|
      \endverbatim
  */


  namespace Impl {

    /**
      \brief See Intrepid2::Basis_HGRAD_QUAD_DEG2_FEM
    */
    template<bool serendipity>
    class Basis_HGRAD_QUAD_DEG2_FEM {
    public:
      typedef struct Quadrilateral<serendipity ? 8 : 9> cell_topology_type;
      /**
        \brief See Intrepid2::Basis_HGRAD_QUAD_DEG2_FEM
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

      template<typename DeviceType,
               typename outputValueValueType, class ...outputValueProperties,
               typename inputPointValueType,  class ...inputPointProperties>
      static void
      getValues( const typename DeviceType::execution_space& space,
                       Kokkos::DynRankView<outputValueValueType,outputValueProperties...> outputValues,
                 const Kokkos::DynRankView<inputPointValueType, inputPointProperties...>  inputPoints,
                 const EOperator operatorType);

      /**
        \brief See Intrepid2::Basis_HGRAD_QUAD_DEG2_FEM
      */
      template<typename outputValueViewType,
               typename inputPointViewType,
               EOperator opType>
      struct Functor {
              outputValueViewType _outputValues;
        const inputPointViewType  _inputPoints;

        KOKKOS_INLINE_FUNCTION
        Functor(      outputValueViewType outputValues_,
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
          case OPERATOR_CURL :
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
                                      opType != OPERATOR_CURL &&
                                      opType != OPERATOR_D1 &&
                                      opType != OPERATOR_D2 &&
                                      opType != OPERATOR_D3 &&
                                      opType != OPERATOR_D4 &&
                                      opType != OPERATOR_MAX,
                                      ">>> ERROR: (Intrepid2::Basis_HGRAD_QUAD_DEG2_FEM::Serial::getValues) operator is not supported");
          }
          }
        }
      };
    };
  }

  template<bool serendipity,
           typename DeviceType,
           typename outputValueType,
           typename pointValueType>
  class Basis_HGRAD_QUAD_DEG2_FEM : public Basis<DeviceType,outputValueType,pointValueType> {
  public:
    using BasisBase = Basis<DeviceType, outputValueType, pointValueType>;
    using typename BasisBase::ExecutionSpace;
    using typename BasisBase::OrdinalTypeArray1DHost;
    using typename BasisBase::OrdinalTypeArray2DHost;
    using typename BasisBase::OrdinalTypeArray3DHost;

    /** \brief Constructor.
     */
    Basis_HGRAD_QUAD_DEG2_FEM();

    using typename BasisBase::OutputViewType;
    using typename BasisBase::PointViewType;
    using typename BasisBase::ScalarViewType;

    using BasisBase::getValues;

    virtual
    void
    getValues( const ExecutionSpace& space,
                     OutputViewType outputValues,
               const PointViewType  inputPoints,
               const EOperator operatorType = OPERATOR_VALUE ) const override {
#ifdef HAVE_INTREPID2_DEBUG
      // Verify arguments
      Intrepid2::getValues_HGRAD_Args(outputValues,
                                      inputPoints,
                                      operatorType,
                                      this->getBaseCellTopology(),
                                      this->getCardinality() );
#endif
      Impl::Basis_HGRAD_QUAD_DEG2_FEM<serendipity>::
          template getValues<DeviceType>(space,
                                outputValues,
                                inputPoints,
                                operatorType);
    }

    virtual
    void
    getDofCoords( ScalarViewType dofCoords ) const override {
#ifdef HAVE_INTREPID2_DEBUG
      // Verify rank of output array.
      INTREPID2_TEST_FOR_EXCEPTION( rank(dofCoords) != 2, std::invalid_argument,
                                    ">>> ERROR: (Intrepid2::Basis_HGRAD_QUAD_DEG2_FEM::getDofCoords) rank = 2 required for dofCoords array");
      // Verify 0th dimension of output array.
      INTREPID2_TEST_FOR_EXCEPTION( static_cast<ordinal_type>(dofCoords.extent(0)) != this->getCardinality(), std::invalid_argument,
                                    ">>> ERROR: (Intrepid2::Basis_HGRAD_QUAD_DEG2_FEM::getDofCoords) mismatch in number of dof and 0th dimension of dofCoords array");
      // Verify 1st dimension of output array.
      INTREPID2_TEST_FOR_EXCEPTION( dofCoords.extent(1) != this->getBaseCellTopology().getDimension(), std::invalid_argument,
                                    ">>> ERROR: (Intrepid2::Basis_HGRAD_QUAD_DEG2_FEM::getDofCoords) incorrect reference cell (1st) dimension in dofCoords array");
#endif
      Kokkos::deep_copy(dofCoords, this->dofCoords_);
    }

    virtual
    void
    getDofCoeffs( ScalarViewType dofCoeffs ) const override {
#ifdef HAVE_INTREPID2_DEBUG
      // Verify rank of output array.
      INTREPID2_TEST_FOR_EXCEPTION( rank(dofCoeffs) != 1, std::invalid_argument,
                                    ">>> ERROR: (Intrepid2::Basis_HGRAD_QUAD_DEG2_FEM::getdofCoeffs) rank = 1 required for dofCoeffs array");
      // Verify 0th dimension of output array.
      INTREPID2_TEST_FOR_EXCEPTION( static_cast<ordinal_type>(dofCoeffs.extent(0)) != this->getCardinality(), std::invalid_argument,
                                    ">>> ERROR: (Intrepid2::Basis_HGRAD_QUAD_DEG2_FEM::getdofCoeffs) mismatch in number of dof and 0th dimension of dofCoeffs array");
#endif
      Kokkos::deep_copy(dofCoeffs, 1.0);
    }

    virtual
    const char*
    getName() const override {
      return serendipity ? "Intrepid2_HGRAD_QUAD_I2_FEM" : "Intrepid2_HGRAD_QUAD_C2_FEM";
    }

    /** \brief returns the basis associated to a subCell.

        The bases of the subCell are the restriction to the subCell
        of the bases of the parent cell.
        \param [in] subCellDim - dimension of subCell
        \param [in] subCellOrd - position of the subCell among of the subCells having the same dimension
        \return pointer to the subCell basis of dimension subCellDim and position subCellOrd
     */
    BasisPtr<DeviceType, outputValueType, pointValueType>
      getSubCellRefBasis(const ordinal_type subCellDim, const ordinal_type subCellOrd) const override{
      if(subCellDim == 1) {
        return Teuchos::rcp(new
            Basis_HGRAD_LINE_C2_FEM<DeviceType, outputValueType, pointValueType>());
      }
      INTREPID2_TEST_FOR_EXCEPTION(true,std::invalid_argument,"Input parameters out of bounds");
    }

    BasisPtr<typename Kokkos::HostSpace::device_type,outputValueType,pointValueType>
    getHostBasis() const override{
      return Teuchos::rcp(new Basis_HGRAD_QUAD_DEG2_FEM<serendipity, typename Kokkos::HostSpace::device_type,outputValueType,pointValueType>());
    }  
  };

  template<typename DeviceType = void, typename outputValueType = double, typename pointValueType = double>
  using Basis_HGRAD_QUAD_C2_FEM = Basis_HGRAD_QUAD_DEG2_FEM<false, DeviceType, outputValueType, pointValueType>;

  template<typename DeviceType = void, typename outputValueType = double, typename pointValueType = double>
  using Basis_HGRAD_QUAD_I2_FEM = Basis_HGRAD_QUAD_DEG2_FEM<true, DeviceType, outputValueType, pointValueType>;

}// namespace Intrepid2

#include "Intrepid2_HGRAD_QUAD_C2_FEMDef.hpp"

#endif
