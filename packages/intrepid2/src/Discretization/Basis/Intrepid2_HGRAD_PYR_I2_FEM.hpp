// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_HGRAD_PYR_I2_FEM.hpp
    \brief  Header file for the Intrepid2::Basis_HGRAD_PYR_I2_FEM class.
    \author Created by P. Bochev, D. Ridzal and M. Perego
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_HGRAD_PYR_I2_SERENDIPITY_FEM_HPP__
#define __INTREPID2_HGRAD_PYR_I2_SERENDIPITY_FEM_HPP__

#include "Intrepid2_Basis.hpp"
#include "Intrepid2_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid2_HGRAD_TRI_C2_FEM.hpp"

namespace Intrepid2 {

  /** \class  Intrepid2::Basis_HGRAD_PYR_I2_FEM
      \brief  Implementation of an H(grad)-compatible FEM basis of degree 2 on a Pyramid cell

      Implements Lagrangian basis of degree 2 on the reference Pyramid cell. The basis has
      cardinality 13 and spans an INCOMPLETE bi-quadratic polynomial space. Basis functions are dual 
      to a unisolvent set of degrees-of-freedom (DoF) defined and enumerated as follows:

      \verbatim
      =================================================================================================
      |         |           degree-of-freedom-tag table                    |                           |
      |   DoF   |----------------------------------------------------------|      DoF definition       |
      | ordinal |  subc dim    | subc ordinal | subc DoF ord |subc num DoF |                           |
      |=========|==============|==============|==============|=============|===========================|
      |    0    |       0      |       0      |       0      |      1      |   L_0(u) = u(-1,-1, 0)    |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |    1    |       0      |       1      |       0      |      1      |   L_1(u) = u( 1,-1, 0)    |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |    2    |       0      |       2      |       0      |      1      |   L_2(u) = u( 1, 1, 0)    |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |    3    |       0      |       3      |       0      |      1      |   L_3(u) = u(-1, 1, 0)    |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |    4    |       0      |       4      |       0      |      1      |   L_4(u) = u( 0, 0, 1)    |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |    5    |       1      |       0      |       0      |      1      |   L_5(u) = u( 0,-1, 0)    |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |    6    |       1      |       1      |       0      |      1      |   L_6(u) = u( 1, 0, 0)    |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |    7    |       1      |       2      |       0      |      1      |   L_7(u) = u( 0, 1, 0)    |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |    8    |       1      |       3      |       0      |      1      |   L_8(u) = u(-1, 0, 0)    |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |    9    |       1      |       4      |       0      |      1      | L_9(u) = u(-1/2,-1/2,1/2) |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |   10    |       1      |       5      |       0      |      1      | L_10(u)= u( 1/2,-1/2,1/2) |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |   11    |       1      |       6      |       0      |      1      | L_11(u)= u( 1/2, 1/2,1/2) |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |   12    |       1      |       7      |       0      |      1      | L_12(u)= u(-1/2, 1/2,1/2) |
      |=========|==============|==============|==============|=============|===========================|
      |   MAX   |  maxScDim=1  |  maxScOrd=7  |  maxDfOrd=0  |      -      |                           |
      |=========|==============|==============|==============|=============|===========================|
      \endverbatim
  */

  namespace Impl {

    /**
      \brief See Intrepid2::Basis_HGRAD_PYR_I2_FEM
    */
    class Basis_HGRAD_PYR_I2_FEM {
    public:
      typedef struct Pyramid<5> cell_topology_type;
      /**
        \brief See Intrepid2::Basis_HGRAD_PYR_I2_FEM
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
      getValues(       Kokkos::DynRankView<outputValueValueType,outputValueProperties...> outputValues,
                 const Kokkos::DynRankView<inputPointValueType, inputPointProperties...>  inputPoints,
                 const EOperator operatorType);

      /**
        \brief See Intrepid2::Basis_HGRAD_PYR_I2_FEM
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
          case OPERATOR_CURL :
          case OPERATOR_D2 : {
            auto       output = Kokkos::subview( _outputValues, Kokkos::ALL(), pt, Kokkos::ALL() );
            const auto input  = Kokkos::subview( _inputPoints,                 pt, Kokkos::ALL() );
            Serial<opType>::getValues( output, input );
            break;
          }
          default: {
            INTREPID2_TEST_FOR_ABORT( opType != OPERATOR_VALUE &&
                                      opType != OPERATOR_GRAD &&
                                      opType != OPERATOR_CURL &&
                                      opType != OPERATOR_D2,
                                      ">>> ERROR: (Intrepid2::Basis_HGRAD_PYR_I2_FEM::Serial::getValues) operator is not supported");
          }
          }
        }
      };

    };
  }
  template<typename DeviceType = void,
           typename outputValueType = double,
           typename pointValueType = double>
  class Basis_HGRAD_PYR_I2_FEM : public Basis<DeviceType,outputValueType,pointValueType> {
  public:
    using OrdinalTypeArray1DHost = typename Basis<DeviceType,outputValueType,pointValueType>::OrdinalTypeArray1DHost;
    using OrdinalTypeArray2DHost = typename Basis<DeviceType,outputValueType,pointValueType>::OrdinalTypeArray2DHost;
    using OrdinalTypeArray3DHost = typename Basis<DeviceType,outputValueType,pointValueType>::OrdinalTypeArray3DHost;

    /** \brief  Constructor.
     */
    Basis_HGRAD_PYR_I2_FEM();

    using OutputViewType = typename Basis<DeviceType,outputValueType,pointValueType>::OutputViewType;
    using PointViewType  = typename Basis<DeviceType,outputValueType,pointValueType>::PointViewType;
    using ScalarViewType = typename Basis<DeviceType,outputValueType,pointValueType>::ScalarViewType;

    using Basis<DeviceType,outputValueType,pointValueType>::getValues;

    virtual
    void
    getValues(       OutputViewType outputValues,
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
      Impl::Basis_HGRAD_PYR_I2_FEM::
        getValues<DeviceType>( outputValues,
                                  inputPoints,
                                  operatorType );
    }

    virtual
    void
    getDofCoords( ScalarViewType dofCoords ) const override {
#ifdef HAVE_INTREPID2_DEBUG
      // Verify rank of output array.
      INTREPID2_TEST_FOR_EXCEPTION( dofCoords.rank() != 2, std::invalid_argument,
                                    ">>> ERROR: (Intrepid2::Basis_HGRAD_PYR_I2_FEM::getDofCoords) rank = 2 required for dofCoords array");
      // Verify 0th dimension of output array.
      INTREPID2_TEST_FOR_EXCEPTION( static_cast<ordinal_type>(dofCoords.extent(0)) != this->getCardinality(), std::invalid_argument,
                                    ">>> ERROR: (Intrepid2::Basis_HGRAD_PYR_I2_FEM::getDofCoords) mismatch in number of dof and 0th dimension of dofCoords array");
      // Verify 1st dimension of output array.
      INTREPID2_TEST_FOR_EXCEPTION( dofCoords.extent(1) != this->getBaseCellTopology().getDimension(), std::invalid_argument,
                                    ">>> ERROR: (Intrepid2::Basis_HGRAD_PYR_I2_FEM::getDofCoords) incorrect reference cell (1st) dimension in dofCoords array");
#endif
      Kokkos::deep_copy(dofCoords, this->dofCoords_);
    }

    virtual
    void
    getDofCoeffs( ScalarViewType dofCoeffs ) const override {
#ifdef HAVE_INTREPID2_DEBUG
      // Verify rank of output array.
      INTREPID2_TEST_FOR_EXCEPTION( dofCoeffs.rank() != 1, std::invalid_argument,
                                    ">>> ERROR: (Intrepid2::Basis_HGRAD_PYR_I2_FEM::getdofCoeffs) rank = 1 required for dofCoeffs array");
      // Verify 0th dimension of output array.
      INTREPID2_TEST_FOR_EXCEPTION( static_cast<ordinal_type>(dofCoeffs.extent(0)) != this->getCardinality(), std::invalid_argument,
                                    ">>> ERROR: (Intrepid2::Basis_HGRAD_PYR_I2_FEM::getdofCoeffs) mismatch in number of dof and 0th dimension of dofCoeffs array");
#endif
      Kokkos::deep_copy(dofCoeffs, 1.0);
    }

    virtual
    const char*
    getName() const override {
      return "Intrepid2_HGRAD_PYR_I2_FEM";
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
            Basis_HGRAD_LINE_C2_FEM<DeviceType,outputValueType,pointValueType>());
      } else if(subCellDim == 2) {
        if(subCellOrd == 0)  //pyramid base
          return Teuchos::rcp(new
            Basis_HGRAD_QUAD_I2_FEM<DeviceType,outputValueType,pointValueType>());
        else
          return Teuchos::rcp(new
            Basis_HGRAD_TRI_C2_FEM<DeviceType,outputValueType,pointValueType>());
      }
      INTREPID2_TEST_FOR_EXCEPTION(true,std::invalid_argument,"Input parameters out of bounds");
    }

    BasisPtr<typename Kokkos::HostSpace::device_type,outputValueType,pointValueType>
    getHostBasis() const override{
      return Teuchos::rcp(new Basis_HGRAD_PYR_I2_FEM<typename Kokkos::HostSpace::device_type,outputValueType,pointValueType>());
    }

  };
}// namespace Intrepid2

#include "Intrepid2_HGRAD_PYR_I2_FEMDef.hpp"

#endif
