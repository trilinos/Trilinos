// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_HGRAD_WEDGE_C2_FEM.hpp
    \brief  Header file for the Intrepid2::Basis_HGRAD_WEDGE_C2_FEM class.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_HGRAD_WEDGE_C2_FEM_HPP__
#define __INTREPID2_HGRAD_WEDGE_C2_FEM_HPP__

#include "Intrepid2_Basis.hpp"
#include "Intrepid2_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid2_HGRAD_TRI_C2_FEM.hpp"

namespace Intrepid2 {

  /** \class  Intrepid2::Basis_HGRAD_WEDGE_DEG2_FEM
      \brief  Implementation of the default H(grad)-compatible FEM basis of degree 2 on Wedge cell

      Implements Lagrangian basis of degree 2 on the reference Wedge cell. 

      When the serendipity template argument is false, the basis has
      cardinality 18 and spans a COMPLETE bi-quadratic polynomial space. 
      Note, Basis_HGRAD_WEDGE_C2_FEM = Basis_HGRAD_WEDGE_DEG2_FEM<false>

      When the serendipity template argument is true, the basis has
      cardinality 15 and spans an ICOMPLETE bi-quadratic polynomial space. 
      Note, Basis_HGRAD_WEDGE_I2_FEM = Basis_HGRAD_WEDGE_DEG2_FEM<true>
      
      Basis functions are dual to a unisolvent set of degrees-of-freedom (DoF) defined and enumerated as follows:

      \verbatim
      =================================================================================================
      |         |           degree-of-freedom-tag table                    |                           |
      |   DoF   |----------------------------------------------------------|      DoF definition       |
      | ordinal |  subc dim    | subc ordinal | subc DoF ord |subc num DoF |                           |
      |=========|==============|==============|==============|=============|===========================|
      |    0    |       0      |       0      |       0      |      1      |   L_0(u) = u( 0, 0,-1)    |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |    1    |       0      |       1      |       0      |      1      |   L_1(u) = u( 1, 0,-1)    |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |    2    |       0      |       2      |       0      |      1      |   L_2(u) = u( 0, 1,-1)    |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |    3    |       0      |       3      |       0      |      1      |   L_3(u) = u( 0, 0, 1)    |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |    4    |       0      |       4      |       0      |      1      |   L_4(u) = u( 1, 0, 1)    |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |    5    |       0      |       5      |       0      |      1      |   L_5(u) = u( 0, 1, 1)    |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |    6    |       1      |       0      |       0      |      1      |   L_6(u) = u(1/2, 0,-1)   |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |    7    |       1      |       1      |       0      |      1      |   L_7(u) = u(1/2,1/2,-1)  |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |    8    |       1      |       2      |       0      |      1      |   L_8(u) = u( 0,1/2,-1)   |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |    9    |       1      |       6      |       0      |      1      |   L_9(u) = u( 0, 0, 0)    |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |   10    |       1      |       7      |       0      |      1      |   L_10(u)= u( 1, 0, 0)    |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |   11    |       1      |       8      |       0      |      1      |   L_11(u)= u( 0, 1, 0)    |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |   12    |       1      |       3      |       0      |      1      |   L_12(u)= u(1/2, 0, 1)   |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |   13    |       1      |       4      |       0      |      1      |   L_13(u)= u(1/2,1/2, 1)  |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |   14    |       1      |       5      |       0      |      1      |   L_14(u)= u( 0,1/2, 1)   |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |   15    |       2      |       0      |       0      |      1      |   L_15(u)= u(1/2, 0, 0)   |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |   16    |       2      |       1      |       0      |      1      |   L_16(u)= u(1/2,1/2, 0)  |
      |---------|--------------|--------------|--------------|-------------|---------------------------|
      |   17    |       2      |       2      |       0      |      1      |   L_17(u)= u( 0,1/2, 0)   |
      |=========|==============|==============|==============|=============|===========================|
      |   MAX   |  maxScDim=2  |  maxScOrd=8  |  maxDfOrd=0  |      -      |                           |
      |=========|==============|==============|==============|=============|===========================|
      \endverbatim

      \remark   Ordering of DoFs follows the node order in Wedge<18> topology. Note that node
      order in this topology does not follow the natural oder of k-subcells where the nodes
      are located, except for nodes 0 to 5 which coincide with the vertices of the base
      Wedge<6> topology. As a result, L_0 to L_5 are associated with nodes 0 to 5, but
      L_6 to L_14 are not associated with edges 0 to 9 in that order. The last three nodes,
      for the non sereendipity basis, are located on 2-subcells (faces) and follow their order. 
      Thus, L_15, L_16 and L17 are associated with faces 0, 1 and 2 in that order.
  */

  namespace Impl {

    /**
      \brief See Intrepid2::Basis_HGRAD_WEDGE_DEG2_FEM
    */
    template<bool serendipity>
    class Basis_HGRAD_WEDGE_DEG2_FEM {
    public:
      typedef struct Wedge<serendipity ? 15 : 18> cell_topology_type;
      /**
        \brief See Intrepid2::Basis_HGRAD_WEDGE_DEG2_FEM
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
        \brief See Intrepid2::Basis_HGRAD_WEDGE_DEG2_FEM
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
                                      opType != OPERATOR_D2 &&
                                      opType != OPERATOR_MAX,
                                      ">>> ERROR: (Intrepid2::Basis_HGRAD_WEDGE_DEG2_FEM::Serial::getValues) operator is not supported");
          }
          }
        }
      };
    };
  }

  template<bool serendipity,
           typename DeviceType = void,
           typename outputValueType = double,
           typename pointValueType = double>
  class Basis_HGRAD_WEDGE_DEG2_FEM : public Basis<DeviceType,outputValueType,pointValueType> {
  public:
    using OrdinalTypeArray1DHost = typename Basis<DeviceType,outputValueType,pointValueType>::OrdinalTypeArray1DHost;
    using OrdinalTypeArray2DHost = typename Basis<DeviceType,outputValueType,pointValueType>::OrdinalTypeArray2DHost;
    using OrdinalTypeArray3DHost = typename Basis<DeviceType,outputValueType,pointValueType>::OrdinalTypeArray3DHost;
    /** \brief  Constructor.
     */
    Basis_HGRAD_WEDGE_DEG2_FEM();

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
      if constexpr (serendipity)
        Impl::Basis_HGRAD_WEDGE_DEG2_FEM<true>::
          getValues<DeviceType>( outputValues,
                                    inputPoints,
                                    operatorType );
      else 
        Impl::Basis_HGRAD_WEDGE_DEG2_FEM<false>::
          getValues<DeviceType>( outputValues,
                                    inputPoints,
                                    operatorType );;
    }

    virtual
    void
    getDofCoords( ScalarViewType dofCoords ) const override {
#ifdef HAVE_INTREPID2_DEBUG
      // Verify rank of output array.
      INTREPID2_TEST_FOR_EXCEPTION( dofCoords.rank() != 2, std::invalid_argument,
                                    ">>> ERROR: (Intrepid2::Basis_HGRAD_WEDGE_DEG2_FEM::getDofCoords) rank = 2 required for dofCoords array");
      // Verify 0th dimension of output array.
      INTREPID2_TEST_FOR_EXCEPTION( static_cast<ordinal_type>(dofCoords.extent(0)) != this->getCardinality(), std::invalid_argument,
                                    ">>> ERROR: (Intrepid2::Basis_HGRAD_WEDGE_DEG2_FEM::getDofCoords) mismatch in number of dof and 0th dimension of dofCoords array");
      // Verify 1st dimension of output array.
      INTREPID2_TEST_FOR_EXCEPTION( dofCoords.extent(1) != this->getBaseCellTopology().getDimension(), std::invalid_argument,
                                    ">>> ERROR: (Intrepid2::Basis_HGRAD_WEDGE_DEG2_FEM::getDofCoords) incorrect reference cell (1st) dimension in dofCoords array");
#endif
      Kokkos::deep_copy(dofCoords, this->dofCoords_);
    }

    virtual
    void
    getDofCoeffs( ScalarViewType dofCoeffs ) const override {
#ifdef HAVE_INTREPID2_DEBUG
      // Verify rank of output array.
      INTREPID2_TEST_FOR_EXCEPTION( dofCoeffs.rank() != 1, std::invalid_argument,
                                    ">>> ERROR: (Intrepid2::Basis_HGRAD_WEDGE_DEG2_FEM::getdofCoeffs) rank = 1 required for dofCoeffs array");
      // Verify 0th dimension of output array.
      INTREPID2_TEST_FOR_EXCEPTION( static_cast<ordinal_type>(dofCoeffs.extent(0)) != this->getCardinality(), std::invalid_argument,
                                    ">>> ERROR: (Intrepid2::Basis_HGRAD_WEDGE_DEG2_FEM::getdofCoeffs) mismatch in number of dof and 0th dimension of dofCoeffs array");
#endif
      Kokkos::deep_copy(dofCoeffs, 1.0);
    }

    virtual
    const char*
    getName() const override {
      return serendipity ? "Intrepid2_HGRAD_WEDGE_I2_FEM" : "Intrepid2_HGRAD_WEDGE_C2_FEM";
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
        if(subCellOrd < 3)  //lateral faces
          return Teuchos::rcp(new Basis_HGRAD_QUAD_DEG2_FEM<serendipity, DeviceType, outputValueType, pointValueType>());
        else
          return Teuchos::rcp(new Basis_HGRAD_TRI_C2_FEM<DeviceType, outputValueType, pointValueType>());
      }
      INTREPID2_TEST_FOR_EXCEPTION(true,std::invalid_argument,"Input parameters out of bounds");
    }

    BasisPtr<typename Kokkos::HostSpace::device_type,outputValueType,pointValueType>
    getHostBasis() const override{
      return Teuchos::rcp(new Basis_HGRAD_WEDGE_DEG2_FEM<serendipity, typename Kokkos::HostSpace::device_type,outputValueType,pointValueType>());
    }
  };

  template<typename DeviceType = void, typename outputValueType = double, typename pointValueType = double>
  using Basis_HGRAD_WEDGE_C2_FEM = Basis_HGRAD_WEDGE_DEG2_FEM<false, DeviceType, outputValueType, pointValueType>;

  template<typename DeviceType = void, typename outputValueType = double, typename pointValueType = double>
  using Basis_HGRAD_WEDGE_I2_FEM = Basis_HGRAD_WEDGE_DEG2_FEM<true, DeviceType, outputValueType, pointValueType>;

}// namespace Intrepid2

#include "Intrepid2_HGRAD_WEDGE_C2_FEMDef.hpp"

#endif
