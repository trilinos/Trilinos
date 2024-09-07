// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_HDIV_HEX_I1_FEM.hpp
    \brief  Header file for the Intrepid2::Basis_HDIV_HEX_I1_FEM class.
    \author Created by P. Bochev and D. Ridzal and K. Petrson.
            Kokkorized by Kyungjoo Kim
 */

#ifndef __INTREPID2_HDIV_HEX_I1_FEM_HPP__
#define __INTREPID2_HDIV_HEX_I1_FEM_HPP__

#include "Intrepid2_Basis.hpp"
#include "Intrepid2_HVOL_C0_FEM.hpp"

namespace Intrepid2 {

  /** \class  Intrepid2::Basis_HDIV_HEX_I1_FEM
      \brief  Implementation of the default H(div)-compatible FEM basis of degree 1 on Hexahedron cell

      Implements Raviart-Thomas basis of degree 1 on the reference Hexahedron cell. The basis has
      cardinality 6 and spans a INCOMPLETE tri-linear polynomial space. Basis functions are dual
      to a unisolvent set of degrees-of-freedom (DoF) defined and enumerated as follows:

      \verbatim
      ===================================================================================================
      |         |           degree-of-freedom-tag table                    |                            |
      |   DoF   |----------------------------------------------------------|       DoF definition       |
      | ordinal |  subc dim    | subc ordinal | subc DoF ord |subc num DoF |                            |
      |=========|==============|==============|==============|=============|============================|
      |    0    |       2      |       0      |       0      |      1      |   L_0(u) = (u.n)(0,-1,0)   |
      |---------|--------------|--------------|--------------|-------------|----------------------------|
      |    1    |       2      |       1      |       0      |      1      |   L_1(u) = (u.n)(1,0,0)    |
      |---------|--------------|--------------|--------------|-------------|----------------------------|
      |    2    |       2      |       2      |       0      |      1      |   L_2(u) = (u.n)(0,1,0)    |
      |---------|--------------|--------------|--------------|-------------|----------------------------|
      |    3    |       2      |       3      |       0      |      1      |   L_3(u) = (u.n)(-1,0,0)   |
      |---------|--------------|--------------|--------------|-------------|----------------------------|
      |    4    |       2      |       4      |       0      |      1      |   L_4(u) = (u.n)(0,0,-1)   |
      |---------|--------------|--------------|--------------|-------------|----------------------------|
      |    5    |       2      |       5      |       0      |      1      |   L_5(u) = (u.n)(0,0,1)    |
      |=========|==============|==============|==============|=============|============================|
      |   MAX   |  maxScDim=2  |  maxScOrd=5  |  maxDfOrd=0  |      -      |                            |
      |=========|==============|==============|==============|=============|============================|
      \endverbatim

      \remarks
      \li     In the DoF functional \f${\bf n}\f$ is a face normal. Direction of face normals
      is determined by the right-hand rule applied to faces oriented by their vertex order
      in the cell topology, from face vertex 0 to last face vertex, whereas their length is
      set equal to face area (see http://mathworld.wolfram.com/Right-HandRule.html for definition
      of right-hand rule). For example, face 1 of all Hexahedron cells has vertex order
      {1,2,6,5} and its right-hand rule normal can be computed, e.g., by the vector product of
      edge tangents to edges {1,2} and {2,6}. On the reference Hexahedron the coordinates of
      face 1 vertices are (1,-1,-1), (1,1,-1), (1,1,1) and (1,-1,1), the edge tangents are
      (0,2,0) and (0,0,2) and the face normal direction is (0,2,0) X (0,0,2) = (4,0,0). In this
      case the normal length already equals face area and no further normalization is needed.

      \li     The length of the face normal equals the face area. As a result, the DoF functional
      is the value of the normal component of a vector field at the face center times the
      face area. The resulting basis is equivalent to a basis defined by using the face
      flux as a DoF functional. Note that the faces of reference Hexahedron<> cells all
      have the same area equal to 4.

  */

  namespace Impl {

    /**
      \brief See Intrepid2::Basis_HDIV_HEX_I1_FEM
    */
    class Basis_HDIV_HEX_I1_FEM {
    public:
      typedef struct Hexahedron<8> cell_topology_type;
      /**
        \brief See Intrepid2::Basis_HDIV_HEX_I1_FEM
      */
      template<EOperator opType>
      struct Serial {
        template<typename OutputViewType,
                 typename inputViewType>
        KOKKOS_INLINE_FUNCTION
        static void
        getValues(       OutputViewType output,
                   const inputViewType  input );
      };

      template<typename DeviceType,
               typename outputValueValueType, class ...outputValueProperties,
               typename inputPointValueType,  class ...inputPointProperties>
      static void
      getValues(       Kokkos::DynRankView<outputValueValueType,outputValueProperties...> outputValues,
                 const Kokkos::DynRankView<inputPointValueType, inputPointProperties...>  inputPoints,
                 const EOperator operatorType);

      /**
        \brief See Intrepid2::Basis_HDIV_HEX_I1_FEM
      */
      template<typename outputValueViewType,
               typename inputPointViewType,
               EOperator opType>
      struct Functor {
              outputValueViewType _outputValues;
        const inputPointViewType  _inputPoints;

        KOKKOS_INLINE_FUNCTION
        Functor( outputValueViewType outputValues_,
                 inputPointViewType  inputPoints_ )
          : _outputValues(outputValues_), _inputPoints(inputPoints_) {}

        KOKKOS_INLINE_FUNCTION
        void operator()(const ordinal_type pt) const {
          switch (opType) {
          case OPERATOR_VALUE : {
            auto       output = Kokkos::subview( _outputValues, Kokkos::ALL(), pt, Kokkos::ALL() );
            const auto input  = Kokkos::subview( _inputPoints,                 pt, Kokkos::ALL() );
            Serial<opType>::getValues( output, input );
            break;
          }
          case OPERATOR_DIV :{
            auto       output = Kokkos::subview( _outputValues, Kokkos::ALL(), pt );
            const auto input  = Kokkos::subview( _inputPoints,                 pt, Kokkos::ALL() );
            Serial<opType>::getValues( output, input );
            break;
          }
          default: {
            INTREPID2_TEST_FOR_ABORT( opType != OPERATOR_VALUE &&
                                      opType != OPERATOR_DIV,
                                      ">>> ERROR: (Intrepid2::Basis_HDIV_HEX_I1_FEM::Serial::getValues) operator is not supported");
          }
          }
        }
      };
    };

  }

  template<typename DeviceType = void,
           typename outputValueType = double,
           typename pointValueType = double>
  class Basis_HDIV_HEX_I1_FEM : public Basis<DeviceType,outputValueType,pointValueType> {
  public:
    using OrdinalTypeArray1DHost = typename Basis<DeviceType,outputValueType,pointValueType>::OrdinalTypeArray1DHost;
    using OrdinalTypeArray2DHost = typename Basis<DeviceType,outputValueType,pointValueType>::OrdinalTypeArray2DHost;
    using OrdinalTypeArray3DHost = typename Basis<DeviceType,outputValueType,pointValueType>::OrdinalTypeArray3DHost;

    template<EOperator opType>
    struct Serial {
      template<typename outputValueValueType, class ...outputValueProperties,
               typename inputPointValueType,  class ...inputPointProperties>
      KOKKOS_INLINE_FUNCTION
      static void
      getValues(       Kokkos::DynRankView<outputValueValueType,outputValueProperties...> outputValues,
                 const Kokkos::DynRankView<inputPointValueType, inputPointProperties...>  inputPoints );

    };


    /** \brief  Constructor.
     */
    Basis_HDIV_HEX_I1_FEM();

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
      Intrepid2::getValues_HDIV_Args(outputValues,
                                     inputPoints,
                                     operatorType,
                                     this->getBaseCellTopology(),
                                     this->getCardinality() );
#endif
      Impl::Basis_HDIV_HEX_I1_FEM::
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
                                    ">>> ERROR: (Intrepid2::Basis_HDIV_HEX_I1_FEM::getDofCoords) rank = 2 required for dofCoords array");
      // Verify 0th dimension of output array.
      INTREPID2_TEST_FOR_EXCEPTION( static_cast<ordinal_type>(dofCoords.extent(0)) != this->basisCardinality_, std::invalid_argument,
                                    ">>> ERROR: (Intrepid2::Basis_HDIV_HEX_I1_FEM::getDofCoords) mismatch in number of dof and 0th dimension of dofCoords array");
      // Verify 1st dimension of output array.
      INTREPID2_TEST_FOR_EXCEPTION( dofCoords.extent(1) != this->getBaseCellTopology().getDimension(), std::invalid_argument,
                                    ">>> ERROR: (Intrepid2::Basis_HDIV_HEX_I1_FEM::getDofCoords) incorrect reference cell (1st) dimension in dofCoords array");
#endif
      Kokkos::deep_copy(dofCoords, this->dofCoords_);
    }

  virtual
  void
  getDofCoeffs( ScalarViewType dofCoeffs ) const override {
#ifdef HAVE_INTREPID2_DEBUG
    // Verify rank of output array.
    INTREPID2_TEST_FOR_EXCEPTION( dofCoeffs.rank() != 2, std::invalid_argument,
        ">>> ERROR: (Intrepid2::Basis_HDIV_HEX_I1_FEM::getDofCoeffs) rank = 2 required for dofCoeffs array");
    // Verify 0th dimension of output array.
    INTREPID2_TEST_FOR_EXCEPTION( static_cast<ordinal_type>(dofCoeffs.extent(0)) != this->getCardinality(), std::invalid_argument,
        ">>> ERROR: (Intrepid2::Basis_HDIV_HEX_I1_FEM::getDofCoeffs) mismatch in number of dof and 0th dimension of dofCoeffs array");
    // Verify 1st dimension of output array.
    INTREPID2_TEST_FOR_EXCEPTION( dofCoeffs.extent(1) != this->getBaseCellTopology().getDimension(), std::invalid_argument,
        ">>> ERROR: (Intrepid2::Basis_HDIV_HEX_I1_FEM::getDofCoeffs) incorrect reference cell (1st) dimension in dofCoeffs array");
#endif
    Kokkos::deep_copy(dofCoeffs, this->dofCoeffs_);
  }

    virtual
    const char*
    getName() const override {
      return "Intrepid2_HDIV_HEX_I1_FEM";
    }

    virtual
    bool
    requireOrientation() const override {
      return true;
    }

    /** \brief returns the basis associated to a subCell.

        The bases of the subCell are the restriction to the subCell of the bases of the parent cell,
        projected along normal to the subCell.

        \param [in] subCellDim - dimension of subCell
        \param [in] subCellOrd - position of the subCell among of the subCells having the same dimension
        \return pointer to the subCell basis of dimension subCellDim and position subCellOrd
     */
    BasisPtr<DeviceType,outputValueType,pointValueType>
    getSubCellRefBasis(const ordinal_type subCellDim, const ordinal_type subCellOrd) const override{

      if(subCellDim == 2) {
        return Teuchos::rcp(new
          Basis_HVOL_C0_FEM<DeviceType,outputValueType,pointValueType>(shards::getCellTopologyData<shards::Quadrilateral<4> >()));
      }
      INTREPID2_TEST_FOR_EXCEPTION(true,std::invalid_argument,"Input parameters out of bounds");
    }
    BasisPtr<typename Kokkos::HostSpace::device_type,outputValueType,pointValueType>
    getHostBasis() const override{
      return Teuchos::rcp(new Basis_HDIV_HEX_I1_FEM<typename Kokkos::HostSpace::device_type,outputValueType,pointValueType>());
    }
  };

}// namespace Intrepid2

#include "Intrepid2_HDIV_HEX_I1_FEMDef.hpp"

#endif
