// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_HCURL_TRI_I1_FEM.hpp
    \brief  Header file for the Intrepid2::Basis_HCURL_TRI_I1_FEM class.
    \author Created by P. Bochev, D. Ridzal and K. Peterson.
            Kokkorized by Kyungjoo Kim
*/

#ifndef INTREPID2_HCURL_TRI_I1_FEM_HPP
#define INTREPID2_HCURL_TRI_I1_FEM_HPP

#include "Intrepid2_Basis.hpp"
#include "Intrepid2_HVOL_C0_FEM.hpp"

namespace Intrepid2 {

  /** \class  Intrepid2::Basis_HCURL_TRI_I1_FEM
      \brief  Implementation of the default H(curl)-compatible FEM basis of degree 1 on Triangle cell

      Implements Nedelec basis of the first kind of degree 1 on the reference Triangle cell.
      The basis has cardinality 3 and spans an INCOMPLETE linear polynomial space. Basis functions
      are dual to a unisolvent set of degrees-of-freedom (DoF) defined and enumerated as follows:

      \verbatim
      ===================================================================================================
      |         |           degree-of-freedom-tag table                    |                            |
      |   DoF   |----------------------------------------------------------|       DoF definition       |
      | ordinal |  subc dim    | subc ordinal | subc DoF ord |subc num DoF |                            |
      |=========|==============|==============|==============|=============|============================|
      |    0    |       1      |       0      |       0      |      1      |  L_0(u) = (u.t)(0.5, 0)    |
      |---------|--------------|--------------|--------------|-------------|----------------------------|
      |    1    |       1      |       1      |       0      |      1      |  L_1(u) = (u.t)(0.5,0.5)   |
      |---------|--------------|--------------|--------------|-------------|----------------------------|
      |    2    |       1      |       2      |       0      |      1      |  L_2(u) = (u.t)( 0,0.5)    |
      |=========|==============|==============|==============|=============|============================|
      |   MAX   |  maxScDim=1  |  maxScOrd=2  |  maxDfOrd=0  |      -      |                            |
      |=========|==============|==============|==============|=============|============================|
      \endverbatim

      \remarks
      \li     In the DoF functional \f${\bf t}\f$ is an edge tangent. Direction of edge
      tangents follows the vertex order of the edges in the cell topology and runs from
      edge vertex 0 to edge vertex 1, whereas their length is set equal to the edge length. For
      example, edge 1 of all Triangle reference cells has vertex order {1,2}, i.e., its tangent
      runs from vertex 1 of the reference Triangle to vertex 2 of that cell. On the reference
      Triangle the coordinates of these vertices are (1,0) and (0,1), respectively. Therefore,
      the tangent to edge 1 is (0,1)-(1,0) = (-1,1). Because its length already equals edge length,
      no further rescaling of the edge tangent is needed.

      \li     The length of the edge tangent equals the edge length. As a result, the DoF functional
      is the value of the tangent component of a vector field at the edge midpoint times the
      edge length. The resulting basis is equivalent to a basis defined by using the edge circulation as a DoF functional. Note that edges
      0 and 2 of reference Triangle<> cells have unit lengths and edge 1 has length Sqrt(2).

  */

  namespace Impl {

    /**
      \brief See Intrepid2::Basis_HCURL_TRI_I1_FEM
    */
    class Basis_HCURL_TRI_I1_FEM {
    public:
      typedef struct Triangle<3> cell_topology_type;
      /**
       \brief See Intrepid2::Basis_HCURL_TRI_I1_FEM
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
       \brief See Intrepid2::Basis_HCURL_TRI_I1_FEM
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
            auto       output = Kokkos::subview( _outputValues, Kokkos::ALL(), pt, Kokkos::ALL() );
            const auto input  = Kokkos::subview( _inputPoints,                 pt, Kokkos::ALL() );
            Serial<opType>::getValues( output, input );
            break;
          }
          case OPERATOR_CURL : {
            auto       output = Kokkos::subview( _outputValues, Kokkos::ALL(), pt );
            const auto input  = Kokkos::subview( _inputPoints,                 pt, Kokkos::ALL() );
            Serial<opType>::getValues( output, input );
            break;
          }
          case OPERATOR_DIV: {
            INTREPID2_TEST_FOR_ABORT( (opType == OPERATOR_DIV),
                                      ">>> ERROR (Basis_HCURL_TRI_I1_FEM): DIV is invalid operator for HCURL Basis Functions");
            break;
          }
          case OPERATOR_GRAD: {
            INTREPID2_TEST_FOR_ABORT( (opType == OPERATOR_GRAD),
                                      ">>> ERROR (Basis_HCURL_TRI_I1_FEM): GRAD is invalid operator for HCURL Basis Functions");
            break;
          }
          default: {
            INTREPID2_TEST_FOR_ABORT( opType != OPERATOR_VALUE &&
                                      opType != OPERATOR_CURL,
                                      ">>> ERROR: (Intrepid2::Basis_HCURL_QUAD_I1_FEM::Serial::getValues) operator is not supported");
          }
          }
        }
      };

    };
  }

  template< typename DeviceType = void,
            typename outputValueType = double,
            typename pointValueType = double >
  class Basis_HCURL_TRI_I1_FEM : public Basis<DeviceType, outputValueType, pointValueType> {
  public:
    using BasisBase = Basis<DeviceType,outputValueType,pointValueType>;
    using typename BasisBase::ExecutionSpace;

    using typename BasisBase::OrdinalTypeArray1DHost;
    using typename BasisBase::OrdinalTypeArray2DHost;
    using typename BasisBase::OrdinalTypeArray3DHost;

    using typename BasisBase::OutputViewType;
    using typename BasisBase::PointViewType ;
    using typename BasisBase::ScalarViewType;

    /** \brief  Constructor.
     */
    Basis_HCURL_TRI_I1_FEM();

    using BasisBase::getValues;

    virtual
    void
    getValues( const ExecutionSpace& space,
                     OutputViewType outputValues,
               const PointViewType  inputPoints,
               const EOperator operatorType = OPERATOR_VALUE ) const override {
#ifdef HAVE_INTREPID2_DEBUG
      // Verify arguments
      Intrepid2::getValues_HCURL_Args(outputValues,
                                      inputPoints,
                                      operatorType,
                                      this->getBaseCellTopology(),
                                      this->getCardinality() );
#endif
      Impl::Basis_HCURL_TRI_I1_FEM::
        getValues<DeviceType>(space,
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
                                    ">>> ERROR: (Intrepid2::Basis_HCURL_TRI_I1_FEM::getDofCoords) rank = 2 required for dofCoords array");
      // Verify 0th dimension of output array.
      INTREPID2_TEST_FOR_EXCEPTION( static_cast<ordinal_type>(dofCoords.extent(0)) != this->basisCardinality_, std::invalid_argument,
                                    ">>> ERROR: (Intrepid2::Basis_HCURL_TRI_I1_FEM::getDofCoords) mismatch in number of dof and 0th dimension of dofCoords array");
      // Verify 1st dimension of output array.
      INTREPID2_TEST_FOR_EXCEPTION( dofCoords.extent(1) != this->getBaseCellTopology().getDimension(), std::invalid_argument,
                                    ">>> ERROR: (Intrepid2::Basis_HCURL_TRI_I1_FEM::getDofCoords) incorrect reference cell (1st) dimension in dofCoords array");
#endif
      Kokkos::deep_copy(dofCoords, this->dofCoords_);
    }

  virtual
  void
  getDofCoeffs( ScalarViewType dofCoeffs ) const override {
#ifdef HAVE_INTREPID2_DEBUG
    // Verify rank of output array.
    INTREPID2_TEST_FOR_EXCEPTION( rank(dofCoeffs) != 2, std::invalid_argument,
        ">>> ERROR: (Intrepid2::Basis_HCURL_TRI_I1_FEM::getDofCoeffs) rank = 2 required for dofCoeffs array");
    // Verify 0th dimension of output array.
    INTREPID2_TEST_FOR_EXCEPTION( static_cast<ordinal_type>(dofCoeffs.extent(0)) != this->getCardinality(), std::invalid_argument,
        ">>> ERROR: (Intrepid2::Basis_HCURL_TRI_I1_FEM::getDofCoeffs) mismatch in number of dof and 0th dimension of dofCoeffs array");
    // Verify 1st dimension of output array.
    INTREPID2_TEST_FOR_EXCEPTION( dofCoeffs.extent(1) != this->getBaseCellTopology().getDimension(), std::invalid_argument,
        ">>> ERROR: (Intrepid2::Basis_HCURL_TRI_I1_FEM::getDofCoeffs) incorrect reference cell (1st) dimension in dofCoeffs array");
#endif
    Kokkos::deep_copy(dofCoeffs, this->dofCoeffs_);
  }


    virtual
    const char*
    getName() const override {
      return "Intrepid2_HCURL_TRI_I1_FEM";
    }

    virtual
    bool
    requireOrientation() const override {
      return true;
    }

    /** \brief returns the basis associated to a subCell.

        The bases of the subCell are the restriction to the subCell of the bases of the parent cell,
        projected to the subCell line.

        \param [in] subCellDim - dimension of subCell
        \param [in] subCellOrd - position of the subCell among of the subCells having the same dimension
        \return pointer to the subCell basis of dimension subCellDim and position subCellOrd
     */
    BasisPtr<DeviceType,outputValueType,pointValueType>
    getSubCellRefBasis(const ordinal_type subCellDim, const ordinal_type subCellOrd) const override{
      if(subCellDim == 1)
        return Teuchos::rcp( new
          Basis_HVOL_C0_FEM<DeviceType,outputValueType,pointValueType>(shards::getCellTopologyData<shards::Line<2> >()));

      INTREPID2_TEST_FOR_EXCEPTION(true,std::invalid_argument,"Input parameters out of bounds");
    }

    BasisPtr<typename Kokkos::HostSpace::device_type,outputValueType,pointValueType>
    getHostBasis() const override{
      return Teuchos::rcp(new Basis_HCURL_TRI_I1_FEM<typename Kokkos::HostSpace::device_type,outputValueType,pointValueType>());
    }
  };

}// namespace Intrepid2

#include "Intrepid2_HCURL_TRI_I1_FEMDef.hpp"

#endif
