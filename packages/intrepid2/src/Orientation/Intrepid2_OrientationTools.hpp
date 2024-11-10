// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_OrientationTools.hpp
    \brief  Header file for the Intrepid2::OrientationTools and Intrepid2::Impl::OrientationTools classes.
    \author Created by Kyungjoo Kim
*/

#ifndef __INTREPID2_ORIENTATIONTOOLS_HPP__
#define __INTREPID2_ORIENTATIONTOOLS_HPP__

#include "Intrepid2_ConfigDefs.hpp"
#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"

#include "Shards_CellTopology.hpp"
#include "Shards_BasicTopologies.hpp"

#include "Intrepid2_PointTools.hpp"

#include "Intrepid2_Basis.hpp"

// -- HGRAD family
#include "Intrepid2_HGRAD_LINE_Cn_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_Cn_FEM.hpp"
#include "Intrepid2_HGRAD_HEX_Cn_FEM.hpp"

#include "Intrepid2_HGRAD_TRI_Cn_FEM.hpp"
#include "Intrepid2_HGRAD_TET_Cn_FEM.hpp"

// -- HCURL family
#include "Intrepid2_HCURL_QUAD_In_FEM.hpp"
#include "Intrepid2_HCURL_HEX_In_FEM.hpp"

#include "Intrepid2_HCURL_TRI_In_FEM.hpp"
#include "Intrepid2_HCURL_TET_In_FEM.hpp"
#include "Intrepid2_HVOL_LINE_Cn_FEM.hpp"

// -- HDIV family
#include "Intrepid2_HDIV_QUAD_In_FEM.hpp"
#include "Intrepid2_HDIV_HEX_In_FEM.hpp"

#include "Intrepid2_HDIV_TRI_In_FEM.hpp"
#include "Intrepid2_HDIV_TET_In_FEM.hpp"
#include "Intrepid2_HVOL_TRI_Cn_FEM.hpp"

// -- Lower order family
#include "Intrepid2_HCURL_QUAD_I1_FEM.hpp"
#include "Intrepid2_HCURL_TRI_I1_FEM.hpp"

#include "Intrepid2_HDIV_QUAD_I1_FEM.hpp"
#include "Intrepid2_HDIV_TRI_I1_FEM.hpp"

#include "Intrepid2_HCURL_HEX_I1_FEM.hpp"
#include "Intrepid2_HCURL_TET_I1_FEM.hpp"
#include "Intrepid2_HCURL_WEDGE_I1_FEM.hpp"

#include "Intrepid2_HDIV_HEX_I1_FEM.hpp"
#include "Intrepid2_HDIV_TET_I1_FEM.hpp"
#include "Intrepid2_HDIV_WEDGE_I1_FEM.hpp"

#include "Intrepid2_HierarchicalBasisFamily.hpp"
#include "Intrepid2_NodalBasisFamily.hpp"

#include "Teuchos_LAPACK.hpp"


namespace Intrepid2 {

  namespace Impl {

    /**
       \brief Tools to compute orientations for degrees-of-freedom
    */ 
    class OrientationTools {
    public:

      // -----------------------------------------------------------------------------
      // Point modification
      //
      //

      /** \brief  Computes modified point for line segment.
          
          \param  ot       [out] - modified point value
          \param  pt       [in]  - input point in [-1.0 , 1.0]
          \param  ort      [in]  - orientation number between 0 and 1
      */
      template<typename ValueType>
      KOKKOS_INLINE_FUNCTION
      static void 
      getModifiedLinePoint(ValueType &ot,
                           const ValueType pt,
                           const ordinal_type ort);
      
      /** \brief  Computes modified point for triangle.
          
          \param  ot0      [out] - modified coordinate 0
          \param  ot1      [out] - modified coordinate 1
          \param  pt0      [out] - input coordinate 0
          \param  pt1      [out] - input coordinate 1
          \param  ort      [in]  - orientation number between 0 and 5
      */
      template<typename ValueType>
      KOKKOS_INLINE_FUNCTION
      static void 
      getModifiedTrianglePoint(ValueType &ot0,
                               ValueType &ot1,
                               const ValueType pt0,
                               const ValueType pt1,
                               const ordinal_type ort);
      
      /** \brief  Computes modified point for quadrilateral.
          
          \param  ot0      [out] - modified coordinate 0
          \param  ot1      [out] - modified coordinate 1
          \param  pt0      [out] - input coordinate 0
          \param  pt1      [out] - input coordinate 1
          \param  ort      [in]  - orientation number between 0 and 7
      */
      template<typename ValueType>
      KOKKOS_INLINE_FUNCTION
      static void 
      getModifiedQuadrilateralPoint(ValueType &ot0,
                                    ValueType &ot1,
                                    const ValueType pt0,
                                    const ValueType pt1,
                                    const ordinal_type ort);

      /** \brief  Computes modified parameterization maps of 1- and 2-subcells with orientation.
          
          \param  outPoints       [out] - rank-2 (P,D2) array with points in 1D or 2D modified domain with orientation
          \param  refPoints       [in]  - rank-2 (P,D2) array with points in 1D or 2D parameter domain
          \param  cellTopo        [in]  - cell topology of the parameterized domain (1- and 2-subcells)
          \param  cellOrt         [in]  - cell orientation number (zero is aligned with shards default configuration
      */
      template<typename outPointViewType,
               typename refPointViewType>
      inline
      static void 
      mapToModifiedReference(outPointViewType outPoints,
                             const refPointViewType refPoints,
                             const shards::CellTopology cellTopo,
                             const ordinal_type cellOrt = 0);

      /** \brief  Computes modified parameterization maps of 1- and 2-subcells with orientation.

          \param  outPoints       [out] - rank-2 (P,D2) array with points in 1D or 2D modified domain with orientation
          \param  refPoints       [in]  - rank-2 (P,D2) array with points in 1D or 2D parameter domain
          \param  cellTopoKey     [in]  - key of the cell topology of the parameterized domain (1- and 2-subcells)
          \param  cellOrt         [in]  - cell orientation number (zero is aligned with shards default configuration
      */
      template<typename outPointViewType,
               typename refPointViewType>
      KOKKOS_INLINE_FUNCTION
      static void
      mapToModifiedReference(outPointViewType outPoints,
                             const refPointViewType refPoints,
                             const unsigned cellTopoKey,
                             const ordinal_type cellOrt = 0);


      /** \brief  Computes Jacobian of orientation map for line segment.

          \param  jacobian    [out] - rank-2 (D,D) array with jacobian of the orientation map
          \param  ort         [in]  - orientation number between 0 and 1
      */
      template<typename JacobianViewType>
      KOKKOS_INLINE_FUNCTION
      static void
      getLineJacobian(JacobianViewType jacobian, const ordinal_type ort);

      /** \brief  Computes Jacobian of orientation map for triangle.

          \param  jacobian    [out] - rank-2 (D,D) array with jacobian of the orientation map
          \param  ort         [in]  - orientation number between 0 and 5
      */
      template<typename JacobianViewType>
      KOKKOS_INLINE_FUNCTION
      static void
      getTriangleJacobian(JacobianViewType jacobian, const ordinal_type ort);

      /** \brief  Computes Jacobian of orientation map for quadrilateral.

          \param  jacobian    [out] - rank-2 (D,D) array with jacobian of the orientation map
          \param  ort         [in]  - orientation number between 0 and 7
      */
      template<typename JacobianViewType>
      KOKKOS_INLINE_FUNCTION
      static void
      getQuadrilateralJacobian(JacobianViewType jacobian, const ordinal_type ort);


      /** \brief  Computes jacobian of the parameterization maps of 1- and 2-subcells with orientation.

          \param  jacobian        [out] - rank-2 (D,D) array with jacobian of the orientation map
          \param  cellTopo        [in]  - cell topology of the parameterized domain (1- and 2-subcells)
          \param  cellOrt         [in]  - cell orientation number (zero is aligned with shards default configuration
      */
      template<typename JacobianViewType>
      inline
      static void
      getJacobianOfOrientationMap(JacobianViewType jacobian,
                                  const shards::CellTopology cellTopo,
                                  const ordinal_type cellOrt);

      /** \brief  Computes jacobian of the parameterization maps of 1- and 2-subcells with orientation.

          \param  jacobian        [out] - rank-2 (D,D) array with jacobian of the orientation map
          \param  cellTopoKey     [in]  - key of the cell topology of the parameterized domain (1- and 2-subcells)
          \param  cellOrt         [in]  - cell orientation number (zero is aligned with shards default configuration
      */
      template<typename JacobianViewType>
      KOKKOS_INLINE_FUNCTION
      static void
      getJacobianOfOrientationMap(JacobianViewType jacobian,
                                  const unsigned cellTopoKey,
                                  const ordinal_type cellOrt);


      /** \brief  Computes the (oriented) subCell tangents
          \param  tangents        [out] - rank-2 (scD,D), tangents of the subcell. scD: subCell dimension, D: parent cell dimension
          \param  subcellParam    [in]  - rank-2 (N, scD+1), subCells parametrization. N:number of subcells, scD: subCell dimension
          \param  cellTopoKey     [in]  - key of the cell topology of the parameterized domain (1- and 2-subcells)
          \param  subCellOrd      [in]  - position of the subCell among subCells of a given dimension
          \param  ort             [in]  - subCell orientation number
      */
      template<typename TanViewType, typename ParamViewType>
      KOKKOS_INLINE_FUNCTION
      static void getRefSubcellTangents(TanViewType tangents,
                                        const ParamViewType subCellParametrization,
                                        const unsigned subcellTopoKey,
                                        const ordinal_type subCellOrd,
                                        const ordinal_type ort);


      /** \brief  Computes the (oriented) tangents and normal of the side subCell
          The normals are only defined for sides (subCells of dimension D-1) and it
          requires the computation of tangents

          \param  tangentsAndNormal [out] - rank-2 (D,D), (D-1) tangents and 1 normal. D: parent cell dimension
          \param  subcellParam      [in]  - rank-2 (N, scD+1), subCells parametrization. N:number of subcells, scD: subCell dimension
          \param  cellTopoKey       [in]  - key of the cell topology of the parameterized domain (1- and 2-subcells)
          \param  subCellOrd        [in]  - position of the subCell among subCells of a given dimension
          \param  ort               [in]  - subCell orientation number
      */
      template<typename TanNormViewType, typename ParamViewType>
      KOKKOS_INLINE_FUNCTION
      static void getRefSideTangentsAndNormal(TanNormViewType tangentsAndNormal,
                                              const ParamViewType subCellParametrization,
                                              const unsigned subcellTopoKey,
                                              const ordinal_type subCellOrd,
                                              const ordinal_type ort);


      /** \brief  Maps points defined on the subCell manifold into the parent Cell accounting for orientation

          \param  cellCoords        [out] - rank-2 (P, D), P points mapped in the parent cell manifold.
          \param  subCellCoords     [in]  - rank-2 (P, scD), P points defined on the subCell manifold, scD: subCell dimension
          \param  subcellParam      [in]  - rank-2 (N, scD+1), subCells parametrization. N:number of subCells, scD: subCell dimension
          \param  cellTopoKey       [in]  - key of the cell topology of the parameterized domain (1- and 2-subCells)
          \param  subCellOrd        [in]  - position of the subCell among subCells of a given dimension
          \param  ort               [in]  - subCell orientation number
      */
      template<typename coordsViewType, typename subcellCoordsViewType, typename ParamViewType>
      KOKKOS_INLINE_FUNCTION
      static void mapSubcellCoordsToRefCell(coordsViewType cellCoords,
                                            const subcellCoordsViewType subCellCoords,
                                            const ParamViewType subcellParametrization,
                                            const unsigned subcellTopoKey,
                                            const ordinal_type subCellOrd,
                                            const ordinal_type ort);

      // -----------------------------------------------------------------------------
      // Coefficient Matrix
      //
      //
      
      /** \brief  Compute orientation matrix for HGRAD basis for a given subcell
          and its reference basis
          
          \param  output      [out]  - rank 2 coefficient matrix
          \param  subcellBasis [in]  - reference subcell basis function
          \param  cellBasis    [in]  - cell basis function
          \param  subcellId    [in]  - subcell Id in the cell topology
          \param  subcellOrt   [in]  - orientation number between 0 and 1
          \param  inverse      [in]  - boolean, when true the inverse of the orintation matrix is computed

      */
      template<typename OutputViewType,
               typename subcellBasisHostType,
               typename cellBasisHostType>
      inline
      static void
      getCoeffMatrix_HGRAD(OutputViewType &output,
                           const subcellBasisHostType& subcellBasis,
                           const cellBasisHostType& cellBasis,
                           const ordinal_type subcellId,
                           const ordinal_type subcellOrt,
                           const bool inverse = false);

      /** \brief  Compute orientation matrix for HCURL basis for a given subcell
          and its reference basis
          
          \param  output      [out]  - rank 2 coefficient matrix
          \param  subcellBasis [in]  - reference subcell basis function
          \param  cellBasis    [in]  - cell basis function
          \param  subcellId    [in]  - subcell Id in the cell topology
          \param  subcellOrt   [in]  - orientation number between 0 and 1
          \param  inverse      [in]  - boolean, when true the inverse of the orintation matrix is computed
      */
      template<typename OutputViewType,
               typename subcellBasisHostType,
               typename cellBasisHostType>
      inline
      static void
      getCoeffMatrix_HCURL(OutputViewType &output,
                           const subcellBasisHostType& subcellBasis,
                           const cellBasisHostType& cellBasis,
                           const ordinal_type subcellId,
                           const ordinal_type subcellOrt,
                           const bool inverse = false);


      /** \brief  Compute orientation matrix for HDIV basis for a given subcell
          and its reference basis
          
          \param  output      [out]  - rank 2 orientation matrix
          \param  subcellBasis [in]  - reference subcell basis
          \param  cellBasis    [in]  - cell basis function
          \param  subcellId    [in]  - subcell Id in the cell topology
          \param  subcellOrt   [in]  - orientation number between 0 and 1
          \param  inverse      [in]  - boolean, when true the inverse of the orintation matrix is computed
      */
      template<typename OutputViewType,
               typename subcellBasisHostType,
               typename cellBasisHostType>
      inline
      static void
      getCoeffMatrix_HDIV(OutputViewType &output,
                          const subcellBasisHostType& subcellBasis,
                          const cellBasisHostType& cellBasis,
                          const ordinal_type subcellId,
                          const ordinal_type subcellOrt,
                           const bool inverse = false);


      /** \brief  Compute orientation matrix for HVOL basis for a given (2D or 1D) cell
          and its reference basis. This method is used only for side basis.
          
          \param  output      [out]  - rank 2 coefficient matrix
          \param  cellBasis    [in]  - cell basis function
          \param  subcellOrt   [in]  - orientation number between 0 and 1
          \param  inverse      [in]  - boolean, when true the inverse of the orintation matrix is computed
      */
      template<typename OutputViewType,
               typename cellBasisHostType>
      inline
      static void
      getCoeffMatrix_HVOL(OutputViewType &output,
                           const cellBasisHostType& cellBasis,
                           const ordinal_type cellOrt,
                           const bool inverse = false);

    };
  }

  /**
     \brief Tools to compute orientations for degrees-of-freedom
  */ 
  template<typename DeviceType>
  class OrientationTools {
  public:

    /** \brief  subcell ordinal, orientation, matrix m x n
     */
    typedef Kokkos::View<double****,DeviceType> CoeffMatrixDataViewType;

    // 
    /** \brief  key :: basis name, order, value :: matrix data view type 
     */
    static std::map<std::pair<std::string,ordinal_type>,CoeffMatrixDataViewType> ortCoeffData;
    static std::map<std::pair<std::string,ordinal_type>,CoeffMatrixDataViewType> ortInvCoeffData;
    
  private:

    template<typename BasisHostType>
    inline 
    static CoeffMatrixDataViewType createCoeffMatrixInternal(const BasisHostType* basis, const bool invTrans = false);
    

    /** \brief  Compute orientation matrix for HGRAD basis
     */
    template<typename BasisHostType>
    inline
    static void init_HGRAD(CoeffMatrixDataViewType matData,
                           BasisHostType const *cellBasis,
                           const bool inverse = false);

    /** \brief  Compute orientation matrix for HCURL basis
     */
    template<typename BasisHostType>
    inline
    static void init_HCURL(CoeffMatrixDataViewType matData,
                           BasisHostType const *cellBasis,
                           const bool inverse = false);

    /** \brief  Compute orientation matrix for HDIV basis
     */
    template<typename BasisHostType>
    inline
    static void init_HDIV(CoeffMatrixDataViewType matData,
                          BasisHostType const *cellBasis,
                           const bool inverse = false);
    
    /** \brief  Compute orientation matrix for HVOL basis
     */
    template<typename BasisHostType>
    inline
    static void init_HVOL(CoeffMatrixDataViewType matData,
                           BasisHostType const *cellBasis,
                           const bool inverse = false);

  public:

    /** \brief  Create coefficient matrix.
        \param  basis      [in]  - basis type
    */
    template<typename BasisType>
    inline 
    static CoeffMatrixDataViewType createCoeffMatrix(const BasisType* basis);

    /** \brief  Create inverse of coefficient matrix.
        \param  basis      [in]  - basis type
    */
    template<typename BasisType>
    inline 
    static CoeffMatrixDataViewType createInvCoeffMatrix(const BasisType* basis);

    /** \brief  Clear coefficient matrix
     */
    inline 
    static void clearCoeffMatrix();

    /** \brief  Compute orientations of cells in a workset
        \param  elemOrts      [out]  - cell orientations
        \param  elemNodes      [in]  - node coordinates
        \param  cellTopo       [in]  - shards cell topology
        \param  isSide         [in]  - boolean, whether the cell is a side
    */
    template<typename elemOrtValueType, class ...elemOrtProperties,
             typename elemNodeValueType, class ...elemNodeProperties>
    inline 
    static void
    getOrientation(Kokkos::DynRankView<elemOrtValueType,elemOrtProperties...> elemOrts,
                   const Kokkos::DynRankView<elemNodeValueType,elemNodeProperties...> elemNodes,
                   const shards::CellTopology cellTopo,
                   bool isSide = false);


    /** \brief  Modify basis due to orientation
        \param  output        [out]  - output array, of shape (C,F,P[,D])
        \param  input          [in]  - input array, of shape (C,F,P[,D]) or (F,P[,D])
        \param  orts           [in]  - orientations, of shape (C)
        \param  basis          [in]  - basis of cardinality F
        \param  transpose      [in]  - boolean, when true the transpose of the orintation matrix is applied
    */
    template<typename outputValueType, class ...outputProperties,
             typename inputValueType,  class ...inputProperties,
             typename OrientationViewType,
             typename BasisType>
    inline
    static void
    modifyBasisByOrientation(Kokkos::DynRankView<outputValueType,outputProperties...> output,
                             const Kokkos::DynRankView<inputValueType, inputProperties...>  input,
                             const OrientationViewType orts,
                             const BasisType * basis,
                             const bool transpose = false);
    
    /** \brief  Modify basis due to orientation, applying the transpose of the operator applied in modifyBasisByOrientation().  If the input provided represents basis coefficents in the global orientation, then this method will appropriately transform them to the local orientation.
        \param  output        [out]  - output array, of shape (C,F,P[,D])
        \param  input          [in]  - input array, of shape (C,F,P[,D]) or (F,P[,D])
        \param  orts           [in]  - orientations, of shape (C)
        \param  basis          [in]  - basis of cardinality F
    */
    template<typename outputValueType, class ...outputProperties,
             typename inputValueType,  class ...inputProperties,
             typename OrientationViewType,
             typename BasisType>
    inline
    static void
    modifyBasisByOrientationTranspose(Kokkos::DynRankView<outputValueType,outputProperties...> output,
                                      const Kokkos::DynRankView<inputValueType, inputProperties...>  input,
                                      const OrientationViewType orts,
                                      const BasisType * basis);


    /** \brief  Modify basis due to orientation, applying the inverse of the operator applied in modifyBasisByOrientation().
        \param  output        [out]  - output array, of shape (C,F,P[,D])
        \param  input          [in]  - input array, of shape (C,F,P[,D]) or (F,P[,D])
        \param  orts           [in]  - orientations, of shape (C)
        \param  basis          [in]  - basis of cardinality F
        \param  transpose      [in]  - boolean, when true the transpose of the inverse of the operatore applied
    */
    template<typename outputValueType, class ...outputProperties,
             typename inputValueType,  class ...inputProperties,
             typename OrientationViewType,
             typename BasisType>
    inline
    static void
    modifyBasisByOrientationInverse(Kokkos::DynRankView<outputValueType,outputProperties...> output,
                                      const Kokkos::DynRankView<inputValueType, inputProperties...>  input,
                                      const OrientationViewType orts,
                                      const BasisType * basis,
                                      const bool transpose = false);
    
    
    /** \brief  Modify an assembled (C,F1,F2) matrix according to orientation of the cells.
        \param  output           [out]  - output array, shape (C,F1,F2)
        \param  input             [in]  - input array, shape (C,F1,F2)
        \param  orts               [in]  - orientations, shape (C)
        \param  basisLeft    [in]  - basis with cardinality F1
        \param  basisRight  [in]  - basis with cardinality F2
    */
    template<typename outputValueType, class ...outputProperties,
             typename inputValueType,  class ...inputProperties,
             typename OrientationViewType,
             typename BasisTypeLeft,
             typename BasisTypeRight>
    inline
    static void
    modifyMatrixByOrientation(Kokkos::DynRankView<outputValueType,outputProperties...> output,
                              const Kokkos::DynRankView<inputValueType, inputProperties...>  input,
                              const OrientationViewType orts,
                              const BasisTypeLeft* basisLeft,
                              const BasisTypeRight* basisRight);
  };
  
  template<typename T> 
  std::map<std::pair<std::string,ordinal_type>, typename OrientationTools<T>::CoeffMatrixDataViewType>
  OrientationTools<T>::ortCoeffData;

  template<typename T> 
  std::map<std::pair<std::string,ordinal_type>, typename OrientationTools<T>::CoeffMatrixDataViewType>
  OrientationTools<T>::ortInvCoeffData;
}

// include templated function definitions
#include "Intrepid2_OrientationToolsDefModifyPoints.hpp"
#include "Intrepid2_OrientationToolsDefCoeffMatrix_HGRAD.hpp"
#include "Intrepid2_OrientationToolsDefCoeffMatrix_HCURL.hpp"
#include "Intrepid2_OrientationToolsDefCoeffMatrix_HDIV.hpp"
#include "Intrepid2_OrientationToolsDefCoeffMatrix_HVOL.hpp"
#include "Intrepid2_OrientationToolsDefMatrixData.hpp"
#include "Intrepid2_OrientationToolsDefModifyBasis.hpp"

#endif
