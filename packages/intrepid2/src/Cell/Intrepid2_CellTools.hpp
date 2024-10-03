// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_CellTools.hpp
    \brief  Header file for the Intrepid2::CellTools class.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_CELLTOOLS_HPP__
#define __INTREPID2_CELLTOOLS_HPP__

#include "Intrepid2_ConfigDefs.hpp"

#include "Shards_CellTopology.hpp"
#include "Shards_BasicTopologies.hpp"

#include "Teuchos_RCP.hpp" 

#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"

#include "Intrepid2_RealSpaceTools.hpp"

#include "Intrepid2_Basis.hpp"

#include "Intrepid2_HGRAD_LINE_C1_FEM.hpp"
#include "Intrepid2_HGRAD_LINE_C2_FEM.hpp"

#include "Intrepid2_HGRAD_TRI_C1_FEM.hpp"
#include "Intrepid2_HGRAD_TRI_C2_FEM.hpp"

#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_C2_FEM.hpp"

#include "Intrepid2_HGRAD_TET_C1_FEM.hpp"
#include "Intrepid2_HGRAD_TET_C2_FEM.hpp"
#include "Intrepid2_HGRAD_TET_COMP12_FEM.hpp"

#include "Intrepid2_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid2_HGRAD_HEX_C2_FEM.hpp"

#include "Intrepid2_HGRAD_WEDGE_C1_FEM.hpp"
#include "Intrepid2_HGRAD_WEDGE_C2_FEM.hpp"

#include "Intrepid2_HGRAD_PYR_C1_FEM.hpp"
#include "Intrepid2_HGRAD_PYR_I2_FEM.hpp"

#include "Intrepid2_Data.hpp"
#include "Intrepid2_CellData.hpp"

namespace Intrepid2 {

  //============================================================================================//
  //                                                                                            //
  //                                          CellTools                                         //
  //                                                                                            //
  //============================================================================================//

  /** \class  Intrepid2::CellTools
      \brief  A stateless class for operations on cell data. Provides methods for:
      \li     computing Jacobians of reference-to-physical frame mappings, their inverses and determinants
      \li     application of the reference-to-physical frame mapping and its inverse
      \li     parametrizations of edges and faces of reference cells needed for edge and face integrals,
      \li     computation of edge and face tangents and face normals on both reference and physical frames
      \li     inclusion tests for point sets in reference and physical cells.
  */


  template<typename DeviceType>
  class CellTools {
    using ExecSpaceType = typename DeviceType::execution_space;
    using MemSpaceType = typename DeviceType::memory_space;
  public:

    /** \brief  Checks if a cell topology has reference cell
        \param  cell              [in]  - cell topology
        \return true if the cell topology has reference cell, false otherwise
    */
    inline
    static bool
    hasReferenceCell( const shards::CellTopology cellTopo ) {
      return RefSubcellParametrization<DeviceType>::isSupported(cellTopo.getKey());
    }

  private:

    /** \brief Generates default HGrad basis based on cell topology 
        \param cellTopo            [in] - cell topology 
     */
    template<typename outputValueType, 
             typename pointValueType>
    static Teuchos::RCP<Basis<DeviceType,outputValueType,pointValueType> >
    createHGradBasis( const shards::CellTopology cellTopo ) {
      Teuchos::RCP<Basis<DeviceType,outputValueType,pointValueType> > r_val;

      switch (cellTopo.getKey()) {
      case shards::Line<2>::key:          r_val = Teuchos::rcp(new Basis_HGRAD_LINE_C1_FEM   <DeviceType,outputValueType,pointValueType>()); break;
      case shards::Line<3>::key:          r_val = Teuchos::rcp(new Basis_HGRAD_LINE_C2_FEM   <DeviceType,outputValueType,pointValueType>()); break;
      case shards::Triangle<3>::key:      r_val = Teuchos::rcp(new Basis_HGRAD_TRI_C1_FEM    <DeviceType,outputValueType,pointValueType>()); break;
      case shards::Quadrilateral<4>::key: r_val = Teuchos::rcp(new Basis_HGRAD_QUAD_C1_FEM   <DeviceType,outputValueType,pointValueType>()); break;
      case shards::Tetrahedron<4>::key:   r_val = Teuchos::rcp(new Basis_HGRAD_TET_C1_FEM    <DeviceType,outputValueType,pointValueType>()); break;
      case shards::Hexahedron<8>::key:    r_val = Teuchos::rcp(new Basis_HGRAD_HEX_C1_FEM    <DeviceType,outputValueType,pointValueType>()); break;
      case shards::Wedge<6>::key:         r_val = Teuchos::rcp(new Basis_HGRAD_WEDGE_C1_FEM  <DeviceType,outputValueType,pointValueType>()); break;
      case shards::Pyramid<5>::key:       r_val = Teuchos::rcp(new Basis_HGRAD_PYR_C1_FEM    <DeviceType,outputValueType,pointValueType>()); break;

      case shards::Triangle<6>::key:      r_val = Teuchos::rcp(new Basis_HGRAD_TRI_C2_FEM    <DeviceType,outputValueType,pointValueType>()); break;
      case shards::Quadrilateral<8>::key: r_val = Teuchos::rcp(new Basis_HGRAD_QUAD_I2_FEM   <DeviceType,outputValueType,pointValueType>()); break;
      case shards::Quadrilateral<9>::key: r_val = Teuchos::rcp(new Basis_HGRAD_QUAD_C2_FEM   <DeviceType,outputValueType,pointValueType>()); break;
      case shards::Tetrahedron<10>::key:  r_val = Teuchos::rcp(new Basis_HGRAD_TET_C2_FEM    <DeviceType,outputValueType,pointValueType>()); break;
      case shards::Tetrahedron<11>::key:  r_val = Teuchos::rcp(new Basis_HGRAD_TET_COMP12_FEM<DeviceType,outputValueType,pointValueType>()); break;
      case shards::Hexahedron<20>::key:   r_val = Teuchos::rcp(new Basis_HGRAD_HEX_I2_FEM    <DeviceType,outputValueType,pointValueType>()); break;
      case shards::Hexahedron<27>::key:   r_val = Teuchos::rcp(new Basis_HGRAD_HEX_C2_FEM    <DeviceType,outputValueType,pointValueType>()); break;
      case shards::Wedge<15>::key:        r_val = Teuchos::rcp(new Basis_HGRAD_WEDGE_I2_FEM  <DeviceType,outputValueType,pointValueType>()); break;
      case shards::Wedge<18>::key:        r_val = Teuchos::rcp(new Basis_HGRAD_WEDGE_C2_FEM  <DeviceType,outputValueType,pointValueType>()); break;
      case shards::Pyramid<13>::key:      r_val = Teuchos::rcp(new Basis_HGRAD_PYR_I2_FEM    <DeviceType,outputValueType,pointValueType>()); break;

      case shards::Beam<2>::key:
      case shards::Beam<3>::key:
      case shards::ShellLine<2>::key:
      case shards::ShellLine<3>::key:
      case shards::ShellTriangle<3>::key:
      case shards::ShellTriangle<6>::key:
      case shards::ShellQuadrilateral<4>::key:
      case shards::ShellQuadrilateral<8>::key:
      case shards::ShellQuadrilateral<9>::key: 
      default: {
        INTREPID2_TEST_FOR_EXCEPTION( true, std::invalid_argument,
                                      ">>> ERROR (Intrepid2::CellTools::createHGradBasis): Cell topology not supported.");
      }
      }
      return r_val;
    }

public:

   /** \brief  Default constructor.
      */
    CellTools() = default;

    /** \brief  Destructor
      */
    ~CellTools() = default;

    //============================================================================================//
    //                                                                                            //
    //                     Jacobian, inverse Jacobian and Jacobian determinant                    //
    //                                                                                            //
    //============================================================================================//

    /** \brief  Computes the Jacobian matrix \e DF of the reference-to-physical frame map \e F.

                There are two use cases:
        \li     Computes Jacobians \f$DF_{c}\f$ of the reference-to-physical map for \b all physical cells
                in a cell workset on a \b single set of reference points stored in a rank-2 (P,D) array;
        \li     Computes Jacobians \f$DF_{c}\f$ of the reference-to-physical map for \b all physical cells
                in a cell workset on \b multiple reference point sets having the same number of points,
                indexed by cell ordinal, and stored in a rank-3 (C,P,D) array;

                For a single point set in a rank-2 array (P,D) returns a rank-4 (C,P,D,D) array
                such that
        \f[
                \mbox{jacobian}(c,p,i,j) = [DF_{c}(\mbox{points}(p))]_{ij} \quad c=0,\ldots, C
        \f]
                For multiple sets of reference points in a rank-3 (C,P,D) array returns
                rank-4 (C,P,D,D) array such that
        \f[
                \mbox{jacobian}(c,p,i,j) = [DF_{c}(\mbox{points}(c,p))]_{ij} \quad c=0,\ldots, C
        \f]

                Requires pointer to HGrad basis that defines reference to physical cell mapping.  
                See Section \ref sec_cell_topology_ref_map_DF for definition of the Jacobian.

                \warning
                The points are not required to be in the reference cell associated with the specified
                cell topology. CellTools provides several inclusion tests methods to check whether
                or not the points are inside a reference cell.

        \param  jacobian          [out] - rank-4 array with dimensions (C,P,D,D) with the Jacobians
        \param  points               [in]  - rank-2/3 array with dimensions (P,D)/(C,P,D) with the evaluation points
        \param  cellWorkset    [in]  - rank-3 container with logical dimensions (C,N,D) with the nodes of the cell workset
        \param  basis                 [in]  - HGrad basis for reference to physical cell mapping
        \param  startCell        [in] - first cell index in cellWorkset for which we should compute the Jacobian; corresponds to the 0 index in Jacobian and/or points container.  Default: 0.
        \param  endCell            [in] - first cell index in cellWorkset that we do not process; endCell - startCell must equal the extent of the Jacobian container in dimension 0.  Default: -1, a special value that indicates the extent of the cellWorkset should be used.
     */
    template<typename JacobianViewType,
             typename PointViewType,
             typename WorksetType,
             typename HGradBasisType>
    static void
    setJacobian(       JacobianViewType             jacobian,
                 const PointViewType                points,
                 const WorksetType                  worksetCell,
                 const Teuchos::RCP<HGradBasisType> basis,
                 const int startCell=0, const int endCell=-1);
    
    /** \brief  Computes the Jacobian matrix \e DF of the reference-to-physical frame map \e F.

                There are two use cases:
        \li     Computes Jacobians \f$DF_{c}\f$ of the reference-to-physical map for \b all physical cells
                in a cell workset on a \b single set of reference points stored in a rank-2 (P,D) array;
        \li     Computes Jacobians \f$DF_{c}\f$ of the reference-to-physical map for \b all physical cells
                in a cell workset on \b multiple reference point sets having the same number of points,
                indexed by cell ordinal, and stored in a rank-3 (C,P,D) array;

                For a single point set in a rank-2 array (P,D) returns a rank-4 (C,P,D,D) array
                such that
        \f[
                \mbox{jacobian}(c,p,i,j) = [DF_{c}(\mbox{points}(p))]_{ij} \quad c=0,\ldots, C
        \f]
                For multiple sets of reference points in a rank-3 (C,P,D) array returns
                rank-4 (C,P,D,D) array such that
        \f[
                \mbox{jacobian}(c,p,i,j) = [DF_{c}(\mbox{points}(c,p))]_{ij} \quad c=0,\ldots, C
        \f]

                Requires pointer to HGrad basis that defines reference to physical cell mapping.
                See Section \ref sec_cell_topology_ref_map_DF for definition of the Jacobian.

                \warning
                The points are not required to be in the reference cell associated with the specified
                cell topology. CellTools provides several inclusion tests methods to check whether
                or not the points are inside a reference cell.

        \param  jacobian          [out] - rank-4 array with dimensions (C,P,D,D) with the Jacobians
        \param  cellWorkset    [in]  - rank-3 container with logical dimensions (C,N,D) with the nodes of the cell workset
        \param  gradients         [in]  - rank-3/4 array with dimensions (N,P,D)/(C,N,P,D) with the gradients of the physical-to-cell mapping
        \param  startCell        [in] - first cell index in cellWorkset for which we should compute the Jacobian; corresponds to the 0 index in Jacobian and/or points container.  Default: 0.
        \param  endCell            [in] - first cell index in cellWorkset that we do not process; endCell - startCell must equal the extent of the Jacobian container in dimension 0.  Default: -1, a special value that indicates the extent of the cellWorkset should be used.
     */
    template<typename JacobianViewType,
             typename BasisGradientsType,
             typename WorksetType>
    static void
    setJacobian(       JacobianViewType   jacobian,
                 const WorksetType        worksetCell,
                 const BasisGradientsType gradients,
                 const int startCell=0, const int endCell=-1);
    
    /** \brief  Computes the Jacobian matrix \e DF of the reference-to-physical frame map \e F.

                There are two use cases:
        \li     Computes Jacobians \f$DF_{c}\f$ of the reference-to-physical map for \b all physical cells
                in a cell workset on a \b single set of reference points stored in a rank-2 (P,D) array;
        \li     Computes Jacobians \f$DF_{c}\f$ of the reference-to-physical map for \b all physical cells
                in a cell workset on \b multiple reference point sets having the same number of points,
                indexed by cell ordinal, and stored in a rank-3 (C,P,D) array;

                For a single point set in a rank-2 array (P,D) returns a rank-4 (C,P,D,D) array
                such that
        \f[
                \mbox{jacobian}(c,p,i,j) = [DF_{c}(\mbox{points}(p))]_{ij} \quad c=0,\ldots, C
        \f]
                For multiple sets of reference points in a rank-3 (C,P,D) array returns
                rank-4 (C,P,D,D) array such that
        \f[
                \mbox{jacobian}(c,p,i,j) = [DF_{c}(\mbox{points}(c,p))]_{ij} \quad c=0,\ldots, C
        \f]

                Requires cell topology with a reference cell. See Section \ref sec_cell_topology_ref_map_DF
                for definition of the Jacobian.

                \warning
                The points are not required to be in the reference cell associated with the specified
                cell topology. CellTools provides several inclusion tests methods to check whether
                or not the points are inside a reference cell.

        \param  jacobian          [out] - rank-4 array with dimensions (C,P,D,D) with the Jacobians
        \param  points            [in]  - rank-2/3 array with dimensions (P,D)/(C,P,D) with the evaluation points
        \param  cellWorkset       [in]  - rank-3 array with dimensions (C,N,D) with the nodes of the cell workset
        \param  cellTopo          [in]  - cell topology of the cells stored in \c cellWorkset
     */

    template<typename JacobianViewType,
             typename PointViewType,
             typename WorksetCellViewType>
    static void
    setJacobian(       JacobianViewType    jacobian,
                 const PointViewType       points,
                 const WorksetCellViewType worksetCell,
                 const shards::CellTopology cellTopo ) {
    using nonConstPointValueType = typename PointViewType::non_const_value_type;
    auto basis = createHGradBasis<nonConstPointValueType,nonConstPointValueType>(cellTopo);
    setJacobian(jacobian, 
                points, 
                worksetCell, 
                basis);
   }

    /** \brief  Computes the inverse of the Jacobian matrix \e DF of the reference-to-physical frame map \e F.

        Returns rank-4 array with dimensions (C,P,D,D) such that
        \f[
        \mbox{jacobianInv}(c,p,*,*) = \mbox{jacobian}^{-1}(c,p,*,*) \quad c = 0,\ldots, C
        \f]

        \param  jacobianInv       [out] - rank-4 array with dimensions (C,P,D,D) with the inverse Jacobians
        \param  jacobian          [in]  - rank-4 array with dimensions (C,P,D,D) with the Jacobians
    */
    template<typename JacobianInvViewType,
             typename JacobianViewType>
    static void
    setJacobianInv(       JacobianInvViewType jacobianInv,
                    const JacobianViewType    jacobian );

    /** \brief  Computes the determinant of the Jacobian matrix \e DF of the reference-to-physical frame map \e F.

        Returns rank-2 array with dimensions (C,P) such that
        \f[
        \mbox{jacobianDet}(c,p) = \mbox{det}(\mbox{jacobian}(c,p,*,*)) \quad c=0,\ldots, C
        \f]

        \param  jacobianDet       [out] - rank-2 array with dimensions (C,P) with Jacobian determinants
        \param  jacobian          [in]  - rank-4 array with dimensions (C,P,D,D) with the Jacobians
    */
    template<typename JacobianDetViewType,
             typename JacobianViewType>
    static void
    setJacobianDet(       JacobianDetViewType jacobianDet,
                    const JacobianViewType    jacobian );

    /** \brief  Allocates and returns a Data container suitable for storing determinants corresponding to the Jacobians in the Data container provided

        \param  jacobian          [in]  - data with shape (C,P,D,D), as returned by CellGeometry::allocateJacobianData()
        \return Data container with shape (C,P)
    */
    template<class PointScalar>
    static Data<PointScalar,DeviceType> allocateJacobianDet( const Data<PointScalar,DeviceType> & jacobian );

    /** \brief  Allocates and returns a Data container suitable for storing inverses corresponding to the Jacobians in the Data container provided

        \param  jacobian          [in]  - data with shape (C,P,D,D), as returned by CellGeometry::allocateJacobianData()
        \return Data container with shape (C,P,D,D)
    */
    template<class PointScalar>
    static Data<PointScalar,DeviceType> allocateJacobianInv( const Data<PointScalar,DeviceType> & jacobian );

    /** \brief  Computes determinants corresponding to the Jacobians in the Data container provided

        \param  jacobianDet   [out]  - data with shape (C,P), as returned by CellTools::allocateJacobianDet()
        \param  jacobian          [in]    - data with shape (C,P,D,D), as returned by CellGeometry::allocateJacobianData()
    */
    template<class PointScalar>
    static void setJacobianDet( Data<PointScalar,DeviceType> & jacobianDet,
                               const Data<PointScalar,DeviceType> & jacobian);
    
    /** \brief  Computes reciprocals of determinants corresponding to the Jacobians in the Data container provided

        \param  jacobianDetInv   [out]  - data with shape (C,P), as returned by CellTools::allocateJacobianDet()
        \param  jacobian                 [in]    - data with shape (C,P,D,D), as returned by CellGeometry::allocateJacobianData()
    */
    template<class PointScalar>
    static void setJacobianDetInv( Data<PointScalar,DeviceType> & jacobianDetInv,
                                  const Data<PointScalar,DeviceType> & jacobian);

    /** \brief  Computes determinants corresponding to the Jacobians in the Data container provided

        \param  jacobianInv   [out]  - data container with shape (C,P,D,D), as returned by CellTools::allocateJacobianInv()
        \param  jacobian          [in]    - data with shape (C,P,D,D), as returned by CellGeometry::allocateJacobianData()
    */
    template<class PointScalar>
    static void setJacobianInv( Data<PointScalar,DeviceType> & jacobianInv,
                               const Data<PointScalar,DeviceType> & jacobian);
    
    /** \brief  Multiplies the Jacobian with shape (C,P,D,D) by the reciprocals of the determinants, with shape (C,P), entrywise.

        \param  jacobianDividedByDet   [out]  - data container with shape (C,P,D,D), as returned by CellTools::allocateJacobianInv()
        \param  jacobian   [in]  - data container with shape (C,P,D,D), as returned by CellTools::allocateJacobianInv()
        \param  jacobianDetInv          [in]    - data with shape (C,P,D,D), as returned by CellGeometry::allocateJacobianData()
    */
    template<class PointScalar>
    static void setJacobianDividedByDet( Data<PointScalar,DeviceType> & jacobianDividedByDet,
                                        const Data<PointScalar,DeviceType> & jacobian,
                                        const Data<PointScalar,DeviceType> & jacobianDetInv);
    
    //============================================================================================//
    //                                                                                            //
    //                     Node information                                                       //
    //                                                                                            //
    //============================================================================================//

    // the node information can be used inside of kokkos functor and needs kokkos inline and
    // exception should be an abort. for now, let's not decorate

    /** \brief  Computes the Cartesian coordinates of reference cell barycenter.

        Requires cell topology with a reference cell.

        \param  cellCenter        [out] - coordinates of the specified reference cell center
        \param  cell              [in]  - cell topology 
    */
    template<typename cellCenterValueType, class ...cellCenterProperties>
    static void
    getReferenceCellCenter( Kokkos::DynRankView<cellCenterValueType,cellCenterProperties...> cellCenter,
                            const shards::CellTopology cell );

    /** \brief  Retrieves the Cartesian coordinates of a reference cell vertex.

        Requires cell topology with a reference cell.

        \param  cellVertex        [out] - coordinates of the specified reference cell vertex
        \param  cell              [in]  - cell topology of the cell
        \param  vertexOrd         [in]  - vertex ordinal
    */
    template<typename cellVertexValueType, class ...cellVertexProperties>
    static void
    getReferenceVertex(       Kokkos::DynRankView<cellVertexValueType,cellVertexProperties...> cellVertex,
                        const shards::CellTopology cell,
                        const ordinal_type         vertexOrd );


    /** \brief  Retrieves the Cartesian coordinates of all vertices of a reference subcell.

        Returns rank-2 array with the Cartesian coordinates of the vertices of the
        specified reference cell subcell. Requires cell topology with a reference cell.

        \param  subcellVertices   [out] - rank-2 (V,D) array with the Cartesian coordinates of the reference subcell vertices
        \param  subcellDim        [in]  - dimension of the subcell; 0 <= subcellDim <= parentCell dimension
        \param  subcellOrd        [in]  - ordinal of the subcell
        \param  parentCell        [in]  - topology of the cell that owns the subcell

        \remark When \c subcellDim = dimension of the \c parentCell this method returns the Cartesian
        coordinates of the vertices of the reference cell itself.
        Note that this requires \e subcellOrd=0.
    */
    template<typename subcellVertexValueType, class ...subcellVertexProperties>
    static void
    getReferenceSubcellVertices(       Kokkos::DynRankView<subcellVertexValueType,subcellVertexProperties...> subcellVertices,
                                 const ordinal_type         subcellDim,
                                 const ordinal_type         subcellOrd,
                                 const shards::CellTopology parentCell );



    /** \brief  Retrieves the Cartesian coordinates of a reference cell node.

        Returns Cartesian coordinates of a reference cell node. Requires cell topology
        with a reference cell. Node coordinates are always returned as an (x,y,z)-triple
        regardlesss of the actual topological cell dimension. The unused coordinates are
        set to zero, e.g., node 0 of Line<2> is returned as {-1,0,0}.

        \remark
        Because the nodes of a cell with a base topology coincide with its vertices, for cells
        with base topology this method is equivalent to Intrepid2::CellTools::getReferenceVertex.

        \param  cellNode          [out] - coordinates of the specified reference vertex
        \param  cell              [in]  - cell topology of the cell
        \param  vertexOrd         [in]  - node ordinal
    */
    template<typename cellNodeValueType, class ...cellNodeProperties>
    static void
    getReferenceNode(       Kokkos::DynRankView<cellNodeValueType,cellNodeProperties...> cellNode,
                      const shards::CellTopology  cell,
                      const ordinal_type          nodeOrd );



    /** \brief  Retrieves the Cartesian coordinates of all nodes of a reference subcell.

        Returns rank-2 array with the Cartesian coordinates of the nodes of the
        specified reference cell subcell. Requires cell topology with a reference cell.

        \param  subcellNodes      [out] - rank-2 (N,D) array with the Cartesian coordinates of the reference subcell nodes
        \param  subcellDim        [in]  - dimension of the subcell; 0 <= subcellDim <= parentCell dimension
        \param  subcellOrd        [in]  - ordinal of the subcell
        \param  parentCell        [in]  - topology of the cell that owns the subcell

        \remark When \c subcellDim = dimension of the \c parentCell this method returns the Cartesian
        coordinates of the nodes of the reference cell itself. Note that this requires \c subcellOrd=0.
    */
    template<typename SubcellNodeViewType>
    static void
    getReferenceSubcellNodes(       SubcellNodeViewType  subcellNodes,
                              const ordinal_type         subcellDim,
                              const ordinal_type         subcellOrd,
                              const shards::CellTopology parentCell );

    /** \brief  Computes constant tangent vectors to edges of 2D or 3D reference cells.

        Returns rank-1 array with dimension (D), D=2 or D=3; such that
        \f[
        {refEdgeTangent}(*) = \hat{\bf t}_i = {\partial\hat{\Phi}_i(t)\over\partial t}\,,
        \f]
        where \f$\hat{\Phi}_i : R =[-1,1]\mapsto \hat{\mathcal E}_i\f$ is the parametrization map
        of the specified reference edge \f$\hat{\mathcal E}_i\f$, given by
        \f[
        \hat{\Phi}_i(t) = \left\{\begin{array}{ll}
        (\hat{x}(t),\hat{y}(t),\hat{z}(t)) & \mbox{for 3D parent cells} \\[1ex]
        (\hat{x}(t),\hat{y}(t))            & \mbox{for 2D parent cells} \\[1ex]
        \end{array}\right.
        \f]
        The length of computed edge tangents is one-half the length of their associated edges:
        \f[
        |\hat{\bf t}_i | = {1\over 2} |\hat{\mathcal E}_i |\,.
        \f]
        Because the edges of all reference cells are always affine images of [-1,1],
        the edge tangent is constant vector field.

        \param  refEdgeTangent    [out] - rank-1 array (D) with the edge tangent; D = cell dimension
        \param  edgeOrd           [in]  - ordinal of the edge whose tangent is computed
        \param  parentCell        [in]  - cell topology of the parent reference cell
    */
    template<typename RefEdgeTangentViewType>
    static void
    getReferenceEdgeTangent(       RefEdgeTangentViewType refEdgeTangent,
                             const ordinal_type           edgeOrd,
                             const shards::CellTopology   parentCell );

    /** \brief  Computes pairs of constant tangent vectors to faces of a 3D reference cells.

        Returns 2 rank-1 arrays with dimension (D), D=3, such that
        \f[
        {refFaceTanU}(*) = \hat{\bf t}_{i,u} = {\partial\hat{\Phi}_i(u,v)\over\partial u} =
        \left({\partial\hat{x}(u,v)\over\partial u},
        {\partial\hat{y}(u,v)\over\partial u},
        {\partial\hat{z}(u,v)\over\partial u} \right) ;
        \f]
        \f[
        {refFaceTanV}(*) = \hat{\bf t}_{i,v} = {\partial\hat{\Phi}_i(u,v)\over \partial v} =
        \left({\partial\hat{x}(u,v)\over\partial v},
        {\partial\hat{y}(u,v)\over\partial v},
        {\partial\hat{z}(u,v)\over\partial v} \right)\,;
        \f]
        where \f$\hat{\Phi}_i: R \mapsto \hat{\mathcal F}_i\f$
        is the parametrization map of the specified reference face \f$\hat{\mathcal F}_i\f$ given by
        \f[
        \hat{\Phi}_i(u,v) =(\hat{x}(u,v),\hat{y}(u,v),\hat{z}(u,v))
        \f]
        and
        \f[
        R = \left\{\begin{array}{rl}
        \{(0,0),(1,0),(0,1)\} & \mbox{if $\hat{\mathcal F}_i$  is Triangle} \\[1ex]
        [-1,1]\times [-1,1] & \mbox{if $\hat{\mathcal F}_i$ is Quadrilateral} \,.
        \end{array}\right.
        \f]
        Because the faces of all reference cells are always affine images of \e R ,
        the coordinate functions \f$\hat{x},\hat{y},\hat{z}\f$ of the parametrization map
        are linear and the face tangents are constant vectors.

        \param  refFaceTanU       [out] - rank-1 array (D) with (constant) tangent in u-direction
        \param  refFaceTanV       [out] - rank-1 array (D) with (constant) tangent in v-direction
        \param  faceOrd           [in]  - ordinal of the face whose tangents are computed
        \param  parentCell        [in]  - cell topology of the parent 3D reference cell
    */
    template<typename RefFaceTanViewType>
    static void
    getReferenceFaceTangents(       RefFaceTanViewType   refFaceTanU,
                                    RefFaceTanViewType   refFaceTanV,
                              const ordinal_type         faceOrd,
                              const shards::CellTopology parentCell );

    /** \brief  Computes constant normal vectors to sides of 2D or 3D reference cells.

        A side is defined as a subcell of dimension one less than that of its parent cell.
        Therefore, sides of 2D cells are 1-subcells (edges) and sides of 3D cells
        are 2-subcells (faces).

        Returns rank-1 array with dimension (D), D = 2 or 3 such that
        \f[
        {refSideNormal}(*) = \hat{\bf n}_i =
        \left\{\begin{array}{rl}
        \displaystyle
        \left({\partial\hat{\Phi}_i(t)\over\partial t}\right)^{\perp}
        & \mbox{for 2D parent cells} \\[2ex]
        \displaystyle
        {\partial\hat{\Phi}_{i}\over\partial u} \times
        {\partial\hat{\Phi}_{i}\over\partial v}   & \mbox{for 3D parent cells}
        \end{array}\right.
        \f]
        where \f$ (u_1,u_2)^\perp = (u_2, -u_1)\f$, and \f$\hat{\Phi}_i: R \mapsto \hat{\mathcal S}_i\f$
        is the parametrization map of the specified reference side \f$\hat{\mathcal S}_i\f$ given by
        \f[
        \hat{\Phi}_i(u,v) =
        \left\{\begin{array}{rl}
        (\hat{x}(t),\hat{y}(t))                   & \mbox{for 2D parent cells} \\[1ex]
        (\hat{x}(u,v),\hat{y}(u,v),\hat{z}(u,v))  & \mbox{for 3D parent cells}
        \end{array}\right.

        \f]
        For sides of 2D cells \e R=[-1,1] and for sides of 3D cells
        \f[
        R = \left\{\begin{array}{rl}
        \{(0,0),(1,0),(0,1)\}   & \mbox{if $\hat{\mathcal S}_i$ is Triangle} \\[1ex]
        [-1,1]\times [-1,1] & \mbox{if $\hat{\mathcal S}_i$ is Quadrilateral} \,.
        \end{array}\right.
        \f]
        For 3D cells the length of computed side normals is proportional to side area:
        \f[
        |\hat{\bf n}_i | = \left\{\begin{array}{rl}
        2 \mbox{Area}(\hat{\mathcal F}_i) & \mbox{if $\hat{\mathcal F}_i$  is Triangle} \\[1ex]
        \mbox{Area}(\hat{\mathcal F}_i) & \mbox{if $\hat{\mathcal F}_i$ is Quadrilateral} \,.
        \end{array}\right.
        \f]
        For 2D cells the length of computed side normals is proportional to side length:
        \f[
        |\hat{\bf n}_i | = {1\over 2} |\hat{\mathcal F}_i |\,.
        \f]
        Because the sides of all reference cells are always affine images of \e R ,
        the coordinate functions \f$\hat{x},\hat{y},\hat{z}\f$ of the parametrization maps
        are linear and the side normal is a constant vector.

        \remark
        - For 3D cells the reference side normal coincides with the face normal computed by
        Intrepid2::CellTools::getReferenceFaceNormal and these two methods are completely equivalent.
        - For 2D cells the reference side normal is defined by \f$\hat{{\bf n}}= \hat{\bf t}^\perp = (t_2,-t_1)\f$
        where \f$\hat{{\bf t}}=(t_1,t_2)\f$ is the tangent vector computed by
        Intrepid2::CellTools::getReferenceEdgeTangent. Therefore, the pair
        \f$(\hat{{\bf n}},\hat{{\bf t}})\f$ is positively oriented.

        \param  refSideNormal     [out] - rank-1 array (D) with (constant) side normal
        \param  sideOrd           [in]  - ordinal of the side whose normal is computed
        \param  parentCell        [in]  - cell topology of the parent reference cell
    */
    template<typename RefSideNormalViewType>
    static void
    getReferenceSideNormal(       RefSideNormalViewType refSideNormal,
                            const ordinal_type          sideOrd,
                            const shards::CellTopology  parentCell );

    /** \brief  Computes constant normal vectors to faces of 3D reference cell.

        Returns rank-1 array with dimension (D), D=3 such that
        \f[
        {refFaceNormal}(*) = \hat{\bf n}_i = {\partial\hat{\Phi}_{i}\over\partial u} \times
        {\partial\hat{\Phi}_{i}\over\partial v}
        \f]
        where \f$\hat{\Phi}_i: R \mapsto \hat{\mathcal F}_i\f$
        is the parametrization map of the specified reference face \f$\hat{\mathcal F}_i\f$ given by
        \f[
        \hat{\Phi}_i(u,v) =(\hat{x}(u,v),\hat{y}(u,v),\hat{z}(u,v))
        \f]
        and
        \f[
        R = \left\{\begin{array}{rl}
        \{(0,0),(1,0),(0,1)\} & \mbox{if ${\mathcal F}$  is Triangle} \\[1ex]
        [-1,1]\times [-1,1] & \mbox{if ${\mathcal F}$ is Quadrilateral} \,.
        \end{array}\right.
        \f]
        The length of computed face normals is proportional to face area:
        \f[
        |\hat{\bf n}_i | = \left\{\begin{array}{rl}
        2 \mbox{Area}(\hat{\mathcal F}_i) & \mbox{if $\hat{\mathcal F}_i$  is Triangle} \\[1ex]
        \mbox{Area}(\hat{\mathcal F}_i) & \mbox{if $\hat{\mathcal F}_i$ is Quadrilateral} \,.
        \end{array}\right.
        \f]
        Because the faces of all reference cells are always affine images of \e R ,
        the coordinate functions \f$\hat{x},\hat{y},\hat{z}\f$ of the parametrization map
        are linear and the face normal is a constant vector.

        \remark
        The method Intrepid2::CellTools::getReferenceFaceTangents computes the reference face tangents
        \f${\partial\hat{\Phi}_{i}/\partial u}\f$ and \f${\partial\hat{\Phi}_{i}/\partial v}\f$.

        \param  refFaceNormal     [out] - rank-1 array (D) with (constant) face normal
        \param  faceOrd           [in]  - ordinal of the face whose normal is computed
        \param  parentCell        [in]  - cell topology of the parent reference cell
    */
    template<typename RefFaceNormalViewType>
    static void
    getReferenceFaceNormal(       RefFaceNormalViewType refFaceNormal,
                            const ordinal_type          faceOrd,
                            const shards::CellTopology  parentCell );

    /** \brief  Computes non-normalized tangent vectors to physical edges in an edge workset
        \f$\{\mathcal{E}_{c,i}\}_{c=0}^{N}\f$; (see \ref sec_cell_topology_subcell_wset for definition of edge worksets).

        For every edge in the workset the tangents are computed at the points
        \f${\bf x}_p = F_c(\hat{\Phi}_i(t_p))\in\mathcal{E}_{c,i}\f$ that are images of points
        from <var>R=[-1,1]</var> on edge \f$\mathcal{E}_{c,i}\f$. Returns rank-3 array with
        dimensions (C,P,D1), D1=2 or D1=3 such that
        \f[
        {edgeTangents}(c,p,d) =
        DF_c(\hat{\Phi}_i(t_p))\cdot {\partial{\hat{\Phi}}_{i}(t_p)\over\partial t}\,; \qquad t_p \in R
        \f]
        In this formula:
        \li     \f$ DF_c \f$ is Jacobian of parent cell \f${\mathcal C}\f$ that owns physical edge \f${\mathcal E}_{c,i}\f$;
        \li     \f$ {\partial{\hat{\Phi}}_{i}/\partial t}\f$ is the (constant) tangent to reference edge
        \f$\hat{\mathcal E}_i\f$; see Intrepid2::CellTools::getReferenceEdgeTangent that has the
        same local ordinal as the edges in the workset;
        \li     \f$ \hat{\Phi}_i R\mapsto\hat{\mathcal E}_i \f$ is parametrization of reference edge \f$\hat{\mathcal E}_i\f$;

        \warning
        \c worksetJacobians must contain the values of \f$DF_c(\hat{\Phi}_i(t_p))\f$,
        where \f$ t_p \in R=[-1,1] \f$, i.e., Jacobians of the parent cells evaluated at points
        that are located on reference edge \f$\hat{\mathcal E}_i\f$ having the same local ordinal as
        the edges in the workset.

        \param  edgeTangents      [out] - rank-3 array (C,P,D1) with tangents on workset edges
        \param  worksetJacobians  [in]  - rank-4 array (C,P,D1,D1) with Jacobians evaluated at ref. edge points
        \param  worksetEdgeOrd    [in]  - edge ordinal, relative to ref. cell, of the edge workset
        \param  parentCell        [in]  - cell topology of the parent reference cell
    */
    template<typename edgeTangentValueType,     class ...edgeTangentProperties,
             typename worksetJacobianValueType, class ...worksetJacobianProperties>
    static void
    getPhysicalEdgeTangents(       Kokkos::DynRankView<edgeTangentValueType,edgeTangentProperties...>         edgeTangents,
                             const Kokkos::DynRankView<worksetJacobianValueType,worksetJacobianProperties...> worksetJacobians,
                             const ordinal_type         worksetEdgeOrd,
                             const shards::CellTopology parentCell );


    /** \brief  Computes non-normalized tangent vectors to physical edges in an edge workset
        \f$\{\mathcal{E}_{c,i}\}_{c=0}^{N}\f$; (see \ref sec_cell_topology_subcell_wset for definition of edge worksets).

        It is similar to the <var>CellTools::getPhysicalEdgeTangents</var> function above, with the difference that the edge ordinal can change from point to point,
        and it is provided by the rank-2 input array <var><b>worksetEdgeOrds</b></var>, indexed by (C,P).

        \param  edgeTangents      [out] - rank-3 array (C,P,D1) with tangents on workset edges
        \param  worksetJacobians  [in]  - rank-4 array (C,P,D1,D1) with Jacobians evaluated at ref. edge points
        \param  worksetEdgeOrds   [in]  - rank-2 array (C,P) with edge ordinals, relative to ref. cell, of the edge workset
        \param  parentCell        [in]  - cell topology of the parent reference cell
    */
    template<typename edgeTangentValueType,     class ...edgeTangentProperties,
             typename worksetJacobianValueType, class ...worksetJacobianProperties,
             typename edgeOrdValueType,         class ...edgeOrdProperties>
    static void
    getPhysicalEdgeTangents(       Kokkos::DynRankView<edgeTangentValueType,edgeTangentProperties...>         edgeTangents,
                             const Kokkos::DynRankView<worksetJacobianValueType,worksetJacobianProperties...> worksetJacobians,
                             const Kokkos::DynRankView<edgeOrdValueType,edgeOrdProperties...>                 worksetEdgeOrds,
                             const shards::CellTopology parentCell );

    /** \brief  Computes non-normalized tangent vector pairs to physical faces in a face workset
        \f$\{\mathcal{F}_{c,i}\}_{c=0}^{N}\f$; (see \ref sec_cell_topology_subcell_wset for definition of face worksets).

        For every face in the workset the tangents are computed at the points
        \f${\bf x}_p = F_c(\hat{\Phi}_i(u_p,v_p))\in\mathcal{F}_{c,i}\f$ that are images of points
        from the parametrization domain \e R  on face \f$\mathcal{F}_{c,i}\f$.
        Returns 2 rank-3 arrays with dimensions (C,P,D), D=3 such that
        \f[
        {faceTanU}(c,p,d) = DF_c(\hat{\Phi}_i(u_p, v_p))\cdot {\partial\hat{\Phi}_i\over\partial u};\qquad
        {faceTanV}(c,p,d) = DF_c(\hat{\Phi}_i(u_p, v_p))\cdot {\partial\hat{\Phi}_{i}\over\partial v}\,;
        \qquad (u_p, v_p) \in R \,.
        \f]
        In this formula:
        \li     \f$ DF_c \f$ is Jacobian of parent cell \f${\mathcal C}\f$ that owns physical face \f${\mathcal F}_{c,i}\f$;
        \li     \f$ {\partial\hat{\Phi}_i/\partial u}, {\partial\hat{\Phi}_i/\partial v}\f$ are the (constant)
        tangents on reference face \f$\hat{\mathcal F}_i\f$; see Intrepid2::CellTools::getReferenceFaceTangents;
        that has the same local ordinal as the faces in the workset;
        \li     \f$ \hat{\Phi}_i : R\mapsto \hat{\mathcal F}_i\f$ is parametrization of reference face \f$\hat{\mathcal F}_i\f$;
        \li     \e R  is the parametrization domain for reference face \f$\hat{\mathcal F}_i\f$:
        \f[
        R =
        \left\{\begin{array}{rl}
        \{(0,0),(1,0),(0,1)\} & \mbox{if $\hat{\mathcal F}_i$ is Triangle} \\[1ex]
        [-1,1]\times [-1,1] & \mbox{if $\hat{\mathcal F}_i$ is Quadrilateral}
        \end{array}\right.
        \f]

        \warning
        \c worksetJacobians must contain the values of \f$DF_c(\hat{\Phi}_i(u_p,v_p))\f$,
        where \f$(u_p,v_p)\in R\f$, i.e., Jacobians of the parent cells evaluated at points
        that are located on reference face \f$\hat{\mathcal F}_i\f$ having the same local ordinal as
        the faces in the workset.

        \param  faceTanU          [out] - rank-3 array (C,P,D), image of ref. face u-tangent at workset faces
        \param  faceTanV          [out] - rank-3 array (C,P,D), image of ref. face u-tangent at workset faces
        \param  worksetJacobians  [in]  - rank-4 array (C,P,D,D) with Jacobians at ref. face points
        \param  worksetFaceOrd    [in]  - face ordinal, relative to ref. cell, of the face workset
        \param  parentCell        [in]  - cell topology of the parent reference cell
    */
    template<typename faceTanValueType,        class ...faceTanProperties,
             typename worksetJacobianValueType, class ...worksetJacobianProperties>
    static void
    getPhysicalFaceTangents(       Kokkos::DynRankView<faceTanValueType,faceTanProperties...> faceTanU,
                                   Kokkos::DynRankView<faceTanValueType,faceTanProperties...> faceTanV,
                             const Kokkos::DynRankView<worksetJacobianValueType,worksetJacobianProperties...> worksetJacobians,
                             const ordinal_type         worksetFaceOrd,
                             const shards::CellTopology parentCell );


    /** \brief  Computes non-normalized tangent vector pairs to physical faces in a face workset
        \f$\{\mathcal{F}_{c,i}\}_{c=0}^{N}\f$; (see \ref sec_cell_topology_subcell_wset for definition of face worksets).

        It is similar to the <var>CellTools::getPhysicalFaceTangents</var> function above, with the difference that the face ordinal can change from point to point,
        and it is provided by the rank-2 input array <var><b>worksetFaceOrds</b></var>, indexed by (C,P).

        \param  faceTanU          [out] - rank-3 array (C,P,D), image of ref. face u-tangent at workset faces
        \param  faceTanV          [out] - rank-3 array (C,P,D), image of ref. face u-tangent at workset faces
        \param  worksetJacobians  [in]  - rank-4 array (C,P,D,D) with Jacobians at ref. face points
        \param  worksetFaceOrds   [in]  - rank-2 array (C,P) with face ordinals, relative to ref. cell, of the face workset
        \param  parentCell        [in]  - cell topology of the parent reference cell
    */
    template<typename faceTanValueType,        class ...faceTanProperties,
             typename worksetJacobianValueType, class ...worksetJacobianProperties,
             typename faceOrdValueType, class ...faceOrdProperties>
    static void
    getPhysicalFaceTangents(       Kokkos::DynRankView<faceTanValueType,faceTanProperties...> faceTanU,
                                   Kokkos::DynRankView<faceTanValueType,faceTanProperties...> faceTanV,
                             const Kokkos::DynRankView<worksetJacobianValueType,worksetJacobianProperties...> worksetJacobians,
                             const Kokkos::DynRankView<faceOrdValueType,faceOrdProperties...>  worksetFaceOrds,
                             const shards::CellTopology parentCell );



    /** \brief  Computes non-normalized normal vectors to physical sides in a side workset
        \f$\{\mathcal{S}_{c,i}\}_{c=0}^{N}\f$.

        For every side in the workset the normals are computed at the points
        \f${\bf x}_p = F_c(\hat{\Phi}_i(P_p))\in\mathcal{S}_{c,i}\f$ that are images of points
        from the parametrization domain \e R  on side \f$\mathcal{S}_{c,i}\f$.
        A side is defined as a subcell of dimension one less than that of its  parent cell.
        Therefore, sides of 2D cells are 1-subcells (edges) and sides of 3D cells are 2-subcells (faces).

        Returns rank-3 array with dimensions (C,P,D), D = 2 or 3, such that
        \f[
        {sideNormals}(c,p,d) =
        \left\{\begin{array}{crl}
        \displaystyle
        \left(DF_c(\hat{\Phi}_i(t_p))\cdot
        {\partial{\hat{\Phi}}_{i}(t_p)\over\partial t}\right)^{\perp} &  t_p\in R
        & \mbox{for 2D parent cells} \\[2ex]
        \displaystyle
        \left( DF_c(\hat{\Phi}_i(u_p, v_p))\cdot {\partial\hat{\Phi}_i\over\partial u}\right) \times
        \left( DF_c(\hat{\Phi}_i(u_p, v_p))\cdot {\partial\hat{\Phi}_i\over\partial v}\right) \,;
        & (u_p, v_p) \in R                          & \mbox{for 3D parent cells}
        \end{array}\right.
        \f]
        In this formula:
        \li     \f$ DF_c \f$ is Jacobian of parent cell \f${\mathcal C}\f$ that owns physical side \f${\mathcal S}_{c,i}\f$;
        \li     For 2D cells: \f$ {\partial{\hat{\Phi}}_{i}/\partial t}\f$ is the (constant) tangent to reference side (edge)
        \f$\hat{\mathcal S}_i\f$; see Intrepid2::CellTools::getReferenceEdgeTangent, that has the
        same local ordinal as the sides in the workset;
        \li     For 3D cells: \f$ {\partial\hat{\Phi}_i/\partial u}, {\partial\hat{\Phi}_i/\partial v}\f$ are the (constant)
        tangents on reference side (face) \f$\hat{\mathcal S}_i\f$; see Intrepid2::CellTools::getReferenceFaceTangents,
        that has the same local ordinal as the sides in the workset;
        \li     \f$ \hat{\Phi}_i : R\mapsto \hat{\mathcal S}_i\f$ is parametrization of reference side \f$\hat{\mathcal S}_i\f$;
        \li     \e R  is the parametrization domain for reference side \f$\hat{\mathcal S}_i\f$. For
        2D parent cells \e R=[-1,1] and for 3D parent cells
        \f[
        R = \left\{\begin{array}{rl}
        \{(0,0),(1,0),(0,1)\} & \mbox{if $\hat{\mathcal S}_i$ is Triangle} \\[1ex]
        [-1,1]\times [-1,1] & \mbox{if $\hat{\mathcal S}_i$ is Quadrilateral}
        \end{array}\right.
        \f]

        \remark
        - For 3D cells the physical side normals coincides with the face normals computed by
        Intrepid2::CellTools::getPhysicalFaceNormals and these two methods are completely equivalent.
        - For 2D cells the physical side normals are defined by \f${\bf n}=(t_2,-t_1)\f$
        where \f${\bf t}=(t_1,t_2)\f$ are the physical edge tangents computed by
        Intrepid2::CellTools::getPhysicalEdgeTangents. Therefore, the pairs \f$({\bf n},{\bf t})\f$ are positively oriented.


        \warning
        \c worksetJacobians must contain the values of \f$DF_c(\hat{\Phi}_i(P_p))\f$,
        where \f$P_p\in R\f$, i.e., Jacobians of the parent cells evaluated at points
        that are located on reference side \f$\hat{\mathcal S}_i\f$ having the same local ordinal as
        the sides in the workset.

        \param  sideNormals       [out] - rank-3 array (C,P,D), normals at workset sides
        \param  worksetJacobians  [in]  - rank-4 array (C,P,D,D) with Jacobians at ref. side points
        \param  worksetSideOrd    [in]  - side ordinal, relative to ref. cell, of the side workset
        \param  parentCell        [in]  - cell topology of the parent reference cell
    */
    template<typename sideNormalValueType,      class ...sideNormalProperties,
             typename worksetJacobianValueType, class ...worksetJacobianProperties>
    static void
    getPhysicalSideNormals(       Kokkos::DynRankView<sideNormalValueType,sideNormalProperties...> sideNormals,
                            const Kokkos::DynRankView<worksetJacobianValueType,worksetJacobianProperties...> worksetJacobians,
                            const ordinal_type         worksetSideOrd,
                            const shards::CellTopology parentCell );


/** \brief  Computes non-normalized normal vectors to physical sides in a side workset
    \f$\{\mathcal{S}_{c,i}\}_{c=0}^{N}\f$.
    It is similar to the <var>CellTools::getPhysicalSideNormals</var> function above, with the difference that the side ordinal can change from point to point,
    and it is provided by the rank-2 input array <var><b>worksetSideOrds</b></var>, indexed by (C,P).

    \param  sideNormals       [out] - rank-3 array (C,P,D), normals at workset sides
    \param  worksetJacobians  [in]  - rank-4 array (C,P,D,D) with Jacobians at ref. side points
    \param  worksetSideOrds   [in]  - rank-2 array (C,P) with side ordinals, relative to ref. cell, of the side workset
    \param  parentCell        [in]  - cell topology of the parent reference cell
*/
    template<typename sideNormalValueType,      class ...sideNormalProperties,
             typename worksetJacobianValueType, class ...worksetJacobianProperties,
             typename edgeOrdValueType,         class ...edgeOrdProperties>
    static void
    getPhysicalSideNormals(       Kokkos::DynRankView<sideNormalValueType,sideNormalProperties...> sideNormals,
                            const Kokkos::DynRankView<worksetJacobianValueType,worksetJacobianProperties...> worksetJacobians,
                            const Kokkos::DynRankView<edgeOrdValueType,edgeOrdProperties...>                 worksetSideOrds,
                            const shards::CellTopology parentCell );

    /** \brief  Computes non-normalized normal vectors to physical faces in a face workset
        \f$\{\mathcal{F}_{c,i}\}_{c=0}^{N}\f$; (see \ref sec_cell_topology_subcell_wset for definition of face worksets).

        For every face in the workset the normals are computed at the points
        \f${\bf x}_p = F_c(\hat{\Phi}_i(u_p,v_p))\in\mathcal{F}_{c,i}\f$ that are images of points
        from the parametrization domain \e R  on face \f$\mathcal{F}_{c,i}\f$.
        Returns rank-3 array with dimensions (C,P,D), D=3, such that
        \f[
        {faceNormals}(c,p,d) =
        \left( DF_c(\hat{\Phi}_i(u_p, v_p))\cdot {\partial\hat{\Phi}_i\over\partial u}\right) \times
        \left( DF_c(\hat{\Phi}_i(u_p, v_p))\cdot {\partial\hat{\Phi}_i\over\partial v}\right) \,;
        \qquad (u_p, v_p) \in R \,.
        \f]
        In this formula:
        \li     \f$ DF_c \f$ is Jacobian of parent cell \f${\mathcal C}\f$ that owns physical face \f${\mathcal F}_{c,i}\f$;
        \li     \f$ {\partial\hat{\Phi}_i/\partial u}, {\partial\hat{\Phi}_i/\partial v}\f$ are the (constant)
        tangents on reference face \f$\hat{\mathcal F}_i\f$; see Intrepid2::CellTools::getReferenceFaceTangents;
        that has the same local ordinal as the faces in the workset;
        \li     \f$ \hat{\Phi}_i : R\mapsto \hat{\mathcal F}_i\f$ is parametrization of reference face \f$\hat{\mathcal F}_i\f$;
        \li     \e R  is the parametrization domain for reference face \f$\hat{\mathcal F}_i\f$:
        \f[
        R = \left\{\begin{array}{rl}
        \{(0,0),(1,0),(0,1)\} & \mbox{if $\hat{\mathcal F}_i$ is Triangle} \\[1ex]
        [-1,1]\times [-1,1] & \mbox{if $\hat{\mathcal F}_i$ is Quadrilateral}
        \end{array}\right.
        \f]

        \warning
        \c worksetJacobians must contain the values of \f$DF_c(\hat{\Phi}_i(u_p,v_p))\f$,
        where \f$(u_p,v_p)\in R\f$, i.e., Jacobians of the parent cells evaluated at points
        that are located on reference face \f$\hat{\mathcal F}_i\f$ having the same local ordinal as
        the faces in the workset.

        \param  faceNormals       [out] - rank-3 array (C,P,D), normals at workset faces
        \param  worksetJacobians  [in]  - rank-4 array (C,P,D,D) with Jacobians at ref. face points
        \param  worksetFaceOrd    [in]  - face ordinal, relative to ref. cell, of the face workset
        \param  parentCell        [in]  - cell topology of the parent reference cell
    */
    template<typename faceNormalValueType,      class ...faceNormalProperties,
             typename worksetJacobianValueType, class ...worksetJacobianProperties>
    static void
    getPhysicalFaceNormals(       Kokkos::DynRankView<faceNormalValueType,faceNormalProperties...> faceNormals,
                            const Kokkos::DynRankView<worksetJacobianValueType,worksetJacobianProperties...> worksetJacobians,
                            const ordinal_type         worksetFaceOrd,
                            const shards::CellTopology parentCell );


    /** \brief  Computes non-normalized normal vectors to physical faces in a face workset
        \f$\{\mathcal{F}_{c,i}\}_{c=0}^{N}\f$; (see \ref sec_cell_topology_subcell_wset for definition of face worksets).
        It is similar to the <var>CellTools::getPhysicalSideNormals</var> function above, with the difference that the side ordinal can change from point to point,
        and it is provided by the rank-2 input array <var><b>worksetSideOrds</b></var>, indexed by (C,P).

        \param  faceNormals       [out] - rank-3 array (C,P,D), normals at workset faces
        \param  worksetJacobians  [in]  - rank-4 array (C,P,D,D) with Jacobians at ref. face points
        \param  worksetFaceOrds   [in]  - rank-2 array (C,P) with face ordinals, relative to ref. cell, of the face workset
        \param  parentCell        [in]  - cell topology of the parent reference cell
    */

    template<typename faceNormalValueType,      class ...faceNormalProperties,
             typename worksetJacobianValueType, class ...worksetJacobianProperties,
             typename faceOrdValueType, class ...faceOrdProperties>
    static void
    getPhysicalFaceNormals(       Kokkos::DynRankView<faceNormalValueType,faceNormalProperties...> faceNormals,
                            const Kokkos::DynRankView<worksetJacobianValueType,worksetJacobianProperties...> worksetJacobians,
                            const Kokkos::DynRankView<faceOrdValueType,faceOrdProperties...>  worksetFaceOrds,
                            const shards::CellTopology parentCell );

    //============================================================================================//
    //                                                                                            //
    //                      Reference-to-physical frame mapping and its inverse                   //
    //                                                                                            //
    //============================================================================================//

    /** \brief  Computes \e F, the reference-to-physical frame map.

        There are 2 use cases:
        \li     Applies \f$ F_{c} \f$ for \b all cells in a cell workset to a \b single point set stored
        in a rank-2 (P,D) array;
        \li     Applies \f$ F_{c} \f$ for \b all cells in a cell workset to \b multiple point sets having
        the same number of points, indexed by cell ordinal, and stored in a rank-3 (C,P,D) array;

        For a single point set in a rank-2 array (P,D) returns a rank-3 (C,P,D) array such that
        \f[
        \mbox{physPoints}(c,p,d) = \Big(F_c(\mbox{refPoint}(p,*)) \Big)_d \quad c=0,\ldots, C
        \f]
        For multiple point sets in a rank-3 (C,P,D) array returns a rank-3 (C,P,D) array such that
        \f[
        \mbox{physPoints}(c,p,d) = \Big(F_c(\mbox{refPoint}(c,p,*)) \Big)_d \quad c=0,\ldots, C
        \f]
        This corresponds to mapping multiple sets of reference points to a matching number of
        physical cells.

        Requires pointer to HGrad basis that defines reference to physical cell mapping.  
        See Section \ref sec_cell_topology_ref_map for definition of the mapping function. 

        \warning
        The array \c refPoints represents an arbitrary set of points in the reference
        frame that are not required to be in the reference cell corresponding to the
        specified cell topology. As a result, the images of these points under a given
        reference-to-physical map are not necessarily contained in the physical cell that
        is the image of the reference cell under that map. CellTools provides several
        inclusion tests methods to check whether or not the points are inside a reference cell.

        \param  physPoints        [out] - rank-3 array with dimensions (C,P,D) with the images of the ref. points
        \param  refPoints          [in]  - rank-3/2 array with dimensions (C,P,D)/(P,D) with points in reference frame
        \param  cellWorkset      [in]  - rank-3 container with logical dimensions (C,N,D) with the nodes of the cell workset
        \param  basis                   [in]  - pointer to HGrad basis used in reference-to-physical cell mapping

    */
    template<typename PhysPointValueType,
             typename RefPointValueType,
             typename WorksetType,
             typename HGradBasisPtrType>
    static void
    mapToPhysicalFrame(       PhysPointValueType physPoints,
                        const RefPointValueType  refPoints,
                        const WorksetType        worksetCell,
                        const HGradBasisPtrType  basis );

    /** \brief  Computes \e F, the reference-to-physical frame map.

        There are 2 use cases:
        \li     Applies \f$ F_{c} \f$ for \b all cells in a cell workset to a \b single point set stored
        in a rank-2 (P,D) array;
        \li     Applies \f$ F_{c} \f$ for \b all cells in a cell workset to \b multiple point sets having
        the same number of points, indexed by cell ordinal, and stored in a rank-3 (C,P,D) array;

        For a single point set in a rank-2 array (P,D) returns a rank-3 (C,P,D) array such that
        \f[
        \mbox{physPoints}(c,p,d) = \Big(F_c(\mbox{refPoint}(p,*)) \Big)_d \quad c=0,\ldots, C
        \f]
        For multiple point sets in a rank-3 (C,P,D) array returns a rank-3 (C,P,D) array such that
        \f[
        \mbox{physPoints}(c,p,d) = \Big(F_c(\mbox{refPoint}(c,p,*)) \Big)_d \quad c=0,\ldots, C
        \f]
        This corresponds to mapping multiple sets of reference points to a matching number of
        physical cells.

        Requires cell topology with a reference cell. See Section \ref sec_cell_topology_ref_map
        for definition of the mapping function. Presently supported cell topologies are

        \li     1D:   \c Line<2>
        \li     2D:   \c Triangle<3>, \c Triangle<6>, \c Quadrilateral<4>, \c Quadrilateral<9>
        \li     3D:   \c Tetrahedron<4>, \c Tetrahedron<10>, \c Hexahedron<8>, \c Hexahedron<27>

        \warning
        The array \c refPoints represents an arbitrary set of points in the reference
        frame that are not required to be in the reference cell corresponding to the
        specified cell topology. As a result, the images of these points under a given
        reference-to-physical map are not necessarily contained in the physical cell that
        is the image of the reference cell under that map. CellTools provides several
        inclusion tests methods to check whether or not the points are inside a reference cell.

        \param  physPoints        [out] - rank-3 array with dimensions (C,P,D) with the images of the ref. points
        \param  refPoints         [in]  - rank-3/2 array with dimensions (C,P,D)/(P,D) with points in reference frame
        \param  cellWorkset       [in]  - rank-3 array with dimensions (C,N,D) with the nodes of the cell workset
        \param  cellTopo          [in]  - cell topology of the cells stored in \c cellWorkset

    */
    template<typename PhysPointViewType,
             typename RefPointViewType,
             typename WorksetCellViewType>
    static void
    mapToPhysicalFrame(       PhysPointViewType    physPoints,
                        const RefPointViewType     refPoints,
                        const WorksetCellViewType  worksetCell,
                        const shards::CellTopology cellTopo ) {
      using nonConstRefPointValueType = typename RefPointViewType::non_const_value_type;
      auto basis = createHGradBasis<nonConstRefPointValueType,nonConstRefPointValueType>(cellTopo);
      mapToPhysicalFrame(physPoints, 
                         refPoints, 
                         worksetCell, 
                         basis);
    }


    /** \brief  Computes parameterization maps of 1- and 2-subcells of reference cells.

        Applies \f$\hat{\Phi}_i\f$, the parametrization map of a subcell \f$\hat{\mathcal{S}}_i\f$
        from a given reference cell, to a set of points in the parametrization domain
        \e R  of \f$\hat{\mathcal{S}}_i\f$. Returns a rank-2 array with dimensions
        (P,D) where for 1-subcells:
        \f[
        {subcellPoints}(p,*) = \hat{\Phi}_i(t_p) \quad\mbox{and}\quad
        \hat{\Phi}_i(t_p) = \left\{
        \begin{array}{ll}
        (\hat{x}(t_p),\hat{y}(t_p),\hat{z}(t_p)) & \mbox{for 3D parent cells}\\[1.5ex]
        (\hat{x}(t_p),\hat{y}(t_p))              & \mbox{for 2D parent cells}
        \end{array} \right.
        \quad t_p \in R = [-1,1] \,;
        \f]
        for 2-subcells:
        \f[
        {subcellPoints}(p,*) = \hat{\Phi}_i(u_p,v_p)\quad\mbox{and}\quad
        \hat{\Phi}_i(u_p,v_p) = (\hat{x}(u_p,v_p), \hat{y}(u_p,v_p), \hat{z}(u_p, v_p))
        \quad (u_p,v_p)\in R
        \f]
        and
        \f[
        R = \left\{\begin{array}{rl}
        \{(0,0),(1,0),(0,1)\} & \mbox{if face is Triangle} \\[1ex]
        [-1,1]\times [-1,1] & \mbox{if face is Quadrilateral}
        \end{array}\right.
        \f]

        \remarks
        \li     Parametrization of 1-subcells is defined for all 2D and 3D cell topologies with reference
        cells, including special 2D and 3D topologies such as shell and beams.
        \li     Parametrization of 2-subcells is defined for all 3D cell topologies with reference cells,
        including special 3D topologies such as shells.

        To map a set of points from a parametrization domain to a physical subcell workset, apply
        Intrepid2::CellTools::mapToPhysicalFrame to the output of this method.  This will effectively
        apply the parametrization map \f$ \Phi_{c,i} =  F_{c}\circ\hat{\Phi}_i \f$
        of each subcell in the workset to \c paramPoints. Here <var>c</var> is
        ordinal of a parent cell, relative to subcell workset, and <var>i</var> is subcell
        ordinal, relative to a reference cell, of the subcell workset. See Section
        \ref sec_cell_topology_subcell_wset for definition of subcell workset and Section
        \ref sec_cell_topology_subcell_map for definition of parametrization maps.

        \param  refSubcellPoints  [out] - rank-2 (P,D1) array with images of parameter space points
        \param  paramPoints       [in]  - rank-2 (P,D2) array with points in 1D or 2D parameter domain
        \param  subcellDim        [in]  - dimension of the subcell where points are mapped to
        \param  subcellOrd        [in]  - subcell ordinal
        \param  parentCell        [in]  - cell topology of the parent cell.
    */
    template<typename refSubcellViewType,
             typename paramPointViewType>
    static void
    mapToReferenceSubcell(       refSubcellViewType refSubcellPoints,
                           const paramPointViewType paramPoints,
                           const ordinal_type subcellDim,
                           const ordinal_type subcellOrd,
                           const shards::CellTopology parentCell );


    /** \brief  Computes parameterization maps of 1- and 2-subcells of reference cells.

        Overload of the previous function (see explanation above) where the subcell parametrization is used instead of 
        passing the parent cell topology.
    */

    template<typename refSubcellViewType, typename paramPointViewType>
    static void
    mapToReferenceSubcell(     refSubcellViewType                                             refSubcellPoints,
                         const paramPointViewType                                             paramPoints,
                         const typename RefSubcellParametrization<DeviceType>::ConstViewType  subcellParametrization,
                         const ordinal_type subcellOrd);

    /** \brief  Computes parameterization maps of 1- and 2-subcells of reference cells.

        Overload of the previous function (see explanation above) where the subcellOrd is a rank-1 array (P).
    */

    template<typename refSubcellViewType, typename paramPointViewType, typename ordViewType>
    static void
    mapToReferenceSubcell(     refSubcellViewType                                             refSubcellPoints,
                         const paramPointViewType                                             paramPoints,
                         const typename RefSubcellParametrization<DeviceType>::ConstViewType  subcellParametrization,
                         const ordViewType                                                    subcellOrd);


    //============================================================================================//
    //                                                                                            //
    //                      Physical-to-reference frame mapping and its inverse                   //
    //                                                                                            //
    //============================================================================================//


    /** \brief  Computes \f$ F^{-1}_{c} \f$, the inverse of the reference-to-physical frame map
        using a default initial guess.

        Applies \f$ F^{-1}_{c} \f$ for \b all cells in a cell workset to \b multiple point sets
        having the same number of points, indexed by cell ordinal, and stored in a rank-3
        (C,P,D) array. Returns a rank-3 (C,P,D) array such that 
        \f[
        \mbox{refPoints}(c,p,d) = \Big(F^{-1}_c(physPoint(c,p,*)) \Big)_d
        \f]

        Requires cell topology with a reference cell. See Section \ref sec_cell_topology_ref_map
        for definition of the mapping function. Presently supported cell topologies are

        \li     1D:   \c Line<2>
        \li     2D:   \c Triangle<3>, \c Triangle<6>, \c Quadrilateral<4>, \c Quadrilateral<9>
        \li     3D:   \c Tetrahedron<4>, \c Tetrahedron<10>, \c Hexahedron<8>, \c Hexahedron<27>

        \warning
        Computation of the inverse map in this method uses default selection of the initial guess
        based on cell topology:
        \li     \c Line topologies: line center (0)
        \li     \c Triangle topologies: the point (1/3, 1/3)
        \li     \c Quadrilateral topologies: the point (0, 0)
        \li     \c Tetrahedron topologies: the point (1/6, 1/6, 1/6)
        \li     \c Hexahedron topologies: the point (0, 0, 0)
        \li     \c Wedge topologies: the point (1/2, 1/2, 0)
        For some cells with extended topologies, these initial guesses may not be good enough
        for Newton's method to converge in the allotted number of iterations. A version of this
        method with user-supplied initial guesses is also available.

        \warning
        The array \c physPoints represents an arbitrary set (or sets) of points in the physical
        frame that are not required to belong in the physical cell (cells) that define(s) the reference
        to physical mapping. As a result, the images of these points in the reference frame
        are not necessarily contained in the reference cell corresponding to the specified
        cell topology.

        \param  refPoints         [out] - rank-3 array with dimensions (C,P,D) with the images of the physical points
        \param  physPoints        [in]  - rank-3 array with dimensions (C,P,D) with points in physical frame
        \param  worksetCell       [in]  - rank-3 array with dimensions (C,N,D) with the nodes of the cell workset
        \param  cellTopo          [in]  - cell topology of the cells stored in \c cellWorkset

    */
    template<typename refPointValueType,    class ...refPointProperties,
             typename physPointValueType,   class ...physPointProperties,
             typename worksetCellValueType, class ...worksetCellProperties>
    static void
    mapToReferenceFrame(       Kokkos::DynRankView<refPointValueType,refPointProperties...>    refPoints,
                         const Kokkos::DynRankView<physPointValueType,physPointProperties...>  physPoints,
                         const Kokkos::DynRankView<worksetCellValueType,worksetCellProperties...>   worksetCell,
                         const shards::CellTopology cellTopo );

    /** \brief  Computation of \f$ F^{-1}_{c} \f$, the inverse of the reference-to-physical frame map
        using user-supplied initial guess.

        Applies \f$ F^{-1}_{c} \f$ for \b all cells in a cell workset to \b multiple point sets
        having the same number of points, indexed by cell ordinal, and stored in a rank-3
        (C,P,D) array. Returns a rank-3 (C,P,D) array such that 
        \f[
        \mbox{refPoints}(c,p,d) = \Big(F^{-1}_c(physPoint(c,p,*)) \Big)_d
        \f]

        Requires pointer to HGrad basis that defines reference to physical cell mapping.  
        See Section \ref sec_cell_topology_ref_map for definition of the mapping function. 

        \warning
        The array \c physPoints represents an arbitrary set (or sets) of points in the physical
        frame that are not required to belong in the physical cell (cells) that define(s) the reference
        to physical mapping. As a result, the images of these points in the reference frame
        are not necessarily contained in the reference cell corresponding to the specified
        cell topology.

        \param  refPoints         [out] - rank-3/2 array with dimensions (C,P,D)/(P,D) with the images of the physical points
        \param  initGuess         [in]  - rank-3/2 array with dimensions (C,P,D)/(P,D) with the initial guesses for each point
        \param  physPoints        [in]  - rank-3/2 array with dimensions (C,P,D)/(P,D) with points in physical frame
        \param  worksetCell       [in]  - rank-3 array with dimensions (C,N,D) with the nodes of the cell workset
        \param  basis             [in]  - pointer to HGrad basis used for reference to physical cell mapping
    */
    template<typename refPointValueType,    class ...refPointProperties,
             typename initGuessValueType,   class ...initGuessProperties,
             typename physPointValueType,   class ...physPointProperties,
             typename worksetCellValueType, class ...worksetCellProperties,
             typename HGradBasisPtrType>
    static void
    mapToReferenceFrameInitGuess(       Kokkos::DynRankView<refPointValueType,refPointProperties...>       refPoints,
                                  const Kokkos::DynRankView<initGuessValueType,initGuessProperties...>     initGuess,
                                  const Kokkos::DynRankView<physPointValueType,physPointProperties...>     physPoints,
                                  const Kokkos::DynRankView<worksetCellValueType,worksetCellProperties...> worksetCell,
                                  const HGradBasisPtrType basis );


    /** \brief  Computation of \f$ F^{-1}_{c} \f$, the inverse of the reference-to-physical frame map
        using user-supplied initial guess.

        Applies \f$ F^{-1}_{c} \f$ for \b all cells in a cell workset to \b multiple point sets
        having the same number of points, indexed by cell ordinal, and stored in a rank-3
        (C,P,D) array. Returns a rank-3 (C,P,D) array such that 
        \f[
        \mbox{refPoints}(c,p,d) = \Big(F^{-1}_c(physPoint(c,p,*)) \Big)_d
        \f]

        Requires cell topology with a reference cell. See Section \ref sec_cell_topology_ref_map
        for definition of the mapping function. Presently supported cell topologies are

        \li     1D:   \c Line<2>
        \li     2D:   \c Triangle<3>, \c Triangle<6>, \c Quadrilateral<4>, \c Quadrilateral<9>
        \li     3D:   \c Tetrahedron<4>, \c Tetrahedron<10>, \c Hexahedron<8>, \c Hexahedron<27>

        \warning
        The array \c physPoints represents an arbitrary set (or sets) of points in the physical
        frame that are not required to belong in the physical cell (cells) that define(s) the reference
        to physical mapping. As a result, the images of these points in the reference frame
        are not necessarily contained in the reference cell corresponding to the specified
        cell topology.

        \param  refPoints         [out] - rank-3/2 array with dimensions (C,P,D)/(P,D) with the images of the physical points
        \param  initGuess         [in]  - rank-3/2 array with dimensions (C,P,D)/(P,D) with the initial guesses for each point
        \param  physPoints        [in]  - rank-3/2 array with dimensions (C,P,D)/(P,D) with points in physical frame
        \param  worksetCell       [in]  - rank-3 array with dimensions (C,N,D) with the nodes of the cell workset
        \param  cellTopo          [in]  - cell topology of the cells stored in \c cellWorkset
    */
    template<typename refPointValueType,    class ...refPointProperties,
             typename initGuessValueType,   class ...initGuessProperties,
             typename physPointValueType,   class ...physPointProperties,
             typename worksetCellValueType, class ...worksetCellProperties>
    static void
    mapToReferenceFrameInitGuess(       Kokkos::DynRankView<refPointValueType,refPointProperties...>       refPoints,
                                  const Kokkos::DynRankView<initGuessValueType,initGuessProperties...>     initGuess,
                                  const Kokkos::DynRankView<physPointValueType,physPointProperties...>     physPoints,
                                  const Kokkos::DynRankView<worksetCellValueType,worksetCellProperties...> worksetCell,
                                  const shards::CellTopology cellTopo ) {
      auto basis = createHGradBasis<refPointValueType,refPointValueType>(cellTopo);
      mapToReferenceFrameInitGuess(refPoints,
                                   initGuess,
                                   physPoints,
                                   worksetCell,
                                   basis);
    }


    //============================================================================================//
    //                                                                                            //
    //                             Control Volume Coordinates                                     //
    //                                                                                            //
    //============================================================================================//

    /** \brief Computes coordinates of sub-control volumes in each primary cell.

      To build the system of equations for the control volume finite element method we
      need to compute geometric data for integration over control volumes. A control
      volume is polygon or polyhedron that surrounds a primary cell node and has
      vertices that include the surrounding primary cells' barycenter, edge midpoints,
      and face midpoints if in 3-d.

      When using element-based assembly of the discrete equations over the primary mesh,
      a single element will contain a piece of each control volume surrounding each of
      the primary cell nodes. This piece of control volume (sub-control volume) is
      always a quadrilateral in 2-d and a hexahedron in 3-d.

      In 2-d the sub-control volumes are defined in the following way:

      \verbatim

       Quadrilateral primary element:

           O________M________O
           |        |        |
           |   3    |   2    |     B = cell barycenter
           |        |        |     O = primary cell nodes
           M________B________M     M = cell edge midpoints
           |        |        |
           |   0    |   1    |     sub-control volumes 0, 1, 2, 3
           |        |        |
           O________M________O


       Triangle primary element:

                    O
                   / \
                  /   \             B = cell barycenter
                 /     \            O = primary cell nodes
                M   2   M           M = cell edge midpoints
               / \     / \
              /   \ B /   \         sub-control volumes 0, 1, 2
             /      |      \
            /   0   |   1   \
           O________M________O

      \endverbatim

      In 3-d the sub-control volumes are defined by the primary cell face
      centers and edge midpoints. The eight sub-control volumes for a
      hexahedron are shown below:

      \verbatim
             O__________E__________O
            /|         /|         /|
           E_|________F_|________E |
          /| |       /| |       /| |
         O_|_|______E_|_|______O | |      O = primary cell nodes
         | | E------|-|-F------|-|-E      B = cell barycenter
         | |/|      | |/|      | |/|      F = cell face centers
         | F-|------|-B-|------|-F |      E = cell edge midpoints
         |/| |      |/| |      |/| |
         E_|_|______F_|_|______E | |
         | | O------|-|-E------|-|-O
         | |/       | |/       | |/
         | E--------|-F--------|-E
         |/         |/         |/
         O__________E__________O

      \endverbatim

    \param subCVCoords     [out] - array containing sub-control volume coordinates
    \param cellCoords       [in] - array containing coordinates of primary cells
    \param primaryCell      [in] - primary cell topology

    */
    template<typename subcvCoordValueType, class ...subcvCoordProperties,
             typename cellCoordValueType,  class ...cellCoordProperties>
    static void 
    getSubcvCoords(       Kokkos::DynRankView<subcvCoordValueType,subcvCoordProperties...> subcvCoords, 
                    const Kokkos::DynRankView<cellCoordValueType,cellCoordProperties...>   cellCoords,
                    const shards::CellTopology primaryCell );
    
    //============================================================================================//
    //                                                                                            //
    //                                        Inclusion tests                                     //
    //                                                                                            //
    //============================================================================================//

    /** \brief  Checks if a point belongs to a reference cell.
        
        Requires cell topology with a reference cell.
        
        \param  point             [in]  - rank-1 view (D) of the point tested for inclusion
        \param  cellTopo          [in]  - cell topology 
        \param  threshold         [in]  - "tightness" of the inclusion test
        \return true if the point is in the closure of the specified reference cell and false otherwise.
    */
    template<typename PointViewType>
    static bool 
    checkPointInclusion( const PointViewType        point,
                         const shards::CellTopology cellTopo,
                         const typename ScalarTraits<typename PointViewType::value_type>::scalar_type thres = 
                               threshold<typename ScalarTraits<typename PointViewType::value_type>::scalar_type>() );



    /** \brief  Checks every point for inclusion in the reference cell of a given topology.
        The points can belong to a global set and stored in a rank-2 view (P,D) ,
        or to multiple sets indexed by a cell ordinal and stored in a rank-3 view (C,P,D).
        The cell topology key is a template argument.
        Requires cell topology with a reference cell.
        \param  inCell            [out] - rank-1 view (P) or rank-2 view (C,P). On return, its entries will be set to 1 or 0 depending on whether points are included in cells 
        \param  point             [in]  - rank-2 view (P,D) or rank-3 view (C,P,D) with reference coordinates of the points tested for inclusion
        \param  threshold         [in]  - "tightness" of the inclusion test
    */
    template<unsigned cellTopologyKey,
             typename OutputViewType,
             typename InputViewType>
    static void checkPointwiseInclusion(       OutputViewType inCell, 
                                         const InputViewType points,
                                         const typename ScalarTraits<typename InputViewType::value_type>::scalar_type thresh =
                                               threshold<typename ScalarTraits<typename InputViewType::value_type>::scalar_type>()); 



    /** \brief  Checks every point in multiple sets indexed by a cell ordinal for inclusion in the reference cell of a given topology.
        Requires cell topology with a reference cell. 

        \param  inRefCell         [out] - rank-2 view (C,P) with results from the pointwise inclusion test
        \param  refPoints         [in]  - rank-3 view (C,P,D)
        \param  cellTopo          [in]  - cell topology
        \param  threshold         [in]  - "tightness" of the inclusion test
    */
    template<typename InCellViewType,                                                       
             typename PointViewType>
    static void checkPointwiseInclusion(       InCellViewType inCell,                     
                                         const PointViewType points,                       
                                         const shards::CellTopology cellTopo,                                                       
                                         const typename ScalarTraits<typename PointViewType::value_type>::scalar_type thres = 
                                               threshold<typename ScalarTraits<typename PointViewType::value_type>::scalar_type>() );

    /** \brief  Checks every points for inclusion in physical cells from a cell workset.
                The points can belong to a global set and stored in a rank-2 (P,D) view,
                or to multiple sets indexed by a cell ordinal and stored in a rank-3 (C,P,D) view
                cells in a cell workset.

        \param  inCell            [out] - rank-2 view (P,D)  with results from the pointwise inclusion test
        \param  points            [in]  - rank-2 view (P,D) or rank-3 view (C,P,D) with the physical points
        \param  cellWorkset       [in]  - rank-3 view with dimensions (C,N,D) with the nodes of the cell workset
        \param  cellTopo          [in]  - cell topology
        \param  threshold         [in]  - tolerance for inclusion tests on the input points
      */
    template<typename inCellValueType, class ...inCellProperties,                                                       
             typename pointValueType, class ...pointProperties,                                                        
             typename cellWorksetValueType, class ...cellWorksetProperties>  
    static void checkPointwiseInclusion(       Kokkos::DynRankView<inCellValueType,inCellProperties...> inCell,                     
                                         const Kokkos::DynRankView<pointValueType,pointProperties...> points,                       
                                         const Kokkos::DynRankView<cellWorksetValueType,cellWorksetProperties...> cellWorkset,      
                                         const shards::CellTopology cellTopo,                                                       
                                         const typename ScalarTraits<pointValueType>::scalar_type thres = 
                                               threshold<typename ScalarTraits<pointValueType>::scalar_type>() );


    // //============================================================================================//
    // //                                                                                            //
    // //                                           Debug                                            //
    // //                                                                                            //
    // //============================================================================================//


    // /** \brief  Prints the reference space coordinates of the vertices of the specified subcell
    //     \param  subcellDim        [in]  - dimension of the subcell where points are mapped to
    //     \param  subcellOrd        [in]  - subcell ordinal
    //     \param  parentCell        [in]  - cell topology of the parent cell.
    // */
    // static void printSubcellVertices(const int subcellDim,
    //                                  const int subcellOrd,
    //                                  const shards::CellTopology & parentCell);



    // /** \brief  Prints the nodes of a subcell from a cell workset

    //   */
    // template<class ArrayCell>
    // static void printWorksetSubcell(const ArrayCell &             cellWorkset,
    //                                 const shards::CellTopology &  parentCell,
    //                                 const int&                    pCellOrd,
    //                                 const int&                    subcellDim,
    //                                 const int&                    subcellOrd,
    //                                 const int&                    fieldWidth = 3);
  };

  //============================================================================================//
  //                                                                                            //
  //                  Validation of input/output arguments for CellTools methods                //
  //                                                                                            //
  //============================================================================================//

  /** \brief  Validates arguments to Intrepid2::CellTools::setJacobian
      \param  jacobian          [in]  - rank-4 (C,P,D,D) array array required
      \param  points            [in]  - rank-2 (P,D) or rank-3 (C,P,D) array required
      \param  cellWorkset       [in]  - rank-3 (C,N,D) array required
      \param  cellTopo          [in]  - cell topology with a reference cell required
  */
  template<typename jacobianViewType,
           typename PointViewType,
           typename worksetCellViewType>
  static void
  CellTools_setJacobianArgs( const jacobianViewType     jacobian,
                             const PointViewType        points,
                             const worksetCellViewType  worksetCell,
                             const shards::CellTopology cellTopo );

  /** \brief  Validates arguments to Intrepid2::CellTools::setJacobianInv
      \param  jacobianInv       [in]  - rank and dimensions must match jacobian array
      \param  jacobian          [in]  - rank-4 (C,P,D,D) array required
  */
  template<typename jacobianInvViewType,
           typename jacobianViewType>
  static void
  CellTools_setJacobianInvArgs( const jacobianInvViewType jacobianInv,
                                const jacobianViewType    jacobian );


  /** \brief  Validates arguments to Intrepid2::CellTools::setJacobianDet
      \param  jacobianDet       [in]  - rank = (jacobian rank - 2) required
      \param  jacobian          [in]  - rank-4 (C,P,D,D) array required
  */
  template<typename jacobianDetViewType,
           typename jacobianViewType>
  static void
  CellTools_setJacobianDetArgs( const jacobianDetViewType jacobianDet,
                                const jacobianViewType    jacobian );


  /** \brief  Validates arguments to Intrepid2::CellTools::mapToPhysicalFrame
      \param  physPoints        [in]  - rank-3 (C,P,D) array array required
      \param  refPoints         [in]  - rank-3 (C,P,D) array or rank-2 (P,D) array required
      \param  cellWorkset       [in]  - rank-3 (C,N,D) array required
      \param  cellTopo          [in]  - cell topology with a reference cell required
  */
  template<typename physPointViewType,
           typename refPointViewType,
           typename worksetCellViewType>
  static void
  CellTools_mapToPhysicalFrameArgs( const physPointViewType    physPoints,
                                    const refPointViewType     refPoints,
                                    const worksetCellViewType  worksetCell,
                                    const shards::CellTopology cellTopo );


  /** \brief  Validates arguments to Intrepid2::CellTools::mapToReferenceFrame with default initial guess.
      \param  physPoints        [in]  - rank-3 (C,P,D) array required
      \param  refPoints         [in]  - rank-3 (C,P,D) array or rank-2 (P,D) array required
      \param  cellWorkset       [in]  - rank-3 (C,N,D) array required
      \param  cellTopo          [in]  - cell topology with a reference cell required
  */
  template<typename refPointViewType, 
           typename physPointViewType, 
           typename worksetCellViewType>
  static void
  CellTools_mapToReferenceFrameArgs( const refPointViewType     refPoints,
                                     const physPointViewType    physPoints,
                                     const worksetCellViewType  worksetCell,
                                     const shards::CellTopology cellTopo );



  /** \brief  Validates arguments to Intrepid2::CellTools::mapToReferenceFrame with user-defined initial guess.
      \param  physPoints        [in]  - rank-3 (C,P,D) array required
      \param  initGuess         [in]  - rank and dimensions must match those of physPoints
      \param  refPoints         [in]  - rank-3 (C,P,D) array required
      \param  cellWorkset       [in]  - rank-3 (C,N,D) array required
      \param  cellTopo          [in]  - cell topology with a reference cell required
  */
  template<typename refPointViewType, 
           typename initGuessViewType, 
           typename physPointViewType, 
           typename worksetCellViewType>
  static void
  CellTools_mapToReferenceFrameInitGuess( const refPointViewType     refPoints,
                                          const initGuessViewType    initGuess,
                                          const physPointViewType    physPoints,
                                          const worksetCellViewType  worksetCell,
                                          const shards::CellTopology cellTopo );

}

#include "Intrepid2_CellToolsDocumentation.hpp"

#include "Intrepid2_CellToolsDefValidateArguments.hpp"
#include "Intrepid2_CellToolsDefNodeInfo.hpp"

#include "Intrepid2_CellToolsDefJacobian.hpp"
#include "Intrepid2_CellToolsDefRefToPhys.hpp"
#include "Intrepid2_CellToolsDefPhysToRef.hpp"

#include "Intrepid2_CellToolsDefControlVolume.hpp"

#include "Intrepid2_CellToolsDefInclusion.hpp"


#endif

