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

/** \file   Intrepid_CellTools.hpp
    \brief  Header file for the Intrepid2::CellTools class.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_CELLTOOLS_HPP__
#define __INTREPID2_CELLTOOLS_HPP__

#include "Intrepid2_ConfigDefs.hpp"

#include "Shards_CellTopology.hpp"
#include "Shards_BasicTopologies.hpp"

#include "Intrepid2_Types.hpp"
#include "Intrepid2_Utils.hpp"

#include "Intrepid2_RealSpaceTools.hpp"

#include "Intrepid2_Basis.hpp"

#include "Intrepid2_HGRAD_LINE_C1_FEM.hpp"
#include "Intrepid2_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid2_HGRAD_HEX_C1_FEM.hpp"

//#include "Intrepid2_HGRAD_TRI_C1_FEM.hpp"
//#include "Intrepid2_HGRAD_TET_C1_FEM.hpp"
//#include "Intrepid2_HGRAD_WEDGE_C1_FEM.hpp"
//#include "Intrepid2_HGRAD_PYR_C1_FEM.hpp"

//#include "Intrepid2_HGRAD_QUAD_C2_FEM.hpp"
//#include "Intrepid2_HGRAD_HEX_C2_FEM.hpp"

//#include "Intrepid2_HGRAD_TRI_C2_FEM.hpp"
//#include "Intrepid2_HGRAD_TET_C2_FEM.hpp"
//#include "Intrepid2_HGRAD_TET_COMP12_FEM.hpp"

//#include "Intrepid2_HGRAD_WEDGE_C2_FEM.hpp"
//#include "Intrepid2_HGRAD_WEDGE_I2_FEM.hpp"

#include "Kokkos_Core.hpp"

namespace Intrepid2 {
  
  //nn  
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
  

  template<typename ExecSpaceType>
  class CellTools {
  public:

    typedef double value_type;

    /** \brief  Checks if a cell topology has reference cell
        \param  cell              [in]  - cell topology
        \return true if the cell topology has reference cell, false oterwise
    */
    KOKKOS_FORCEINLINE_FUNCTION
    static bool 
    hasReferenceCell( const shards::CellTopology cellTopo ) {
      bool r_val = false;
      switch ( cellTopo.getKey() ) {
      case shards::Line<2>::key:
      case shards::Line<3>::key:
      case shards::ShellLine<2>::key:
      case shards::ShellLine<3>::key:
      case shards::Beam<2>::key:
      case shards::Beam<3>::key:

      // case shards::Triangle<3>::key:
      // case shards::Triangle<4>::key:
      // case shards::Triangle<6>::key:
      // case shards::ShellTriangle<3>::key:
      // case shards::ShellTriangle<6>::key:

      case shards::Quadrilateral<4>::key:
      case shards::Quadrilateral<8>::key:
      case shards::Quadrilateral<9>::key:
        //case shards::ShellQuadrilateral<4>::key:
        //case shards::ShellQuadrilateral<8>::key:
        //case shards::ShellQuadrilateral<9>::key:

      // case shards::Tetrahedron<4>::key:
      // case shards::Tetrahedron<8>::key:
      // case shards::Tetrahedron<10>::key:
      // case shards::Tetrahedron<11>::key:

      case shards::Hexahedron<8>::key:
      case shards::Hexahedron<20>::key:
      case shards::Hexahedron<27>::key:

      // case shards::Pyramid<5>::key:
      // case shards::Pyramid<13>::key:
      // case shards::Pyramid<14>::key:

      // case shards::Wedge<6>::key:
      // case shards::Wedge<15>::key:
      // case shards::Wedge<18>::key:
        r_val = true;
      }
      return r_val;
    }

  private:

    // reference nodes initialized 
    typedef Kokkos::DynRankView<const value_type,Kokkos::LayoutRight,Kokkos::HostSpace> referenceNodeDataViewHostType;
    struct ReferenceNodeDataStatic {
      value_type line[2][3], line_3[3][3];
      value_type triangle[3][3], triangle_4[4][3], triangle_6[6][3];
      value_type quadrilateral[4][3], quadrilateral_8[8][3], quadrilateral_9[9][3];
      value_type tetrahedron[4][3], tetrahedron_8[8][3], tetrahedron_10[10][3], tetrahedron_11[10][3];
      value_type hexahedron[8][3], hexahedron_20[20][3], hexahedron_27[27][3];
      value_type pyramid[5][3], pyramid_13[13][3], pyramid_14[14][3];
      value_type wedge[6][3], wedge_15[15][3], wedge_18[18][3];
    };

    typedef Kokkos::DynRankView<value_type,Kokkos::LayoutRight,ExecSpaceType> referenceNodeDataViewType;
    struct ReferenceNodeData {
      referenceNodeDataViewType line, line_3;
      referenceNodeDataViewType triangle, triangle_4, triangle_6;
      referenceNodeDataViewType quadrilateral, quadrilateral_8, quadrilateral_9;
      referenceNodeDataViewType tetrahedron, tetrahedron_8, tetrahedron_10, tetrahedron_11;
      referenceNodeDataViewType hexahedron, hexahedron_20, hexahedron_27;
      referenceNodeDataViewType pyramid, pyramid_13, pyramid_14;
      referenceNodeDataViewType wedge, wedge_15, wedge_18;
    };

    static  const ReferenceNodeDataStatic refNodeDataStatic_;
    static        ReferenceNodeData       refNodeData_;

    static bool isReferenceNodeDataSet_;
    static void setReferenceNodeData();

    //============================================================================================//
    //                                                                                            //
    //          Parametrization coefficients of edges and faces of reference cells                //
    //                                                                                            //
    //============================================================================================//

    // static variables
    typedef Kokkos::DynRankView<value_type,ExecSpaceType> subcellParamViewType;
    struct SubcellParamData {
      subcellParamViewType dummy; 
      subcellParamViewType lineEdges;  // edge maps for 2d non-standard cells; shell line and beam
      subcellParamViewType triEdges, quadEdges; // edge maps for 2d standard cells
      subcellParamViewType shellTriEdges, shellQuadEdges; // edge maps for 3d non-standard cells; shell tri and quad
      subcellParamViewType tetEdges, hexEdges, pyrEdges, wedgeEdges; // edge maps for 3d standard cells
      subcellParamViewType shellTriFaces, shellQuadFaces; // face maps for 3d non-standard cells
      subcellParamViewType tetFaces, hexFaces, pyrFaces, wedgeFaces; // face maps for 3d standard cells
    };

    static SubcellParamData subcellParamData_;

    static bool isSubcellParametrizationSet_;    
    static void setSubcellParametrization();

    // // hgrad functions
    // struct CachedHgradBasis {
    //   struct c1 {
    //     static const Basis_HGRAD_LINE_C1_FEM<ExecSpaceType>     line;
    //     static const Basis_HGRAD_TRI_C1_FEM<ExecSpaceType>      tri;
    //     static const Basis_HGRAD_QUAD_C1_FEM<ExecSpaceType>     quad;
    //     static const Basis_HGRAD_TET_C1_FEM<ExecSpaceType>      tet;
    //     static const Basis_HGRAD_HEX_C1_FEM<ExecSpaceType>      hex;
    //     static const Basis_HGRAD_WEDGE_C1_FEM<ExecSpaceType>    wedge;
    //     static const Basis_HGRAD_PYR_C1_FEM<ExecSpaceType>      pyr;
    //   };
    //   // struct c2 {
    //   //   static const Basis_HGRAD_TRI_C2_FEM<ExecSpaceType>      tri;
    //   //   static const Basis_HGRAD_QUAD_C2_FEM<ExecSpaceType>     quad;
    //   //   static const Basis_HGRAD_TET_C2_FEM<ExecSpaceType>      tet;
    //   //   static const Basis_HGRAD_HEX_C2_FEM<ExecSpaceType>      hex;
    //   //   static const Basis_HGRAD_WEDGE_C2_FEM<ExecSpaceType>    wedge;
    //   // };  
    // };
    // static Basis<ExecSpaceType>* 
    // getHgradBasis( const shards::CellTopology cellTopo );



    /** \brief  Returns array with the coefficients of the parametrization maps for the edges or faces
        of a reference cell topology. 
        
        See CellTools<Scalar>::setSubcellParametrization and Section \ref sec_cell_topology_subcell_map 
        more information about parametrization maps.
        
        \param  subcellDim        [in]  - dimension of subcells whose parametrization map is returned
        \param  parentCell        [in]  - topology of the reference cell owning the subcells
        
        \return FieldContainer<double> with the coefficients of the parametrization map for all subcells
        of the specified dimension. 
    */
    static void
    getSubcellParametrization( /**/  subcellParamViewType &subcellParam,
                               const ordinal_type          subcellDim, 
                               const shards::CellTopology  parentCell );
    
    /** \brief  Defines orientation-preserving parametrizations of reference edges and faces of cell 
        topologies with reference cells. 
        
        Given an edge {V0, V1} of some reference cell, its parametrization is a mapping from
        [-1,1] onto the edge. Parametrization of a triangular face {V0,V1,V2} is mapping from
        the standard 2-simplex {(0,0,0), (1,0,0), (0,1,0)}, embedded in 3D onto that face. 
        Parametrization of a quadrilateral face {V0,V1,V2,V3} is mapping from the standard 
        2-cube {(-1,-1,0),(1,-1,0),(1,1,0),(-1,1,0)}, embedded in 3D, onto that face.  
        
        This method computes the coefficients of edge and face parametrization maps and stores
        them in static arrays owned by CellTools<Scalar>::getSubcellParametrization method. 
        All mappings are affine and orientation-preserving, i.e., they preserve the tangent
        and normal directions implied by the vertex order of the edge or the face relative to
        the reference cell:
        
        \li     the tangent on [-1,1] from -1 in the direction of 1 is mapped to a tangent on edge {V0,V1}
        from V0 in the direction of V1  (the forward direction of the edge determined by its 
        start and end vertices)
        
        \li     the normal in the direction of (0,0,1) to the standard 2-simplex {(0,0,0),(1,0,0),(0,1,0)} 
        and the standard 2-cube {(-1,-1,0),(1,-1,0),(1,1,0),(-1,1,0)} is mapped to a normal
        on {V0,V1,V2} and {V0,V1,V2,V3}, determined according to the right-hand rule 
        (see http://mathworld.wolfram.com/Right-HandRule.html for definition of right-hand rule
        and Section \ref Section sec_cell_topology_subcell_map for further details).
        
        Because faces of all reference cells supported in Intrepid are affine images of either
        the standard 2-simplex or the standard 2-cube, the coordinate functions of the respective
        parmetrization maps are linear polynomials in the parameter variables (u,v), i.e., they
        are of the form \c F_i(u,v)=C_0(i)+C_1(i)u+C_2(i)v;  \c 0<=i<3 (face parametrizations
        are supported only for 3D cells, thus parametrization maps have 3 coordinate functions).   
        As a result, application of these maps is independent of the face type which is convenient
        for cells such as Wedge or Pyramid that have both types of faces. Also, coefficients of
        coordinate functions for all faces can be stored together in the same array.
        
        \param  subcellParametrization [out]  - array with the coefficients of the parametrization map
        \param  subcellDim             [in]   - dimension of the subcells being parametrized (1 or 2)
        \param  parentCell             [in]   - topology of the parent cell owning the subcells.
    */
    static void 
    setSubcellParametrization( /**/  subcellParamViewType &subcellParam,
                               const ordinal_type          subcellDim,
                               const shards::CellTopology  parentCell );
    
public:
  
    CellTools() = default;
    ~CellTools() = default;
    
    //============================================================================================//
    //                                                                                            //
    //                     Jacobian, inverse Jacobian and Jacobian determinant                    //
    //                                                                                            //
    //============================================================================================//
    
    // /** \brief  Computes the Jacobian matrix \e DF of the reference-to-physical frame map \e F.
      
    //             There are three use cases:
    //     \li     Computes Jacobians \f$DF_{c}\f$ of the reference-to-physical map for a \b specified physical
    //             cell \f${\mathcal C}\f$ from a cell workset on a \b single set of reference points stored
    //             in a rank-2 (P,D) array;
    //     \li     Computes Jacobians \f$DF_{c}\f$ of the reference-to-physical map for \b all physical cells
    //             in a cell workset on a \b single set of reference points stored in a rank-2 (P,D) array;
    //     \li     Computes Jacobians \f$DF_{c}\f$ of the reference-to-physical map for \b all physical cells
    //             in a cell workset on \b multiple reference point sets having the same number of points, 
    //             indexed by cell ordinal, and stored in a rank-3 (C,P,D) array;
      
    //             For a single point set in a rank-2 array (P,D) and \c whichCell set to a valid cell ordinal
    //             relative to \c cellWorkset returns a rank-3 (P,D,D) array such that
    //     \f[ 
    //             \mbox{jacobian}(p,i,j)   = [DF_{c}(\mbox{points}(p))]_{ij} \quad \mbox{for $0\le c < C$ - fixed} 
    //     \f]
    //             For a single point set in a rank-2 array (P,D) and \c whichCell=-1 (default value) returns
    //             a rank-4 (C,P,D,D) array such that
    //     \f[
    //             \mbox{jacobian}(c,p,i,j) = [DF_{c}(\mbox{points}(p))]_{ij} \quad c=0,\ldots, C
    //     \f]
    //             For multiple sets of reference points in a rank-3 (C,P,D) array returns 
    //             rank-4 (C,P,D,D) array such that
    //     \f[ 
    //             \mbox{jacobian}(c,p,i,j) = [DF_{c}(\mbox{points}(c,p))]_{ij} \quad c=0,\ldots, C
    //     \f]
    //             This setting requires the default value \c whichCell=-1. 
      
    //             Requires cell topology with a reference cell. See Section \ref sec_cell_topology_ref_map_DF 
    //             for definition of the Jacobian. 
      
    //             The default \c whichCell = -1 forces computation of all cell Jacobians and 
    //             requiers rank-4 output array. Computation of single cell Jacobians is forced by  
    //             selecting a valid cell ordinal \c whichCell and requires rank-3 output array.
      
    //             \warning
    //             The points are not required to be in the reference cell associated with the specified 
    //             cell topology. CellTools provides several inclusion tests methods to check whether 
    //             or not the points are inside a reference cell.
      
    //     \param  jacobian          [out] - rank-4/3 array with dimensions (C,P,D,D)/(P,D,D) with the Jacobians
    //     \param  points            [in]  - rank-2/3 array with dimensions (P,D)/(C,P,D) with the evaluation points 
    //     \param  cellWorkset       [in]  - rank-3 array with dimensions (C,N,D) with the nodes of the cell workset
    //     \param  cellTopo          [in]  - cell topology of the cells stored in \c cellWorkset
    //     \param  whichCell         [in]  - cell ordinal (for single cell Jacobian computation); default is -1      
    //  */
    // template<typename jacobianValueType, class ...jacobianProperties,
    //          typename pointValueType,    class ...pointProperties,
    //          typename worksetCellValueType,  class ...worksetCellProperties,
    //          typename HgradBasisType>
    // static void 
    // setJacobian( /**/  Kokkos::DynRankView<jacobianValueType,...jacobianProperties> jacobian,
    //              const Kokkos::DynRankView<pointValueType,...pointProperties>       points,
    //              const Kokkos::DynRankView<worksetCellValueType,...worksetCellProperties>   worksetCell,
    //              const shards::CellTopology cellTopo,
    //              const ordinal_type         whichCell = -1 );
    
    // /** \brief  Computes the inverse of the Jacobian matrix \e DF of the reference-to-physical frame map \e F.
        
    //     Returns rank-4 or rank-3 array with dimensions (C,P,D,D) or (P,D,D) such that 
    //     \f[ 
    //     \mbox{jacobianInv}(c,p,*,*) = \mbox{jacobian}^{-1}(c,p,*,*) \quad c = 0,\ldots, C
    //     \quad\mbox{or}\quad
    //     \mbox{jacobianInv}(p,*,*)   = \mbox{jacobian}^{-1}(p,*,*) 
    //     \f]
        
    //     \param  jacobianInv       [out] - rank-4/3 array with dimensions (C,P,D,D)/(P,D,D) with the inverse Jacobians
    //     \param  jacobian          [in]  - rank-4/3 array with dimensions (C,P,D,D)/(P,D,D) with the Jacobians
    // */
    // template<typename jacobianInvValueType, class ...jacobianInvProperties,
    //          typename jacobianValueType,    class ...jacobianProperties>
    // static void 
    // setJacobianInv( /**/  Kokkos::DynRankView<jacobianInvValueType,jacobianInvProperties...> jacobianInv,
    //                 const Kokkos::DynRankView<jacobianValueType,jacobianProperties...>       jacobian );
    
    // /** \brief  Computes the determinant of the Jacobian matrix \e DF of the reference-to-physical frame map \e F.
        
    //     Returns rank-2 or rank-1 array with dimensions (C,P)/(P) such that 
    //     \f[ 
    //     \mbox{jacobianDet}(c,p) = \mbox{det}(\mbox{jacobian}(c,p,*,*)) \quad c=0,\ldots, C 
    //     \quad\mbox{or}\quad
    //     \mbox{jacobianDet}(p)   = \mbox{det}(\mbox{jacobian}(p,*,*)) 
    //     \f]
        
    //     \param  jacobianDet       [out] - rank-2/1 array with dimensions (C,P)/(P) with Jacobian determinants
    //     \param  jacobian          [in]  - rank-4/3 array with dimensions (C,P,D,D)/(P,D,D) with the Jacobians
    // */
    // template<typename jacobianDetValueType, class ...jacobianDetProperties,
    //          typename jacobianValueType,    class ...jacobianProperties>
    // static void 
    // setJacobianDet( /**/  Kokkos::DynRankView<jacobianDetValueType,jacobianDetProperties...>  jacobianDet,
    //                 const Kokkos::DynRankView<jacobianValueType,jacobianProperties...>        jacobian );
    

    //============================================================================================//
    //                                                                                            //
    //                     Node information                                                       //
    //                                                                                            //
    //============================================================================================//

    // the node information can be used inside of kokkos functor and needs kokkos inline and 
    // exception should be an abort. for now, let's not decorate 

    /** \brief  Retrieves the Cartesian coordinates of a reference cell vertex. 
      
                Requires cell topology with a reference cell. Vertex coordinates are always returned 
                as an (x,y,z)-triple regardlesss of the actual topological cell dimension. The unused 
                coordinates are set to zero, e.g., vertex 0 of Line<2> is returned as {-1,0,0}.                
                
        \param  cell              [in]  - cell topology of the cell
        \param  vertexOrd         [in]  - vertex ordinal 
        \return pointer to array with the (x,y,z) coordinates of the specified reference vertex
    */
    template<typename cellVertexValueType, class ...cellVertexProperties>
    static void 
    getReferenceVertex( /**/  Kokkos::DynRankView<cellVertexValueType,cellVertexProperties...> cellVertices, 
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
    getReferenceSubcellVertices( /**/  Kokkos::DynRankView<subcellVertexValueType,subcellVertexProperties...> subcellVertices,
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
               with base topology this method is equivalent to CellTools<Scalar>::getReferenceVertex.
      
        \param  cell              [in]  - cell topology of the cell
        \param  vertexOrd         [in]  - node ordinal 
        \return pointer to array with the (x,y,z) coordinates of the specified reference vertex
    */
    template<typename cellNodeValueType, class ...cellNodeProperties>
    static void 
    getReferenceNode( /**/  Kokkos::DynRankView<cellNodeValueType,cellNodeProperties...> cellNodes,
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
    template<typename subcellNodeValueType, class ...subcellNodeProperties>
    static void
    getReferenceSubcellNodes( /**/  Kokkos::DynRankView<subcellNodeValueType,subcellNodeProperties...> subcellNodes,
                              const ordinal_type         subcellDim,
                              const ordinal_type         subcellOrd,
                              const shards::CellTopology parentCell );

    // /** \brief  Computes constant tangent vectors to edges of 2D or 3D reference cells. 
        
    //     Returns rank-1 array with dimension (D), D=2 or D=3; such that
    //     \f[
    //     {refEdgeTangent}(*) = \hat{\bf t}_i = {\partial\hat{\Phi}_i(t)\over\partial t}\,,
    //     \f]
    //     where \f$\hat{\Phi}_i : R =[-1,1]\mapsto \hat{\mathcal E}_i\f$ is the parametrization map
    //     of the specified reference edge \f$\hat{\mathcal E}_i\f$, given by
    //     \f[
    //     \hat{\Phi}_i(t) = \left\{\begin{array}{ll}
    //     (\hat{x}(t),\hat{y}(t),\hat{z}(t)) & \mbox{for 3D parent cells} \\[1ex]
    //     (\hat{x}(t),\hat{y}(t))            & \mbox{for 2D parent cells} \\[1ex]
    //     \end{array}\right.
    //     \f]
    //     The length of computed edge tangents is one-half the length of their associated edges:
    //     \f[
    //     |\hat{\bf t}_i | = {1\over 2} |\hat{\mathcal E}_i |\,.
    //     \f]
    //     Because the edges of all reference cells are always affine images of [-1,1],
    //     the edge tangent is constant vector field. 
        
    //     \param  refEdgeTangent    [out] - rank-1 array (D) with the edge tangent; D = cell dimension
    //     \param  edgeOrd           [in]  - ordinal of the edge whose tangent is computed
    //     \param  parentCell        [in]  - cell topology of the parent reference cell
    // */
    // template<typename refEdgeTangentValueType, class ...refEdgeTangentProperties>
    // static void 
    // getReferenceEdgeTangent( /**/ Kokkos::DynRankView<refEdgeTangentValueType,refEdgeTangentProperties...> refEdgeTangent,
    //                          const ordinal_type         edgeOrd,
    //                          const shards::CellTopology parentCell );

    // /** \brief  Computes pairs of constant tangent vectors to faces of a 3D reference cells. 
        
    //     Returns 2 rank-1 arrays with dimension (D), D=3, such that       
    //     \f[
    //     {refFaceTanU}(*) = \hat{\bf t}_{i,u} = {\partial\hat{\Phi}_i(u,v)\over\partial u} = 
    //     \left({\partial\hat{x}(u,v)\over\partial u}, 
    //     {\partial\hat{y}(u,v)\over\partial u},
    //     {\partial\hat{z}(u,v)\over\partial u} \right) ;
    //     \f]
    //     \f[
    //     {refFaceTanV}(*) = \hat{\bf t}_{i,v} = {\partial\hat{\Phi}_i(u,v)\over \partial v} = 
    //     \left({\partial\hat{x}(u,v)\over\partial v}, 
    //     {\partial\hat{y}(u,v)\over\partial v},
    //     {\partial\hat{z}(u,v)\over\partial v} \right)\,;
    //     \f]
    //     where \f$\hat{\Phi}_i: R \mapsto \hat{\mathcal F}_i\f$  
    //     is the parametrization map of the specified reference face \f$\hat{\mathcal F}_i\f$ given by
    //     \f[
    //     \hat{\Phi}_i(u,v) =(\hat{x}(u,v),\hat{y}(u,v),\hat{z}(u,v))
    //     \f]
    //     and 
    //     \f[
    //     R = \left\{\begin{array}{rl} 
    //     \{(0,0),(1,0),(0,1)\} & \mbox{if $\hat{\mathcal F}_i$  is Triangle} \\[1ex]
    //     [-1,1]\times [-1,1] & \mbox{if $\hat{\mathcal F}_i$ is Quadrilateral} \,.
    //     \end{array}\right.
    //     \f]
    //     Because the faces of all reference cells are always affine images of \e R , 
    //     the coordinate functions \f$\hat{x},\hat{y},\hat{z}\f$ of the parametrization map 
    //     are linear and the face tangents are constant vectors.  
        
    //     \param  refFaceTanU       [out] - rank-1 array (D) with (constant) tangent in u-direction
    //     \param  refFaceTanV       [out] - rank-1 array (D) with (constant) tangent in v-direction
    //     \param  faceOrd           [in]  - ordinal of the face whose tangents are computed
    //     \param  parentCell        [in]  - cell topology of the parent 3D reference cell
    // */
    // template<typename refFaceTanUValueType,...refFaceTanUProperties,
    //          typename refFaceTanVValueType,...refFaceTanVProperties>
    // static void 
    // getReferenceFaceTangents( /**/  Kokkos::DynRankView<refFaceTanUValueType,refFaceTanUProperties...> refFaceTanU,
    //                           /**/  Kokkos::DynRankView<refFaceTanVValueType,refFaceTanVProperties...> refFaceTanV,
    //                           const ordinal_type         faceOrd,
    //                           const shards::CellTopology parentCell );
    
    // /** \brief  Computes constant normal vectors to sides of 2D or 3D reference cells. 
        
    //     A side is defined as a subcell of dimension one less than that of its parent cell. 
    //     Therefore, sides of 2D cells are 1-subcells (edges) and sides of 3D cells
    //     are 2-subcells (faces).
        
    //     Returns rank-1 array with dimension (D), D = 2 or 3 such that
    //     \f[
    //     {refSideNormal}(*) = \hat{\bf n}_i =
    //     \left\{\begin{array}{rl} 
    //     \displaystyle
    //     \left({\partial\hat{\Phi}_i(t)\over\partial t}\right)^{\perp} 
    //     & \mbox{for 2D parent cells} \\[2ex]
    //     \displaystyle
    //     {\partial\hat{\Phi}_{i}\over\partial u} \times 
    //     {\partial\hat{\Phi}_{i}\over\partial v}   & \mbox{for 3D parent cells} 
    //     \end{array}\right.
    //     \f]
    //     where \f$ (u_1,u_2)^\perp = (u_2, -u_1)\f$, and \f$\hat{\Phi}_i: R \mapsto \hat{\mathcal S}_i\f$ 
    //     is the parametrization map of the specified reference side \f$\hat{\mathcal S}_i\f$ given by
    //     \f[
    //     \hat{\Phi}_i(u,v) = 
    //     \left\{\begin{array}{rl}
    //     (\hat{x}(t),\hat{y}(t))                   & \mbox{for 2D parent cells} \\[1ex]
    //     (\hat{x}(u,v),\hat{y}(u,v),\hat{z}(u,v))  & \mbox{for 3D parent cells}
    //     \end{array}\right.
        
    //     \f]
    //     For sides of 2D cells \e R=[-1,1] and for sides of 3D cells 
    //     \f[
    //     R = \left\{\begin{array}{rl} 
    //     \{(0,0),(1,0),(0,1)\}   & \mbox{if $\hat{\mathcal S}_i$ is Triangle} \\[1ex]
    //     [-1,1]\times [-1,1] & \mbox{if $\hat{\mathcal S}_i$ is Quadrilateral} \,.
    //     \end{array}\right.
    //     \f]
    //     For 3D cells the length of computed side normals is proportional to side area:
    //     \f[
    //     |\hat{\bf n}_i | = \left\{\begin{array}{rl} 
    //     2 \mbox{Area}(\hat{\mathcal F}_i) & \mbox{if $\hat{\mathcal F}_i$  is Triangle} \\[1ex]
    //     \mbox{Area}(\hat{\mathcal F}_i) & \mbox{if $\hat{\mathcal F}_i$ is Quadrilateral} \,.
    //     \end{array}\right.
    //     \f]
    //     For 2D cells the length of computed side normals is proportional to side length:
    //     \f[
    //     |\hat{\bf n}_i | = {1\over 2} |\hat{\mathcal F}_i |\,.
    //     \f]
    //     Because the sides of all reference cells are always affine images of \e R , 
    //     the coordinate functions \f$\hat{x},\hat{y},\hat{z}\f$ of the parametrization maps 
    //     are linear and the side normal is a constant vector.  
        
    //     \remark
    //     - For 3D cells the reference side normal coincides with the face normal computed by
    //     CellTools<Scalar>::getReferenceFaceNormal and these two methods are completely equivalent.
    //     - For 2D cells the reference side normal is defined by \f$\hat{{\bf n}}= \hat{\bf t}^\perp = (t_2,-t_1)\f$
    //     where \f$\hat{{\bf t}}=(t_1,t_2)\f$ is the tangent vector computed by 
    //     CellTools<Scalar>::getReferenceEdgeTangent. Therefore, the pair 
    //     \f$(\hat{{\bf n}},\hat{{\bf t}})\f$ is positively oriented.
        
    //     \param  refSideNormal     [out] - rank-1 array (D) with (constant) side normal
    //     \param  sideOrd           [in]  - ordinal of the side whose normal is computed
    //     \param  parentCell        [in]  - cell topology of the parent reference cell
    // */
    // template<typename refSideNormalValueType, class ...refSideNormalProperties>
    // static void
    // getReferenceSideNormal( /**/  Kokkos::DynRankView<refSideNormalValueType,refSideNormalProperties...> refSideNormal,
    //                         const ordinal_type         sideOrd,
    //                         const shards::CellTopology parentCell );

    // /** \brief  Computes constant normal vectors to faces of 3D reference cell. 
        
    //     Returns rank-1 array with dimension (D), D=3 such that
    //     \f[
    //     {refFaceNormal}(*) = \hat{\bf n}_i = {\partial\hat{\Phi}_{i}\over\partial u} \times 
    //     {\partial\hat{\Phi}_{i}\over\partial v}
    //     \f]
    //     where \f$\hat{\Phi}_i: R \mapsto \hat{\mathcal F}_i\f$  
    //     is the parametrization map of the specified reference face \f$\hat{\mathcal F}_i\f$ given by
    //     \f[
    //     \hat{\Phi}_i(u,v) =(\hat{x}(u,v),\hat{y}(u,v),\hat{z}(u,v))
    //     \f]
    //     and 
    //     \f[
    //     R = \left\{\begin{array}{rl} 
    //     \{(0,0),(1,0),(0,1)\} & \mbox{if ${\mathcal F}$  is Triangle} \\[1ex]
    //     [-1,1]\times [-1,1] & \mbox{if ${\mathcal F}$ is Quadrilateral} \,.
    //     \end{array}\right.
    //     \f]
    //     The length of computed face normals is proportional to face area:
    //     \f[
    //     |\hat{\bf n}_i | = \left\{\begin{array}{rl} 
    //     2 \mbox{Area}(\hat{\mathcal F}_i) & \mbox{if $\hat{\mathcal F}_i$  is Triangle} \\[1ex]
    //     \mbox{Area}(\hat{\mathcal F}_i) & \mbox{if $\hat{\mathcal F}_i$ is Quadrilateral} \,.
    //     \end{array}\right.
    //     \f]
    //     Because the faces of all reference cells are always affine images of \e R , 
    //     the coordinate functions \f$\hat{x},\hat{y},\hat{z}\f$ of the parametrization map 
    //     are linear and the face normal is a constant vector.  
        
    //     \remark
    //     The method CellTools::getReferenceFaceTangents computes the reference face tangents
    //     \f${\partial\hat{\Phi}_{i}/\partial u}\f$ and \f${\partial\hat{\Phi}_{i}/\partial v}\f$.
        
    //     \param  refFaceNormal     [out] - rank-1 array (D) with (constant) face normal
    //     \param  faceOrd           [in]  - ordinal of the face whose normal is computed
    //     \param  parentCell        [in]  - cell topology of the parent reference cell
    // */
    // template<typename refFaceNormalValueType, class ...refFaceNormalProperties>>
    // static void 
    // getReferenceFaceNormal( /**/  Kokkos::DynRankView<refFaceNormalValueType,refFaceNormalProperties...> refFaceNormal,
    //                         const ordinal_type         faceOrd,
    //                         const shards::CellTopology parentCell );
    
    // /** \brief  Computes non-normalized tangent vectors to physical edges in an edge workset 
    //     \f$\{\mathcal{E}_{c,i}\}_{c=0}^{N}\f$; (see \ref sec_cell_topology_subcell_wset for definition of edge worksets). 
        
    //     For every edge in the workset the tangents are computed at the points 
    //     \f${\bf x}_p = F_c(\hat{\Phi}_i(t_p))\in\mathcal{E}_{c,i}\f$ that are images of points
    //     from <var>R=[-1,1]</var> on edge \f$\mathcal{E}_{c,i}\f$. Returns rank-3 array with 
    //     dimensions (C,P,D1), D1=2 or D1=3 such that 
    //     \f[
    //     {edgeTangents}(c,p,d) = 
    //     DF_c(\hat{\Phi}_i(t_p))\cdot {\partial{\hat{\Phi}}_{i}(t_p)\over\partial t}\,; \qquad t_p \in R
    //     \f]
    //     In this formula: 
    //     \li     \f$ DF_c \f$ is Jacobian of parent cell \f${\mathcal C}\f$ that owns physical edge \f${\mathcal E}_{c,i}\f$;
    //     \li     \f$ {\partial{\hat{\Phi}}_{i}/\partial t}\f$ is the (constant) tangent to reference edge
    //     \f$\hat{\mathcal E}_i\f$; see CellTools<Scalar>::getReferenceEdgeTangent that has the 
    //     same local ordinal as the edges in the workset;
    //     \li     \f$ \hat{\Phi}_i R\mapsto\hat{\mathcal E}_i \f$ is parametrization of reference edge \f$\hat{\mathcal E}_i\f$;
        
    //     \warning
    //     \c worksetJacobians must contain the values of \f$DF_c(\hat{\Phi}_i(t_p))\f$, 
    //     where \f$ t_p \in R=[-1,1] \f$, i.e., Jacobians of the parent cells evaluated at points 
    //     that are located on reference edge \f$\hat{\mathcal E}_i\f$ having the same local ordinal as
    //     the edges in the workset.
        
    //     \param  edgeTangents      [out] - rank-3 array (C,P,D1) with tangents on workset edges
    //     \param  worksetJacobians  [in]  - rank-4 array (C,P,D1,D1) with Jacobians evaluated at ref. edge points
    //     \param  worksetEdgeOrd    [in]  - edge ordinal, relative to ref. cell, of the edge workset
    //     \param  parentCell        [in]  - cell topology of the parent reference cell
    // */
    // template<typename edgeTangentValueType,     class ...edgeTangentProperties,
    //          typename worksetJacobianValueType, class ...worksetJacobianProperties>
    // static void
    // getPhysicalEdgeTangents( /**/  Kokkos::DynRankView<edgeTangentValueType,edgeTangentProperties...>         edgeTangents,
    //                          const Kokkos::DynRankView<worksetJacobianValueType,worksetJacobianProperties...> worksetJacobians,
    //                          const ordinal_type         worksetEdgeOrd,
    //                          const shards::CellTopology parentCell );
    
    // /** \brief  Computes non-normalized tangent vector pairs to physical faces in a face workset 
    //     \f$\{\mathcal{F}_{c,i}\}_{c=0}^{N}\f$; (see \ref sec_cell_topology_subcell_wset for definition of face worksets). 
        
    //     For every face in the workset the tangents are computed at the points 
    //     \f${\bf x}_p = F_c(\hat{\Phi}_i(u_p,v_p))\in\mathcal{F}_{c,i}\f$ that are images of points
    //     from the parametrization domain \e R  on face \f$\mathcal{F}_{c,i}\f$. 
    //     Returns 2 rank-3 arrays with dimensions (C,P,D), D=3 such that
    //     \f[
    //     {faceTanU}(c,p,d) = DF_c(\hat{\Phi}_i(u_p, v_p))\cdot {\partial\hat{\Phi}_i\over\partial u};\qquad
    //     {faceTanV}(c,p,d) = DF_c(\hat{\Phi}_i(u_p, v_p))\cdot {\partial\hat{\Phi}_{i}\over\partial v}\,;
    //     \qquad (u_p, v_p) \in R \,.
    //     \f]
    //     In this formula:
    //     \li     \f$ DF_c \f$ is Jacobian of parent cell \f${\mathcal C}\f$ that owns physical face \f${\mathcal F}_{c,i}\f$;
    //     \li     \f$ {\partial\hat{\Phi}_i/\partial u}, {\partial\hat{\Phi}_i/\partial v}\f$ are the (constant) 
    //     tangents on reference face \f$\hat{\mathcal F}_i\f$; see CellTools<Scalar>::getReferenceFaceTangents; 
    //     that has the same local ordinal as the faces in the workset;
    //     \li     \f$ \hat{\Phi}_i : R\mapsto \hat{\mathcal F}_i\f$ is parametrization of reference face \f$\hat{\mathcal F}_i\f$;
    //     \li     \e R  is the parametrization domain for reference face \f$\hat{\mathcal F}_i\f$:
    //     \f[
    //     R = 
    //     \left\{\begin{array}{rl} 
    //     \{(0,0),(1,0),(0,1)\} & \mbox{if $\hat{\mathcal F}_i$ is Triangle} \\[1ex]
    //     [-1,1]\times [-1,1] & \mbox{if $\hat{\mathcal F}_i$ is Quadrilateral}
    //     \end{array}\right.
    //     \f]
        
    //     \warning
    //     \c worksetJacobians must contain the values of \f$DF_c(\hat{\Phi}_i(u_p,v_p))\f$, 
    //     where \f$(u_p,v_p)\in R\f$, i.e., Jacobians of the parent cells evaluated at points 
    //     that are located on reference face \f$\hat{\mathcal F}_i\f$ having the same local ordinal as
    //     the faces in the workset.
        
    //     \param  faceTanU          [out] - rank-3 array (C,P,D), image of ref. face u-tangent at workset faces
    //     \param  faceTanV          [out] - rank-3 array (C,P,D), image of ref. face u-tangent at workset faces
    //     \param  worksetJacobians  [in]  - rank-4 array (C,P,D,D) with Jacobians at ref. face points
    //     \param  worksetFaceOrd    [in]  - face ordinal, relative to ref. cell, of the face workset
    //     \param  parentCell        [in]  - cell topology of the parent reference cell
    // */
    // template<class ArrayFaceTangentU, class ArrayFaceTangentV, class ArrayJac>
    // static void 
    // getPhysicalFaceTangents(ArrayFaceTangentU &           faceTanU,
    //                         ArrayFaceTangentV &           faceTanV,
    //                         const ArrayJac &              worksetJacobians,
    //                         const int &                   worksetFaceOrd,
    //                         const shards::CellTopology &  parentCell);
    
    // /** \brief  Computes non-normalized normal vectors to physical sides in a side workset 
    //     \f$\{\mathcal{S}_{c,i}\}_{c=0}^{N}\f$. 
        
    //     For every side in the workset the normals are computed at the points 
    //     \f${\bf x}_p = F_c(\hat{\Phi}_i(P_p))\in\mathcal{S}_{c,i}\f$ that are images of points 
    //     from the parametrization domain \e R  on side \f$\mathcal{S}_{c,i}\f$.   
    //     A side is defined as a subcell of dimension one less than that of its  parent cell. 
    //     Therefore, sides of 2D cells are 1-subcells (edges) and sides of 3D cells are 2-subcells (faces).
        
    //     Returns rank-3 array with dimensions (C,P,D), D = 2 or 3, such that 
    //     \f[
    //     {sideNormals}(c,p,d) = 
    //     \left\{\begin{array}{crl} 
    //     \displaystyle
    //     \left(DF_c(\hat{\Phi}_i(t_p))\cdot 
    //     {\partial{\hat{\Phi}}_{i}(t_p)\over\partial t}\right)^{\perp} &  t_p\in R 
    //     & \mbox{for 2D parent cells} \\[2ex]
    //     \displaystyle
    //     \left( DF_c(\hat{\Phi}_i(u_p, v_p))\cdot {\partial\hat{\Phi}_i\over\partial u}\right) \times
    //     \left( DF_c(\hat{\Phi}_i(u_p, v_p))\cdot {\partial\hat{\Phi}_i\over\partial v}\right) \,;
    //     & (u_p, v_p) \in R                          & \mbox{for 3D parent cells}
    //     \end{array}\right.
    //     \f]
    //     In this formula:
    //     \li     \f$ DF_c \f$ is Jacobian of parent cell \f${\mathcal C}\f$ that owns physical side \f${\mathcal S}_{c,i}\f$;
    //     \li     For 2D cells: \f$ {\partial{\hat{\Phi}}_{i}/\partial t}\f$ is the (constant) tangent to reference side (edge)
    //     \f$\hat{\mathcal S}_i\f$; see CellTools<Scalar>::getReferenceEdgeTangent, that has the 
    //     same local ordinal as the sides in the workset;
    //     \li     For 3D cells: \f$ {\partial\hat{\Phi}_i/\partial u}, {\partial\hat{\Phi}_i/\partial v}\f$ are the (constant) 
    //     tangents on reference side (face) \f$\hat{\mathcal S}_i\f$; see CellTools<Scalar>::getReferenceFaceTangents,
    //     that has the same local ordinal as the sides in the workset;
    //     \li     \f$ \hat{\Phi}_i : R\mapsto \hat{\mathcal S}_i\f$ is parametrization of reference side \f$\hat{\mathcal S}_i\f$;
    //     \li     \e R  is the parametrization domain for reference side \f$\hat{\mathcal S}_i\f$. For
    //     2D parent cells \e R=[-1,1] and for 3D parent cells
    //     \f[
    //     R = \left\{\begin{array}{rl} 
    //     \{(0,0),(1,0),(0,1)\} & \mbox{if $\hat{\mathcal S}_i$ is Triangle} \\[1ex]
    //     [-1,1]\times [-1,1] & \mbox{if $\hat{\mathcal S}_i$ is Quadrilateral}
    //     \end{array}\right.
    //     \f]
        
    //     \remark
    //     - For 3D cells the physical side normals coincides with the face normals computed by
    //     CellTools<Scalar>::getPhysicalFaceNormals and these two methods are completely equivalent.
    //     - For 2D cells the physical side normals are defined by \f${\bf n}=(t_2,-t_1)\f$
    //     where \f${\bf t}=(t_1,t_2)\f$ are the physical edge tangents computed by 
    //     CellTools<Scalar>::getPhysicalEdgeTangents. Therefore, the pairs \f$({\bf n},{\bf t})\f$ are positively oriented.
        
        
    //     \warning
    //     \c worksetJacobians must contain the values of \f$DF_c(\hat{\Phi}_i(P_p))\f$, 
    //     where \f$P_p\in R\f$, i.e., Jacobians of the parent cells evaluated at points 
    //     that are located on reference side \f$\hat{\mathcal S}_i\f$ having the same local ordinal as
    //     the sides in the workset.
        
    //     \param  sideNormals       [out] - rank-3 array (C,P,D), normals at workset sides
    //     \param  worksetJacobians  [in]  - rank-4 array (C,P,D,D) with Jacobians at ref. side points
    //     \param  worksetSideOrd    [in]  - side ordinal, relative to ref. cell, of the side workset
    //     \param  parentCell        [in]  - cell topology of the parent reference cell
    // */
    // template<class ArraySideNormal, class ArrayJac>
    // static void 
    // getPhysicalSideNormals(ArraySideNormal &             sideNormals,
    //                        const ArrayJac &              worksetJacobians,
    //                        const int &                   worksetSideOrd,
    //                        const shards::CellTopology &  parentCell);
    
    // /** \brief  Computes non-normalized normal vectors to physical faces in a face workset 
    //     \f$\{\mathcal{F}_{c,i}\}_{c=0}^{N}\f$; (see \ref sec_cell_topology_subcell_wset for definition of face worksets). 
        
    //     For every face in the workset the normals are computed at the points 
    //     \f${\bf x}_p = F_c(\hat{\Phi}_i(u_p,v_p))\in\mathcal{F}_{c,i}\f$ that are images of points
    //     from the parametrization domain \e R  on face \f$\mathcal{F}_{c,i}\f$.  
    //     Returns rank-3 array with dimensions (C,P,D), D=3, such that 
    //     \f[
    //     {faceNormals}(c,p,d) = 
    //     \left( DF_c(\hat{\Phi}_i(u_p, v_p))\cdot {\partial\hat{\Phi}_i\over\partial u}\right) \times
    //     \left( DF_c(\hat{\Phi}_i(u_p, v_p))\cdot {\partial\hat{\Phi}_i\over\partial v}\right) \,;
    //     \qquad (u_p, v_p) \in R \,.
    //     \f]
    //     In this formula:
    //     \li     \f$ DF_c \f$ is Jacobian of parent cell \f${\mathcal C}\f$ that owns physical face \f${\mathcal F}_{c,i}\f$;
    //     \li     \f$ {\partial\hat{\Phi}_i/\partial u}, {\partial\hat{\Phi}_i/\partial v}\f$ are the (constant) 
    //     tangents on reference face \f$\hat{\mathcal F}_i\f$; see CellTools<Scalar>::getReferenceFaceTangents;
    //     that has the same local ordinal as the faces in the workset;
    //     \li     \f$ \hat{\Phi}_i : R\mapsto \hat{\mathcal F}_i\f$ is parametrization of reference face \f$\hat{\mathcal F}_i\f$;
    //     \li     \e R  is the parametrization domain for reference face \f$\hat{\mathcal F}_i\f$:
    //     \f[
    //     R = \left\{\begin{array}{rl} 
    //     \{(0,0),(1,0),(0,1)\} & \mbox{if $\hat{\mathcal F}_i$ is Triangle} \\[1ex]
    //     [-1,1]\times [-1,1] & \mbox{if $\hat{\mathcal F}_i$ is Quadrilateral}
    //     \end{array}\right.
    //     \f]
        
    //     \warning
    //     \c worksetJacobians must contain the values of \f$DF_c(\hat{\Phi}_i(u_p,v_p))\f$, 
    //     where \f$(u_p,v_p)\in R\f$, i.e., Jacobians of the parent cells evaluated at points 
    //     that are located on reference face \f$\hat{\mathcal F}_i\f$ having the same local ordinal as
    //     the faces in the workset.
        
    //     \param  faceNormals       [out] - rank-3 array (C,P,D), normals at workset faces
    //     \param  worksetJacobians  [in]  - rank-4 array (C,P,D,D) with Jacobians at ref. face points
    //     \param  worksetFaceOrd    [in]  - face ordinal, relative to ref. cell, of the face workset
    //     \param  parentCell        [in]  - cell topology of the parent reference cell
    // */
    // template<class ArrayFaceNormal, class ArrayJac>
    // static void 
    // getPhysicalFaceNormals(ArrayFaceNormal &             faceNormals,
    //                        const ArrayJac &              worksetJacobians,
    //                        const int &                   worksetFaceOrd,
    //                        const shards::CellTopology &  parentCell);
    
    
    
    



    //============================================================================================//
    //                                                                                            //
    //                      Reference-to-physical frame mapping and its inverse                   //
    //                                                                                            //
    //============================================================================================//
    
    /** \brief  Computes \e F, the reference-to-physical frame map.
        
        There are 3 use cases:
        \li     Applies \f$ F_{c} \f$ for a \b specified physical cell \f${\mathcal C}\f$ from a cell 
        workset to a \b single point set stored in a rank-2 (P,D) array;
        \li     Applies \f$ F_{c} \f$ for \b all cells in a cell workset to a \b single point set stored 
        in a rank-2 (P,D) array;
        \li     Applies \f$ F_{c} \f$ for \b all cells in a cell workset to \b multiple point sets having  
        the same number of points, indexed by cell ordinal, and stored in a rank-3 (C,P,D) array;
        
        For a single point set in a rank-2 array (P,D) and \c whichCell set to a valid cell ordinal
        relative to \c cellWorkset returns a rank-2 (P,D) array such that
        \f[  
        \mbox{physPoints}(p,d)   = \Big(F_c(\mbox{refPoint}(p,*)) \Big)_d \quad \mbox{for $0\le c < C$ - fixed}    
        \f]
        For a single point set in a rank-2 array (P,D) and \c whichCell=-1 (default value) returns
        a rank-3 (C,P,D) array such that
        \f[
        \mbox{physPoints}(c,p,d) = \Big(F_c(\mbox{refPoint}(p,*)) \Big)_d \quad c=0,\ldots, C 
        \f]
        For multiple point sets in a rank-3 (C,P,D) array returns a rank-3 (C,P,D) array such that 
        \f[  
        \mbox{physPoints}(c,p,d) = \Big(F_c(\mbox{refPoint}(c,p,*)) \Big)_d \quad c=0,\ldots, C 
        \f]
        This corresponds to mapping multiple sets of reference points to a matching number of 
        physical cells and requires the default value \c whichCell=-1.
        
        Requires cell topology with a reference cell. See Section \ref sec_cell_topology_ref_map
        for definition of the mapping function. Presently supported cell topologies are
        
        \li     1D:   \c Line<2>
        \li     2D:   \c Triangle<3>, \c Triangle<6>, \c Quadrilateral<4>, \c Quadrilateral<9>
        \li     3D:   \c Tetrahedron<4>, \c Tetrahedron<10>, \c Hexahedron<8>, \c Hexahedron<27>
        
        The default \c whichCell = -1 requires rank-3 output array and
        forces application of all reference-to-physical frame mappings corresponding to the
        cells stored in \c cellWorkset. Application of a single mapping is forced by selecting 
        a valid cell ordinal \c whichCell and requires rank-2 input/output point arrays.  
        
        \warning
        The array \c refPoints represents an arbitrary set of points in the reference
        frame that are not required to be in the reference cell corresponding to the
        specified cell topology. As a result, the images of these points under a given 
        reference-to-physical map are not necessarily contained in the physical cell that
        is the image of the reference cell under that map. CellTools provides several 
        inclusion tests methods to check whether or not the points are inside a reference cell.
        
        \param  physPoints        [out] - rank-3/2 array with dimensions (C,P,D)/(P,D) with the images of the ref. points
        \param  refPoints         [in]  - rank-3/2 array with dimensions (C,P,D)/(P,D) with points in reference frame
        \param  cellWorkset       [in]  - rank-3 array with dimensions (C,N,D) with the nodes of the cell workset
        \param  cellTopo          [in]  - cell topology of the cells stored in \c cellWorkset
        \param  whichCell         [in]  - ordinal of the cell that defines the reference-to-physical map; default is -1 
        
        \todo   Implement method for non-standard (shell, beam, etc) topologies.
    */
    template<typename physPointValueType, class ...physPointProperties,
             typename refPointValueType,  class ...refPointProperties,
             typename worksetCellValueType,   class ...worksetCellProperties>
    static void 
    mapToPhysicalFrame( /**/  Kokkos::DynRankView<physPointValueType,physPointProperties...>  physPoints,
                        const Kokkos::DynRankView<refPointValueType,refPointProperties...>    refPoints,
                        const Kokkos::DynRankView<worksetCellValueType,worksetCellProperties...>      worksetCell,
                        const shards::CellTopology cellTopo );


    //============================================================================================//
    //                                                                                            //
    //                      Reference-to-physical frame mapping and its inverse                   //
    //                                                                                            //
    //============================================================================================//


    /** \brief  Computes \f$ F^{-1}_{c} \f$, the inverse of the reference-to-physical frame map
        using a default initial guess. 
        
        There are two use cases: 
        \li     Applies \f$ F^{-1}_{c} \f$ for a \b specified physical cell \f${\mathcal C}\f$ from a 
        cell workset to a \b single set of points stored in a rank-2 (P,D) array;
        \li     Applies \f$ F^{-1}_{c} \f$ for \b all cells in a cell workset to \b multiple point sets
        having the same number of points, indexed by cell ordinal, and stored in a rank-3 
        (C,P,D) array (default mode).
        
        For a single point set in a rank-2 array (P,D) returns a rank-2 (P,D) array such that      
        \f[            
        \mbox{refPoints}(p,d) = \Big(F^{-1}_c(physPoint(p,*)) \Big)_d         
        \f]
        The \c whichCell argument selects the physical cell and is required to be a valid cell
        ordinal for \c cellWorkset array.
        
        For multiple point sets in a rank-3 (C,P,D) array returns a rank-3 (C,P,D) array such that 
        \f[            
        \mbox{refPoints}(c,p,d) = \Big(F^{-1}_c(physPoint(c,p,*)) \Big)_d         
        \f]
        The default value \e whichCell=-1 selects this mode.  
        
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
        
        \param  refPoints         [out] - rank-3/2 array with dimensions (C,P,D)/(P,D) with the images of the physical points
        \param  physPoints        [in]  - rank-3/2 array with dimensions (C,P,D)/(P,D) with points in physical frame
        \param  cellWorkset       [in]  - rank-3 array with dimensions (C,N,D) with the nodes of the cell workset
        \param  whichCell         [in]  - ordinal of the cell that defines the reference-to-physical map; default is -1
        \param  cellTopo          [in]  - cell topology of the cells stored in \c cellWorkset   
        
        \todo   Implement method for non-standard (shell, beam, etc) topologies.
    */
    template<typename refPointValueType,  class ...refPointProperties,
             typename physPointValueType, class ...physPointProperties,
             typename worksetCellValueType,   class ...worksetCellProperties>
    static void 
    mapToReferenceFrame( /**/  Kokkos::DynRankView<refPointValueType,refPointProperties...>    refPoints,
                         const Kokkos::DynRankView<physPointValueType,physPointProperties...>  physPoints,
                         const Kokkos::DynRankView<worksetCellValueType,worksetCellProperties...>      worksetCell,
                         const shards::CellTopology cellTopo );
    
    /** \brief  Computation of \f$ F^{-1}_{c} \f$, the inverse of the reference-to-physical frame map
        using user-supplied initial guess. 
        
        There are two use cases: 
        \li       Applies \f$ F^{-1}_{c} \f$ for a \b specified physical cell \f${\mathcal C}\f$ from a 
        cell workset to a \b single set of points stored in a rank-2 (P,D) array;
        \li       Applies \f$ F^{-1}_{c} \f$ for \b all cells in a cell workset to \b multiple point sets
        having the same number of points, indexed by cell ordinal, and stored in a rank-3 
        (C,P,D) array (default mode).
        
        For a single point set in a rank-2 array (P,D) returns a rank-2 (P,D) array such that      
        \f[            
        \mbox{refPoints}(p,d) = \Big(F^{-1}_c(physPoint(p,*)) \Big)_d         
        \f]
        The \c whichCell argument selects the physical cell and is required to be a valid cell
        ordinal for \c cellWorkset array.
        
        For multiple point sets in a rank-3 (C,P,D) array returns a rank-3 (C,P,D) array such that 
        \f[            
        \mbox{refPoints}(c,p,d) = \Big(F^{-1}_c(physPoint(c,p,*)) \Big)_d         
        \f]
        The default value \c whichCell=-1 selects this mode.  
        
        Requires cell topology with a reference cell. See Section \ref sec_cell_topology_ref_map
        for definition of the mapping function. Presently supported cell topologies are
        
        \li       1D:   \c Line<2>
        \li       2D:   \c Triangle<3>, \c Triangle<6>, \c Quadrilateral<4>, \c Quadrilateral<9>
        \li       3D:   \c Tetrahedron<4>, \c Tetrahedron<10>, \c Hexahedron<8>, \c Hexahedron<27>
        
        \warning 
        The array \c physPoints represents an arbitrary set (or sets) of points in the physical
        frame that are not required to belong in the physical cell (cells) that define(s) the reference
        to physical mapping. As a result, the images of these points in the reference frame
        are not necessarily contained in the reference cell corresponding to the specified
        cell topology. 
        
        \param  refPoints         [out] - rank-3/2 array with dimensions (C,P,D)/(P,D) with the images of the physical points
        \param  initGuess         [in]  - rank-3/2 array with dimensions (C,P,D)/(P,D) with the initial guesses for each point
        \param  physPoints        [in]  - rank-3/2 array with dimensions (C,P,D)/(P,D) with points in physical frame
        \param  cellWorkset       [in]  - rank-3 array with dimensions (C,N,D) with the nodes of the cell workset
        \param  whichCell         [in]  - ordinal of the cell that defines the reference-to-physical map; default is -1
        \param  cellTopo          [in]  - cell topology of the cells stored in \c cellWorkset      
    */
    template<typename refPointValueType,  class ...refPointProperties,
             typename initGuessValueType, class ...initGuessProperties,
             typename physPointValueType, class ...physPointProperties,
             typename worksetCellValueType,   class ...worksetCellProperties>
    static void
    mapToReferenceFrameInitGuess( /**/  Kokkos::DynRankView<refPointValueType,refPointProperties...>    refPoints,
                                  const Kokkos::DynRankView<initGuessValueType,initGuessProperties...>  initGuess,
                                  const Kokkos::DynRankView<physPointValueType,physPointProperties...>  physPoints,
                                  const Kokkos::DynRankView<worksetCellValueType,worksetCellProperties...>      worksetCell,
                                  const shards::CellTopology cellTopo );
    
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
        CellTools<Scalar>::mapToPhysicalFrame to the output of this method.  This will effectively 
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
    template<typename refSubcellPointValueType, class ...refSubcellPointProperties,
             typename paramPointValueType, class ...paramPointProperties>
    static void
    mapToReferenceSubcell( /**/  Kokkos::DynRankView<refSubcellPointValueType,refSubcellPointProperties...> refSubcellPoints,
                           const Kokkos::DynRankView<paramPointValueType,paramPointProperties...>           paramPoints,
                           const ordinal_type subcellDim,
                           const ordinal_type subcellOrd,
                           const shards::CellTopology parentCell );
           
    // //============================================================================================//
    // //                                                                                            //
    // //                                        Inclusion tests                                     //
    // //                                                                                            //
    // //============================================================================================//
    
    // /** \brief  Checks if a point belongs to a reference cell. 
      
    //             Requires cell topology with a reference cell.
      
    //     \param  point             [in]  - spatial coordinates of the point tested for inclusion
    //     \param  pointDim          [in]  - spatial dimension of that point
    //     \param  cellTopo          [in]  - cell topology of the cells stored in \c cellWorkset
    //     \param  threshold         [in]  - "tightness" of the inclusion test
    //     \return 1 if the point is in the closure of the specified reference cell and 0 otherwise.
    // */ 
    // static int checkPointInclusion(const Scalar*                 point,
    //                                const int                     pointDim,
    //                                const shards::CellTopology &  cellTopo,
    //                                const double &                threshold = INTREPID2_THRESHOLD);
    
    
    
    // /** \brief  Checks if a set of points belongs to a reference cell. 
      
    //             Requires cell topology with a reference cell. See Intrepid2::CellTools::checkPointwiseInclusion 
    //             for admissible ranks and dimensions of the input point array.
      
    //     \param  points            [in]  - rank-1, 2 or 3 array (point, vector of points, matrix of points)
    //     \param  cellTopo          [in]  - cell topology of the cells stored in \c cellWorkset
    //     \param  threshold         [in]  - "tightness" of the inclusion test
    //     \return 1 if all points are in the closure of the specified reference cell
    //             0 if at least one point is outside the closure of the reference cell 
                
    // */
    // template<class ArrayPoint>
    // static int checkPointsetInclusion(const ArrayPoint &            points,
    //                                   const shards::CellTopology &  cellTopo, 
    //                                   const double &                threshold = INTREPID2_THRESHOLD);
    
    
    
    // /** \brief  Checks every point in a set for inclusion in a reference cell.  
      
    //             Requires cell topology with a reference cell. Admissible ranks and dimensions of the 
    //             input point array and the corresponding rank and dimension of the output array are as follows:
    //             \verbatim
    //             |-------------------|-------------|-------------|-------------|
    //             |  rank: (in)/(out) |    1/1      |     2/1     |    3/2      |
    //             |-------------------|-------------|-------------|-------------|
    //             |  points    (in)   |     (D)     |    (I, D)   |  (I, J, D)  |
    //             |-------------------|-------------|-------------|-------------|
    //             |  inRefCell (out)  |     (1)     |    (I)      |  (I, J)     |
    //             |------------------ |-------------|-------------|-------------|
    //             \endverbatim
    //             Example: if \c points is rank-3 array with dimensions (I, J, D), then
    //     \f[
    //            \mbox{inRefCell}(i,j) = 
    //              \left\{\begin{array}{rl} 
    //                 1 & \mbox{if $points(i,j,*)\in\hat{\mathcal{C}}$} \\[2ex]
    //                 0 & \mbox{if $points(i,j,*)\notin\hat{\mathcal{C}}$} 
    //              \end{array}\right.
    //       \f]
    //     \param  inRefCell         [out] - rank-1 or 2 array with results from the pointwise inclusion test
    //     \param  refPoints         [in]  - rank-1,2 or 3 array (point, vector of points, matrix of points)
    //     \param  cellTopo          [in]  - cell topology of the cells stored in \c cellWorkset
    //     \param  threshold         [in]  - "tightness" of the inclusion test
    // */
    // template<class ArrayIncl, class ArrayPoint>
    // static void checkPointwiseInclusion(ArrayIncl &                   inRefCell,
    //                                     const ArrayPoint &            points,
    //                                     const shards::CellTopology &  cellTopo, 
    //                                     const double &                threshold = INTREPID2_THRESHOLD);
    
    
    
    // /** \brief  Checks every point in a set or multiple sets for inclusion in physical cells from a cell workset.
      
    //             There are two use cases:
    //     \li     Checks every point from a \b single point set, stored in a rank-2 array (P,D) for inclusion
    //             in a \b specified physical cell \f${\mathcal C}\f$ from a cell workset.
    //     \li     Checks every point from \b multiple point sets indexed by a cell ordinal, and stored in a rank-3
    //             (C,P,D) array, for inclusion in the physical cell having the same cell ordinal, for \b all 
    //             cells in a cell workset. 
      
    //             For a single point set in a rank-2 array (P,D) and \c whichCell set to a valid cell ordinal
    //             relative to \c cellWorkset returns a rank-1 (P) array such that
    //     \f[
    //             \mbox{inCell}(p) = 
    //               \left\{\begin{array}{rl} 
    //                  1 & \mbox{if $points(p,*)\in {\mathcal{C}}$} \\ [2ex]
    //                  0 & \mbox{if $points(p,*)\notin {\mathcal{C}}$} 
    //               \end{array}\right.
    //     \f]
    //             For multiple point sets in a rank-3 array (C,P,D) and \c whichCell=-1 (default value)
    //             returns a rank-2 (C,P) array such that
    //     \f[
    //             \mbox{inCell}(c,p) = 
    //               \left\{\begin{array}{rl} 
    //                   1 & \mbox{if $points(c,p,*)\in {\mathcal{C}}$} \\ [2ex]
    //                   0 & \mbox{if $points(c,p,*)\notin {\mathcal{C}}$} 
    //             \end{array}\right.
    //     \f]

    //     \param  inCell            [out] - rank-1  array with results from the pointwise inclusion test
    //     \param  points            [in]  - rank-2 array with dimensions (P,D) with the physical points
    //     \param  cellWorkset       [in]  - rank-3 array with dimensions (C,N,D) with the nodes of the cell workset
    //     \param  cellTopo          [in]  - cell topology of the cells stored in \c cellWorkset
    //     \param  whichCell         [in]  - ordinal of the cell used in the inclusion test
    //     \param  threshold         [in]  - tolerance for inclusion tests on the input points
    //   */
    // template<class ArrayIncl, class ArrayPoint, class ArrayCell>
    // static void checkPointwiseInclusion(ArrayIncl &                   inCell,
    //                                     const ArrayPoint &            points,
    //                                     const ArrayCell &             cellWorkset,
    //                                     const shards::CellTopology &  cell,
    //                                     const int &                   whichCell = -1, 
    //                                     const double &                threshold = INTREPID2_THRESHOLD);
    
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

    // //============================================================================================//
    // //                                                                                            //
    // //                             Control Volume Coordinates                                     //
    // //                                                                                            //
    // //============================================================================================//

    // /** \brief Computes coordinates of sub-control volumes in each primary cell.

    //   To build the system of equations for the control volume finite element method we
    //   need to compute geometric data for integration over control volumes. A control
    //   volume is polygon or polyhedron that surrounds a primary cell node and has
    //   vertices that include the surrounding primary cells' barycenter, edge midpoints,
    //   and face midpoints if in 3-d.

    //   When using element-based assembly of the discrete equations over the primary mesh,
    //   a single element will contain a piece of each control volume surrounding each of
    //   the primary cell nodes. This piece of control volume (sub-control volume) is
    //   always a quadrilateral in 2-d and a hexahedron in 3-d.

    //   In 2-d the sub-control volumes are defined in the following way:

    //   \verbatim

    //    Quadrilateral primary element:

    //        O________M________O
    //        |        |        |
    //        |   3    |   2    |     B = cell barycenter
    //        |        |        |     O = primary cell nodes
    //        M________B________M     M = cell edge midpoints
    //        |        |        |
    //        |   0    |   1    |     sub-control volumes 0, 1, 2, 3
    //        |        |        |
    //        O________M________O


    //    Triangle primary element:

    //                 O
    //                / \
    //               /   \             B = cell barycenter
    //              /     \            O = primary cell nodes
    //             M   2   M           M = cell edge midpoints
    //            / \     / \
    //           /   \ B /   \         sub-control volumes 0, 1, 2
    //          /      |      \
    //         /   0   |   1   \
    //        O________M________O

    //   \endverbatim

    //   In 3-d the sub-control volumes are defined by the primary cell face
    //   centers and edge midpoints. The eight sub-control volumes for a
    //   hexahedron are shown below:

    //   \verbatim
    //          O__________E__________O
    //         /|         /|         /|
    //        E_|________F_|________E |
    //       /| |       /| |       /| |
    //      O_|_|______E_|_|______O | |      O = primary cell nodes
    //      | | E------|-|-F------|-|-E      B = cell barycenter
    //      | |/|      | |/|      | |/|      F = cell face centers
    //      | F-|------|-B-|------|-F |      E = cell edge midpoints
    //      |/| |      |/| |      |/| |
    //      E_|_|______F_|_|______E | |
    //      | | O------|-|-E------|-|-O
    //      | |/       | |/       | |/
    //      | E--------|-F--------|-E
    //      |/         |/         |/
    //      O__________E__________O

    //   \endverbatim

    // \param subCVCoords     [out] - array containing sub-control volume coordinates
    // \param cellCoords       [in] - array containing coordinates of primary cells
    // \param primaryCell      [in] - primary cell topology

    // */
    
    // template<class ArrayCVCoord, class ArrayCellCoord>
    // static void getSubCVCoords(ArrayCVCoord & subCVcoords, const ArrayCellCoord & cellCoords,
    //                            const shards::CellTopology& primaryCell);

    // /** \brief Compute cell barycenters.
        
    //     \param barycenter      [out] - array containing cell baycenters
    //     \param cellCoords       [in] - array containing cell coordinates
        
    // */
    // template<class ArrayCent, class ArrayCellCoord>
    // static void getBarycenter(ArrayCent & barycenter, const ArrayCellCoord & cellCoords);

  };

  //============================================================================================//
  //                                                                                            //
  //                  Validation of input/output arguments for CellTools methods                //
  //                                                                                            //
  //============================================================================================//
  
  /** \brief  Validates arguments to Intrepid2::CellTools::setJacobian
      \param  jacobian          [in]  - rank-4 (C,P,D,D) array or rank-3 (P,D,D) array required
      \param  points            [in]  - rank-2 (P,D) or rank-3 (C,P,D) array required
      \param  cellWorkset       [in]  - rank-3 (C,N,D) array required
      \param  whichCell         [in]  - default = -1 or 0 <= whichCell < C required
      \param  cellTopo          [in]  - cell topology with a reference cell required
  */
  template<typename jacobianViewType,
           typename pointViewType,
           typename worksetCellViewType>
  static void
  CellTools_setJacobianArgs( const jacobianViewType     jacobian,
                             const pointViewType        points,
                             const worksetCellViewType  worksetCell,
                             const shards::CellTopology cellTopo );

  /** \brief  Validates arguments to Intrepid2::CellTools::setJacobianInv
      \param  jacobianInv       [in]  - rank and dimensions must match jacobian array
      \param  jacobian          [in]  - rank-4 (C,P,D,D) array or rank-3 (P,D,D) array required
  */
  template<typename jacobianInvViewType,
           typename jacobianViewType>
  static void
  CellTools_setJacobianInvArgs( const jacobianInvViewType jacobianInv,
                                const jacobianViewType    jacobian );
  
  
  /** \brief  Validates arguments to Intrepid2::CellTools::setJacobianDet
      \param  jacobianDet       [in]  - rank = (jacobian rank - 2) required
      \param  jacobian          [in]  - rank-4 (C,P,D,D) array or rank-3 (P,D,D) array required
  */
  template<typename jacobianDetViewType,
           typename jacobianViewType>
  static void
  CellTools_setJacobianDetArgs( const jacobianDetViewType jacobianDet,
                                const jacobianViewType    jacobian );

  
  /** \brief  Validates arguments to Intrepid2::CellTools::mapToPhysicalFrame
      \param  physPoints        [in]  - rank-3 (C,P,D) array or rank-2 (P,D) array required
      \param  refPoints         [in]  - rank-3 (C,P,D) array or rank-2 (P,D) array required
      \param  cellWorkset       [in]  - rank-3 (C,N,D) array required
      \param  whichCell         [in]  - default = -1 or 0 <= whichCell < C required
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
  
  
  // /** \brief  Validates arguments to Intrepid2::CellTools::mapToReferenceFrame with default initial guess.
  //     \param  physPoints        [in]  - rank-3 (C,P,D) array or rank-2 (P,D) array required
  //     \param  refPoints         [in]  - rank-3 (C,P,D) array or rank-2 (P,D) array required
  //     \param  cellWorkset       [in]  - rank-3 (C,N,D) array required
  //     \param  whichCell         [in]  - default = -1 or 0 <= whichCell < C required
  //     \param  cellTopo          [in]  - cell topology with a reference cell required
  // */
  // template<class ArrayRefPoint, class ArrayPhysPoint, class ArrayCell>
  // static void 
  // validateArguments_mapToReferenceFrame(const ArrayRefPoint  &        refPoints,
  //                                       const ArrayPhysPoint &        physPoints,
  //                                       const ArrayCell      &        cellWorkset,
  //                                       const shards::CellTopology &  cellTopo,
  //                                       const int&                    whichCell);
  
  
  
  // /** \brief  Validates arguments to Intrepid2::CellTools::mapToReferenceFrame with user-defined initial guess.
  //     \param  physPoints        [in]  - rank-3 (C,P,D) array or rank-2 (P,D) array required
  //     \param  initGuess         [in]  - rank and dimensions must match those of physPoints
  //     \param  refPoints         [in]  - rank-3 (C,P,D) array or rank-2 (P,D) array required
  //     \param  cellWorkset       [in]  - rank-3 (C,N,D) array required
  //     \param  whichCell         [in]  - default = -1 or 0 <= whichCell < C required
  //     \param  cellTopo          [in]  - cell topology with a reference cell required
  // */
  // template<class ArrayRefPoint, class ArrayInitGuess, class ArrayPhysPoint, class ArrayCell>
  // static void 
  // validateArguments_mapToReferenceFrame(const ArrayRefPoint  &        refPoints,
  //                                       const ArrayInitGuess &        initGuess,
  //                                       const ArrayPhysPoint &        physPoints,
  //                                       const ArrayCell      &        cellWorkset,
  //                                       const shards::CellTopology &  cellTopo,
  //                                       const int&                    whichCell);
  
  
  
  // /** \brief  Validates arguments to Intrepid2::CellTools::checkPointwiseInclusion
  //     \param  inCell            [out] - rank-1  (P) array required
  //     \param  physPoints        [in]  - rank-2  (P,D) array required
  //     \param  cellWorkset       [in]  - rank-3  (C,N,D) array required
  //     \param  whichCell         [in]  - 0 <= whichCell < C required
  //     \param  cellTopo          [in]  - cell topology with a reference cell required
  // */
  // template<class ArrayIncl, class ArrayPoint, class ArrayCell>
  // static void
  // validateArguments_checkPointwiseInclusion(ArrayIncl        &            inCell,
  //                                           const ArrayPoint &            physPoints,
  //                                           const ArrayCell  &            cellWorkset,
  //                                           const int &                   whichCell,
  //                                           const shards::CellTopology &  cell);
  
  

} 

#include "Intrepid2_CellToolsDocumentation.hpp"

#include "Intrepid2_CellToolsDefValidateArguments.hpp"
#include "Intrepid2_CellToolsDefNodeInfo.hpp"

#include "Intrepid2_CellToolsDefParametrization.hpp"
//#include "Intrepid2_CellToolsDefJacobian.hpp"

#include "Intrepid2_CellToolsDefRefToPhys.hpp"

// not yet converted ...

// #include "Intrepid2_CellToolsDefPhysToRef.hpp"
// #include "Intrepid2_CellToolsDefInclusion.hpp"

// #include "Intrepid2_CellToolsDefDebug.hpp"
// #include "Intrepid2_CellToolsDefControlVolume.hpp"

#endif

