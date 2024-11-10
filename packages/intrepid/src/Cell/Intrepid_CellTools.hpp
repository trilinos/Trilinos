#ifndef INTREPID_CELTOOLS_HPP
#define INTREPID_CELTOOLS_HPP

// @HEADER
// ************************************************************************
//
//                           Intrepid Package
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
// Questions? Contact Pavel Bochev  (pbboche@sandia.gov)
//                    Denis Ridzal  (dridzal@sandia.gov), or
//                    Kara Peterson (kjpeter@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_CellTools.hpp
    \brief  Header file for the Intrepid::CellTools class.
    \author Created by P. Bochev and D. Ridzal.
*/

#ifndef INTREPID_CELLTOOLS_HPP
#define INTREPID_CELLTOOLS_HPP


#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_RealSpaceTools.hpp"
#include "Intrepid_ConfigDefs.hpp"
#include "Intrepid_Types.hpp"
#include "Intrepid_Utils.hpp"
#include "Intrepid_Basis.hpp"
#include "Intrepid_HGRAD_TRI_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_TET_C1_FEM.hpp"
#include "Intrepid_HGRAD_WEDGE_C1_FEM.hpp"
#include "Intrepid_HGRAD_PYR_C1_FEM.hpp"
#include "Intrepid_HGRAD_HEX_C1_FEM.hpp"

#include "Intrepid_HGRAD_LINE_C1_FEM.hpp"

#include "Intrepid_HGRAD_TRI_C2_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_HGRAD_TET_C2_FEM.hpp"
#include "Intrepid_HGRAD_TET_COMP12_FEM.hpp"
#include "Intrepid_HGRAD_WEDGE_C2_FEM.hpp"
#include "Intrepid_HGRAD_WEDGE_I2_FEM.hpp"
#include "Intrepid_HGRAD_HEX_C2_FEM.hpp"
#include "Intrepid_HGRAD_HEX_I2_FEM.hpp"
#include "Intrepid_HGRAD_PYR_I2_FEM.hpp"

#include "Shards_CellTopology.hpp"
#include "Shards_BasicTopologies.hpp"

#include "Teuchos_Assert.hpp"
#include "Teuchos_RCP.hpp"

#include <Intrepid_Rank.hpp>

namespace Intrepid {
  
//nn  
  //============================================================================================//
  //                                                                                            //
  //                                          CellTools                                         //
  //                                                                                            //
  //============================================================================================//
  
/** \class  Intrepid::CellTools
    \brief  A stateless class for operations on cell data. Provides methods for:
    \li     computing Jacobians of reference-to-physical frame mappings, their inverses and determinants    
    \li     application of the reference-to-physical frame mapping and its inverse
    \li     parametrizations of edges and faces of reference cells needed for edge and face integrals,
    \li     computation of edge and face tangents and face normals on both reference and physical frames
    \li     inclusion tests for point sets in reference and physical cells.  
*/
 

template<class Scalar>
class CellTools {
private:
  
  //============================================================================================//
  //                                                                                            //
  //          Parametrization coefficients of edges and faces of reference cells                //
  //                                                                                            //
  //============================================================================================//
  
  
  /** \brief  Returns array with the coefficients of the parametrization maps for the edges or faces
              of a reference cell topology. 
  
              See CellTools<Scalar>::setSubcellParametrization and Section \ref sec_cell_topology_subcell_map 
              more information about parametrization maps.
   
      \param  subcellDim        [in]  - dimension of subcells whose parametrization map is returned
      \param  parentCell        [in]  - topology of the reference cell owning the subcells
  
      \return FieldContainer<double> with the coefficients of the parametrization map for all subcells
              of the specified dimension. 
    */
  static const FieldContainer<double>& getSubcellParametrization(const int                   subcellDim, 
                                                                 const shards::CellTopology& parentCell);
  
  
  
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
  
  static void setSubcellParametrization(FieldContainer<double>&     subcellParametrization,
                                        const int                   subcellDim,
                                        const shards::CellTopology& parentCell);
  
  //============================================================================================//
  //                                                                                            //
  //                  Validation of input/output arguments for CellTools methods                //
  //                                                                                            //
  //============================================================================================//
  
  /** \brief  Validates arguments to Intrepid::CellTools::setJacobian
      \param  jacobian          [in]  - rank-4 (C,P,D,D) array or rank-3 (P,D,D) array required
      \param  points            [in]  - rank-2 (P,D) or rank-3 (C,P,D) array required
      \param  cellWorkset       [in]  - rank-3 (C,N,D) array required
      \param  whichCell         [in]  - default = -1 or 0 <= whichCell < C required
      \param  cellTopo          [in]  - cell topology with a reference cell required
  */
  template<class ArrayJac, class ArrayPoint, class ArrayCell>
  static void validateArguments_setJacobian(const ArrayJac    &          jacobian,
                                            const ArrayPoint  &          points,
                                            const ArrayCell   &          cellWorkset,
                                            const int &                  whichCell,
                                            const shards::CellTopology & cellTopo);
  
  
  
  /** \brief  Validates arguments to Intrepid::CellTools::setJacobianInv
      \param  jacobianInv       [in]  - rank and dimensions must match jacobian array
      \param  jacobian          [in]  - rank-4 (C,P,D,D) array or rank-3 (P,D,D) array required
  */
  template<class ArrayJacInv, class ArrayJac>
  static void validateArguments_setJacobianInv(const ArrayJacInv &  jacobianInv,
                                               const ArrayJac    &  jacobian);
  
  
  
  /** \brief  Validates arguments to Intrepid::CellTools::setJacobianDet
      \param  jacobianDet       [in]  - rank = (jacobian rank - 2) required
      \param  jacobian          [in]  - rank-4 (C,P,D,D) array or rank-3 (P,D,D) array required
    */
  template<class ArrayJacDet, class ArrayJac>
  static void validateArguments_setJacobianDetArgs(const ArrayJacDet &  jacobianDet,
                                                   const ArrayJac    &  jacobian);
  
  
  
  /** \brief  Validates arguments to Intrepid::CellTools::mapToPhysicalFrame
      \param  physPoints        [in]  - rank-3 (C,P,D) array or rank-2 (P,D) array required
      \param  refPoints         [in]  - rank-3 (C,P,D) array or rank-2 (P,D) array required
      \param  cellWorkset       [in]  - rank-3 (C,N,D) array required
      \param  whichCell         [in]  - default = -1 or 0 <= whichCell < C required
      \param  cellTopo          [in]  - cell topology with a reference cell required
    */
  template<class ArrayPhysPoint, class ArrayRefPoint, class ArrayCell>
  static void validateArguments_mapToPhysicalFrame(const ArrayPhysPoint &        physPoints,
                                                   const ArrayRefPoint  &        refPoints,
                                                   const ArrayCell      &        cellWorkset,
                                                   const shards::CellTopology &  cellTopo,
                                                   const int&                    whichCell);
  
  
  
  /** \brief  Validates arguments to Intrepid::CellTools::mapToReferenceFrame with default initial guess.
      \param  physPoints        [in]  - rank-3 (C,P,D) array or rank-2 (P,D) array required
      \param  refPoints         [in]  - rank-3 (C,P,D) array or rank-2 (P,D) array required
      \param  cellWorkset       [in]  - rank-3 (C,N,D) array required
      \param  whichCell         [in]  - default = -1 or 0 <= whichCell < C required
      \param  cellTopo          [in]  - cell topology with a reference cell required
    */
  template<class ArrayRefPoint, class ArrayPhysPoint, class ArrayCell>
  static void validateArguments_mapToReferenceFrame(const ArrayRefPoint  &        refPoints,
                                                    const ArrayPhysPoint &        physPoints,
                                                    const ArrayCell      &        cellWorkset,
                                                    const shards::CellTopology &  cellTopo,
                                                    const int&                    whichCell);

  
  
  /** \brief  Validates arguments to Intrepid::CellTools::mapToReferenceFrame with user-defined initial guess.
      \param  physPoints        [in]  - rank-3 (C,P,D) array or rank-2 (P,D) array required
      \param  initGuess         [in]  - rank and dimensions must match those of physPoints
      \param  refPoints         [in]  - rank-3 (C,P,D) array or rank-2 (P,D) array required
      \param  cellWorkset       [in]  - rank-3 (C,N,D) array required
      \param  whichCell         [in]  - default = -1 or 0 <= whichCell < C required
      \param  cellTopo          [in]  - cell topology with a reference cell required
    */
  template<class ArrayRefPoint, class ArrayInitGuess, class ArrayPhysPoint, class ArrayCell>
  static void validateArguments_mapToReferenceFrame(const ArrayRefPoint  &        refPoints,
                                                    const ArrayInitGuess &        initGuess,
                                                    const ArrayPhysPoint &        physPoints,
                                                    const ArrayCell      &        cellWorkset,
                                                    const shards::CellTopology &  cellTopo,
                                                    const int&                    whichCell);
  
  
  
  /** \brief  Validates arguments to Intrepid::CellTools::checkPointwiseInclusion
      \param  inCell            [out] - rank-1  (P) array required
      \param  physPoints        [in]  - rank-2  (P,D) array required
      \param  cellWorkset       [in]  - rank-3  (C,N,D) array required
      \param  whichCell         [in]  - 0 <= whichCell < C required
      \param  cellTopo          [in]  - cell topology with a reference cell required
    */
  template<class ArrayIncl, class ArrayPoint, class ArrayCell>
  static void validateArguments_checkPointwiseInclusion(ArrayIncl        &            inCell,
                                                        const ArrayPoint &            physPoints,
                                                        const ArrayCell  &            cellWorkset,
                                                        const int &                   whichCell,
                                                        const shards::CellTopology &  cell);
public:
  
    /** \brief  Default constructor.
      */
    CellTools(){ };
      
      
    /** \brief  Destructor
      */
    ~CellTools(){ };
    
    //============================================================================================//
    //                                                                                            //
    //                     Jacobian, inverse Jacobian and Jacobian determinant                    //
    //                                                                                            //
    //============================================================================================//
    
    /** \brief  Computes the Jacobian matrix \e DF of the reference-to-physical frame map \e F.
      
                There are three use cases:
        \li     Computes Jacobians \f$DF_{c}\f$ of the reference-to-physical map for a \b specified physical
                cell \f${\mathcal C}\f$ from a cell workset on a \b single set of reference points stored
                in a rank-2 (P,D) array;
        \li     Computes Jacobians \f$DF_{c}\f$ of the reference-to-physical map for \b all physical cells
                in a cell workset on a \b single set of reference points stored in a rank-2 (P,D) array;
        \li     Computes Jacobians \f$DF_{c}\f$ of the reference-to-physical map for \b all physical cells
                in a cell workset on \b multiple reference point sets having the same number of points, 
                indexed by cell ordinal, and stored in a rank-3 (C,P,D) array;
      
                For a single point set in a rank-2 array (P,D) and \c whichCell set to a valid cell ordinal
                relative to \c cellWorkset returns a rank-3 (P,D,D) array such that
        \f[ 
                \mbox{jacobian}(p,i,j)   = [DF_{c}(\mbox{points}(p))]_{ij} \quad \mbox{for $0\le c < C$ - fixed} 
        \f]
                For a single point set in a rank-2 array (P,D) and \c whichCell=-1 (default value) returns
                a rank-4 (C,P,D,D) array such that
        \f[
                \mbox{jacobian}(c,p,i,j) = [DF_{c}(\mbox{points}(p))]_{ij} \quad c=0,\ldots, C
        \f]
                For multiple sets of reference points in a rank-3 (C,P,D) array returns 
                rank-4 (C,P,D,D) array such that
        \f[ 
                \mbox{jacobian}(c,p,i,j) = [DF_{c}(\mbox{points}(c,p))]_{ij} \quad c=0,\ldots, C
        \f]
                This setting requires the default value \c whichCell=-1. 
      
                Requires cell topology with a reference cell. See Section \ref sec_cell_topology_ref_map_DF 
                for definition of the Jacobian. 
      
                The default \c whichCell = -1 forces computation of all cell Jacobians and 
                requiers rank-4 output array. Computation of single cell Jacobians is forced by  
                selecting a valid cell ordinal \c whichCell and requires rank-3 output array.
      
                \warning
                The points are not required to be in the reference cell associated with the specified 
                cell topology. CellTools provides several inclusion tests methods to check whether 
                or not the points are inside a reference cell.
      
        \param  jacobian          [out] - rank-4/3 array with dimensions (C,P,D,D)/(P,D,D) with the Jacobians
        \param  points            [in]  - rank-2/3 array with dimensions (P,D)/(C,P,D) with the evaluation points 
        \param  cellWorkset       [in]  - rank-3 array with dimensions (C,N,D) with the nodes of the cell workset
        \param  cellTopo          [in]  - cell topology of the cells stored in \c cellWorkset
        \param  whichCell         [in]  - cell ordinal (for single cell Jacobian computation); default is -1      
     */

    template<class ArrayJac, class ArrayPoint, class ArrayCell, bool typecheck>
	struct setJacobianTempSpec;
	
	
	/*template<class ArrayJac, class ArrayPoint, class ArrayCell>
	static void setJacobianTemp(ArrayJac &               jacobian,
                            const ArrayPoint &           points,
                            const ArrayCell  &           cellWorkset,
                            const shards::CellTopology & cellTopo,
                            const int &                  whichCell = -1);*/
    template<class ArrayJac, class ArrayPoint, class ArrayCell>
    static void setJacobian(ArrayJac &                   jacobian,
                            const ArrayPoint &           points,
                            const ArrayCell  &           cellWorkset,
                            const shards::CellTopology & cellTopo,
                            const int &                  whichCell = -1);
    
    template<class ArrayJac, class ArrayPoint, class ArrayCell>
    static void setJacobian(ArrayJac &                   jacobian,
                            const ArrayPoint &           points,
                            const ArrayCell  &           cellWorkset,
                            const Teuchos::RCP< Basis< Scalar, FieldContainer<Scalar> > > HGRAD_Basis,
                            const int &                  whichCell = -1);
  
    
    /** \brief  Computes the inverse of the Jacobian matrix \e DF of the reference-to-physical frame map \e F.
      
                Returns rank-4 or rank-3 array with dimensions (C,P,D,D) or (P,D,D) such that 
        \f[ 
                \mbox{jacobianInv}(c,p,*,*) = \mbox{jacobian}^{-1}(c,p,*,*) \quad c = 0,\ldots, C
                \quad\mbox{or}\quad
                \mbox{jacobianInv}(p,*,*)   = \mbox{jacobian}^{-1}(p,*,*) 
        \f]
        
        \param  jacobianInv       [out] - rank-4/3 array with dimensions (C,P,D,D)/(P,D,D) with the inverse Jacobians
        \param  jacobian          [in]  - rank-4/3 array with dimensions (C,P,D,D)/(P,D,D) with the Jacobians
    */
    template<class ArrayJacInv, class ArrayJac>
    static void setJacobianInv(ArrayJacInv &     jacobianInv,
                               const ArrayJac &  jacobian);
                               
 /*   template<class ArrayJacInv, class ArrayJac>
    static void setJacobianInvTemp(ArrayJacInv &     jacobianInv,
                               const ArrayJac &  jacobian);   
    */
    
    /** \brief  Computes the determinant of the Jacobian matrix \e DF of the reference-to-physical frame map \e F.
      
                Returns rank-2 or rank-1 array with dimensions (C,P)/(P) such that 
        \f[ 
                \mbox{jacobianDet}(c,p) = \mbox{det}(\mbox{jacobian}(c,p,*,*)) \quad c=0,\ldots, C 
                \quad\mbox{or}\quad
                \mbox{jacobianDet}(p)   = \mbox{det}(\mbox{jacobian}(p,*,*)) 
          \f]
      
        \param  jacobianDet       [out] - rank-2/1 array with dimensions (C,P)/(P) with Jacobian determinants
        \param  jacobian          [in]  - rank-4/3 array with dimensions (C,P,D,D)/(P,D,D) with the Jacobians
      */
    template<class ArrayJacDet, class ArrayJac>
    static void setJacobianDet(ArrayJacDet &     jacobianDet,
                               const ArrayJac &  jacobian);
                               
  /*  template<class ArrayJacDet, class ArrayJac>
    static void setJacobianDetTemp(ArrayJacDet &     jacobianDet,
                               const ArrayJac &  jacobian);*/
    
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
    template<class ArrayPhysPoint, class ArrayRefPoint, class ArrayCell>
    static void mapToPhysicalFrame(ArrayPhysPoint      &         physPoints,
                                   const ArrayRefPoint &         refPoints,
                                   const ArrayCell     &         cellWorkset,
                                   const shards::CellTopology &  cellTopo,
                                   const int &                   whichCell = -1);
    
    template<class ArrayPhysPoint, class ArrayRefPoint, class ArrayCell>
    static void mapToPhysicalFrame(ArrayPhysPoint      &         physPoints,
                                   const ArrayRefPoint &         refPoints,
                                   const ArrayCell     &         cellWorkset,
                                   const Teuchos::RCP< Basis< Scalar, FieldContainer<Scalar> > > HGRAD_Basis,
                                   const int &                   whichCell = -1);
                                   
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
 /*       template<class ArrayPhysPoint, class ArrayRefPoint, class ArrayCell>
    static void mapToPhysicalFrameTemp(ArrayPhysPoint      &         physPoints,
                                   const ArrayRefPoint &         refPoints,
                                   const ArrayCell     &         cellWorkset,
                                   const shards::CellTopology &  cellTopo,
                                   const int &                   whichCell = -1);
                */                   
	    template<class ArrayPhysPoint, class ArrayRefPoint, class ArrayCell, int refRank,int phyptsrank>
	struct mapToPhysicalFrameTempSpec;
	
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
    template<class ArrayRefPoint, class ArrayPhysPoint, class ArrayCell>
    static void mapToReferenceFrame(ArrayRefPoint        &        refPoints,
                                    const ArrayPhysPoint &        physPoints,
                                    const ArrayCell      &        cellWorkset,
                                    const shards::CellTopology &  cellTopo,
                                    const int &                   whichCell = -1);
    
    
    template<class ArrayRefPoint, class ArrayPhysPoint, class ArrayCell>
    static void mapToReferenceFrame(ArrayRefPoint        &        refPoints,
                                    const ArrayPhysPoint &        physPoints,
                                    const ArrayCell      &        cellWorkset,
                                    const Teuchos::RCP< Basis< Scalar, FieldContainer<Scalar> > > HGRAD_Basis,
                                    const int &                   whichCell = -1);
    
    
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
    template<class ArrayRefPoint, class ArrayInitGuess, class ArrayPhysPoint, class ArrayCell>
    static void mapToReferenceFrameInitGuess(ArrayRefPoint        &        refPoints,
                                             const ArrayInitGuess &        initGuess,
                                             const ArrayPhysPoint &        physPoints,
                                             const ArrayCell      &        cellWorkset,
                                             const shards::CellTopology &  cellTopo,
                                             const int &                   whichCell = -1);
    
    template<class ArrayRefPoint, class ArrayInitGuess, class ArrayPhysPoint, class ArrayCell>
    static void mapToReferenceFrameInitGuess(ArrayRefPoint        &        refPoints,
                                             const ArrayInitGuess &        initGuess,
                                             const ArrayPhysPoint &        physPoints,
                                             const ArrayCell      &        cellWorkset,
                                             const Teuchos::RCP< Basis< Scalar, FieldContainer<Scalar> > > HGRAD_Basis,
                                             const int &                   whichCell = -1);
    

    
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
    template<class ArraySubcellPoint, class ArrayParamPoint>
    static void mapToReferenceSubcell(ArraySubcellPoint     &       refSubcellPoints,
                                      const ArrayParamPoint &       paramPoints,
                                      const int                     subcellDim,
                                      const int                     subcellOrd,
                                      const shards::CellTopology &  parentCell);
    
    
    
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
    template<class ArrayEdgeTangent>
    static void getReferenceEdgeTangent(ArrayEdgeTangent &            refEdgeTangent,
                                        const int &                   edgeOrd,
                                        const shards::CellTopology &  parentCell);

    
    
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
    template<class ArrayFaceTangentU, class ArrayFaceTangentV>
    static void getReferenceFaceTangents(ArrayFaceTangentU &           refFaceTanU,
                                         ArrayFaceTangentV &           refFaceTanV,
                                         const int &                   faceOrd,
                                         const shards::CellTopology &  parentCell);
    
    
    
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
                CellTools<Scalar>::getReferenceFaceNormal and these two methods are completely equivalent.
              - For 2D cells the reference side normal is defined by \f$\hat{{\bf n}}= \hat{\bf t}^\perp = (t_2,-t_1)\f$
                where \f$\hat{{\bf t}}=(t_1,t_2)\f$ is the tangent vector computed by 
                CellTools<Scalar>::getReferenceEdgeTangent. Therefore, the pair 
                \f$(\hat{{\bf n}},\hat{{\bf t}})\f$ is positively oriented.

        \param  refSideNormal     [out] - rank-1 array (D) with (constant) side normal
        \param  sideOrd           [in]  - ordinal of the side whose normal is computed
        \param  parentCell        [in]  - cell topology of the parent reference cell
*/
    template<class ArraySideNormal>
    static void getReferenceSideNormal(ArraySideNormal &             refSideNormal,
                                       const int &                   sideOrd,
                                       const shards::CellTopology &  parentCell);

    
    
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
                The method CellTools::getReferenceFaceTangents computes the reference face tangents
                \f${\partial\hat{\Phi}_{i}/\partial u}\f$ and \f${\partial\hat{\Phi}_{i}/\partial v}\f$.
 
        \param  refFaceNormal     [out] - rank-1 array (D) with (constant) face normal
        \param  faceOrd           [in]  - ordinal of the face whose normal is computed
        \param  parentCell        [in]  - cell topology of the parent reference cell
      */
    template<class ArrayFaceNormal>
    static void getReferenceFaceNormal(ArrayFaceNormal &             refFaceNormal,
                                       const int &                   faceOrd,
                                       const shards::CellTopology &  parentCell);
    
    
    
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
                \f$\hat{\mathcal E}_i\f$; see CellTools<Scalar>::getReferenceEdgeTangent that has the 
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
    template<class ArrayEdgeTangent, class ArrayJac>
    static void getPhysicalEdgeTangents(ArrayEdgeTangent &            edgeTangents,
                                        const ArrayJac &              worksetJacobians,
                                        const int &                   worksetEdgeOrd,
                                        const shards::CellTopology &  parentCell);
    
  /*  template<class ArrayEdgeTangent, class ArrayJac>
    static void getPhysicalEdgeTangentsTemp(ArrayEdgeTangent &            edgeTangents,
                                        const ArrayJac &              worksetJacobians,
                                        const int &                   worksetEdgeOrd,
                                        const shards::CellTopology &  parentCell);    */
    
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
                tangents on reference face \f$\hat{\mathcal F}_i\f$; see CellTools<Scalar>::getReferenceFaceTangents; 
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
    template<class ArrayFaceTangentU, class ArrayFaceTangentV, class ArrayJac>
    static void getPhysicalFaceTangents(ArrayFaceTangentU &           faceTanU,
                                        ArrayFaceTangentV &           faceTanV,
                                        const ArrayJac &              worksetJacobians,
                                        const int &                   worksetFaceOrd,
                                        const shards::CellTopology &  parentCell);
    
 /*   template<class ArrayFaceTangentU, class ArrayFaceTangentV, class ArrayJac>
    static void getPhysicalFaceTangentsTemp(ArrayFaceTangentU &           faceTanU,
                                        ArrayFaceTangentV &           faceTanV,
                                        const ArrayJac &              worksetJacobians,
                                        const int &                   worksetFaceOrd,
                                        const shards::CellTopology &  parentCell);    
    */
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
                \f$\hat{\mathcal S}_i\f$; see CellTools<Scalar>::getReferenceEdgeTangent, that has the 
                same local ordinal as the sides in the workset;
        \li     For 3D cells: \f$ {\partial\hat{\Phi}_i/\partial u}, {\partial\hat{\Phi}_i/\partial v}\f$ are the (constant) 
                tangents on reference side (face) \f$\hat{\mathcal S}_i\f$; see CellTools<Scalar>::getReferenceFaceTangents,
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
                CellTools<Scalar>::getPhysicalFaceNormals and these two methods are completely equivalent.
              - For 2D cells the physical side normals are defined by \f${\bf n}=(t_2,-t_1)\f$
                where \f${\bf t}=(t_1,t_2)\f$ are the physical edge tangents computed by 
                CellTools<Scalar>::getPhysicalEdgeTangents. Therefore, the pairs \f$({\bf n},{\bf t})\f$ are positively oriented.


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
    template<class ArraySideNormal, class ArrayJac>
    static void getPhysicalSideNormals(ArraySideNormal &             sideNormals,
                                       const ArrayJac &              worksetJacobians,
                                       const int &                   worksetSideOrd,
                                       const shards::CellTopology &  parentCell);
    
    
    
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
                tangents on reference face \f$\hat{\mathcal F}_i\f$; see CellTools<Scalar>::getReferenceFaceTangents;
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
    template<class ArrayFaceNormal, class ArrayJac>
    static void getPhysicalFaceNormals(ArrayFaceNormal &             faceNormals,
                                       const ArrayJac &              worksetJacobians,
                                       const int &                   worksetFaceOrd,
                                       const shards::CellTopology &  parentCell);

  /*  template<class ArrayFaceNormal, class ArrayJac>
    static void getPhysicalFaceNormalsTemp(ArrayFaceNormal &             faceNormals,
                                       const ArrayJac &              worksetJacobians,
                                       const int &                   worksetFaceOrd,
                                       const shards::CellTopology &  parentCell);    
    */
    //============================================================================================//
    //                                                                                            //
    //                                        Inclusion tests                                     //
    //                                                                                            //
    //============================================================================================//
    
    /** \brief  Checks if a point belongs to a reference cell. 
      
                Requires cell topology with a reference cell.
      
        \param  point             [in]  - spatial coordinates of the point tested for inclusion
        \param  pointDim          [in]  - spatial dimension of that point
        \param  cellTopo          [in]  - cell topology of the cells stored in \c cellWorkset
        \param  threshold         [in]  - "tightness" of the inclusion test
        \return 1 if the point is in the closure of the specified reference cell and 0 otherwise.
    */ 
    static int checkPointInclusion(const Scalar*                 point,
                                   const int                     pointDim,
                                   const shards::CellTopology &  cellTopo,
                                   const double &                threshold = INTREPID_THRESHOLD);
    
    
    
    /** \brief  Checks if a set of points belongs to a reference cell. 
      
                Requires cell topology with a reference cell. See Intrepid::CellTools::checkPointwiseInclusion 
                for admissible ranks and dimensions of the input point array.
      
        \param  points            [in]  - rank-1, 2 or 3 array (point, vector of points, matrix of points)
        \param  cellTopo          [in]  - cell topology of the cells stored in \c cellWorkset
        \param  threshold         [in]  - "tightness" of the inclusion test
        \return 1 if all points are in the closure of the specified reference cell
                0 if at least one point is outside the closure of the reference cell 
                
    */
    template<class ArrayPoint>
    static int checkPointsetInclusion(const ArrayPoint &            points,
                                      const shards::CellTopology &  cellTopo, 
                                      const double &                threshold = INTREPID_THRESHOLD);
    
    
    
    /** \brief  Checks every point in a set for inclusion in a reference cell.  
      
                Requires cell topology with a reference cell. Admissible ranks and dimensions of the 
                input point array and the corresponding rank and dimension of the output array are as follows:
                \verbatim
                |-------------------|-------------|-------------|-------------|
                |  rank: (in)/(out) |    1/1      |     2/1     |    3/2      |
                |-------------------|-------------|-------------|-------------|
                |  points    (in)   |     (D)     |    (I, D)   |  (I, J, D)  |
                |-------------------|-------------|-------------|-------------|
                |  inRefCell (out)  |     (1)     |    (I)      |  (I, J)     |
                |------------------ |-------------|-------------|-------------|
                \endverbatim
                Example: if \c points is rank-3 array with dimensions (I, J, D), then
        \f[
               \mbox{inRefCell}(i,j) = 
                 \left\{\begin{array}{rl} 
                    1 & \mbox{if $points(i,j,*)\in\hat{\mathcal{C}}$} \\[2ex]
                    0 & \mbox{if $points(i,j,*)\notin\hat{\mathcal{C}}$} 
                 \end{array}\right.
          \f]
        \param  inRefCell         [out] - rank-1 or 2 array with results from the pointwise inclusion test
        \param  refPoints         [in]  - rank-1,2 or 3 array (point, vector of points, matrix of points)
        \param  cellTopo          [in]  - cell topology of the cells stored in \c cellWorkset
        \param  threshold         [in]  - "tightness" of the inclusion test
    */
    template<class ArrayIncl, class ArrayPoint>
    static void checkPointwiseInclusion(ArrayIncl &                   inRefCell,
                                        const ArrayPoint &            points,
                                        const shards::CellTopology &  cellTopo, 
                                        const double &                threshold = INTREPID_THRESHOLD);
    
    
    
    /** \brief  Checks every point in a set or multiple sets for inclusion in physical cells from a cell workset.
      
                There are two use cases:
        \li     Checks every point from a \b single point set, stored in a rank-2 array (P,D) for inclusion
                in a \b specified physical cell \f${\mathcal C}\f$ from a cell workset.
        \li     Checks every point from \b multiple point sets indexed by a cell ordinal, and stored in a rank-3
                (C,P,D) array, for inclusion in the physical cell having the same cell ordinal, for \b all 
                cells in a cell workset. 
      
                For a single point set in a rank-2 array (P,D) and \c whichCell set to a valid cell ordinal
                relative to \c cellWorkset returns a rank-1 (P) array such that
        \f[
                \mbox{inCell}(p) = 
                  \left\{\begin{array}{rl} 
                     1 & \mbox{if $points(p,*)\in {\mathcal{C}}$} \\ [2ex]
                     0 & \mbox{if $points(p,*)\notin {\mathcal{C}}$} 
                  \end{array}\right.
        \f]
                For multiple point sets in a rank-3 array (C,P,D) and \c whichCell=-1 (default value)
                returns a rank-2 (C,P) array such that
        \f[
                \mbox{inCell}(c,p) = 
                  \left\{\begin{array}{rl} 
                      1 & \mbox{if $points(c,p,*)\in {\mathcal{C}}$} \\ [2ex]
                      0 & \mbox{if $points(c,p,*)\notin {\mathcal{C}}$} 
                \end{array}\right.
        \f]

        \param  inCell            [out] - rank-1  array with results from the pointwise inclusion test
        \param  points            [in]  - rank-2 array with dimensions (P,D) with the physical points
        \param  cellWorkset       [in]  - rank-3 array with dimensions (C,N,D) with the nodes of the cell workset
        \param  cellTopo          [in]  - cell topology of the cells stored in \c cellWorkset
        \param  whichCell         [in]  - ordinal of the cell used in the inclusion test
        \param  threshold         [in]  - tolerance for inclusion tests on the input points
      */
    template<class ArrayIncl, class ArrayPoint, class ArrayCell>
    static void checkPointwiseInclusion(ArrayIncl &                   inCell,
                                        const ArrayPoint &            points,
                                        const ArrayCell &             cellWorkset,
                                        const shards::CellTopology &  cell,
                                        const int &                   whichCell = -1, 
                                        const double &                threshold = INTREPID_THRESHOLD);
    
    
    
    /** \brief  Retrieves the Cartesian coordinates of a reference cell vertex. 
      
                Requires cell topology with a reference cell. Vertex coordinates are always returned 
                as an (x,y,z)-triple regardlesss of the actual topological cell dimension. The unused 
                coordinates are set to zero, e.g., vertex 0 of Line<2> is returned as {-1,0,0}.                
                
        \param  cell              [in]  - cell topology of the cell
        \param  vertexOrd         [in]  - vertex ordinal 
        \return pointer to array with the (x,y,z) coordinates of the specified reference vertex
      */
    static const double* getReferenceVertex(const shards::CellTopology& cell,
                                            const int                   vertexOrd);
    
    
    
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
    template<class ArraySubcellVert>
    static void getReferenceSubcellVertices(ArraySubcellVert &          subcellVertices,
                                            const int                   subcellDim,
                                            const int                   subcellOrd,
                                            const shards::CellTopology& parentCell);
    

    
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
    static const double* getReferenceNode(const shards::CellTopology& cell,
                                          const int                   nodeOrd);
    
    
    
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
    template<class ArraySubcellNode>
    static void getReferenceSubcellNodes(ArraySubcellNode&          subcellNodes,
                                        const int                   subcellDim,
                                        const int                   subcellOrd,
                                        const shards::CellTopology& parentCell);
    
   
    
    /** \brief  Checks if a cell topology has reference cell
        \param  cell              [in]  - cell topology
        \return 1 if the cell topology has reference cell, 
                0 oterwise
      */
    static int hasReferenceCell(const shards::CellTopology &  cellTopo);
    
    
    
    //============================================================================================//
    //                                                                                            //
    //                                           Debug                                            //
    //                                                                                            //
    //============================================================================================//
    
    
    /** \brief  Prints the reference space coordinates of the vertices of the specified subcell
        \param  subcellDim        [in]  - dimension of the subcell where points are mapped to
        \param  subcellOrd        [in]  - subcell ordinal
        \param  parentCell        [in]  - cell topology of the parent cell.
    */
    static void printSubcellVertices(const int subcellDim,
                                     const int subcellOrd,
                                     const shards::CellTopology & parentCell);
    
    
    
    /** \brief  Prints the nodes of a subcell from a cell workset 
      
      */
    template<class ArrayCell>
    static void printWorksetSubcell(const ArrayCell &             cellWorkset,
                                    const shards::CellTopology &  parentCell,
                                    const int&                    pCellOrd,
                                    const int&                    subcellDim,
                                    const int&                    subcellOrd,
                                    const int&                    fieldWidth = 3);

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
  template<class ArrayCVCoord, class ArrayCellCoord>
  static void getSubCVCoords(ArrayCVCoord & subCVcoords, const ArrayCellCoord & cellCoords,
                             const shards::CellTopology& primaryCell);

  /** \brief Compute cell barycenters. 

    \param barycenter      [out] - array containing cell baycenters
    \param cellCoords       [in] - array containing cell coordinates 

 */
  template<class ArrayCent, class ArrayCellCoord>
  static void getBarycenter(ArrayCent & barycenter, const ArrayCellCoord & cellCoords);

    
  }; // class CellTools

} // namespace Intrepid

// include templated function definitions
#include "Intrepid_CellToolsDef.hpp"

#endif

/***************************************************************************************************
 **                                                                                               **
 **                           D O C U M E N T A T I O N   P A G E S                               **
 **                                                                                               **
 **************************************************************************************************/

/**  
\page cell_tools_page                 Cell tools

<b>Table of contents </b>

-   \ref cell_topology_sec
-   \ref cell_topology_ref_cells
  - \ref sec_cell_topology_ref_map
  - \ref sec_cell_topology_ref_map_DF
  - \ref sec_cell_topology_subcell_map
  - \ref sec_cell_topology_subcell_wset

\section cell_topology_sec            Cell topologies

The range of admissible cell shapes in Intrepid is restricted to <var>d</var>-dimensional 
<strong>polytopes</strong>, <var>d=1,2,3</var>. A polytope is defined by a set of its vertices 
\f$ \{ {\bf v}_0,\ldots {\bf v}_V\} \f$ and a <strong>base topology</strong>  <var>BT</var> that 
defines how these verices are connected into <var>k</var>-dimensional, <var>k < d</var> facets 
(k-subcells) of that polytope.  

The base topology of any polytope can be extended by augmenting the set of its vertices by 
an additional set of points \f$\{ {\bf p}_0,\ldots {\bf p}_P\} \f$. The <strong>extended topology</strong>
<var>ET</var> is defined by specifying the connectivity of the set
\f$ \{ {\bf v}_0,\ldots {\bf v}_V\}\cup \{ {\bf p}_0,\ldots {\bf p}_P\} \f$
relative to the subcells specified by its base topology <var>BT</var>.
The vertices and the extra points are collectively referred to as <strong>nodes</strong>. Thus,
a polytope with <strong>extended topology</strong> <var>ET</var> is defined by a set of nodes 
\f$\{{\bf q}_0,\ldots,{\bf q}_N\}\f$, where \f$N = V + P\f$, and a connectivity rule for these nodes.

Intrepid requires any cell to have a valid base topology. The nodes of the cell should always be 
ordered by listing its vertices <strong>first</strong>, i.e.,
\f[ 
  \{{\bf q}_0,\ldots,{\bf q}_N\} = \{ {\bf v}_0,\ldots {\bf v}_V,{\bf p}_0,\ldots {\bf p}_P\} 
  \f]
To manage cell topologies Intrepid uses the Shards package http://trilinos.sandia.gov/packages/shards .
Shards provides definitions for a standard set of base and extended cell topologies plus tools to
construct custom, user defined cell topologies, such as arbitrary polyhedral cells. For further
details see Shards documentation. 



\section cell_topology_ref_cells      Reference cells

For some cell topologies there exist simple, e.g., polynomial, mappings that allow to obtain any 
cell having that topology as an image of a single "standard" cell. We refer to such standard cells
as <strong>reference</strong> cells. 

Just like in the general case, a reference cell with a base topology <var>BT</var> is defined by a 
set of vertices, and a reference cell with extended topology <var>ET</var> is defined by a set of 
nodes that include the original vertices and some additional points. 

The actual vertex and node coordinates for the reference cells can be chosen arbitrarily; however, 
once selected they should not be changed because in many cases, e.g., in finite element reconstructions,
all calculations are done on a reference cell and then transformed to physical cells by an appropriate
pullback (see Section \ref sec_pullbacks).

In Intrepid base and extended reference cell topologies are defined using the following selections
of vertex and node coordinates:

\verbatim
|=======================|==============================================================================|
| Topology family       |    reference cell vertices/additional nodes defining extended topology       |
|=======================|==============================================================================|
| Line<2>               |                                                                              |
| Beam<2>               | {(-1,0,0),(1,0,0)}                                                           |
| ShellLine<2>          |                                                                              |
|-----------------------|------------------------------------------------------------------------------|
| Line<3>               |                                                                              |
| Beam<3>               | {(0,0,0)}                                                                    |
| ShellLine<3>          |                                                                              |
|=======================|==============================================================================|
| Triangle<3>           | {(0,0,0),(1,0,0),(0,1,0)}                                                    |
| ShellTriangle<3>      |                                                                              |
|-----------------------|------------------------------------------------------------------------------|
| Triangle<4>           | {(1/3,1/3,0)}                                                                |
|.......................|..............................................................................|
| Triangle<6>           | {(1/2,0,0),(1/2,1/2,0),(0,1/2,0)}                                            |
| ShellTriangle<6>      |                                                                              |
|=======================|==============================================================================|
| Quadrilateral<4>      | {(-1,-1,0),(1,-1,0), (1,1,0),(-1,1,0)}                                       |
| ShellQuadrilateral<4> |                                                                              |
|-----------------------|------------------------------------------------------------------------------|
| Quadrilateral<8>      | {(0,-1,0),(1,0,0),(0,1,0),(-1,0,0)}                                          |
| ShellQuadrilateral<8> |                                                                              |
|.......................|..............................................................................|
| Quadrilateral<9>      | {(0,-1,0),(1,0,0),(0,1,0),(-1,0,0),(0,0,0)}                                  |
| ShellQuadrilateral<9> |                                                                              |
|=======================|==============================================================================|
| Tetrahedron<4>        | {(0,0,0),(1,0,0),(0,1,0),(0,0,1)}                                            |
|-----------------------|------------------------------------------------------------------------------|
| Tetrahedron<8>        | {(1/2,0,0),(1/2,1/2,0),(0,1/2,0),(1/3,1/3,1/3)}                              |
| Tetrahedron<10>       | {(1/2,0,0),(1/2,1/2,0),(0,1/2,0),(0,0,1/2),(1/2,0,1/2),(0,1/2,1/2)}          |
|=======================|==============================================================================|
| Pyramid<5>            | {(-1,-1,0),(1,-1,0),(1,1,0),(-1,1,0),(0,0,1)}                                |
|-----------------------|------------------------------------------------------------------------------|
| Pyramid<13>           | {(0,-1,0),(1,0,0),(0,1,0),(-1,0,0), 1/2((-1,-1,1),(1,-1,1),(1,1,1),(-1,1,1))}|
| Pyramid<14>           | all of the above and (0,0,0)                                                 |
|=======================|==============================================================================|
| Wedge<6>              | {(0,0,-1),(1,0,-1),(0,1,-1),(0,0,1),(1,0,1),(0,1,1)}                         |
|-----------------------|------------------------------------------------------------------------------|
| Wedge<15>             | {(1/2,0,-1),(1/2,1/2,-1),(0,1/2,-1), (0,0,0),(1,0,0),(0,1,0),                |  
|                       |  (1/2,0, 1),(1/2,1/2, 1),(0,1/2, 1)                                          |
|.......................|..............................................................................|
| Wedge<18>             | All of the above plus {(1/2,0,0),(1/2,1/2,0),(0,1/2,0)}                      |
|=======================|==============================================================================|
| Hexahedron<8>         | {(-1,-1,-1),(1,-1,-1),(1,1,-1),(-1,1,-1),(-1,-1,1),(1,-1,1),(1,1,1),(-1,1,1)}|
|-----------------------|------------------------------------------------------------------------------|
| Hexahedron<20>        | {(0,-1,-1),(1,0,-1),(0,1,-1),(-1,0,-1), (0,-1,0),(1,0,0),(0,1,0),(-1,0,0),   |
|                       |  (0,-1, 1),(1,0, 1),(0,1, 1),(-1,0, 1) }                                     |
|.......................|..............................................................................|
| Hexahedron<27>        | All of the above plus center point and face midpoints:                       |
|                       | {(0,0,0), (0,0,-1),(0,0,1), (-1,0,0),(1,0,0), (0,-1,0),(0,1,0)}              |
|=======================|==============================================================================|
\endverbatim

Finite element reconstruction methods based on pullbacks (see Section \ref sec_pullbacks) are 
restricted to the above cell topologies.  


  
\subsection sec_cell_topology_ref_map      Reference-to-physical cell mapping

The mapping that takes a given reference cell to a physical cell with the same topology is defined 
using a nodal Lagrangian basis corresponding to the nodes of the reference cell. In other words, the
mapping is constructed using basis functions that are dual to the nodes of the reference cell.
Implementation details are as follows.

Assume that \f$ \hat{\mathcal{C}} \f$ is a reference cell with topology <var>T</var> and nodes
\f$\{\hat{{\bf q}}_0,\ldots,\hat{{\bf q}}_{N}\}\f$, and that \f$ \{\hat{\phi}_i\}_{i=0}^{N} \f$ is
the Lagrangian basis dual to these nodes, i.e., \f$ \hat{\phi}_i( \hat{{\bf q}}_j)=\delta_{ij} \f$.
A physical cell \f$\mathcal{C}\f$ with the same topology <var>T</var> as \f$\hat{\mathcal{C}}\f$ is  
then defined as the image of \f$ \hat{\mathcal{C}} \f$ under the mapping
\f[               
  F_\mathcal{C}(\hat{\bf x}) = \sum_{m=0}^{N}  {\bf q}_m(\mathcal{C})  \hat{\phi}_m(\hat{\bf x})          
\f]
where \f$\{{\bf q}_0(\mathcal{C}),\ldots,{\bf q}_N(\mathcal{C})\}\f$ is the set of <strong>physical nodes</strong>
that defines \f$\mathcal{C}\f$. The number of physical nodes is required to match the number of reference  
nodes in the specified cell topology <var>T</var>. The <var>i</var>-th coordinate function of the 
reference-to-physical mapping is given by
\f[        
  \big(F_{\mathcal{C}}(\hat{\bf x})\big)_i = 
              \sum_{m=0}^{N}  \big({\bf q}_m(\mathcal{C})\big)_i\hat{\phi}_m(\hat{\bf x})           
\f]
where \f$ \big({\bf q}_m(\mathcal{C})\big)_i \f$ is the <var>i</var>-th spatial coordinate of the <var>m</var>-th node.

For simplicity, unless there's a chance for confusion, the cell symbol will be ommitted from the
designations of physical points and reference-to-physical maps, i.e., we shall simply write
\f$F(\hat{\bf x})\ \mbox{and}\ {\bf q}_m\f$.

\par Summary
  
\li      \f$F(\hat{\bf x})\f$: implemented in Intrepid::CellTools::mapToPhysicalFrame 
\li      \f$F^{-1}({\bf x})\f$: implemented in Intrepid::CellTools::mapToReferenceFrame 
  
\warning Intrepid::CellTools does not check for non-degeneracy of the physical cell obtained from a
given set of physical nodes. As a result, <var>F</var> is not guaranteed to be a diffeomorphism,
i.e., it may not have a continuously differentiable inverse. In this case some 
Intrepid::CellTools methods, such as Intrepid::CellTools::setJacobianInv, 
and Intrepid::CellTools::mapToReferenceFrame will fail.
    

  
\subsection sec_cell_topology_ref_map_DF   Jacobian of the reference-to-physical cell mapping

Intrepid follows the convention that the rows of the Jacobian are the transposed gradients of the
coordinate functions of the mapping, i.e.,
\f[                    
      DF_{ij}(\hat{{\bf x}}) = \frac{\partial F_i(\hat{{\bf x}})}{\partial\hat{{\bf x}}_j}                    
\f]
In light of the definition of <var>F</var> in Section \ref sec_cell_topology_ref_map, it follows that
\f[  
      DF_{ij}(\hat{{\bf x}}) = \sum_{m=0}^{N}
                ({\bf q}_m)_i\frac{\partial\hat{\phi}_m(\hat{{\bf x}})}{\partial\hat{{\bf x}}_j} \,.       
\f]

\par Summary
  
\li     \f$DF_{ij}(\hat{{\bf x}}) \f$: implemented in Intrepid::CellTools::setJacobian
\li     \f$DF_{ij}^{-1}(\hat{{\bf x}}) \f$: implemented in Intrepid::CellTools::setJacobianInv
\li     \f$\mbox{det} DF_{ij}(\hat{{\bf x}}) \f$: implemented in  Intrepid::CellTools::setJacobianDet                      


  
\subsection sec_cell_topology_subcell_map  Parametrization of physical 1- and 2-subcells
  
Parametrization of a given physical k-subcell \f$\mathcal{S}_i\f$, k=1,2, is a map from a 
k-dimensional parametrization domain \e R to that subcell:
\f[
    \Phi : R \mapsto \mathcal{S} \,.
\f]
Parametrization domains play role similar to that of reference cells in the sense that they allow 
computation of line and surface integrals on 1- and 2-subcells (edges and faces) to be reduced to   
computation of integrals on \e R . 

Parametrization maps are supported for 1- and 2-subcells (edges and faces) that belong to physical 
cells with reference cells. The reason is that these maps are defined by the composition of  the
parametrization maps for reference edges and faces with the mapping \e F defined in \ref sec_cell_topology_ref_map.
As a result, parametrization of a given physical k-subcell requires selection of a 
<strong>parent cell</strong> that contains the subcell. 

\remark  
Because a given k-subcell may belong to more than one physical cell, its parent cell is not unique. For a single
k-subcell the choice of a parent cell is not important, however, when dealing with subcell worksets
parent cells must all have the same topology (see \ref sec_cell_topology_subcell_wset for details about
subcell worksets).
  
  
Implementation of subcell parametrization is as follows. Assume that \f$\mathcal{S}_i\f$ is a k-subcell   
with parent cell \f$\mathcal{C}\f$; \f$\mathcal{\hat{C}}\f$ is the associated reference cell and 
\e i is the local ordinal of the subcell relative to the reference cell. To this physical subcell
corresponds a reference subcell \f$\hat{\mathcal{S}}_i\f$ having the same local ordinal. Parametrization of
the reference k-subcell is a map from the k-dimensional parametrization domain \e R to that subcell:
\f[
    \hat{\Phi} : R \mapsto \hat{\mathcal{S}}_i \,.
\f]
Parametrization of \f$\mathcal{S}_i\f$ is then defined as
\f[
    \Phi = F\circ \hat{\Phi} \,,
\f]
where \e F is the reference-to-physical mapping between the parent cell and its reference cell. 
  

A 1-subcell (edge) always has \c Line<N> topology and so, the parametrization domain for edges is the standard 1-cube:
\f[
      R = [-1,1]\,.
\f]
On the other hand, faces of reference cells can have \c Triangle<N> and/or \c Quadrilateral<N> topologies. 
Thus, the parametrization domain for a 2-subcell depends on its topology and is either the standard 
2-simplex or the standard 2-cube:
\f[
      R = \left\{\begin{array}{rl} 
          \{(0,0),(1,0),(0,1)\} & \mbox{if the face is Triangle} \\[1ex]
            [-1,1]\times [-1,1] & \mbox{if the face is Quadrilateral}
          \end{array}\right.
\f]
\par Summary

-    \f$ \Phi : R \mapsto {\mathcal{S}}_i \f$ requires two steps: 
  -# Intrepid::CellTools::mapToReferenceSubcell to apply \f$\hat{\Phi}: R \mapsto \hat{\mathcal{S}}_i\f$;
  -# Intrepid::CellTools::mapToPhysicalFrame to apply \e F



\subsection sec_cell_topology_subcell_wset Subcell worksets

A subcell workset comprises of 1- or 2-subcells and associated parent cells that satisfy the 
following conditions

\li     all subcells have the same cell topology;
\li     all parent cells have the same cell topology;
\li     The parent cell topology has a reference cell;
\li     relative to that reference cell, all subcells in the workset have the same local ordinal

Therefore, a subcell workset is defined by 

-#      collecting a set of 1- or 2-subcells having the same topology
-#      selecting a parent cell for every subcell in such a way that 
  -#      all parent cells have the same cell topology
  -#      all subcells in the workset have the same local ordinal relative to the parent cell topology

  Obviously, a subcell can have multiple parent cells. For example, in a mesh consisiting of Triangle<3> cells, every
  edge is shared by 2 triangle cells. To define an edge workset we can use either one of the two traingles
  sharing the cell.

  Suppose now that the mesh comprises of Triangle<3> and Quadrilateral<4> cells and we want to define an
  edge workset. Let's say the first few edges in our workset happen to be shared by 2 triangles and so 
  we select one of them as the parent cell. Now suppose the next edge is shared by a traingle and a 
  quadrilateral. Because all parent cells in the workset must have the same cell topology we cannot use
  the quadrilateral as a parent cell and so we choose the triangle. Finally suppose that one of the candidate 
  edges for our workset is shared by 2 quadrilaterals. Because of the requirement that all parent cells
  have the same topology, we will have to reject this edge because it does not posses a potential parent cell
  with the same topology as the rest of the edges in our workset.

  A subcell workset is denoted by \f$ \{\mathcal{S}_{c,i}\}_{c=0}^N \f$, where

  \li <var>c</var> is parent cell ordinal;
  \li <var>i</var> is the local subcell ordinal (relative to the topology of the parent cell) shared
  by all subcells in the workset. 

  */
#endif

#if defined(Intrepid_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Intrepid package is deprecated"
#endif
#endif

