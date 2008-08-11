// @HEADER
// ************************************************************************
//
//                           Intrepid Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Pavel Bochev (pbboche@sandia.gov) or
//                    Denis Ridzal (dridzal@sandia.gov).
//
// ************************************************************************
// @HEADER

/** \file   Intrepid_MultiCell.hpp
\brief  Header file for the Intrepid::MultiCell class.
\author Created by P. Bochev and D. Ridzal.
*/

#ifndef INTREPID_MULTICELL_HPP
#define INTREPID_MULTICELL_HPP

#include "Intrepid_ConfigDefs.hpp"
#include "Intrepid_RealSpace.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_TestForException.hpp"

namespace Intrepid {
  
  /** \relates Intrepid::MultiCell
  \brief Provides additional shape points to compute a higher degree chart for a given physical 
  cell. Can only be used with cells that have a reference cell. 
  
  The default chart of a physical cell has degree 1 and is computed using only its vertices. Such a
  chart may not be adequate for cells that have curved sides. To approximate such cells Intrepid 
  provides charts of degree 2 and 3. Their computation requires additional shape points that live
  on some or all of the subcells of the cell. The precise allocation of these points depends on the
  degree of the chart and the cell shape. To use the optional higher degree charts the user must 
  populate the data of this struct with a valid set of shape points as follows:
  
  shapePoints_[0] -> 1D, 2D, 3D: array of edge shape points ordered by local edge ID
  shapePoints_[1] -> 3D: array of face shape points ordered by local face ID
           -> 2D: array of cell shape points
  shapePoints_[2] -> 3D: array of cell shape points
  
  The valid shape shapePoints_ sets for charts of degree 2 and 3 are
  \verbatim
  |==========|==========================================|===========================================|
  | cell type|            chart of degree 2             |             chart of degree 3             |
  |==========|==========================================|===========================================|
  | EDGE     |  shapePoints_[0] -> 1 pt.                |  shapePoints_[0] -> 2 points              |
  |==========|==========================================|===========================================|
  | QUAD     |  shapePoints_[0] -> 4 pts. (1 per edge)  |  shapePoints_[0] -> 8 pts.  (2 per edge)  |
  |          |  shapePoints_[1] -> 1 pt                 |  shapePoints_[1] -> 4 pts.                |
  |==========|==========================================|===========================================|
  |          |  shapePoints_[0] -> 12 pts. (1 per edge) | shapePoints_[0] -> 24 pts.  (2 per edge)  |
  | HEX      |  shapePoints_[1] -> 6 pts.  (1 per face) | shapePoints_[1] -> 24 pts.  (4 per face)  |
  |          |  shapePoints_[2] -> 1 pt.                | shapePoints_[2] -> 8 pts.                 |
  |==========|==========================================|===========================================|
  |  TRI     |  shapePoints_[0] -> 3 pts.  (1 per edge) |  shapePoints_[0] -> 6 pts.  (2 per edge)  |
  |          |                                          |  shapePoints_[1] -> 1 pt.                 |
  |==========|==========================================|===========================================|
  |  TET     |  shapePoints_[0] -> 6 pts.  (1 per edge) |  shapePoints_[0] -> 12 pts. (2 per edge)  |
  |          |                                          |  shapePoints_[1] -> 4 pts.  (1 per face)  |
  |==========|==========================================|===========================================|
  | PYRAMID  |
  |==========|==========================================|===========================================|
  | TRIPRISM |
  |==========|==========================================|===========================================|
  \endverbatim
  
  To define the chart of a cell Intrepid uses the standard nodal basis functions of appropriate 
  degree for the particular cell shape. However, for efficiency, instead of using 
  */
  template<class Scalar>
  struct ShapePoints {
    Teuchos::Array< Teuchos::Array< Point<Scalar> > > shapePoints_;
    int chartDegree_;
  };
  
  
/** \class Intrepid::MultiCell
    \brief A MultiCell (batch or group of cells) object is used to communicate cell information from a 
    (global) mesh object to a (local) interpolation/reconstruction operator. MultiCell interface
    provides methods to access global data by local cell connectivity templates and to obtain
    geometrical information about the cells such as cell Jacobians, cell charts (transformations)
    end etc. 
  
    A MultiCell object allows the user to aggregate cells that use the exact same reconstruction
    operator into groups. This may improve performance when computing local operators for the cells
    in the batch.
  
    \todo  Carefully define all terms (cell, canonical cell, subcell, dimension, orientation), and 
    describe the Intrepid::MultiCell interface. Also discuss creation of custom cell shapes.
*/
template<class Scalar>
class MultiCell {
  private:
    
    /** \brief Stores all canonical cell templates (see Intrepid_CellTemplates.hpp for their definition)
    */
    static const ConnMapTemplate connMapCanonical_[CELL_CANONICAL_MAX][3];
    
    
    /** \brief Provides storage for user-defined connectivity templates.
    */
    static ConnMapTemplate connMapCustom_[CELL_MAX - CELL_CANONICAL_MAX - 1][3];
    
    
    /** \brief Names of canonical and custom cell templates (corresponds to the enumeration in ECell)
    */
    static const char *cellNames_[];
    
    
    /** \brief Number of cells in the multicell.
    */
    int numCells_;
    
    
    /** \brief Type of the generating cell.
    */
    ECell myCellType_;
    
    
    /** \brief Rank-two array containing the vertices of the cells in the MultiCell as Point objects.
      The leading dimension equals the number of cells in the MultiCell; the number of columns equals
      the number of vertices of the generating cell type. In sum, vertices_[i][j] is the j-th vertex 
      of the i-th cell in the MultiCell.
    */
    Teuchos::Array< Teuchos::Array<Point<Scalar> > >  vertices_;
    
    
    /** \brief Rank-two array containing the edge signs of the cells in the MultiCell. The leading 
      dimension equals the number of cells in the MultiCell and the number of columns equals the 
      number of edges (1-subcells) of the generating cell type. In sum, edgeSigns_[i][j] is the sign
      of the j-th edge of the i-th cell. Upon construction of a MultiCell the size of this array is 
      zero thereby indicating that edge signs are undefined. Edge signs can only be defined in 2D 
      and 3D where they have the following meaning:
      
      \verbatim
      edgeSigns_[i][j] = +1 if local edge direction coincides with the global edge direction
      edgeSigns_[i][j] = -1 if local edge direction is opposite to the global edge direction
      \endverbatim
      
      Local edge direction is defined by the vertex order in the cell template of the generating cell 
      type. For example, local directions of the edges in a QUAD cell are defined by the order of 
      their vertex pairs {{0,1}, {1,2}, {2,3}, {3,0}}, i.e., direction of edge 0 is from vertex 0 to
      vertex 1; direction of edge 1 is  from vertex 1 to vertex 2 and so on. Setting edge signs is
      user's responsibility because the MultiCell is not aware of the global mesh structures and 
      orientation choices in the user code.
    */
    Teuchos::Array< Teuchos::Array<short> > edgeSigns_;
    
    
    /** \brief Rank-two array containing the face signs of the cells in the MultiCell. The leading
      dimension equals the number of cells in the MultiCell and the number of columns equals the 
      number of "faces" (2-subcells in 3D and 1-subcells in 2D) of the generating cell type. In sum,
      faceSigns_[i][j] is the sign of the j-th face of the i-th cell. Upon construction of a MultiCell 
      the size of this array is zero thereby indicating that face signs are undefined. Face signs can 
      only be defined in 2D and 3D. While in 2D edges and "faces" are the same 1-subcells of the
      generating cell, separate "face" and edge signs are needed to distinguish between orientation
      of 1-subcells by the unit tangent (edge signs) vs. orientation by the unit normal ("face" signs)
      Face signs are defined as follows:
      
      \verbatim
      faceSigns_[i][j] = +1 if local unit normal coincides with the global unit normal
      faceSigns_[i][j] = -1 if local unit normal coincides with the global unit normal
      \endverbatim
      
      Local unit normals in 3D are defined using the right hand rule and the vertex order of the faces
      in the cell template of the generating cell type.  For example, the local unit normals of the
      faces in a PYRAMID are defined by the order of their vertices in the PYRAMID cell template:
      {{0,3,2,1}, {0,1,4}, {1,2,4}, {2,3,4}, {3,0,4}} and application of the right hand rule, i.e., 
      local orientation on all faces is provided by the outer unit normal. In 2D face signs are set
      by comparing local and global unit normals to 1-subcells. Setting face signs is
      user's responsibility because the MultiCell is not aware of the global mesh structures and 
      orientation choices in the user code.
    */
    Teuchos::Array< Teuchos::Array<short> > faceSigns_;
    
    
    /** \brief Rank-two array containing the edge tags of the cells in the MultiCell. The leading 
      dimension equals the number of cells in the MultiCell and the number of columns equals the 
      number of edges (1-subcells) of the generating cell type. In sum, edgeTags_[i][j] is the tag
      of the j-th edge of the i-th cell. Upon construction of a MultiCell the size of this array is 
      zero thereby indicating that edge tags are undefined. Intrepid uses tags to tell other methods 
      which edges they should work on. The defaults are 0 (skip this edge) and 1 (work on this edge). 
      Tags can be used to set computation of, e.g., Neumann data, or edge integrals along some interfaces.
    */
    Teuchos::Array< Teuchos::Array<short> > edgeTags_;
    
      
    /** \brief Rank-two array containing the face tags of the cells in the MultiCell. The leading
      dimension equals the number of cells in the MultiCell and the number of columns equals the 
      number of "faces" (2-subcells in 3D and 1-subcells in 2D) of the generating cell type. In sum,
      faceTags_[i][j] is the tag of the j-th face of the i-th cell. Upon construction of a MultiCell 
      the size of this array is zero thereby indicating that face signs are undefined. Intrepid uses 
      tags to tell other methods which faces they should work on. The defaults are 0 (skip this face) 
      and 1 (work on this face). While in 2D "faces" and edges are the same 1-subcells, separate
      face and edge tags may be needed to distingusih between, e.g., boundaries where normal and 
      tangential components of a vector field are specified.
    */
    Teuchos::Array< Teuchos::Array<short> > faceTags_;
    
    
    /** \brief Rank-one array of charts, i.e., mappings between the cells in the MultiCell and their 
      standard reference cell. Upon construction of a MultiCell the size of this array is zero thereby 
      indicating that the atlas is undefined. An atlas can be populated with charts if and only if
      the generating cell type has a reference cell, i.e., it is one of CELL_EDGE, CELL_QUAD, CELL_TRI, 
      CELL_TET, CELL_HEX, CELL_PYRAMID, or CELL_TRIPRISM.
    */
    Teuchos::Array< ChartTemplate<Scalar> >  atlas_;
    
    //--------------------------------------------------------------------------------------------//
    //                                                                                            //
    //      Private members of the MultiCell class which are used by LocalField classes to store  //
    //      information needed to compute integrals in operator and functional methods            //                                                                                      //
    //--------------------------------------------------------------------------------------------//
    
    
    /** \brief A four-dimensional array of Matrices representing cell Jacobians \f$ DF \f$  
      evaluated at cubature point sets assigned to the subcells of each cell. The leading dimension
      equals the number of subcells in the MultiCell. This array is filled on demand by 
      initializeMeasure, when Jacobian values and the associated subcell measure values at a
      particular cubature point set are requested by LocalField methods. 
      Specifically, <var>jacobianMat_[cellId][dimIndex][subcellId][cubPt]</var>  stores the cell
      Jacobian matrix evaluated at a cubature point:
      
      \li <var>cellId</var>     is the cell Id relative to the MultiCell;
      \li <var>dimIndex</var>   is cell dimension - subcell dimension;
      \li <var>subcellId</var>  is the local order of the subcell (remember, subcells of a given
                                dimension are implicitely ordered in the cell template) 
      \li <var>cubPt</var>      is the order of the cubature point from the set of points
                                assigned to that subcell.
      */
    Teuchos::Array<Teuchos::Array<Teuchos::Array<Teuchos::Array<Matrix<Scalar> > > > > jacobianMat_;
    
    
    /** \brief A four-dimensional array of Matrices which has the same dimensions as 
      <var>jacobianMat_</var> and stores the inverse transposes \f$ DF^{-T}\f$ of the 
      Jacobians in that array. Initially, <var>jacobianTInv</var> is empty and is filled  
      on demand by getJacobianTInv. Specifically, 
      <var>jacobianTInv_[cellId][dimIndex][subcellId][cubPt]</var>  stores the inverse 
      transpose of the cell Jacobian matrix evaluated at a cubature point:
      
      \li <var>cellId</var>     is the cell Id relative to the MultiCell;
      \li <var>dimIndex</var>   is cell dimension - subcell dimension;
      \li <var>subcellId</var>  is the local order of the subcell (remember, subcells of a given
                                dimension are implicitely ordered in the cell template)                                                        
      \li <var>cubPt</var>      is the order of the cubature point from the set of points
                                assigned to that subcell.
      */
    Teuchos::Array<Teuchos::Array<Teuchos::Array<Teuchos::Array<Matrix<Scalar> > > > > jacobianTInv_;
    
    
    /** \brief A four-dimensional array of Scalars which has the same dimension as 
      <var>jacobianMat_</var> and stores values of the subcell measure function at the cubature 
      point set on the subcell. This array is filled on demand by initializeMeasure, when Jacobian 
      values and subcell measure values for a specific subcell are requested by LocalField methods. 
      Definition of subcell measure function depends on the generating cell dimension and the subcell 
      dimension. In 3,2 and 1 dimensions, the contents of <var>measure_[cellId]</var> are as follows:
      
      3D:
      \li [0][0][cp]      is volume function ( = Jacobian determinant) 
      \li [1][faceId][cp] is surface area function of face with <var>faceId</var>
      \li [2][edgeId][cp] is arclength function of edge with <var>edgeId</var>
      
      2D:
      \li [0][0][cp]      is area function ( = Jacobian determinant) 
      \li [1][edgeId][cp] is arclength function of edge with <var>edgeId</var>
      
      1D:
      \li [0][0][cp]      is arclength function ( = Jacobian determinant) 

      all evaluated at cubature point <var>cp</var> from a cubature point set located on the subcell.
    */
    Teuchos::Array<Teuchos::Array<Teuchos::Array<Teuchos::Array<Scalar> > > >  measure_;
    
    
    /** \brief A four-dimensional array of Scalars which has the same dimension as 
      <var>jacobianMat_</var> and stores values of the subcell measure function at the cubature 
      point set on the subcell, <strong>times the associated cubature weight</strong>. This array is 
      filled on demand by initializeMeasure, when Jacobian values and subcell measure values for a 
      specific subcell are requested by LocalField methods. The structure of <var>weightedMeasure_</var>
      is the same as that of <var>measure_[cellId]</var>.
       */    
    Teuchos::Array<Teuchos::Array<Teuchos::Array<Teuchos::Array<Scalar> > > >  weightedMeasure_;

    
    /** \brief Disable default constructor by making it private.
    */
    MultiCell();
    
  public:
      
      /** \brief Virtual destructor
      */
      virtual ~MultiCell(){ };
      
    /** \brief Creates a multicell from a list of vertices and a generating cell type.
      
        \param numCells [in]            - Number of cells in this multicell.
        \param generatingCellType [in]  - Generating cell type: can be a canonical (TET, HEX, etc.). 
                                          or a custom, i.e. user-defined type.
        \param vertices [in]            - Physical coordinates of the vertices for all cells of the
                                          multicell.
      
      \remarks Physical coordinates must be stored in an interleaved format, e.g., for a two-cell
      multicell in 3D consisting of cells \f$A\f$ and \f$B\f$, the order is:
      \f[ \{
        \bf{v}^{A,1}_1, \bf{v}^{A,1}_2, \bf{v}^{A,1}_3, 
        \bf{v}^{A,2}_1, \bf{v}^{A,2}_2, \bf{v}^{A,2}_3,\ldots, 
        \bf{v}^{A,n}_1, \bf{v}^{A,n}_2, \bf{v}^{A,n}_3,
        \bf{v}^{B,1}_1, \bf{v}^{B,2}_2, \bf{v}^{B,1}_3, 
        \bf{v}^{B,2}_1, \bf{v}^{B,2}_2, \bf{v}^{B,2}_3,\ldots, 
        \bf{v}^{B,n}_1, \bf{v}^{B,n}_2, \bf{v}^{B,n}_3
        \}\f],\n\n
      where \f$\bf{v}^{A,p}_d\f$ is the <var>d</var>th coordinate of <var>p</var>th point in cell
      <var>A</var> and \f$n\f$ is the number of vertices in the cell.
    */
    MultiCell(const int      numCells,
              const ECell    generatingCellType,
              const Scalar*  vertices);
    
    
    /** \brief Creates a MultiCell from an array of vertices and a generating cell type. 
      
      \param generatingCellType [in]  - Generating cell type: can be a canonical (TET, HEX, etc.). 
                                        or a custom, i.e. user-defined type.
      \param vertices [in]            - Physical coordinates of the vertices for all cells of the
                                        multicell.
      
      \remarks Physical coordinates must be stored in a rank-3 array with dimensions <var>(C,P,D)</var> where
      \li <var>C</var> is the number of cells
      \li <var>P</var> is the number of nodes per cell
      \li <var>D</var> is the number of coordinates, i.e., the space dimension
      Therefore, \f$ x^{c,p}_d =  vertices(c,p,d) \f$ where \f$c\f$ is the cell number, \f$p\f$ is 
      the vertex number and \f$d\f$ is the coordinate of the point. 
    */
    template<class ArrayType>
    MultiCell(const ECell       generatingCellType,
              const ArrayType & vertices);
    
    
    //--------------------------------------------------------------------------------------------//
    //                                                                                            //
    //          Public methods of the MultiCell class which are defined for all cell types        //
    //                                                                                            //
    //--------------------------------------------------------------------------------------------//
    
    
    /** \brief Sets edge signs for 2D and 3D generating cells. 
      \warning MultiCell cannot validate correctness of the sign data in the input argument because
      it is not aware of the orientations rules applied to the global edges. The user is responsible
      for setting correct edge signs that are consistent with Intrepid's local edge orientations.
      
      \param edgeSigns [in]         - Edge signs, per each cell. For a two-cell multicell  
                                      consisting of cells \f$A\f$ and \f$B\f$, the input sequence is:\n\n
                                      \f$\{s^A_1, s^A_2, ..., s^A_m, s^B_1, s^B_2, ..., s^B_m\}\f$,\n\n
                                      where \f$m\f$ is the number of edges of the generating cell, and 
                                      \f$s^X_j \in \{1,-1,0\}\f$.
    */
    void setEdgeSigns(const short* edgeSigns);
    
    
    /** \brief Sets face signs for 2D and 3D generating cells. In 3D faces are 2-dimensional subcells, 
      i.e., true faces. In 2D the "faces" are 1-dimensional subcells, i.e., edges.
      \warning MultiCell cannot validate correctness of the sign data in the input argument because
      it is not aware of the orientations rules applied to the global faces. The user is responsible
      for setting correct face signs that are consistent with Intrepid's local face orientations.
       
      \param faceSigns [in]         - face signs, per each cell. For a two-cell multicell  
                                      consisting of cells \f$A\f$ and \f$B\f$, the input sequence is:\n\n
                                      \f$\{s^A_1, s^A_2, ..., s^A_m, s^B_1, s^B_2, ..., s^B_m\}\f$,\n\n
                                      where \f$m\f$ is the number of faces of the generating cell, and 
                                      \f$s^X_j \in \{1,-1,0\}\f$.
    */
    void setFaceSigns(const short* faceSigns);
    
    
    /** \brief Sets edge tags. 
      
      \param edgeTags  [in]         - Edge tags, per each cell. For a two-cell multicell  
                                      consisting of cells \f$A\f$ and \f$B\f$, the input sequence is:\n\n
                                      \f$\{t^A_1, t^A_2, ..., t^A_m, t^B_1, t^B_2, ..., t^B_m\}\f$,\n\n
                                      where \f$m\f$ is the number of edges in the generating cell, and 
                                      \f$t^X_j \in \{0,1\}\f$.
    */
    void setEdgeTags(const short* edgeTags);
    
    
    /** \brief Sets face tags. In 3D faces are 2-dimensional subcells, i.e., true faces. In 2D the
      "faces" are 1-dimensional subcells, i.e., edges. 
      
      \param faceTags  [in]         - Face tags, per each cell. For a two-cell multicell  
                                      consisting of cells \f$A\f$ and \f$B\f$, the input sequence is:\n\n
                                      \f$\{t^A_1, t^A_2, ..., t^A_m, t^B_1, t^B_2, ..., t^B_m\}\f$,\n\n
                                      where \f$m\f$ is the number of faces in the generating cell, and 
                                      \f$s^A_j \in \{0,1\}\f$. 
    */
    void setFaceTags(const short* faceTags);
        
    
    /** \brief A static function (can be called without prior object instantiation) to set the 
        definition of a custom cell type inside the static data member <var>connMapCustom_</var>.
      
        \param customCellType [in]      - Cell type: must be a custom type (POLY0 - POLY9).
        \param customCellTemplate [in]  - An array of 3 (three) ConnMapTemplate structs.
      
      See detailed explanation below.
      
      The <var>customCellTemplate_</var> should contain valid topological information about a user 
      defined cell. For example, the <var>customCellTemplate_</var> for a prism with a triangular
      base would look as follows:
      \code
      customMapTemplate[3] =
        {   // prism with triangular base (wedge)
          { // triprism->1cell                    //      DEFINITIONS:
            3,                                    //----> topological dimension of the cell
            9,                                    //----> number of subcells that are 1-cells 
            {2,2,2,2,2,2,2,2,2},                  //----> number of nodes per subcell
            {CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,
             CELL_EDGE,CELL_EDGE,CELL_EDGE},      //----> canonical or custom types of subcells
            {{0,1}, {1,2}, {2,0}, {0,3},
            {1,4}, {2,5}, {3,4}, {4,5}, {5,3}}    //----> local node connectivities for each subcell
          },
          { // triprism->2cell                    //      MORE CONCRETELY:
            3,                                    //----> a prism is a 3D object
            5,                                    //----> a wedge prism contains five faces
            {4,4,4,3,3},                          //----> number of nodes per face
            {CELL_QUAD,CELL_QUAD,CELL_QUAD,       //----> the faces are three quads and two triangles
             CELL_TRI,CELL_TRI},           
            {{0,1,4,3}, {1,2,5,4}, {2,0,3,5},
            {0,1,2}, {3,4,5}}                     //----> local node connectivities for each face
          },
          { // triprism->3cell                    //      MORE CONCRETELY:
            3,                                    //----> again, a prism is a 3D object
            1,                                    //----> a prism consists of one 3-cell
            {6},                                  //----> a prism has six nodes
            {CELL_TRIPRISM},                      //----> the only 3-cell is ... a prism
            {{0,1,2,3,4,5}}                       //----> local node numbers
          }
        };  // end prism
    \endcode 
    */
    static void setConnMapCustom(const ECell customCellType, 
                                 const ConnMapTemplate customCellTemplate[]);
    
    
    //--------------------------------------------------------------------------------------------//
    //                                                                                            //
    //                      Accessors operating on a specific MultiCell instance                  //
    //                                                                                            //
    //--------------------------------------------------------------------------------------------//
    
    
    /** \brief returns the number of cells in the multicell.
      */
    int getMyNumCells() const;
    
    
    /** \brief Returns the type of the generating cell of the MultiCell object.
    */
    ECell getMyCellType() const;
    
    
    /** \brief Returns the name of the generating cell of the MultiCell object.
    */
    const char* getMyCellName() const;
    
    
    /** \brief Returns the topological dimension of the generating cell of the MultiCell object.
    */
    int getMyCellDim() const;
    
    
    /** \brief Returns number of subcells of specified dimension for the generating cell of the MultiCell object.
      
      \param subcellDim [in] - dimension of the subcells whose number we want to find      
      */
    int getMyCellNumSubcells(const int subcellDim) const;
    
    
    /** \brief Returns the subcell type of a particular dimension, based on its local ID, for the 
        generating cell (remember: subcells are implicitly indexed in the cell template, starting with 0).
      
        \param subcellDim [in] - dimension of the target subcell    
        \param subcellId  [in] - local subcell ID  \n
      
      For example, if the generating cell type is CELL_PYRAMID, then
      \code
        getMySubcellType(2,0) = CELL_QUAD
      \endcode
      because the first 2-subcell in the template is the base of the cell, and
      \code
        getMySubcellType(2,1) = CELL_TRI
      \endcode
      because the second 2-subcell is a triangular side face, etc.
    */
    ECell getMySubcellType(const int subcellDim,
                           const int subcellId) const;
    
    
    /** \brief Returns local node IDs of a subcell of a particular dimension, based on its local ID, for the 
      generating cell (remember: subcells are implicitly indexed in the cell template, starting with 0).
      
      \param subcellNodeIDs [out] - ordered list of the local node IDs that are vertices of the subcell 
      \param subcellDim     [in]  - dimension of the target subcell    
      \param subcellID      [in]  - local subcell ID 

      For example, if the generating cell type is CELL_PYRAMID, then

      \code
        getMySubcellVertexIDs(subcellNodeIDs,2,0)
      \endcode
      loads <tt>{0,3,2,1}</tt> into subcellNodeIDs (the nodes of the first 2-subcell (face) of the pyramid
      in the cell template).
     
      \code 
        getMySubcellVertexIDs(subcellNodeIDs,1,2)
      \endcode
      loads <tt>{2,3}</tt> into subcellNodeIDs (the nodes of the third 1-subcell (edge) of the pyramid
      in the cell template), etc.
    */
    void getMySubcellVertexIDs(Teuchos::Array<int> & subcellNodeIDs,
                             const int subcellDim,
                             const int subcellID) const;
    
    
    /** \brief Returns reference to the data member with the edge signs for the specified cell in 
      the MultiCell (the cells stored in the MultiCell are implicitely indexed by the order 
      in which they were provided by the user).
      
      \param cellID            [in]  - cell ID relative to the MultiCell object
    */
    const Teuchos::Array<short> & getCellEdgeSigns(const int cellID) const;
    
    
    /** \brief Returns reference to the data member with the face signs for the specified cell in 
      the MultiCell (the cells stored in the MultiCell are implicitely indexed by the order 
      in which they were provided by the user).
      
      \param cellID            [in]  - cell ID relative to the MultiCell object
      */
    const Teuchos::Array<short> & getCellFaceSigns(const int cellID) const;
    
    
    /** \brief Returns reference to the data member with the edge tags for the specified cell in 
      the MultiCell (the cells stored in the MultiCell are implicitely indexed by the order 
      in which they were provided by the user).
      
      \param cellID            [in]  - cell ID relative to the MultiCell object
      */
    const Teuchos::Array<short> & getCellEdgeTags(const int cellID) const;
    
    
    /** \brief Returns reference to the data member with the face tags for the specified cell in 
      the MultiCell (the cells stored in the MultiCell are implicitely indexed by the order 
      in which they were provided by the user).
      
      \param cellID            [in]  - cell ID relative to the MultiCell object
      */
    const Teuchos::Array<short> & getCellFaceTags(const int cellID) const;
    
    
    /** \brief Returns coordinates of a single vertex in a cell as a Point object.
      
        \param cellID   [in]  - cell ID relative to the multicell
        \param vertexID [in]  - vertex ID relative to the generating cell template       
    */
    const Point<Scalar> & getCellVertex(const int cellID,
                                        const int vertexID) const;
    
    
    /** \brief Returns coordinates of all vertices in a cell as a vector of Point objects.
      
        \param cellID  [in]  - cell ID relative to the MultiCell.
    */
    const Teuchos::Array< Point<Scalar> > & getCellVertices(const int cellID) const;
    
    
    /** \brief Overloaded [] operator; returns the vertices of the cell with the specified
        cellID (relative to the MultiCell). NOTE: this allows us to use the [][] operator as well.
      
      \param cellID  [in]  - cell ID relative to the MultiCell.
    */
    const Teuchos::Array< Point<Scalar> > & operator [] (const int cellID) const;
    
    
    //--------------------------------------------------------------------------------------------//
    //                                                                                            //
    //       Static accessor functions: can be called without prior MultiCell instantiation       //
    //                                                                                            //
    //--------------------------------------------------------------------------------------------//
    
    /** \brief Returns cell type based on its name.
      
        \param cellName [in] - the name of the cell defined in <var>cellNames_</var>
    */
    static ECell getCellType(const char* cellName);
    
    
    /** \brief Returns cell name based on its type. 
      
        \param cellType [in] - target cell type 
    */
    static const char* getCellName(const ECell cellType);
    
    
    /** \brief Returns topological dimension of a cell type.
      
        \param cellType [in] - target cell type 
      */
    static int getCellDim(const ECell cellType);
    
    
    /** \brief Returns number of subcells of a particular dimension for a given cell type.
      
        \param parentCellType   [in] - parent cell type
        \param subcellDim       [in] - dimension of the subcells whose number we want to get 
    */
    static int getCellNumSubcells(const ECell parentCellType, 
                                  const int   subcellDim);
    
    
    /** \brief Returns type of subcell of a particular dimension, based on its index, for a given cell type.
      
        \param parentCellType   [in] - parent cell type
        \param subcellDim       [in] - dimension of the subcell whose type we want to get        
        \param subcellID        [in] - local ID of the desired subcell relative to the parent cell \n
                                       (remember: subcells are implicitly indexed in the cell template, starting with 0)       
    */
    static ECell getSubcellType(const ECell parentCellType,
                                const int subcellDim,
                                const int subcellID);
    
    
    /** \brief Returns local node IDs of subcell of a particular dimension, based on its index, 
               for a given cell type.
      
        \param subcellNodeIDs   [out] - ordered list of local node IDs of the desired subcell
        \param parentCellType   [in]  - parent cell type
        \param subcellDim       [in]  - dimension of the subcells whose type we want to get        
        \param subcellID        [in]  - local ID of the desired subcell, relative to the parent cell \n
                                        (remember: they are implicitly indexed in the cell template, starting with 0)
    */
    static void getSubcellVertexIDs(Teuchos::Array<int> & subcellNodeIDs,
                                  const ECell parentCellType,
                                  const int subcellDim,
                                  const int subcellID);
    
    //--------------------------------------------------------------------------------------------//
    //                                                                                            //
    //      Public methods of the MultiCell class which require the generating cell to have       //
    //      a reference cell, i.e., to be EDGE, TRI, QUAD, TET, HEX, TRIPRISM, or PYRAMID         //
    //                                                                                            //
    //--------------------------------------------------------------------------------------------//
    
    
    /** \brief Sets a default chart into an existing atlas. 
      The atlas must have size > 0, i.e., STATUS_DEFINED, for this method to work. Admissible 
      generating cell types are EDGE, TRI, QUAD, TET, HEX, TRIPRISM, and PYRAMID, i.e., cells that 
      have reference cells.
       
      \param cellID           [in]  - ID of the cell whose chart we want to set
    */
    void setChart(const int cellID);
    
    /** \brief Sets a chart of degree 1, 2 or 3 into the MultiCell atlas. 
      The atlas must have size > 0, i.e., STATUS_DEFINED, for this method to work. If 
      <var>shapePoints</var> has size zero, i.e., there are no additional shape points provided, the 
      method sets the default chart of degree 1. To define a chart of degree 2 or 3, the user must 
      provide a <var>shapePoints</var> argument with properly filled additional shape points. The 
      number and allocation of these points to the subcells must be consistent with the generating 
      cell type and the desired chart degree. Admissible generating cell types are EDGE, TRI, QUAD, 
      TET, HEX, TRIPRISM, and PYRAMID, i.e., cells that have reference cells.
      
      \param cellID           [in]  - ID of the cell whose chart we want to set
      \param shapePoints      [in]  - a set of additional shape points to compute higher degree chart.
    */
    void setChart(const int cellID, 
                  const ShapePoints<Scalar>& shapePoints);
    
    
    /** \brief Fills an undefined atlas with default cell charts, or overwrites the cell charts
      in an existing atlas by the default charts. Admissible generating cell types are EDGE, TRI, 
      QUAD, TET, HEX, TRIPRISM, and PYRAMID, i.e., cells that have reference cells.
    */
    virtual void setAtlas();

    
    /** \brief Fills an undefined atlas with cell charts based on the provided array of shapePoints, 
      or overwrites the cell charts in an existing atlas by these charts. The user is responsible for
      populating the array with <var>ShapePoints</var> structs consistent with the desired chart types  
      for each cell. The shape point sets can correspond to any valid chart degree, i.e., the charts 
      of the individual cells are not required to be of  the same degree. If a <var>ShapePoints</var>
      set of a cell has zero size, the default chart will be set for that cell. Admissible generating
      cell types  are EDGE, TRI, QUAD, TET, HEX, TRIPRISM, and PYRAMID, i.e., cells that have 
      reference cells. 
      
      \param shapePoints      [in]  - array of shape point structs for each cell in the MultiCell.
    */
    void setAtlas(const Teuchos::Array< ShapePoints<Scalar> >& shapePoints);
    
    
    /** \brief Returns the status of the <var>atlas_</var> data member (DEFINED or UNDEFINED)
    */
    EStatus getAtlasStatus() const;
    
    
    /** \brief Returns the status of the <var>atlas_</var> data member as a string.
    */
    std::string getAtlasStatusName() const;
    
    
    /** \brief Returns the mapping type (affine, non-affine) of a cell. Will throw an exception if
      the atlas has not been defined
      
      \param cellId            [in] - cell Id relative to the multicell.
      \return
        - type of the maping for this cell
    */
    const EMapping getCellMappingType(const int cellId) const;
    
    
    /** \brief Returns a Matrix object representing the Jacobian matrix of the maping specified by
      the chart of the cell, evaluated at a point in the reference frame. Admissible generating
      cell types  are EDGE, TRI, QUAD, TET, HEX, TRIPRISM, and PYRAMID, i.e., cells that have 
      reference cells.
      
      \warning Not to be confused with the Jacobian determinant, the determinant of this matrix.
      
      \param cellId    [in]   - cell Id relative to the MultiCell
      \param refPoint  [in]   - point in the reference frame where Jacobian is evaluated
      \param threshold [in]   - "tightness" in the inclusion test carried out inside the method
      */
    Matrix<Scalar> jacobian(const int            cellId, 
                            const Point<Scalar>& refPoint,
                            const double         threshold = INTREPID_THRESHOLD) const;
    
    
    /** \brief Resizes the internal storage arrays <var>jacobianMat_</var>, <var>measure_</var> and
      <var>weightedMeasure_</var> for the values of \f$DF\f$, the subcell measure function, and the
      subcell measure function multiplied by the cubature weight. Fills these arrays with data for
      all cells in the MultiCell as follows:
      
      \li Fills <var>jacobianMat_[cellId][subCellDim][subcellId]</var> with an array of <var>Matrix</var> 
      objects representing the Jacobian of the cell with the <var>cellId</var> evaluated at a cubature 
      point set assigned to the subcell specified by its subcell dimension <var>subcellDim</var>
      and local order <var>subcellId</var>. For cells with AFFINE mappings stores <strong>only one value</strong>, 
      for cells with NONAFFINE mappings stores all values.
      
      \li Fills <var>measure_[cellId][subCellDim][subcellId]</var> with an array of scalars representing 
      the values of the measure function of the specified subcell, evaluated on a cubature point set
      assigned to that subcell. For cells with AFFINE mappings stores <strong>only one value</strong>, 
      for cells with NONAFFINE mappings stores all values
      
      \li  Fills <var>weightedMeasure_[cellId][subCellDim][subcellId]</var> with an array of scalars  
      representing the values of the subcell measure function defined above, multiplied by the cubature 
      weights. Because different cubature points can have different weights, this array <strong>always
      stores as many values as there are cub. points</strong>, regardless of whether the cell is AFFINE or
      NONAFFINE.
      
      This method must be called before using any of getJacobian, getJacobianTInv, getMeasure and
      getWeightedMeasure methods.
        
      \param subcellDim       [in]     - Dimension of subcells at whose cubature points the values are computed.
      \param subcellId        [in]     - Subcell Id.
      \param cubPoints        [in]     - Array of cubature points assigned to the specified subcell
      \param cubWeights       [in]     - Array of cubature weights assigned to the points
      
      \warning <var>cubPoints</var> and <var>cubWeights</var> must belong to the same cubature rule!
      */
    void initializeMeasures(const int                              subcellDim,
                            const int                              subcellId,
                            const Teuchos::Array<Point<Scalar> > & cubPoints,
                            const Teuchos::Array<Scalar>         & cubWeights);
    
    
    /**\brief Deletes the values of \f$DF\f$, \f$DF^{-T}\f$, and measure functions at the cubature point
      set assigned to the subcell with dimension <var>subcellDim</var> and local order <var>subcellId</var>
      for all cells in the MultiCell.
      
      This method should be called <strong>before</strong> initializeMeasures any time the cubature
      point set assigned to the specified subcell has been changed.
      
      \param subcellDim       [in]     - Dimension of subcells at whose cubature points the values are computed.
      \param subcellId        [in]     - Subcell id.      
      */
    void deleteSubcellMeasures(const int  subcellDim,
                               const int  subcellId);
    
    
    /**\brief Deletes the values of \f$DF\f$, \f$DF^{-T}\f$, and measure functions for all cells in the
      MultiCell. This method should be called before initializeMeasures, if the set of cubature rules 
      for the generating cell has been completely changed.
      */
    void deleteAllMeasures();
    
    
    /** \brief Returns reference to <var>jacobianMat_[cellId][subCellDim][subcellId]</var> - an array 
      of Intrepid::Matrix objects representing the Jacobian of the cell with the <var>cellId</var>, 
      evaluated at a set of cubature points located on a subcell with dimension <var>subcell</var> 
      and local <var>subcellId</var>.
      
      \warning initializeMeasures must be called before using this method.
      
      \param cellId           [in]     - Cell Id relative to the MultiCell.
      \param subcellDim       [in]     - Dimension of subcells at whose cubature points the values are computed.
      \param subcellId        [in]     - Subcell id.
      
      \return
        - reference to an array containing \f$ DF\f$ evaluated at cubature points 
      
      \warning If the mapping of the specified cell is AFFINE, the return array contains 
      <strong>only one value</strong>!      
      */
    const Teuchos::Array<Matrix<Scalar> > & getJacobian(const int  cellId,
                                                        const int  subcellDim,
                                                        const int  subcellId);
    
    
    /** \brief Returns reference to <var>jacobianTInv_[cellId][subCellDim][subcellId]</var>
      - an array of Intrepid::Matrix objects representing the inverse transpose of the 
      Jacobian of the cell with with the <var>cellId</var>, evaluated at a set of cubature points 
      located on a subcell with dimension <var>subcell</var> and local <var>subcellId</var>.
      
      \warning initializeMeasures must be called before using this method.
      
      \param cellId           [in]     - Cell Id relative to the MultiCell.
      \param subcellDim       [in]     - Dimension of subcells at whose cubature points the values are computed.
      \param subcellId        [in]     - Subcell id.
      
      \return
      - reference to an array containing \f$ DF^{-T}\f$ evaluated at cubature points 
      
      \warning If the mapping of the specified cell is AFFINE, the return array contains 
      <strong>only one value</strong>!
    */
    const Teuchos::Array<Matrix<Scalar> > & getJacobianTInv(const int  cellId,
                                                            const int  subcellDim,
                                                            const int  subcellId);
    
    
    /** \brief Returns reference to <var>measure_[cellId][subCellDim][subcellId]</var> - an array of 
      <var>Scalar</var> values representing the measure function of the subcell with dimension 
      <var>subcell</var> and local <var>subcellId</var>, evaluated on a cubature point set assigned 
      to that subcell. 
            
      \warning initializeMeasures must be called before using this method.
      
      \param cellId           [in]     - Input operator (primitive).
      \param subcellDim       [in]     - Dimension of subcells at whose cubature points the values are computed.
      \param subcellId        [in]     - Subcell id.
      
      \return
        - reference to an array containing subcell measure function evaluated at cubature points
      
      \warning If the mapping of the specified cell is AFFINE, the return array contains 
      <strong>only one value</strong>!
    */
    const Teuchos::Array<Scalar> & getMeasure(const int  cellId,
                                              const int  subcellDim,
                                              const int  subcellId);
    
    /** \brief Returns reference to <var>weightedMeasure_[cellId][subCellDim][subcellId]</var> - an 
      array of <var>Scalar</var> values representing the measure function of the subcell with dimension 
      <var>subcell</var> and local <var>subcellId</var>, evaluated on a cubature point set assigned 
      to that subcell, and multiplied by the cubature weight assigned to each cubature point.
      
      \warning initializeMeasures must be called before using this method.
 
      \param cellId           [in]     - Input operator (primitive).
      \param subcellDim       [in]     - Dimension of subcells at whose cubature points the values are computed.
      \param subcellId        [in]     - Subcell id.
      
      \return
      - reference to an array containing subcell measure function values at cubature points
        times the cubature weight.
      
      \warning Because different cubature points can have different weights, return array <strong>always
      stores as many values as there are cub. points</strong>, regardless of whether the cell is AFFINE or
      NONAFFINE.
    */
    const Teuchos::Array<Scalar> & getWeightedMeasure(const int  cellId,
                                                      const int  subcellDim,
                                                      const int  subcellId);
    
    
    /** \brief Returns a Point object with FRAME_PHYSICAL representing the image of a point from
      the generating reference cell in the specified physical cell (relative to the MultiCell).
      The method uses the chart of the specified cell and requires <var>atlas_</var> to have 
      STATUS_DEFINED. Admissible generating cell types  are EDGE, TRI, QUAD, TET, HEX, TRIPRISM, and 
      PYRAMID, i.e., cells that have reference cells.
      
        \param  cellID [in]    - cell ID relative to the MultiCell
        \param  refPoint [in]  - point in the reference cell, must have FRAME_REFERENCE type
        \param  threshold [in] - "tightness" of the inclusion test carried out inside the method
    */
    Point<Scalar> mapToPhysicalCell(const int cellID, 
                                    const Point<Scalar>& refPoint,
                                    const double threshold = INTREPID_THRESHOLD) const;

    
    /** \brief Returns a Point object with FRAME_REFERENCE representing the image of a point from
      the specified physical cell (relative to the MultiCell) in the generating reference cell.
      The method uses the chart of the specified cell and requires <var>atlas_</var> to have 
      STATUS_DEFINED. Admissible generating cell types  are EDGE, TRI, QUAD, TET, HEX, TRIPRISM, and 
      PYRAMID, i.e., cells that have reference cells. This method uses an internal dynamic threshold value.

      
        \param cellID [in]     - cell ID relative to the MultiCell
        \param physPoint [in]  - physical point, must have FRAME_PHYSICAL type
    */
    Point<Scalar> mapToReferenceCell(const int cellID, 
                                     const Point<Scalar>& physPoint) const;
    
    //--------------------------------------------------------------------------------------------//
    //                                                                                            //
    //                                        Inclusion tests                                     //
    //                                                                                            //
    //--------------------------------------------------------------------------------------------//
    
    /** \brief Returns true if the Point argument is in the closure of the specified reference cell.
      
      Admissible cell types  are EDGE, TRI, QUAD, TET, HEX, TRIPRISM, and PYRAMID, i.e., 
      cells that have reference cells. Dimension of the Point argument must match the cell dimension 
      and its coordinate frame must be of FRAME_REFERENCE kind. This is a static function, i.e.,
      it can be called without instantiation of a MultiCell object.
      
      \param cellType  [in]  - cell type that has reference cell to be tested against
      \param refPoint  [in]  - point in FRAME_REFERENCE coordinates that is being tested
      \param threshold [in]  - "tightness" of the inclusion test
    */
    static bool inReferenceCell(const ECell cellType, 
                                const Point<Scalar>& refPoint,
                                const double threshold = INTREPID_THRESHOLD);
    
    
    /** \brief Returns true if the Point argument is in the closure of the cell with the specified 
      cellID (relative to the MultiCell).
      
      Dimension of the Point argument must match the cell dimension and its coordinate
      frame must be of FRAME_PHYSICAL kind. This method checks inclusion in physical
      space and so it requires a properly instantiated MultiCell object that holds
      valid vertex data for at least one cell. 
      
      \warning Inclusion tests for cells without reference cells, or whose charts are 
      non-affine can be slow!
      
      \param cellID [in]     - cell ID relative to the MultiCell
      \param physPoint [in]  - point in FRAME_PHYSICAL coordinates
      \param threshold [in]  - "tightness" of the inclusion test
      */
    bool inPhysicalCell(const int cellID, 
                        const Point<Scalar>& physPoint,
                        const double threshold = INTREPID_THRESHOLD) const;
    
    //--------------------------------------------------------------------------------------------//
    //                                                                                            //
    //                                         Debugging                                          //
    //                                                                                            //
    //--------------------------------------------------------------------------------------------//    
    
    /** \brief Prints multicell info to <var>os</var> stream.
    */
    void printMyInfo(std::ostream& os) const;
    
  }; // class MultiCell
  
/** \relates Intrepid::MultiCell
    \brief   Overloaded << operator.
*/
template<class Scalar>
std::ostream& operator << (std::ostream& os, const MultiCell<Scalar>& base);

} // namespace Intrepid


// include (templated) inline functions
#include "Intrepid_MultiCell.icpp"

// include templated function definitions
#include "Intrepid_MultiCellDef.hpp"

#endif

//==============================================================================
//    D O X Y G E N        D O C U M E N T A T I O N:   P A G E S             //   
//==============================================================================
/*!
 \page chart_atlas Charts and Atlasses in Intrepid
 
 Pullback is defined for canonical cells that have a standard (reference) cell. It is
 a function R^n -> R^n, where n=ambient_dimension, that maps the standard cell to
 a cell of the same type in the physical space. Therefore, to define the chart, 
 the cell dimension must match the ambient dimension. For example, it is OK to ask
 for the chart of a TRI cell if ambient_dimension = 2, but we cannot get a chart
 for a TRI cell if the ambient_dimension = 3. In this case, the TRI cell is a subcell
 of a higher dimensional cell (e.g. a TET) and its chart can be obtained by restricting
 the chart of its parent cell. 
 
 This function computes the standard charts of canonical cells, i.e., the 
 chart is a polynomial function. Pullback coefficients are stored in a 
 Pullback struct. Definition of coefficients and the storage convention is
 as follows (\f$v_{ij}\f$ denotes the \f$j\f$-th coordinate of the \f$i\f$-th vertex, except  
              in 1D where we simply write \f$v_i\f$):

 \code
 ----------------------------------------------------------------------------
 
 EDGE:  
 F: [-1,1] -> [v0,v1]
 
 F0(x) = x*(v1-v0)/2 + (v1+v0)/2
 
 Storage:  
 
 F[0][0] -> (v1-v0)/2 
 F[0][1] -> (v1+v0)/2
 ----------------------------------------------------------------------------
 TRI:    
 F: [(0,0),(1,0),(0,1)] -> [v0,v1,v2]
 
 F0(x,y) = (-v00 + v10)*x + (-v00 + v20)*y + v00
 F1(x,y) = (-v01 + v11)*x + (-v01 + v21)*y + v01
 
 Storage:  F[i][j]; i=0,1 is the coordinate function index:
 
 F[i][0] -> (-v0i + v1i)                x coefficient
 F[i][1] -> (-v0i + v2i)                y coefficient
 F[i][2] ->   v0i                       free term
 ----------------------------------------------------------------------------
 QUAD:
 F: [(-1,-1),(1,-1),(1,1),(-1,1)] -> [v0,v1,v2,v3]
 
 F0(x,y) = xy*( v00 - v10 + v20 - v30)/4 + 
 x* (-v00 + v10 + v20 - v30)/4 + 
 y* (-v00 - v10 + v20 + v30)/4 + 
 ( v00 + v10 + v20 + v30)/4
 F1(x,y) = xy*( v01 - v11 + v21 - v31)/4 + 
 x* (-v01 + v11 + v21 - v31)/4 + 
 y* (-v01 - v11 + v21 + v31)/4 + 
 ( v01 + v11 + v21 + v31)/4
 
 Storage:  F[i][j]; i=0,1 is the coordinate function index:
 
 F[i][0] -> ( v0i - v1i + v2i - v3i)/4      xy coefficient
 F[i][1] -> (-v0i + v1i + v2i - v3i)/4      x  coefficient
 F[i][2] -> (-v0i + v1i + v2i - v3i)/4      y  coefficient
 F[i][3] -> ( v0i + v1i + v2i + v3i)/4      free term
 ----------------------------------------------------------------------------
 TET:
 F:[(0,0,0),(1,0,0),(0,1,0),(0,0,1)] -> [v0,v1,v2,v3]
 
 F0(x,y,z) = x*(-v00 + v10) + y*(-v00 + v20) + z*(-v00 + v30) + v00
 F1(x,y,z) = x*(-v01 + v11) + y*(-v01 + v21) + z*(-v01 + v31) + v01
 F2(x,y,z) = x*(-v02 + v12) + y*(-v02 + v22) + z*(-v02 + v32) + v02
 
 Storage:  F[i][j]; i=0,1,2 is the coordinate function index:
 
 F[i][0] -> (-v0i + v1i)                x coefficient
 F[i][1] -> (-v0i + v2i)                y coefficient
 F[i][2] -> (-v0i + v3i)                z coefficient
 F[i][3] -> (-v0i + v3i)                free term
 ----------------------------------------------------------------------------
 HEX: 
 F:[(-1,-1,-1),(1,-1,-1),(1,1,-1),(-1,1,-1) 
   (-1,-1, 1),(1,-1, 1),(1,1, 1),(-1,1, 1)] -> [v0,v1,v2,v3,v4,v5,v6,v7]
 
 Fi(x,y,z) = xyz*(-v0i + v1i - v2i + v3i + v4i - v5i + v6i - v7i)/8 +
 xy*( v0i - v1i + v2i - v3i + v4i - v5i + v6i - v7i)/8 +
 xz*( v0i - v1i - v2i + v3i - v4i + v5i + v6i - v7i)/8 +
 yz*( v0i + v1i - v2i - v3i - v4i - v5i + v6i + v7i)/8 +
 x*(-v0i + v1i + v2i - v3i - v4i + v5i + v6i - v7i)/8 +
 y*(-v0i - v1i + v2i + v3i - v4i - v5i + v6i + v7i)/8 +
 z*(-v0i - v1i - v2i - v3i + v4i + v5i + v6i + v7i)/8 +
 ( v0i + v1i + v2i + v3i + v4i + v5i + v6i + v7i)/8
 
 Storage: F[i][j]; i=0,1,2 is the coordinate function index:
 
 F[i][0] -> (-v0i + v1i - v2i + v3i + v4i - v5i + v6i - v7i)/8
 F[i][1] -> ( v0i - v1i + v2i - v3i + v4i - v5i + v6i - v7i)/8
 F[i][2] -> ( v0i - v1i - v2i + v3i - v4i + v5i + v6i - v7i)/8
 F[i][3] -> ( v0i + v1i - v2i - v3i - v4i - v5i + v6i + v7i)/8
 F[i][4] -> (-v0i + v1i + v2i - v3i - v4i + v5i + v6i - v7i)/8
 F[i][5] -> (-v0i - v1i + v2i + v3i - v4i - v5i + v6i + v7i)/8
 F[i][6] -> (-v0i - v1i - v2i - v3i + v4i + v5i + v6i + v7i)/8
 F[i][7] -> ( v0i + v1i + v2i + v3i + v4i + v5i + v6i + v7i)/8	
 
 ----------------------------------------------------------------------------
 PRISM: 
 ----------------------------------------------------------------------------
 PYRAMID: 
 
 ----------------------------------------------------------------------------
 ----------------------------------------------------------------------------
 \endcode
*/
