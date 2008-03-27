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

namespace Intrepid {
  
/** \class Intrepid::MultiCell
    \brief A MultiCell (batch or group of cells) object is used to communicate cell information from a 
           (global) mesh object to a (local) interpolation/reconstruction operator.

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
    
    
    /** \brief Dimension of the ambient space. Equals the topological dimension of the generating cell.
    */
    int ambientDim_;
    
    
    /** \brief Type of the generating cell.
    */
    ECell myCellType_;
    
    
    /** \brief Two-dimensional array whose leading dimension equals the number of cells in the MultiCell
               and the number of columns equals the number of vertices of the generating cell type.

        The i-th row  stores the vertices of the i-th cell as Point objects, i.e., vertices[i][j] is the j-th   
        vertex of the i-th cell in the MultiCell.
    */
    Teuchos::Array< Teuchos::Array< Point<Scalar> > >  vertices_;
    
    
    /** \brief Two-dimensional array whose leading dimension equals the number of cells in the MultiCell
               and the number of columns equals the number of edges (1-subcells) of the generating cell type. 

      The i-th row stores the signs of the edges of the i-th cell, i.e., edgeSigns_[i][j] is the sign
      of the j-th edge of the i-th cell. Edge signs are defined as follows:
      
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
    
    
    /** \brief STATUS_UNDEFINED if MultiCell was constructed without providing edge sign data.
      */
    EStatus edgeSignsStatus_;
    
    
    /** \brief Two-dimensional array whose leading dimension equals the number of cells in the MultiCell
               and the number of columns equals the number of faces (2-subcells) of the generating cell type. 

      The i-th row stores the signs of the faces of the i-th cell, i.e., faceSigns_[i][j] is the sign of 
      the j-th face of the i-th cell. Face signs are defined as follows:
      
      \verbatim
      faceSigns_[i][j] = +1 if local unit normal coincides with the global unit normal
      faceSigns_[i][j] = -1 if local unit normal coincides with the global unit normal
      \endverbatim
      
      Local unit normals are defined using the right hand rule and the vertex order of the faces
      in the cell template of the generating cell type.  For example, the local unit normals of the
      faces in a PYRAMID are defined by the order of their vertices in the PYRAMID cell template:
      {{0,3,2,1}, {0,1,4}, {1,2,4}, {2,3,4}, {3,0,4}} and application of the right hand rule, i.e., 
      local orientation on all faces is provided by the outer unit normal. Setting face signs is
      user's responsibility because the MultiCell is not aware of the global mesh structures and 
      orientation choices in the user code.
    */
    Teuchos::Array< Teuchos::Array<short> > faceSigns_;
    
    
    /** \brief STATUS_UNDEFINED if MultiCell was constructed without providing face sign data.
    */
    EStatus faceSignsStatus_;
    
    
    /** \brief A one-dimensional array of charts, i.e., mappings between the cells in the MultiCell 
               and their standard reference cell.

      Can be initialized if and only if the generating cell type has a chart, i.e., it is one of
      CELL_EDGE, CELL_QUAD, CELL_TRI, CELL_TET, CELL_HEX, CELL_PYRAMID, or CELL_TRIPRISM.
      By default, <var>atlas_<var> is not initialized by the ctor because not all reconstruction
      methods require this information. 
    */
    Teuchos::Array< ChartTemplate<Scalar> >  atlas_;
    
    
    /** \brief Default value is STATUS_UNDEFINED. Changed to STATUS_DEFINED by Intrepid::setAtlas()
    */
    EStatus atlasStatus_;
    
    
    /** \brief Disable default constructor.
    */
    MultiCell();
    
  public:
      
    /** \brief Use this constructor if subcell signs are not needed by reconstruction interface.
      
        \param numCells [in]            - Number of cells in this multicell.
        \param generatingCellType [in]  - Generating cell type: can be a canonical (TET, HEX, etc.). 
                                          or a custom, i.e. user-defined type.
        \param vertices [in]            - Physical coordinates of the vertices for all cells of the
                                          multicell, in an interleaved format, e.g., for a 
                                          two-cell multicell consisting of cells \f$A\f$ and \f$B\f$, 
                                          the order is:\n\n
                                          \f$\{x^A_1, y^A_1, z^A_1, x^A_2, y^A_2, z^A_2, ..., x^A_n, y^A_n, z^A_n,
                                               x^B_1, y^B_1, z^B_1, x^B_2, y^B_2, z^B_2, ..., x^B_n, y^B_n, z^B_n\}\f$,\n\n
                                          where \f$n\f$ is the number of nodes (points) in the cell.
    */
    MultiCell(const int      numCells,
              const ECell    generatingCellType,
              const Scalar*  vertices);
    
    
    /** \brief Use this constructor if the reconstruction interface needs edge OR face signs.
      
        \param numCells   [in]           - Number of cells in this multicell.
        \param generatingCellType [in]   - Generating cell type: can be a canonical (TET, HEX, etc.). 
                                            or a custom, i.e. user-defined type.
        \param vertices   [in]           - Physical coordinates of the vertices for all cells of the
                                           multicell in an interleaved format (see above).      
        \param subcellSigns [in]         - Edge or face signs, per each cell. For a two-cell multicell  
                                           consisting of cells \f$A\f$ and \f$B\f$, the input sequence is:\n\n
                                           \f$\{s^A_1, s^A_2, ..., s^A_m, s^B_1, s^B_2, ..., s^B_m\}\f$,\n\n
                                           where \f$m\f$ is the number of edges/faces per cell, and 
                                           \f$s^X_j \in \{1,-1,0\}\f$.
        \param subcellDim [in]           - dimension of the subcell type for which signs are provided (1 or 2)
      
      \warning Constructor cannot check correctness of the signs in <var>subcellSigns</var> because 
      MultiCell does not have access to global mesh data. The user is responsible for setting the
      correct subcell signs.
    */
    MultiCell(const int      numCells,
              const ECell    generatingCellType,
              const Scalar*  vertices,
              const short*   subcellSigns,
              const int      subcellDim);
    
    
    /** \brief Use this constructor if the reconstruction interface needs edge AND face signs
      
        \param numCells   [in]          - Number of cells in this multicell.
        \param generatingCellType [in]  - Generating cell type: can be a canonical (TET, HEX, etc.). 
                                          or a custom, i.e. user-defined type.
        \param vertices   [in]          - Physical coordinates of the vertices for all cells of the
                                          multicell in an interleaved format (see above).
        \param edgeSigns [in]           - Edge signs, per each cell. For a two-cell multicell consisting
                                          of cells \f$A\f$ and \f$B\f$, the input sequence is:\n\n
                                          \f$\{s^A_1, s^A_2, ..., s^A_m, s^B_1, s^B_2, ..., s^B_m\}\f$,\n\n
                                          where \f$m\f$ is the number of edges per cell, and 
                                          \f$s^X_j \in \{1,-1,0\}\f$.
      \param faceSigns [in]             - Face signs, per each cell. For a two-cell multicell consisting 
                                          of cells \f$A\f$ and \f$B\f$, the input sequence is:\n\n
                                          \f$\{s^A_1, s^A_2, ..., s^A_m, s^B_1, s^B_2, ..., s^B_m\}\f$,\n\n
                                          where \f$m\f$ is the number of faces per cell, and 
                                          \f$s^X_j \in \{1,-1,0\}\f$.
      
      \warning Constructor cannot check correctness of the signs in <var>edgeSigns</var> and 
      <var>faceSigns</var> because MultiCell does not have access to global mesh data. The user is 
      responsible for setting the correct edge and face signs.
    */
    MultiCell(const int      numCells,
              const ECell    generatingCellType,
              const Scalar*  vertices,
              const short*   edgeSigns,
              const short*   faceSigns);
    
    
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
    //                      Accessors operating on a specific MultiCell instance                  //
    //--------------------------------------------------------------------------------------------//
    
    
    /** \brief Returns the ambient dimension of the MultiCell object.
    */
    int getMyAmbientDim() const;
    
    
    /** \brief Returns the type of the generating cell of the MultiCell object.
    */
    ECell getMyCellType() const;
    
    
    /** \brief Returns the name of the generating cell of the MultiCell object.
    */
    const char* getMyCellName() const;
    
    
    /** \brief Returns the topological dimension of the generating cell of the MultiCell object.
    */
    int getMyTopologicalDim() const;
    
    
    /** \brief Returns the number of nodes of the generating cell of the MultiCell object.
    */
    int getMyNumNodes() const;
    
    
    /** \brief Returns number of subcells of specified dimension for the generating cell of the MultiCell object.
        \param subcellDim [in] - dimension of the subcells whose number we want to find      
    */
    int getMyNumSubcells(const int subcellDim) const;
    
    
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
      
      \param subcellDim     [in]  - dimension of the target subcell    
      \param subcellID      [in]  - local subcell ID 
      \param subcellNodeIDs [out] - ordered list of the local node IDs that are vertices of the subcell 

      For example, if the generating cell type is CELL_PYRAMID, then

      \code
        getMySubcellNodeIDs(2,0,subcellNodeIDs)
      \endcode
      loads <tt>{0,3,2,1}</tt> into subcellNodeIDs (the nodes of the first 2-subcell (face) of the pyramid
      in the cell template).
     
      \code 
        getMySubcellNodeIDs(1,2,subcellNodeIDs)
      \endcode
      loads <tt>{2,3}</tt> into subcellNodeIDs (the nodes of the third 1-subcell (edge) of the pyramid
      in the cell template), etc.
    */
    void getMySubcellNodeIDs(const int subcellDim,
                             const int subcellID,
                             Teuchos::Array<int> & subcellNodeIDs) const;
    
    
    /** \brief Returns reference to the data member that holds the signs of 1 or 2 subcells of the 
               specified cell in the MultiCell (remember, the cells stored in the MultiCell are implicitely 
               indexed by the order in which they were provided by the user).

        \param cellID            [in]  - cell ID relative to the MultiCell object
        \param subcellDim        [in]  - dimension of the subcells whose orientations we want
    */
    const Teuchos::Array<short> & getMySubcellSigns(const int cellID,
                                                    const int subcellDim) const;
    
    
    /** \brief Returns the coordinates of a vertex as a Point object.
      
        \param cellID   [in]  - cell ID relative to the multicell
        \param vertexID [in]  - vertex ID relative to the generating cell template       
    */
    const Point<Scalar> & getVertex(const int cellID,
                                    const int vertexID) const;
    
    
    /** \brief Returns the coordinates of all vertices in a cell as a vector of Point objects.
      
        \param cellID  [in]  - cell ID relative to the MultiCell.
    */
    const Teuchos::Array< Point<Scalar> > & getCell(const int cellID) const;
    
    
    /** \brief Overloaded [] operator; returns the vertices of the cell with the specified
               cellID (relative to the MultiCell). NOTE: this allows us to use the [][] operator as well.
      
      \param cellID  [in]  - cell ID relative to the MultiCell.
    */
    const Teuchos::Array< Point<Scalar> > & operator [] (const int cellID) const;
    
    
    //-------------------------------------------------------------------------------------//
    //   Static accessor functions: can be called without prior MultiCell instantiation    //
    //-------------------------------------------------------------------------------------//
    
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
    static int getTopologicalDim(const ECell cellType);
    
    
    /** \brief Returns number of nodes per cell for a specified cell type.

        \param cellType [in] - target cell type 
    */
    static int getNumNodes(const ECell cellType);
    
    
    /** \brief Returns number of subcells of a particular dimension for a given cell type.
      
        \param parentCellType   [in] - parent cell type
        \param subcellDim       [in] - dimension of the subcells whose number we want to get 
    */
    static int getNumSubcells(const ECell parentCellType, 
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
    
    
    /** \brief Returns node IDs of subcell of a particular dimension, based on its index, 
               for a given cell type.
      
        \param parentCellType   [in]  - parent cell type
        \param subcellDim       [in]  - dimension of the subcells whose type we want to get        
        \param subcellID        [in]  - local ID of the desired subcell, relative to the parent cell \n
                                        (remember: they are implicitly indexed in the cell template, starting with 0)
        \param subcellNodeIDs   [out] - ordered list of local node IDs of the desired subcell
    */
    static void getSubcellNodeIDs(const ECell parentCellType,
                                  const int subcellDim,
                                  const int subcellID,
                                  Teuchos::Array<int> & subcellNodeIDs);
    
    //-------------------------------------------------------------------------------------//
    //                         Charts and atlas                                            //
    //-------------------------------------------------------------------------------------//
    
    /** \brief Computes charts for all cells in the MultiCell and stores them in <var>atlas_</var>.

      The chart defines the mapping that takes a reference cell to an ambient space cell.
      Admissible generating cell types for this method are EDGE, TRI, QUAD, TET, HEX, PRISM, PYRAMID.
    */
    void setAtlas();

    /** \brief Returns the status of the <var>atlas_</var> data member (DEFINED or UNDEFINED)
    */
    EStatus getAtlasStatus() const;
    
    
    /** \brief Returns the status of the <var>atlas_</var> data member as a string.
    */
    std::string getAtlasStatusName() const;
    
    
    /** \brief Returns reference to <var>atlas_</var> that contains the chart of the specified cell  
      
        \param cellID    [in]   - cell ID relative to the MultiCell
    */
    const ChartTemplate<Scalar> & getChart(const int cellID) const;
    
    
    /** \brief Returns Matrix object containing the Jacobian matrix of the chart mapping for the
               specified cell, evaluated at a point in the reference frame. 
        
      \warning Not to be confused with the Jacobian determinant, the determinant of this matrix.
       
        \param cellID [in]     - cell ID relative to the MultiCell
        \param refPoint [in]   - point in the reference frame where Jacobian is evaluated
        \param threshold [in]  - "tightness" in the inclusion test carried out inside the method
    */
    Matrix<Scalar> jacobian(const int cellID, 
                            const Point<Scalar>& refPoint,
                            const double threshold = INTREPID_THRESHOLD) const;

    /** \brief  Maps a point from a reference cell to a target physical cell (relative to
                the MultiCell).
 
       The status of <var>atlas_</var> must be STATUS_DEFINED for this method
       to work because it uses the reference-to-physical mapping defined in the cell's chart. 
       Result is returned as a Point object with FRAME_PHYSICAL.
      
        \param  cellID [in]    - cell ID relative to the MultiCell
        \param  refPoint [in]  - point in the reference cell, must have FRAME_REFERENCE type
        \param  threshold [in] - "tightness" of the inclusion test carried out inside the method
    */
    Point<Scalar> mapToPhysicalCell(const int cellID, 
                                    const Point<Scalar>& refPoint,
                                    const double threshold = INTREPID_THRESHOLD) const;

    /** \brief Maps a point from a physical cell, relative to the MultiCell, to its reference cell.

      The status of <var>atlas_</var> must be STATUS_DEFINED for this method to work
      because it uses Newton's method to invert the reference-to-physical mapping defined
      in the cell's chart.  Result is returned as a Point object with FRAME_PHYSICAL.

        \param cellID [in]     - cell ID relative to the MultiCell
        \param physPoint [in]  - physical point, must have FRAME_PHYSICAL type
    */
    Point<Scalar> mapToReferenceCell(const int cellID, const Point<Scalar>& physPoint) const;
    
    //-------------------------------------------------------------------------------------//
    //                               Inclusion tests                                       //
    //-------------------------------------------------------------------------------------//
    
    /** \brief Checks if the Point argument belongs to the reference cell of the specified type.

      Valid cell type range is CELL_EDGE to CELL_TRIPRISM, i.e., cell types that have
      reference cells. Dimension of the Point argument must match the cell dimension and
      its coordinate frame must be of FRAME_REFERENCE kind. This is a static function, i.e.,
      it can be called without instantiation of a MultiCell object.
      
        \param cellType  [in]  - cell type that has reference cell to be tested against
        \param refPoint  [in]  - point in FRAME_REFERENCE coordinates that is being tested
        \param threshold [in]  - "tightness" of the inclusion test
    */
    static EFailCode insideReferenceCell(const ECell cellType, 
                                         const Point<Scalar>& refPoint,
                                         const double threshold = INTREPID_THRESHOLD);
    
    /** \brief Checks if the Point argument belongs to the cell with the specified cellID
               (relative to the MultiCell).
   
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
    EFailCode insidePhysicalCell(const int cellID, 
                                 const Point<Scalar>& physPoint,
                                 const double threshold = INTREPID_THRESHOLD) const;
    
    //-------------------------------------------------------------------------------------//
    //                                     Debugging                                       //
    //-------------------------------------------------------------------------------------//    
    
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
