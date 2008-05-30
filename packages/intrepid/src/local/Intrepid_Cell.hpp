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

/** \file   Intrepid_Cell.hpp
\brief  Header file for the Intrepid::Cell class.
\author Created by P. Bochev and D. Ridzal.
*/

#ifndef INTREPID_CELL_HPP
#define INTREPID_CELL_HPP

#include "Intrepid_ConfigDefs.hpp"
#include "Intrepid_MultiCell.hpp"

namespace Intrepid {
  
/** \class Intrepid::Cell
    \brief The Cell object (a MultiCell with only one cell) is used to communicate
           single-cell information from a (global) mesh object to a (local)
           interpolation/reconstruction operator.
*/
template<class Scalar>
class Cell : public Intrepid::MultiCell<Scalar> {
  private:
    
    /** \brief Disable default constructor.
    */
    Cell();
    
  public:
        
    /** \brief Creates a Cell from a list of vertices and a generating cell type.
      
        \param generatingCellType [in]  - Generating cell type: can be a canonical (TET, HEX, etc.). 
                                          or a custom, i.e. user-defined type.
        \param vertices [in]            - Physical coordinates of the vertices, in an
                                          interleaved format, e.g., for a cell \f$A\f$, 
                                          the order is:\n\n
                                          \f$\{x^A_1, y^A_1, z^A_1, x^A_2, y^A_2, z^A_2, ..., x^A_n, y^A_n, z^A_n\}\f$,\n\n
                                          where \f$n\f$ is the number of nodes (points) in the cell.
    */
    Cell(const ECell    generatingCellType,
         const Scalar*  vertices);
    
    
  /** \brief Returns reference to the data member with the edge signs of the cell.
  */
  const Teuchos::Array<short> & getCellEdgeSigns() const;
  
  
  /** \brief Returns reference to the data member with the face signs of the cell.
    */
  const Teuchos::Array<short> & getCellFaceSigns() const;
  
  
  /** \brief Returns reference to the data member with the edge tags of the cell.
    */
  const Teuchos::Array<short> & getCellEdgeTags() const;
  
  
  /** \brief Returns reference to the data member with the face tags of the cell.
    */
  const Teuchos::Array<short> & getCellFaceTags() const;
  

  /** \brief Returns the coordinates of a vertex as a Point object.
      
        \param vertexID [in]  - vertex ID relative to the generating cell template       
    */
    const Point<Scalar> & getCellVertex(const int vertexID) const;
    
    
    /** \brief Returns the coordinates of all vertices in a cell as an array of Point objects.
    */
    const Teuchos::Array< Point<Scalar> > & getCellVertices() const;
    
    //-------------------------------------------------------------------------------------//
    //                         Charts and atlas                                            //
    //-------------------------------------------------------------------------------------//
    
    /** \brief Sets the default chart into the cell atlas. The atlas must have size > 0,
      i.e., STATUS_DEFINED, for this method to work. Requires the generating cell to 
      have a reference cell, i.e.,  admissible cell types are EDGE, TRI, QUAD, TET, HEX, TRIPRISM, 
      and PYRAMID.
    */
    void setChart();
    
    
    /** \brief Sets a chart of degree 1, 2 or 3 into the cell atlas. The atlas must have size > 0,
      i.e., STATUS_DEFINED, for this method to work. If the argument has size zero, i.e., there are
      no additional shape points provided, the method sets the default chart of degree 1. To define 
      a chart of degree 1, 2 or 3, the user must provide a shapePoints argument with properly filled 
      additional shape points. The number and allocation of these points to the subcells must be 
      consistent with the generating cell type and the desired chart degree. This method requires 
      the generating cell to have a reference cell, i.e., admissible cell types are EDGE, TRI, QUAD, 
      TET, HEX, TRIPRISM, and PYRAMID.
      
      \param shapePoints      [in]  - a set of additional shape points to compute higher degree chart.
    */
    void setChart(const ShapePoints<Scalar>& shapePoints);
    
    
    /** \brief Fills an undefined atlas with the default cell chart, or overwrites the cell chart
      in an existing atlas by the default chart. Requires the generating cell to have a reference
      cell, i.e.,  admissible cell types are EDGE, TRI, QUAD, TET, HEX, TRIPRISM, and PYRAMID.
    */
    virtual void setAtlas();
    
    
    /** \brief Fills an undefined atlas with a cell chart based on the provided shapePoints argument, 
      or overwrites the cell chart in an existing atlas by that chart. The user is responsible for
      populating the argument with shape points consistent with the desired chart type. For argument
      of size 0 (no shape points) sets the default chart. Requires the generating cell to have a 
      reference cell, i.e.,  admissible cell types are EDGE, TRI, QUAD, TET, HEX, TRIPRISM, and PYRAMID.
      
      \param shapePoints      [in]  - a set of additional shape points to compute higher degree chart.
    */
    void setAtlas(const ShapePoints<Scalar> & shapePoints);
    
    
    /** \brief Returns the mapping type (affine, non-affine) of a cell. Will throw an exception if
      the atlas has not been defined
      
      \param cellId            [in] - cell Id relative to the multicell.
      \return
      - type of the maping for this cell
    */
    const EMapping getCellMappingType() const;    
        
    
    /** \brief Returns Matrix object containing the Jacobian matrix of the chart mapping for
      this cell, evaluated at a point in the reference frame. 
        
      \warning Not to be confused with the Jacobian determinant, the determinant of this matrix.
       
        \param refPoint [in]   - point in the reference frame where Jacobian is evaluated
        \param threshold [in]  - "tightness" in the inclusion test carried out inside the method
    */
    Matrix<Scalar> jacobian(const Point<Scalar>& refPoint,
                            const double threshold = INTREPID_THRESHOLD) const;

    /** \brief  Returns a Point object with FRAME_PHYSICAL representing the image of a point from
      the generating reference cell in the physical cell.
      The method uses the chart of the physical cell and requires <var>atlas_</var> to have 
      STATUS_DEFINED. Admissible generating cell types  are EDGE, TRI, QUAD, TET, HEX, TRIPRISM, and 
      PYRAMID, i.e., cells that have reference cells.
      
      \param  refPoint  [in]  - point in the reference cell, must have FRAME_REFERENCE type
      \param  threshold [in] - "tightness" of the inclusion test carried out inside the method
      */
    Point<Scalar> mapToPhysicalCell(const Point<Scalar>&  refPoint,
                                    const double          threshold = INTREPID_THRESHOLD) const;

    /** \brief Returns a Point object with FRAME_REFERENCE representing the image of a point from
      the physical cell in the generating reference cell.
      The method uses the chart of the physical cell and requires <var>atlas_</var> to have 
      STATUS_DEFINED. Admissible generating cell types  are EDGE, TRI, QUAD, TET, HEX, TRIPRISM, and 
      PYRAMID, i.e., cells that have reference cells. This method uses an internal dynamic threshold value.

        \param physPoint [in]  - physical point, must have FRAME_PHYSICAL type
    */
    Point<Scalar> mapToReferenceCell(const Point<Scalar>& physPoint) const;
    
    //-------------------------------------------------------------------------------------//
    //                               Inclusion tests                                       //
    //-------------------------------------------------------------------------------------//
    
    /** \brief Returns true if the Point argument belongs to the cell.
   
      Dimension of the Point argument must match the cell dimension and its coordinate
      frame must be of FRAME_PHYSICAL kind. This method checks inclusion in physical
      space and so it requires a properly instantiated Cell object that holds
      valid vertex data. 
      
      \warning Inclusion tests for cells without reference cells, or whose charts are 
               non-affine can be slow!
       
        \param physPoint [in]  - point in FRAME_PHYSICAL coordinates
        \param threshold [in]  - "tightness" of the inclusion test
    */
    bool inPhysicalCell(const Point<Scalar>& physPoint,
                                 const double threshold = INTREPID_THRESHOLD) const;
    
}; // class Cell
  
} // namespace Intrepid


// include templated function definitions
#include "Intrepid_CellDef.hpp"

#endif
