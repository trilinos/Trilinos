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
    //Cell();
    
  public:
      
    /** \brief Use this constructor if subcell signs are not needed by reconstruction interface.
      
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
    
    
    /** \brief Use this constructor if the reconstruction interface needs edge OR face signs.
      
        \param generatingCellType [in]   - Generating cell type: can be a canonical (TET, HEX, etc.). 
                                            or a custom, i.e. user-defined type.
        \param vertices   [in]           - Physical coordinates of the vertices an interleaved format (see above).      
        \param subcellSigns [in]         - Edge or face signs. For the cell \f$A\f$, the input sequence is:\n\n
                                           \f$\{s^A_1, s^A_2, ..., s^A_m\}\f$,\n\n
                                           where \f$m\f$ is the number of edges/faces per cell, and 
                                           \f$s^A_j \in \{1,-1,0\}\f$.
        \param subcellDim [in]           - dimension of the subcell type for which signs are provided (1 or 2)
      
      \warning Constructor cannot check correctness of the signs in <var>subcellSigns</var> because 
      Cell does not have access to global mesh data. The user is responsible for setting the
      correct subcell signs.
    */
    Cell(const ECell    generatingCellType,
         const Scalar*  vertices,
         const short*   subcellSigns,
         const int      subcellDim);
    
    
    /** \brief Use this constructor if the reconstruction interface needs edge AND face signs
      
        \param generatingCellType [in]  - Generating cell type: can be a canonical (TET, HEX, etc.). 
                                          or a custom, i.e. user-defined type.
        \param vertices   [in]          - Physical coordinates of the vertices in interleaved format (see above).
        \param edgeSigns [in]           - Edge signs. For the cell \f$A\f$ the input sequence is:\n\n
                                          \f$\{s^A_1, s^A_2, ..., s^A_m\}\f$,\n\n
                                          where \f$m\f$ is the number of edges per cell, and 
                                          \f$s^A_j \in \{1,-1,0\}\f$.
      \param faceSigns [in]             - Face signs. For the cell \f$A\f$ the input sequence is:\n\n
                                          \f$\{s^A_1, s^A_2, ..., s^A_m\}\f$,\n\n
                                          where \f$m\f$ is the number of faces per cell, and 
                                          \f$s^A_j \in \{1,-1,0\}\f$.
      
      \warning Constructor cannot check correctness of the signs in <var>edgeSigns</var> and 
      <var>faceSigns</var> because Cell does not have access to global mesh data. The user is 
      responsible for setting the correct edge and face signs.
    */
    Cell(const ECell    generatingCellType,
         const Scalar*  vertices,
         const short*   edgeSigns,
         const short*   faceSigns);
    

    /** \brief Returns reference to the data member that holds the signs of 1 or 2 subcells of the cell.

        \param subcellDim        [in]  - dimension of the subcells whose orientations we want
    */
    const Teuchos::Array<short> & getMySubcellSigns(const int subcellDim) const;
    
    
    /** \brief Returns the coordinates of a vertex as a Point object.
      
        \param vertexID [in]  - vertex ID relative to the generating cell template       
    */
    const Point<Scalar> & getVertex(const int vertexID) const;
    
    
    /** \brief Returns the coordinates of all vertices in a cell as an array of Point objects.
    */
    const Teuchos::Array< Point<Scalar> > & getCell() const;
    
    //-------------------------------------------------------------------------------------//
    //                         Charts and atlas                                            //
    //-------------------------------------------------------------------------------------//
    
    /** \brief Returns reference to <var>atlas_</var> that contains the chart of this cell.  
    */
    const ChartTemplate<Scalar> & getChart() const;
    
    
    /** \brief Returns Matrix object containing the Jacobian matrix of the chart mapping for
               this cell, evaluated at a point in the reference frame. 
        
      \warning Not to be confused with the Jacobian determinant, the determinant of this matrix.
       
        \param refPoint [in]   - point in the reference frame where Jacobian is evaluated
        \param threshold [in]  - "tightness" in the inclusion test carried out inside the method
    */
    Matrix<Scalar> jacobian(const Point<Scalar>& refPoint,
                            const double threshold = INTREPID_THRESHOLD) const;

    /** \brief  Maps a point from a reference cell to a target physical cell.
 
       The status of <var>atlas_</var> must be STATUS_DEFINED for this method
       to work because it uses the reference-to-physical mapping defined in the cell's chart. 
       Result is returned as a Point object with FRAME_PHYSICAL.
      
        \param  refPoint [in]  - point in the reference cell, must have FRAME_REFERENCE type
        \param  threshold [in] - "tightness" of the inclusion test carried out inside the method
    */
    Point<Scalar> mapToPhysicalCell(const Point<Scalar>& refPoint,
                                    const double threshold = INTREPID_THRESHOLD) const;

    /** \brief Maps a point from a physical cell to its reference cell.

      The status of <var>atlas_</var> must be STATUS_DEFINED for this method to work
      because it uses Newton's method to invert the reference-to-physical mapping defined
      in the cell's chart.  Result is returned as a Point object with FRAME_PHYSICAL.

        \param physPoint [in]  - physical point, must have FRAME_PHYSICAL type
    */
    Point<Scalar> mapToReferenceCell(const Point<Scalar>& physPoint) const;
    
    //-------------------------------------------------------------------------------------//
    //                               Inclusion tests                                       //
    //-------------------------------------------------------------------------------------//
    
    /** \brief Checks if the Point argument belongs to the cell.
   
      Dimension of the Point argument must match the cell dimension and its coordinate
      frame must be of FRAME_PHYSICAL kind. This method checks inclusion in physical
      space and so it requires a properly instantiated Cell object that holds
      valid vertex data. 
      
      \warning Inclusion tests for cells without reference cells, or whose charts are 
               non-affine can be slow!
       
        \param physPoint [in]  - point in FRAME_PHYSICAL coordinates
        \param threshold [in]  - "tightness" of the inclusion test
    */
    EFailCode insidePhysicalCell(const Point<Scalar>& physPoint,
                                 const double threshold = INTREPID_THRESHOLD) const;
    
}; // class Cell
  
} // namespace Intrepid


// include templated function definitions
#include "Intrepid_CellDef.hpp"

#endif
