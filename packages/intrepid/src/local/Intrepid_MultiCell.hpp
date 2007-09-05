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
    \author Created by P. Bochev, D. Ridzal, and D. Day.
*/


#ifndef INTREPID_MULTICELL_HPP
#define INTREPID_MULTICELL_HPP

#include "Intrepid_ConfigDefs.hpp"
#include "Intrepid_RealSpace.hpp"


namespace Intrepid {

/** \class Intrepid::MultiCell
    \brief A MultiCell (batch or group of cells) object is used to communicate
           cell information from a (global) mesh object to a (local)
           interpolation/reconstruction operator.

    \todo  Carefully define all terms (cell, canonical cell, subcell, dimension,
           orientation), and describe the Intrepid::MultiCell interface. Also
           discuss creation of custom cell shapes.
*/
template<class Scalar>
class MultiCell {
private:

  /** \brief Stores connectivity and other information about the canonical cell templates.
  */
  static const ConnMapTemplate conn_map_canonical[MAXCANONICAL][3];

  /** \brief Stores connectivity and other information about user-defined cell templates.
  */
  static ConnMapTemplate conn_map_custom[MAXTYPE-MAXCANONICAL-1][3];

  /** \brief Names of canonical cell templates.
  */
  static const char *cell_type_names[];

  /** \brief Number of cells in the multicell.
  */
  int num_cells;

  /** \brief Ambient dimension of the multicell.
  */
  int ambient_dimension;

  /** \brief Type of the multicell.
  */
  CellType my_cell_type;

  /** \brief Node coordinates.
  */
  std::vector< std::vector< Point<Scalar> > >  node_coords;

  /** \brief Edge orientations.
  */
  std::vector< std::vector<short> > edge_orients;

  /** \brief Face orientations.
  */
  std::vector< std::vector<short> > face_orients;

  /** \brief Array of pullback structs for each cell of the multicell. 
    Can be used to precompute and store pullback coefficients for each cell
    */
  std::vector< PullbackTemplate<Scalar> >  pullback;
  
  /** \brief Set to UNDEFINED by the ctor, indicates whether or not setPullbackMethod was applied
    to this instance
    */
  StatusType pullback_status;
  
  /** \brief UNDEFINED if MultiCell was constructed without providing orientation data
   */
  StatusType orients_status;

  /** \brief Disable default constructor.
  */
  MultiCell();

public:

  /** \brief Use this constructor if subcell orientations are immaterial or will be computed by Intrepid.

      \param num_cells_ [in]          - Number of cells in this multicell.
      \param ambient_dimension_ [in]  - Ambient dimension, i.e. dimension of a geometric point.
      \param cell_type_ [in]          - Generating cell type: can be a canonical (TET, HEX, etc.). or a custom,
                                        i.e. user-defined type.
      \param node_coords_ [in]        - Physical coordinates of nodes (points) for all cells of the
                                        multicell. Should be specified in the interleaved format, e.g.
                                        for a two-cell multicell consisting of cells \f$A\f$ and \f$B\f$, the order is:\n
                                        \f$\{x^A_1, y^A_1, z^A_1, x^A_2, y^A_2, z^A_2, ..., x^A_n, y^A_n, z^A_n,
                                             x^B_1, y^B_1, z^B_1, x^B_2, y^B_2, z^B_2, ..., x^B_n, y^B_n, z^B_n\}\f$,\n
                                        where \f$n\f$ is the number of nodes (points) in the cell.
  */
  MultiCell(const int      num_cells_,
            const int      ambient_dimension_,
            const CellType cell_type_,
            const Scalar*  node_coords_);

  /** \brief Use this constructor if the interpolation/reconstruction interface needs
             subcell orientations.

      \param num_cells_ [in]          - Number of cells in this multicell.
      \param ambient_dimension_ [in]  - Ambient dimension, i.e. dimension of a geometric point.
      \param cell_type_ [in]          - Generating cell type: can be a canonical (TET, HEX, etc.). or a custom,
                                        i.e. user-defined type.
      \param node_coords_ [in]        - Physical coordinates of nodes (points) for all cells of the
                                        multicell. Should be specified in the interleaved format, e.g.
                                        for a two-cell multicell consisting of cells \f$A\f$ and \f$B\f$, the input
                                        sequence is:\n
                                        \f$\{x^A_1, y^A_1, z^A_1, x^A_2, y^A_2, z^A_2, ..., x^A_n, y^A_n, z^A_n,
                                             x^B_1, y^B_1, z^B_1, x^B_2, y^B_2, z^B_2, ..., x^B_n, y^B_n, z^B_n\}\f$,\n
                                        where \f$n\f$ is the number of nodes (points) in the cell.
      \param edge_orients_ [in]       - Edge orientations, per each cell. For a two-cell multicell consisting of cells
                                        \f$A\f$ and \f$B\f$, the input sequence is:\n
                                        \f$\{o^A_1, o^A_2, ..., o^A_m, o^B_1, o^B_2, ..., o^B_m\}\f$,\n
                                         where \f$m\f$ is the number of edges per cell, and \f$o^X_j \in \{1,-1,0\}\f$.
      \param face_orients_ [in]       - Face orientations, per each cell. For a two-cell multicell consisting of cells
                                        \f$A\f$ and \f$B\f$, the input sequence is:\n
                                        \f$\{o^A_1, o^A_2, ..., o^A_m, o^B_1, o^B_2, ..., o^B_m\}\f$,\n
                                         where \f$m\f$ is the number of faces per cell, and \f$o^X_j \in \{1,-1,0\}\f$.
										 
      \warning Constructor does not check consistency of the orientations in <var>edge_orient</var>
               and <var>face_orient</var>. It is responsibility of the user to ensure that the orientation values
               corresponding to the same edge or face are always the same.
  */
  MultiCell(const int      num_cells_,
            const int      ambient_dimension_,
            const CellType cell_type_,
            const Scalar*  node_coords_,
            const short*   edge_orients_,
            const short*   face_orients_);

  /** \brief A static function (can be called without prior object instantiation)
             that sets the definition of a custom cell type inside the static data
             member <var>conn_map_custom</var>.

      \param cell_type_ [in]         - Cell type: must be a custom type.
      \param conn_map_template_ [in] - An array of 3 (three) ConnMapTemplate objects.
                                       See detailed explanation below.

             The <var>conn_map_template_</var> variable contains the topological information about
             the custom cell template. For example, the <var>conn_map_template_</var> for a prism with
             a triangular base would look as follows:
  \verbatim
  conn_map_template[3] =
  {   // prism with triangular base (wedge)
    { // prism->1cell                        DEFINITIONS:
      3,                                     ----> topological dimension of the cell
      9,                                     ----> number of subcells that are 1-cells 
      {2,2,2,2,2,2,2,2,2},                   ----> number of nodes per subcell
      {EDGE,EDGE,EDGE,EDGE,
       EDGE,EDGE,EDGE,EDGE,EDGE},            ----> canonical or custom types of subcells
      {{0,1}, {1,2}, {2,0}, {0,3},
       {1,4}, {2,5}, {3,4}, {4,5}, {5,3}}    ----> local node connectivities for each subcell
    },
    { // prism->2cell                        MORE CONCRETELY:
      3,                                     ----> a prism is a 3D object
      5,                                     ----> a wedge prism contains five faces
      {4,4,4,3,3},                           ----> number of nodes per face
      {QUAD,QUAD,QUAD,
       TRI,TRI},                             ----> the faces are three quads and two triangles
      {{0,1,4,3}, {1,2,5,4}, {2,0,3,5},
       {0,1,2}, {3,4,5}}                     ----> local node connectivities for each face
    },
    { // prism->3cell                        MORE CONCRETELY:
      3,                                     ----> again, a prism is a 3D object
      1,                                     ----> a prism consists of one 3-cell
      {6},                                   ----> a prism has six nodes
      {PRISM},                               ----> the only 3-cell is ... a prism
      {{0,1,2,3,4,5}}                        ----> local node numbers
    }
  };  // end prism
  \endverbatim 
  */
  static void setConnMapCustom(const CellType        cell_type_,
                               const ConnMapTemplate conn_map_template_[]);


  //-------------------------------------------------------------------------------------//
  //       Accessors operating on the multicell instance -> using my_cell_type           //
  //-------------------------------------------------------------------------------------//

  /** \brief Returns ambient dimension of this multicell instance
  */
  int getMyAmbientDimension() const;


  /** \brief Returns the cell type of this multicell instance.
  */
  CellType getMyType() const;


  /** \brief Returns the cell name of this multicell instance.
  */
  const char* getMyName() const;


  /** \brief Returns topological dimension of the generating cell type.
  */
  int getMyTopoDimension() const;


  /** \brief Returns number of nodes of the generating cell for this instance of multicell.
  */
  int getMyNumNodes() const;


  /** \brief Returns number of subcells of a particular dimension for the generating cell of this multicell instance.

      \param target_dim_ [in] - target dimension       
  */
  int getMyNumSubcells(const int target_dim_) const;


  /** \brief Returns type of subcell of a particular dimension, based on its index, for the generating cell.

      \param target_dim_    [in] - target dimension       
      \param subcell_index_ [in] - index of the desired subcell \n
                                   (remember: they are implicitly indexed in the cell template, starting with 0)       
  */
  CellType getMySubcellType(const int target_dim_,
                            const int subcell_index_) const;


  /** \brief Returns nodes of subcell of a particular dimension, based on its index, for the generating cell.

      \param target_dim_        [in]  - target dimension       
      \param subcell_index_     [in]  - index of the desired subcell \n
                                       (remember: they are implicitly indexed in the cell template, starting with 0)
      \param subcell_node_conn_ [out] - node connectivities of the desired subcell
  */
  void getMySubcellNodes(const int target_dim_,
                         const int subcell_index_,
                         std::vector<int> & subcell_node_conn_) const;


  /** \brief Returns orientations of 1 or 2 subcells of a specified cell 
    
    \param cell_ID             [in]  - target cell ID 
    \param subcell_dim         [in]  - dimension of the subcells whose orientations we want
    \param subcell_orient      [out] - array with the orientation values for the subcells
    */
  void getMyOrientations(const int cell_ID,
                         const int subcell_dim,
                         std::vector<short> & subcell_orient) const;
  
  //-------------------------------------------------------------------------------------//
  //                Accessors operating on any cell type specified by type_              //
  //-------------------------------------------------------------------------------------//

  /** \brief Returns cell type based on its name.
  */
  static CellType getType(const char* name_);


  /** \brief Returns cell name based on its type. This is a static function, i.e.
             it can be accessed without object instantiation.
  */
  static const char* getName(const CellType type_);


  /** \brief Returns topological dimension of a cell type.
  */
  static int getTopoDimension(const CellType type_);


  /** \brief Returns number of nodes per cell for a specified cell type.
  */
  static int getNumNodes(const CellType type_);


  /** \brief Returns number of subcells of a particular dimension for a given cell type.

      \param type_       [in] - source cell type
      \param target_dim_ [in] - target dimension
  */
  static int getNumSubcells(const CellType type_, 
                            const int target_dim_);


  /** \brief Returns type of subcell of a particular dimension, based on its index, for a given cell type.

      \param type_          [in] - source cell type
      \param target_dim_    [in] - target dimension       
      \param subcell_index_ [in] - index of the desired subcell \n
                                   (remember: they are implicitly indexed in the cell template, starting with 0)       
  */
  static CellType getSubcellType(const CellType type_,
                                 const int target_dim_,
                                 const int subcell_index_);


  /** \brief Returns nodes of subcell of a particular dimension, based on its index, for a given cell type.

      \param type_              [in]  - source cell type
      \param target_dim_        [in]  - target dimension       
      \param subcell_index_     [in]  - index of the desired subcell \n
                                       (remember: they are implicitly indexed in the cell template, starting with 0)
      \param subcell_node_conn_ [out] - node connectivities of the desired subcell
  */
  static void getSubcellNodes(const CellType type_,
                              const int target_dim_,
                              const int subcell_index_,
                              std::vector<int> & subcell_node_conn_);


  /** \brief Returns the coordinates of a node as a Point object.

      \param cell_id_  [in]  - cell ID within this multicell
      \param point_id_ [in]  - local point ID within the cell       
  */
  const Point<Scalar> & getPoint(const int cell_id_,
                                 const int point_id_) const;
  
  
  /** \brief Returns the coordinates of all nodes in a cell as a vector of Point objects.

      \param cell_id_  [in]  - cell ID within this multicell
  */
  const std::vector< Point<Scalar> > & getCell(const int cell_id_) const;


  /** \brief   Overloaded [] operator; NOTE: this allows us to use
               the [][] operator as well.
  */
  const std::vector< Point<Scalar> > & operator [] (const int cell_id_) const;

  //-------------------------------------------------------------------------------------//
  //               Pullback and other stuff..                                            //
  //-------------------------------------------------------------------------------------//
  
  /** \brief Returns the coefficients of the pullback for the cell with cell_id
    
    \param cell_id_  [in]  - cell ID within this multicell
    
    Pullback is defined for canonical cells that have a standard (reference) cell. It is
    a function R^n -> R^n, where n=ambient_dimension, that maps the standard cell to
    a cell of the same type in the physical space. Therefore, to define the pullback, 
    the cell dimension must match the ambient dimension. For example, it is OK to ask
    for the pullback of a TRI cell if ambient_dimension = 2, but we cannot get a pullback
    for a TRI cell if the ambient_dimension = 3. In this case, the TRI cell is a subcell
    of a higher dimensional cell (e.g. a TET) and its pullback can be obtained by restricting
    the pullback of its parent cell. 
    
    This function computes the standard pullbacks of canonical cells, i.e., the 
    pullback is a polynomial function. Pullback coefficients are stored in a 
    Pullback struct. Definition of coefficients and the storage convention is
    as folloows: (vij denotes the j-th coordinate of the i-th vertex, except  
    in 1D where we simply write vi)
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
    */

  /** \brief Computes the standard pullback coefficients for the given cell 
    and stores them in pullback privite member. Can only be used with canonical cell
    shapes (EDGE, TRI, QUAD, TET, HEX, PRISM, PYRAMID). 
    */
  void setPullback();
  
  
  /** \brief Returns the status of the pullback data member (DEFINED or UNDEFINED)
    */
  StatusType getPullbackStatus() const;
  
  
  /** \brief Returns the status of the pullback data member (defined or 
    undefined) as a string
    */
  const char* getPullbackInfo() const;
  
  
  /** \brief Computes the  pullback coefficients for the given cell.
    Can only be used with canonical cell shapes (EDGE, TRI, QUAD, TET, HEX,
    PRISM, PYRAMID). In Intrepid, pullback is defined as the map that takes a 
    reference cell to a specific ambient space cell.
    
    \param cell_id_    [in] - The number of the cell in the MultiCell instance
    \param pullback_   [in] - struct defining the pullback for this cell 
    
    */
  void getPullback(const int                  cell_id_, 
                   PullbackTemplate<Scalar>&  pullback_) const;
  
  
  /** \brief Evaluates the Jacobian of the pullback using a given pullback 
    template. The Point argument must have PointType = REFERENCE.
    This method provides the basic functionality needed for the next method.
    Can only be used with canonical cell shapes (EDGE, TRI, QUAD, TET, HEX,
    PRISM, PYRAMID). 
    
    In Intrepid the Jacobian of a pullback is a matrix function whose rows
    are the transposed gradients of the coordinate functions of the pullback. 
    For example, in 2D the Jacobian is DF = { {F0_x, F0_y},{F1_x,F1_y}}.

    \param cell_id_    [in] - The number of the cell in the MultiCell instance
    \param argument_   [in] - Point in its reference cell with REFERENCE type
    \param pullback_   [in] - struct defining the pullback for this cell 
    
    */
  LinearMap<Scalar> Jacobian(const int                        cell_id_,
                             const Point<Scalar>&             argument_,
                             const PullbackTemplate<Scalar>&  pullback_) const;

  
  /** \brief Evaluates the Jacobian of the pullback for the given cell and point.
    The Point argument must have PointType = REFERENCE. This method will work 
     with or without precomputed pullbacks:
    
    1) if setPullback method has been used and pullback_status = DEFINED,
       the method uses pullback data from the pullback private data member.

    2) if pullback_status == UNDEFINED, the method calls getPullback to 
       compute the pullback coefficients on the fly. This is inefficient
       if done inside a loop.

    \param cell_id_    [in] - The number of the cell in the MultiCell instance
    \param argument_   [in] - Point in its reference cell with REFERENCE type

*/  
  LinearMap<Scalar> Jacobian(const int            cell_id_,
                             const Point<Scalar>& argument_) const;
  
  
  /** \brief Finds the image of a point from the reference cell in the 
    ambient cell with the specified cell_id. Since pullback maps REFERENCE
    points to AMBIENT points, the argument must have PointType = REFERENCE.
    The return point has PointType = AMBIENT.
    
    The method does check if the preimage_ is inside the reference 
    cell. This method provides the basic pullback functionality used by
    the next pullback method. 
    
    \param cell_id_    [in] - The number of the cell in the MultiCell instance
    \param preimage_   [in] - Point in its reference cell with REFERENCE type
    \param pullback_   [in] - struct defining the pullback for the cell 
    
    */  
  Point<Scalar> Pullback(const int                       cell_id_, 
                         const Point<Scalar>&            preimage_,
                         const PullbackTemplate<Scalar>& pullback_) const;
                                              
  
  /** \brief Finds the image of a point from the reference cell in the 
    ambient cell with the specified cell_id. Since pullback maps REFERENCE
    points to AMBIENT points, the argument must have PointType = REFERENCE.
    The return point has PointType = AMBIENT.
    
    The method does check if preimage_ is inside the reference 
    cell. This method works with or without precomputed pullbacks:
    
    1) if setPullback method has been used and pullback_status = DEFINED,
       the method uses pullback data from the pullback private data member

    2) if pullback_status == UNDEFINED, the method calls getPullback to 
       compute the pullback coefficients on the fly. This is inefficient
       if done inside a loop.

    \param cell_id_    [in] - The number of the cell in the MultiCell instance
    \param preimage_   [in] - Point in its reference cell with REFERENCE type

    */
  Point<Scalar> Pullback(const int             cell_id_, 
                         const Point<Scalar>&  preimage_) const;
    
  
  /** \brief Finds the preimage of a point from the ambient cell with the 
    given cell_id in its reference cell. Reverses the action of Pullback methods. 
    Since inverse of a pullback maps AMBIENT points to REFERENCE points, the 
    argument must have PointType = AMBIENT. The return point has 
    PointType = REFERENCE. If the PullBack type is NON_AFFINE this method uses 
    Newton iteration and can be expensive. It is used by the next pullback method. 
    
    \param cell_id_    [in] - The number of the ambient cell in the MultiCell instance
    \param image_      [in] - Point in the ambient cell with AMBIENT type
    \param pullback_   [in] - struct defining the pullback for this cell 
        
    */  
  Point<Scalar> InversePullback(const int                       cell_id_, 
                                const Point<Scalar>&            image_,
                                const PullbackTemplate<Scalar>& pullback_) const;
  
  
  /**  \brief Finds the preimage of a point from the ambient  cell with the
    given cell_id in its reference cell. Reverses the action of Pullback methods. 
    Since inverse of a pullback maps AMBIENT points to REFERENCE points, the 
    argument must have PointType = AMBIENT. The return point has 
    PointType = REFERENCE. This method calls the previous InversePullback method
    and works with or without precomputed pullbacks:
    
    1) if setPullback method has been used and pullback_status = DEFINED,
       the method uses pullback data from the pullback private data member

    2) if pullback_status == UNDEFINED, the method calls getPullback to 
       compute the pullback coefficients on the fly. This is inefficient
       if done inside a loop.

    \param cell_id_    [in] - The number of the cell in the MultiCell instance
    \param image_      [in] - Point in the ambient cell with AMBIENT type
    \param pullback_   [in] - struct defining the pullback for this cell 

    */  
  Point<Scalar> InversePullback(const int               cell_id_, 
                                const Point<Scalar>&    image_) const;
  

  /** \brief Prints multicell info to <var>os</var> stream.
  */
  void printMyInfo(std::ostream& os) const;

}; // class MultiCell

/** \relates Intrepid::MultiCell
    \brief   Overloaded << operator.
*/
template<class Scalar>
std::ostream& operator << (std::ostream& os,
                           const MultiCell<Scalar>& base);

} // namespace Intrepid

// include (templated) inline functions
#include "Intrepid_MultiCell.icpp"

// include templated function definitions
#include "Intrepid_MultiCellDef.hpp"

#endif
