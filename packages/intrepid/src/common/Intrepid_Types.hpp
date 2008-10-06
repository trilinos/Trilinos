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

/** \file
  \brief  Contains definitions of custom data types in Intrepid.
  \author Created by P. Bochev and D. Ridzal.
*/

#ifndef INTREPID_INTREPID_TYPES_HPP
#define INTREPID_INTREPID_TYPES_HPP

#include <Teuchos_ScalarTraits.hpp>


/** \def   INTREPID_MAX_ADJ_CELLS
  \brief The maximum number of 1,2, or 3-cells adjacent to any base 0,1,2, or 3-cell.
*/
#define INTREPID_MAX_ADJ_CELLS  24

/** \def   INTREPID_MAX_CELL_NODES
\brief The maximum number of (topological) nodes in a cell.
*/
#define INTREPID_MAX_CELL_NODES 32


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//  replaces INTREPID_MAX_CELL_NODES                         
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
/** \def   INTREPID_MAX_CELL_VERTICES
\brief The maximum number of vertices in a cell.
*/
#define INTREPID_MAX_CELL_VERTICES 32

/** \def   INTREPID_MAX_ORDER
  \brief The maximum reconstruction order.
*/
#define INTREPID_MAX_ORDER 10

/** \def   INTREPID_MAX_INTEGRATION_POINTS
  \brief The maximum number of integration points for direct cubature rules.
*/
#define INTREPID_MAX_INTEGRATION_POINTS 1001

/** \def   INTREPID_MAX_CUBATURE_DEGREE_EDGE
  \brief The maximum degree of the polynomial that can be integrated exactly by
         a direct edge rule.
*/
#define INTREPID_MAX_CUBATURE_DEGREE_EDGE 61

/** \def   INTREPID_MAX_CUBATURE_DEGREE_TRI
  \brief The maximum degree of the polynomial that can be integrated exactly by
         a direct triangle rule.
*/
#define INTREPID_MAX_CUBATURE_DEGREE_TRI 20

/** \def   INTREPID_MAX_CUBATURE_DEGREE_TET
  \brief The maximum degree of the polynomial that can be integrated exactly by
         a direct tetrahedron rule.
*/
#define INTREPID_MAX_CUBATURE_DEGREE_TET 20

/** \def   INTREPID_MAX_CUBATURE_DEGREE_PYRAMID
  \brief The maximum degree of the polynomial that can be integrated exactly by
         a direct pyramid rule.
*/
#define INTREPID_MAX_CUBATURE_DEGREE_PYRAMID 4

/** \def   INTREPID_MAX_DIMENSION
  \brief The maximum ambient space dimension.
*/
#define INTREPID_MAX_DIMENSION 3

/** \def   INTREPID_MAX_MAPPING_COEFF
  \brief  The maximum allowed number of coefficients in polynomial mappings that take reference cells
          to their ambient space images. The current setting allows to accomodate standard mappings for 
          EDGE, TRI, TET, PRISM, PYRAMID, QUAD, HEX cells. 
*/
#define INTREPID_MAX_MAPPING_COEFF 8

/** \def INTREPID_MAX_NEWTON 
  \brief Maximum number of Newton iterations to use when inverting non-affine maps between reference
         and physical elements
*/
#define INTREPID_MAX_NEWTON 15

/** \def INTREPID_MAX_DERIVATIVE
  \brief Maximum order of derivatives allowed in intrepid
*/
#define INTREPID_MAX_DERIVATIVE 10

namespace Intrepid {
  
  /* 
  Define global platform-dependent constants for various reference cell inclusion tests
  */
  static const double INTREPID_EPSILON   = std::abs(Teuchos::ScalarTraits<double>::eps());
    
  // Used in tests for inclusion of a Point in a reference cell
  static const double INTREPID_THRESHOLD = 10.0 * INTREPID_EPSILON;
  
  // Used as tolerance in e.g. Newton's method to invert non-affine mappings
  static const double INTREPID_TOL       = 10.0* INTREPID_THRESHOLD;
  
  // Used as tolerance in testing FIAT elements
  static const double INTREPID_FIAT_TOL  = 1000.0 * INTREPID_TOL;

  /** \enum  Intrepid::EStatus
      \brief Indicates the status of an object.
  */
  enum EStatus{
    STATUS_UNDEFINED=0,
    STATUS_DEFINED,
    STATUS_MAX
  };

  std::string EStatusToString(EStatus status) {
    std::string retString;
    switch(status) {
      case STATUS_UNDEFINED: retString = "Undefined";   break;
      case STATUS_DEFINED:   retString = "Defined";     break;
      case STATUS_MAX:       retString = "Max. Status"; break;
      default:               retString = "INVALID EStatus";
    }
    return retString;
  }
  

  
  /** \enum  Intrepid::EFrame
      \brief Enumeration of coordinate frames (reference/ambient) for geometrical entities (cells, points).
  */
  enum EFrame{
    FRAME_PHYSICAL=0,
    FRAME_REFERENCE
  };

  std::string EFrameToString(EFrame frame) {
    std::string retString;
    switch(frame) {
      case FRAME_PHYSICAL:  retString = "Physical";  break;
      case FRAME_REFERENCE: retString = "Reference"; break;
      default:              retString = "INVALID EFrame";
    }
    return retString;
  }
  

  
  /** \enum  Intrepid::ECoordinates
      \brief Enumeration of coordinate systems for geometrical entities (cells, points).
  */
  enum ECoordinates{
    COORDINATES_CARTESIAN=0,
    COORDINATES_POLAR,
    COORDINATES_CYLINDRICAL,
    COORDINATES_SPHERICAL,
    COORDINATES_MAX
  };

  std::string ECoordinatesToString(ECoordinates coords) {
    std::string retString;
    switch(coords) {
      case COORDINATES_CARTESIAN:   retString = "Cartesian";            break;
      case COORDINATES_POLAR:       retString = "Polar";                break;
      case COORDINATES_CYLINDRICAL: retString = "Cylindrical";          break;
      case COORDINATES_SPHERICAL:   retString = "Spherical";            break;
      case COORDINATES_MAX:         retString = "Max. Coordinates";     break;
      default:                      retString = "INVALID ECoordinates";
    }
    return retString;
  }


  
  /** \enum  Intrepid::ENorm
      \brief Enumeration of norm types for vectors and functions
  */
  enum ENorm{
    NORM_ONE = 0,
    NORM_TWO,
    NORM_INF,
    NORM_FRO,    // Frobenius matrix norm
    NORM_MAX
  };

  std::string ENormToString(ENorm norm) {
    std::string retString;
    switch(norm) {
      case NORM_ONE:   retString = "1-Norm";         break;
      case NORM_TWO:   retString = "2-Norm";         break;
      case NORM_INF:   retString = "Infinity Norm";  break;
      case NORM_FRO:   retString = "Frobenius Norm"; break;
      case NORM_MAX:   retString = "Max. Norm";      break;
      default:         retString = "INVALID ENorm";
    }
    return retString;
  }


  
  /** \enum  Intrepid::EOperator
    \brief Enumeration of primitive operators available in Intrepid. Primitive operators act on
    reconstructed functions or basis functions. Pairs of primitive operators are used to 
    specify what kind of local weak operator should be constructed.
  */
  enum EOperator{
    OPERATOR_VALUE = 0,
    OPERATOR_GRAD,      // 1
    OPERATOR_CURL,      // 2
    OPERATOR_DIV,       // 3
    OPERATOR_D1,        // 4
    OPERATOR_D2,        // 5
    OPERATOR_D3,        // 6
    OPERATOR_D4,        // 7
    OPERATOR_D5,        // 8
    OPERATOR_D6,        // 9
    OPERATOR_D7,        // 10
    OPERATOR_D8,        // 11
    OPERATOR_D9,        // 12
    OPERATOR_D10,       // 13
    OPERATOR_MAX        // 14
  };
  
  std::string EOperatorToString(EOperator op) {
    std::string retString;
    switch(op) {
      case OPERATOR_VALUE: retString = "Value";         break;
      case OPERATOR_GRAD:  retString = "Grad";          break;
      case OPERATOR_CURL:  retString = "Curl";          break;
      case OPERATOR_DIV:   retString = "Div";           break;
      case OPERATOR_D1:    retString = "D1";            break;
      case OPERATOR_D2:    retString = "D2";            break;
      case OPERATOR_D3:    retString = "D3";            break;
      case OPERATOR_D4:    retString = "D4";            break;
      case OPERATOR_D5:    retString = "D5";            break;
      case OPERATOR_D6:    retString = "D6";            break;
      case OPERATOR_D7:    retString = "D7";            break;
      case OPERATOR_D8:    retString = "D8";            break;
      case OPERATOR_D9:    retString = "D9";            break;
      case OPERATOR_D10:   retString = "D10";           break;
      case OPERATOR_MAX:   retString = "Max. Operator"; break;
      default:             retString = "INVALID EOperator";
    }
    return retString;
  }
  
  inline EOperator & operator++(EOperator &type) {
    return type = static_cast<EOperator>(type+1);
  }
    
  inline EOperator operator++(EOperator &type, int) {
    EOperator oldval = type;
    ++type;
    return oldval;
  }
   
  inline EOperator & operator--(EOperator &type) {
    return type = static_cast<EOperator>(type-1);
  }
    
  inline EOperator operator--(EOperator &type, int) {
    EOperator oldval = type;
    --type;
    return oldval;
  }

    
  
  /** \enum    Intrepid::ECell
      \brief   Enumeration of admissible cells in Intrepid. A canonical cell is one for which Intrepid 
               provides a cell template. A fixed number of enumerations is provided for user-defined cells.
               For summary of polygon types and names see http://mathworld.wolfram.com/Polygon.html .

    \warning The order of the enumeration must be exactly the same as the order of the cell
             templates defined in MultiCell<Scalar>::connMapCanonical_, Intrepid_CellTemplates. If the
             order of two enumerations is changed, the order of the associated cell template definitions in that 
             file also must be changed!
  */
  enum ECell{
    CELL_NODE = 0,       // 0-simplex, i.e. node                        value = 0
    CELL_EDGE,           // 1-simplex, i.e. edge                        value = 1
    CELL_TRI,            // 2-simplex, i.e. triangular cell             value = 2
    CELL_QUAD,           // quadrilateral cell                          value = 3
    CELL_TET,            // 3-simplex, i.e. tetrahedral cell            value = 4
    CELL_HEX,            // hexahedral cell                             value = 5
    CELL_PYRAMID,        // pyramid cell                                value = 6
    CELL_PENTAGON,       // polygon with 5 sides                        value = 7
    CELL_HEXAGON,        // polygon with 6 sides                        value = 8
    CELL_HEPTAGON,       // polygon with 7 sides                        value = 9
    CELL_OCTAGON,        // polygon with 8 sides                        value = 10
    CELL_NONAGON,        // polygon with 9 sides                        value = 11
    CELL_DECAGON,        // polygon with 10 sides                       value = 12
    CELL_TRIPRISM,       // prismatic cell with a triangle base         value = 13
    CELL_PENTAPRISM,     // prismatic polyhedron with a pentagon base   value = 14
    CELL_HEXAPRISM,      // prismatic polyhedron with a hexagon base    value = 15
    CELL_CANONICAL_MAX,  // used as the maximum number of canonical types (current value = 16)
    CELL_POLY0,          // user defined cell 
    CELL_POLY1,          // user defined cell 
    CELL_POLY2,          // user defined cell 
    CELL_POLY3,          // user defined cell 
    CELL_POLY4,          // user defined cell 
    CELL_POLY5,          // user defined cell 
    CELL_POLY6,          // user defined cell 
    CELL_POLY7,          // user defined cell 
    CELL_POLY8,          // user defined cell 
    CELL_POLY9,          // user defined cell 
    CELL_MAX             // placeholder for looping over all types        (current value = 27)
  };

  inline ECell & operator++(ECell &type) {
    return type = static_cast<ECell>(type+1);
  }
  
  inline ECell operator++(ECell &type, int) {
    ECell oldval = type;
    ++type;
    return oldval;
  }
  
  inline ECell & operator--(ECell &type) {
    return type = static_cast<ECell>(type-1);
  }
  
  inline ECell operator--(ECell &type, int) {
    ECell oldval = type;
    --type;
    return oldval;
  }


  
  /** \struct Intrepid::ConnMapTemplate
      \brief  Relational (with respect to a 1, 2, or 3-cell) connectivity map template for an arbitrary 
              MultiCell.  Three of ConnMapTemplate objects are involved in fully describing the topological 
              connectivity information of a cell. For example, the topological definition of a prism with a 
              triangular base looks as the following 3-array of ConnMapTemplate objects:
 
    \verbatim
    conn_map_template[3] =
  { // prism with triangular base (wedge)
    { // prism->1cell                        DEFINITIONS:
       3,                                  ----> topological dimension of the cell
       9,                                  ----> number of subcells that are 1-cells
      {2,2,2,2,2,2,2,2,2},                 ----> number of nodes per subcell
      {CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,
       CELL_EDGE,CELL_EDGE,CELL_EDGE,CELL_EDGE,
       CELL_EDGE},                         ----> canonical or custom types of subcells
      {{0,1}, {1,2}, {2,0}, {0,3},
       {1,4}, {2,5}, {3,4}, {4,5}, {5,3}}  ----> local node connectivities for each subcell
    },
    { // prism->2cell                        MORE CONCRETELY:
      3,                                   ----> a prism is a 3D object
      5,                                   ----> a wedge prism contains five faces
     {4,4,4,3,3},                          ----> number of nodes per face
     {CELL_QUAD,CELL_QUAD,CELL_QUAD,
      CELL_TRI,CELL_TRI},                  ----> the faces are three quads and two triangles
     {{0,1,4,3}, {1,2,5,4}, {2,0,3,5},
      {0,1,2}, {3,4,5}}                    ----> local node connectivities for each face
    },
    { // prism->3cell                        MORE CONCRETELY:
      3,                                   ----> again, a prism is a 3D object
      1,                                   ----> a prism consists of one 3-cell
      {6},                                 ----> a prism has six nodes
      {CELL_TRIPRISM},                        ----> the only 3-cell is ... a prism
      {{0,1,2,3,4,5}}                      ----> local node numbers
     }
  };  // end prism
    \endverbatim

    Also see data member documentation.
  */
  struct ConnMapTemplate {
    /** \brief Topological dimension of this cell.
    */
    int topologicalDim_;
    
    /** \brief Number of subcells of a given dimension.
    */
    int numSubcells_;
    
    /** \brief Number of nodes for each subcell of a given dimension in this cell.
      */
    int numNodesPerSubcell_[INTREPID_MAX_ADJ_CELLS];
    
    /** \brief Type of each subcell.
      */
    ECell subcellType_[INTREPID_MAX_ADJ_CELLS];
    
    /** \brief Node connectivity of each subcell.
      */
    int nodeList_[INTREPID_MAX_ADJ_CELLS][INTREPID_MAX_CELL_NODES];
  };

  

  /** \enum  Intrepid::EMapping
      \brief Enumeration of the admissible mappings in Intrepid.
  */
  enum EMapping {
    MAPPING_AFFINE = 0,
    MAPPING_NON_AFFINE,
    MAPPING_MAX
  };

  std::string EMappingToString(EMapping mapping) {
    std::string retString;
    switch(mapping) {
      case MAPPING_AFFINE:     retString = "Affine";       break;
      case MAPPING_NON_AFFINE: retString = "Non-Affine";   break;
      case MAPPING_MAX:        retString = "Max. Mapping"; break;
      default:                 retString = "INVALID EMapping";
    }
    return retString;
  }


  
  /** \struct Intrepid::ChartTemplate
      \brief  A struct to store information about the reference cell shape and the polynomial mapping
              that takes this cell to its ambient space image. The name chart alludes to the mathematical
              notion of a chart as a pair {K,F} consisting of an open set and a one-to-one mapping.
  */
  template <class Scalar>
    struct ChartTemplate {
      ECell          refCellType_;
      Scalar         mapping_[INTREPID_MAX_DIMENSION][INTREPID_MAX_MAPPING_COEFF];
      EMapping       mappingType_;
    };
  


  /** \struct Intrepid::CanonicalDofTemplate
      \brief  Array of shorts giving the number of DOF per k-dimensional subcell for the available Intrepid
              reconstruction operators. 
  */
  struct CanonicalDofTemplate {    // For use with uniform DoF allocations, i.e.,
    short numDofPerSubcell_[4];     // number of DoFs per k-subcell is same for every k-subcell
  };
  


  /** \struct Intrepid::LocalDofTag
      \brief  A data type that allows to associate a local dofId (assigned using Intrepid's canonical 
              local dof order) with a global dofId. For a fixed local dof, the data is:
              \li tag[0] - dimension of the subcell associated with the local dof
              \li tag[1] - the local id of the subcell (defined in file Intrepid_CellTemplates)
              \li tag[2] - the order of the dof relative to the subcell
              \li tag[3] - total number of dofs assigned to this subcell
  */
  struct LocalDofTag {
    int tag_[4];
  };


  
  /** \enum  Intrepid::CubatureType
      \brief Enumerates canonical (default) cubature rules in Intrepid.
  */
  enum ECubature
  {
    //
    //  GAUSS cubatures on the reference EDGE cell 
    //
    CUBATURE_GAUSS_0,
    CUBATURE_GAUSS_1,             
    CUBATURE_GAUSS_2,
    CUBATURE_GAUSS_3,
    CUBATURE_GAUSS_4,
    CUBATURE_GAUSS_5,
    CUBATURE_GAUSS_6,
    CUBATURE_GAUSS_7,
    CUBATURE_GAUSS_8,
    CUBATURE_GAUSS_9,
    CUBATURE_GAUSS_10,
    CUBATURE_GAUSS_11,
    CUBATURE_GAUSS_12,
    CUBATURE_GAUSS_13,
    CUBATURE_GAUSS_14,
    CUBATURE_GAUSS_15,
    CUBATURE_GAUSS_16,
    CUBATURE_GAUSS_17,
    CUBATURE_GAUSS_18,
    CUBATURE_GAUSS_19,
    CUBATURE_GAUSS_20,
    CUBATURE_GAUSS_21,             
    CUBATURE_GAUSS_22,
    CUBATURE_GAUSS_23,
    CUBATURE_GAUSS_24,
    CUBATURE_GAUSS_25,
    CUBATURE_GAUSS_26,
    CUBATURE_GAUSS_27,
    CUBATURE_GAUSS_28,
    CUBATURE_GAUSS_29,
    CUBATURE_GAUSS_30,
    CUBATURE_GAUSS_31,
    CUBATURE_GAUSS_32,
    CUBATURE_GAUSS_33,
    CUBATURE_GAUSS_34,
    CUBATURE_GAUSS_35,
    CUBATURE_GAUSS_36,
    CUBATURE_GAUSS_37,
    CUBATURE_GAUSS_38,
    CUBATURE_GAUSS_39,
    CUBATURE_GAUSS_40,
    CUBATURE_GAUSS_41,             
    CUBATURE_GAUSS_42,
    CUBATURE_GAUSS_43,
    CUBATURE_GAUSS_44,
    CUBATURE_GAUSS_45,
    CUBATURE_GAUSS_46,
    CUBATURE_GAUSS_47,
    CUBATURE_GAUSS_48,
    CUBATURE_GAUSS_49,
    CUBATURE_GAUSS_50,
    CUBATURE_GAUSS_51,
    CUBATURE_GAUSS_52,
    CUBATURE_GAUSS_53,
    CUBATURE_GAUSS_54,
    CUBATURE_GAUSS_55,
    CUBATURE_GAUSS_56,
    CUBATURE_GAUSS_57,
    CUBATURE_GAUSS_58,
    CUBATURE_GAUSS_59,
    CUBATURE_GAUSS_60,
    CUBATURE_GAUSS_61,
    //
    // 2D cubatures on the reference TRI
    //
    CUBATURE_TRI_0,
    CUBATURE_TRI_1,                  
    CUBATURE_TRI_2,                
    CUBATURE_TRI_3,                
    CUBATURE_TRI_4,                
    CUBATURE_TRI_5,                
    CUBATURE_TRI_6,                
    CUBATURE_TRI_7,
    CUBATURE_TRI_8,
    CUBATURE_TRI_9,
    CUBATURE_TRI_10,
    CUBATURE_TRI_11,
    CUBATURE_TRI_12,
    CUBATURE_TRI_13,
    CUBATURE_TRI_14,
    CUBATURE_TRI_15,
    CUBATURE_TRI_16,
    CUBATURE_TRI_17,
    CUBATURE_TRI_18,
    CUBATURE_TRI_19,
    CUBATURE_TRI_20,
    //
    // 3D cubatures on the reference TET
    //
    CUBATURE_TET_0,
    CUBATURE_TET_1,               
    CUBATURE_TET_2,               
    CUBATURE_TET_3,
    CUBATURE_TET_4,
    CUBATURE_TET_5,
    CUBATURE_TET_6,
    CUBATURE_TET_7,
    CUBATURE_TET_8,
    CUBATURE_TET_9,
    CUBATURE_TET_10,
    CUBATURE_TET_11,
    CUBATURE_TET_12,
    CUBATURE_TET_13,
    CUBATURE_TET_14,
    CUBATURE_TET_15,
    CUBATURE_TET_16,
    CUBATURE_TET_17,
    CUBATURE_TET_18,
    CUBATURE_TET_19,
    CUBATURE_TET_20,
    //
    // 3D cubatures on the reference PYRAMID
    //
    CUBATURE_PYRAMID_0,
    CUBATURE_PYRAMID_1,
    CUBATURE_PYRAMID_2,
    CUBATURE_PYRAMID_3,
    CUBATURE_PYRAMID_4,
    CUBATURE_MAX         // Maximum number of integration rules in Intrepid
  };
  
  

  /** \struct Intrepid::CubatureTemplate
      \brief  Template for the cubature rules used by Intrepid. Cubature template consists of cubature 
              points and cubature weights. Intrepid provides a collection of cubature templates for the 
              canonical cell shapes. The templates are defined in reference coordinates using a standard 
              reference cell for each canonical cell type. Cubature points are always specified by a triple
              of (X,Y,Z) coordinates even if the cell dimension is less than 3. The unused dimensions should
              be padded by zeroes.
    
              For example, a set of Gauss rules on [-1,1]  (the reference CELL_EDGE cell) looks as 
              the following array of CubatureTemplate objects:
    
    \verbatim
    cubature_rule[4] =
  {                                         // Collection of Gauss rules on [-1,1]
  {CELL_EDGE,                               ----> type of the cell where the rule is defined
    1,                                      ----> number of points in the rule
  {{0.0,0.0,0.0}},                          ----> X,Y,Z coordinates of the cubature points  
  {0.5}                                     ----> the cubature weight
  },
  {CELL_EDGE,
    2,
  {{-sqrt(1.0/3.0),0.0,0.0},
  {+sqrt(1.0/3.0),0.0,0.0}},
  {1.0,1.0}
  },
  {CELL_EDGE,
    3,
  {{-sqrt(3.0/5.0),0.0,0.0},
  {0.0,0.0,0.0},
  {+sqrt(3.0/5.0),0.0,0.0}},
  {5.0/9.0, 8.0/9.0,5.0/9.0}
  },
  {CELL_EDGE,
    4,
  {{-sqrt((3.0+4.0*sqrt(0.3))/7.0),0.0,0.0},
  {-sqrt((3.0-4.0*sqrt(0.3))/7.0),0.0,0.0},
  {+sqrt((3.0-4.0*sqrt(0.3))/7.0),0.0,0.0},
  {+sqrt((3.0+4.0*sqrt(0.3))/7.0),0.0,0.0}},
    //
  {0.5-sqrt(10.0/3.0)/12.0,
	0.5+sqrt(10.0/3.0)/12.0,
	0.5+sqrt(10.0/3.0)/12.0,
	0.5-sqrt(10.0/3.0)/12.0}
    
  };  // end Gauss
    \endverbatim

      Also see data member documentation.
  */
  struct CubatureTemplate {
    
    /** \brief Type of the cubature stored in the template (for verification).
    */
    ECubature         cubatureType_;
    
    /** \brief Type of the cell where the cubature is defined.
    */
    ECell             cellType_;
    
    /** \brief Number of cubature points stored in the template.
    */  
    int               numPoints_;
    
    /** \brief Array with the (X,Y,Z) coordinates of the cubature points.
    */
    double            points_[INTREPID_MAX_INTEGRATION_POINTS][INTREPID_MAX_DIMENSION];
    
    /** \brief Array with the associated cubature weights.
    */
    double            weights_[INTREPID_MAX_INTEGRATION_POINTS];
    
  };
  


  /** \enum  Intrepid::EField
    \brief Enumeration of the admissible field types in Intrepid.
  */
  enum EField
  {
    FIELD_FORM_0 = 0,
    FIELD_FORM_1,
    FIELD_FORM_2,
    FIELD_FORM_3,
    FIELD_VECTOR,
    FIELD_TENSOR,
    FIELD_MAX
  };

  std::string EFieldToString(EField field) {
    std::string retString;
    switch(field) {
      case FIELD_FORM_0: retString = "Form 0";     break;
      case FIELD_FORM_1: retString = "Form 1";     break;
      case FIELD_FORM_2: retString = "Form 2";     break;
      case FIELD_FORM_3: retString = "Form 3";     break;
      case FIELD_VECTOR: retString = "Vector";     break;
      case FIELD_TENSOR: retString = "Tensor";     break;
      case FIELD_MAX:    retString = "Max. Field"; break;
      default:           retString = "INVALID EField";
    }
    return retString;
  }
  
  inline EField & operator++(EField &type) {
    return type = static_cast<EField>(type+1);
  }
  
  inline EField operator++(EField &type, int) {
    EField oldval = type;
    ++type;
    return oldval;
  }
  
  inline EField & operator--(EField &type) {
    return type = static_cast<EField>(type-1);
  }
  
  inline EField operator--(EField &type, int) {
    EField oldval = type;
    --type;
    return oldval;
  }
  

  
  /** \enum  Intrepid::EReconstructionSpace
      \brief Enumeration of the available reconstruction spaces. Intrepid allows three basic kinds
             of reconstruction spaces for each cell type, although not all three kinds have to be defined
             for each cell type.
    
    \arg COMPLETE   on simplicial cells is complete polynomial space, on HEX and QUAD cells is 
                    tensor product space whose degree is the same in each coordinate direction
    \arg INCOMPLETE used for div and curl conforming elements of the 1st and 2nd kind
    \arg BROKEN     used for reconstructions that may subdivide a cell into subcells, i.e., a 
                    piecewise polynomial space on a cell.
                                                                                                
  */
  enum EReconstructionSpace
  {
    RECONSTRUCTION_SPACE_COMPLETE = 0,        // value = 0
    RECONSTRUCTION_SPACE_INCOMPLETE,          // value = 1
    RECONSTRUCTION_SPACE_BROKEN,              // value = 2
    RECONSTRUCTION_SPACE_MAX                  // value = 3
  };

  std::string EReconstructionSpaceToString(EReconstructionSpace space) {
    std::string retString;
    switch(space) {
      case RECONSTRUCTION_SPACE_COMPLETE:   retString = "Complete";        break;
      case RECONSTRUCTION_SPACE_INCOMPLETE: retString = "Incomplete";      break;
      case RECONSTRUCTION_SPACE_BROKEN:     retString = "Broken";          break;
      case RECONSTRUCTION_SPACE_MAX:        retString = "Max. Rec. Space"; break;
      default:                              retString = "INVALID EReconstructionSpace";
    }
    return retString;
  }

  
 
  /** \enum  Intrepid::EBasis
      \brief Enumeration of basis types in Intrepid.
  */
  enum EBasis
  {
    BASIS_FEM_DEFAULT = 0,                // value = 0
    BASIS_FEM_HIERARCHICAL,               // value = 1  
    BASIS_FEM_FIAT,                       // value = 2
    BASIS_FVD_DEFAULT,                    // value = 3
    BASIS_FVD_COVOLUME,                   // value = 4
    BASIS_FVD_MIMETIC,                    // value = 5
    BASIS_MAX                             // value = 6
  };

  std::string EBasisToString(EBasis basis) {
    std::string retString;
    switch(basis) {
      case BASIS_FEM_DEFAULT:      retString = "FEM Default";        break;
      case BASIS_FEM_HIERARCHICAL: retString = "FEM Hierarchical";   break;
      case BASIS_FEM_FIAT:         retString = "FEM FIAT";           break;
      case BASIS_FVD_DEFAULT:      retString = "FVD Default";        break;
      case BASIS_FVD_COVOLUME:     retString = "FVD Covolume";       break;
      case BASIS_FVD_MIMETIC:      retString = "FVD Mimetic";        break;
      case BASIS_MAX:              retString = "Max. Basis";         break;
      default:                     retString = "INVALID EBasis";
    }
    return retString;
  }

  
  
  /** \enum  Intrepid::EIntegrationDomain
      \brief Enumeration of integration domains.
  */
  enum EIntegrationDomain
  {
    INTEGRATION_DOMAIN_CELL = 0,             
    INTEGRATION_DOMAIN_FACE,                 
    INTEGRATION_DOMAIN_EDGE,
    INTEGRATION_DOMAIN_MAX                     
  };

  std::string EIntegrationDomainToString(EIntegrationDomain domain) {
    std::string retString;
    switch(domain) {
      case INTEGRATION_DOMAIN_CELL:    retString = "Cell";        break;
      case INTEGRATION_DOMAIN_FACE:    retString = "Face";        break;
      case INTEGRATION_DOMAIN_EDGE:    retString = "Edge";        break;
      case INTEGRATION_DOMAIN_MAX:     retString = "Max. Domain"; break;
      default:                         retString = "INVALID EIntegrationDomain";
    }
    return retString;
  }



  /** \todo Remove when sure not needed!
    \enum  Intrepid::EDataFormat
      \brief Data format used when, e.g., calling getOperator
             and getFunctional methods of the Intrepid::LocalField interface.
  */
  /*
  enum EDataFormat
  {
    DATA_SCALAR = 0,             
    DATA_VECTOR,
    DATA_TENSOR,
    DATA_FORMAT_MAX
  };

  std::string EDataFormatToString(EDataFormat dFormat) {
    std::string retString;
    switch(dFormat) {
      case DATA_SCALAR:          retString = "Scalar";              break;
      case DATA_VECTOR:          retString = "Vector";              break;
      case DATA_TENSOR:          retString = "Tensor";              break;
      case DATA_FORMAT_MAX:      retString = "Max. Data Format";    break;
      default:                   retString = "INVALID EDataFormat";
    }
    return retString;
  }
   */



  /** \enum  Intrepid::ECompEngine
      \brief Specifies how operators and functionals are computed internally
             (COMP_MANUAL = native C++ implementation, COMP_BLAS = BLAS implementation, etc.).
  */
  enum ECompEngine
  {
    COMP_CPP = 0,             
    COMP_BLAS,
    COMP_ENGINE_MAX
  };

  std::string ECompEngineToString(ECompEngine cEngine) {
    std::string retString;
    switch(cEngine) {
      case COMP_CPP:             retString = "Native C++";           break;
      case COMP_BLAS:            retString = "BLAS";                 break;
      case COMP_ENGINE_MAX:      retString = "Max. Comp. Engine";    break;
      default:                   retString = "INVALID ECompEngine";
    }
    return retString;
  }
  
  inline ECompEngine & operator++(ECompEngine &type) {
    return type = static_cast<ECompEngine>(type+1);
  }

  inline ECompEngine operator++(ECompEngine &type, int) {
    ECompEngine oldval = type;
    ++type;
    return oldval;
  }

  inline ECompEngine & operator--(ECompEngine &type) {
    return type = static_cast<ECompEngine>(type-1);
  }

  inline ECompEngine operator--(ECompEngine &type, int) {
    ECompEngine oldval = type;
    --type;
    return oldval;
  }


  
  /** \enum  Intrepid::EFailCode
      \brief Enumeration of failure codes in Intrepid.
  */
  enum EFailCode
  {
    FAIL_CODE_SUCCESS = 0,
    FAIL_CODE_INDEX_OUT_OF_RANGE,
    FAIL_CODE_TYPE_OUT_OF_RANGE,
    FAIL_CODE_MEMORY_ALLOCATION_FAILED,
    FAIL_CODE_ENTITY_NOT_FOUND,
    FAIL_CODE_MULTIPLE_ENTITIES_FOUND,
    FAIL_CODE_TAG_NOT_FOUND,
    FAIL_CODE_FILE_DOES_NOT_EXIST,
    FAIL_CODE_FILE_WRITE_ERROR,
    FAIL_CODE_NOT_IMPLEMENTED,
    FAIL_CODE_ALREADY_ALLOCATED,
    FAIL_CODE_NOT_IN_REF_CELL,
    FAIL_CODE_ERROR
  };
  


  /** \enum  Intrepid::ETagType
      \brief Enumeration of the admissible tag types in the Intrepid mesh manager.
  */
  enum ETagType
  {
    TAG_TYPE_BIT = 0,
    TAG_TYPE_SPARSE,
    TAG_TYPE_DENSE,
    TAG_TYPE_MESH,
    TAG_TYPE_LAST
  };


  
  /** \enum  Intrepid::EDataType
      \brief Enumeration of the admissible data types in the Intrepid mesh manager.
  */
  enum EDataType
  {
    DATA_TYPE_OPAQUE = 0,
    DATA_TYPE_INTEGER,
    DATA_TYPE_DOUBLE,
    DATA_TYPE_BIT,
    DATA_TYPE_HANDLE
  };
  


  /** \enum  Intrepid::EStorage
      \brief Enumeration of the admissible storage modes used to encode information in Intrepid.
  */
  enum EStorage
  {
    STORAGE_INTERLEAVED=0,
    STORAGE_BLOCKED
  };


  
  /** \enum  Intrepid::ESetOp
      \brief Enumeration of the set operations in the Intrepid mesh manager.
  */
  enum ESetOp
  {
    SET_OP_INTERSECT = 0,
    SET_OP_UNION
  };


  
  /** \typedef Intrepid::CellHandle
      \brief   CellHandle is a 32-bit integer handle.
  */
  typedef unsigned int CellHandle;


  
  /** \typedef Intrepid::Tag
      \brief   A Tag can be anything, so we define it as a void**.
  */
  typedef void** Tag;
  
  
  //===========================================================================================
  //===========================================================================================
  //===========================================================================================
  //===========================================================================================
  //===========================================================================================
  // new types added after July 20
  //===========================================================================================
  //===========================================================================================
  //===========================================================================================
  //===========================================================================================
  //===========================================================================================
  //===========================================================================================
  
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //  replaces ConnTemplate                       
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  /** \struct Intrepid::TopologyTemplate
    \brief  Relational (with respect to a 1, 2, or 3-cell) topology template to describe arbitrary
    cell types. Four TopologyTemplate structs are used to fully describe topological connectivity of
    any given (parent) cell type in up to three space dimensions. For example, the topological   
    definition of a prism with a triangular base looks as the following 4-array of TopologyTemplate 
    objects:
    
    \verbatim
    triPrismTemplate[4] =
    { // prism with triangular base
      { // triprism -> 0 cell                     DEFINITIONS:
        3,                                    ----> topological dimension of the parent cell
        6,                                    ----> number of subcells that are 0-cells
        {CELL_NODE, CELL_NODE, CELL_NODE,     ----> canonical or custom types of subcells
         CELL_NODE, CELL_NODE, CELL_NODE},    
        {{0},{1},{2},{3},{4},{5}}             ----> local vertex connectivities for each subcell
      },
      { // triprism -> 1cell                       
        3,                                    ----> topological dimension of the parent cell
        9,                                    ----> number of subcells that are 1-cells
        {CELL_EDGE, CELL_EDGE, CELL_EDGE,     ----> canonical or custom types of subcells
         CELL_EDGE, CELL_EDGE, CELL_EDGE, 
         CELL_EDGE, CELL_EDGE, CELL_EDGE},                         
        {{0,1}, {1,2}, {2,0}, {0,3},          ----> local vertex connectivities for each subcell
         {1,4}, {2,5}, {3,4}, {4,5}, {5,3}}   
      },
        { // triprism -> 2cell                     MORE CONCRETELY:
          3,                                  ----> a prism is a 3D object
          5,                                  ----> a wedge prism contains five faces
          {CELL_QUAD,CELL_QUAD,CELL_QUAD,     ----> the faces are three quads and two triangles
           CELL_TRI, CELL_TRI},                  
          {{0,1,4,3}, {1,2,5,4}, {2,0,3,5},   ----> local node connectivities for each face
           {0,1,2}, {3,4,5}}                    
        },
          { // triprism -> 3cell                   MORE CONCRETELY:
            3,                                ----> again, a prism is a 3D object
            1,                                ----> a prism consists of one 3-cell
            {CELL_TRIPRISM},                  ----> the only 3-cell is ... a prism
            {{0,1,2,3,4,5}}                   ----> local vertex numbers
          }
    };  // end prism
  \endverbatim
    
    Also see data member documentation.
    */
  struct TopologyTemplate {
    /** \brief Topological dimension of the parent cell.
    */
    int topologicalDim_;
    
    /** \brief Number of subcells of a given dimension in the parent cell.
    */
    int numSubcells_;
    
    /** \brief Type of each subcell.
      */
    ECell subcellType_[INTREPID_MAX_ADJ_CELLS];
    
    /** \brief Node connectivity of each subcell.
      */
    int vertexList_[INTREPID_MAX_ADJ_CELLS][INTREPID_MAX_CELL_VERTICES];
  };
  
  
  
  
  
  
  
  
  
  
  
  
  
} //namespace Intrepid
#endif
