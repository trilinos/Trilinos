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
    \brief Contains definitions of custom data types in Intrepid.
    \author Created by P. Bochev, D. Ridzal, and D. Day.
*/


#ifndef INTREPID_INTREPID_TYPES_HPP
#define INTREPID_INTREPID_TYPES_HPP

/** \def   INTREPID_MAX_ADJ_CELLS
    \brief The maximum number of 1,2, or 3-cells adjacent to any base 0,1,2, or 3-cell.
*/
#define INTREPID_MAX_ADJ_CELLS  24

/** \def   INTREPID_MAX_CELL_NODES
\brief The maximum number of (topological) nodes in a cell.
*/
#define INTREPID_MAX_CELL_NODES 32

/** \def   INTREPID_MAX_ORDER
\brief The maximum reconstruction order.
*/
#define INTREPID_MAX_ORDER 10

/** \def   INTREPID_MAX_INTEGRATION_POINTS
\brief The maximum number of integration points.
*/
#define INTREPID_MAX_INTEGRATION_POINTS 1001

/** \def   INTREPID_MAX_DIMENSION
\brief The maximum ambient space dimension.
*/
#define INTREPID_MAX_DIMENSION 3


/** \def   INTREPID_MAX_PULLBACK_COEFF
\brief The maximum number of coefficients in a pullback map that takes reference
cell to its physical image. The present value is set to accomodate standard pullbacks
for EDGE, TRI, TET, PRISM, PYRAMID, QUAD, HEX
*/
#define INTREPID_MAX_PULLBACK_COEFF 8

namespace Intrepid {

/** \enum  Intrepid::StatusType
  \brief To indicate the status of an object.
  */
enum StatusType{
  UNDEFINED=0,
  DEFINED,
  MAX_STATUS_TYPE
};

static const char* StatusNames[] = {
  "UNDEFINED",
  "DEFINED",
  "MAX_STATUS_TYPE"
};

/** \enum Intrepid::PointType
\brief Indicates the space (reference or ambient) where a point is located
*/
enum PointType{
  AMBIENT=0,
  REFERENCE
};

static const char* PointNames[]={
  "  AMBIENT",
  "REFERENCE"
};

/** \enum Intrepid::OperatorType
\brief A list of differential and other operators available in Intrepid
*/

enum OperatorType{
  MASS = 0,
  GRAD,
  CURL,
  DIV,
  DIV_GRAD,
  CURL_CURL,
  GRAD_DIV,
  MAX_CANONICAL_OPERATORS
};

/* These will be disabled until we actually need them,
   in order to prevent compiler warnings.
static const char* OperatorNames[]={
  "Mass",
  "Grad",
  "Curl",
  "Div",
  "DivGrad",
  "CurlCurl",
  "GradDiv",
  "Max. Canonical Operator"
};
*/

/** \enum  Intrepid::CellType
    \brief Enumerates types of cells in Intrepid.

           In addition to the canonical cell types, there is a fixed
           number of generic polygonal and polyhedral cells. They are
           intended as placeholders for user-defined (custom) cell shapes.
*/
enum CellType
{
  NODE = 0,       // 0-simplex, i.e. node
  EDGE,           // 1-simplex, i.e. edge
  TRI,            // 2-simplex, i.e. triangular cell
  QUAD,           // quadrilateral cell
  TET,            // 3-simplex, i.e. tetrahedral cell
  HEX,            // hexahedral cell
  PYRAMID,        // pyramid cell
  PRISM,          // prism (wedge) cell
  POLYGON,        // general (unknown) polygon
  POLYHEDRON,     // general (unknown) polyhedron
  CELLSET,        // set of arbitrary cells
  MAXCANONICAL,   // used as the maximum number of canonical types
  POLYGON1,       // general polygonal cell (user-defined)
  POLYGON2,       // general polygonal cell (user-defined)
  POLYGON3,       // general polygonal cell (user-defined)
  POLYGON4,       // general polygonal cell (user-defined)
  POLYGON5,       // general polygonal cell (user-defined)
  POLYHEDRON1,    // general polyhedral cell (user-defined)
  POLYHEDRON2,    // general polyhedral cell (user-defined)
  POLYHEDRON3,    // general polyhedral cell (user-defined)
  POLYHEDRON4,    // general polyhedral cell (user-defined)
  POLYHEDRON5,    // general polyhedral cell (user-defined)
  MAXTYPE         // placeholder for looping over all types
};

inline CellType & operator++(CellType &type_) {
  return type_ = static_cast<CellType>(type_+1);
}

inline CellType operator++(CellType &type_, int) {
  CellType oldval = type_;
  ++type_;
  return oldval;
}

inline CellType & operator--(CellType &type_) {
  return type_ = static_cast<CellType>(type_-1);
}

inline CellType operator--(CellType &type_, int) {
  CellType oldval = type_;
  --type_;
  return oldval;
}

/** \struct Intrepid::ConnMapTemplate
    \brief  Relational (with respect to a 1, 2, or 3-cell) connectivity map template for
            an arbitrary MultiCell.

            Three of ConnMapTemplate objects are involved in fully describing the topological
            connectivity information of a cell.
            For example, the topological definition of a prism with a triangular base looks
            as the following 3-array of ConnMapTemplate objects:
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

            Also see data member documentation.
*/
struct ConnMapTemplate {
  /** \brief Topological dimension of this cell.
  */
  int topo_dimension;

  /** \brief Number of subcells of a given dimension.
  */
  int num_subcells;

  /** \brief Number of nodes for each subcell of a given dimension in this cell.
  */
  int num_nodes_per_subcell[INTREPID_MAX_ADJ_CELLS];

  /** \brief Type of each subcell.
  */
  CellType type_subcell[INTREPID_MAX_ADJ_CELLS];

  /** \brief Node connectivity of each subcell.
  */
  int node_conn[INTREPID_MAX_ADJ_CELLS][INTREPID_MAX_CELL_NODES];
};

/** \enum Intrepid::PullbackType
\brief Describes the type of pullback map (affine or non-affine)
*/
enum PullbackType {
  AFFINE = 0,
  NON_AFFINE
};

/* These will be disabled until we actually need them,
   in order to prevent compiler warnings.
static const char* PullbackNames[] = {
  "AFFINE",
  "NON_AFFINE"
};
*/

/** \struct Intrepid::PullbackTemplate
\brief  Stores coefficients of the pullback for the canonical cell types. 
*/
template <class Scalar>
struct PullbackTemplate {
  PullbackType  pullback_type;
  CellType      cell_type;
  int           cell_id;
  Scalar        F[INTREPID_MAX_DIMENSION][INTREPID_MAX_PULLBACK_COEFF];
};


/** \struct Intrepid::CanonicalDofTemplate
\brief  Array of shorts giving the number of DOF per k-dimensional subcell for the available Intrepid
        reconstruction operators. 
*/
struct CanonicalDofTemplate {    // For use with uniform DoF allocations, i.e.,
  short num_dof_per_subcell[4];  // number of DoFs per k-subcell is same for every k-subcell
};


/** \struct Intrepid::LocalDofTag
\brief A data type that allows to associate a local dof_id (assigned using 
 Intrepid's canonical local dof order) with a global dof_id. For a fixed local
 dof, the data is:
 tag[0] - dimension of the subcell associated with the local dof
 tag[1] - the local id of the subcell (defined in Intrepid_CellTemplates)
 tag[2] - the order of the dof relative to the subcell (if more than 1 per subcell)
*/
struct LocalDofTag {
  short tag[3];
};

/** \enum  Intrepid::CubatureType
\brief Enumerates types of cubature templates in Intrepid.

*/
enum CubatureType
{
  DEFAULT = 0,          // contains no data, just tells to use whatever is the default for the reconstruction
  //
  //  GAUSS cubatures on the reference EDGE cell 
  //
  GAUSS_0,
  GAUSS_1,             
  GAUSS_2,
  GAUSS_3,
  GAUSS_4,
  GAUSS_5,
  GAUSS_6,
  GAUSS_7,
  GAUSS_8,
  GAUSS_9,
  GAUSS_10,
  GAUSS_11,
  GAUSS_12,
  GAUSS_13,
  GAUSS_14,
  GAUSS_15,
  GAUSS_16,
  GAUSS_17,
  GAUSS_18,
  GAUSS_19,
  GAUSS_20,
  //
  // 2D cubatures on the reference TRI
  //
  TRI_0,
  TRI_1,                // 1-point  
  TRI_2,                // 3-point rule
  TRI_3,                // 7-point rule
  TRI_4,                // 6-point rule
  TRI_5,                // 7-point rule
  TRI_6,                // 12-point rule
  TRI_7,
  TRI_8,
  TRI_9,
  TRI_10,
  TRI_11,
  TRI_12,
  TRI_13,
  TRI_14,
  TRI_15,
  TRI_16,
  TRI_17,
  TRI_18,
  TRI_19,
  TRI_20,
  
  //
  // 3D cubatures on the reference TET
  //
  TET_0,
  TET_1,               // 1-point rule
  TET_2,               // 4-point rule
  TET_3,
  TET_4,
  TET_5,
  TET_6,
  TET_7,
  TET_8,
  TET_9,
  TET_10,
  TET_11,
  TET_12,
  TET_13,
  TET_14,
  TET_15,
  TET_16,
  TET_17,
  TET_18,
  TET_19,
  TET_20,
  //
  // 3D cubatures on the reference PYRAMID
  //
  PYRAMID_1,
  PYRAMID_2,
  PYRAMID_3,
  PYRAMID_4,
  //
  // 3D cubatures on the reference PRISM
  //
  PRISM_1,
  PRISM_2,
  PRISM_3,
  PRISM_4,
  //
  MAX_LOWORDER,         // Maximum number of simple integration rules on standard cell shapes (value = 21)
  //
  // High-order cubature templates on standard cell shapes (EDGE, TRI, TET, PYRAMID, PRISM)
  //
  GAUSS_P,              //  Higher order 1D Gauss rules on an EDGE
  TRI_P,                // Higher order 2D rules on a TRI
  TET_P,                // Higher order 3D rules on a TET
  PYRAMID_P,            // Higher order 3D rules on a PYRAMID
  PRISM_P,              // Higher order 3D rules on a PRYSM
  //
  MAX_HIGHORDER,        // Maximum number of higher order integration rules on standard cell shapes (value = 27)
  //
  // Cubature templates for FVD/CUSTOM reconstructions on standard cell shapes 
  //
  GAUSS_CUSTOM,         // Use for FVD and CUSTOM reconstrruction types on EDGE, QUAD and HEX
  TRI_CUSTOM,           // Use for FVD and CUSTOM reconstrruction types on TRI
  TET_CUSTOM,           // Use for FVD and CUSTOM reconstrruction types on TET
  PYRAMID_CUSTOM,       // Use for FVD and CUSTOM reconstrruction types on PYRAMID
  PRISM_CUSTOM,         // Use for FVD and CUSTOM reconstrruction types on PRISM
  //
  // Cubature templates for FVD/CUSTOM reconstruction on POLYGON and POLYHEDRON cell families
  //
  POLYGON_RULE,         // generic 2D integration rule for a POLYGON
  POLYHEDRON_RULE,      // Generic 3D integration rule for POLYHEDRON
  //
  MAX_CUSTOMRULE        // Maximum number of generic integration rules for custom cell shapes (value = 30)
};


/** \struct Intrepid::CubatureTemplate
\brief  Template for the cubature rules used by Intrepid. 

Cubature template consists of cubature points and cubature weights. Intrepid 
provides a collection of cubature templates for the canonical cell shapes. 
The templates are defined in reference coordinates using a standard reference 
cell for each canonical cell type. Cubature points in physical coordinates
are obtained by using the MultiCell::getPullback method.

For example, a set of Gauss rules on [-1,1]  (the reference EDGE cell) looks as 
the following array of CubatureTemplate objects:
\verbatim
cubature_rule[4] =
{   // Collection of Gauss rules on [-1,1]
  {EDGE,
    1,
    {{0.0,0.0,0.0}},
    {0.5}
    },
  {EDGE,
    2,
    {{-sqrt(1.0/3.0),0.0,0.0},
    {+sqrt(1.0/3.0),0.0,0.0}},
    {1.0,1.0}
    },
  {EDGE,
    3,
    {{-sqrt(3.0/5.0),0.0,0.0},
	{0.0,0.0,0.0},
	{+sqrt(3.0/5.0),0.0,0.0}},
    {5.0/9.0, 8.0/9.0,5.0/9.0}
    },
  {EDGE,
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
Rules for general cells are defined separately because they depend on the 
reconstruction operator.  Also see data member documentation.
*/
struct CubatureTemplate {
  
  /** \brief  Type of the cubature  template (for verification of correctness).
   */
  CubatureType      template_type;
  
  /** \brief The admiissible cell type for this cubature template.
  */
  CellType          admissible_cell;
  
  /** \brief Number of cubature points in the template.
  */  
  int               num_points;
  
  /** \brief Array of cubature points
    */
  double            points[INTREPID_MAX_INTEGRATION_POINTS][INTREPID_MAX_DIMENSION];
    
  /** \brief Array containing the weights for the cubature template.
    */
  double            weights[INTREPID_MAX_INTEGRATION_POINTS];
  
};
  


/** \enum  Intrepid::ReconstructionType
    \brief Enumerates the reconstruction types available for p-Forms.
*/
enum ReconstructionType
{
  FEM = 0,
  FVD,
  CUSTOM,
  MAX_RECONSTRUCTION_TYPES
};

static const char* ReconstructionTypeNames[] = {
  "FEM",
  "FVD",
  "CUSTOM",
  "MAX_RECONSTRUCTION_TYPES"
};

/** \enum  Intrepid::ReconstructionFlavor_0Form
    \brief Enumerates the flavor of the reconstruction types available for 0-Forms.
*/
enum ReconstructionFlavor_0Form
{
  LAGRANGE_0FORM = 0,
  HIERARCHICAL_0FORM,
  MAX_FEM_FLAVORS_0FORM,    // Last flavor for FEM reconstruction
  COVOLUME_0FORM,
  MIMETIC_0FORM,
  CUSTOM_0FORM,
  MAX_FVD_FLAVORS_0FORM     // Last flavor for FVD reconstruction
};


/** \enum  Intrepid::ReconstructionFlavor_1Form
    \brief Enumerates the flavor of the reconstruction types available for 1-Forms.
*/
enum ReconstructionFlavor_1Form
{
  TYPE_A_1FORM = 0,
  TYPE_B_1FORM,
  HIERARCHICAL_1FORM,
  MAX_FEM_FLAVORS_1FORM,    // Last flavor for FEM reconstruction
  COVOLUME_1FORM,
  MIMETIC_1FORM,
  CUSTOM_1FORM,
  MAX_FVD_FLAVORS_1FORM     // Last flavor for FVD reconstruction
};


/** \enum  Intrepid::ReconstructionFlavor_2Form
    \brief Enumerates the flavor of the reconstruction types available for 2-Forms.
*/
enum ReconstructionFlavor_2Form{
  TYPE_A_2FORM = 0,
  TYPE_B_2FORM,
  HIERARCHICAL_2FORM,
  MAX_FEM_FLAVORS_2FORM,    // Last flavor for FEM reconstruction
  COVOLUME_2FORM,
  MIMETIC_2FORM,
  CUSTOM_2FORM,
  MAX_FVD_FLAVORS_2FORM     // Last flavor for FVD reconstruction	
};


/** \enum  Intrepid::ReconstructionFlavor_3Form
    \brief Enumerates the flavor of the reconstruction types available for 3-Forms.
*/
enum ReconstructionFlavor_3Form{
  LAGRANGE_3FORM = 0,
  HIERARCHICAL_3FORM,
  MAX_FEM_FLAVORS_3FORM,    // Last flavor for FEM reconstruction
  COVOLUME_3FORM,
  MIMETIC_3FORM,
  CUSTOM_3FORM,
  MAX_FVD_FLAVORS_3FORM     // Last flavor for FVD reconstruction	
};


/** \enum  Intrepid::ErrorCode
    \brief Enumerates error codes in Intrepid.
*/
enum ErrorCode
{
  SUCCESS = 0,
  INDEX_OUT_OF_RANGE,
  TYPE_OUT_OF_RANGE,
  MEMORY_ALLOCATION_FAILED,
  ENTITY_NOT_FOUND,
  MULTIPLE_ENTITIES_FOUND,
  TAG_NOT_FOUND,
  FILE_DOES_NOT_EXIST,
  FILE_WRITE_ERROR,
  NOT_IMPLEMENTED,
  ALREADY_ALLOCATED,
  FAILURE
};


/** \enum  Intrepid::TagType
    \brief Enumerates tag types in the Intrepid mesh manager.
*/
enum TagType
{
  TAG_BIT = 0,
  TAG_SPARSE,
  TAG_DENSE,
  TAG_MESH,
  TAG_LAST
};


/** \enum  Intrepid::DataType
    \brief Enumerates data types in the Intrepid mesh manager.
*/
enum DataType
{
  TYPE_OPAQUE = 0,
  TYPE_INTEGER,
  TYPE_DOUBLE,
  TYPE_BIT,
  TYPE_HANDLE
};

/** \enum  Intrepid::StorageType
\brief Enumerates storage types used to pass various pieces of information in Intrepid.
*/
enum StorageType
{
  INTERLEAVED=0,
  BLOCKED
};


/** \enum  Intrepid::SetOp
    \brief Enumerates set operations in the Intrepid mesh manager.
*/
enum SetOp
{
  INTERSECT = 0,
  UNION
};


/** \typedef Intrepid::CellHandle
    \brief   CellHandle is a 32-bit integer handle.
*/
typedef unsigned int CellHandle;


/** \typedef Intrepid::Tag
    \brief   A Tag can be anything, so we define it as a void**.
*/
typedef void** Tag;

} //namespace Intrepid
#endif
