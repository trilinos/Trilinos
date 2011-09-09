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

/** \file
    \brief  Contains definitions of custom data types in Intrepid.
    \author Created by P. Bochev and D. Ridzal.
 */

#ifndef INTREPID_INTREPID_TYPES_HPP
#define INTREPID_INTREPID_TYPES_HPP

#ifdef  HAVE_INTREPID_DEBUG
#define INTREPID_VALIDATE( A )  A
#else
#define INTREPID_VALIDATE( A ) /* empty */
#endif

#include <Teuchos_ScalarTraits.hpp>

/** \def    INTREPID_MAX_ORDER
    \brief  The maximum reconstruction order.
 */
#define INTREPID_MAX_ORDER 10

/** \def    INTREPID_MAX_INTEGRATION_POINTS
    \brief  The maximum number of integration points for direct cubature rules.
 */
#define INTREPID_MAX_INTEGRATION_POINTS 1001

/** \def    INTREPID_MAX_CUBATURE_DEGREE_EDGE
    \brief  The maximum degree of the polynomial that can be integrated exactly by
            a direct edge rule.
 */
#define INTREPID_MAX_CUBATURE_DEGREE_EDGE 61

/** \def    INTREPID_MAX_CUBATURE_DEGREE_TRI
    \brief  The maximum degree of the polynomial that can be integrated exactly by
            a direct triangle rule.
 */
#define INTREPID_MAX_CUBATURE_DEGREE_TRI 20

/** \def    INTREPID_MAX_CUBATURE_DEGREE_TET
    \brief  The maximum degree of the polynomial that can be integrated exactly by
            a direct tetrahedron rule.
 */
#define INTREPID_MAX_CUBATURE_DEGREE_TET 20

/** \def    INTREPID_MAX_CUBATURE_DEGREE_PYRAMID
    \brief  The maximum degree of the polynomial that can be integrated exactly by
            a direct pyramid rule.
 */
#define INTREPID_MAX_CUBATURE_DEGREE_PYRAMID 4

/** \def    INTREPID_MAX_DIMENSION
    \brief  The maximum ambient space dimension.
 */
#define INTREPID_MAX_DIMENSION 3

/** \def    INTREPID_MAX_NEWTON 
    \brief  Maximum number of Newton iterations used internally in methods such as computing the
            action of the inverse reference to physical cell map.
 */
#define INTREPID_MAX_NEWTON 15

/** \def    INTREPID_MAX_DERIVATIVE
    \brief  Maximum order of derivatives allowed in intrepid
*/
#define INTREPID_MAX_DERIVATIVE 10

namespace Intrepid {
  
  /** \brief  Platform-dependent machine epsilon. 
   */
  static const double INTREPID_EPSILON   = std::abs(Teuchos::ScalarTraits<double>::eps());
    
  /** \brief  Tolerance for various cell inclusion tests
   */
  static const double INTREPID_THRESHOLD = 10.0 * INTREPID_EPSILON;
  
  /** \brief  General purpose tolerance in, e.g., internal Newton's method to invert ref to phys maps
   */
  static const double INTREPID_TOL       = 10.0* INTREPID_THRESHOLD;
  
  /** \enum   Intrepid::ECoordinates
      \brief  Enumeration of coordinate systems for geometrical entities (cells, points).
   */
  enum ECoordinates{
    COORDINATES_CARTESIAN=0,
    COORDINATES_POLAR,
    COORDINATES_CYLINDRICAL,
    COORDINATES_SPHERICAL,
    COORDINATES_MAX
  };

  inline std::string ECoordinatesToString(ECoordinates coords) {
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

  /** \brief  Verifies validity of a Coordinate enum.
    
      \param  coordinateType      [in]  - enum of the coordinate system
      \return 1 if the argument is valid coordinate system; 0 otherwise
    */
  inline int isValidCoordinate(ECoordinates coordinateType){
    return( ( coordinateType == COORDINATES_CARTESIAN)   ||
            ( coordinateType == COORDINATES_POLAR)       ||
            ( coordinateType == COORDINATES_CYLINDRICAL) ||
            ( coordinateType == COORDINATES_SPHERICAL) );
  }
  
  
  
  /** \enum   Intrepid::ENorm
      \brief  Enumeration of norm types for vectors and functions
   */
  enum ENorm{
    NORM_ONE = 0,
    NORM_TWO,
    NORM_INF,
    NORM_FRO,    // Frobenius matrix norm
    NORM_MAX
  };

  inline std::string ENormToString(ENorm norm) {
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
  
  /** \brief  Verifies validity of a Norm enum.
    
      \param  normType      [in]  - enum of the norm
      \return 1 if the argument is valid norm; 0 otherwise
    */
  inline int isValidNorm(ENorm normType){
    return( (normType == NORM_ONE) ||
            (normType == NORM_TWO) ||
            (normType == NORM_INF) ||
            (normType == NORM_FRO) ||
            (normType == NORM_MAX) );
  }


  
  /** \enum   Intrepid::EOperator
      \brief  Enumeration of primitive operators available in Intrepid. Primitive operators act on
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
  
  inline std::string EOperatorToString(EOperator op) {
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

  /** \brief  Verifies validity of an operator enum.
    
      \param  operatorType      [in]  - enum of the operator
      \return 1 if the argument is valid operator; 0 otherwise
    */
  inline int isValidOperator(const EOperator operatorType){
    return ( (operatorType == OPERATOR_VALUE) ||
             (operatorType == OPERATOR_GRAD)  || 
             (operatorType == OPERATOR_CURL)  || 
             (operatorType == OPERATOR_DIV)   ||
             (operatorType == OPERATOR_D1)    || 
             (operatorType == OPERATOR_D2)    || 
             (operatorType == OPERATOR_D3)    || 
             (operatorType == OPERATOR_D4)    || 
             (operatorType == OPERATOR_D5)    || 
             (operatorType == OPERATOR_D6)    || 
             (operatorType == OPERATOR_D7)    || 
             (operatorType == OPERATOR_D8)    || 
             (operatorType == OPERATOR_D9)    || 
             (operatorType == OPERATOR_D10) );
  }
  
  
  /** \enum   Intrepid::FunctionSpace
      \brief  Enumeration of the admissible function space types in Intrepid.
    */
  enum EFunctionSpace
    {
      FUNCTION_SPACE_HGRAD = 0,
      FUNCTION_SPACE_HCURL,
      FUNCTION_SPACE_HDIV,
      FUNCTION_SPACE_HVOL,
      FUNCTION_SPACE_VECTOR_HGRAD,
      FUNCTION_SPACE_TENSOR_HGRAD,
      FUNCTION_SPACE_MAX
    };
  
  inline std::string EFunctionSpaceToString(EFunctionSpace space) {
    std::string retString;
    switch(space) {
      case FUNCTION_SPACE_HGRAD:        retString = "H(grad)";     break;
      case FUNCTION_SPACE_HCURL:        retString = "H(curl)";     break;
      case FUNCTION_SPACE_HDIV:         retString = "H(div)";     break;
      case FUNCTION_SPACE_HVOL:         retString = "H(vol)";     break;
      case FUNCTION_SPACE_VECTOR_HGRAD: retString = "Vector H(grad)";     break;
      case FUNCTION_SPACE_TENSOR_HGRAD: retString = "Tensor H(grad)";     break;
      case FUNCTION_SPACE_MAX:          retString = "Max. Function space"; break;
      default:                          retString = "INVALID EFunctionSpace";
    }
    return retString;
  }
  
  /** \brief  Verifies validity of a function space enum
    
      \param  spaceType      [in]  - enum of the function space
      \return 1 if the argument is valid function space; 0 otherwise
    */
  inline int isValidFunctionSpace(const EFunctionSpace spaceType){
    return ( (spaceType == FUNCTION_SPACE_HGRAD) ||
             (spaceType == FUNCTION_SPACE_HCURL) || 
             (spaceType == FUNCTION_SPACE_HDIV)  ||
             (spaceType == FUNCTION_SPACE_HVOL)  ||
             (spaceType == FUNCTION_SPACE_VECTOR_HGRAD) || 
             (spaceType == FUNCTION_SPACE_TENSOR_HGRAD) );
  }
  
  
  
  /** \enum   Intrepid::EDiscreteSpace
      \brief  Enumeration of the discrete spaces used to define bases for function spaces. 
              Intrepid allows up to three basic kinds of discrete spaces for each cell type.
    
      \arg    COMPLETE     complete polynomial or tensor product space
      \arg    INCOMPLETE   incomplete polynomial or tensor product space, such as used in RT elements
      \arg    BROKEN       piecewise smooth, with respect to a cell, polynomial space
    */
  enum EDiscreteSpace
    {
      DISCRETE_SPACE_COMPLETE = 0,        // value = 0
      DISCRETE_SPACE_INCOMPLETE,          // value = 1
      DISCRETE_SPACE_BROKEN,              // value = 2
      DISCRETE_SPACE_MAX                  // value = 3
    };
  
  inline std::string EDiscreteSpaceToString(EDiscreteSpace space) {
    std::string retString;
    switch(space) {
      case DISCRETE_SPACE_COMPLETE:   retString = "Complete";        break;
      case DISCRETE_SPACE_INCOMPLETE: retString = "Incomplete";      break;
      case DISCRETE_SPACE_BROKEN:     retString = "Broken";          break;
      case DISCRETE_SPACE_MAX:        retString = "Max. Rec. Space"; break;
      default:                              retString = "INVALID EDiscreteSpace";
    }
    return retString;
  }
  
  /** \brief  Verifies validity of a discrete space enum
    
      \param  spaceType      [in]  - enum of the function space
      \return 1 if the argument is valid discrete space; 0 otherwise
    */
  inline int isValidDiscreteSpace(const EDiscreteSpace spaceType){
    return ( (spaceType == DISCRETE_SPACE_COMPLETE) ||
             (spaceType == DISCRETE_SPACE_INCOMPLETE) || 
             (spaceType ==DISCRETE_SPACE_BROKEN) );
  }

  /** \enum   Intrepid::EPointType
      \brief  Enumeration of types of point distributions in Intrepid
    */
  enum EPointType
    {
      POINTTYPE_EQUISPACED = 0,             // value = 0
      POINTTYPE_SPECTRAL,
      POINTTYPE_SPECTRAL_OPEN,
      POINTTYPE_WARPBLEND 
    };
  
  inline std::string EPointTypeToString(EPointType pointType) {
    std::string retString;
    switch (pointType) {
    case POINTTYPE_EQUISPACED:
      retString = "Equispaced Points";
      break;
    case POINTTYPE_WARPBLEND:
      retString = "WarpBlend Points";
      break;
    case POINTTYPE_SPECTRAL:
      retString = "Spectral Points";
      break;
    case POINTTYPE_SPECTRAL_OPEN:
      retString = "Open Spectral Points";
      break;
    }
    return retString;
  }
  
  /** \brief Verifies validity of a point type enum
      \param pointType      [in] - enum of the point type
      \return 1 if the argument is a valid point type; 0 otherwise
   */

  inline int isValidPointType( const EPointType pointType ) {
    return ( (pointType == POINTTYPE_EQUISPACED ) ||
           (pointType == POINTTYPE_WARPBLEND ) );
  }

  /** \enum   Intrepid::EBasis
      \brief  Enumeration of basis types for discrete spaces in Intrepid.
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
  
  inline std::string EBasisToString(EBasis basis) {
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
  
  /** \brief  Verifies validity of a basis enum
    
      \param  basisType      [in]  - enum of the basis
      \return 1 if the argument is valid discrete space; 0 otherwise
    */
  inline int isValidBasis(const EBasis basisType){
    return ( (basisType == BASIS_FEM_DEFAULT) ||
             (basisType == BASIS_FEM_HIERARCHICAL) ||
             (basisType == BASIS_FEM_FIAT) ||
             (basisType == BASIS_FVD_DEFAULT) ||
             (basisType == BASIS_FVD_COVOLUME) ||
             (basisType == BASIS_FVD_MIMETIC) );
  }

  
  /** \struct Intrepid::CubatureTemplate
      \brief  Template for the cubature rules used by Intrepid. Cubature template consists of  
              cubature points and cubature weights. Intrepid provides a collection of cubature  
              templates for most standard cell topologies. The templates are defined in reference  
              coordinates using a standard reference cell for each canonical cell type. Cubature 
              points are always specified by a triple of (X,Y,Z) coordinates even if the cell 
              dimension is less than 3. The unused dimensions should be padded by zeroes.
    
              For example, a set of Gauss rules on [-1,1] looks as the following array of CubatureTemplate structs:
    
  \verbatim
  cubature_rule[4] =
  {                                         // Collection of Gauss rules on [-1,1]
    {
      1,                                      ----> number of points in the rule
      {{0.0,0.0,0.0}},                        ----> X,Y,Z coordinates of the cubature points  
      {0.5}                                   ----> the cubature weight
    },
    {
      2,
      {{-sqrt(1.0/3.0),0.0,0.0},
       {+sqrt(1.0/3.0),0.0,0.0}},
      {1.0,1.0}
    },
    {
      3,
      {{-sqrt(3.0/5.0),0.0,0.0},
       {0.0,0.0,0.0},
       {+sqrt(3.0/5.0),0.0,0.0}},
      {5.0/9.0, 8.0/9.0,5.0/9.0}
    },
    {
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
    }  
  };  // end Gauss
  \endverbatim

      Also see data member documentation.
  */
  struct CubatureTemplate {

    /** \brief  Number of cubature points stored in the template.
     */  
    int               numPoints_;
    
    /** \brief  Array with the (X,Y,Z) coordinates of the cubature points.
     */
    double            points_[INTREPID_MAX_INTEGRATION_POINTS][INTREPID_MAX_DIMENSION];
    
    /** \brief  Array with the associated cubature weights.
     */
    double            weights_[INTREPID_MAX_INTEGRATION_POINTS];
    
  };



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

  inline std::string ECompEngineToString(ECompEngine cEngine) {
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
  
  
  /** \brief  Verifies validity of a computational engine enum
    
      \param  compEngType    [in]  - enum of the computational engine
      \return 1 if the argument is valid computational engine; 0 otherwise
   */
  inline int isValidCompEngine(const ECompEngine compEngType){
    return ( (compEngType == COMP_CPP) ||
             (compEngType == COMP_BLAS) );
  }

} //namespace Intrepid

/*! \mainpage INTREPID Documentation (Development Version)
 
  \image html intrepid.png
  \image latex intrepid.jpg "Reconnaissance balloon ``Intrepid''" width=1in

  \section intro_sec Introduction

  %Intrepid is a library of interoperable tools for compatible discretizations of
  Partial Differential Equations (PDEs). Included with the Trilinos 10.0 release
  is the &quot;<em>expert version</em>&quot; of %Intrepid. This version is intended
  primarily for application developers who want to reuse large parts of their existing
  code frameworks such as I/O, data structures, assembly routines, etc. while gaining
  access to advanced discretization capabilities provided by %Intrepid. In such cases
  the bulk of the data is owned and managed by the user rather than by %Intrepid.
  To avoid unnecessary and possibly expensive copying of data to and from %Intrepid,
  the expert version of the package comprises of mostly stateless classes operating on
  user-owned data. Virtually all numerical data required by PDE codes can be represented
  as a multi-dimensional array of scalar values. For this reason, and to enhance
  interoprability, %Intrepid classes are templated on generic multi-dimensional arrays.
  The <a href="http://trilinos.sandia.gov/packages/shards/">Shards</a> package provides
  an implementation of a multi-dimensional array that can be used for that purpose, or
  users can write their own multi-dimensional arrays as long as a minimal interface
  is supported.
 
  \section overview_sec Overview

  Current release of %Intrepid includes the following features:
  \li Default finite element basis functions for <em>H(grad)</em>, <em>H(curl)</em>,
       <em>H(div)</em> and <em>L2</em> spaces of orders up to 2 on standard cell
       topologies in 1D, 2D and 3D</li>
  \li High-order (up to 10) basis functions for <em>H(grad)</em>, <em>H(curl)</em>,
       <em>H(div)</em> and <em>L2</em> spaces on select cell topologies</li>
  \li Pullbacks (transformations) from reference coordinate frame of <em>H(grad)</em>,
       <em>H(curl)</em>, <em>H(div)</em> and <em>L2</em> fields</li>
  \li Pullbacks of gradient, curl and divergence of <em>H(grad)</em>, <em>H(curl)</em>,
       <em>H(div)</em> fields</li>
  \li Cubature rules of orders up to 20 on most standard 1D, 2D and 3D cell topologies</li>
  \li Implementation of multi-diumensional arrays and algebraic operations on them</li>
  \li Examples showing solution of basic 2nd order elliptic boundary value problems
       (Poisson, div-curl, and curl-curl systems) using %Intrepid</li>

  \section quickstart_sec Quick Start

  Familiarity with with the following concepts, objects, and tools is required:
  \li <a href="http://trilinos.sandia.gov/packages/shards/">Shards</a> cell topologies,
  \li numerical integration / Intrepid::Cubature,
  \li discrete (e.g. finite element) bases / Intrepid::Basis / \ref basis_page,
  \li multi-dimensional arrays / Intrepid::FieldContainer / \ref md_array_page,
  \li cell mappings and transformations / Intrepid::CellTools / \ref cell_tools_page, and
  \li function mappings (pullbacks) / Intrepid::FunctionSpaceTools / \ref function_space_tools_page.

  The following example demonstrates, in 7 steps, the computation of finite element
  stiffness matrices on a set of tetrahedral cells using a piecewise linear basis
  and an appropriate integration rule.

  \subsection topo_qs_sec Step 1: Select a cell topology
  
  \code
      shards::CellTopology cellType = shards::getCellTopologyData< shards::Tetrahedron<> >(); // cell type: tetrahedron
      int spaceDim = cellType->getDimension();                                                // retrieve spatial dimension
      int numNodes = cellType->getNodeCount();                                                // retrieve number of 0-cells (nodes)
  \endcode

  We additionally set the number of computational cells \c numCells.


  \subsection integration_qs_sec Step 2: Select integration (cubature) rule

  \code
      DefaultCubatureFactory<double> cubFactory;                                              // create cubature factory
      int cubDegree = 2;                                                                      // set cubature degree, e.g. 2
      Teuchos::RCP<Cubature<double> > myCub = cubFactory.create(cellType, cubDegree);         // create default cubature
      int numCubPoints = myCub->getNumPoints();                                               // retrieve number of cubature points
  \endcode


  \subsection bases_qs_sec Step 3: Select discrete basis
  
  \code
      Basis_HGRAD_TET_C1_FEM<double, FieldContainer<double> > tetBasis;                       // create tet basis
      int numFields = tetBasis.getCardinality();                                              // get basis cardinality
  \endcode


  \subsection mdarray_qs_sec Step 4: Format multi-dimensional arrays

  \code
      FieldContainer<double> cub_points(numCubPoints, spaceDim);
      FieldContainer<double> cub_weights(numCubPoints);

      FieldContainer<double> cell_nodes(numCells, numNodes, spaceDim);

      FieldContainer<double> jacobian(numCells, numCubPoints, spaceDim, spaceDim);
      FieldContainer<double> jacobian_inv(numCells, numCubPoints, spaceDim, spaceDim);
      FieldContainer<double> jacobian_det(numCells, numCubPoints);
      FieldContainer<double> weighted_measure(numCells, numCubPoints);

      FieldContainer<double> grad_at_cub_points(numFields, numCubPoints, spaceDim);
      FieldContainer<double> transformed_grad_at_cub_points(numCells, numFields, numCubPoints, spaceDim);
      FieldContainer<double> weighted_transformed_grad_at_cub_points(numCells, numFields, numCubPoints, spaceDim);
      FieldContainer<double> stiffness_matrices(numCells, numFields, numFields);
  \endcode

  We assume that the array \c cell_nodes is filled with nodes defining a set of computational (physical) cells.

  \subsection tabulate_qs_sec Step 5: Evaluate differential operator applied to basis at cubature points 

  \code
      myCub->getCubature(cub_points, cub_weights);                                          // retrieve cubature points and weights
      tetBasis.getValues(grad_at_cub_points, cub_points, OPERATOR_GRAD);                    // evaluate grad operator at cubature points
  \endcode


  \subsection ct_qs_sec Step 6: Apply cell tools

  \code
      CellTools<double>::setJacobian(jacobian, cub_points, cell_nodes, cellType);           // compute cell Jacobians
      CellTools<double>::setJacobianInv(jacobian_inv, jacobian);                            // compute inverses of cell Jacobians
      CellTools<double>::setJacobianDet(jacobian_det, jacobian);                            // compute determinants of cell Jacobians
  \endcode


  \subsection fst_qs_sec Step 7: Apply function space tools

  \code
      FunctionSpaceTools::computeCellMeasure<double>(weighted_measure,                      // compute weighted cell measure
                                                     jacobian_det,
                                                     cub_weights);
      FunctionSpaceTools::HGRADtransformGRAD<double>(transformed_grad_at_cub_points,        // transform reference gradients into physical space
                                                     jacobian_inv,
                                                     grad_at_cub_points);
      FunctionSpaceTools::multiplyMeasure<double>(weighted_transformed_grad_at_cub_points,  // multiply with weighted measure
                                                  weighted_measure,
                                                  transformed_grad_at_cub_points);
      FunctionSpaceTools::integrate<double>(stiffness_matrices,                             // compute stiffness matrices
                                            transformed_grad_at_cub_points,
                                            weighted_transformed_grad_at_cub_points,
                                            COMP_CPP);
  \endcode

  The computed (local) stiffness matrices can now be used in the assembly of a (global)
  discrete differential operator, e.g. a discrete Laplacian.

  \subsection doen_qs_sec Done!


  \section devplans_sec Development Plans

  The next release of %Intrepid is expected to support Finite Difference and Finite Volume
  discretizations on standard and non-standard (polygon and polyhedron) cell topologies.
  A &quot;<em>user-friendly</em>&quot; version for rapid development of PDE codes is also
  under development.
*/
  

#endif
