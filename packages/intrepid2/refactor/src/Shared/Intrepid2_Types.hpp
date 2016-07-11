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

/** \file
    \brief  Contains definitions of custom data types in Intrepid.
    \author Created by P. Bochev and D. Ridzal and Kyungjoo Kim.
*/

#ifndef __INTREPID2_TYPES_HPP__
#define __INTREPID2_TYPES_HPP__

#include <Kokkos_Core.hpp>
#include <Kokkos_DynRankView.hpp>

#include <stdexcept>

namespace Intrepid2 {

  // use ordinal_type and size_type everywhere (no index type)
  typedef int    ordinal_type;
  typedef size_t size_type;

  template<typename ValueType>
  KOKKOS_FORCEINLINE_FUNCTION
  ValueType epsilon() {
    return 0;
  }

  template<>
  KOKKOS_FORCEINLINE_FUNCTION
  double epsilon<double>() {
    typedef union {
      long long i64;
      double d64;
    } dbl_64;

    dbl_64 s;
    s.d64 = 1;
    s.i64++;
    return (s.i64 < 0 ? 1 - s.d64 : s.d64 - 1);
  }

  template<>
  KOKKOS_FORCEINLINE_FUNCTION
  float epsilon<float>() {
    typedef union {
      int i32;
      float f32;
    } flt_32;

    flt_32 s;
    s.f32 = 1;
    s.f32++;
    return (s.i32 < 0 ? 1 - s.f32 : s.f32 - 1);
  }

  KOKKOS_FORCEINLINE_FUNCTION
  double epsilon() {
    return epsilon<double>();
  }

  KOKKOS_FORCEINLINE_FUNCTION
  double tolerence() {
    return 100.0*epsilon();
  }

  KOKKOS_FORCEINLINE_FUNCTION
  double threshold() {
    return 10.0*epsilon();
  }

  /// define constants
  class Parameters {
  public:
    static constexpr unsigned int MaxNumPtsPerBasisEval= 16;      /// The maximum number of points to eval in serial mode
    static constexpr unsigned int MaxOrder             = 10;      /// The maximum reconstruction order.
    static constexpr unsigned int MaxIntegrationPoints = 1001;    /// The maximum number of integration points for direct cubature rules.
    static constexpr unsigned int MaxCubatureDegreeEdge= 20;      /// The maximum degree of the polynomial that can be integrated exactly by a direct edge rule.
    static constexpr unsigned int MaxCubatureDegreeTri = 20;      /// The maximum degree of the polynomial that can be integrated exactly by a direct triangle rule.
    static constexpr unsigned int MaxCubatureDegreeTet = 20;      /// The maximum degree of the polynomial that can be integrated exactly by a direct tetrahedron rule.
    static constexpr unsigned int MaxCubatureDegreePyr = 11;      /// The maximum degree of the polynomial that can be integrated exactly by a direct pyramid rule.
    static constexpr unsigned int MaxDimension         = 3;       /// The maximum ambient space dimension.
    static constexpr unsigned int MaxNewton            = 15;      /// Maximum number of Newton iterations used internally in methods such as computing the action of the inverse reference to physical cell map.
    static constexpr unsigned int MaxDerivative        = 10;      /// Maximum order of derivatives allowed in intrepid.

    // we do not want to use hard-wired epsilon, threshold and tolerence. 
    // static constexpr double Epsilon   = 1.0e-16; 
    // static constexpr double Threshold = 1.0e-15;
    // static constexpr double Tolerence = 1.0e-14;
  };
  //  const double Parameters::Epsilon   =       epsilon<double>();   /// Platform-dependent machine epsilon.
  //  const double Parameters::Threshold =  10.0*epsilon<double>();   /// Tolerance for various cell inclusion tests
  //  const double Parameters::Tolerence = 100.0*epsilon<double>();   /// General purpose tolerance in, e.g., internal Newton's method to invert ref to phys maps

  // ===================================================================
  // Enum classes
  // - Enum, String (char*) helper, valid
  // - can be used on device and inside of kernel for debugging purpose
  // - let's decorate kokkos inline

  /** \enum  PolyType
      \brief Enumeration of polynomial type; used in polylib
  */
  enum EPolyType {
    POLYTYPE_GAUSS=0,
    POLYTYPE_GAUSS_RADAU_LEFT,
    POLYTYPE_GAUSS_RADAU_RIGHT,
    POLYTYPE_GAUSS_LOBATTO,
    POLYTYPE_MAX
  };

  KOKKOS_INLINE_FUNCTION
  const char* EPolyTypeToString(const EPolyType polytype) {
    switch(polytype) {
    case POLYTYPE_GAUSS:                return "Gauss";
    case POLYTYPE_GAUSS_RADAU_LEFT:     return "GaussRadauLeft";
    case POLYTYPE_GAUSS_RADAU_RIGHT:    return "GaussRadauRight";
    case POLYTYPE_GAUSS_LOBATTO:        return "GaussRadauLobatto";
    case POLYTYPE_MAX:                  return "Max PolyType";
    default:                            return "INVALID EPolyType";
    }
    return "Error";
  }
  
  /** \brief  Verifies validity of a PolyType enum

      \param  polytype      [in]  - enum of the coordinate system
      \return 1 if the argument is valid poly type; 0 otherwise
  */
  KOKKOS_FORCEINLINE_FUNCTION
  bool isValidPolyType(const EPolyType polytype){
    return( polytype == POLYTYPE_GAUSS ||
            polytype == POLYTYPE_GAUSS_RADAU_LEFT ||
            polytype == POLYTYPE_GAUSS_RADAU_RIGHT ||
            polytype == POLYTYPE_GAUSS_LOBATTO );
  }


  /** \enum   Intrepid2::ECoordinates
      \brief  Enumeration of coordinate systems for geometrical entities (cells, points).
  */
  enum ECoordinates{
    COORDINATES_CARTESIAN=0,
    COORDINATES_POLAR,
    COORDINATES_CYLINDRICAL,
    COORDINATES_SPHERICAL,
    COORDINATES_MAX
  };

  KOKKOS_INLINE_FUNCTION
  const char* ECoordinatesToString(const ECoordinates coords) {
    switch(coords) {
    case COORDINATES_CARTESIAN:   return "Cartesian";
    case COORDINATES_POLAR:       return "Polar";
    case COORDINATES_CYLINDRICAL: return "Cylindrical";
    case COORDINATES_SPHERICAL:   return "Spherical";
    case COORDINATES_MAX:         return "Max. Coordinates";
    default:                      return "INVALID ECoordinates";
    }
    return "Error";
  }

  /** \brief  Verifies validity of a Coordinate enum.

      \param  coordinateType      [in]  - enum of the coordinate system
      \return 1 if the argument is valid coordinate system; 0 otherwise
  */
  KOKKOS_FORCEINLINE_FUNCTION
  bool isValidCoordinate(const ECoordinates coordinateType){
    return( coordinateType == COORDINATES_CARTESIAN   ||
            coordinateType == COORDINATES_POLAR       ||
            coordinateType == COORDINATES_CYLINDRICAL ||
            coordinateType == COORDINATES_SPHERICAL );
  }

  /** \enum   Intrepid2::ENorm
      \brief  Enumeration of norm types for vectors and functions
  */
  enum ENorm{
    NORM_ONE = 0,
    NORM_TWO,
    NORM_INF,
    NORM_FRO,    // Frobenius matrix norm
    NORM_MAX
  };

  KOKKOS_INLINE_FUNCTION
  const char* ENormToString(const ENorm norm) {
    switch(norm) {
    case NORM_ONE:   return "1-Norm";
    case NORM_TWO:   return "2-Norm";
    case NORM_INF:   return "Infinity Norm";
    case NORM_FRO:   return "Frobenius Norm";
    case NORM_MAX:   return "Max. Norm";
    default:         return "INVALID ENorm";
    }
    return "Error";
  }

  /** \brief  Verifies validity of a Norm enum.

      \param  normType      [in]  - enum of the norm
      \return 1 if the argument is valid norm; 0 otherwise
  */
  KOKKOS_FORCEINLINE_FUNCTION
  bool isValidNorm(const ENorm normType){
    return( normType == NORM_ONE ||
            normType == NORM_TWO ||
            normType == NORM_INF ||
            normType == NORM_FRO ||
            normType == NORM_MAX );
  }

  /** \enum   Intrepid2::EOperator
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
    OPERATOR_Dn,        // 14
    OPERATOR_MAX = OPERATOR_Dn // 14
  };

  KOKKOS_INLINE_FUNCTION
  const char* EOperatorToString(const EOperator op) {
    switch(op) {
    case OPERATOR_VALUE: return "Value";
    case OPERATOR_GRAD:  return "Grad";
    case OPERATOR_CURL:  return "Curl";
    case OPERATOR_DIV:   return "Div";
    case OPERATOR_D1:    return "D1";
    case OPERATOR_D2:    return "D2";
    case OPERATOR_D3:    return "D3";
    case OPERATOR_D4:    return "D4";
    case OPERATOR_D5:    return "D5";
    case OPERATOR_D6:    return "D6";
    case OPERATOR_D7:    return "D7";
    case OPERATOR_D8:    return "D8";
    case OPERATOR_D9:    return "D9";
    case OPERATOR_D10:   return "D10";
    case OPERATOR_MAX:   return "Dn Operator";
    default:             return "INVALID EOperator";
    }
    return "Error";
  }

  /** \brief  Verifies validity of an operator enum.

      \param  operatorType      [in]  - enum of the operator
      \return 1 if the argument is valid operator; 0 otherwise
  */
  KOKKOS_FORCEINLINE_FUNCTION
  bool isValidOperator(const EOperator operatorType){
    return ( operatorType == OPERATOR_VALUE ||
             operatorType == OPERATOR_GRAD  ||
             operatorType == OPERATOR_CURL  ||
             operatorType == OPERATOR_DIV   ||
             operatorType == OPERATOR_D1    ||
             operatorType == OPERATOR_D2    ||
             operatorType == OPERATOR_D3    ||
             operatorType == OPERATOR_D4    ||
             operatorType == OPERATOR_D5    ||
             operatorType == OPERATOR_D6    ||
             operatorType == OPERATOR_D7    ||
             operatorType == OPERATOR_D8    ||
             operatorType == OPERATOR_D9    ||
             operatorType == OPERATOR_D10 );
  }


  /** \enum   Intrepid2::FunctionSpace
      \brief  Enumeration of the admissible function space types in Intrepid.
  */
  enum EFunctionSpace {
    FUNCTION_SPACE_HGRAD = 0,
    FUNCTION_SPACE_HCURL,
    FUNCTION_SPACE_HDIV,
    FUNCTION_SPACE_HVOL,
    FUNCTION_SPACE_VECTOR_HGRAD,
    FUNCTION_SPACE_TENSOR_HGRAD,
    FUNCTION_SPACE_MAX
  };

  KOKKOS_INLINE_FUNCTION
  const char* EFunctionSpaceToString(const EFunctionSpace space) {
    switch(space) {
    case FUNCTION_SPACE_HGRAD:        return "H(grad)";
    case FUNCTION_SPACE_HCURL:        return "H(curl)";
    case FUNCTION_SPACE_HDIV:         return "H(div)";
    case FUNCTION_SPACE_HVOL:         return "H(vol)";
    case FUNCTION_SPACE_VECTOR_HGRAD: return "Vector H(grad)";
    case FUNCTION_SPACE_TENSOR_HGRAD: return "Tensor H(grad)";
    case FUNCTION_SPACE_MAX:          return "Max. Function space";
    default:                          return "INVALID EFunctionSpace";
    }
    return "Error";
  }

  /** \brief  Verifies validity of a function space enum

      \param  spaceType      [in]  - enum of the function space
      \return 1 if the argument is valid function space; 0 otherwise
  */
  KOKKOS_FORCEINLINE_FUNCTION
  bool isValidFunctionSpace(const EFunctionSpace spaceType){
    return ( spaceType == FUNCTION_SPACE_HGRAD ||
             spaceType == FUNCTION_SPACE_HCURL ||
             spaceType == FUNCTION_SPACE_HDIV  ||
             spaceType == FUNCTION_SPACE_HVOL  ||
             spaceType == FUNCTION_SPACE_VECTOR_HGRAD ||
             spaceType == FUNCTION_SPACE_TENSOR_HGRAD );
  }

  /** \enum   Intrepid2::EDiscreteSpace
      \brief  Enumeration of the discrete spaces used to define bases for function spaces.
      Intrepid allows up to three basic kinds of discrete spaces for each cell type.

      \arg    COMPLETE     complete polynomial or tensor product space
      \arg    INCOMPLETE   incomplete polynomial or tensor product space, such as used in RT elements
      \arg    BROKEN       piecewise smooth, with respect to a cell, polynomial space
  */
  enum EDiscreteSpace {
    DISCRETE_SPACE_COMPLETE = 0,        // value = 0
    DISCRETE_SPACE_INCOMPLETE,          // value = 1
    DISCRETE_SPACE_BROKEN,              // value = 2
    DISCRETE_SPACE_MAX                  // value = 3
  };

  KOKKOS_INLINE_FUNCTION
  const char* EDiscreteSpaceToString(const EDiscreteSpace space) {
    switch(space) {
    case DISCRETE_SPACE_COMPLETE:   return "Complete";
    case DISCRETE_SPACE_INCOMPLETE: return "Incomplete";
    case DISCRETE_SPACE_BROKEN:     return "Broken";
    case DISCRETE_SPACE_MAX:        return "Max. Rec. Space";
    default:                        return "INVALID EDiscreteSpace";
    }
    return "Error";
  }

  /** \brief  Verifies validity of a discrete space enum

      \param  spaceType      [in]  - enum of the function space
      \return 1 if the argument is valid discrete space; 0 otherwise
  */
  KOKKOS_FORCEINLINE_FUNCTION
  bool isValidDiscreteSpace(const EDiscreteSpace spaceType){
    return ( spaceType == DISCRETE_SPACE_COMPLETE   ||
             spaceType == DISCRETE_SPACE_INCOMPLETE ||
             spaceType == DISCRETE_SPACE_BROKEN );
  }

  /** \enum   Intrepid2::EPointType
      \brief  Enumeration of types of point distributions in Intrepid
  */
  enum EPointType {
    POINTTYPE_EQUISPACED = 0,             // value = 0
    POINTTYPE_WARPBLEND,
    POINTTYPE_GAUSS
  };

  KOKKOS_INLINE_FUNCTION
  const char* EPointTypeToString(const EPointType pointType) {
    switch (pointType) {
    case POINTTYPE_EQUISPACED:     return "Equispaced Points";
    case POINTTYPE_WARPBLEND:      return "WarpBlend Points";
    case POINTTYPE_GAUSS:          return "Gauss Points";
    default:                       return "INVALID EPointType";
    }
    return "Error";
  }

  /** \brief Verifies validity of a point type enum
      \param pointType      [in] - enum of the point type
      \return 1 if the argument is a valid point type; 0 otherwise
  */
  KOKKOS_FORCEINLINE_FUNCTION
  bool isValidPointType(const EPointType pointType) {
    return ( pointType == POINTTYPE_EQUISPACED ||
             pointType == POINTTYPE_WARPBLEND ||
             pointType == POINTTYPE_GAUSS );
  }

  /** \enum   Intrepid2::EBasis
      \brief  Enumeration of basis types for discrete spaces in Intrepid.
  */
  enum EBasis {
    BASIS_FEM_DEFAULT = 0,                // value = 0
    BASIS_FEM_HIERARCHICAL,               // value = 1
    BASIS_FEM_FIAT,                       // value = 2
    BASIS_FVD_DEFAULT,                    // value = 3
    BASIS_FVD_COVOLUME,                   // value = 4
    BASIS_FVD_MIMETIC,                    // value = 5
    BASIS_MAX                             // value = 6
  };

  KOKKOS_INLINE_FUNCTION
  const char* EBasisToString(const EBasis basis) {
    switch(basis) {
    case BASIS_FEM_DEFAULT:      return "FEM Default";
    case BASIS_FEM_HIERARCHICAL: return "FEM Hierarchical";
    case BASIS_FEM_FIAT:         return "FEM FIAT";
    case BASIS_FVD_DEFAULT:      return "FVD Default";
    case BASIS_FVD_COVOLUME:     return "FVD Covolume";
    case BASIS_FVD_MIMETIC:      return "FVD Mimetic";
    case BASIS_MAX:              return "Max. Basis";
    default:                     return "INVALID EBasis";
    }
    return "Error";
  }

  /** \brief  Verifies validity of a basis enum

      \param  basisType      [in]  - enum of the basis
      \return 1 if the argument is valid discrete space; 0 otherwise
  */
  KOKKOS_FORCEINLINE_FUNCTION
  bool isValidBasis(const EBasis basisType){
    return ( basisType == BASIS_FEM_DEFAULT      ||
             basisType == BASIS_FEM_HIERARCHICAL ||
             basisType == BASIS_FEM_FIAT         ||
             basisType == BASIS_FVD_DEFAULT      ||
             basisType == BASIS_FVD_COVOLUME     ||
             basisType == BASIS_FVD_MIMETIC );
  }

  // /** \enum  Intrepid2::ECompEngine
  //     \brief Specifies how operators and functionals are computed internally
  //     (COMP_MANUAL = native C++ implementation, COMP_BLAS = BLAS implementation, etc.).
  // */
  // enum ECompEngine {
  //   COMP_CPP = 0,
  //   COMP_BLAS,
  //   COMP_ENGINE_MAX
  // };

  // KOKKOS_INLINE_FUNCTION
  // const char* ECompEngineToString(const ECompEngine cEngine) {
  //   switch(cEngine) {
  //   case COMP_CPP:             return "Native C++";
  //   case COMP_BLAS:            return "BLAS";
  //   case COMP_ENGINE_MAX:      return "Max. Comp. Engine";
  //   default:                   return "INVALID ECompEngine";
  //   }
  //   return "Error";
  // }


  // /** \brief  Verifies validity of a computational engine enum

  //     \param  compEngType    [in]  - enum of the computational engine
  //     \return 1 if the argument is valid computational engine; 0 otherwise
  // */
  // KOKKOS_FORCEINLINE_FUNCTION
  // bool isValidCompEngine(const ECompEngine compEngType){
  //   //at the moment COMP_BLAS is not a valid CompEngine.
  //   return (compEngType == COMP_CPP);
  // }


} //namespace Intrepid2

#endif
