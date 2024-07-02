// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/**
   \file   Amesos2_NewSolver_TypeMap.hpp
   \author John Doe <jd@sandia.gov>
   \date   
   
   \brief Provides definition of NewSolver types as well as
          conversions and type traits.  For the purpose of
          demonstration, we assume that NewSolver has defined its own
          complex data-types called `complex' and `doublecomplex'.
*/

#ifndef AMESOS2_NEWSOLVER_TYPEMAP_HPP
#define AMESOS2_NEWSOLVER_TYPEMAP_HPP

#ifdef HAVE_TEUCHOS_COMPLEX
#include <complex>
#endif

#include <Teuchos_as.hpp>

#include "Amesos2_TypeMap.hpp"

// If necessary, include headers from the NewSolver library that
// declare NewSolver-specific type definitions.
namespace New_Solver {
  extern "C" {

    namespace C {
      // include solver-specific single-precision complex data type definitions
    }

    namespace Z {
      // include solver-specific double-precision complex data type definitions
    }

    // include other solver data type definitions as needed

  } // end extern "C"
} // end namespace New_Solver


/* ==================== Conversion ====================
 *
 * Define here, in the Teuchos namespace, any conversions between commonly
 * used date types and the solver-specific data types.  Use template
 * specializations of the Teuchos::ValueTypeConversionTraits class.
 */
namespace Teuchos {

  /**
   * \defgroup new_solver_conversion Conversion definitions for NewSolver types.
   * 
   * Define specializations of Teuchos::as<> for the NewSolver types.
   *
   * These specializations are meant to work with any complex data type that
   * implements the same interface as the STL complex type.
   *
   * @{
   */
  template <typename TypeFrom>
  class ValueTypeConversionTraits<New_Solver::C::complex, TypeFrom>
  {
  public:
    static New_Solver::C::complex convert( const TypeFrom t )
    {                           // adapt conversion as necessary
      New_Solver::C::complex ret;
      ret.r = Teuchos::as<float>(t.real());
      ret.i = Teuchos::as<float>(t.imag());
      return( ret );
    }
  
    static New_Solver::C::complex safeConvert( const TypeFrom t )
    {                           // adapt conversion as necessary
      New_Solver::C::complex ret;
      ret.r = Teuchos::as<float>(t.real());
      ret.i = Teuchos::as<float>(t.imag());
      return( ret );
    }
  };


  template <typename TypeFrom>
  class ValueTypeConversionTraits<New_Solver::Z::doublecomplex, TypeFrom>
  {
  public:
    static New_Solver::Z::doublecomplex convert( const TypeFrom t )
    {                           // adapt conversion as necessary
      New_Solver::Z::doublecomplex ret;
      ret.r = Teuchos::as<double>(t.real());
      ret.i = Teuchos::as<double>(t.imag());
      return( ret );
    }
  
    static New_Solver::Z::doublecomplex safeConvert( const TypeFrom t )
    {                           // adapt conversion as necessary
      New_Solver::Z::doublecomplex ret;
      ret.r = Teuchos::as<double>(t.real());
      ret.i = Teuchos::as<double>(t.imag());
      return( ret );
    }
  };


  // Also convert *from* New_Solver types
  template <typename TypeTo>
  class ValueTypeConversionTraits<TypeTo, New_Solver::C::complex>
  {
  public:
    static TypeTo convert( const New_Solver::C::complex t )
    {                           // adapt conversion as necessary
      typedef typename TypeTo::value_type value_type;
      value_type ret_r = Teuchos::as<value_type>( t.r );
      value_type ret_i = Teuchos::as<value_type>( t.i );
      return ( TypeTo( ret_r, ret_i ) );
    }
  
    static TypeTo safeConvert( const New_Solver::C::complex t )
    {                           // adapt conversion as necessary
      typedef typename TypeTo::value_type value_type;
      value_type ret_r = Teuchos::as<value_type>( t.r );
      value_type ret_i = Teuchos::as<value_type>( t.i );
      return ( TypeTo( ret_r, ret_i ) );
    }
  };


  template <typename TypeTo>
  class ValueTypeConversionTraits<TypeTo, New_Solver::Z::doublecomplex>
  {
  public:
    static TypeTo convert( const New_Solver::Z::doublecomplex t )
    {
      typedef typename TypeTo::value_type value_type;
      value_type ret_r = Teuchos::as<value_type>( t.r );
      value_type ret_i = Teuchos::as<value_type>( t.i );
      return ( TypeTo( ret_r, ret_i ) );
    }
  
    // No special checks for safe Convert
    static TypeTo safeConvert( const New_Solver::Z::doublecomplex t )
    {
      typedef typename TypeTo::value_type value_type;
      value_type ret_r = Teuchos::as<value_type>( t.r );
      value_type ret_i = Teuchos::as<value_type>( t.i );
      return ( TypeTo( ret_r, ret_i ) );
    }
  };

  //@}  End Conversion group

} // end namespace Teuchos


namespace Amesos2 {

  // forward declaration due to circular reference
  template <class, class> class NewSolver;

  /* Specialize the Amesos2::TypeMap struct for NewSolver types.
   *
   * Additional nested types may be added without harm.  For an example, look at
   * Amesos2_Superlu_TypeMap.hpp
   */
  template <>
  struct TypeMap<NewSolver,float>
  {
    typedef float type;
    typedef float magnitude_type;
  };


  template <>
  struct TypeMap<NewSolver,double>
  {
    typedef double type;
    typedef double magnitude_type;
  };

#ifdef HAVE_TEUCHOS_COMPLEX

  /*
   * We map the std complex types to the appropriate NewSolver complex
   * types.
   */

  template <>
  struct TypeMap<NewSolver,std::complex<float> >
  {
    typedef New_Solver::C::complex type;
    typedef float magnitude_type;
  };


  template <>
  struct TypeMap<NewSolver,std::complex<double> >
  {
    typedef New_Solver::Z::doublecomplex type;
    typedef double magnitude_type;
  };


  template <>
  struct TypeMap<NewSolver,New_Solver::C::complex>
  {
    typedef New_Solver::C::complex type;
    typedef float magnitude_type;
  };


  template <>
  struct TypeMap<NewSolver,New_Solver::Z::doublecomplex>
  {
    typedef New_Solver::Z::doublecomplex type;
    typedef double magnitude_type;
  };
#endif	// HAVE_TEUCHOS_COMPLEX

} // end namespace Amesos2

#endif  // AMESOS2_NEWSOLVER_TYPEMAP_HPP
