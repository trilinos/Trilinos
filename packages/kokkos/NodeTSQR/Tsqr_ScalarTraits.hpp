//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2009) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#ifndef __TSQR_Tsqr_ScalarTraits_hpp
#define __TSQR_Tsqr_ScalarTraits_hpp

#include <Tsqr_ConfigDefs.hpp>

#include <cmath> // std::abs
#include <complex>


namespace TSQR {

  /// \class ScalarTraits
  ///
  /// \brief Map from Scalar type to its arithmetic properties
  ///
  /// ScalarTraits dispatches from a Scalar data type, to its
  /// arithmetic properties.  These include the type of its absolute
  /// value / magnitude, zero, one, \f$\pi\f$, and functions for
  /// computing its absolute value resp. conjugate.
  ///
  /// \note Models for Scalar: float, double, std::complex<float>,
  ///   std::complex<double> (which correspond to the four data types
  ///   S,D,C,Z supported by the BLAS and LAPACK).  If ScalarTraits<
  ///   Scalar >::is_complex, then Scalar should have at least two
  ///   different constructors: a default constructor, and a
  ///   constructor of the form "Scalar (const magnitude_type&, const
  ///   magnitude_type& = magnitude_type(0))".  Otherwise, Scalar
  ///   should have a default constructor, and a "Scalar (const
  ///   magnitude_type&)" constructor.  magnitude_type should follow
  ///   the latter model.
  template< class Scalar >
  class ScalarTraits {
  public:
    /// \brief Was this traits class specialized for Scalar?
    ///
    /// Whether we've specialized this traits class for the particular
    /// Scalar type.  If you're writing your own specialization, you
    /// should set this to true.
    static const bool is_specialized = false;

    //! Whether Scalar represents a complex number.
    static const bool is_complex = false;
    /// If Scalar is complex, this is the type of its magnitude, and
    /// the type of its real and complex parts.  (That means, if the
    /// real and complex parts can be negative, magnitude_type is
    /// allowed to be negative as well, even though magnitudes
    /// themselves (as returned by abs(), see below) are nonnegative.)
    typedef Scalar magnitude_type;
    ///
    /// The arithmetic identity for the given Scalar data type.
    static Scalar zero();
    ///
    /// The multiplicative identity for the given Scalar data type.
    static Scalar one();

    /// The value of \f$\pi\f$ (ratio of a circle's circumference to
    /// its diameter) for magnitude_type.  
    static magnitude_type pi();

    //! Complex conjugate of z, in case is_complex == true, else just z
    inline static Scalar conj (const Scalar& z);

    //! Absolute value of z
    inline static magnitude_type abs (const Scalar& z);
  };

  template<>
  class ScalarTraits< double > {
  public:
    static const bool is_specialized = true;
    static const bool is_complex = false;
    typedef double magnitude_type;

    static double zero() { return 0.0; }
    static double one() { return 1.0; }
    static magnitude_type pi() {
      // In double precision, 17 digits suffice.  We include 20 just
      // for good measure.  Hopefully the C++ compiler won't do a
      // stupid thing with them.
      return 3.14159265358979323846;
    }

    inline static double conj (const double& z) { return z; }
    inline static magnitude_type abs (const double& z) {
      return std::abs (z);
    }
  };

  template<>
  class ScalarTraits< float > {
  public:
    static const bool is_specialized = true;
    static const bool is_complex = false;
    typedef float magnitude_type;

    static float zero() { return 0.0; }
    static float one() { return 1.0; }
    static magnitude_type pi() {
      return 3.14159265358979323846;
    }

    inline static float conj (const float& z) { return z; }
    inline static magnitude_type abs (const float& z) {
      return std::abs (z);
    }
  };

  template<>
  class ScalarTraits< std::complex< double > > {
  public:
    static const bool is_specialized = true;
    static const bool is_complex = true;
    typedef double magnitude_type;

    static std::complex<double> zero() { return std::complex<double>(0.0, 0.0); }
    static std::complex<double> one()  { return std::complex<double>(1.0, 0.0); }
    static magnitude_type pi() { return ScalarTraits< magnitude_type >::pi(); }

    inline static std::complex<double> conj (const std::complex<double>& z) {
      return std::conj (z);
    }
    inline static magnitude_type abs (const std::complex<double>& z) {
      return std::abs (z);
    }
  };

  template<>
  class ScalarTraits< std::complex< float > > {
  public:
    static const bool is_specialized = true;
    static const bool is_complex = true;
    typedef float magnitude_type;

    static std::complex<float> zero() { return std::complex<float>(0.0, 0.0); }
    static std::complex<float> one()  { return std::complex<float>(1.0, 0.0); }
    static magnitude_type pi() { return ScalarTraits< magnitude_type >::pi(); }

    inline static std::complex<float> conj (const std::complex<float>& z) {
      return std::conj (z);
    }
    inline static magnitude_type abs (const std::complex<float>& z) {
      return std::abs (z);
    }
  };


} // namespace TSQR

#endif // __TSQR_Tsqr_ScalarTraits_hpp
