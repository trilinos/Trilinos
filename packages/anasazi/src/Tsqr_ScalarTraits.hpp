#ifndef __TSQR_Tsqr_ScalarTraits_hpp
#define __TSQR_Tsqr_ScalarTraits_hpp

#include <cmath> // std::abs
#include <complex>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {

  /// \class ScalarTraits
  ///
  /// \brief A traits class to dispatch to the correct arithmetic properties.
  ///
  /// A traits class to dispatch to the correct arithmetic properties
  /// for a given Scalar data type.
  template< class Scalar >
  class ScalarTraits {
  public:
    /// Whether we've specialized this traits class for the particular Scalar type
    static const bool is_specialized = false;
    /// Whether Scalar represents a complex number
    static const bool is_complex = false;
    typedef Scalar magnitude_type;

    /// The arithmetic identity for the given Scalar data type
    static Scalar zero();
    /// The multiplicative identity for the given Scalar data type
    static Scalar one();
    /// Complex conjugate of z, in case is_complex == true, else just z
    inline static Scalar conj (const Scalar& z);
    /// Absolute value of z
    inline static magnitude_type abs (const Scalar& z);
  };

  template<>
  class ScalarTraits< std::complex< double > > {
  public:
    static const bool is_specialized = true;
    static const bool is_complex = true;
    typedef double magnitude_type;

    static std::complex<double> zero() { return std::complex<double>(0.0, 0.0); }
    static std::complex<double> one()  { return std::complex<double>(1.0, 0.0); }

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

    inline static std::complex<float> conj (const std::complex<float>& z) {
      return std::conj (z);
    }
    inline static magnitude_type abs (const std::complex<float>& z) {
      return std::abs (z);
    }
  };

  template<>
  class ScalarTraits< double > {
  public:
    static const bool is_specialized = true;
    static const bool is_complex = false;
    typedef double magnitude_type;

    static double zero() { return 0.0; }
    static double one() { return 1.0; }

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

    inline static float conj (const float& z) { return z; }
    inline static magnitude_type abs (const float& z) {
      return std::abs (z);
    }
  };

} // namespace TSQR

#endif // __TSQR_Tsqr_ScalarTraits_hpp
