// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// \file default_blas_rot.cpp
/// \brief Test Teuchos' generic ROTG implementation.
/// \author Mark Hoemmen
///
/// This program tests Teuchos::DefaultBLASImpl::ROTG against
/// Teuchos::BLAS::ROTG for ScalarType=double.  This is a reasonable
/// correctness test, assuming that Teuchos::BLAS invokes the BLAS
/// library.
///
/// TODO: Test for ScalarType=std::complex<double>.  This would
/// require constructing some complex-valued tests.

#include <Teuchos_as.hpp>
#include <Teuchos_BLAS.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>

// Anonymous namespace to avoid name collisions.
namespace {

  /// \class GivensTester
  /// \brief Test Givens rotations in Teuchos::DefaultBLASImpl.
  /// \author Mark Hoemmen
  template<class OrdinalType, class ScalarType>
  class GivensTester {
  public:
    typedef ScalarType scalar_type;
    typedef typename Teuchos::ScalarTraits<ScalarType>::magnitudeType MagnitudeType;
    typedef MagnitudeType magnitude_type;

  private:
    //! Stream to which to write test results.
    std::ostream& out_;

    /// \brief Error bound for computed cosines and sines.
    ///
    /// We compute this once, in the constructor, to save a little
    /// overhead.
    MagnitudeType trigErrorBound_;

    typedef Teuchos::ScalarTraits<ScalarType> STS;
    typedef Teuchos::ScalarTraits<MagnitudeType> STM;

    //! Compute $\sqrt{|a|^2 + |b|^2}$ without overflow.
    MagnitudeType norm2 (const ScalarType a, const ScalarType& b)
    {
      const MagnitudeType scale = STS::magnitude(a) + STS::magnitude(b);

      if (scale == STM::zero()) {
	return STM::zero();
      } else {
	const ScalarType a_scaled = a / scale;
	const ScalarType b_scaled = b / scale;
	return scale * STM::squareroot (a_scaled*a_scaled + b_scaled*b_scaled);
      }
    }

    /// \brief Maximum allowed absolute difference between two
    ///   different BLAS implementations' computed cosines or sines.
    ///
    /// See Lemma 19.7 in Nicholas J. Higham's "Accuracy and Stability
    /// of Numerical Algorithms," second edition (SIAM).  This method
    /// returns twice that quantity, because we are comparing two
    /// different computations, rather than comparing one computation
    /// with the exact value.  In the worst case, the two computations
    /// might be on opposite sides of the exact value.
    static MagnitudeType trigErrorBound ()
    {
      // NOTE (mfh 12 Oct 2011) I'm not sure if this error bound holds
      // for complex arithmetic.  Does STS report the right machine
      // precision for complex numbers?
      const MagnitudeType u = STS::eps();

      // In Higham's notation, $\gamma_k = ku / (1 - ku)$.
      return 2 * (4*u) / (1 - 4*u);
    }

  public:
    /// \brief Constructor.
    ///
    /// \param out [out] Stream to which to write test results.
    GivensTester (std::ostream& out) :
      out_ (out), trigErrorBound_ (trigErrorBound())
    {}

    //! Test the results of ROTG on [a; b].
    bool compare (ScalarType a, ScalarType b)
    {
      using std::endl;
      typedef Teuchos::DefaultBLASImpl<int, ScalarType> generic_blas_type;
      typedef Teuchos::BLAS<int, ScalarType> library_blas_type;

      out_ << "Comparing Givens rotations for [a; b] = ["
	   << a << "; " << b << "]" << endl;

      generic_blas_type genericBlas;
      library_blas_type libraryBlas;

      // ROTG overwrites its input arguments.
      ScalarType a_generic = a;
      ScalarType b_generic = b;
      MagnitudeType c_generic;
      ScalarType s_generic;
      genericBlas.ROTG (&a_generic, &b_generic, &c_generic, &s_generic);

      ScalarType a_library = a;
      ScalarType b_library = b;
      MagnitudeType c_library;
      ScalarType s_library;
      libraryBlas.ROTG (&a_library, &b_library, &c_library, &s_library);

      out_ << "-- DefaultBLASImpl results: a,b,c,s = "
	   << a_generic << ", " << b_generic << ", "
	   << c_generic << ", " << s_generic << endl;
      out_ << "-- (Library) BLAS results: a,b,c,s = "
	   << a_library << ", " << b_library << ", "
	   << c_library << ", " << s_library << endl;

      bool success = true; // Innocent until proven guilty.

      // Test the difference between the computed cosines.
      out_ << "-- |c_generic - c_library| = "
	   << STS::magnitude(c_generic - c_library) << endl;
      if (STS::magnitude(c_generic - c_library) > trigErrorBound_) {
	success = false;
	out_ << "---- Difference exceeded error bound " << trigErrorBound_ << endl;
      }

      // Test the difference between the computed sines.
      out_ << "-- |s_generic - s_library| = "
	   << STS::magnitude(s_generic - s_library) << endl;
      if (STS::magnitude(s_generic - s_library) > trigErrorBound_) {
	success = false;
	out_ << "---- Difference exceeded error bound " << trigErrorBound_ << endl;
      }

      // Test the forward error of applying the Givens rotation.
      // Remember that ROTG applies the rotation to its input
      // arguments [a; b], overwriting them with the resulting [r; z].
      //
      // See Higham's Lemma 19.8.
      const MagnitudeType inputNorm = norm2 (a, b);
      const MagnitudeType outputDiffNorm =
	norm2 (a_generic - a_library, b_generic - b_library);

      out_ << "-- ||[a; b]||_2 = " << inputNorm << endl;
      out_ << "-- ||[a_generic - a_library; b_generic - b_library]||_2 = "
	   << outputDiffNorm << endl;

      // Multiply by a fudge factor of the base, just in case the
      // forward error bound wasn't computed accurately.  Also
      // multiply by 2, since we don't know the exact result.  The
      // latter is because the two computed results could be on either
      // side of the exact result: sqrt((2 * x_diff)^2 + (2 *
      // y_diff)^2) = sqrt(4) * sqrt(x_diff^2 + y_diff^2).
      const MagnitudeType two = STM::one() + STM::one();
      const MagnitudeType fwdErrorBound =
	2 * STS::base() * STM::squareroot(two) * (6*STS::eps() / (1 - 6*STS::eps()));

      if (outputDiffNorm > fwdErrorBound * inputNorm) {
	success = false;
	out_ << "---- Forward error exceeded relative error bound "
	     << fwdErrorBound << endl;
      }
      return success;
    }

    /// \brief Run \c compare() for several different [a; b] values.
    ///
    /// \warning FIXME This only tests reals at the moment.
    bool test ()
    {
      using Teuchos::as;

      const ScalarType zero = STS::zero();
      const ScalarType one = STS::one();
      //const ScalarType two = one + one;
      //const ScalarType four = two + two;

      bool success = true; // Innocent until proven guilty.

      // First test the corner cases: [\pm 1, 0] and [0, \pm 1].
      success = success && compare (one, zero);
      success = success && compare (zero, one);
      success = success && compare (-one, zero);
      success = success && compare (zero, -one);

      // Test a range of other values.
      {
	const ScalarType incr = one / as<ScalarType> (10);
	for (int k = -30; k < 30; ++k) {
	  const ScalarType a = as<ScalarType> (k) * incr;
	  const ScalarType b = one - as<ScalarType> (k) * incr;
	  success = success && compare (a, b);
	}
      }

      //
      // Try some big values just to see whether ROTG correctly
      // scales its inputs to avoid overflow.
      //
      //success = success && compare (STS::rmax() / four, STS::rmax() / four);
      //success = success && compare (-STS::rmax() / four, STS::rmax() / four);
      //success = success && compare (STS::rmax() / four, -STS::rmax() / four);

      //success = success && compare (STS::rmax() / two, STS::rmax() / two);
      //success = success && compare (-STS::rmax() / two, STS::rmax() / two);
      //success = success && compare (-STS::rmax() / two, -STS::rmax() / two);

      //success = success && compare (STS::rmax() / two, zero);
      //success = success && compare (zero, STS::rmax() / two);
      //success = success && compare (-STS::rmax() / two, zero);
      //success = success && compare (zero, -STS::rmax() / two);

      //success = success && compare (STS::rmax() / two, one);
      //success = success && compare (one, STS::rmax() / two);
      //success = success && compare (-STS::rmax() / two, one);
      //success = success && compare (one, -STS::rmax() / two);

      return success;
    }
  };
} // namespace (anonymous)


int
main (int argc, char *argv[])
{
  using std::endl;
  using Teuchos::CommandLineProcessor;

  Teuchos::GlobalMPISession mpiSession (&argc, &argv, NULL);
  const int myRank = mpiSession.getRank();
  Teuchos::oblackholestream blackHole;
  std::ostream& out = (myRank == 0) ? std::cout : blackHole;

  bool verbose = true;
  CommandLineProcessor cmdp(false,true);
  cmdp.setOption ("verbose", "quiet", &verbose, "Print messages and results.");

  // Parse the command-line arguments.
  {
    const CommandLineProcessor::EParseCommandLineReturn parseResult =
      cmdp.parse (argc,argv);
    // If the caller asks us to print the documentation, let the
    // "test" pass trivially.
    if (parseResult == CommandLineProcessor::PARSE_HELP_PRINTED) {
      out << "End Result: TEST PASSED" << endl;
      return EXIT_SUCCESS;
    }
    TEUCHOS_TEST_FOR_EXCEPTION(parseResult != CommandLineProcessor::PARSE_SUCCESSFUL,
		       std::invalid_argument,
		       "Failed to parse command-line arguments");
  }

  // Only let the tester print if in verbose mode.
  GivensTester<int, double> tester (verbose ? out : blackHole);
  const bool success = tester.test ();

  if (success) {
    out << "End Result: TEST PASSED" << endl;
    return EXIT_SUCCESS;
  } else {
    out << "End Result: TEST FAILED" << endl;
    return EXIT_FAILURE;
  }
}

