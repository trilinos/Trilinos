// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER


#include "Teuchos_OrdinalTraits.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_TypeTraits.hpp"
#include "Teuchos_Version.hpp"

#include "Teuchos_TestingHelpers.hpp"

namespace {

//
// Output of the ordinal that will avoid printing non-asci chars.
//
// The issue is that if you just print a raw char then if it prints non-asci
// character, it can not be printed in some cases and causes trouble with
// CTest/CDash.
//

// For general types, just print the type
template <class T>
T outputOrdinal(const T &t)
{
  return t;
}

// For char, print the int value (avoid non-ansi chars)
int outputOrdinal(const char& t)
{
  return t;
}


//
// Type outputting
//

template <class T>
void TYPE_CHAIN_A(Teuchos::FancyOStream &out) {
  using Teuchos::TypeTraits::is_same;
  typedef typename Teuchos::ScalarTraits<T>::doublePrecision double_type;

  T b;
  // double_type d; // unused
  out << Teuchos::typeName (b);
  if (! is_same<T, double_type>::value) {
    out << " -> ";
    TYPE_CHAIN_A<double_type>(out);
  }
}

template <class T>
void TYPE_CHAIN_D(Teuchos::FancyOStream &out) {
  using Teuchos::TypeTraits::is_same;
  typedef typename Teuchos::ScalarTraits<T>::halfPrecision half_type;

  T b;
  // half_type d; // unused
  out << Teuchos::typeName (b);

  if (! is_same<T, half_type>::value) {
    out << " -> ";
    TYPE_CHAIN_D<half_type>(out);
  }
}


inline std::string passfail(bool result)
{
  return ( result ? "passed" : "failed" );
}


template<class Scalar>
bool testScalarTraits(
  Teuchos::FancyOStream &out
  )
{

  typedef Teuchos::ScalarTraits<Scalar> ST;

  bool success = true;

  out << "\nTesting: " << Teuchos::TypeNameTraits<ST>::name() << " ...\n";

  Teuchos::OSTab tab(out);

  const Scalar nan = ST::nan();

  out << "Type chain (ascending) : "; TYPE_CHAIN_A<Scalar>(out); out << "\n";
  out << "Type chain (descending): "; TYPE_CHAIN_D<Scalar>(out); out << "\n";

  TEUCHOS_TEST_EQUALITY_CONST(ST::isnaninf(1.0), false, out, success);

  out << "\nTesting that squareroot(NaN) == NaN! ...\n";
  {
    const Scalar sqrtNan = ST::squareroot(nan);
    out << "squareroot("<<nan<<") = " << sqrtNan << "\n";
    TEUCHOS_TEST_EQUALITY_CONST(ST::isnaninf(sqrtNan), true, out, success);
  }

  out << "\nTesting that squareroot(-NaN) == NaN! ...\n";
  {
    const Scalar negNan = -nan;
    const Scalar sqrtNegNan = ST::squareroot(negNan);
    out << "squareroot("<<negNan<<") = " << sqrtNegNan << "\n";
    TEUCHOS_TEST_EQUALITY_CONST(ST::isnaninf(sqrtNegNan), true, out, success);
  }

  if (ST::isComplex == false)
  {
    out << "\nTesting that squareroot(-1) == NaN! ...\n";
    {
      const Scalar negOne = -ST::one();
      const Scalar sqrtNegOne = ST::squareroot(negOne);
      out << "squareroot("<<negOne<<") = " << sqrtNegOne << "\n";
      TEUCHOS_TEST_EQUALITY_CONST(ST::isnaninf(sqrtNegOne), true, out, success);
    }
  }

#ifdef HAVE_NUMERIC_LIMITS

  out << "\nTesting that squareroot(quiet_NaN) == NaN! ...\n";
  {
    const Scalar nan = std::numeric_limits<Scalar>::quiet_NaN();
    const Scalar sqrtNan = ST::squareroot(nan);
    TEUCHOS_TEST_EQUALITY_CONST(ST::isnaninf(sqrtNan), true, out, success);
  }

  out << "\nTesting that squareroot(signaling_NaN) == NaN! ...\n";
  {
    const Scalar nan = std::numeric_limits<Scalar>::signaling_NaN();
    const Scalar sqrtNan = ST::squareroot(nan);
    TEUCHOS_TEST_EQUALITY_CONST(ST::isnaninf(sqrtNan), true, out, success);
  }

  out << "\nTesting that squareroot(inf) == NaN! ...\n";
  {
    const Scalar inf = std::numeric_limits<Scalar>::infinity();
    const Scalar sqrtInf = ST::squareroot(inf);
    TEUCHOS_TEST_EQUALITY_CONST(ST::isnaninf(sqrtInf), true, out, success);
  }

#endif // HAVE_NUMERIC_LIMITS

  return success;

}


template<class Ordinal>
bool testOrdinalTraits(
  Teuchos::FancyOStream &out
  )
{

  typedef Teuchos::OrdinalTraits<Ordinal> OT;

  bool success = true;
  bool result;

  out << "\nTesting: " << Teuchos::TypeNameTraits<OT>::name() << " ...\n";

  Teuchos::OSTab tab(out);

  const Ordinal zero = OT::zero();
  const Ordinal one  = OT::one();
  const Ordinal max  = OT::max();
  const Ordinal invalid  = OT::invalid();
  out << "\nmax() == " << outputOrdinal(max) << "\n";
  out << "\ninvalid() == " << outputOrdinal(invalid) << "\n";

  out << "\nTesting that zero() * one() == zero() ...\n";
  {
    const Ordinal zto = zero*one;
    result = (zto == zero);
    if (!result) success = false;
    out
      << "zero*one = " << outputOrdinal(zto) << " == " << outputOrdinal(zero) << " : "
      << passfail(result) << "\n";
  }

  out << "\nTesting that one() * one() == one() ...\n";
  {
    const Ordinal oto = one*one;
    result = (oto == one);
    if (!result) success = false;
    out
      << "one*one = " << outputOrdinal(oto) << " == " << outputOrdinal(one) << " : "
      << passfail(result) << "\n";
  }

  out << "\nTesting that one() + zero() == zero() + one() == one() ...\n";
  {
    const Ordinal opz = one+zero;
    const Ordinal zpo = zero+one;
    result = (opz == one) && (zpo == one);
    if (!result) success = false;
    out
      << "one+zero = " << outputOrdinal(opz) << " == zero+one = "
      << outputOrdinal(zpo) << " == " << outputOrdinal(one) << " : "
      << passfail(result) << "\n";
  }

  out << "\nTesting that one() - one() == zero() ...\n";
  {
    const Ordinal omo = one-one;
    result = (omo == zero);
    if (!result) success = false;
    out
      << "one-one = " << outputOrdinal(omo) << " == " << outputOrdinal(zero) << " : "
      << passfail(result) << "\n";
  }

  out << "\nTesting that zero() < one() <= max() ...\n";
  {
    result = (zero < one) && (one <= max) && (zero < max);
    if (!result) success = false;
    out
      << "(zero < one) = " << (zero < one) << " == "
      << "(one <= max) = " << (one <= max) << " == "
      << "(zero < max) = " << (zero < max) << " == "
      << true << " : "
      << passfail(result) << "\n";
  }

  out << "\nTesting that invalid() not in [zero(),max()]...\n";
  {
    result = !( (invalid > zero || invalid==zero) && (invalid <= max) );
    if (!result) success = false;
    out
      << "invalid in [zero,max] == false : " << passfail(result) << "\n";
  }

  return success;

}


} // namespace


int main( int argc, char* argv[] ) {

  using Teuchos::CommandLineProcessor;

        bool success = true;
  bool result;

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);
  //const int procRank = Teuchos::GlobalMPISession::getRank();

  Teuchos::RCP<Teuchos::FancyOStream>
    out = Teuchos::VerboseObjectBase::getDefaultOStream();

        try {

                // Read options from the commandline
    CommandLineProcessor  clp(false); // Don't throw exceptions
                CommandLineProcessor::EParseCommandLineReturn parse_return = clp.parse(argc,argv);
                if( parse_return != CommandLineProcessor::PARSE_SUCCESSFUL ) {
                        *out << "\nEnd Result: TEST FAILED" << std::endl;
                        return parse_return;
                }

    result = testScalarTraits<float>(*out);
    if(!result) success = false;

    result = testScalarTraits<double>(*out);
    if(!result) success = false;

    result = testOrdinalTraits<char>(*out);
    if(!result) success = false;

    result = testOrdinalTraits<short int>(*out);
    if(!result) success = false;

    result = testOrdinalTraits<int>(*out);
    if(!result) success = false;

    result = testOrdinalTraits<long int>(*out);
    if(!result) success = false;

    result = testOrdinalTraits<size_t>(*out);
    if(!result) success = false;

    result = testOrdinalTraits<long long int>(*out);
    if(!result) success = false;

#ifdef HAVE_TEUCHOS___INT64
    result = testOrdinalTraits<__int64>(*out);
    if(!result) success = false;

    result = testOrdinalTraits<unsigned __int64>(*out);
    if(!result) success = false;
#endif

// #ifdef HAVE_TEUCHOS_COMPLEX
//     result = testScalarTraits<std::complex<double> >(*out);
//     if(!result) success = false;
// #endif

#ifdef HAVE_TEUCHOSCORE_QUADMATH
    result = testScalarTraits<__float128>(*out);
    if(!result) success = false;
#endif // HAVE_TEUCHOSCORE_QUADMATH

#ifdef HAVE_TEUCHOS_QD
    result = testScalarTraits<dd_real>(*out);
    if(!result) success = false;
    result = testScalarTraits<qd_real>(*out);
    if(!result) success = false;
#endif

        }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(true,std::cerr,success);

  if(success)
    *out << "\nEnd Result: TEST PASSED\n" << std::endl;
  else
    *out << "\nEnd Result: TEST FAILED\n" << std::endl;

  return ( success ? 0 : 1 );

}
