// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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
// ***********************************************************************
// @HEADER


#include "Teuchos_OrdinalTraits.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#include "Teuchos_Version.hpp"


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
// Type ouptutting
//

template <class T>
void TYPE_CHAIN_A(Teuchos::FancyOStream &out) {
  T b; typename Teuchos::ScalarTraits<T>::doublePrecision d;
  out << Teuchos::typeName(b);
  if (typeid(b) != typeid(d)) {
    out << " -> ";
    TYPE_CHAIN_A<typename Teuchos::ScalarTraits<T>::doublePrecision>(out);
  }
}

template <class T>
void TYPE_CHAIN_D(Teuchos::FancyOStream &out) {
  T b; typename Teuchos::ScalarTraits<T>::halfPrecision d;
  out << Teuchos::typeName(b);
  if (typeid(b) != typeid(d)) {
    out << " -> ";
    TYPE_CHAIN_D<typename Teuchos::ScalarTraits<T>::halfPrecision>(out);
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
  bool result;

  out << "\nTesting: " << Teuchos::TypeNameTraits<ST>::name() << " ...\n";

  Teuchos::OSTab tab(out);

  const Scalar nan = ST::nan();

  out << "Type chain (ascending) : "; TYPE_CHAIN_A<Scalar>(out); out << "\n";
  out << "Type chain (descending): "; TYPE_CHAIN_D<Scalar>(out); out << "\n";

  out << "\nTesting that squareroot(NaN) == NaN! ...\n";
  {
    const Scalar sqrtNan = ST::squareroot(nan);
    result = (ST::isnaninf(sqrtNan));
    if (!result) success = false;
    out
      << "squareroot("<<nan<<") = " << sqrtNan << " == " << nan << " : "
      << passfail(result) << "\n";
  }

  out << "\nTesting that squareroot(-NaN) == NaN! ...\n";
  {
    const Scalar negNan = -nan;
    const Scalar sqrtNegNan = ST::squareroot(negNan);
    result = (ST::isnaninf(sqrtNegNan));
    if (!result) success = false;
    out
      << "squareroot("<<negNan<<") = " << sqrtNegNan << " == " << nan << " : "
      << passfail(result) << "\n";
  }
    
  if (ST::isComplex == false)
  {
    out << "\nTesting that squareroot(-1) == NaN! ...\n";
    {
      const Scalar negOne = -ST::one();
      const Scalar sqrtNegOne = ST::squareroot(negOne);
      result = (ST::isnaninf(sqrtNegOne));
      if (!result) success = false;
      out
        << "squareroot("<<negOne<<") = " << sqrtNegOne << " == " << nan << " : "
        << passfail(result) << "\n";
    }
  }

#ifdef HAVE_NUMERIC_LIMITS
    
  out << "\nTesting that squareroot(quiet_NaN) == NaN! ...\n";
  {
    const Scalar nan = std::numeric_limits<Scalar>::quiet_NaN();
    const Scalar sqrtNan = ST::squareroot(nan);
    result = (ST::isnaninf(sqrtNan));
    if (!result) success = false;
    out
      << "squareroot("<<nan<<") = " << sqrtNan << " == " << nan << " : "
      << passfail(result) << "\n";
  }
    
  out << "\nTesting that squareroot(signaling_NaN) == NaN! ...\n";
  {
    const Scalar nan = std::numeric_limits<Scalar>::signaling_NaN();
    const Scalar sqrtNan = ST::squareroot(nan);
    result = (ST::isnaninf(sqrtNan));
    if (!result) success = false;
    out
      << "squareroot("<<nan<<") = " << sqrtNan << " == " << nan << " : "
      << passfail(result) << "\n";
  }
    
  out << "\nTesting that squareroot(inf) == NaN! ...\n";
  {
    const Scalar inf = std::numeric_limits<Scalar>::infinity();
    const Scalar sqrtInf = ST::squareroot(inf);
    result = (ST::isnaninf(sqrtInf));
    if (!result) success = false;
    out
      << "squareroot("<<inf<<") = " << sqrtInf << " == " << nan << " : "
      << passfail(result) << "\n";
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

#ifdef HAVE_TEUCHOS_LONG_LONG_INT
    result = testOrdinalTraits<long long int>(*out);
    if(!result) success = false;
#endif


// #ifdef HAVE_TEUCHOS_COMPLEX
//     result = testScalarTraits<std::complex<double> >(*out);
//     if(!result) success = false;
// #endif
    
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
