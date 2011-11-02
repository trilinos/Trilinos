// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_TestingHelpers.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Teuchos_Array.hpp"
#include "Sacado.hpp"
#include "Sacado_Tay_CacheTaylor.hpp"
#include "Sacado_mpl_apply.hpp"
#include "Sacado_Random.hpp"

template <typename TayType>
bool testSerialization(const Teuchos::Array<TayType>& x, 
		       const std::string& tag,
		       Teuchos::FancyOStream& out) {
  
  typedef int Ordinal;
  typedef Teuchos::SerializationTraits<Ordinal,TayType> SerT;
  
  // Serialize
  Ordinal count = x.size();
  Ordinal bytes = SerT::fromCountToIndirectBytes(count, &x[0]);
  char *charBuffer = new char[bytes];
  SerT::serialize(count, &x[0], bytes, charBuffer);
  
  // Deserialize
  Ordinal count2 = SerT::fromIndirectBytesToCount(bytes, charBuffer);
  Teuchos::Array<TayType> x2(count2);
  SerT::deserialize(bytes, charBuffer, count2, &x2[0]);

  delete [] charBuffer;
  
  // Check counts match
  bool success = (count == count2);
  out << tag << " serialize/deserialize count test";
  if (success)
    out << " passed";
  else
    out << " failed";
  out << ":  \n\tExpected:  " << count << ", \n\tGot:       " << count2 << "." 
      << std::endl;
  
  // Check coefficients match
  for (Ordinal i=0; i<count; i++) {
    bool success2 = Sacado::IsEqual<TayType>::eval(x[i], x2[i]);
    out << tag << " serialize/deserialize fad test " << i;
    if (success2)
      out << " passed";
    else
	out << " failed";
    out << ":  \n\tExpected:  " << x[i] << ", \n\tGot:       " << x2[i] 
	<< "." << std::endl;
    success = success && success2;
  }

  return success;
}

#define TAY_SERIALIZATION_TESTS(TayType, TAY)				\
  TEUCHOS_UNIT_TEST( TAY##_Serialization, Uniform ) {			\
    int n = 7;								\
    int p = 5;								\
    Teuchos::Array<TayType> x(n);					\
    for (int i=0; i<n; i++) {						\
      x[i] = TayType(p, rnd.number());					\
      for (int j=0; j<=p; j++)						\
	x[i].fastAccessCoeff(j) = rnd.number();				\
    }									\
    success = testSerialization(x, std::string(#TAY) + " Uniform", out); \
  }									\
									\
  TEUCHOS_UNIT_TEST( TAY##_Serialization, Empty ) {			\
    int n = 7;								\
    Teuchos::Array<TayType> x(n);					\
    for (int i=0; i<n; i++) {						\
      x[i] = TayType(rnd.number());					\
    }									\
    success = testSerialization(x, std::string(#TAY) + " Empty", out);	\
  }									\
									\
  TEUCHOS_UNIT_TEST( TAY##_Serialization, Mixed ) {			\
    int n = 6;								\
    int p[] = { 5, 0, 8, 8, 3, 0 };					\
    Teuchos::Array<TayType> x(n);					\
    for (int i=0; i<n; i++) {						\
      x[i] = TayType(p[i], rnd.number());				\
      for (int j=0; j<=p[i]; j++)					\
	x[i].fastAccessCoeff(j) = rnd.number();				\
    }									\
    success = testSerialization(x, std::string(#TAY) + " Mixed", out);	\
  }									\
									\
  TEUCHOS_UNIT_TEST( TAY##_Serialization, NestedUniform ) {		\
    typedef Sacado::mpl::apply<TayType,TayType>::type TayTayType;	\
    int n = 7;								\
    int p1 = 5;								\
    int p2 = 3;								\
    Teuchos::Array<TayTayType> x(n);					\
    for (int i=0; i<n; i++) {						\
      TayType f(p1, rnd.number());					\
      for (int k=0; k<=p1; k++)						\
	f.fastAccessCoeff(k) = rnd.number();				\
      x[i] = TayTayType(p2, f);						\
      for (int j=0; j<=p2; j++) {					\
	TayType g(p1, rnd.number());					\
	for (int k=0; k<=p1; k++)					\
	  g.fastAccessCoeff(k) = rnd.number();				\
	x[i].fastAccessCoeff(j) = g;					\
      }									\
    }									\
    success =								\
      testSerialization(x, std::string(#TAY) + " Nested Uniform", out); \
  }									\
									\
  TEUCHOS_UNIT_TEST( TAY##_Serialization, NestedEmptyInner ) {		\
    typedef Sacado::mpl::apply<TayType,TayType>::type TayTayType;	\
    int n = 7;								\
    int p1 = 5;								\
    int p2 = 3;								\
    Teuchos::Array<TayTayType> x(n);					\
    for (int i=0; i<n; i++) {						\
      TayType f(p1, rnd.number());					\
      for (int k=0; k<=p1; k++)						\
	f.fastAccessCoeff(k) = rnd.number();				\
      x[i] = TayTayType(p2, f);						\
      for (int j=0; j<=p2; j++)						\
	x[i].fastAccessCoeff(j) =  rnd.number();			\
    }									\
    success =								\
      testSerialization(						\
	x, std::string(#TAY) + " Nested Empty Inner", out);		\
  }									\
									\
  TEUCHOS_UNIT_TEST( TAY##_Serialization, NestedEmptyOuter ) {		\
    typedef Sacado::mpl::apply<TayType,TayType>::type TayTayType;	\
    int n = 7;								\
    int p1 = 5;								\
    Teuchos::Array<TayTayType> x(n);					\
    for (int i=0; i<n; i++) {						\
      TayType f(p1, rnd.number());					\
      for (int k=0; k<=p1; k++)						\
	f.fastAccessCoeff(k) = rnd.number();				\
      x[i] = TayTayType(f);						\
    }									\
    success =								\
      testSerialization(						\
	x, std::string(#TAY) + " Nested Empty Outer", out);		\
  }									\
									\
  TEUCHOS_UNIT_TEST( TAY##_Serialization, NestedEmptyAll ) {		\
    typedef Sacado::mpl::apply<TayType,TayType>::type TayTayType;	\
    int n = 7;								\
    Teuchos::Array<TayTayType> x(n);					\
    for (int i=0; i<n; i++) {						\
      x[i] = rnd.number();						\
    }									\
    success =								\
      testSerialization(						\
	x, std::string(#TAY) + " Nested Empty All", out);		\
  }									\

Sacado::Random<double> rnd;
TAY_SERIALIZATION_TESTS(Sacado::Tay::Taylor<double>, Taylor)
TAY_SERIALIZATION_TESTS(Sacado::Tay::CacheTaylor<double>, CacheTaylor)

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
