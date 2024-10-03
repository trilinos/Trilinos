// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_TestingHelpers.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Teuchos_Array.hpp"
#include "Sacado_No_Kokkos.hpp"
#include "Sacado_Tay_CacheTaylor.hpp"
#include "Sacado_mpl_apply.hpp"
#include "Sacado_Random.hpp"

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::ValueTypeSerializer;

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
    out << tag << " serialize/deserialize taylor test " << i;
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

template <typename TayType, typename Serializer>
bool testSerializationObject(const Serializer& serializer,
			     Teuchos::Array<TayType>& x, 
			     const std::string& tag,
			     Teuchos::FancyOStream& out) {
  
  typedef int Ordinal;
  
  // Serialize
  Ordinal count = x.size();
  Ordinal bytes = serializer.fromCountToIndirectBytes(count, &x[0]);
  char *charBuffer = new char[bytes];
  serializer.serialize(count, &x[0], bytes, charBuffer);

  // Expand x to serializer size
  Ordinal sz = serializer.getSerializerSize();
  for (Ordinal i=0; i<count; i++)
    x[i].resize(sz-1, true);
  
  // Deserialize
  Ordinal count2 = serializer.fromIndirectBytesToCount(bytes, charBuffer);
  Teuchos::Array<TayType> x2(count2);
  serializer.deserialize(bytes, charBuffer, count2, &x2[0]);

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
    out << tag << " serialize/deserialize taylor test " << i;
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

template <typename TayType, typename Serializer>
bool testNestedSerializationObject(const Serializer& serializer,
				   Teuchos::Array<TayType>& x, 
				   const std::string& tag,
				   Teuchos::FancyOStream& out) {
  
  typedef int Ordinal;
  
  // Serialize
  Ordinal count = x.size();
  Ordinal bytes = serializer.fromCountToIndirectBytes(count, &x[0]);
  char *charBuffer = new char[bytes];
  serializer.serialize(count, &x[0], bytes, charBuffer);

  // Expand x to serializer size
  Ordinal sz = serializer.getSerializerSize();
  typedef typename Serializer::value_serializer_type VST;
  RCP<const VST> vs = serializer.getValueSerializer();
  Ordinal sz_inner = vs->getSerializerSize();
  for (Ordinal i=0; i<count; i++) {
    x[i].resize(sz-1, true);
    for (Ordinal j=0; j<sz; j++)
      x[i].fastAccessCoeff(j).resize(sz_inner-1, true);
    x[i].val().resize(sz_inner-1, true);
  }
  
  // Deserialize
  Ordinal count2 = serializer.fromIndirectBytesToCount(bytes, charBuffer);
  Teuchos::Array<TayType> x2(count2);
  serializer.deserialize(bytes, charBuffer, count2, &x2[0]);

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
    out << tag << " serialize/deserialize taylor test " << i;
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
    ValueTypeSerializer<int,TayType> tts(				\
      rcp(new ValueTypeSerializer<int,double>), p+1);			\
    Teuchos::Array<TayType> x(n);					\
    for (int i=0; i<n; i++) {						\
      x[i] = TayType(p, rnd.number());					\
      for (int j=0; j<=p; j++)						\
	x[i].fastAccessCoeff(j) = rnd.number();				\
    }									\
    bool success1 = testSerialization(					\
      x, std::string(#TAY) + " Uniform", out);				\
    bool success2 = testSerializationObject(				\
      tts, x, std::string(#TAY) + " Uniform TTS", out);			\
    success = success1 && success2;					\
  }									\
									\
  TEUCHOS_UNIT_TEST( TAY##_Serialization, Empty ) {			\
    int n = 7;								\
    ValueTypeSerializer<int,TayType> tts(				\
      rcp(new ValueTypeSerializer<int,double>), 7);			\
    Teuchos::Array<TayType> x(n);					\
    for (int i=0; i<n; i++) {						\
      x[i] = TayType(rnd.number());					\
    }									\
    bool success1 = testSerialization(					\
      x, std::string(#TAY) + " Empty", out);				\
    bool success2 = testSerializationObject(				\
      tts, x, std::string(#TAY) + " Empty TTS", out);			\
    success = success1 && success2;					\
  }									\
									\
  TEUCHOS_UNIT_TEST( TAY##_Serialization, Mixed ) {			\
    int n = 6;								\
    int p[] = { 5, 0, 8, 8, 3, 0 };					\
    ValueTypeSerializer<int,TayType> tts(				\
      rcp(new ValueTypeSerializer<int,double>), 9);			\
    Teuchos::Array<TayType> x(n);					\
    for (int i=0; i<n; i++) {						\
      x[i] = TayType(p[i], rnd.number());				\
      for (int j=0; j<=p[i]; j++)					\
	x[i].fastAccessCoeff(j) = rnd.number();				\
    }									\
    bool success1 = testSerialization(					\
      x, std::string(#TAY) + " Mixed", out);				\
    bool success2 = testSerializationObject(				\
      tts, x, std::string(#TAY) + " Mixed TTS", out);			\
    success = success1 && success2;					\
  }									\
									\
  TEUCHOS_UNIT_TEST( TAY##_Serialization, NestedUniform ) {		\
    typedef Sacado::mpl::apply<TayType,TayType>::type TayTayType;	\
    int n = 7;								\
    int p1 = 5;								\
    int p2 = 3;								\
    RCP< ValueTypeSerializer<int,TayType> > tts =			\
      rcp(new ValueTypeSerializer<int,TayType>(				\
	    rcp(new ValueTypeSerializer<int,double>), p1+1));		\
     ValueTypeSerializer<int,TayTayType> ttts(tts, p2+1);		\
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
    bool success1 = testSerialization(					\
      x, std::string(#TAY) + " Nested Uniform", out);			\
    bool success2 = testNestedSerializationObject(			\
      ttts, x, std::string(#TAY) + " Nested Uniform TTS", out);		\
    success = success1 && success2;					\
  }									\
									\
  TEUCHOS_UNIT_TEST( TAY##_Serialization, NestedEmptyInner ) {		\
    typedef Sacado::mpl::apply<TayType,TayType>::type TayTayType;	\
    int n = 7;								\
    int p1 = 5;								\
    int p2 = 3;								\
    RCP< ValueTypeSerializer<int,TayType> > tts =			\
      rcp(new ValueTypeSerializer<int,TayType>(				\
	    rcp(new ValueTypeSerializer<int,double>), p1+1));		\
     ValueTypeSerializer<int,TayTayType> ttts(tts, p2+1);		\
    Teuchos::Array<TayTayType> x(n);					\
    for (int i=0; i<n; i++) {						\
      TayType f(p1, rnd.number());					\
      for (int k=0; k<=p1; k++)						\
	f.fastAccessCoeff(k) = rnd.number();				\
      x[i] = TayTayType(p2, f);						\
      for (int j=0; j<=p2; j++)						\
	x[i].fastAccessCoeff(j) =  rnd.number();			\
    }									\
    bool success1 = testSerialization(					\
      x, std::string(#TAY) + " Nested Empty Inner", out);		\
    bool success2 = testNestedSerializationObject(			\
      ttts, x, std::string(#TAY) + " Nested Empty Inner TTS", out);	\
    success = success1 && success2;					\
  }									\
									\
  TEUCHOS_UNIT_TEST( TAY##_Serialization, NestedEmptyOuter ) {		\
    typedef Sacado::mpl::apply<TayType,TayType>::type TayTayType;	\
    int n = 7;								\
    int p1 = 5;								\
    RCP< ValueTypeSerializer<int,TayType> > tts =			\
      rcp(new ValueTypeSerializer<int,TayType>(				\
	    rcp(new ValueTypeSerializer<int,double>), p1+1));		\
    ValueTypeSerializer<int,TayTayType> ttts(tts, 5);			\
    Teuchos::Array<TayTayType> x(n);					\
    for (int i=0; i<n; i++) {						\
      TayType f(p1, rnd.number());					\
      for (int k=0; k<=p1; k++)						\
	f.fastAccessCoeff(k) = rnd.number();				\
      x[i] = TayTayType(f);						\
    }									\
    bool success1 = testSerialization(					\
      x, std::string(#TAY) + " Nested Empty Outer", out);		\
    bool success2 = testNestedSerializationObject(			\
      ttts, x, std::string(#TAY) + " Nested Empty Outer TTS", out);	\
    success = success1 && success2;					\
  }									\
									\
  TEUCHOS_UNIT_TEST( TAY##_Serialization, NestedEmptyAll ) {		\
    typedef Sacado::mpl::apply<TayType,TayType>::type TayTayType;	\
    int n = 7;								\
    RCP< ValueTypeSerializer<int,TayType> > tts =			\
      rcp(new ValueTypeSerializer<int,TayType>(				\
	    rcp(new ValueTypeSerializer<int,double>), 5));		\
    ValueTypeSerializer<int,TayTayType> ttts(tts, 5);			\
    Teuchos::Array<TayTayType> x(n);					\
    for (int i=0; i<n; i++) {						\
      x[i] = rnd.number();						\
    }									\
    bool success1 = testSerialization(					\
      x, std::string(#TAY) + " Nested Empty All", out);			\
    bool success2 = testNestedSerializationObject(			\
      ttts, x, std::string(#TAY) + " Nested Empty All TTS", out);	\
    success = success1 && success2;					\
  }									\

Sacado::Random<double> rnd;
TAY_SERIALIZATION_TESTS(Sacado::Tay::Taylor<double>, Taylor)
TAY_SERIALIZATION_TESTS(Sacado::Tay::CacheTaylor<double>, CacheTaylor)

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
