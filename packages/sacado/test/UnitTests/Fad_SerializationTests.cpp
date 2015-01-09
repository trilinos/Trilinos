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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
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
#include "Sacado_Fad_SimpleFad.hpp"
#include "Sacado_CacheFad_DFad.hpp"
#include "Sacado_CacheFad_SFad.hpp"
#include "Sacado_CacheFad_SLFad.hpp"
#include "Sacado_mpl_apply.hpp"
#include "Sacado_Random.hpp"

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::ValueTypeSerializer;

template <typename FadType>
bool testSerialization(const Teuchos::Array<FadType>& x, 
		       const std::string& tag,
		       Teuchos::FancyOStream& out) {
  
  typedef int Ordinal;
  typedef Teuchos::SerializationTraits<Ordinal,FadType> SerT;
  
  // Serialize
  Ordinal count = x.size();
  Ordinal bytes = SerT::fromCountToIndirectBytes(count, &x[0]);
  char *charBuffer = new char[bytes];
  SerT::serialize(count, &x[0], bytes, charBuffer);
  
  // Deserialize
  Ordinal count2 = SerT::fromIndirectBytesToCount(bytes, charBuffer);
  Teuchos::Array<FadType> x2(count2);
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
  
  // Check Fads match
  for (Ordinal i=0; i<count; i++) {
    bool success2 = Sacado::IsEqual<FadType>::eval(x[i], x2[i]);
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

template <typename FadType, typename Serializer>
bool testSerializationObject(const Serializer& serializer,
			     Teuchos::Array<FadType>& x, 
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
    x[i].expand(sz);
  
  // Deserialize
  Ordinal count2 = serializer.fromIndirectBytesToCount(bytes, charBuffer);
  Teuchos::Array<FadType> x2(count2);
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
  
  // Check Fads match
  for (Ordinal i=0; i<count; i++) {
    bool success2 = Sacado::IsEqual<FadType>::eval(x[i], x2[i]);
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

template <typename FadType, typename Serializer>
bool testNestedSerializationObject(const Serializer& serializer,
				   Teuchos::Array<FadType>& x, 
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
    x[i].expand(sz);
    for (Ordinal j=0; j<sz; j++)
      x[i].fastAccessDx(j).expand(sz_inner);
    x[i].val().expand(sz_inner);
  }
  
  // Deserialize
  Ordinal count2 = serializer.fromIndirectBytesToCount(bytes, charBuffer);
  Teuchos::Array<FadType> x2(count2);
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
  
  // Check Fads match
  for (Ordinal i=0; i<count; i++) {
    bool success2 = Sacado::IsEqual<FadType>::eval(x[i], x2[i]);
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

#define FAD_SERIALIZATION_TESTS(FadType, FAD)				\
  TEUCHOS_UNIT_TEST( FAD##_Serialization, FadUniform ) {		\
    int n = 7;								\
    int p = 5;								\
    ValueTypeSerializer<int,FadType> fts(				\
      rcp(new ValueTypeSerializer<int,double>), p);			\
    Teuchos::Array<FadType> x(n);					\
    for (int i=0; i<n; i++) {						\
      x[i] = FadType(p, rnd.number());					\
      for (int j=0; j<p; j++)						\
	x[i].fastAccessDx(j) = rnd.number();				\
    }									\
    bool success1 = testSerialization(					\
      x, std::string(#FAD) + " Uniform", out);				\
    bool success2 = testSerializationObject(				\
      fts, x, std::string(#FAD) + " Uniform FTS", out);			\
    success = success1 && success2;					\
  }									\
									\
  TEUCHOS_UNIT_TEST( FAD##_Serialization, FadEmpty ) {			\
    int n = 7;								\
    ValueTypeSerializer<int,FadType> fts(				\
      rcp(new ValueTypeSerializer<int,double>), 5);			\
    Teuchos::Array<FadType> x(n);					\
    for (int i=0; i<n; i++) {						\
      x[i] = FadType(rnd.number());					\
    }									\
    bool success1 = testSerialization(x, std::string(			\
					#FAD) + " Empty", out);		\
    bool success2 = testSerializationObject(				\
      fts, x, std::string(#FAD) + " Empty FTS", out);			\
    success = success1 && success2;					\
  }									\
									\
  TEUCHOS_UNIT_TEST( FAD##_Serialization, FadMixed ) {			\
    int n = 6;								\
    int p[] = { 5, 0, 8, 8, 3, 0 };					\
    ValueTypeSerializer<int,FadType> fts(				\
      rcp(new ValueTypeSerializer<int,double>), 8);			\
    Teuchos::Array<FadType> x(n);					\
    for (int i=0; i<n; i++) {						\
      x[i] = FadType(p[i], rnd.number());				\
      for (int j=0; j<p[i]; j++)					\
	x[i].fastAccessDx(j) = rnd.number();				\
    }									\
    bool success1 = testSerialization(					\
      x, std::string(#FAD) + " Mixed", out);				\
    bool success2 = testSerializationObject(				\
      fts, x, std::string(#FAD) + " Mixed FTS", out);			\
    success = success1 && success2;					\
  }									\
									\
  TEUCHOS_UNIT_TEST( FAD##_Serialization, FadFadUniform ) {		\
    typedef Sacado::mpl::apply<FadType,FadType>::type FadFadType;	\
    int n = 7;								\
    int p1 = 5;								\
    int p2 = 3;								\
    RCP< ValueTypeSerializer<int,FadType> > fts =			\
      rcp(new ValueTypeSerializer<int,FadType>(				\
	    rcp(new ValueTypeSerializer<int,double>), p1));		\
    ValueTypeSerializer<int,FadFadType> ffts(fts, p2);			\
    Teuchos::Array<FadFadType> x(n);					\
    for (int i=0; i<n; i++) {						\
      FadType f(p1, rnd.number());					\
      for (int k=0; k<p1; k++)						\
	f.fastAccessDx(k) = rnd.number();				\
      x[i] = FadFadType(p2, f);						\
      for (int j=0; j<p2; j++) {					\
	FadType g(p1, rnd.number());					\
	for (int k=0; k<p1; k++)					\
	  g.fastAccessDx(k) = rnd.number();				\
	x[i].fastAccessDx(j) = g;					\
      }									\
    }									\
    bool success1 = testSerialization(					\
      x, std::string(#FAD) + " Nested Uniform", out);			\
    bool success2 =							\
      testNestedSerializationObject(					\
	ffts, x, std::string(#FAD) + " Nested Uniform", out);		\
    success = success1 && success2;					\
  }									\
									\
  TEUCHOS_UNIT_TEST( FAD##_Serialization, FadFadEmptyInner ) {		\
    typedef Sacado::mpl::apply<FadType,FadType>::type FadFadType;	\
    int n = 7;								\
    int p1 = 5;								\
    int p2 = 3;								\
    RCP< ValueTypeSerializer<int,FadType> > fts =			\
      rcp(new ValueTypeSerializer<int,FadType>(				\
	    rcp(new ValueTypeSerializer<int,double>), p1));		\
    ValueTypeSerializer<int,FadFadType> ffts(fts, p2);			\
    Teuchos::Array<FadFadType> x(n);					\
    for (int i=0; i<n; i++) {						\
      FadType f(p1, rnd.number());					\
      for (int k=0; k<p1; k++)						\
	f.fastAccessDx(k) = rnd.number();				\
      x[i] = FadFadType(p2, f);						\
      for (int j=0; j<p2; j++)						\
	x[i].fastAccessDx(j) =  rnd.number();				\
    }									\
    bool success1 = testSerialization(					\
      x, std::string(#FAD) + " Nested Empty Inner", out);		\
    bool success2 = testNestedSerializationObject(			\
      ffts, x, std::string(#FAD) + " Nested Empty Inner FTS", out);	\
    success = success1 && success2;					\
  }									\
									\
  TEUCHOS_UNIT_TEST( FAD##_Serialization, FadFadEmptyOuter ) {		\
    typedef Sacado::mpl::apply<FadType,FadType>::type FadFadType;	\
    int n = 7;								\
    int p1 = 5;								\
    RCP< ValueTypeSerializer<int,FadType> > fts =			\
      rcp(new ValueTypeSerializer<int,FadType>(				\
	    rcp(new ValueTypeSerializer<int,double>), p1));		\
    ValueTypeSerializer<int,FadFadType> ffts(fts, 5);			\
    Teuchos::Array<FadFadType> x(n);					\
    for (int i=0; i<n; i++) {						\
      FadType f(p1, rnd.number());					\
      for (int k=0; k<p1; k++)						\
	f.fastAccessDx(k) = rnd.number();				\
      x[i] = FadFadType(f);						\
    }									\
    bool success1 =testSerialization(					\
	x, std::string(#FAD) + " Nested Empty Outer", out);		\
    bool success2 =testNestedSerializationObject(			\
      ffts, x, std::string(#FAD) + " Nested Empty Outer FTS", out);	\
    success = success1 && success2;					\
  }									\
									\
  TEUCHOS_UNIT_TEST( FAD##_Serialization, FadFadEmptyAll ) {		\
    typedef Sacado::mpl::apply<FadType,FadType>::type FadFadType;	\
    int n = 7;								\
    RCP< ValueTypeSerializer<int,FadType> > fts =			\
      rcp(new ValueTypeSerializer<int,FadType>(				\
	    rcp(new ValueTypeSerializer<int,double>), 5));		\
    ValueTypeSerializer<int,FadFadType> ffts(fts, 5);			\
    Teuchos::Array<FadFadType> x(n);					\
    for (int i=0; i<n; i++) {						\
      x[i] = rnd.number();						\
    }									\
    bool success1 = testSerialization(					\
	x, std::string(#FAD) + " Nested Empty All", out);		\
    bool success2 = testNestedSerializationObject(			\
      ffts, x, std::string(#FAD) + " Nested Empty All FTS", out);	\
    success = success1 && success2;					\
  }									\

#define SFAD_SERIALIZATION_TESTS(FadType, FAD)				\
  TEUCHOS_UNIT_TEST( FAD##_Serialization, FadUniform ) {		\
    int n = 7;								\
    int p = 5;								\
    ValueTypeSerializer<int,FadType> fts(				\
      rcp(new ValueTypeSerializer<int,double>), p);			\
    Teuchos::Array<FadType> x(n);					\
    for (int i=0; i<n; i++) {						\
      x[i] = FadType(p, rnd.number());					\
      for (int j=0; j<p; j++)						\
	x[i].fastAccessDx(j) = rnd.number();				\
    }									\
    bool success1 = testSerialization(					\
      x, std::string(#FAD) + " Uniform", out);				\
    bool success2 = testSerializationObject(				\
      fts, x, std::string(#FAD) + " Uniform FTS", out);			\
    success = success1 && success2;					\
  }									\
									\
  TEUCHOS_UNIT_TEST( FAD##_Serialization, FadFadUniform ) {		\
    typedef Sacado::mpl::apply<FadType,FadType>::type FadFadType;	\
    int n = 7;								\
    int p1 = 5;								\
    int p2 = 5;								\
    RCP< ValueTypeSerializer<int,FadType> > fts =			\
      rcp(new ValueTypeSerializer<int,FadType>(				\
	    rcp(new ValueTypeSerializer<int,double>), p1));		\
    ValueTypeSerializer<int,FadFadType> ffts(fts, p2);			\
    Teuchos::Array<FadFadType> x(n);					\
    for (int i=0; i<n; i++) {						\
      FadType f(p1, rnd.number());					\
      for (int k=0; k<p1; k++)						\
	f.fastAccessDx(k) = rnd.number();				\
      x[i] = FadFadType(p2, f);						\
      for (int j=0; j<p2; j++) {					\
	FadType g(p1, rnd.number());					\
	for (int k=0; k<p1; k++)					\
	  g.fastAccessDx(k) = rnd.number();				\
	x[i].fastAccessDx(j) = g;					\
      }									\
    }									\
    bool success1 = testSerialization(					\
      x, std::string(#FAD) + " Nested Uniform", out);			\
    bool success2 = testNestedSerializationObject(			\
      ffts, x, std::string(#FAD) + " Nested Uniform FTS", out);		\
    success = success1 && success2;					\
  }									\


Sacado::Random<double> rnd;
FAD_SERIALIZATION_TESTS(Sacado::Fad::DFad<double>, Fad_DFad)
FAD_SERIALIZATION_TESTS(Sacado::ELRFad::DFad<double>, ELRFad_DFad)
FAD_SERIALIZATION_TESTS(Sacado::ELRCacheFad::DFad<double>, ELRCacheFad_DFad)
FAD_SERIALIZATION_TESTS(Sacado::CacheFad::DFad<double>, CacheFad_DFad)

typedef Sacado::Fad::SLFad<double,10> Fad_SLFadType;
typedef Sacado::ELRFad::SLFad<double,10> ELRFad_SLFadType;
typedef Sacado::ELRCacheFad::SLFad<double,10> ELRCacheFad_SLFadType;
typedef Sacado::CacheFad::SLFad<double,10> CacheFad_SLFadType;
FAD_SERIALIZATION_TESTS(Fad_SLFadType, Fad_SLFad)
FAD_SERIALIZATION_TESTS(ELRFad_SLFadType, ELRFad_SLFad)
FAD_SERIALIZATION_TESTS(ELRCacheFad_SLFadType, ELRCacheFad_SLFad)
FAD_SERIALIZATION_TESTS(CacheFad_SLFadType, CacheFad_SLFad)

typedef Sacado::Fad::SFad<double,5> Fad_SFadType;
typedef Sacado::ELRFad::SFad<double,5> ELRFad_SFadType;
typedef Sacado::ELRCacheFad::SFad<double,5> ELRCacheFad_SFadType;
typedef Sacado::CacheFad::SFad<double,5> CacheFad_SFadType;
SFAD_SERIALIZATION_TESTS(Fad_SFadType, Fad_SFad)
SFAD_SERIALIZATION_TESTS(ELRFad_SFadType, ELRFad_SFad)
SFAD_SERIALIZATION_TESTS(ELRCacheFad_SFadType, ELRCacheFad_SFad)
SFAD_SERIALIZATION_TESTS(CacheFad_SFadType, CacheFad_SFad)

FAD_SERIALIZATION_TESTS(Sacado::Fad::DMFad<double>, Fad_DMFad)
//typedef Sacado::LFad::LogicalSparse<double,int> Fad_LSType;
//FAD_SERIALIZATION_TESTS(Fad_LSType, LFad_LS)

// DVFad, LFad, Flop

template <>
Sacado::Fad::MemPool* Sacado::Fad::MemPoolStorage<double>::defaultPool_ = NULL;
template <>
Sacado::Fad::MemPool* Sacado::Fad::MemPoolStorage< Sacado::Fad::DMFad<double> >::defaultPool_ = NULL;

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  Sacado::Fad::MemPoolManager<double> poolManager(100);
  Sacado::Fad::MemPool *pool = poolManager.getMemoryPool(10);
  Sacado::Fad::DMFad<double>::setDefaultPool(pool);

  Sacado::Fad::MemPoolManager< Sacado::Fad::DMFad<double> > poolManager2(100);
  Sacado::Fad::MemPool *pool2 = poolManager2.getMemoryPool(5);
  Sacado::Fad::DMFad< Sacado::Fad::DMFad<double> >::setDefaultPool(pool2);

  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
