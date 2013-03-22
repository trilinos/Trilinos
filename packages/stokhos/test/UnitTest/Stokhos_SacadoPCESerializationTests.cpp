// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_TestingHelpers.hpp"
#include "Teuchos_UnitTestRepository.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "Teuchos_Array.hpp"
#include "Stokhos_Sacado.hpp"
#include "Sacado_Fad_DFad.hpp"
#include "Sacado_mpl_apply.hpp"
#include "Sacado_Random.hpp"

using Teuchos::RCP;
using Teuchos::rcp;

// Common setup for unit tests
template <typename PCEType, typename FadType>
struct UnitTestSetup {
  RCP<const Stokhos::CompletePolynomialBasis<int,double> > basis;
  RCP<Stokhos::Sparse3Tensor<int,double> > Cijk;
  RCP< Stokhos::AlgebraicOrthogPolyExpansion<int,double> > exp;

  typedef Teuchos::ValueTypeSerializer<int, PCEType> PCESerializerT;
  RCP<PCESerializerT> pce_serializer;

  typedef typename Sacado::mpl::apply<FadType,PCEType>::type FadPCEType;
  typedef Teuchos::ValueTypeSerializer<int, FadPCEType> FadPCESerializerT;
  RCP<FadPCESerializerT> fad_pce_serializer;
  int sz;

  UnitTestSetup() {
    const int d = 2;
    const int p = 7;
      
    // Create product basis
    Teuchos::Array< RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > bases(d);
    for (int i=0; i<d; i++)
      bases[i] = 
	rcp(new Stokhos::LegendreBasis<int,double>(p));
    basis =
      rcp(new Stokhos::CompletePolynomialBasis<int,double>(bases));
    
    // Triple product tensor
    Cijk = basis->computeTripleProductTensor();
      
    // Expansion
    exp = 
      rcp(new Stokhos::AlgebraicOrthogPolyExpansion<int,double>(basis, Cijk));

    // Serializers
    pce_serializer = 
      rcp(new PCESerializerT(
	    exp,
	    rcp(new Teuchos::ValueTypeSerializer<int,double>())));
    fad_pce_serializer = rcp(new FadPCESerializerT(pce_serializer, 8));

    sz = basis->size();
  }
};

template <typename PCEType>
bool testSerialization(const Teuchos::Array<PCEType>& x, 
		       const std::string& tag,
		       Teuchos::FancyOStream& out) {
  
  typedef int Ordinal;
  typedef Teuchos::SerializationTraits<Ordinal,PCEType> SerT;
  
  // Serialize
  Ordinal count = x.size();
  Ordinal bytes = SerT::fromCountToIndirectBytes(count, &x[0]);
  char *charBuffer = new char[bytes];
  SerT::serialize(count, &x[0], bytes, charBuffer);
  Ordinal count2 = SerT::fromIndirectBytesToCount(bytes, charBuffer);

  // Check counts match
  bool success = (count == count2);
  out << tag << " serialize/deserialize count test";
  if (success)
    out << " passed";
  else
    out << " failed";
  out << ":  \n\tExpected:  " << count << ", \n\tGot:       " << count2 << "." 
      << std::endl;

  // Deserialize
  Teuchos::Array<PCEType> x2(count2);
  for (Ordinal i=0; i<count2; i++)
    x2[i].reset(x[i].expansion());
  SerT::deserialize(bytes, charBuffer, count2, &x2[0]);

  delete [] charBuffer;
  
  // Check coefficients match
  for (Ordinal i=0; i<count; i++) {
    bool success2 = Sacado::IsEqual<PCEType>::eval(x[i], x2[i]);
    out << tag << " serialize/deserialize pce test " << i;
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

template <typename PCEType, typename Serializer>
bool testSerialization(Teuchos::Array<PCEType>& x, 
		       const Serializer& serializer,
		       const std::string& tag,
		       Teuchos::FancyOStream& out) {
  
  typedef int Ordinal;
  
  // Serialize
  Ordinal count = x.size();
  Ordinal bytes = serializer.fromCountToIndirectBytes(count, &x[0]);
  char *charBuffer = new char[bytes];
  serializer.serialize(count, &x[0], bytes, charBuffer);

  // Reset x to given expansion
  for (Ordinal i=0; i<count; i++)
    x[i].reset(serializer.getSerializerExpansion());
  
  // Deserialize
  Ordinal count2 = serializer.fromIndirectBytesToCount(bytes, charBuffer);
  Teuchos::Array<PCEType> x2(count2);
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
    bool success2 = Sacado::IsEqual<PCEType>::eval(x[i], x2[i]);
    out << tag << " serialize/deserialize pce test " << i;
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

template <typename PCEType, typename Serializer>
bool testNestedSerialization(Teuchos::Array<PCEType>& x, 
			     const Serializer& serializer,
			     const std::string& tag,
			     Teuchos::FancyOStream& out) {
  
  typedef int Ordinal;
  
  // Serialize
  Ordinal count = x.size();
  Ordinal bytes = serializer.fromCountToIndirectBytes(count, &x[0]);
  char *charBuffer = new char[bytes];
  serializer.serialize(count, &x[0], bytes, charBuffer);

  // Reset x to given expansion
  Ordinal sz = serializer.getSerializerSize();
  typedef typename Serializer::value_serializer_type VST;
  RCP<const VST> vs = serializer.getValueSerializer();
  for (Ordinal i=0; i<count; i++) {
    x[i].expand(sz);
    for (Ordinal j=0; j<sz; j++)
      x[i].fastAccessDx(j).reset(vs->getSerializerExpansion());
    x[i].val().reset(vs->getSerializerExpansion());
  }
  
  // Deserialize
  Ordinal count2 = serializer.fromIndirectBytesToCount(bytes, charBuffer);
  Teuchos::Array<PCEType> x2(count2);
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
    bool success2 = Sacado::IsEqual<PCEType>::eval(x[i], x2[i]);
    out << tag << " serialize/deserialize pce test " << i;
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

#define PCE_SERIALIZATION_TESTS(PCEType, FadType, PCE)			\
  TEUCHOS_UNIT_TEST( PCE##_Serialization, Uniform ) {			\
    int n = 7;								\
    Teuchos::Array<PCEType> x(n);					\
    for (int i=0; i<n; i++) {						\
      x[i] = PCEType(setup.exp);					\
      for (int j=0; j<setup.sz; j++)					\
	x[i].fastAccessCoeff(j) = rnd.number();				\
    }									\
    bool success1 = testSerialization(					\
      x, std::string(#PCE) + " Uniform", out);				\
    bool success2 = testSerialization(					\
      x, *setup.pce_serializer, std::string(#PCE) + " Uniform PTS", out); \
    success = success1 && success2;					\
  }									\
									\
  TEUCHOS_UNIT_TEST( PCE##_Serialization, Empty ) {			\
    int n = 7;								\
    Teuchos::Array<PCEType> x(n);					\
    for (int i=0; i<n; i++) {						\
      x[i] = rnd.number();						\
    }									\
    bool success1 = testSerialization(					\
      x, std::string(#PCE) + " Empty", out);				\
    bool success2 = testSerialization(					\
      x, *setup.pce_serializer, std::string(#PCE) + " Empty PTS", out);	\
    success = success1 && success2;					\
  }									\
									\
  TEUCHOS_UNIT_TEST( PCE##_Serialization, Mixed ) {			\
    int n = 6;								\
    int p[] = { 5, 0, 8, 8, 3, 0 };					\
    Teuchos::Array<PCEType> x(n);					\
    for (int i=0; i<n; i++) {						\
      x[i] = PCEType(setup.exp, p[i]);					\
      for (int j=0; j<p[i]; j++)					\
	x[i].fastAccessCoeff(j) = rnd.number();				\
    }									\
    bool success1 = testSerialization(					\
      x, std::string(#PCE) + " Mixed", out);				\
    bool success2 = testSerialization(					\
      x, *setup.pce_serializer, std::string(#PCE) + " Mixed PTS", out);	\
    success = success1 && success2;					\
  }									\
									\
  TEUCHOS_UNIT_TEST( PCE##_Serialization, FadPCEUniform ) {		\
    typedef Sacado::mpl::apply<FadType,PCEType>::type FadPCEType;	\
    int n = 7;								\
    int p = 3;								\
    Teuchos::Array<FadPCEType> x(n);					\
    for (int i=0; i<n; i++) {						\
      PCEType f(setup.exp);						\
      for (int k=0; k<setup.sz; k++)					\
	f.fastAccessCoeff(k) = rnd.number();				\
      x[i] = FadPCEType(p, f);						\
      for (int j=0; j<p; j++) {						\
	PCEType g(setup.exp);						\
	for (int k=0; k<setup.sz; k++)					\
	  g.fastAccessCoeff(k) = rnd.number();				\
	x[i].fastAccessDx(j) = g;					\
      }									\
    }									\
    success =								\
      testNestedSerialization(x, *setup.fad_pce_serializer,		\
			std::string(#PCE) + " Nested Uniform", out);	\
  }									\
  TEUCHOS_UNIT_TEST( PCE##_Serialization, FadPCEEmptyInner ) {		\
    typedef Sacado::mpl::apply<FadType,PCEType>::type FadPCEType;	\
    int n = 7;								\
    int p = 3;								\
    Teuchos::Array<FadPCEType> x(n);					\
    for (int i=0; i<n; i++) {						\
      PCEType f(setup.exp);						\
      for (int k=0; k<setup.sz; k++)					\
	f.fastAccessCoeff(k) = rnd.number();				\
      x[i] = FadPCEType(p, f);						\
      for (int j=0; j<p; j++)						\
	x[i].fastAccessDx(j) =  rnd.number();				\
    }									\
    success =								\
      testNestedSerialization(						\
	x, *setup.fad_pce_serializer,					\
	std::string(#PCE) + " Nested Empty Inner", out);		\
  }									\
  TEUCHOS_UNIT_TEST( PCE##_Serialization, FadPCEEmptyOuter ) {		\
    typedef Sacado::mpl::apply<FadType,PCEType>::type FadPCEType;	\
    int n = 7;								\
    Teuchos::Array<FadPCEType> x(n);					\
    for (int i=0; i<n; i++) {						\
      PCEType f(setup.exp);						\
      for (int k=0; k<setup.sz; k++)					\
	f.fastAccessCoeff(k) = rnd.number();				\
      x[i] = FadPCEType(f);						\
    }									\
    success =								\
      testNestedSerialization(						\
	x, *setup.fad_pce_serializer,					\
	std::string(#PCE) + " Nested Empty Outer", out);		\
  }									\
  TEUCHOS_UNIT_TEST( PCE##_Serialization, FadPCEEmptyAll ) {		\
    typedef Sacado::mpl::apply<FadType,PCEType>::type FadPCEType;	\
    int n = 7;								\
    Teuchos::Array<FadPCEType> x(n);					\
    for (int i=0; i<n; i++) {						\
      x[i] = rnd.number();						\
    }									\
    success =								\
      testNestedSerialization(						\
	x, *setup.fad_pce_serializer,					\
	std::string(#PCE) + " Nested Empty All", out);			\
  }

typedef Stokhos::StandardStorage<int,double> storage_type;
typedef Sacado::Fad::DFad<double> fad_type;
namespace PCETest {
  Sacado::Random<double> rnd;
  typedef Sacado::PCE::OrthogPoly<double,storage_type> pce_type;
  UnitTestSetup<pce_type, fad_type> setup;
  PCE_SERIALIZATION_TESTS(pce_type, fad_type, OrthogPoly)
}

namespace ETPCETest {
  Sacado::Random<double> rnd;
  typedef Sacado::ETPCE::OrthogPoly<double,storage_type> pce_type;
  UnitTestSetup<pce_type, fad_type> setup;
  PCE_SERIALIZATION_TESTS(pce_type, fad_type, ETOrthogPoly)
}

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
