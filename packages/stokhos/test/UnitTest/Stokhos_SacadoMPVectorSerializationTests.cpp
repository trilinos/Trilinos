// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
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
#include "Stokhos_Sacado_Kokkos_MP_Vector.hpp"
#include "Sacado_Fad_DFad.hpp"
#include "Sacado_mpl_apply.hpp"
#include "Sacado_Random.hpp"

#include <Kokkos_Core.hpp>

//
// Currently this doesn't test:
//   * the device
//   * threaded storage (needs the device)
//   * strided storage with non-trivial stride
//

using Teuchos::RCP;
using Teuchos::rcp;

// Common setup for unit tests
template <typename VecType, typename FadType>
struct UnitTestSetup {
  int sz;

  typedef Teuchos::ValueTypeSerializer<int, VecType> VecSerializerT;
  RCP<VecSerializerT> vec_serializer;

  typedef typename Sacado::mpl::apply<FadType,VecType>::type FadVecType;
  typedef Teuchos::ValueTypeSerializer<int, FadVecType> FadVecSerializerT;
  RCP<FadVecSerializerT> fad_vec_serializer;

  UnitTestSetup() {
    sz = 8;

    // Serializers
    vec_serializer =
      rcp(new VecSerializerT(
            rcp(new Teuchos::ValueTypeSerializer<int,double>()), sz));
    fad_vec_serializer = rcp(new FadVecSerializerT(vec_serializer, 5));
  }
};

template <typename VecType>
bool testSerialization(const Teuchos::Array<VecType>& x,
                       const std::string& tag,
                       Teuchos::FancyOStream& out) {

  typedef int Ordinal;
  typedef Teuchos::SerializationTraits<Ordinal,VecType> SerT;

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
  Teuchos::Array<VecType> x2(count2);
  for (Ordinal i=0; i<count2; i++)
    x2[i].reset(x[i].size());
  SerT::deserialize(bytes, charBuffer, count2, &x2[0]);

  delete [] charBuffer;

  // Check coefficients match
  for (Ordinal i=0; i<count; i++) {
    bool success2 = Sacado::IsEqual<VecType>::eval(x[i], x2[i]);
    out << tag << " serialize/deserialize vec test " << i;
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

template <typename VecType, typename Serializer>
bool testSerialization(Teuchos::Array<VecType>& x,
                       const Serializer& serializer,
                       const std::string& tag,
                       Teuchos::FancyOStream& out) {

  typedef int Ordinal;

  // Serialize
  Ordinal count = x.size();
  Ordinal bytes = serializer.fromCountToIndirectBytes(count, &x[0]);
  char *charBuffer = new char[bytes];
  serializer.serialize(count, &x[0], bytes, charBuffer);

  // Reset x to given size
  for (Ordinal i=0; i<count; i++)
    x[i].reset(serializer.getSerializerSize());

  // Deserialize
  Ordinal count2 = serializer.fromIndirectBytesToCount(bytes, charBuffer);
  Teuchos::Array<VecType> x2(count2);
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
    bool success2 = Sacado::IsEqual<VecType>::eval(x[i], x2[i]);
    out << tag << " serialize/deserialize vec test " << i;
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

template <typename VecType, typename Serializer>
bool testNestedSerialization(Teuchos::Array<VecType>& x,
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
      x[i].fastAccessDx(j).reset(vs->getSerializerSize());
    x[i].val().reset(vs->getSerializerSize());
  }

  // Deserialize
  Ordinal count2 = serializer.fromIndirectBytesToCount(bytes, charBuffer);
  Teuchos::Array<VecType> x2(count2);
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
    bool success2 = Sacado::IsEqual<VecType>::eval(x[i], x2[i]);
    out << tag << " serialize/deserialize vec test " << i;
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

#define VEC_SERIALIZATION_TESTS(VecType, FadType, Vec)                  \
  TEUCHOS_UNIT_TEST( Vec##_Serialization, Uniform ) {                   \
    int n = 7;                                                          \
    Teuchos::Array<VecType> x(n);                                       \
    for (int i=0; i<n; i++) {                                           \
      x[i] = VecType(setup.sz, 0.0);                                    \
      for (int j=0; j<setup.sz; j++)                                    \
        x[i].fastAccessCoeff(j) = rnd.number();                         \
    }                                                                   \
    bool success1 = testSerialization(                                  \
      x, std::string(#Vec) + " Uniform", out);                          \
    bool success2 = testSerialization(                                  \
      x, *setup.vec_serializer, std::string(#Vec) + " Uniform PTS", out); \
    success = success1 && success2;                                     \
  }                                                                     \
                                                                        \
  TEUCHOS_UNIT_TEST( Vec##_Serialization, Empty ) {                     \
    int n = 7;                                                          \
    Teuchos::Array<VecType> x(n);                                       \
    for (int i=0; i<n; i++) {                                           \
      x[i] = VecType(1, 0.0);                                           \
      x[i].val() = rnd.number();                                        \
    }                                                                   \
    bool success1 = testSerialization(                                  \
      x, std::string(#Vec) + " Empty", out);                            \
    bool success2 = testSerialization(                                  \
      x, *setup.vec_serializer, std::string(#Vec) + " Empty PTS", out); \
    success = success1 && success2;                                     \
  }                                                                     \
                                                                        \
  TEUCHOS_UNIT_TEST( Vec##_Serialization, Mixed ) {                     \
    int n = 6;                                                          \
    int p[] = { 5, 0, 8, 8, 3, 0 };                                     \
    Teuchos::Array<VecType> x(n);                                       \
    for (int i=0; i<n; i++) {                                           \
      x[i] = VecType(p[i], 0.0);                                        \
      for (int j=0; j<p[i]; j++)                                        \
        x[i].fastAccessCoeff(j) = rnd.number();                         \
    }                                                                   \
    bool success1 = testSerialization(                                  \
      x, std::string(#Vec) + " Mixed", out);                            \
    bool success2 = testSerialization(                                  \
      x, *setup.vec_serializer, std::string(#Vec) + " Mixed PTS", out); \
    success = success1 && success2;                                     \
  }                                                                     \
                                                                        \
  TEUCHOS_UNIT_TEST( Vec##_Serialization, FadVecUniform ) {             \
    typedef Sacado::mpl::apply<FadType,VecType>::type FadVecType;       \
    int n = 7;                                                          \
    int p = 3;                                                          \
    Teuchos::Array<FadVecType> x(n);                                    \
    for (int i=0; i<n; i++) {                                           \
      VecType f(setup.sz, 0.0);                                         \
      for (int k=0; k<setup.sz; k++)                                    \
        f.fastAccessCoeff(k) = rnd.number();                            \
      x[i] = FadVecType(p, f);                                          \
      for (int j=0; j<p; j++) {                                         \
        VecType g(setup.sz, 0.0);                                       \
        for (int k=0; k<setup.sz; k++)                                  \
          g.fastAccessCoeff(k) = rnd.number();                          \
        x[i].fastAccessDx(j) = g;                                       \
      }                                                                 \
    }                                                                   \
    success =                                                           \
      testNestedSerialization(x, *setup.fad_vec_serializer,             \
                              std::string(#Vec) + " Nested Uniform", out); \
  }                                                                     \
  TEUCHOS_UNIT_TEST( Vec##_Serialization, FadVecEmptyInner ) {          \
    typedef Sacado::mpl::apply<FadType,VecType>::type FadVecType;       \
    int n = 7;                                                          \
    int p = 3;                                                          \
    Teuchos::Array<FadVecType> x(n);                                    \
    for (int i=0; i<n; i++) {                                           \
      VecType f(setup.sz, 0.0);                                         \
      for (int k=0; k<setup.sz; k++)                                    \
        f.fastAccessCoeff(k) = rnd.number();                            \
      x[i] = FadVecType(p, f);                                          \
      for (int j=0; j<p; j++)                                           \
        x[i].fastAccessDx(j) =  rnd.number();                           \
    }                                                                   \
    success =                                                           \
      testNestedSerialization(                                          \
        x, *setup.fad_vec_serializer,                                   \
        std::string(#Vec) + " Nested Empty Inner", out);                \
  }                                                                     \
  TEUCHOS_UNIT_TEST( Vec##_Serialization, FadVecEmptyOuter ) {          \
    typedef Sacado::mpl::apply<FadType,VecType>::type FadVecType;       \
    int n = 7;                                                          \
    Teuchos::Array<FadVecType> x(n);                                    \
    for (int i=0; i<n; i++) {                                           \
      VecType f(setup.sz, 0.0);                                         \
      for (int k=0; k<setup.sz; k++)                                    \
        f.fastAccessCoeff(k) = rnd.number();                            \
      x[i] = FadVecType(f);                                             \
    }                                                                   \
    success =                                                           \
      testNestedSerialization(                                          \
        x, *setup.fad_vec_serializer,                                   \
        std::string(#Vec) + " Nested Empty Outer", out);                \
  }                                                                     \
  TEUCHOS_UNIT_TEST( Vec##_Serialization, FadVecEmptyAll ) {            \
    typedef Sacado::mpl::apply<FadType,VecType>::type FadVecType;       \
    int n = 7;                                                          \
    Teuchos::Array<FadVecType> x(n);                                    \
    for (int i=0; i<n; i++) {                                           \
      x[i] = rnd.number();                                              \
    }                                                                   \
    success =                                                           \
      testNestedSerialization(                                          \
        x, *setup.fad_vec_serializer,                                   \
        std::string(#Vec) + " Nested Empty All", out);                  \
  }

namespace DynamicVecTest {
  Sacado::Random<double> rnd;
  typedef Kokkos::Threads device_type;
  typedef Stokhos::DynamicStorage<int,double,device_type> storage_type;
  typedef Sacado::Fad::DFad<double> fad_type;
  typedef Sacado::MP::Vector<storage_type> vec_type;
  UnitTestSetup<vec_type, fad_type> setup;
  VEC_SERIALIZATION_TESTS(vec_type, fad_type, DynamicVector)
}

namespace DynamicStridedVecTest {
  Sacado::Random<double> rnd;
  typedef Kokkos::Threads device_type;
  typedef Stokhos::DynamicStridedStorage<int,double,device_type> storage_type;
  typedef Sacado::Fad::DFad<double> fad_type;
  typedef Sacado::MP::Vector<storage_type> vec_type;
  UnitTestSetup<vec_type, fad_type> setup;
  VEC_SERIALIZATION_TESTS(vec_type, fad_type, DynamicStridedVector)
}

namespace StaticVecTest {
  Sacado::Random<double> rnd;
  typedef Kokkos::Threads device_type;
  typedef Stokhos::StaticStorage<int,double,8,device_type> storage_type;
  typedef Sacado::Fad::DFad<double> fad_type;
  typedef Sacado::MP::Vector<storage_type> vec_type;
  UnitTestSetup<vec_type, fad_type> setup;
  VEC_SERIALIZATION_TESTS(vec_type, fad_type, StaticVector)
}

namespace StaticFixedVecTest {
  Sacado::Random<double> rnd;
  typedef Kokkos::Threads device_type;
  typedef Stokhos::StaticFixedStorage<int,double,8,device_type> storage_type;
  typedef Sacado::Fad::DFad<double> fad_type;
  typedef Sacado::MP::Vector<storage_type> vec_type;
  UnitTestSetup<vec_type, fad_type> setup;
  VEC_SERIALIZATION_TESTS(vec_type, fad_type, StaticFixedVector)
}

int main( int argc, char* argv[] ) {
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
