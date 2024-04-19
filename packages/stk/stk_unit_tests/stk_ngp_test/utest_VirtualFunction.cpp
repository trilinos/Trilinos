// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include <stk_ngp_test/ngp_test.hpp>

namespace ngp
{
class NgpBase
{
 public:
  KOKKOS_DEFAULTED_FUNCTION NgpBase() = default;
  KOKKOS_FUNCTION virtual ~NgpBase() {}

  virtual void host_function() = 0;
};

template <typename T>
class NgpDerived : public NgpBase
{
 public:
  KOKKOS_FUNCTION
  NgpDerived() : NgpBase() {}

  KOKKOS_FUNCTION
  ~NgpDerived() {}

  virtual void host_function() {}
};

struct SimpleStruct {
  KOKKOS_FUNCTION
  void print() {
#if KOKKOS_VERSION < 40200
    printf("Printing from A located at %p\n", static_cast<void*>(this));
#else
    Kokkos::printf("Printing from A located at %p\n", static_cast<void*>(this));
#endif
  }
};

struct BaseStruct {
  virtual void set_i(const int) = 0;
  KOKKOS_FUNCTION
  virtual void print() {
#if KOKKOS_VERSION < 40200
    printf("Printing from base located at %p\n", static_cast<void*>(this));
#else
    Kokkos::printf("Printing from base located at %p\n", static_cast<void*>(this));
#endif
  }
};

struct ChildStruct : public BaseStruct {
  int i;
  virtual void set_i(const int _i) { i = _i; }
  KOKKOS_FUNCTION
  virtual void print() {
#if KOKKOS_VERSION < 40200
    printf("Printing from child located at %p with i %i\n", static_cast<void*>(this), i); }
#else
    Kokkos::printf("Printing from child located at %p with i %i\n", static_cast<void*>(this), i); }
#endif
};

}  // namespace ngp

namespace
{
void test_device_class()
{
  int constructionFinished = 0;
  Kokkos::parallel_reduce(
      stk::ngp::DeviceRangePolicy(0, 1),
      KOKKOS_LAMBDA(const unsigned& i, int& localFinished) {
        ngp::NgpDerived<int> derivedClass;
        localFinished = 1;
      },
      constructionFinished);

  EXPECT_EQ(1, constructionFinished);
}

TEST(NgpDevice, virtualFunction)
{
  test_device_class();
}

void testSimpleFunction()
{
  ngp::SimpleStruct* devicePointer = static_cast<ngp::SimpleStruct*>(Kokkos::kokkos_malloc(sizeof(ngp::SimpleStruct)));
  Kokkos::parallel_for(
      stk::ngp::DeviceRangePolicy(0, 1), KOKKOS_LAMBDA(const int) { new (devicePointer) ngp::SimpleStruct(); });

  Kokkos::fence();

  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1), [devicePointer] KOKKOS_FUNCTION(const int) { devicePointer->print(); });

  Kokkos::fence();
  Kokkos::kokkos_free(devicePointer);
}

TEST(PlacementNew, simple)
{
  testSimpleFunction();
}

void testVirtualFunction()
{
  ngp::BaseStruct* devicePointer = static_cast<ngp::BaseStruct*>(Kokkos::kokkos_malloc(sizeof(ngp::ChildStruct)));
  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1), [devicePointer] KOKKOS_FUNCTION(const int) { new (devicePointer) ngp::ChildStruct(); });

  Kokkos::fence();

  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1), [devicePointer] KOKKOS_FUNCTION(const int) { devicePointer->print(); });

  Kokkos::fence();
  Kokkos::kokkos_free(devicePointer);
}

TEST(PlacementNew, Virtual)
{
  testVirtualFunction();
}

void testVirtualFunctionWithCopyConstructor()
{
  ngp::ChildStruct hostObject;
  hostObject.set_i(3);
  const ngp::ChildStruct& hostReference = hostObject;
  ngp::BaseStruct* devicePointer = static_cast<ngp::BaseStruct*>(Kokkos::kokkos_malloc(sizeof(ngp::ChildStruct)));
  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1), [devicePointer, hostReference] KOKKOS_FUNCTION(
                              const int) { new (devicePointer) ngp::ChildStruct(hostReference); });

  Kokkos::fence();

  Kokkos::parallel_for(stk::ngp::DeviceRangePolicy(0, 1), [devicePointer] KOKKOS_FUNCTION(const int) { devicePointer->print(); });

  Kokkos::fence();
  Kokkos::kokkos_free(devicePointer);
}

TEST(PlacementNew, VirtualCopy)
{
  testVirtualFunctionWithCopyConstructor();
}

}  // namespace
