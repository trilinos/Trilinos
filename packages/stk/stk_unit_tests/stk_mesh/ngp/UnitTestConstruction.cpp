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

#include <gtest/gtest.h>
#include <stk_mesh/base/NgpField.hpp>
#include <stk_ngp_test/ngp_test.hpp>

namespace {

#ifdef STK_USE_DEVICE_MESH
void test_device_field_default_constructor()
{
  int constructionFinished = 0;
  Kokkos::parallel_reduce(stk::ngp::DeviceRangePolicy(0, 1), KOKKOS_LAMBDA(const unsigned& /*i*/, int& localFinished) {
                            stk::mesh::DeviceField<double> deviceField;
                            NGP_EXPECT_EQ(stk::topology::INVALID_RANK, deviceField.get_rank());
                            localFinished = 1;
                          }, constructionFinished);

  EXPECT_EQ(1, constructionFinished);
}

TEST(NgpDeviceConstruction, deviceField)
{
  test_device_field_default_constructor();
}

TEST(NgpDeviceConstruction, deviceField_onHost)
{
  stk::mesh::DeviceField<double> deviceField;
  EXPECT_EQ(stk::topology::INVALID_RANK, deviceField.get_rank());
}
#endif

TEST(NgpDeviceConstruction, hostField)
{
  stk::mesh::HostField<double> hostField;
  EXPECT_EQ(stk::topology::INVALID_RANK, hostField.get_rank());
}

#ifdef STK_USE_DEVICE_MESH
struct MimicNaluWindKernelBase
{
  KOKKOS_DEFAULTED_FUNCTION MimicNaluWindKernelBase() = default;
  KOKKOS_DEFAULTED_FUNCTION MimicNaluWindKernelBase(const MimicNaluWindKernelBase&) = default;
  KOKKOS_FUNCTION virtual ~MimicNaluWindKernelBase()
  {
  }

  KOKKOS_FUNCTION virtual unsigned get_num() const { return 0; }
};

struct MimicNaluWindKernel : public MimicNaluWindKernelBase
{
  KOKKOS_FUNCTION MimicNaluWindKernel()
  : ngpField(), num(0)
  {
  }

  KOKKOS_FUNCTION MimicNaluWindKernel(const MimicNaluWindKernel& src)
  : ngpField(src.ngpField), num(src.num)
  {
  }

  KOKKOS_FUNCTION ~MimicNaluWindKernel()
  {
  }

  KOKKOS_FUNCTION unsigned get_num() const override { return num; }

  stk::mesh::NgpField<double> ngpField;
  unsigned num = 0;
};

void test_ngp_field_placement_new()
{
  MimicNaluWindKernel hostObj;
  hostObj.num = 42;

  std::string debugName("MimicNaluWindKernel");
  MimicNaluWindKernel* devicePtr = static_cast<MimicNaluWindKernel*>(Kokkos::kokkos_malloc<stk::ngp::MemSpace>(debugName, sizeof(MimicNaluWindKernel)));

  int constructionFinished = 0;
  Kokkos::parallel_reduce(stk::ngp::DeviceRangePolicy(0, 1), KOKKOS_LAMBDA(const unsigned& /*i*/, int& localFinished) {
    new (devicePtr) MimicNaluWindKernel(hostObj);
    localFinished = 1;
  }, constructionFinished);
  EXPECT_EQ(1, constructionFinished);

  int numFromDevice = 0;
  Kokkos::parallel_reduce(stk::ngp::DeviceRangePolicy(0, 1), KOKKOS_LAMBDA(const unsigned& /*i*/, int& localNum) {
    localNum = devicePtr->get_num();
  }, numFromDevice);
  EXPECT_EQ(42, numFromDevice);
  Kokkos::kokkos_free<stk::ngp::MemSpace>(devicePtr);
}

TEST(NgpDeviceConstruction, structWithNgpField)
{
  test_ngp_field_placement_new();
}
#endif

}

