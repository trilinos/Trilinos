// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// ************************************************************************
// @HEADER

#include <TpetraKernels_Blas1_MV_UnitTests.hpp>

#ifdef KOKKOS_HAVE_PTHREAD

namespace KokkosBlas {
namespace Impl {

#define TPETRAKERNELS_BLAS1_MV_DEVICE Kokkos::Device<Kokkos::Threads, Kokkos::HostSpace>

template<>
bool
testOverScalarsAndLayouts<TPETRAKERNELS_BLAS1_MV_DEVICE> (std::ostream& out,
                                                          const char deviceName[],
                                                          const int numCols,
                                                          const bool oneCol,
                                                          const bool testComplex)
{
  using std::endl;
  typedef TPETRAKERNELS_BLAS1_MV_DEVICE device_type;
  bool success = true;

  if (device_type::execution_space::is_initialized ()) {
    out << endl << "Testing Device = " << deviceName << endl;
    bool curSuccess = true;

    curSuccess = testOverLayouts<double, device_type> (out, numCols, oneCol);
    success = success && curSuccess;
    curSuccess = testOverLayouts<float, device_type> (out, numCols, oneCol);
    success = success && curSuccess;
    curSuccess = testOverLayouts<int, device_type> (out, numCols, oneCol);
    success = success && curSuccess;

    if (testComplex) {
      curSuccess = testOverLayouts<Kokkos::complex<float>, device_type> (out, numCols, oneCol);
      success = success && curSuccess;
      curSuccess = testOverLayouts<Kokkos::complex<double>, device_type> (out, numCols, oneCol);
      success = success && curSuccess;
    }
  }
  else {
    out << "Device = " << deviceName << " NOT initialized; skipping test" << endl;
  }

  return success;
}

#undef TPETRAKERNELS_BLAS1_MV_DEVICE

} // namespace Impl
} // namespace KokkosBlas

#endif // KOKKOS_HAVE_PTHREAD
