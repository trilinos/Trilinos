/*
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
// ************************************************************************
// @HEADER
*/

#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_Details_minMaxNumRowEntries.hpp"
#include "Kokkos_Core.hpp"

namespace { // (anonymous)

  //
  // UNIT TESTS
  //

  template<class LocalOrdinal, class RowOffset>
  void
  testMinMaxNumRowEntries(bool& success, Teuchos::FancyOStream& out)
  {
    using LO = LocalOrdinal;
    using OT = RowOffset;
    using Tpetra::Details::minMaxNumRowEntries;
    Kokkos::HostSpace host;

    {
      const LO lclNumRows = 0;

      Kokkos::View<OT*> ptr1("ptr1", 0); // valid for this to be empty
      auto result = minMaxNumRowEntries("minMax", ptr1, lclNumRows);

      TEST_EQUALITY( result.first, LO(0) );
      TEST_EQUALITY( result.second, LO(0) );

      Kokkos::View<OT*> ptr2("ptr2", lclNumRows+1);
      result = minMaxNumRowEntries("minMax", ptr1, lclNumRows);

      TEST_EQUALITY( result.first, LO(0) );
      TEST_EQUALITY( result.second, LO(0) );
    }

    {
      const LO lclNumRows = 1;

      Kokkos::View<OT*> ptr("ptr", lclNumRows+1);
      auto result = minMaxNumRowEntries("minMax", ptr, lclNumRows);

      TEST_EQUALITY( result.first, LO(0) );
      TEST_EQUALITY( result.second, LO(0) );

      auto ptr_h = Kokkos::create_mirror_view(host, ptr);
      ptr_h[0] = OT(0);
      ptr_h[1] = OT(5);
      Kokkos::deep_copy(ptr, ptr_h);

      result = minMaxNumRowEntries("minMax", ptr, lclNumRows);
      TEST_EQUALITY( result.first, LO(5) );
      TEST_EQUALITY( result.second, LO(5) );
    }

    {
      const LO lclNumRows = 2;

      Kokkos::View<OT*> ptr("ptr", lclNumRows+1);
      auto result = minMaxNumRowEntries("minMax", ptr, lclNumRows);

      TEST_EQUALITY( result.first, LO(0) );
      TEST_EQUALITY( result.second, LO(0) );

      auto ptr_h = Kokkos::create_mirror_view(host, ptr);

      ptr_h[0] = OT(0);
      ptr_h[1] = OT(5);
      ptr_h[2] = OT(5);
      Kokkos::deep_copy(ptr, ptr_h);

      result = minMaxNumRowEntries("minMax", ptr, lclNumRows);
      TEST_EQUALITY( result.first, LO(0) );
      TEST_EQUALITY( result.second, LO(5) );

      ptr_h[0] = OT(0);
      ptr_h[1] = OT(6);
      ptr_h[2] = OT(11);
      Kokkos::deep_copy(ptr, ptr_h);

      result = minMaxNumRowEntries("minMax", ptr, lclNumRows);
      TEST_EQUALITY( result.first, LO(5) );
      TEST_EQUALITY( result.second, LO(6) );
    }
  }

  TEUCHOS_UNIT_TEST( TpetraUtils, MinMaxNumRowEntries )
  {
    testMinMaxNumRowEntries<int, size_t>(success, out);
    testMinMaxNumRowEntries<long long, size_t>(success, out);
    testMinMaxNumRowEntries<int, ptrdiff_t>(success, out);
    testMinMaxNumRowEntries<int, unsigned int>(success, out);
  }

} // namespace (anonymous)

int
main(int argc, char* argv[])
{
  int errCode = 0;
  {
    Kokkos::ScopeGuard kokkosScope(argc, argv);
    errCode = Teuchos::UnitTestRepository::
      runUnitTestsFromMain(argc, argv);
  }
  return errCode;
}
