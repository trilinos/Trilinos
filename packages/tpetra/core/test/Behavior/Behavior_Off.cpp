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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER
*/
#include <iostream>
#include <string>
#include <cstdlib>
#include <Tpetra_Core.hpp>
#include <Tpetra_TestingUtilities.hpp>
#include <Tpetra_Details_Behavior.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_CommHelpers.hpp>

/*
 * Note: The tests for the Behavior class are scattered across several files,
 * rather than being confined to a single file.  The reason is that the Behavior
 * class is instantiated one time only and once environment variables are read,
 * they are cached for future use.  Therefore, to test several values of an
 * environment variable, several tests need to be created (one for each distinct
 * value of the environment variable).
*/

namespace {

TEUCHOS_STATIC_SETUP()
{
  setenv("TPETRA_DEBUG", "OFF", 1);
  setenv("TPETRA_VERBOSE", "OFF", 1);
  setenv("TPETRA_ASSUME_CUDA_AWARE_MPI", "OFF", 1);
}

TEUCHOS_UNIT_TEST(Behavior, Off)
{

  // TPETRA_DEBUG was set globally in TEUCHOS_STATIC_SETUP to OFF, so any query
  // on TPETRA_DEBUG should evaluate to false, including named variants.
  bool dbg = Tpetra::Details::Behavior::debug();
  TEUCHOS_TEST_ASSERT(!dbg, out, success);
  bool dbg_named = Tpetra::Details::Behavior::debug("Named");
  TEUCHOS_TEST_ASSERT(!dbg_named, out, success);

  // TPETRA_VERBOSE was set globally in TEUCHOS_STATIC_SETUP to OFF, so any query
  // on TPETRA_VERBOSE should evaluate to false, including named variants.
  bool verb = Tpetra::Details::Behavior::verbose();
  TEUCHOS_TEST_ASSERT(!verb, out, success);
  bool verb_named = Tpetra::Details::Behavior::verbose("Named");
  TEUCHOS_TEST_ASSERT(!verb_named, out, success);

  // TPETRA_ASSUME_CUDA_AWARE_MPI was set globally in TEUCHOS_STATIC_SETUP to OFF,
  // so any query on TPETRA_ASSUME_CUDA_AWARE_MPI should evaluate to false.
  bool cuda_aware_mpi = Tpetra::Details::Behavior::assumeMpiIsCudaAware();
  TEUCHOS_TEST_ASSERT(!cuda_aware_mpi, out, success);
}
} // namespace (anonymous)
