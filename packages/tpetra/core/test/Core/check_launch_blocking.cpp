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

#include "Tpetra_Core.hpp"
#include "Kokkos_Core.hpp"
#include <cstring>

int main(int argc, char** argv)
{
  bool success = true;
  //These architecture macros are defined in KokkosCore_config.h (included by Core)
#if defined(KOKKOS_ARCH_KEPLER) || defined(KOKKOS_ARCH_MAXWELL)
  constexpr bool initializeShouldThrow = true;
#else
  constexpr bool initializeShouldThrow = false;
#endif
  bool threw = false;
  try
  {
    Tpetra::initialize(&argc, &argv);
  }
  catch(std::exception& e)
  {
    //Make sure the error's message is the one about launch blocking
    //(it mentions CUDA_LAUNCH_BLOCKING verbatim)
    //Otherwise it might be an unrelated exception, in which case this test should fail.
    if(!strstr(e.what(), "CUDA_LAUNCH_BLOCKING"))
    {
      std::cerr << "TEST FAILED: Tpetra::initialize() threw unrelated exception.\n";
      success = false;
    }
    else
    {
      std::cout << "Tpetra::initialize threw exception because CUDA_LAUNCH_BLOCKING required but not set.\n";
    }
    threw = true;
  }
  if(!threw)
  {
    //Initialization succeeded, so clean up
    Tpetra::finalize();
  }
  else
  {
    //Tpetra wasn't fully initialized, but Kokkos and MPI might have been.
    //Finalize those to avoid extra error messages when test terminates.
    if(Kokkos::is_initialized())
    {
      //Tpetra not fully initialized, but Kokkos was
      Kokkos::finalize();
    }
#ifdef HAVE_TPETRACORE_MPI
    int mpiInitialized = 0;
    MPI_Initialized(&mpiInitialized);
    if(mpiInitialized)
    {
      //Tpetra not fully initialized, but MPI was
      MPI_Finalize();
    }
#endif
  }
  if(threw && !initializeShouldThrow)
  {
    std::cerr << "TEST FAILED: Tpetra::initialize() threw an exception when it shouldn't have\n";
    std::cerr << "(CUDA arch is >= Pascal, so launch blocking not needed).\n";
    success = false;
  }
  else if(!threw && initializeShouldThrow)
  {
    std::cerr << "TEST FAILED: Tpetra::initialize() did not catch CUDA_LAUNCH_BLOCKING being unset.\n";
    success = false;
  }
  if(success)
    std::cout << "End Result: TEST PASSED\n";
  return 0;
}

