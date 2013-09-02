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

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#include "sacado_mpvector_example.hpp"

int main(int argc, char **argv)
{
  try {

    Teuchos::GlobalMPISession mpiSession(&argc, &argv);

    // Setup command line options
    Teuchos::CommandLineProcessor CLP;
    CLP.setDocString(
      "This example explores operator overloading on CUDA.\n");

    Storage_Method storage_method = STATIC_FIXED;
    CLP.setOption("storage_method", &storage_method, 
		  num_storage_method, storage_method_values, 
		  storage_method_names, "Vector storage method");

    int n = 100;
    CLP.setOption("num_loop", &n, "Number of loop iterations");

    int sz = MaxSize;
    CLP.setOption("num_vec", &sz, "Number of vector components");

    int nblocks = 400;
    CLP.setOption("num_blocks", &nblocks, "Number of thread blocks");

    int nthreads = 64;
    CLP.setOption("num_threads", &nthreads, "Number of threads per block");

    bool reset = true;
    CLP.setOption("reset", "no-reset", &reset, "Reset initial vectors");

    bool print = false;
    CLP.setOption("print", "quiet", &print, "Print values");

    bool test_host = true;
    CLP.setOption("host", "no-host", &test_host, "Test host");

#ifdef HAVE_KOKKOSCLASSIC_CUDA
    bool test_cuda = true;
    CLP.setOption("cuda", "no-cuda", &test_cuda, "Test CUDA");
#else
    bool test_cuda = false;
#endif

    CLP.parse( argc, argv );

    std::cout << "Summary of command line options:" << std::endl
	      << "\tstorage_method = " 
	      << storage_method_names[storage_method] << std::endl
	      << "\tnum_loop    = " << n << std::endl
	      << "\tnum_vec     = " << sz << std::endl
	      << "\tnum_blocks  = " << nblocks << std::endl
	      << "\tnum_threads = " << nthreads << std::endl
	      << "\treset       = " << reset << std::endl
	      << "\tprint       = " << print << std::endl
	      << "\thost        = " << test_host << std::endl
	      << "\tcuda        = " << test_cuda << std::endl << std::endl;

#ifdef HAVE_KOKKOSCLASSIC_CUDA
    if (test_cuda) {
      bool status = MPVectorExample<MaxSize,Kokkos::Cuda>::run(
	storage_method, n, sz, nblocks, nthreads, reset, print);

      if (status)
	std::cout << "CUDA Test Passed!" << std::endl;
      else
	std::cout << "CUDA Test Failed!" << std::endl;
    }
#endif

     if (test_host) {
       bool status = MPVectorExample<MaxSize,Kokkos::Threads>::run(
	 storage_method, n, sz, nblocks, nthreads, reset, print);
       
       if (status)
	 std::cout << "Host Test Passed!" << std::endl;
       else
	 std::cout << "Host Test Failed!" << std::endl;
    }

    Teuchos::TimeMonitor::summarize(std::cout);
    Teuchos::TimeMonitor::zeroOutTimers();

  }
  catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
  catch (std::string& s) {
    std::cout << s << std::endl;
  }
  catch (char *s) {
    std::cout << s << std::endl;
  }
  catch (...) {
    std::cout << "Caught unknown exception!" <<std::endl;
  }
}
