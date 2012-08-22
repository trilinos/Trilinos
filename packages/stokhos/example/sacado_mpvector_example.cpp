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

#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include "sacado_mpvector_example.hpp"

int main(int argc, char **argv)
{
  try {

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
      bool status = MPVectorExample<MaxSize,KokkosArray::Cuda>::run(
	storage_method, n, sz, nblocks, nthreads, reset, print);

      if (status)
	std::cout << "CUDA Test Passed!" << std::endl;
      else
	std::cout << "CUDA Test Failed!" << std::endl;
    }
#endif

     if (test_host) {
       bool status = MPVectorExample<MaxSize,KokkosArray::Host>::run(
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
