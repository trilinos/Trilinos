/*
//@HEADER
// ************************************************************************
//
//   Kokkos: Manycore Performance-Portable Multidimensional Arrays
//              Copyright (2012) Sandia Corporation
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
// Questions? Contact  H. Carter Edwards (hcedwar@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>

#if !defined(KOKKOS_HAVE_CUDA) || defined(__CUDACC__)
//----------------------------------------------------------------------------

namespace Test {

namespace Impl {

  char** init_kokkos_args(bool do_threads,bool do_numa,bool do_device,bool do_other, int& nargs) {
    nargs = (do_threads?1:0) +
            (do_numa?1:0) +
            (do_device?1:0) +
            (do_other?4:0);
    char** args_kokkos = new char*[nargs];
    for(int i = 0; i < nargs; i++)
      args_kokkos[i] = new char[20];

    int threads_idx = do_other?1:0;
    int numa_idx = (do_other?3:0) + (do_threads?1:0);
    int device_idx = (do_other?3:0) + (do_threads?1:0) + (do_numa?1:0);


    if(do_threads) {
      int nthreads = 3;

#ifdef _OPENMP
      if(omp_get_max_threads() < 3)
        nthreads = omp_get_max_threads();
#endif

      if(Kokkos::hwloc::available())  {
        if(Kokkos::hwloc::get_available_threads_per_core()<3)
            nthreads =   Kokkos::hwloc::get_available_threads_per_core()
                       * Kokkos::hwloc::get_available_numa_count();
      }

      sprintf(args_kokkos[threads_idx],"--threads=%i",nthreads);
    }

    if(do_numa) {
      int numa = 1;
      if(Kokkos::hwloc::available())
        numa = Kokkos::hwloc::get_available_numa_count();
      sprintf(args_kokkos[numa_idx],"--numa=%i",numa);
    }

    if(do_device) {
      sprintf(args_kokkos[device_idx],"--device=%i",0);
    }

    if(do_other) {
      sprintf(args_kokkos[0],"--dummyarg=1");
      sprintf(args_kokkos[threads_idx+(do_threads?1:0)],"--dummy2arg");
      sprintf(args_kokkos[threads_idx+(do_threads?1:0)+1],"dummy3arg");
      sprintf(args_kokkos[device_idx+(do_device?1:0)],"dummy4arg=1");
    }

    return args_kokkos;
  }

  Kokkos::InitArguments init_initstruct(bool do_threads, bool do_numa, bool do_device) {
    Kokkos::InitArguments args;

    if(do_threads) {
      int nthreads = 3;

#ifdef _OPENMP
      if(omp_get_max_threads() < 3)
        nthreads = omp_get_max_threads();
#endif

      if(Kokkos::hwloc::available())  {
        if(Kokkos::hwloc::get_available_threads_per_core()<3)
            nthreads =   Kokkos::hwloc::get_available_threads_per_core()
                       * Kokkos::hwloc::get_available_numa_count();
      }

      args.nthreads = nthreads;
    }

    if(do_numa) {
      int numa = 1;
      if(Kokkos::hwloc::available())
        numa = Kokkos::hwloc::get_available_numa_count();
      args.nnuma = numa;
    }

    if(do_device) {
      args.device = 0;
    }

    return args;
  }

  //ToDo: Add check whether correct number of threads are actually started
  void test_no_arguments() {
    Kokkos::initialize();
    Kokkos::finalize();
  }

  void test_commandline_args(int nargs, char** args) {
    Kokkos::initialize(nargs,args);
    Kokkos::finalize();
  }

  void test_initstruct_args(const Kokkos::InitArguments& args) {
    Kokkos::initialize(args);
    Kokkos::finalize();
  }
}

class defaultdevicetypeinit : public ::testing::Test {
protected:
  static void SetUpTestCase()
  {
  }

  static void TearDownTestCase()
  {
  }
};


TEST_F( defaultdevicetypeinit, no_args) {
  Impl::test_no_arguments();
}

TEST_F( defaultdevicetypeinit, commandline_args_empty) {
  int nargs = 0;
  char** args = Impl::init_kokkos_args(false,false,false,false,nargs);
  Impl::test_commandline_args(nargs,args);
  for(int i = 0; i < nargs; i++)
    delete [] args[i];
  delete [] args;
}

TEST_F( defaultdevicetypeinit, commandline_args_other) {
  int nargs = 0;
  char** args = Impl::init_kokkos_args(false,false,false,true,nargs);
  Impl::test_commandline_args(nargs,args);
  for(int i = 0; i < nargs; i++)
    delete [] args[i];
  delete [] args;
}

TEST_F( defaultdevicetypeinit, commandline_args_nthreads) {
  int nargs = 0;
  char** args = Impl::init_kokkos_args(true,false,false,false,nargs);
  Impl::test_commandline_args(nargs,args);
  for(int i = 0; i < nargs; i++)
    delete [] args[i];
  delete [] args;
}

TEST_F( defaultdevicetypeinit, commandline_args_nthreads_numa) {
  int nargs = 0;
  char** args = Impl::init_kokkos_args(true,true,false,false,nargs);
  Impl::test_commandline_args(nargs,args);
  for(int i = 0; i < nargs; i++)
    delete [] args[i];
  delete [] args;
}

TEST_F( defaultdevicetypeinit, commandline_args_nthreads_numa_device) {
  int nargs = 0;
  char** args = Impl::init_kokkos_args(true,true,true,false,nargs);
  Impl::test_commandline_args(nargs,args);
  for(int i = 0; i < nargs; i++)
    delete [] args[i];
  delete [] args;
}

TEST_F( defaultdevicetypeinit, commandline_args_nthreads_device) {
  int nargs = 0;
  char** args = Impl::init_kokkos_args(true,false,true,false,nargs);
  Impl::test_commandline_args(nargs,args);
  for(int i = 0; i < nargs; i++)
    delete [] args[i];
  delete [] args;
}

TEST_F( defaultdevicetypeinit, commandline_args_numa_device) {
  int nargs = 0;
  char** args = Impl::init_kokkos_args(false,true,true,false,nargs);
  Impl::test_commandline_args(nargs,args);
  for(int i = 0; i < nargs; i++)
    delete [] args[i];
  delete [] args;
}

TEST_F( defaultdevicetypeinit, commandline_args_device) {
  int nargs = 0;
  char** args = Impl::init_kokkos_args(false,false,true,false,nargs);
  Impl::test_commandline_args(nargs,args);
  for(int i = 0; i < nargs; i++)
    delete [] args[i];
  delete [] args;
}

TEST_F( defaultdevicetypeinit, commandline_args_nthreads_numa_device_other) {
  int nargs = 0;
  char** args = Impl::init_kokkos_args(true,true,true,true,nargs);
  Impl::test_commandline_args(nargs,args);
  for(int i = 0; i < nargs; i++)
    delete [] args[i];
  delete [] args;
}

TEST_F( defaultdevicetypeinit, initstruct_default) {
  Kokkos::InitArguments args;
  Impl::test_initstruct_args(args);
}

TEST_F( defaultdevicetypeinit, initstruct_nthreads) {
  Kokkos::InitArguments args = Impl::init_initstruct(true,false,false);
  Impl::test_initstruct_args(args);
}

TEST_F( defaultdevicetypeinit, initstruct_nthreads_numa) {
  Kokkos::InitArguments args = Impl::init_initstruct(true,true,false);
  Impl::test_initstruct_args(args);
}

TEST_F( defaultdevicetypeinit, initstruct_device) {
  Kokkos::InitArguments args = Impl::init_initstruct(false,false,true);
  Impl::test_initstruct_args(args);
}

TEST_F( defaultdevicetypeinit, initstruct_nthreads_device) {
  Kokkos::InitArguments args = Impl::init_initstruct(true,false,true);
  Impl::test_initstruct_args(args);
}


TEST_F( defaultdevicetypeinit, initstruct_nthreads_numa_device) {
  Kokkos::InitArguments args = Impl::init_initstruct(true,true,true);
  Impl::test_initstruct_args(args);
}



} // namespace test

#endif
