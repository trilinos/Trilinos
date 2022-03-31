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

#include "Tpetra_ConfigDefs.hpp"
#include "Teuchos_LocalTestingHelpers.hpp"
#include "Teuchos_FancyOStream.hpp"
#include <stdio.h>
#include <mpi.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_View.hpp>
#include <unistd.h>
#include <map>
#include <string>
#include <list>
#include "ApiTest.h"



class DeepCopyTester {
public:
  DeepCopyTester() {}
  void run(int argc, char *argv[])  {
    log("Tpetra Kokkos DeepCopy Regression test start");

    std::string cudaSync("cudaDeviceSynchronize");
    std::string cudaMemcpy("cudaMemcpy");
    std::string cudaMemcpyAsync("cudaMemcpyAsync");
    ApiTest *counter = ApiTest::getInstance();
    //initialize
#ifdef HAVE_TPETRACORE_MPI
    MPI_Init(&argc, &argv);
#endif
    Kokkos::initialize(argc, argv);
    {
      isConsistent = true;
      const int N = 100;
      Kokkos::View<int*, Kokkos::LayoutLeft, Kokkos::HostSpace> a ("a", N);
      Kokkos::View<int*, Kokkos::LayoutLeft, Kokkos::HostSpace> b ("b", N);
      Kokkos::View<int*, Kokkos::LayoutLeft, Kokkos::CudaSpace> c ("c", N);
      Kokkos::View<int*, Kokkos::LayoutLeft, Kokkos::CudaSpace> d ("d", N);

      counter->map_zero();
      counter->setExpectation(cudaSync, 2);
      OnHost2Arg(a, b);
      if (!counter->testExpectations()) {
        log("OnHost2Arg()",counter);
	isConsistent = false;
      }

      counter->map_zero();
      counter->setExpectation(cudaSync, 2);
      counter->setExpectation(cudaMemcpy, 1);
      OnDevice2Arg(c, d);
      if (!counter->testExpectations()) {
        log("OnDevice2Arg()",counter);
	isConsistent = false;
      }

      counter->map_zero();
      counter->setExpectation(cudaSync, 2);
      counter->setExpectation(cudaMemcpy, 1);
      HostToDevice2Arg(a, c);
      if (!counter->testExpectations()) {
        log("HostToDevice2Arg()",counter);
	isConsistent = false;
      }

      counter->map_zero();
      counter->setExpectation(cudaSync, 2);
      counter->setExpectation(cudaMemcpy, 1);
      DeviceToHost2Arg(c, a);
      if (!counter->testExpectations()) {
	log("DeviceToHost2Arg()",counter);
	isConsistent = false;
      }

      counter->map_zero();
      OnHost3Arg(a, b);
      if (!counter->testExpectations()) {
	log("OnHost3Arg()",counter);
	isConsistent = false;
      }

      counter->map_zero();
      counter->setExpectation(cudaMemcpyAsync, 1);
      OnDevice3Arg(c, d);
      if (!counter->testExpectations()) {
	log("OnDevice3Arg()",counter);
	isConsistent = false;
      }

      counter->map_zero();
      counter->setExpectation(cudaMemcpyAsync, 1);
      HostToDevice3Arg(a, c);
      if (!counter->testExpectations()) {
	log("HostToDevice3Arg()",counter);
	isConsistent = false;
      }

      counter->map_zero();
      counter->setExpectation(cudaMemcpyAsync, 1);
      DeviceToHost3Arg(c, a);
      if (!counter->testExpectations()) {
	log("DeviceToHost3Arg()",counter);
	isConsistent = false;
      }
    }
    Kokkos::finalize();
#ifdef HAVE_TPETRACORE_MPI
    MPI_Finalize();
#endif
    return;
  }

  bool getConsistency() {
    return isConsistent;
  }

  void OnHost2Arg(Kokkos::View<int*, Kokkos::LayoutLeft, Kokkos::HostSpace> &a,
		  Kokkos::View<int*, Kokkos::LayoutLeft, Kokkos::HostSpace> &b) {
    /*  using range_policy = Kokkos::RangePolicy<Kokkos::HostSpace>;

	Kokkos::parallel_for("onHost initialize", range_policy(0,N), KOKKOS_LAMBDA(const int &i) {
	a[i] = 2;
	});
    */
    Kokkos::deep_copy(b, a);
  }

  void OnDevice2Arg(Kokkos::View<int*, Kokkos::LayoutLeft, Kokkos::CudaSpace> &a,
		    Kokkos::View<int*, Kokkos::LayoutLeft, Kokkos::CudaSpace> &b) {
    using range_policy = Kokkos::RangePolicy<Kokkos::Cuda>;
    /*
      Kokkos::parallel_for("onDevice initialize", range_policy(0,N), KOKKOS_LAMBDA(const int &i) {
      a[i] = 2;
      });
    */
    Kokkos::deep_copy(b, a);
  }

  void HostToDevice2Arg(Kokkos::View<int*, Kokkos::LayoutLeft, Kokkos::HostSpace> &a,
			Kokkos::View<int*, Kokkos::LayoutLeft, Kokkos::CudaSpace> &b) {
    /*
      Kokkos::parallel_for("host to device initialize", N, KOKKOS_LAMBDA(const int &i) {
      a[i] = 2;
      });
    */
    Kokkos::deep_copy(b, a);
  }

  void DeviceToHost2Arg(Kokkos::View<int*, Kokkos::LayoutLeft, Kokkos::CudaSpace> &a,
			Kokkos::View<int*, Kokkos::LayoutLeft, Kokkos::HostSpace> &b) {
    /*
      Kokkos::parallel_for("device to host initialize", N, KOKKOS_LAMBDA(const int &i) {
      a[i] = 2;
      });
    */
    Kokkos::deep_copy(b, a);
  }

  void OnHost3Arg(Kokkos::View<int*, Kokkos::LayoutLeft, Kokkos::HostSpace> &a,
		  Kokkos::View<int*, Kokkos::LayoutLeft, Kokkos::HostSpace> &b) {
    Kokkos::DefaultExecutionSpace mySpace;
    /*
      using range_policy = Kokkos::RangePolicy<Kokkos::Host>;
      
      Kokkos::parallel_for("onHost initialize", range_policy(0,N), KOKKOS_LAMBDA(const int &i) {
      a[i] = 2;
      });
    */
    Kokkos::deep_copy(mySpace, b, a);
  }

  void OnDevice3Arg(Kokkos::View<int*, Kokkos::LayoutLeft, Kokkos::CudaSpace> &a,
		    Kokkos::View<int*, Kokkos::LayoutLeft, Kokkos::CudaSpace> &b) {
    Kokkos::DefaultExecutionSpace mySpace;
    using range_policy = Kokkos::RangePolicy<Kokkos::Cuda>;
    
    /*
      Kokkos::parallel_for("onDevice initialize", range_policy(0,N), KOKKOS_LAMBDA(const int &i) {
      a[i] = 2;
      });
    */
    Kokkos::deep_copy(mySpace, b, a);
  }

  void HostToDevice3Arg(Kokkos::View<int*, Kokkos::LayoutLeft, Kokkos::HostSpace> &a,
			Kokkos::View<int*, Kokkos::LayoutLeft, Kokkos::CudaSpace> &b) {
    Kokkos::DefaultExecutionSpace mySpace;
    /*
      Kokkos::parallel_for("host to device initialize", N, KOKKOS_LAMBDA(const int &i) {
      a[i] = 2;
      });
    */
    Kokkos::deep_copy(mySpace, b, a);
  }

  void DeviceToHost3Arg(Kokkos::View<int*, Kokkos::LayoutLeft, Kokkos::CudaSpace> &a,
			Kokkos::View<int*, Kokkos::LayoutLeft, Kokkos::HostSpace> &b) {
    Kokkos::DefaultExecutionSpace mySpace;
    /*
      Kokkos::parallel_for("device to host initialize", N, KOKKOS_LAMBDA(const int &i) {
      a[i] = 2;
      });
    */
    Kokkos::deep_copy(mySpace, b, a);
  }

private:
  bool isConsistent;
  void log(const std::string& msg, ApiTest *counter=nullptr) {
    std::cout << msg << std::endl;
    if(counter) counter->printAll();
  }
};

int main(int argc, char *argv[]) {
  bool success=true;

  std::ostream &out = std::cout; 

  
  DeepCopyTester tester;
  tester.run(argc,argv);
  TEST_EQUALITY_CONST(tester.getConsistency(),true);

  if (success)
    out << "\nEnd Result: TEST PASSED" << std::endl;
  else
    out << "\nEnd Result: TEST FAILED" << std::endl;


  return (success ? 0 : 1);
}

