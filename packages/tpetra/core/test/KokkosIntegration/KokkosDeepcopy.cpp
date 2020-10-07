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

#include <stdio.h>
#include <mpi.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_View.hpp>
#include <unistd.h>

//strategy: kokkos::deepcopy a view, count the intercepts.  For right now, just report them when mpi_finalize happens.
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

int main(int argc, char *argv[]) {
  //initialize
  MPI_Init(&argc, &argv);
  Kokkos::initialize(argc, argv);
  {
    const int N = 100;
    Kokkos::View<int*, Kokkos::LayoutLeft, Kokkos::HostSpace> a ("a", N);
    Kokkos::View<int*, Kokkos::LayoutLeft, Kokkos::HostSpace> b ("b", N);
    Kokkos::View<int*, Kokkos::LayoutLeft, Kokkos::CudaSpace> c ("c", N);
    Kokkos::View<int*, Kokkos::LayoutLeft, Kokkos::CudaSpace> d ("d", N);
  
    OnHost2Arg(a, b);
    OnDevice2Arg(c, d);
    HostToDevice2Arg(a, c);
    DeviceToHost2Arg(c, a);
    OnHost3Arg(a, b);
    OnDevice3Arg(c, d);
    HostToDevice3Arg(a, c);
    DeviceToHost3Arg(c, a);
  }
  fprintf(stderr, "summary statistics\n");
  Kokkos::finalize();
  MPI_Finalize();

  return 0;
}

