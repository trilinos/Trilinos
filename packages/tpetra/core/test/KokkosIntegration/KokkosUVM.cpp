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

#include "Kokkos_Core.hpp"


class UVMTester {
public:
  UVMTester()
    : n(10),
      factorDevice(2),
      factorHost(5),
      factorCombined(factorDevice*factorHost),
      view("UVM View", n)
  { }

  void run() {
    log("Tpetra Kokkos UVM test start");

    fillOnHost();
    modifyOnDevice();
    modifyOnHost();

    if (isConsistent()) log("+++ Result is consistent +++");
    else log("!!! Result is inconsistent !!!");
  }

  void fillOnHost() {
    Kokkos::parallel_for(
        "fill on host",
        Kokkos::RangePolicy<Kokkos::Serial>(0,n),
        KOKKOS_LAMBDA(const int i) {
          view(i) = i;
        });
  }

  void modifyOnDevice() {
    Kokkos::parallel_for(
        "modify on device",
        Kokkos::RangePolicy<Kokkos::Cuda>(0,n),
        KOKKOS_LAMBDA(const int i) {
          view(i) *= factorDevice;
        });
    Kokkos::fence();
  }

  void modifyOnHost() {
    Kokkos::parallel_for(
        "modify on host",
        Kokkos::RangePolicy<Kokkos::Serial>(0,n),
        KOKKOS_LAMBDA(const int i) {
          view(i) *= factorHost;
        });
  }

private:
  bool isConsistent() {
    bool correct = true;
    for (int i=0; i<n && correct; i++) {
      correct &= (view(i) == factorCombined*i);
    }
    return correct;
  }

  void log(const std::string& msg) {
    std::cout << msg << std::endl;
  }

  const int n;
  const int factorDevice;
  const int factorHost;
  const int factorCombined;
  Kokkos::View<int*, Kokkos::CudaUVMSpace> view;
};


int main(int argc, char *argv[]) {
  Kokkos::ScopeGuard kokkos(argc, argv);

  UVMTester tester;
  tester.run();

  return 0;
}
