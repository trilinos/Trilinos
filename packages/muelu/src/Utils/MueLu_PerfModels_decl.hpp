// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_PERFMODELS_HPP
#define MUELU_PERFMODELS_HPP

#include <iostream>
#include <vector>
#include <Kokkos_Core.hpp>
#include <Kokkos_Vector.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_StackedTimer.hpp>


namespace Teuchos { class Time; }
namespace Teuchos { template <typename Ordinal> class Comm; }

namespace MueLu {

  namespace PerfModels {
    using namespace Teuchos;

    /* Single Node tests based upon UVA's STREAM benchmark for measuring memory
     * bandwith and computation rate. These processes compute either the addition
     * of two vectors or the multiplication of dense matrices of any given size.
     * Many iterations occur which then return a vector containing the individual
     * lengths of time per iteration.
     *
     * See further here:
     *    - https://www.cs.virginia.edu/stream/ref.html
     *    - https://github.com/UoB-HPC/BabelStream
     */
    template<class T>
    Kokkos::vector<double> singleNodeVectorAdditionTest(int& KERNEL_REPEATS, int& VECTOR_SIZE) {
      Kokkos::View<T*> a("a", VECTOR_SIZE);
      Kokkos::View<T*> b("b", VECTOR_SIZE);
      Kokkos::View<T*> c("c", VECTOR_SIZE);
      Kokkos::vector<double> test_times(KERNEL_REPEATS);


      Kokkos::parallel_for(VECTOR_SIZE, KOKKOS_LAMBDA (const size_t i) {
        a(i) = 1.0/(i+1);
        b(i) = a(i);
      });

      for(int i = 0; i < KERNEL_REPEATS; i++) {
        clock_t start = clock();

        Kokkos::parallel_for(VECTOR_SIZE, KOKKOS_LAMBDA (const size_t j) { //Vector Addition
            c(j) = a(j) + b(j);
        });

        Kokkos::fence();
        clock_t end = clock();
        double diffs = (end - start)/(double)CLOCKS_PER_SEC;
        test_times[i] = diffs;
      }

      return test_times;
    }


    /* A latency test between two processes based upon the MVAPICH OSU Micro-Benchmarks.
     * The sender process sends a message and then waits for confirmation of reception.
     * Many iterations occur with various message sizes and the average latency values
     * are returned within a map. Utilizes blocking send and recieve.
     *
     * See further: https://mvapich.cse.ohio-state.edu/benchmarks/
     */
    std::map<int,double> pingpong_test(int& KERNEL_REPEATS, int& MAX_SIZE, RCP<const Teuchos::Comm<int>> comm);

  } //namespace PerfModels

} //namespace MueLu

#endif //ifndef MUELU_PERFMODELS_HPP
