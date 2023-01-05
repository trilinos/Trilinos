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
#include <cstdio>
#include <mpi.h>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <ctime>
#include <chrono>
#include <Teuchos_StackedTimer.hpp>

#include "MueLu_PerfModels.hpp"

namespace MueLu {
  namespace PerfModels {

    int MTRX_MAX_SIZE = 800;
    int MTRX_MAX = pow(2,15);
    int MTRX_MIN = 1;


    std::vector<double> singleNodeDenseMatrixMultiplicationTest(int& KERNEL_REPEATS, int row1, int col1, int col2) {

      double **a, **b, **c;
      int i,j;
      std::vector<double> times(KERNEL_REPEATS); // Modify number of tests ran

      a = new double*[row1];
      b = new double*[col1];
      c = new double*[row1];

      // Create random seeds based on the smallest tick period of the node.
      long long t1 = std::chrono::high_resolution_clock::now().time_since_epoch().count();
      srand((size_t) t1); //updates the random number quicker then every second, used to rely on srand(time(nullptr)).

      // Populate the matrices with random values
        for(i = 0; i < row1; i++) {
          a[i] = new double[col1];
          c[i] = new double[col2];

          for(j = 0; j < col1; j++) {
            a[i][j] = rand() % MTRX_MAX + MTRX_MIN;
          }
        }

        for(i = 0; i < col1; i++) {
          b[i] = new double[col2];
          for(j = 0; j < col2; j++) {
            b[i][j] = rand() % MTRX_MAX + MTRX_MIN;
          }
        }

      // TIMED BENCHMARK BEGINS
      for(size_t t = 0; t < times.size(); t++) {
        clock_t start = clock();

        for(i = 0; i < row1; ++i) {
          for(j = 0; j < col2; ++j) {
            for(int k = 0; k < col1; ++k) {
                c[i][j] += a[i][k] * b[k][j];
            }
          }
        }

        clock_t end = clock();
        double diffs = (end - start)/(double)CLOCKS_PER_SEC;
        times[t] = diffs;
      }

  //    double avg = std::accumulate(times.begin(), times.end(), 0.0) / times.size();

      return times;


    }


    std::map<int,double> pingpong_test(int& KERNEL_REPEATS, int& MAX_SIZE, RCP<const Teuchos::Comm<int>> comm) {
      using Teuchos::BaseTimer;
      int rank, nproc, msg_length, i, j;
      std::vector<int> msg_arr((int)pow(2,MAX_SIZE)); // values from 0,1,2... to 2^15. Sizes of each buffer send
      char  *s_buf, *r_buf;  // Send & recieve buffers
      double t_avg;
      std::vector<double> time_array(KERNEL_REPEATS); //Stores the times for every single kernel repetition. Reset with each repeat.

      RCP<CommRequest<Ordinal>> request;
      RCP<CommStatus<Ordinal>> status;
      BaseTimer timer;

      std::map<int, double> time_map;

      rank = comm->getRank();
      nproc = comm->getSize();

      msg_arr[0] = 0;
      for(i = 0; i < MAX_SIZE; i++) {
        msg_arr[i+1] = (int) pow(2,i);
      }

      if (nproc < 2) {
        if (rank == 0) printf("This benchmark should be run on at least two processes");
        exit(EXIT_FAILURE);
      }

      //Allocating memory for the buffers.
      MPI_Alloc_mem( (int) pow(2,MAX_SIZE), MPI_INFO_NULL, &r_buf);
      MPI_Alloc_mem( (int) pow(2,MAX_SIZE), MPI_INFO_NULL, &s_buf);

      // Populate send buffer
      for(i = 0; i < (int) pow(2,MAX_SIZE); i++)
        s_buf[i] = 1;

      //Send and recieve.
      for(msg_length = 0; msg_length < MAX_SIZE + 1 ; msg_length++) {
        comm->barrier();

        for(j = 0; j < KERNEL_REPEATS; j++) {
          timer.start();

          if(rank == 1) {
            comm->send(msg_arr[msg_length], s_buf, 0);
          }

          else if(rank == 0){
            comm->receive(1, msg_arr[msg_length],r_buf);
          }

          timer.stop();
          time_array[j] = timer.accumulatedTime() * 1.0e6; // Formmated in microseconds (us)
          timer.reset();
        }

        t_avg = accumulate(time_array.begin(), time_array.end(), 0.0) / KERNEL_REPEATS;

        time_map.insert(std::pair<int, double>(msg_arr[msg_length],t_avg));

      }

    return time_map;
    }


  } //namespace PerfModels

} //namespace MueLu
