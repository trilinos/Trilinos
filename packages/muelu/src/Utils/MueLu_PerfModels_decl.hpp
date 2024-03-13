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

#include "MueLu_ConfigDefs.hpp"
#include "Xpetra_Import_fwd.hpp"

#include <vector>
#include <ostream>
#include <Teuchos_DefaultComm.hpp>

#include "MueLu_PerfModels_fwd.hpp"

namespace MueLu {

template <class Scalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class PerfModels {
 public:
  PerfModels();

  /* Single Node tests based upon the STREAM benchmark for measuring memory
   * bandwith and computation rate. These processes compute either the addition
   * of two vectors or the multiplication of dense matrices of any given size.
   * Many iterations occur which then return a vector containing the individual
   * lengths of time per iteration.
   *
   * See further here:
   *    - https://www.cs.virginia.edu/stream/ref.html
   *    - https://github.com/UoB-HPC/BabelStream
   */

  /* This version is for table interpolation and works on chars, so the LOG_MAX_SIZE is for bytes */
  void stream_vector_make_table(int KERNEL_REPEATS, int LOG_MAX_SIZE = 20);
  bool has_stream_vector_table() const { return stream_sizes_.size() > 0; }

  /* Lookup in the stream_vector table */
  double stream_vector_copy_lookup(int SIZE_IN_BYTES);
  double stream_vector_add_lookup(int SIZE_IN_BYTES);
  double latency_corrected_stream_vector_copy_lookup(int SIZE_IN_BYTES);
  double latency_corrected_stream_vector_add_lookup(int SIZE_IN_BYTES);

  // Uses the faster of the tables.  The time is then divided by the number of memory transactions
  // per element in the kernel (e.g. 2 for COPY and 3 for ADD).
  double stream_vector_lookup(int SIZE_IN_BYTES);
  double latency_corrected_stream_vector_lookup(int SIZE_IN_BYTES);

  /* Print table */
  void print_stream_vector_table(std::ostream &out, const std::string &prefix = "");
  void print_latency_corrected_stream_vector_table(std::ostream &out, const std::string &prefix = "");

  /* A latency test between two processes based upon the MVAPICH OSU Micro-Benchmarks.
   * The sender process sends a message and then waits for confirmation of reception.
   * Many iterations occur with various message sizes and the average latency values
   * are returned within a map. Utilizes blocking send and recieve.
   *
   * See further: https://mvapich.cse.ohio-state.edu/benchmarks/
   */
  void pingpong_make_table(int KERNEL_REPEATS, int LOG_MAX_SIZE, const RCP<const Teuchos::Comm<int> > &comm);
  bool has_pingpong_table() const { return pingpong_sizes_.size() > 0; }

  /* Lookup in the pingpong_vector table */
  double pingpong_host_lookup(int SIZE_IN_BYTES);
  double pingpong_device_lookup(int SIZE_IN_BYTES);

  /* Print table */
  void print_pingpong_table(std::ostream &out, const std::string &prefix = "");

  /* A halo-exchange based ping-pong, inspired by halo-mode in MPPTEST from ANL.
   * Here we use exactly the communication pattern specified in the import object
   * and send messages accordingly.  We vary the size in bytes sent per message,
   * which should capture max-rate effects to some degree.
   *
   * See further: https://www.mcs.anl.gov/research/projects/mpi/mpptest/
   */
  void halopong_make_table(int KERNEL_REPEATS, int LOG_MAX_SIZE, const RCP<const Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > &import);
  bool has_halopong_table() const { return halopong_sizes_.size() > 0; }

  /* Lookup in the halopong_vector table */
  double halopong_host_lookup(int SIZE_IN_BYTES_PER_MESSAGE);
  double halopong_device_lookup(int SIZE_IN_BYTES_PER_MESSAGE);

  /* Print table */
  void print_halopong_table(std::ostream &out, const std::string &prefix = "");

  /* Estimate launch latency based on the cost of submitting an empty Kokkos::parallel_for.
   * This necessary to correct the memory bandwidth costs for models on high latency platforms,
   * e.g., GPUS.
   */
  void launch_latency_make_table(int KERNEL_REPEATS);
  bool has_launch_latency_table() const { return launch_and_wait_latency_ > 0; }

  /* Lookup launch latency */
  double launch_latency_lookup();

  /* Print table */
  void print_launch_latency_table(std::ostream &out, const std::string &prefix = "");

 private:
  void print_stream_vector_table_impl(std::ostream &out, bool use_latency_correction, const std::string &prefix);

  std::vector<int> stream_sizes_;
  std::vector<double> stream_copy_times_;
  std::vector<double> stream_add_times_;
  std::vector<double> latency_corrected_stream_copy_times_;
  std::vector<double> latency_corrected_stream_add_times_;

  std::vector<int> pingpong_sizes_;
  std::vector<double> pingpong_host_times_;
  std::vector<double> pingpong_device_times_;

  std::vector<int> halopong_sizes_;
  std::vector<double> halopong_host_times_;
  std::vector<double> halopong_device_times_;

  double launch_and_wait_latency_;

};  // class PerfModels

}  // namespace MueLu

#endif  // ifndef MUELU_PERFMODELS_HPP
