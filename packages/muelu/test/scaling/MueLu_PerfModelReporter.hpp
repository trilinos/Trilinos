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
#ifndef MUELU_PERFMODEL_REPORTER_HPP
#define MUELU_PERFMODEL_REPORTER_HPP
#include <vector>

#include <Teuchos_RCP.hpp>
#include <Xpetra_Matrix.hpp>

#include "MueLu_PerfModels.hpp"

namespace MueLu {

// =========================================================================
// Performance Routines
// =========================================================================
// Report bandwidth in GB / sec

/* @brief Generate performance model reports for SPMV of a given matrix
   Input:
     A       - Matrix to evaluate
     nrepeat - # times to repeat the tests in model construction
     timer_names - Names of Teuchos timers which are SPMVs to compare agains the model
     verbose - Enable verbose output
*/
template <class Matrix>
void report_spmv_performance_models(const Teuchos::RCP<const Matrix> &A, int nrepeat, const std::vector<const char *> &timer_names, Teuchos::RCP<Teuchos::TimeMonitor> &myTimeMonitor, const std::string prefix = "", bool verbose = false) {
  using Teuchos::RCP;
  using Teuchos::rcp;
  const RCP<const Teuchos::Comm<int> > comm = A->getMap()->getComm();
  using SC                                  = typename Matrix::scalar_type;
  using LO                                  = typename Matrix::local_ordinal_type;
  using GO                                  = typename Matrix::global_ordinal_type;
  using NO                                  = typename Matrix::node_type;

  const double GB = 1024.0 * 1024.0 * 1024.0;

  // Conversion function as a lambda
  auto convert_time_to_bandwidth_gbs = [=](double time, int num_calls, double memory_per_call_bytes) {
    double time_per_call = time / num_calls;

    return memory_per_call_bytes / GB / time_per_call;
  };

  // NOTE: We've hardwired this to size_t for the rowptr.  This really should really get read out of a typedef,
  // if Tpetra actually had one
  using rowptr_type = size_t;

  // Make a new model if we need one
  MueLu::PerfModels<SC, LO, GO, NO> PM;

  int rank  = comm->getRank();
  int nproc = comm->getSize();
  int m     = static_cast<int>(A->getLocalNumRows());
  int n     = static_cast<int>(A->getColMap()->getLocalNumElements());
  int nnz   = static_cast<int>(A->getLocalMatrixHost().graph.entries.extent(0));

  // Generate Lookup Tables
  int v_log_max = ceil(log(nnz) / log(2)) + 1;
  PM.stream_vector_make_table(nrepeat, v_log_max);

  int m_log_max = 15;
  PM.pingpong_make_table(nrepeat, m_log_max, comm);

  std::string rank_prefix = prefix + std::to_string(rank) + std::string(": ");

  if (A->hasCrsGraph()) {
    auto importer = A->getCrsGraph()->getImporter();
    if (!importer.is_null()) {
      size_t recv_size   = importer->getRemoteLIDs().size() * sizeof(SC);
      size_t send_size   = importer->getExportLIDs().size() * sizeof(SC);
      int local_log_max  = ceil(log(std::max(send_size, recv_size)) / log(2)) + 1;
      int global_log_max = local_log_max;
      Teuchos::reduceAll<int>(*comm, Teuchos::REDUCE_MAX, 1, &local_log_max, &global_log_max);
      PM.halopong_make_table(nrepeat, global_log_max, importer);
    }
  }

  if (verbose) {
    std::cout << prefix << "********************************************************" << std::endl;
    std::cout << prefix << "Performance model results on " << nproc << " ranks" << std::endl;
    if (rank == 0) {
      std::cout << rank_prefix << " ****** Launch Latency Table ******" << std::endl;
      PM.print_launch_latency_table(std::cout, rank_prefix);
      std::cout << rank_prefix << "****** Stream Table ******" << std::endl;
      PM.print_stream_vector_table(std::cout, rank_prefix);
      std::cout << rank_prefix << "****** Latency Corrected Stream Table ******" << std::endl;
      PM.print_latency_corrected_stream_vector_table(std::cout, rank_prefix);
      std::cout << rank_prefix << "****** Pingpong Table ******" << std::endl;
      PM.print_pingpong_table(std::cout, rank_prefix);
      std::cout << rank_prefix << "****** Halopong Table ******" << std::endl;
      PM.print_halopong_table(std::cout, rank_prefix);
    }
  }

  // For convenience
  const int NUM_TIMERS                    = 6;
  std::string SPMV_test_names[NUM_TIMERS] = {"colind", "rowptr", "vals", "x", "y", "all"};
  std::vector<int> SPMV_num_objects(NUM_TIMERS), SPMV_object_size(NUM_TIMERS), SPMV_corrected(NUM_TIMERS);

  // Composite model: Use latency correction
  SPMV_num_objects[0] = nnz;
  SPMV_object_size[0] = sizeof(LO);
  SPMV_corrected[0]   = 1;  // colind
  SPMV_num_objects[1] = (m + 1);
  SPMV_object_size[1] = sizeof(rowptr_type);
  SPMV_corrected[1]   = 1;  // rowptr
  SPMV_num_objects[2] = nnz;
  SPMV_object_size[2] = sizeof(SC);
  SPMV_corrected[2]   = 1;  // vals
  SPMV_num_objects[3] = n;
  SPMV_object_size[3] = sizeof(SC);
  SPMV_corrected[3]   = 1;  // x
  SPMV_num_objects[4] = m;
  SPMV_object_size[4] = sizeof(SC);
  SPMV_corrected[4]   = 1;  // y

  // All-Model: Do not use latency correction
  SPMV_object_size[5] = 1;
  SPMV_num_objects[5] = (m + 1) * sizeof(rowptr_type) + nnz * sizeof(LO) + nnz * sizeof(SC) +
                        n * sizeof(SC) + m * sizeof(SC);
  SPMV_corrected[5] = 0;

  std::vector<double> gb_per_sec(NUM_TIMERS);
  if (verbose && rank == 0)
    std::cout << rank_prefix << "****** Local Time Model Results ******" << std::endl;
  for (int i = 0; i < NUM_TIMERS; i++) {
    double avg_time;

    // Table interpolation - Take the faster of the two
    int size_in_bytes = SPMV_object_size[i] * SPMV_num_objects[i];
    if (SPMV_corrected[i] == 1)
      avg_time = PM.latency_corrected_stream_vector_lookup(size_in_bytes);
    else
      avg_time = PM.stream_vector_lookup(size_in_bytes);

    // The lookup divides by transactions-per-element already
    double memory_traffic = (double)SPMV_object_size[i] * (double)SPMV_num_objects[i];
    gb_per_sec[i]         = convert_time_to_bandwidth_gbs(avg_time, 1, memory_traffic);

    if (verbose) {
      std::cout << rank_prefix << "Local: " << SPMV_test_names[i] << " # Scalars = " << memory_traffic / sizeof(SC) << " time per call = " << avg_time * 1e6 << " us. GB/sec = " << gb_per_sec[i] << std::endl;
    }
  }

  // Get the latency info
  double avg_latency = PM.launch_latency_lookup();

  /***************************************************************************/
  // *** Calculate SPMV minimum time (composite) ***
  // Model:
  // rowptr = One read per row
  // colind = One read per entry
  // values = One read per entry
  // x      = One read per entry in values array (but we assume the cache will work its magic here)
  // y      = One write per row

  long unsigned int spmv_memory_bytes[NUM_TIMERS] = {
      (m + 1) * sizeof(rowptr_type),  // rowptr
      nnz * sizeof(LO),               // colind
      nnz * sizeof(SC),               // values
      n * sizeof(SC),                 // x
      m * sizeof(SC),                 // y
      0};
  for (int i = 0; i < NUM_TIMERS - 1; i++)
    spmv_memory_bytes[NUM_TIMERS - 1] += spmv_memory_bytes[i];

  double minimum_local_composite_time = avg_latency;
  for (int i = 0; i < NUM_TIMERS - 1; i++)
    minimum_local_composite_time += spmv_memory_bytes[i] / (GB * gb_per_sec[i]);

  double minimum_local_all_time = spmv_memory_bytes[NUM_TIMERS - 1] / (GB * gb_per_sec[NUM_TIMERS - 1]);

  if (verbose) {
    std::cout << rank_prefix << "Local: composite      = " << minimum_local_composite_time * 1e6 << " us.\n"
              << rank_prefix << "Local: all            = " << minimum_local_all_time * 1e6 << " us.\n";
  }

  /***************************************************************************/
  // *** Calculate Remote part of the SPMV ***
  double time_pack_unpack_outofplace = 0.0;
  double time_pack_unpack_inplace    = 0.0;
  double time_communicate_ping       = 0.0;
  double time_communicate_halo       = 0.0;
  // Note: We'll assume that each of the permutes, remotes and exports is a unified
  // memory transaction, even though that's not strictly speaking correct.
  if (A->hasCrsGraph()) {
    auto importer = A->getCrsGraph()->getImporter();
    if (!importer.is_null()) {
      // Sames [pack] - 1 read SC, 1 write SC
      // NOTE: Only if you're out-of-place
      size_t num_sames = importer->getNumSameIDs();
      double same_time = (num_sames == 0) ? 0.0 : 2.0 * PM.latency_corrected_stream_vector_lookup(num_sames * sizeof(SC)) + avg_latency;

      // Permutes [pack] - 2 reads LO [to, from] , 1 read SC [values], 1 write SC [values]
      size_t num_permutes = importer->getNumPermuteIDs();
      double permute_time = (num_permutes == 0) ? 0.0 : 2.0 * PM.latency_corrected_stream_vector_lookup(num_permutes * sizeof(LO)) + 2.0 * PM.latency_corrected_stream_vector_lookup(num_permutes * sizeof(SC)) + avg_latency;

      // Exports [pack] - 1 read LO [exportLIDs], 1 read SC [values], 1 write SC [buffer]
      // This is what Epetra does at least
      size_t num_exports = importer->getNumExportIDs();
      double export_time = (num_exports == 0) ? 0.0 : PM.latency_corrected_stream_vector_lookup(num_exports * sizeof(LO)) + 2.0 * PM.latency_corrected_stream_vector_lookup(num_exports * sizeof(SC)) + avg_latency;

      // Remotes [unpack] - 1 read LO [remoteLIDs],  1 read SC [buffer], 1 write SC [values]
      // NOTE: Only if you're out of place
      size_t num_remotes = importer->getNumRemoteIDs();
      double remote_time = (num_remotes == 0) ? 0.0 : PM.latency_corrected_stream_vector_lookup(num_remotes * sizeof(LO)) + 2.0 * PM.latency_corrected_stream_vector_lookup(num_remotes * sizeof(SC)) + avg_latency;

      // Total pack / unpack time
      time_pack_unpack_outofplace = same_time + permute_time + export_time + remote_time;
      time_pack_unpack_inplace    = permute_time + export_time;

      // We now need to get the size of each message for the ping-pong costs.
      double send_time         = 0.0;
      double recv_time         = 0.0;
      double halo_time         = 0.0;
      size_t total_send_length = 0, total_recv_length = 0;
      double avg_size_per_msg                                 = 0.0;
      RCP<const Xpetra::TpetraImport<LO, GO, NO> > t_importer = Teuchos::rcp_dynamic_cast<const Xpetra::TpetraImport<LO, GO, NO> >(importer);
      if (!t_importer.is_null()) {
        RCP<const Tpetra::Import<LO, GO, NO> > tt_i   = t_importer->getTpetra_Import();
        Tpetra::Distributor &distor                   = tt_i->getDistributor();
        Teuchos::ArrayView<const size_t> recv_lengths = distor.getLengthsFrom();
        Teuchos::ArrayView<const size_t> send_lengths = distor.getLengthsTo();

        for (int i = 0; i < (int)send_lengths.size(); i++) {
          send_time += PM.pingpong_device_lookup(send_lengths[i] * sizeof(SC));
          total_send_length = send_lengths[i] * sizeof(SC);
        }

        for (int i = 0; i < (int)recv_lengths.size(); i++) {
          recv_time += PM.pingpong_device_lookup(recv_lengths[i] * sizeof(SC));
          total_recv_length = recv_lengths[i] * sizeof(SC);
        }

        avg_size_per_msg = (double)total_send_length / (2.0 * send_lengths.size()) + (double)total_recv_length / (2.0 * recv_lengths.size());
        halo_time        = PM.halopong_device_lookup(avg_size_per_msg);
      }

      // NOTE: For now we'll do comm time as the larger of send/recv.  Not sure this is
      // really the optimal thing to do, but we'll start here.
      time_communicate_ping = std::max(send_time, recv_time);
      time_communicate_halo = halo_time;

      if (verbose) {
        std::cout << rank_prefix << "****** Remote Time Model Results ******" << std::endl;
        std::cout << rank_prefix << "Remote: same      = " << same_time * 1e6 << " us.\n"
                  << rank_prefix << "Remote: permutes  = " << permute_time * 1e6 << " us.\n"
                  << rank_prefix << "Remote: exports   = " << export_time * 1e6 << " us.\n"
                  << rank_prefix << "Remote: remotes   = " << remote_time * 1e6 << " us.\n"
                  << rank_prefix << "Remote: sends len = " << total_send_length << " time = " << send_time * 1e6 << " us.\n"
                  << rank_prefix << "Remote: recvs len = " << total_recv_length << " time = " << recv_time * 1e6 << " us.\n"
                  << rank_prefix << "Remote: pack/unpack OOP = " << time_pack_unpack_outofplace * 1e6 << " us.\n"
                  << rank_prefix << "Remote: pack/unpack INP = " << time_pack_unpack_inplace * 1e6 << " us.\n"
                  << rank_prefix << "Remote: ping avg  = " << (size_t)avg_size_per_msg << " time  = " << time_communicate_ping * 1e6 << " us.\n"
                  << rank_prefix << "Remote: halo avg  = " << (size_t)avg_size_per_msg << " time  = " << time_communicate_halo * 1e6 << " us.\n"
                  << std::endl;
      }
    }
  }
  double minimum_time_in_place_ping     = time_communicate_ping + time_pack_unpack_inplace;
  double minimum_time_out_of_place_ping = time_communicate_ping + time_pack_unpack_outofplace;
  double minimum_time_in_place_halo     = time_communicate_halo + time_pack_unpack_inplace;
  double minimum_time_out_of_place_halo = time_communicate_halo + time_pack_unpack_outofplace;

  /***************************************************************************/

  // Get global average/max time sums
  constexpr int NUM_TIMES = 12;
  double avg_times[NUM_TIMES];
  double max_times[NUM_TIMES];
  {
    double comp      = minimum_local_composite_time;
    double alls      = minimum_local_all_time;
    double p_inplace = minimum_time_in_place_ping;
    double p_ooplace = minimum_time_out_of_place_ping;
    double h_inplace = minimum_time_in_place_halo;
    double h_ooplace = minimum_time_out_of_place_halo;

    double all_times_local[NUM_TIMES] = {
        comp, comp + p_inplace, comp + p_ooplace, comp + h_inplace, comp + h_ooplace,  // 0-4
        alls, alls + p_inplace, alls + p_ooplace, alls + h_inplace, alls + h_ooplace,  // 5-9
        time_communicate_ping, time_communicate_halo};                                 // 10,11

    if (nproc > 1) {
      Teuchos::reduceAll(*comm, Teuchos::REDUCE_SUM, NUM_TIMES, &all_times_local[0], &avg_times[0]);
      Teuchos::reduceAll(*comm, Teuchos::REDUCE_MAX, NUM_TIMES, &all_times_local[0], &max_times[0]);
      for (int i = 0; i < NUM_TIMES; i++)
        avg_times[i] /= nproc;
    } else {
      for (int i = 0; i < NUM_TIMES; i++)
        avg_times[i] = max_times[i] = all_times_local[i];
    }

    if (rank == 0)
      std::cout << "\n\n========================================================\n"
                << prefix << "Minimum time model (composite) : " << avg_times[0] * 1e6 << " us.\n"
                << prefix << "Minimum time model (all)       : " << avg_times[5] * 1e6 << " us.\n"
                << prefix << "Pack/unpack in-place           : " << (avg_times[1] - avg_times[0] - avg_times[10]) * 1e6 << "us.\n"  // comp+ping_pack_inplace - comp - ping
                << prefix << "Pack/unpack out-of-place       : " << (avg_times[2] - avg_times[0] - avg_times[10]) * 1e6 << "us.\n"  // comp+ping_pack_ooplace - comp - ping
                << prefix << "Communication time (ping)      : " << avg_times[10] * 1e6 << "us.\n"
                << prefix << "Communication time (halo)      : " << avg_times[11] * 1e6 << "us." << std::endl;
  }

  // Iterate through all of the "MV" timers

  if (rank == 0) {
    if (!myTimeMonitor.is_null() && timer_names.size() > 0) {
      const std::string l[NUM_TIMES] = {"Comp", "Comp+ping+inplace", "Comp+ping+ooplace", "Comp+halo+inplace", "Comp+halo+ooplace",
                                        "All", "All+ping+inplace", "All+ping+ooplace", "All+halo+inplace", "All+halo+ooplace"};
      const std::string div          = {"-------------------"};
      printf("%-60s %20s %20s %20s %20s %20s %20s %20s %20s %20s %20s\n", "Timer", l[0].c_str(), l[1].c_str(), l[2].c_str(), l[3].c_str(), l[4].c_str(),
             l[5].c_str(), l[6].c_str(), l[7].c_str(), l[8].c_str(), l[9].c_str());
      printf("%-60s %20s %20s %20s %20s %20s %20s %20s %20s %20s %20s\n", "-----", div.c_str(), div.c_str(), div.c_str(), div.c_str(), div.c_str(),
             div.c_str(), div.c_str(), div.c_str(), div.c_str(), div.c_str());
      for (int i = 0; i < (int)timer_names.size(); i++) {
        Teuchos::RCP<Teuchos::Time> t = myTimeMonitor->lookupCounter(timer_names[i]);
        if (!t.is_null()) {
          double time_per_call = t->totalElapsedTime() / t->numCalls();
          printf("%-60s %20.2f %20.2f %20.2f %20.2f %20.2f %20.2f %20.2f %20.2f %20.2f %20.2f\n", timer_names[i],
                 max_times[0] / time_per_call, max_times[1] / time_per_call, max_times[2] / time_per_call, max_times[3] / time_per_call, max_times[4] / time_per_call,
                 max_times[5] / time_per_call, max_times[6] / time_per_call, max_times[7] / time_per_call, max_times[8] / time_per_call, max_times[9] / time_per_call);
        }
      }
    } else {
      std::cout << prefix << "Note: Minimum time model individual timers only work with stacked timers off." << std::endl;
    }
  }
}

}  // namespace MueLu

#endif
