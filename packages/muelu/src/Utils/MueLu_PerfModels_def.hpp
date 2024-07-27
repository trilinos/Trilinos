// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "MueLu_PerfModels_decl.hpp"

#include <cstdio>
#include <cmath>
#include <numeric>
#include <utility>
#include <chrono>
#include <iomanip>
#include <Teuchos_ScalarTraits.hpp>
#include <Kokkos_ArithTraits.hpp>
#include <Xpetra_Import.hpp>
#if defined(HAVE_MUELU_TPETRA) && defined(HAVE_MPI)
#include <Xpetra_TpetraImport.hpp>
#include <Tpetra_Import.hpp>
#include <Tpetra_Distributor.hpp>
#include <mpi.h>
#endif

#ifdef HAVE_MPI
#include <mpi.h>
#endif

namespace MueLu {

namespace PerfDetails {
template <class Scalar, class Node>
double stream_vector_add(int KERNEL_REPEATS, int VECTOR_SIZE) {
  // PerfDetails' STREAM routines need to be instantiatiated on impl_scalar_type, not Scalar
  using impl_scalar_type = typename Kokkos::ArithTraits<Scalar>::val_type;

  using exec_space   = typename Node::execution_space;
  using memory_space = typename Node::memory_space;
  using range_policy = Kokkos::RangePolicy<exec_space>;

  Kokkos::View<impl_scalar_type *, memory_space> a("a", VECTOR_SIZE);
  Kokkos::View<impl_scalar_type *, memory_space> b("b", VECTOR_SIZE);
  Kokkos::View<impl_scalar_type *, memory_space> c("c", VECTOR_SIZE);
  double total_test_time = 0.0;

  impl_scalar_type ONE = Teuchos::ScalarTraits<impl_scalar_type>::one();

  Kokkos::parallel_for(
      "stream/fill", range_policy(0, VECTOR_SIZE), KOKKOS_LAMBDA(const size_t i) {
        a(i) = ONE * (double)i;
        b(i) = a(i);
      });
  exec_space().fence();

  using clock = std::chrono::high_resolution_clock;

  clock::time_point start, stop;

  for (int i = 0; i < KERNEL_REPEATS; i++) {
    start = clock::now();
    Kokkos::parallel_for(
        "stream/add", range_policy(0, VECTOR_SIZE), KOKKOS_LAMBDA(const size_t j) {  // Vector Addition
          c(j) = a(j) + b(j);
        });

    exec_space().fence();
    stop                = clock::now();
    double my_test_time = std::chrono::duration<double>(stop - start).count();
    total_test_time += my_test_time;
  }

  return total_test_time / KERNEL_REPEATS;
}

template <class Scalar, class Node>
double stream_vector_copy(int KERNEL_REPEATS, int VECTOR_SIZE) {
  // PerfDetails' STREAM routines need to be instantiatiated on impl_scalar_type, not Scalar
  using impl_scalar_type = typename Kokkos::ArithTraits<Scalar>::val_type;

  using exec_space   = typename Node::execution_space;
  using memory_space = typename Node::memory_space;
  using range_policy = Kokkos::RangePolicy<exec_space>;

  Kokkos::View<impl_scalar_type *, memory_space> a("a", VECTOR_SIZE);
  Kokkos::View<impl_scalar_type *, memory_space> b("b", VECTOR_SIZE);
  double total_test_time = 0.0;

  impl_scalar_type ONE = Teuchos::ScalarTraits<impl_scalar_type>::one();

  Kokkos::parallel_for(
      "stream/fill", range_policy(0, VECTOR_SIZE), KOKKOS_LAMBDA(const size_t i) {
        a(i) = ONE;
      });
  exec_space().fence();

  using clock = std::chrono::high_resolution_clock;
  clock::time_point start, stop;

  for (int i = 0; i < KERNEL_REPEATS; i++) {
    start = clock::now();
    Kokkos::parallel_for(
        "stream/copy", range_policy(0, VECTOR_SIZE), KOKKOS_LAMBDA(const size_t j) {  // Vector Addition
          b(j) = a(j);
        });

    exec_space().fence();
    stop                = clock::now();
    double my_test_time = std::chrono::duration<double>(stop - start).count();
    total_test_time += my_test_time;
  }

  return total_test_time / KERNEL_REPEATS;
}

double table_lookup(const std::vector<int> &x, const std::vector<double> &y, int value) {
  // If there's no table, nan
  if (x.size() == 0) return Teuchos::ScalarTraits<double>::nan();

  // NOTE:  This should probably be a binary search, but this isn't performance sensitive, so we'll go simple
  int N  = (int)x.size();
  int hi = 0;
  for (; hi < N; hi++) {
    if (x[hi] > value)
      break;
  }

  if (hi == 0) {
    // Lower end (return the min time)
    // printf("Lower end: %d < %d\n",value,x[0]);
    return y[0];
  } else if (hi == N) {
    // Higher end (extrapolate from the last two points)
    // printf("Upper end: %d > %d\n",value,x[N-1]);
    hi           = N - 1;
    int run      = x[hi] - x[hi - 1];
    double rise  = y[hi] - y[hi - 1];
    double slope = rise / run;
    int diff     = value - x[hi - 1];

    return y[hi - 1] + slope * diff;
  } else {
    // Interpolate
    // printf("Middle: %d < %d < %d\n",x[hi-1],value,x[hi]);
    int run      = x[hi] - x[hi - 1];
    double rise  = y[hi] - y[hi - 1];
    double slope = rise / run;
    int diff     = value - x[hi - 1];

    return y[hi - 1] + slope * diff;
  }
}

// Report bandwidth in GB / sec
const double GB = 1024.0 * 1024.0 * 1024.0;
double convert_time_to_bandwidth_gbs(double time, int num_calls, double memory_per_call_bytes) {
  double time_per_call = time / num_calls;
  return memory_per_call_bytes / GB / time_per_call;
}

template <class exec_space, class memory_space>
void pingpong_basic(int KERNEL_REPEATS, int MAX_SIZE, const Teuchos::Comm<int> &comm, std::vector<int> &sizes, std::vector<double> &times) {
#ifdef HAVE_MPI
  int rank  = comm.getRank();
  int nproc = comm.getSize();

  if (nproc < 2) return;

  const int buff_size = (int)pow(2, MAX_SIZE);

  sizes.resize(MAX_SIZE + 1);
  times.resize(MAX_SIZE + 1);

  // Allocate memory for the buffers (and fill send)
  Kokkos::View<char *, memory_space> r_buf("recv", buff_size), s_buf("send", buff_size);
  Kokkos::deep_copy(s_buf, 1);

  // Send and recieve.
  //  NOTE:  Do consectutive pair buddies here for simplicity.  We should be smart later
  int odd   = rank % 2;
  int buddy = odd ? rank - 1 : rank + 1;

  for (int i = 0; i < MAX_SIZE + 1; i++) {
    int msg_size = (int)pow(2, i);
    comm.barrier();

    double t0 = MPI_Wtime();
    for (int j = 0; j < KERNEL_REPEATS; j++) {
      if (buddy < nproc) {
        if (odd) {
          comm.send(msg_size, (char *)s_buf.data(), buddy);
          comm.receive(buddy, msg_size, (char *)r_buf.data());
        } else {
          comm.receive(buddy, msg_size, (char *)r_buf.data());
          comm.send(msg_size, (char *)s_buf.data(), buddy);
        }
      }
    }

    double time_per_call = (MPI_Wtime() - t0) / (2.0 * KERNEL_REPEATS);
    sizes[i]             = msg_size;
    times[i]             = time_per_call;
  }
#else
  return;
#endif
}

template <class exec_space, class memory_space, class LocalOrdinal, class GlobalOrdinal, class Node>
void halopong_basic(int KERNEL_REPEATS, int MAX_SIZE, const RCP<const Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > &import, std::vector<int> &sizes, std::vector<double> &times) {
  int nproc = import->getSourceMap()->getComm()->getSize();
  if (nproc < 2) return;
#if defined(HAVE_MUELU_TPETRA) && defined(HAVE_MPI)
  // NOTE: We need to get the distributer here, which means we need Tpetra, since Xpetra does
  // not have a distributor interface
  using x_import_type                     = Xpetra::TpetraImport<LocalOrdinal, GlobalOrdinal, Node>;
  RCP<const x_import_type> Ximport        = Teuchos::rcp_dynamic_cast<const x_import_type>(import);
  RCP<const Teuchos::MpiComm<int> > mcomm = Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >(import->getSourceMap()->getComm());
  MPI_Comm communicator                   = *mcomm->getRawMpiComm();

  if (Ximport.is_null() || mcomm.is_null()) return;
  auto Timport = Ximport->getTpetra_Import();
  auto distor  = Timport->getDistributor();

  // Distributor innards
  Teuchos::ArrayView<const int> procsFrom = distor.getProcsFrom();
  Teuchos::ArrayView<const int> procsTo   = distor.getProcsTo();
  int num_recvs                           = (int)distor.getNumReceives();
  int num_sends                           = (int)distor.getNumSends();

  const int buff_size_per_msg = (int)pow(2, MAX_SIZE);
  sizes.resize(MAX_SIZE + 1);
  times.resize(MAX_SIZE + 1);

  // Allocate memory for the buffers (and fill send)
  Kokkos::View<char *, memory_space> f_recv_buf("forward_recv", buff_size_per_msg * num_recvs), f_send_buf("forward_send", buff_size_per_msg * num_sends);
  Kokkos::View<char *, memory_space> r_recv_buf("reverse_recv", buff_size_per_msg * num_sends), r_send_buf("reverse_send", buff_size_per_msg * num_recvs);
  Kokkos::deep_copy(f_send_buf, 1);
  Kokkos::deep_copy(r_send_buf, 1);

  std::vector<MPI_Request> requests(num_sends + num_recvs);
  std::vector<MPI_Status> status(num_sends + num_recvs);

  for (int i = 0; i < MAX_SIZE + 1; i++) {
    int msg_size = (int)pow(2, i);

    MPI_Barrier(communicator);

    double t0 = MPI_Wtime();
    for (int j = 0; j < KERNEL_REPEATS; j++) {
      int ct = 0;
      // Recv/Send the forward messsages
      for (int r = 0; r < num_recvs; r++) {
        const int tag = 1000 + j;
        MPI_Irecv(f_recv_buf.data() + msg_size * r, msg_size, MPI_CHAR, procsFrom[r], tag, communicator, &requests[ct]);
        ct++;
      }
      for (int s = 0; s < num_sends; s++) {
        const int tag = 1000 + j;
        MPI_Isend(f_send_buf.data() + msg_size * s, msg_size, MPI_CHAR, procsTo[s], tag, communicator, &requests[ct]);
        ct++;
      }
      // Wait for the forward messsages
      MPI_Waitall(ct, requests.data(), status.data());

      ct = 0;
      // Recv/Send the reverse messsages
      for (int r = 0; r < num_sends; r++) {
        const int tag = 2000 + j;
        MPI_Irecv(r_recv_buf.data() + msg_size * r, msg_size, MPI_CHAR, procsTo[r], tag, communicator, &requests[ct]);
        ct++;
      }
      for (int s = 0; s < num_recvs; s++) {
        const int tag = 2000 + j;
        MPI_Isend(r_send_buf.data() + msg_size * s, msg_size, MPI_CHAR, procsFrom[s], tag, communicator, &requests[ct]);
        ct++;
      }
      // Wait for the reverse messsages
      MPI_Waitall(ct, requests.data(), status.data());
    }

    double time_per_call = (MPI_Wtime() - t0) / (2.0 * KERNEL_REPEATS);
    sizes[i]             = msg_size;
    times[i]             = time_per_call;
  }

#endif
}

}  // end namespace PerfDetails

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
PerfModels<Scalar, LocalOrdinal, GlobalOrdinal, Node>::PerfModels()
  : launch_and_wait_latency_(-1.0) {}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
PerfModels<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~PerfModels() {}

/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void PerfModels<Scalar, LocalOrdinal, GlobalOrdinal, Node>::stream_vector_make_table(int KERNEL_REPEATS, int LOG_MAX_SIZE) {
  // We need launch/waits latency estimates for corrected stream
  launch_latency_make_table(KERNEL_REPEATS);
  double latency = launch_latency_lookup();

  if (LOG_MAX_SIZE < 2)
    LOG_MAX_SIZE = 20;

  stream_sizes_.resize(LOG_MAX_SIZE + 1);
  stream_copy_times_.resize(LOG_MAX_SIZE + 1);
  stream_add_times_.resize(LOG_MAX_SIZE + 1);
  latency_corrected_stream_copy_times_.resize(LOG_MAX_SIZE + 1);
  latency_corrected_stream_add_times_.resize(LOG_MAX_SIZE + 1);

  for (int i = 0; i < LOG_MAX_SIZE + 1; i++) {
    int size      = (int)pow(2, i);
    double c_time = PerfDetails::stream_vector_copy<Scalar, Node>(KERNEL_REPEATS, size);
    double a_time = PerfDetails::stream_vector_add<Scalar, Node>(KERNEL_REPEATS, size);

    stream_sizes_[i] = size;

    // Correct for the difference in memory transactions per element
    stream_copy_times_[i] = c_time / 2.0;
    stream_add_times_[i]  = a_time / 3.0;

    // Correct for launch latency too.  We'll note that sometimes the latency estimate
    // is higher than the actual copy/add time estimate.  If so, we don't correct
    latency_corrected_stream_copy_times_[i] = (c_time - latency <= 0.0) ? c_time / 2.0 : ((c_time - latency) / 2.0);
    latency_corrected_stream_add_times_[i]  = (a_time - latency <= 0.0) ? a_time / 3.0 : ((a_time - latency) / 3.0);
  }
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
double
PerfModels<Scalar, LocalOrdinal, GlobalOrdinal, Node>::stream_vector_copy_lookup(int SIZE_IN_BYTES) {
  return PerfDetails::table_lookup(stream_sizes_, stream_copy_times_, SIZE_IN_BYTES / sizeof(Scalar));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
double
PerfModels<Scalar, LocalOrdinal, GlobalOrdinal, Node>::stream_vector_add_lookup(int SIZE_IN_BYTES) {
  return PerfDetails::table_lookup(stream_sizes_, stream_add_times_, SIZE_IN_BYTES / sizeof(Scalar));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
double
PerfModels<Scalar, LocalOrdinal, GlobalOrdinal, Node>::stream_vector_lookup(int SIZE_IN_BYTES) {
  return std::min(stream_vector_copy_lookup(SIZE_IN_BYTES), stream_vector_add_lookup(SIZE_IN_BYTES));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
double
PerfModels<Scalar, LocalOrdinal, GlobalOrdinal, Node>::latency_corrected_stream_vector_copy_lookup(int SIZE_IN_BYTES) {
  return PerfDetails::table_lookup(stream_sizes_, latency_corrected_stream_copy_times_, SIZE_IN_BYTES / sizeof(Scalar));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
double
PerfModels<Scalar, LocalOrdinal, GlobalOrdinal, Node>::latency_corrected_stream_vector_add_lookup(int SIZE_IN_BYTES) {
  return PerfDetails::table_lookup(stream_sizes_, latency_corrected_stream_add_times_, SIZE_IN_BYTES / sizeof(Scalar));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
double
PerfModels<Scalar, LocalOrdinal, GlobalOrdinal, Node>::latency_corrected_stream_vector_lookup(int SIZE_IN_BYTES) {
  return std::min(latency_corrected_stream_vector_copy_lookup(SIZE_IN_BYTES), latency_corrected_stream_vector_add_lookup(SIZE_IN_BYTES));
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void PerfModels<Scalar, LocalOrdinal, GlobalOrdinal, Node>::print_stream_vector_table(std::ostream &out, const std::string &prefix) {
  print_stream_vector_table_impl(out, false, prefix);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void PerfModels<Scalar, LocalOrdinal, GlobalOrdinal, Node>::print_latency_corrected_stream_vector_table(std::ostream &out, const std::string &prefix) {
  print_stream_vector_table_impl(out, true, prefix);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void PerfModels<Scalar, LocalOrdinal, GlobalOrdinal, Node>::print_stream_vector_table_impl(std::ostream &out, bool use_latency_correction, const std::string &prefix) {
  using namespace std;
  std::ios old_format(NULL);
  old_format.copyfmt(out);

  out << prefix
      << setw(20) << "Length in Scalars" << setw(1) << " "
      << setw(20) << "COPY (us)" << setw(1) << " "
      << setw(20) << "ADD (us)" << setw(1) << " "
      << setw(20) << "COPY (GB/s)" << setw(1) << " "
      << setw(20) << "ADD (GB/s)" << std::endl;

  out << prefix
      << setw(20) << "-----------------" << setw(1) << " "
      << setw(20) << "---------" << setw(1) << " "
      << setw(20) << "--------" << setw(1) << " "
      << setw(20) << "-----------" << setw(1) << " "
      << setw(20) << "----------" << std::endl;

  for (int i = 0; i < (int)stream_sizes_.size(); i++) {
    int size      = stream_sizes_[i];
    double c_time = use_latency_correction ? latency_corrected_stream_copy_times_[i] : stream_copy_times_[i];
    double a_time = use_latency_correction ? latency_corrected_stream_add_times_[i] : stream_add_times_[i];
    // We've already corrected for the transactions per element difference
    double c_bw = PerfDetails::convert_time_to_bandwidth_gbs(c_time, 1, size * sizeof(Scalar));
    double a_bw = PerfDetails::convert_time_to_bandwidth_gbs(a_time, 1, size * sizeof(Scalar));

    out << prefix
        << setw(20) << size << setw(1) << " "
        << setw(20) << fixed << setprecision(4) << (c_time * 1e6) << setw(1) << " "
        << setw(20) << fixed << setprecision(4) << (a_time * 1e6) << setw(1) << " "
        << setw(20) << fixed << setprecision(4) << c_bw << setw(1) << " "
        << setw(20) << fixed << setprecision(4) << a_bw << std::endl;
  }

  out.copyfmt(old_format);
}

/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void PerfModels<Scalar, LocalOrdinal, GlobalOrdinal, Node>::pingpong_make_table(int KERNEL_REPEATS, int LOG_MAX_SIZE, const RCP<const Teuchos::Comm<int> > &comm) {
  PerfDetails::pingpong_basic<Kokkos::HostSpace::execution_space, Kokkos::HostSpace::memory_space>(KERNEL_REPEATS, LOG_MAX_SIZE, *comm, pingpong_sizes_, pingpong_host_times_);

  PerfDetails::pingpong_basic<typename Node::execution_space, typename Node::memory_space>(KERNEL_REPEATS, LOG_MAX_SIZE, *comm, pingpong_sizes_, pingpong_device_times_);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
double
PerfModels<Scalar, LocalOrdinal, GlobalOrdinal, Node>::pingpong_host_lookup(int SIZE_IN_BYTES) {
  return PerfDetails::table_lookup(pingpong_sizes_, pingpong_host_times_, SIZE_IN_BYTES);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
double
PerfModels<Scalar, LocalOrdinal, GlobalOrdinal, Node>::pingpong_device_lookup(int SIZE_IN_BYTES) {
  return PerfDetails::table_lookup(pingpong_sizes_, pingpong_device_times_, SIZE_IN_BYTES);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void PerfModels<Scalar, LocalOrdinal, GlobalOrdinal, Node>::print_pingpong_table(std::ostream &out, const std::string &prefix) {
  if (pingpong_sizes_.size() == 0) return;

  using namespace std;
  std::ios old_format(NULL);
  old_format.copyfmt(out);

  out << prefix
      << setw(20) << "Message Size" << setw(1) << " "
      << setw(20) << "Host (us)" << setw(1) << " "
      << setw(20) << "Device (us)" << std::endl;

  out << prefix
      << setw(20) << "------------" << setw(1) << " "
      << setw(20) << "---------" << setw(1) << " "
      << setw(20) << "-----------" << std::endl;

  for (int i = 0; i < (int)pingpong_sizes_.size(); i++) {
    int size      = pingpong_sizes_[i];
    double h_time = pingpong_host_times_[i];
    double d_time = pingpong_device_times_[i];

    out << prefix
        << setw(20) << size << setw(1) << " "
        << setw(20) << fixed << setprecision(4) << (h_time * 1e6) << setw(1) << " "
        << setw(20) << fixed << setprecision(4) << (d_time * 1e6) << setw(1) << std::endl;
  }

  out.copyfmt(old_format);
}

/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void PerfModels<Scalar, LocalOrdinal, GlobalOrdinal, Node>::halopong_make_table(int KERNEL_REPEATS, int LOG_MAX_SIZE, const RCP<const Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> > &import) {
  PerfDetails::halopong_basic<Kokkos::HostSpace::execution_space, Kokkos::HostSpace::memory_space>(KERNEL_REPEATS, LOG_MAX_SIZE, import, halopong_sizes_, halopong_host_times_);

  PerfDetails::halopong_basic<typename Node::execution_space, typename Node::memory_space>(KERNEL_REPEATS, LOG_MAX_SIZE, import, halopong_sizes_, halopong_device_times_);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
double
PerfModels<Scalar, LocalOrdinal, GlobalOrdinal, Node>::halopong_host_lookup(int SIZE_IN_BYTES) {
  return PerfDetails::table_lookup(halopong_sizes_, halopong_host_times_, SIZE_IN_BYTES);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
double
PerfModels<Scalar, LocalOrdinal, GlobalOrdinal, Node>::halopong_device_lookup(int SIZE_IN_BYTES) {
  return PerfDetails::table_lookup(halopong_sizes_, halopong_device_times_, SIZE_IN_BYTES);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void PerfModels<Scalar, LocalOrdinal, GlobalOrdinal, Node>::print_halopong_table(std::ostream &out, const std::string &prefix) {
  if (halopong_sizes_.size() == 0) return;

  using namespace std;
  std::ios old_format(NULL);
  old_format.copyfmt(out);

  out << prefix
      << setw(20) << "Message Size" << setw(1) << " "
      << setw(20) << "Host (us)" << setw(1) << " "
      << setw(20) << "Device (us)" << std::endl;

  out << prefix
      << setw(20) << "------------" << setw(1) << " "
      << setw(20) << "---------" << setw(1) << " "
      << setw(20) << "-----------" << std::endl;

  for (int i = 0; i < (int)halopong_sizes_.size(); i++) {
    int size      = halopong_sizes_[i];
    double h_time = halopong_host_times_[i];
    double d_time = halopong_device_times_[i];

    out << prefix
        << setw(20) << size << setw(1) << " "
        << setw(20) << fixed << setprecision(4) << (h_time * 1e6) << setw(1) << " "
        << setw(20) << fixed << setprecision(4) << (d_time * 1e6) << setw(1) << std::endl;
  }

  out.copyfmt(old_format);
}

/****************************************************************************************/
/****************************************************************************************/
/****************************************************************************************/

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void PerfModels<Scalar, LocalOrdinal, GlobalOrdinal, Node>::launch_latency_make_table(int KERNEL_REPEATS) {
  using exec_space   = typename Node::execution_space;
  using range_policy = Kokkos::RangePolicy<exec_space>;
  using clock        = std::chrono::high_resolution_clock;

  double total_test_time = 0;
  clock::time_point start, stop;
  for (int i = 0; i < KERNEL_REPEATS; i++) {
    start = clock::now();
    Kokkos::parallel_for(
        "empty kernel", range_policy(0, 1), KOKKOS_LAMBDA(const size_t j) {
          ;
        });
    exec_space().fence();
    stop                = clock::now();
    double my_test_time = std::chrono::duration<double>(stop - start).count();
    total_test_time += my_test_time;
  }

  launch_and_wait_latency_ = total_test_time / KERNEL_REPEATS;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
double
PerfModels<Scalar, LocalOrdinal, GlobalOrdinal, Node>::launch_latency_lookup() {
  return launch_and_wait_latency_;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void PerfModels<Scalar, LocalOrdinal, GlobalOrdinal, Node>::print_launch_latency_table(std::ostream &out, const std::string &prefix) {
  using namespace std;
  std::ios old_format(NULL);
  old_format.copyfmt(out);

  out << prefix
      << setw(20) << "Launch+Wait Latency (us)" << setw(1) << " "
      << setw(20) << fixed << setprecision(4) << (launch_and_wait_latency_ * 1e6) << std::endl;

  out.copyfmt(old_format);
}

}  // namespace MueLu
