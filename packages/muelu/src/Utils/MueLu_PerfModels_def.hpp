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
#include <cmath>
#include <numeric>
#include <utility>
#include <chrono>
#include <iomanip>
#include <Teuchos_ScalarTraits.hpp>
#include "MueLu_PerfModels_decl.hpp"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

namespace MueLu {

  namespace PerfDetails {
    template<class Scalar,class Node>
    double stream_vector_add(int KERNEL_REPEATS, int VECTOR_SIZE) {      
      using exec_space   = typename Node::execution_space;
      using memory_space = typename Node::memory_space;
      using range_policy = Kokkos::RangePolicy<exec_space>;

      Kokkos::View<Scalar*,memory_space> a("a", VECTOR_SIZE);
      Kokkos::View<Scalar*,memory_space> b("b", VECTOR_SIZE);
      Kokkos::View<Scalar*,memory_space> c("c", VECTOR_SIZE);
      double total_test_time = 0.0;

      Scalar ONE = Teuchos::ScalarTraits<Scalar>::one();

      Kokkos::parallel_for("stream/fill",range_policy(0,VECTOR_SIZE), KOKKOS_LAMBDA (const size_t i) {
          a(i) = ONE * (double)i;
          b(i) = a(i);
        });
      exec_space().fence();

      using clock = std::chrono::high_resolution_clock;

      clock::time_point start, stop;

      for(int i = 0; i < KERNEL_REPEATS; i++) {       
        start = clock::now();
        Kokkos::parallel_for("stream/add",range_policy(0,VECTOR_SIZE), KOKKOS_LAMBDA (const size_t j) { //Vector Addition
            c(j) = a(j) + b(j);
        });

        exec_space().fence();
        stop = clock::now();
        double my_test_time = std::chrono::duration<double>(stop - start).count();
        total_test_time += my_test_time;
      }

      return total_test_time / KERNEL_REPEATS;
    }

    template<class Scalar,class Node>
    double stream_vector_copy(int KERNEL_REPEATS, int VECTOR_SIZE) {      
      using exec_space   = typename Node::execution_space;
      using memory_space = typename Node::memory_space;
      using range_policy = Kokkos::RangePolicy<exec_space>;

      Kokkos::View<Scalar*,memory_space> a("a", VECTOR_SIZE);
      Kokkos::View<Scalar*,memory_space> b("b", VECTOR_SIZE);
      double total_test_time = 0.0;

      Scalar ONE = Teuchos::ScalarTraits<Scalar>::one();

      Kokkos::parallel_for("stream/fill",range_policy(0,VECTOR_SIZE), KOKKOS_LAMBDA (const size_t i) {
          a(i) = ONE;
        });
      exec_space().fence();

      using clock = std::chrono::high_resolution_clock;
      clock::time_point start, stop;

      for(int i = 0; i < KERNEL_REPEATS; i++) {       
        start = clock::now();
        Kokkos::parallel_for("stream/copy",range_policy(0,VECTOR_SIZE), KOKKOS_LAMBDA (const size_t j) { //Vector Addition
            b(j) = a(j);
        });

        exec_space().fence();
        stop = clock::now();
        double my_test_time = std::chrono::duration<double>(stop - start).count();
        total_test_time += my_test_time;
      }

      return total_test_time / KERNEL_REPEATS;
    }


 
    double table_lookup(const std::vector<int> & x, const std::vector<double> & y, int value) {
      // If there's no table, nan
      if(x.size() == 0) return Teuchos::ScalarTraits<double>::nan();
      
      // NOTE:  This should probably be a binary search, but this isn't performance sensitive, so we'll go simple
      int N = (int) x.size();
      int hi = 0;
      for(  ; hi < N; hi++) {
        if (x[hi] > value) 
          break;
      }
      
      if(hi == 0) {
        // Lower end (return the min time)
        //printf("Lower end: %d < %d\n",value,x[0]);
        return y[0];
      }
      else if (hi == N) {
        // Higher end (extrapolate from the last two points)
        //printf("Upper end: %d > %d\n",value,x[N-1]);
        hi = N-1;
        int run     = x[hi] - x[hi-1];
        double rise = y[hi] - y[hi-1];
        double slope = rise / run;     
        int diff = value - x[hi-1];
        
        return y[hi-1] + slope * diff;
      }
      else {
        // Interpolate
        //printf("Middle: %d < %d < %d\n",x[hi-1],value,x[hi]);
        int run     = x[hi] - x[hi-1];
        double rise = y[hi] - y[hi-1];
        double slope = rise / run;     
        int diff = value - x[hi-1];
        
        return y[hi-1] + slope * diff;
      }
    }

    // Report bandwidth in GB / sec
    const double GB = 1024.0 * 1024.0 * 1024.0;
    double convert_time_to_bandwidth_gbs(double time, int num_calls, double memory_per_call_bytes) {
      double time_per_call = time / num_calls;      
      return memory_per_call_bytes / GB / time_per_call;
    }


    template <class exec_space, class memory_space>
    void pingpong_basic(int KERNEL_REPEATS, int MAX_SIZE,const Teuchos::Comm<int> &comm, std::vector<int> & sizes, std::vector<double> & times) {
      int rank = comm.getRank();
      int nproc = comm.getSize();
      
      if(nproc < 2) return;
      
#ifdef HAVE_MPI
      using range_policy = Kokkos::RangePolicy<exec_space>;
      const int buff_size = (int) pow(2,MAX_SIZE);
      
      sizes.resize(MAX_SIZE+1);
      times.resize(MAX_SIZE+1);
            
      // Allocate memory for the buffers (and fill send)
      Kokkos::View<char*,memory_space> r_buf("recv",buff_size), s_buf("send",buff_size);
      Kokkos::deep_copy(s_buf,1);
      
      //Send and recieve.   
      // NOTE:  Do consectutive pair buddies here for simplicity.  We should be smart later
      int odd = rank % 2;
      int buddy = odd ? rank - 1 : rank + 1;
      
      for(int i = 0; i < MAX_SIZE + 1 ;i ++) {
        int msg_size = (int) pow(2,i);
        comm.barrier();

        double t0 = MPI_Wtime();      
        for(int j = 0; j < KERNEL_REPEATS; j++) {
          if (buddy < nproc) {
            if (odd) {
              comm.send(msg_size, (char*)s_buf.data(), buddy);
              comm.receive(buddy, msg_size, (char*)r_buf.data());
            }
            else {
              comm.receive(buddy, msg_size,(char*)r_buf.data());
              comm.send(msg_size, (char*)s_buf.data(), buddy);
            }
          }
        }

        double time_per_call = (MPI_Wtime() - t0) / (2.0 * KERNEL_REPEATS);
        sizes[i] = msg_size;
        times[i] = time_per_call;
      }
#endif
    }



  }// end namespace PerfDetails


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  PerfModels<Scalar, LocalOrdinal, GlobalOrdinal, Node>::PerfModels():launch_and_wait_latency_(-1.0){}


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void 
  PerfModels<Scalar, LocalOrdinal, GlobalOrdinal, Node>::stream_vector_make_table(int KERNEL_REPEATS,  int LOG_MAX_SIZE) {
    // We need launch/waits latency estimates for corrected stream
    launch_latency_make_table(KERNEL_REPEATS);
    double latency = launch_latency_lookup();
    
    if(LOG_MAX_SIZE < 2)
      LOG_MAX_SIZE=20;

    stream_sizes_.resize(LOG_MAX_SIZE+1);    
    stream_copy_times_.resize(LOG_MAX_SIZE+1);    
    stream_add_times_.resize(LOG_MAX_SIZE+1);    
    latency_corrected_stream_copy_times_.resize(LOG_MAX_SIZE+1);    
    latency_corrected_stream_add_times_.resize(LOG_MAX_SIZE+1);    
        
    for(int i=0; i<LOG_MAX_SIZE+1; i++) {
      int size = (int) pow(2,i);
      double c_time = PerfDetails::stream_vector_copy<Scalar,Node>(KERNEL_REPEATS,size);
      double a_time = PerfDetails::stream_vector_add<Scalar,Node>(KERNEL_REPEATS,size);

      stream_sizes_[i] = size;

      // Correct for the difference in memory transactions per element
      stream_copy_times_[i] = c_time / 2.0;
      stream_add_times_[i] = a_time / 3.0;

      // Correct for launch latency too.  We'll note that sometimes the latency estimate
      // is higher than the actual copy/add time estimate.  If so, we don't correct
      latency_corrected_stream_copy_times_[i] = (c_time - latency <= 0.0) ? c_time / 2.0 : ( (c_time-latency)/2.0 );
      latency_corrected_stream_add_times_[i]  = (a_time - latency <= 0.0) ? a_time / 3.0 : ( (a_time-latency)/3.0 );


    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  double
  PerfModels<Scalar, LocalOrdinal, GlobalOrdinal, Node>::stream_vector_copy_lookup(int SIZE_IN_BYTES) {
    return PerfDetails::table_lookup(stream_sizes_,stream_copy_times_,SIZE_IN_BYTES/sizeof(Scalar));
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  double
  PerfModels<Scalar, LocalOrdinal, GlobalOrdinal, Node>::stream_vector_add_lookup(int SIZE_IN_BYTES) {
    return PerfDetails::table_lookup(stream_sizes_,stream_add_times_,SIZE_IN_BYTES/sizeof(Scalar));
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  double
  PerfModels<Scalar, LocalOrdinal, GlobalOrdinal, Node>::stream_vector_lookup(int SIZE_IN_BYTES) {
    return std::min(stream_vector_copy_lookup(SIZE_IN_BYTES),stream_vector_add_lookup(SIZE_IN_BYTES));
  }
  

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  double
  PerfModels<Scalar, LocalOrdinal, GlobalOrdinal, Node>::latency_corrected_stream_vector_copy_lookup(int SIZE_IN_BYTES) {
    return PerfDetails::table_lookup(stream_sizes_,latency_corrected_stream_copy_times_,SIZE_IN_BYTES/sizeof(Scalar));
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  double
  PerfModels<Scalar, LocalOrdinal, GlobalOrdinal, Node>::latency_corrected_stream_vector_add_lookup(int SIZE_IN_BYTES) {
    return PerfDetails::table_lookup(stream_sizes_,latency_corrected_stream_add_times_,SIZE_IN_BYTES/sizeof(Scalar));
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  double
  PerfModels<Scalar, LocalOrdinal, GlobalOrdinal, Node>::latency_corrected_stream_vector_lookup(int SIZE_IN_BYTES) {
    return std::min(latency_corrected_stream_vector_copy_lookup(SIZE_IN_BYTES),latency_corrected_stream_vector_add_lookup(SIZE_IN_BYTES));
  }
 

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  PerfModels<Scalar, LocalOrdinal, GlobalOrdinal, Node>::print_stream_vector_table(std::ostream & out) {
    print_stream_vector_table_impl(out,false);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  PerfModels<Scalar, LocalOrdinal, GlobalOrdinal, Node>::print_latency_corrected_stream_vector_table(std::ostream & out) {
    print_stream_vector_table_impl(out,true);
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  PerfModels<Scalar, LocalOrdinal, GlobalOrdinal, Node>::print_stream_vector_table_impl(std::ostream & out,bool use_latency_correction) {
    using namespace std;
    std::ios old_format(NULL);
    old_format.copyfmt(out);

    out << setw(20) << "Length in Scalars" << setw(1) << " "
        << setw(20) << "COPY (us)" << setw(1) << " "
        << setw(20) << "ADD (us)" << setw(1) << " "
        << setw(20) << "COPY (GB/s)" << setw(1) << " "
        << setw(20) << "ADD (GB/s)" << std::endl;

    out << setw(20) << "-----------------" << setw(1) << " "
        << setw(20) << "---------" << setw(1) << " "
        << setw(20) << "--------" << setw(1) << " "
        << setw(20) << "-----------" << setw(1) << " "
        << setw(20) << "----------" << std::endl;
    

    for(int i=0; i<(int)stream_sizes_.size(); i++) {
      int size = stream_sizes_[i];
      double c_time = use_latency_correction ? latency_corrected_stream_copy_times_[i] : stream_copy_times_[i];
      double a_time = use_latency_correction ? latency_corrected_stream_add_times_[i]  : stream_add_times_[i];
      // We've already corrected for the transactions per element difference
      double c_bw = PerfDetails::convert_time_to_bandwidth_gbs(c_time,1,size*sizeof(Scalar));
      double a_bw = PerfDetails::convert_time_to_bandwidth_gbs(a_time,1,size*sizeof(Scalar));


      out << setw(20) << size << setw(1) << " "
          << setw(20) << fixed << setprecision(4) << (c_time*1e6) << setw(1) << " "
          << setw(20) << fixed << setprecision(4) << (a_time*1e6) << setw(1) << " "
          << setw(20) << fixed << setprecision(4) << c_bw << setw(1) << " "
          << setw(20) << fixed << setprecision(4) << a_bw << std::endl;        
    }

    out.copyfmt(old_format);
  }


 

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
 PerfModels<Scalar, LocalOrdinal, GlobalOrdinal, Node>::pingpong_make_table(int KERNEL_REPEATS, int LOG_MAX_SIZE, const RCP<const Teuchos::Comm<int> > &comm) {

    PerfDetails::pingpong_basic<Kokkos::HostSpace::execution_space,Kokkos::HostSpace::memory_space>(KERNEL_REPEATS,LOG_MAX_SIZE,*comm,pingpong_sizes_,pingpong_host_times_);

    PerfDetails::pingpong_basic<typename Node::execution_space,typename Node::memory_space>(KERNEL_REPEATS,LOG_MAX_SIZE,*comm,pingpong_sizes_,pingpong_device_times_);

 }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  double
  PerfModels<Scalar, LocalOrdinal, GlobalOrdinal, Node>::pingpong_host_lookup(int SIZE_IN_BYTES) {
    return PerfDetails::table_lookup(pingpong_sizes_,pingpong_host_times_,SIZE_IN_BYTES);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  double
  PerfModels<Scalar, LocalOrdinal, GlobalOrdinal, Node>::pingpong_device_lookup(int SIZE_IN_BYTES) {
    return PerfDetails::table_lookup(pingpong_sizes_,pingpong_device_times_,SIZE_IN_BYTES);
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  PerfModels<Scalar, LocalOrdinal, GlobalOrdinal, Node>::print_pingpong_table(std::ostream & out) {
    if(pingpong_sizes_.size() == 0) return;

    using namespace std;
    std::ios old_format(NULL);
    old_format.copyfmt(out);

    out << setw(20) << "Message Size" << setw(1) << " "
        << setw(20) << "Host (us)" << setw(1) << " "
        << setw(20) << "Device (us)" << std::endl;

    out << setw(20) << "------------" << setw(1) << " "
        << setw(20) << "---------" << setw(1) << " "
        << setw(20) << "-----------" << std::endl;
    

    for(int i=0; i<(int)pingpong_sizes_.size(); i++) {
      int size = pingpong_sizes_[i];
      double h_time = pingpong_host_times_[i];
      double d_time = pingpong_device_times_[i];


      out << setw(20) << size << setw(1) << " "
          << setw(20) << fixed << setprecision(4) << (h_time*1e6) << setw(1) << " "
          << setw(20) << fixed << setprecision(4) << (d_time*1e6) << setw(1) << std::endl;
    }

    out.copyfmt(old_format);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  PerfModels<Scalar, LocalOrdinal, GlobalOrdinal, Node>::launch_latency_make_table(int KERNEL_REPEATS) {
    using exec_space   = typename Node::execution_space;
    using range_policy = Kokkos::RangePolicy<exec_space>;
    using clock = std::chrono::high_resolution_clock;

    double total_test_time = 0;
    clock::time_point start, stop;
    for(int i = 0; i < KERNEL_REPEATS; i++) {       
      start = clock::now();
      Kokkos::parallel_for("empty kernel",range_policy(0,1), KOKKOS_LAMBDA (const size_t j) {
          ;
        });
      exec_space().fence();
      stop = clock::now();
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
  void
  PerfModels<Scalar, LocalOrdinal, GlobalOrdinal, Node>::print_launch_latency_table(std::ostream & out) {
    using namespace std;
    std::ios old_format(NULL);
    old_format.copyfmt(out);

    out << setw(20) << "Launch+Wait Latency (us)" << setw(1) << " " 
        << setw(20) << fixed << setprecision(4) << (launch_and_wait_latency_*1e6) << std::endl;

    out.copyfmt(old_format);
  }


} //namespace MueLu
