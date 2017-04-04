// Copyright(C) 1999-2010
// Sandia Corporation. Under the terms of Contract
// DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
// certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef IOSS_Ioss_ParallelUtils_h
#define IOSS_Ioss_ParallelUtils_h

#include <Ioss_CodeTypes.h> // for Int64Vector, IntVector
#include <Ioss_Utils.h>
#include <assert.h>
#include <stddef.h> // for size_t
#include <string>   // for string
#include <vector>   // for vector

#ifdef HAVE_MPI
#include <mpi.h>
#endif

namespace Ioss {

  class ParallelUtils
  {
  public:
    explicit ParallelUtils(MPI_Comm the_communicator);
    ~ParallelUtils() = default;

    // Assignment operator
    // Copy constructor

    enum MinMax { DO_MAX, DO_MIN, DO_SUM };

    /*!
     * Returns 'true' if 'name' is defined in the environment.
     * The value of the environment variable is returned in 'value'.
     * getenv system call is only done on processor 0.
     * If '!sync_parallel', then don't push to other processors.
     */
    bool get_environment(const std::string &name, std::string &value, bool sync_parallel) const;

    /*!
     * Returns 'true' if 'name' is defined in the environment.
     * The value of the environment variable is converted to an
     * integer via the atoi library call and returned in 'value'.
     * No checking is done to ensure that the environment variable
     * points to a valid integer.
     * getenv system call is only done on processor 0.
     * If '!sync_parallel', then don't push to other processors.
     */
    bool get_environment(const std::string &name, int &value, bool sync_parallel) const;

    /*!
     * Returns 'true' if 'name' is defined in the environment no
     * matter what the value. Returns false otherwise.
     * getenv system call is only done on processor 0.
     * If '!sync_parallel', then don't push to other processors.
     */
    bool get_environment(const std::string &name, bool sync_parallel) const;

    std::string decode_filename(const std::string &filename, bool is_parallel) const;

    MPI_Comm communicator() const { return communicator_; }
    int      parallel_size() const;
    int      parallel_rank() const;

    /*!
     * Global OR of attribute strings, the processors which have no
     * knowledge of the value should initialize to '0' and the
     * processors with knowledge set the appropriate values.
     */
    void attribute_reduction(const int length, char buffer[]) const;

    /*! Return min, max, average memory used by any process */
    void memory_stats(int64_t &min, int64_t &max, int64_t &avg) const;

    /*! Return high-water-mark min, max, average memory used by any process */
    /* May be inaccurate unless system maintains this information */
    void hwm_memory_stats(int64_t &min, int64_t &max, int64_t &avg) const;

    /*! Vector 'local_counts' contains the number of objects
     * local to this processor.  On exit, global_counts
     * contains the total number of objects on all processors.
     * Assumes that ordering is the same on all processors
     */
    void global_count(const IntVector &local_counts, IntVector &global_counts) const;
    void global_count(const Int64Vector &local_counts, Int64Vector &global_counts) const;

    template <typename T> T global_minmax(T local_minmax, MinMax which) const;
    template <typename T>
    void global_array_minmax(std::vector<T> &local_minmax, MinMax which) const;
    template <typename T>
    void global_array_minmax(T *local_minmax, size_t count, MinMax which) const;

    template <typename T> void gather(T my_value, std::vector<T> &result) const;
    template <typename T> void all_gather(T my_value, std::vector<T> &result) const;
    template <typename T> void gather(std::vector<T> &my_values, std::vector<T> &result) const;

    void progress(const std::string &output) const;
    
  private:
    MPI_Comm communicator_;
  };

  inline int power_2(int count)
  {
    // Return the power of two which is equal to or greater than 'count'
    // count = 15 -> returns 16
    // count = 16 -> returns 16
    // count = 17 -> returns 32

    // Use brute force...
    int pow2 = 1;
    while (pow2 < count) {
      pow2 *= 2;
    }
    return pow2;
  }

#ifdef HAVE_MPI
  inline MPI_Datatype mpi_type(double /*dummy*/) { return MPI_DOUBLE; }
  inline MPI_Datatype mpi_type(int /*dummy*/) { return MPI_INT; }
  inline MPI_Datatype mpi_type(int64_t /*dummy*/) { return MPI_LONG_LONG_INT; }
  inline MPI_Datatype mpi_type(unsigned int /*dummy*/) { return MPI_UNSIGNED; }

  template <typename T>
  int MY_Alltoallv64(std::vector<T> &sendbuf, const std::vector<int64_t> &sendcounts,
                     const std::vector<int64_t> &senddisp, std::vector<T> &recvbuf,
                     const std::vector<int64_t> &recvcounts, const std::vector<int64_t> &recvdisp,
                     MPI_Comm comm)
  {
    int processor_count = 0;
    int my_processor    = 0;
    MPI_Comm_size(comm, &processor_count);
    MPI_Comm_rank(comm, &my_processor);

    // Verify that all 'counts' can fit in an integer. Symmetric
    // communication, so recvcounts are sendcounts on another processor.
    for (int i = 0; i < processor_count; i++) {
      int snd_cnt = (int)sendcounts[i];
      if ((int64_t)snd_cnt != sendcounts[i]) {
        std::ostringstream errmsg;
        errmsg << "ERROR: The number of items that must be communicated via MPI calls from\n"
               << "       processor " << my_processor << " to processor " << i << " is "
               << sendcounts[i] << "\n       which exceeds the storage capacity of the integers "
                                   "used by MPI functions.\n";
        std::cerr << errmsg.str();
        exit(EXIT_FAILURE);
      }
    }

    size_t pow_2 = power_2(processor_count);

    for (size_t i = 1; i < pow_2; i++) {
      MPI_Status status;

      int    tag           = 24713;
      size_t exchange_proc = i ^ my_processor;
      if (exchange_proc < (size_t)processor_count) {
        int snd_cnt =
            (int)sendcounts[exchange_proc]; // Converts from int64_t to int as needed by mpi
        int rcv_cnt = (int)recvcounts[exchange_proc];
        if ((size_t)my_processor < exchange_proc) {
          MPI_Send(&sendbuf[senddisp[exchange_proc]], snd_cnt, mpi_type(T(0)), exchange_proc, tag,
                   comm);
          MPI_Recv(&recvbuf[recvdisp[exchange_proc]], rcv_cnt, mpi_type(T(0)), exchange_proc, tag,
                   comm, &status);
        }
        else {
          MPI_Recv(&recvbuf[recvdisp[exchange_proc]], rcv_cnt, mpi_type(T(0)), exchange_proc, tag,
                   comm, &status);
          MPI_Send(&sendbuf[senddisp[exchange_proc]], snd_cnt, mpi_type(T(0)), exchange_proc, tag,
                   comm);
        }
      }
    }

    // Take care of this processor's data movement...
    std::copy(&sendbuf[senddisp[my_processor]],
              &sendbuf[senddisp[my_processor] + sendcounts[my_processor]],
              &recvbuf[recvdisp[my_processor]]);
    return 0;
  }

  template <typename T>
  int MY_Alltoallv(std::vector<T> &sendbuf, const std::vector<int64_t> &sendcnts,
                   const std::vector<int64_t> &senddisp, std::vector<T> &recvbuf,
                   const std::vector<int64_t> &recvcnts, const std::vector<int64_t> &recvdisp,
                   MPI_Comm comm)
  {
    // Wrapper to handle case where send/recv counts and displacements are 64-bit integers.
    // Two cases:
    // 1) They are of type 64-bit integers, but only storing data in the 32-bit integer range.
    //    -- if (sendcnts[#proc-1] + senddisp[#proc-1] < 2^31, then we are ok
    // 2) They are of type 64-bit integers, and storing data in the 64-bit integer range.
    //    -- call special alltoallv which does point-to-point sends
#if 1
    int processor_count = 0;
    MPI_Comm_size(comm, &processor_count);
    size_t max_comm = sendcnts[processor_count - 1] + senddisp[processor_count - 1];
    size_t one      = 1;
    if (max_comm < one << 31) {
      // count and displacement data in range, need to copy to integer vector.
      std::vector<int> send_cnt(sendcnts.begin(), sendcnts.end());
      std::vector<int> send_dis(senddisp.begin(), senddisp.end());
      std::vector<int> recv_cnt(recvcnts.begin(), recvcnts.end());
      std::vector<int> recv_dis(recvdisp.begin(), recvdisp.end());
      return MPI_Alltoallv(TOPTR(sendbuf), TOPTR(send_cnt), TOPTR(send_dis), mpi_type(T(0)),
                           TOPTR(recvbuf), TOPTR(recv_cnt), TOPTR(recv_dis), mpi_type(T(0)), comm);
    }
    else {
#endif
      // Same as if each processor sent a message to every other process with:
      //     MPI_Send(sendbuf+senddisp[i]*sizeof(sendtype),sendcnts[i], sendtype, i, tag, comm);
      // And received a message from each processor with a call to:
      //     MPI_Recv(recvbuf+recvdisp[i]*sizeof(recvtype),recvcnts[i], recvtype, i, tag, comm);
      return MY_Alltoallv64(sendbuf, sendcnts, senddisp, recvbuf, recvcnts, recvdisp, comm);
#if 1
    }
#endif
  }

  template <typename T>
  int MY_Alltoallv(std::vector<T> &sendbuf, const std::vector<int> &sendcnts,
                   const std::vector<int> &senddisp, std::vector<T> &recvbuf,
                   const std::vector<int> &recvcnts, const std::vector<int> &recvdisp,
                   MPI_Comm comm)
  {
    return MPI_Alltoallv(TOPTR(sendbuf), (int *)TOPTR(sendcnts), (int *)TOPTR(senddisp),
                         mpi_type(T(0)), TOPTR(recvbuf), (int *)TOPTR(recvcnts),
                         (int *)TOPTR(recvdisp), mpi_type(T(0)), comm);
  }
#endif
}
#endif
