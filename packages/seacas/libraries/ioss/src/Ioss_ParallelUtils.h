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

#include <Ioss_CodeTypes.h>
#include <string>
#include <vector>

namespace Ioss {

  typedef std::vector<int> IntVector;

  class ParallelUtils {
  public:

    explicit ParallelUtils(MPI_Comm communicator);
    ~ParallelUtils() {};

    // Assignment operator
    // Copy constructor
    
    enum MinMax {DO_MAX, DO_MIN};

    /*! 
     * Returns 'true' if 'name' is defined in the environment.
     * The value of the environment variable is returned in 'value'. 
     * getenv system call is only done on processor 0.
     * If '!sync_parallel', then don't push to other processors.
     */
    bool get_environment(const std::string &name, std::string &value,
			 bool sync_parallel) const;

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

    MPI_Comm communicator() const {return communicator_;}
    int parallel_size() const;
    int parallel_rank() const;

    /*!
     * Global OR of attribute strings, the processors which have no
     * knowledge of the value should initialize to '0' and the
     * processors with knowledge set the appropriate values.
     */
    void attribute_reduction( const int length , char buffer[]) const;

    /*! Vector 'local_counts' contains the number of objects
     * local to this processor.  On exit, global_counts
     * contains the total number of objects on all processors.
     * Assumes that ordering is the same on all processors
     */
    void global_count(const IntVector &local_counts, IntVector &global_counts) const;

    int  global_minmax(int  local_minmax, MinMax which) const;
    double global_minmax(double local_minmax, MinMax which) const;
    unsigned int  global_minmax(unsigned int local_minmax, MinMax which) const;
    void global_array_minmax(int *local_minmax, size_t count, MinMax which) const;
    void global_array_minmax(unsigned int *local_minmax, size_t count, MinMax which) const;
    void global_array_minmax(double *local_minmax, size_t count, MinMax which) const;
    void gather(int my_value, std::vector<int> &result) const;
    void gather(std::vector<int> &my_values, std::vector<int> &result) const;

  private:
    MPI_Comm communicator_;
  };
}
#endif
