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

#include <assert.h>
#include <Ioss_ParallelUtils.h>
#include <Ioss_Utils.h>
#include <stddef.h>
#include <algorithm>
#include <cstdlib>
#include <string>
#include <cstring>
#include <vector>

#include <Ioss_CodeTypes.h>
#ifdef HAVE_MPI
#include <assert.h>
#include <cstring>
#include <Ioss_SerializeIO.h>
#include <mpi.h>
#include <Ioss_SerializeIO.h>

namespace {
  MPI_Datatype mpi_type(double /*dummy*/)  {return MPI_DOUBLE;}
  MPI_Datatype mpi_type(int /*dummy*/)     {return MPI_INT;}
  MPI_Datatype mpi_type(unsigned int /*dummy*/)     {return MPI_UNSIGNED;}
  MPI_Datatype mpi_type(int64_t /*dummy*/) {return MPI_LONG_LONG;}
}
#endif

Ioss::ParallelUtils::ParallelUtils(MPI_Comm the_communicator)
  : communicator_(the_communicator)
{}
  
bool Ioss::ParallelUtils::get_environment(const std::string &name, std::string &value, bool sync_parallel) const
{
#ifdef HAVE_MPI
  char *result_string = NULL;
  char *broadcast_string = NULL;
  int string_length = 0;

  int rank = parallel_rank();
  if (rank == 0) {
    result_string = std::getenv(name.c_str());
    string_length = result_string ? (int)std::strlen(result_string) : 0;
  }

  if (sync_parallel && parallel_size() > 1) {
    MPI_Bcast(&string_length, 1, MPI_INT, 0, communicator_);

    if (string_length > 0) {
      broadcast_string = new char[string_length+1];
      if (rank == 0) {
	std::strncpy(broadcast_string, result_string, (size_t)string_length+1);
      }
      MPI_Bcast(broadcast_string, string_length+1, MPI_CHAR, 0, communicator_);
      result_string = broadcast_string;
      value = std::string(result_string);
      delete [] broadcast_string;
    } else {
      value = std::string("");
    }
  } else {
    if (rank == 0) {
      if (string_length > 0)
	value = std::string(result_string);
      else
	value = std::string("");
    }
  }
  return string_length > 0;
#else
  char *result_string = std::getenv(name.c_str());
  if (result_string != NULL) {
    value = std::string(result_string);
  } else {
    value = std::string("");
  }
  return (result_string != NULL);
#endif
}

bool Ioss::ParallelUtils::get_environment(const std::string &name, int &value, bool sync_parallel) const
{
  std::string str_value;
  bool success = get_environment(name, str_value, sync_parallel);
  if (success) {
    value = std::atoi(str_value.c_str());
  }
  return success;
}

bool Ioss::ParallelUtils::get_environment(const std::string &name, bool sync_parallel) const
{
  // Return true if 'name' defined, no matter what the value.
  // Return false if 'name' not defined.
#ifdef HAVE_MPI
  char *result_string = NULL;
  int string_length = 0;

  int rank = Ioss::ParallelUtils::parallel_rank();
  if (rank == 0) {
    result_string = std::getenv(name.c_str());
    string_length = result_string ? (int)std::strlen(result_string) : 0;
  }

  if (sync_parallel && parallel_size() > 1)
    MPI_Bcast(&string_length, 1, MPI_INT, 0, communicator_);

  return string_length > 0;
#else
  char *result_string = std::getenv(name.c_str());
  return (result_string != NULL);
#endif
}

std::string Ioss::ParallelUtils::decode_filename(const std::string &filename, bool is_parallel) const
{
  std::string decoded_filename(filename);

  if (is_parallel) {
    // Running in parallel, assume nemesis and decode what the filename
    // should be for this processor.
    int processor = parallel_rank();
    int num_processors = parallel_size();

    decoded_filename = Ioss::Utils::decode_filename(filename, processor, num_processors);
  }
  return decoded_filename;
}

int Ioss::ParallelUtils::parallel_size() const
{
  int my_size = 1;
#ifdef HAVE_MPI
  if (communicator_ != MPI_COMM_NULL) {
    MPI_Comm_size(communicator_, &my_size);
  }
#endif
  return my_size;
}

int Ioss::ParallelUtils::parallel_rank() const
{
  int my_rank = 0;
#ifdef HAVE_MPI
  if (communicator_ != MPI_COMM_NULL) {
    MPI_Comm_rank(communicator_, &my_rank);
  }
#endif
  return my_rank;
}

void Ioss::ParallelUtils::attribute_reduction( const int length , char buffer[]) const
{
#ifdef HAVE_MPI
  if ( 1 < parallel_size() ) {
    assert( sizeof(char) == 1 );

    char * const recv_buf = new char[ length ];

    std::memset( recv_buf , 0 , length );

    const int success =
      MPI_Allreduce( buffer , recv_buf , length , MPI_BYTE , MPI_BOR , communicator_);

    if ( MPI_SUCCESS != success ) {
      std::ostringstream errmsg;
      errmsg << "Ioss::ParallelUtils::attribute_reduction - MPI_Allreduce failed";
      IOSS_ERROR(errmsg);
    }

    std::memcpy( buffer , recv_buf , length );

    delete[] recv_buf ;
  }
#endif
}

void Ioss::ParallelUtils::global_count(const IntVector &local_counts, IntVector &global_counts) const
{
  // Vector 'local_counts' contains the number of objects
  // local to this processor.  On exit, global_counts
  // contains the total number of objects on all processors.
  // Assumes that ordering is the same on all processors
  global_counts.resize(local_counts.size());
#ifdef HAVE_MPI
  if (local_counts.size() > 0 && parallel_size() > 1) {
    if (Ioss::SerializeIO::isEnabled() && Ioss::SerializeIO::inBarrier()) {
      std::ostringstream errmsg;
      errmsg << "Attempting mpi while in barrier owned by " << Ioss::SerializeIO::getOwner();
      IOSS_ERROR(errmsg);
    }
    const int success = MPI_Allreduce((void*)&local_counts[0], &global_counts[0],
				      static_cast<int>(local_counts.size()),
				      MPI_INT, MPI_SUM, communicator_);
    if (success !=  MPI_SUCCESS) {
      std::ostringstream errmsg;
      errmsg  << "Ioss::ParallelUtils::global_count - MPI_Allreduce failed";
      IOSS_ERROR(errmsg);
    }
  } else {
    // Serial run, just copy local to global...
    std::copy(local_counts.begin(), local_counts.end(), global_counts.begin());
  }
#else
  std::copy(local_counts.begin(), local_counts.end(), global_counts.begin());
#endif
}

void Ioss::ParallelUtils::global_count(const Int64Vector &local_counts, Int64Vector &global_counts) const
{
  // Vector 'local_counts' contains the number of objects
  // local to this processor.  On exit, global_counts
  // contains the total number of objects on all processors.
  // Assumes that ordering is the same on all processors
  global_counts.resize(local_counts.size());
#ifdef HAVE_MPI
  if (local_counts.size() > 0 && parallel_size() > 1) {
    if (Ioss::SerializeIO::isEnabled() && Ioss::SerializeIO::inBarrier()) {
      std::ostringstream errmsg;
      errmsg << "Attempting mpi while in barrier owned by " << Ioss::SerializeIO::getOwner();
      IOSS_ERROR(errmsg);
    }
    const int success = MPI_Allreduce((void*)&local_counts[0], &global_counts[0],
				      static_cast<int>(local_counts.size()),
				      MPI_LONG_LONG_INT, MPI_SUM, communicator_);
    if (success !=  MPI_SUCCESS) {
      std::ostringstream errmsg;
      errmsg  << "Ioss::ParallelUtils::global_count - MPI_Allreduce failed";
      IOSS_ERROR(errmsg);
    }
  } else {
    // Serial run, just copy local to global...
    std::copy(local_counts.begin(), local_counts.end(), global_counts.begin());
  }
#else
  std::copy(local_counts.begin(), local_counts.end(), global_counts.begin());
#endif
}

template int  Ioss::ParallelUtils::global_minmax(int, Ioss::ParallelUtils::MinMax which) const;
template unsigned int Ioss::ParallelUtils::global_minmax(unsigned int,Ioss::ParallelUtils::MinMax which) const;
template int64_t Ioss::ParallelUtils::global_minmax(int64_t, Ioss::ParallelUtils::MinMax which) const;
template double  Ioss::ParallelUtils::global_minmax(double,  Ioss::ParallelUtils::MinMax which) const;

template <typename T>
T Ioss::ParallelUtils::global_minmax(T local_minmax, Ioss::ParallelUtils::MinMax which) const
{
  T minmax = local_minmax;

#ifdef HAVE_MPI
  if (parallel_size() > 1) {
    if (Ioss::SerializeIO::isEnabled() && Ioss::SerializeIO::inBarrier()) {
      std::ostringstream errmsg;
      errmsg << "Attempting mpi while in barrier owned by " << Ioss::SerializeIO::getOwner();
      IOSS_ERROR(errmsg);
    }
    static T inbuf[1], outbuf[1];
    inbuf[0] = local_minmax;

    MPI_Op oper;
    if (which == DO_MAX)
      oper = MPI_MAX;
    else if (which == DO_MIN)
      oper = MPI_MIN;
    else if (which == DO_SUM)
      oper = MPI_SUM;

    const int success = MPI_Allreduce((void*)&inbuf[0], &outbuf[0], 1,
				      mpi_type(T(1)), oper, communicator_);
    if (success !=  MPI_SUCCESS) {
      std::ostringstream errmsg;
      errmsg << "Ioss::ParallelUtils::global_minmax - MPI_Allreduce failed";
      IOSS_ERROR(errmsg);
    }
    minmax = outbuf[0];
  }
#endif
  return minmax;
}

template <typename T>
void Ioss::ParallelUtils::global_array_minmax(T *local_minmax, size_t count, Ioss::ParallelUtils::MinMax which) const
{
#ifdef HAVE_MPI
  if (parallel_size() > 1 && count > 0) {
    if (Ioss::SerializeIO::isEnabled() && Ioss::SerializeIO::inBarrier()) {
      std::ostringstream errmsg;
      errmsg << "Attempting mpi while in barrier owned by " << Ioss::SerializeIO::getOwner();
      IOSS_ERROR(errmsg);
    }

    std::vector<T> maxout(count);

    MPI_Op oper;
    if (which == DO_MAX)
      oper = MPI_MAX;
    else if (which == DO_MIN)
      oper = MPI_MIN;
    else if (which == DO_SUM)
      oper = MPI_SUM;

    const int success = MPI_Allreduce((void*)&local_minmax[0], &maxout[0],
				      static_cast<int>(count),
				      mpi_type(T(1)), oper, communicator_);
    if (success !=  MPI_SUCCESS) {
      std::ostringstream errmsg;
      errmsg << "Ioss::ParallelUtils::global_array_minmax - MPI_Allreduce failed";
      IOSS_ERROR(errmsg);
    }
    // Now copy back into passed in array...
    for (size_t i=0; i < count; i++) {
      local_minmax[i] = maxout[i];
    }
  }
#endif
}

template void Ioss::ParallelUtils::global_array_minmax(std::vector<int>&,     MinMax) const;
template void Ioss::ParallelUtils::global_array_minmax(std::vector<int64_t>&, MinMax) const;
template void Ioss::ParallelUtils::global_array_minmax(std::vector<double>&,  MinMax) const;

template <typename T>
void Ioss::ParallelUtils::global_array_minmax(std::vector<T> &local_minmax,  MinMax which) const
{
  if (!local_minmax.empty())
    global_array_minmax(&local_minmax[0], local_minmax.size(), which);
}

template void Ioss::ParallelUtils::gather(int, std::vector<int>&) const;
template void Ioss::ParallelUtils::gather(int64_t, std::vector<int64_t>&) const;

template <typename T>
void Ioss::ParallelUtils::gather(T my_value, std::vector<T> &result) const
{
  if (parallel_rank() == 0) {
    result.resize(parallel_size());
  }
#ifdef HAVE_MPI
  if (parallel_size() > 1) {
    const int success = MPI_Gather((void*)&my_value,  1, mpi_type(T(1)),
				   (void*)&result[0], 1, mpi_type(T(1)),
				   0, communicator_);
    if (success !=  MPI_SUCCESS) {
      std::ostringstream errmsg;
      errmsg << "Ioss::ParallelUtils::gather - MPI_Gather failed";
      IOSS_ERROR(errmsg);
    }
  } else {
    result[0] = my_value;
  }
#else
  result[0] = my_value;
#endif
}

template void Ioss::ParallelUtils::gather(std::vector<int> &my_values, std::vector<int> &result) const;
template void Ioss::ParallelUtils::gather(std::vector<int64_t> &my_values, std::vector<int64_t> &result) const;
template <typename T>
void Ioss::ParallelUtils::gather(std::vector<T> &my_values, std::vector<T> &result) const
{
  size_t count = my_values.size();
  if (parallel_rank() == 0) {
    result.resize(count * parallel_size());
  }
#ifdef HAVE_MPI
  if (parallel_size() > 1) {
    const int success = MPI_Gather((void*)TOPTR(my_values), count, mpi_type(T(1)),
				   (void*)TOPTR(result),    count, mpi_type(T(1)),
				   0, communicator_);
    if (success !=  MPI_SUCCESS) {
      std::ostringstream errmsg;
      errmsg << "Ioss::ParallelUtils::gather - MPI_Gather failed";
      IOSS_ERROR(errmsg);
    }
  } else {
    std::copy(my_values.begin(), my_values.end(), result.begin());
  }
#else
  std::copy(my_values.begin(), my_values.end(), result.begin());
#endif
}

