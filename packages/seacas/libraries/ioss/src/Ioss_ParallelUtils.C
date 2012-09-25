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

double Ioss::ParallelUtils::global_minmax(double local_minmax, Ioss::ParallelUtils::MinMax which) const
{
  double minmax = local_minmax;

#ifdef HAVE_MPI
  if (parallel_size() > 1) {
    if (Ioss::SerializeIO::isEnabled() && Ioss::SerializeIO::inBarrier()) {
      std::ostringstream errmsg;
      errmsg << "Attempting mpi while in barrier owned by " << Ioss::SerializeIO::getOwner();
      IOSS_ERROR(errmsg);
    }
    static double inbuf[1], outbuf[1];
    inbuf[0] = local_minmax;

    MPI_Op oper;
    if (which == DO_MAX)
      oper = MPI_MAX;
    else
      oper = MPI_MIN;

    const int success = MPI_Allreduce((void*)&inbuf[0], &outbuf[0], 1,
				      MPI_DOUBLE, oper, communicator_);
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

int Ioss::ParallelUtils::global_minmax(int local_minmax, Ioss::ParallelUtils::MinMax which) const
{
  int minmax = local_minmax;

#ifdef HAVE_MPI
  if (parallel_size() > 1) {
    if (Ioss::SerializeIO::isEnabled() && Ioss::SerializeIO::inBarrier()) {
      std::ostringstream errmsg;
      errmsg << "Attempting mpi while in barrier owned by " << Ioss::SerializeIO::getOwner();
      IOSS_ERROR(errmsg);
    }
    static int inbuf[1], outbuf[1];
    inbuf[0] = local_minmax;

    MPI_Op oper;
    if (which == DO_MAX)
      oper = MPI_MAX;
    else
      oper = MPI_MIN;

    const int success = MPI_Allreduce((void*)&inbuf[0], &outbuf[0], 1,
				      MPI_INT, oper, communicator_);
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

unsigned int Ioss::ParallelUtils::global_minmax(unsigned int local_minmax, Ioss::ParallelUtils::MinMax which) const
{
  unsigned int minmax = local_minmax;

#ifdef HAVE_MPI
  if (parallel_size() > 1) {
    if (Ioss::SerializeIO::isEnabled() && Ioss::SerializeIO::inBarrier()) {
      std::ostringstream errmsg;
      errmsg << "Attempting mpi while in barrier owned by " << Ioss::SerializeIO::getOwner();
      IOSS_ERROR(errmsg);
    }
    static unsigned int inbuf[1], outbuf[1];
    inbuf[0] = local_minmax;

    MPI_Op oper;
    if (which == DO_MAX)
      oper = MPI_MAX;
    else
      oper = MPI_MIN;

    const int success = MPI_Allreduce((void*)&inbuf[0], &outbuf[0], 1,
				      MPI_UNSIGNED, oper, communicator_);
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

void Ioss::ParallelUtils::global_array_minmax(std::vector<int64_t> &local_minmax,  MinMax which) const
{
  if (!local_minmax.empty())
    global_array_minmax(&local_minmax[0], local_minmax.size(), which);
}

void Ioss::ParallelUtils::global_array_minmax(std::vector<int> &local_minmax,  MinMax which) const
{
  if (!local_minmax.empty())
    global_array_minmax(&local_minmax[0], local_minmax.size(), which);
}

void Ioss::ParallelUtils::global_array_minmax(std::vector<double> &local_minmax,  MinMax which) const
{
  if (!local_minmax.empty())
    global_array_minmax(&local_minmax[0], local_minmax.size(), which);
}

void Ioss::ParallelUtils::global_array_minmax(std::vector<unsigned int> &local_minmax,  MinMax which) const
{
  if (!local_minmax.empty())
    global_array_minmax(&local_minmax[0], local_minmax.size(), which);
}

void Ioss::ParallelUtils::global_array_minmax(int *local_minmax, size_t count, Ioss::ParallelUtils::MinMax which) const
{
#ifdef HAVE_MPI
  if (parallel_size() > 1 && count > 0) {
    if (Ioss::SerializeIO::isEnabled() && Ioss::SerializeIO::inBarrier()) {
      std::ostringstream errmsg;
      errmsg << "Attempting mpi while in barrier owned by " << Ioss::SerializeIO::getOwner();
      IOSS_ERROR(errmsg);
    }

    IntVector maxout(count);

    MPI_Op oper;
    if (which == DO_MAX)
      oper = MPI_MAX;
    else
      oper = MPI_MIN;

    const int success = MPI_Allreduce((void*)&local_minmax[0], &maxout[0],
				      static_cast<int>(count),
				      MPI_INT, oper, communicator_);
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

void Ioss::ParallelUtils::global_array_minmax(int64_t *local_minmax, size_t count, Ioss::ParallelUtils::MinMax which) const
{
#ifdef HAVE_MPI
  if (parallel_size() > 1 && count > 0) {
    if (Ioss::SerializeIO::isEnabled() && Ioss::SerializeIO::inBarrier()) {
      std::ostringstream errmsg;
      errmsg << "Attempting mpi while in barrier owned by " << Ioss::SerializeIO::getOwner();
      IOSS_ERROR(errmsg);
    }

    Int64Vector maxout(count);

    MPI_Op oper;
    if (which == DO_MAX)
      oper = MPI_MAX;
    else
      oper = MPI_MIN;

    const int success = MPI_Allreduce((void*)&local_minmax[0], &maxout[0],
				      static_cast<int64_t>(count),
				      MPI_LONG_LONG_INT, oper, communicator_);
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

void Ioss::ParallelUtils::global_array_minmax(unsigned int *local_minmax, size_t count,
				      Ioss::ParallelUtils::MinMax which) const
{
#ifdef HAVE_MPI
  if (parallel_size() > 1 && count > 0) {
    if (Ioss::SerializeIO::isEnabled() && Ioss::SerializeIO::inBarrier()) {
      std::ostringstream errmsg;
      errmsg << "Attempting mpi while in barrier owned by " << Ioss::SerializeIO::getOwner();
      IOSS_ERROR(errmsg);
    }

    std::vector<unsigned int> maxout(count);

    MPI_Op oper;
    if (which == DO_MAX)
      oper = MPI_MAX;
    else
      oper = MPI_MIN;

    const int success = MPI_Allreduce((void*)&local_minmax[0], &maxout[0],
				      static_cast<int>(count),
				      MPI_UNSIGNED, oper, communicator_);
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

void Ioss::ParallelUtils::global_array_minmax(double *local_minmax, size_t count, Ioss::ParallelUtils::MinMax which) const
{
#ifdef HAVE_MPI
  if (parallel_size() > 1 && count > 0) {
    if (Ioss::SerializeIO::isEnabled() && Ioss::SerializeIO::inBarrier()) {
      std::ostringstream errmsg;
      errmsg << "Attempting mpi while in barrier owned by " << Ioss::SerializeIO::getOwner();
      IOSS_ERROR(errmsg);
    }

    std::vector<double> maxout(count);

    MPI_Op oper;
    if (which == DO_MAX)
      oper = MPI_MAX;
    else
      oper = MPI_MIN;

    const int success = MPI_Allreduce((void*)&local_minmax[0], &maxout[0],
				      static_cast<int>(count),
				      MPI_DOUBLE, oper, communicator_);
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

void Ioss::ParallelUtils::gather(int my_value, std::vector<int> &result) const
{
  if (parallel_rank() == 0) {
    result.resize(parallel_size());
  }
#ifdef HAVE_MPI
  if (parallel_size() > 1) {
    const int success = MPI_Gather((void*)&my_value,  1, MPI_INT,
				   (void*)&result[0], 1, MPI_INT,
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

void Ioss::ParallelUtils::gather(int64_t my_value, std::vector<int64_t> &result) const
{
  if (parallel_rank() == 0) {
    result.resize(parallel_size());
  }
#ifdef HAVE_MPI
  if (parallel_size() > 1) {
    const int success = MPI_Gather((void*)&my_value,  1, MPI_LONG_LONG_INT,
				   (void*)&result[0], 1, MPI_LONG_LONG_INT,
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

void Ioss::ParallelUtils::gather(std::vector<int> &my_values, std::vector<int> &result) const
{
  size_t count = my_values.size();
  if (parallel_rank() == 0) {
    result.resize(count * parallel_size());
  }
#ifdef HAVE_MPI
  if (parallel_size() > 1) {
    const int success = MPI_Gather((void*)TOPTR(my_values),  count, MPI_INT,
				   (void*)TOPTR(result), count, MPI_INT,
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

void Ioss::ParallelUtils::gather(std::vector<int64_t> &my_values, std::vector<int64_t> &result) const
{
  size_t count = my_values.size();
  if (parallel_rank() == 0) {
    result.resize(count * parallel_size());
  }
#ifdef HAVE_MPI
  if (parallel_size() > 1) {
    const int success = MPI_Gather((void*)TOPTR(my_values),  count, MPI_LONG_LONG_INT,
				   (void*)TOPTR(result), count, MPI_LONG_LONG_INT,
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


