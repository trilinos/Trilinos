/*--------------------------------------------------------------------*/
/*    Copyright 2008 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <Ioss_ParallelUtils.h>
#include <Ioss_Utils.h>

#include <sstream>
#include <cstring>
#include <cstdlib>
#include <assert.h>

#include <algorithm>

#include <string>
#ifndef NO_MPI
#include <mpi.h>
#endif

#include <Ioss_SerializeIO.h>

#include <iostream>
#include <fstream>

Ioss::ParallelUtils::ParallelUtils(MPI_Comm the_communicator)
  : communicator_(the_communicator)
{}
  
bool Ioss::ParallelUtils::get_environment(const std::string &name, std::string &value, bool sync_parallel) const
{
#ifndef NO_MPI
  char *result_string = NULL;
  char *broadcast_string = NULL;
  int string_length = 0;

  int rank = Ioss::ParallelUtils::parallel_rank();
  if (rank == 0) {
    result_string = std::getenv(name.c_str());
    string_length = result_string ? (int)std::strlen(result_string) : 0;
  }

  if (sync_parallel) {
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
#ifndef NO_MPI
  char *result_string = NULL;
  int string_length = 0;

  int rank = Ioss::ParallelUtils::parallel_rank();
  if (rank == 0) {
    result_string = std::getenv(name.c_str());
    string_length = result_string ? (int)std::strlen(result_string) : 0;
  }

  if (sync_parallel)
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
#ifndef NO_MPI
  MPI_Comm_size(communicator_, &my_size);
#endif
  return my_size;
}

int Ioss::ParallelUtils::parallel_rank() const
{
  int my_rank = 0;
#ifndef NO_MPI
  MPI_Comm_rank(communicator_, &my_rank);
#endif
  return my_rank;
}

void Ioss::ParallelUtils::attribute_reduction( const int length , char buffer[]) const
{
#ifndef NO_MPI
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
#ifndef NO_MPI
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

double Ioss::ParallelUtils::global_minmax(double local_minmax, Ioss::ParallelUtils::MinMax which) const
{
  double minmax = local_minmax;

#ifndef NO_MPI
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

#ifndef NO_MPI
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

#ifndef NO_MPI
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

void Ioss::ParallelUtils::global_array_minmax(int *local_minmax, size_t count, Ioss::ParallelUtils::MinMax which) const
{
#ifndef NO_MPI
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

void Ioss::ParallelUtils::global_array_minmax(unsigned int *local_minmax, size_t count,
				      Ioss::ParallelUtils::MinMax which) const
{
#ifndef NO_MPI
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
#ifndef NO_MPI
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
#ifndef NO_MPI
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

void Ioss::ParallelUtils::gather(std::vector<int> &my_values, std::vector<int> &result) const
{
  size_t count = my_values.size();
  if (parallel_rank() == 0) {
    result.resize(count * parallel_size());
  }
#ifndef NO_MPI
  if (parallel_size() > 1) {
    const int success = MPI_Gather((void*)&my_values[0],  count, MPI_INT,
				   (void*)&result[0], count, MPI_INT,
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


