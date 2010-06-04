/*--------------------------------------------------------------------*/
/*    Copyright 2008 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef SIERRA_Ioss_ParallelUtils_h
#define SIERRA_Ioss_ParallelUtils_h

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
