/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <string>

#include <stk_util/parallel/BroadcastArg.hpp>

namespace stk_classic {

BroadcastArg::BroadcastArg(
  stk_classic::ParallelMachine  parallel_machine,
  int                   argc,
  char **               argv)
{
#ifdef STK_HAS_MPI
  int rank = stk_classic::parallel_machine_rank(parallel_machine);
#else
  int rank = 0;
#endif

  size_t buffer_length = 0;
  char * buffer = 0;
  
// Populate m_argc, m_buffer and buffer_length on rank 0 or !STK_HAS_MPI
  if (rank == 0) {
    m_argc = argc;

    std::string s;
    for (int i = 0; i < argc; ++i) {
      s += argv[i];
      s += '\0';
    }
    
    buffer_length = s.size();
    buffer = new char[buffer_length];
    
    std::copy(s.begin(), s.end(), buffer);
  }

// if STK_HAS_MPI, broadcast m_argc, buffer and buffer_length to processors
#ifdef STK_HAS_MPI
  if (rank == 0) {
    int lengths_buffer[2];
    lengths_buffer[0] = m_argc;
    lengths_buffer[1] = buffer_length;
    
    MPI_Bcast(lengths_buffer, 2, MPI_INT, 0, parallel_machine);

    MPI_Bcast(buffer, buffer_length, MPI_BYTE, 0, parallel_machine);
  }
  else {
    int lengths_buffer[2];
    MPI_Bcast(lengths_buffer, 2, MPI_INT, 0, parallel_machine);

    m_argc = lengths_buffer[0];
    buffer_length = lengths_buffer[1];
    buffer = new char[buffer_length];
    
    MPI_Bcast(buffer, buffer_length, MPI_BYTE, 0, parallel_machine); 
  }
#endif
    
// Populate the m_argv
  m_argv = new char *[m_argc];
  
// argv[0] will always point to buffer, so argv[0] needs to be deleted by the destructor
  char *c = &buffer[0];
  for (int i = 0; i < argc; ++i) {
    m_argv[i] = c;
    while (*c)
      ++c;
    ++c;
  }
}


BroadcastArg::~BroadcastArg()
{
  delete[] m_argv[0];
  delete[] m_argv;
}

} // namespace stk_classic
