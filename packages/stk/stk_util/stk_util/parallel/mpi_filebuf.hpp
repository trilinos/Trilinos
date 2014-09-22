// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Governement retains certain rights in this software.
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
// 
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
// 


#ifndef STK_UTIL_PARALLEL_mpi_filebuf_hpp
#define STK_UTIL_PARALLEL_mpi_filebuf_hpp

#include <stk_util/stk_config.h>

#if defined( STK_HAS_MPI )
#include <mpi.h>                        // for MPI_Comm, etc
#endif
#include <stddef.h>                     // for size_t
#include <cstdio>                       // for NULL, EOF, FILE
#include <ios>                          // for streambuf, ios_base, etc
#include <string>                       // for string, streamsize

//: Specialize the ANSI Standard C++ streambuf class
//: for a parallel file buffer.  The actual file is
//: only touched by the root processor.
//
//  READ MODE: The file is read on the root processor and
//  broadcast one buffer at a time to the remaining processors.
//
//  WRITE MODE: Each processor has a buffer that is locally
//  filled.  When the buffer is full on the root processor the
//  buffer is written to the output file.  When the buffer is
//  full on any other processor the size of the buffer is doubled.
//  The 'mpi_filebuf::flush' method gathers all buffers on the
//  root processor and writes the buffers to the output file.
//
//  GLOBAL: Calls to the 'open', 'flush', 'close', destructor,
//  and 'underflow' methods are global; these calls must be
//  made on all processors.  The 'underflow' method is called
//  by the 'istream' that uses the 'mpi_filebuf' object when
//  ever the input buffer is empty.  Thus reading from an 'istream'
//  that uses an 'mpi_filebuf' must be globally consistent.

class mpi_filebuf : public std::streambuf {
public:

  //: Construct an MPI-parallel input/output file buffer
  mpi_filebuf(bool aprepro=false, const std::string &aprepro_defines=std::string());

  //: GLOBAL: Open a file.
  // The file name is only significant on the root processsor.
  // May only be opened as ios::in, ios::out, ios:app.
  mpi_filebuf * open(
	  MPI_Comm       communicator ,       /* All processors */
    const int            root_processor ,     /* All processors */
    const std::ios_base::openmode file_mode , /* All processors */
    const char * const   file_name = NULL );  /* Root processor */

  //: GLOBAL: Close the file.
  // If output mode then flush the output.
  mpi_filebuf * close();

  //: GLOBAL: Flush the buffered output to the file.
  // Sends all buffers to the root processor,
  // write to the file, and flushes the file.
  mpi_filebuf * flush();

  //: GLOBAL: Destructor
  //  Close and then reclaim memory.
  virtual ~mpi_filebuf();

  //: Query if open, a local operations
  int is_open() const ;

  //: When the file buffer is in the 'closed' state set the buffer length,
  //  The input argument must be consistent on all processors; however,
  //  this condition is not checked until the next 'open' operation.
  mpi_filebuf * set_buffer_length( const size_t buffer_length );

  //: Query the current buffer
  void get_buffer( const char * & , size_t & ) const ;

  //: Query wall-clock time spent communicating.
  double wtime() const ;

protected:

  //: Called to refill the input buffer
  virtual int underflow();

  //: Called when output buffer is filled
  virtual int overflow( int c = EOF );

  //: Sync is a no-op
  virtual int sync();

  //: Setbuf is a no-op
  virtual std::streambuf * setbuf( char * s , std::streamsize n );

private:

  mpi_filebuf( const mpi_filebuf & ); // Not allowed
  mpi_filebuf & operator = ( const mpi_filebuf & ); // Not allowed

  MPI_Comm	comm ;            // Communicator
  int	        comm_root ;       // Rank of root processor
  std::FILE *   comm_root_fp ;    // Root processor's file
  int	        comm_output ;     // Output file
  char *	comm_buffer ;     // local buffer
  size_t	comm_buffer_len ; // length of buffer
  double	comm_time ;       // wall-time spent communicating

  bool    use_aprepro;        // If true, process through aprepro (input only)
  char *  aprepro_buffer;     // Buffer holding results of aprepro processing input file (root only)
  size_t  aprepro_buffer_len; // Length of aprepro buffer.
  size_t  aprepro_buffer_ptr; // Pointer to current location of data returned from aprepro_buffer.
  const std::string aprepro_defines;
};

/*--------------------------------------------------------------------*/

inline int  mpi_filebuf::is_open() const { return NULL != comm_buffer ; }

/* The SUN has the 'streambuf::pptr()' as a non-const method,
   which violates the ISO/ANSI standard specification.
   Therefore, must cast away the const. */

inline void mpi_filebuf::get_buffer( const char * & b , size_t & n ) const
{ b = comm_buffer ; n = (const_cast<mpi_filebuf*>(this))->pptr() - comm_buffer ; }

inline double mpi_filebuf::wtime() const { return comm_time ; }

#endif // STK_UTIL_PARALLEL_mpi_filebuf_hpp
