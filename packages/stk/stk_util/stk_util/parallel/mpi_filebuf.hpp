// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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


#ifndef STK_UTIL_PARALLEL_MPI_FILEBUF_HPP
#define STK_UTIL_PARALLEL_MPI_FILEBUF_HPP

#include "stk_util/parallel/Parallel.hpp"  // for MPI_Comm, ompi_communicator_t
#include <cstdio>                          // for size_t, EOF, FILE
#include <ios>                             // for basic_streambuf<>::char_type, ios_base, ios_ba...
#include <streambuf>
#include <string>                          // for string

namespace SEAMS { class Aprepro; }

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
    const char * const   file_name = nullptr );  /* Root processor */

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

  //: Return number of aprepro parse errors detected
  int aprepro_parse_error_count() const {return aprepro_parsing_error_count;}
  int aprepro_parse_warning_count() const {return aprepro_parsing_warning_count;}
  const std::string& aprepro_parse_errors() const {return aprepro_errors;}
  const std::string& aprepro_parse_warnings() const {return aprepro_warnings;}
  const std::string& aprepro_parse_info() const {return aprepro_info;}
protected:

  //: Called to refill the input buffer
  virtual int underflow() override;

  //: Called when output buffer is filled
  virtual int overflow( int c = EOF ) override;

  //: Sync is a no-op
  virtual int sync() override;

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
  std::string aprepro_warnings;  // (root only)
  std::string aprepro_errors; // (root only)
  std::string aprepro_info; // (root only)
  size_t  aprepro_buffer_len; // Length of aprepro buffer.
  size_t  aprepro_buffer_ptr; // Pointer to current location of data returned from aprepro_buffer.
  int     aprepro_parsing_error_count; 
  int     aprepro_parsing_warning_count; 
  const std::string aprepro_defines;


};

/*--------------------------------------------------------------------*/

inline int  mpi_filebuf::is_open() const { return nullptr != comm_buffer ; }

inline void mpi_filebuf::get_buffer( const char * & b , size_t & n ) const
{ b = comm_buffer ; n = (const_cast<mpi_filebuf*>(this))->pptr() - comm_buffer ; }

inline double mpi_filebuf::wtime() const { return comm_time ; }

#if !defined(NOT_HAVE_STK_SEACASAPREPRO_LIB)
namespace stk
{
  void add_aprepro_defines(SEAMS::Aprepro &aprepro, const std::string &defines);
}
#endif

#endif // STK_UTIL_PARALLEL_MPI_FILEBUF_HPP
