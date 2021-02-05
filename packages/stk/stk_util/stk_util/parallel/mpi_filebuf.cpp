/*
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
 */

#include "stk_util/parallel/mpi_filebuf.hpp"
#include "stk_util/stk_config.h"  // for STK_HAS_MPI
#include <cassert>                // for assert
#include <fstream>                // for fstream
#include <cstddef>                // for size_t
#include <cstdlib>                // for free, malloc, realloc
#include <cstring>                // for memcpy
#include <iostream>               // for operator<<, basic_ostream, cerr, endl, basic_istream
#include <vector>                 // for vector

#if !defined(NOT_HAVE_STK_SEACASAPREPRO_LIB)
#include "apr_tokenize.h"         // for tokenize
#include "aprepro.h"              // for Aprepro, aprepro_options
#endif

enum { buffer_default_length = 4096 };
enum { buffer_putback_length =   16 };

// --------------------------------------------------------------------
namespace {

  double mpi_wall_time()
  {
#if defined(STK_HAS_MPI)
    return MPI_Wtime();
#else
    return 0.0;
#endif
  }
}
// --------------------------------------------------------------------

mpi_filebuf::mpi_filebuf(bool useAprepro, const std::string &apreproDefines)
  : std::streambuf(),
    comm( MPI_COMM_NULL ),
    comm_root( -1 ),
    comm_root_fp( nullptr ),
    comm_output( 0 ),
    comm_buffer( nullptr ),
    comm_buffer_len( buffer_default_length ),
    comm_time(0.0),
    use_aprepro(useAprepro),
    aprepro_buffer(nullptr),
    aprepro_buffer_len(0),
    aprepro_buffer_ptr(0),
    aprepro_parsing_error_count(0),
    aprepro_parsing_warning_count(0),
    aprepro_defines(apreproDefines)
{}

mpi_filebuf::~mpi_filebuf()
{
  close();
}

/*--------------------------------------------------------------------*/

mpi_filebuf * mpi_filebuf::set_buffer_length( const size_t len )
{
  // If already open then abort
  if ( nullptr != comm_buffer ) return nullptr ;

  // Wait and verify upon the attempt to open
  comm_buffer_len = static_cast<size_t>(buffer_putback_length) < len ? len : static_cast<size_t>(buffer_putback_length);

  return this ;
}

/*--------------------------------------------------------------------*/

#if defined( STK_HAS_MPI )
mpi_filebuf * mpi_filebuf::open(
	MPI_Comm       communicator ,
  const int            root_processor ,
  const std::ios_base::openmode file_mode ,
  const char * const   file_name )
{
#if !defined(NOT_HAVE_STK_SEACASAPREPRO_LIB)
  const double start_time = mpi_wall_time();

  // If already open then abort
  if ( nullptr != comm_buffer ) return nullptr ;

  const int mode =
    ( std::ios::in  == file_mode ) ? 'r' : (
    ( std::ios::out == file_mode ) ? 'w' : (
    ( std::ios::app == file_mode ) ? 'a' : -1 ) );


  int rank = 0 ;
  int err = 0 ;
  int local = 0, global = 0 ;
  int data[3] ;

  // Broadcast the selected root processor and 'C' file mode

  data[0] = root_processor ;
  data[1] = mode ;
  data[2] = comm_buffer_len ;

  if ( MPI_SUCCESS != ( err = MPI_Bcast(data,3,MPI_INT,0,communicator) ) )
    MPI_Abort( communicator , err );

  // Verify that all processors have the same root, mode, and buffer length:

  local = data[0] != root_processor || data[1] != mode || data[2] != static_cast<signed>(comm_buffer_len) ;

  if ( MPI_SUCCESS != ( err =
       MPI_Allreduce(&local,&global,1,MPI_INT,MPI_BOR,communicator) ) )
    MPI_Abort( communicator , err );

  if ( global ) {
    comm_time += mpi_wall_time() - start_time ;
    return nullptr ;
  }

  //--------------------------------------------------------------------
  // Root processor and mode are consistent.
  // All processors try to allocate buffers and the
  // root processor tries to open the file.

  if ( MPI_SUCCESS != ( err =  MPI_Comm_rank( communicator , &rank ) ) )
    MPI_Abort( communicator , err );

  char * const tmp_buf = static_cast<char*>(std::malloc( comm_buffer_len ));
  std::FILE *       tmp_fp  = nullptr ;

  local = tmp_buf == nullptr ; // Failed allocation ?

  if ( root_processor == rank && ! local ) {
    tmp_fp = std::fopen( file_name , ( ( ( mode == 'r' ) ? "r" :
				    ( mode == 'w' ) ? "w" : "a" ) ) );
    local = nullptr == tmp_fp ;
  }

  if ( MPI_SUCCESS != ( err =
       MPI_Allreduce(&local,&global,1,MPI_INT,MPI_BOR,communicator) ) )
    MPI_Abort( communicator , err );

  if ( global ) {
    if ( nullptr != tmp_buf ) std::free(   tmp_buf ); // Deallocate
    if ( nullptr != tmp_fp  ) std::fclose( tmp_fp );  // Close the file
    comm_time += mpi_wall_time() - start_time ;
    return nullptr ;
  }

  // If input and use_aprepro, parse the file and store parsed results
  // into buffer on root_processor.
  if (use_aprepro && !comm_output) {

    std::stringstream errStream;
    std::stringstream warnStream;
    std::stringstream infoStream;

    if (root_processor == rank) {
      // Note that file is double-opened.  Aprepro uses an std::fstream
      std::fstream infile(file_name, std::fstream::in);
      if (!infile.good()) {
	if ( nullptr != tmp_buf ) std::free(   tmp_buf ); // Deallocate
	if ( nullptr != tmp_fp  ) std::fclose( tmp_fp );  // Close the file
	std::cerr << "APREPRO: Could not open file: " << file_name << std::endl;
	return nullptr;
      }

      SEAMS::Aprepro aprepro;

      aprepro.set_error_streams(&errStream, &warnStream, &infoStream);

      stk::add_aprepro_defines(aprepro, aprepro_defines);
      aprepro.ap_options.require_defined=true;


      bool result = aprepro.parse_stream(infile);

      if (result) {
	// Get size of buffer needed to store the parsed data...
	std::string tmp = aprepro.parsing_results().str();
	aprepro.clear_results();
	aprepro_buffer_len = tmp.size();
	aprepro_buffer_ptr = 0; // At beginning of buffer...
	aprepro_buffer = static_cast<char*>(std::malloc(aprepro_buffer_len));
	std::memcpy(aprepro_buffer, tmp.data(), aprepro_buffer_len);
      }
      aprepro_parsing_error_count = aprepro.get_error_count();
      aprepro_parsing_warning_count = aprepro.get_warning_count();
    }
    err = MPI_Bcast(&aprepro_parsing_error_count, 1, MPI_INT, root_processor, communicator );
    err = MPI_Bcast(&aprepro_parsing_warning_count, 1, MPI_INT, root_processor, communicator );

    aprepro_warnings = warnStream.str();
    aprepro_errors   = errStream.str();
    aprepro_info     = infoStream.str();     

    if(aprepro_parsing_error_count > 0 && aprepro_errors.empty()) {
        aprepro_errors = "Fatal aprepro errors encountered";
    }


    if (err != MPI_SUCCESS)
      MPI_Abort(communicator,err);
  }

  //--------------------------------------------------------------------
  // All memory allocated and root processor opened the file
  // Update the internal members accordingly.

  comm         = communicator ;
  comm_root    = root_processor ;
  comm_root_fp = tmp_fp ;
  comm_buffer  = tmp_buf ;
  comm_output  = mode != 'r' ;

  // If output then set up put-buffer
  if ( comm_output ) setp( comm_buffer, comm_buffer + comm_buffer_len );

  comm_time += mpi_wall_time() - start_time ;

#endif
  return this ;
}
#else
mpi_filebuf * mpi_filebuf::open(
				      MPI_Comm       communicator,
				const int            root_processor ,
				const std::ios_base::openmode file_mode ,
				const char * const   file_name )
{
  // If already open then abort
  if ( nullptr != comm_buffer ) return nullptr ;

  const int mode =
    ( std::ios::in  == file_mode ) ? 'r' : (
    ( std::ios::out == file_mode ) ? 'w' : (
    ( std::ios::app == file_mode ) ? 'a' : -1 ) );

  char * const tmp_buf = static_cast<char*>(std::malloc( comm_buffer_len ));
  std::FILE *       tmp_fp  = std::fopen( file_name , ( ( ( mode == 'r' ) ? "r" :
							  ( mode == 'w' ) ? "w" : "a" ) ) );
#if !defined(NOT_HAVE_STK_SEACASAPREPRO_LIB)
  // If input and use_aprepro, parse the file and store parsed results
  // into buffer on root_processor.
  if (use_aprepro && !comm_output) {
    // Note that file is double-opened.  Aprepro uses an std::fstream
    std::fstream infile(file_name, std::fstream::in);
    if (!infile.good()) {
      if ( nullptr != tmp_buf ) std::free(   tmp_buf ); // Deallocate
      if ( nullptr != tmp_fp  ) std::fclose( tmp_fp );  // Close the file
      std::cerr << "APREPRO: Could not open file: " << file_name << std::endl;
      return nullptr;
    }

    SEAMS::Aprepro aprepro;

    stk::add_aprepro_defines(aprepro, aprepro_defines);

    aprepro.ap_options.require_defined=true;

    bool result = aprepro.parse_stream(infile);
    if (result) {
      // Get size of buffer needed to store the parsed data...
      std::string tmp = aprepro.parsing_results().str();
      aprepro.clear_results();
      aprepro_buffer_len = tmp.size();
      aprepro_buffer_ptr = 0; // At beginning of buffer...
      aprepro_buffer = static_cast<char*>(std::malloc(aprepro_buffer_len));
      std::memcpy(aprepro_buffer, tmp.data(), aprepro_buffer_len);
    }
    aprepro_parsing_error_count = aprepro.get_error_count();
  }
#endif

  //--------------------------------------------------------------------
  // All memory allocated and root processor opened the file
  // Update the internal members accordingly.

  comm         = communicator ;
  comm_root    = root_processor ;
  comm_root_fp = tmp_fp ;
  comm_buffer  = tmp_buf ;
  comm_output  = mode != 'r' ;

  // If output then set up put-buffer
  if ( comm_output ) setp( comm_buffer, comm_buffer + comm_buffer_len );

  return this ;
}
#endif
/*--------------------------------------------------------------------*/

mpi_filebuf * mpi_filebuf::close()
{
  mpi_filebuf * tmp = nullptr ;

  if ( nullptr != comm_buffer ) {

    flush(); // Flush the buffers

    if ( nullptr != comm_root_fp ) std::fclose( comm_root_fp ); // Close the file

    std::free( comm_buffer ); // Free the buffer

    if ( comm_output ) setp(nullptr,nullptr);
    else               setg(nullptr,nullptr,nullptr);

    // Reset the members:

    comm         = MPI_COMM_NULL ;
    comm_root    = -1 ;
    comm_root_fp = nullptr ;
    comm_output  = 0 ;
    comm_buffer  = nullptr ;

    if (aprepro_buffer != nullptr) {
      std::free(aprepro_buffer);
      aprepro_buffer = nullptr;
      aprepro_buffer_ptr = 0;
      aprepro_buffer_len = 0;
    }
    
    tmp = this ;
  }

  return tmp ;
}

/*--------------------------------------------------------------------*/
/* Underflow, a global call.
   Read more data from the root processor's file and
   broadcast it to all processors.
*/

int mpi_filebuf::underflow()
{
  const double start_time = mpi_wall_time();

  if ( nullptr != comm_buffer && ! comm_output &&
       (gptr() == nullptr || gptr() >= egptr()) ) { // valid get buffer
    // Length of the buffer, consistent on all processors
    // Entire buffer is offset to accomodate putbacks
    const size_t size = comm_buffer_len - buffer_putback_length ;
    char * const buf  = comm_buffer     + buffer_putback_length ;

    int nread = 0;
    if (comm_root_fp != nullptr) {
      if (use_aprepro) {
	// Copy from current location in aprepro_buffer into comm_buffer
	nread = size;
	if (aprepro_buffer_ptr + size > aprepro_buffer_len) {
	  nread = aprepro_buffer_len - aprepro_buffer_ptr;
	}
	std::memcpy(buf, aprepro_buffer+aprepro_buffer_ptr, nread);
	aprepro_buffer_ptr += nread;
      } else {
	// Root processor reads from the file and broadcasts the result
	nread = std::fread(buf,1,size,comm_root_fp);
      }
    }

#if defined( STK_HAS_MPI )
    int err = MPI_Bcast(&nread, 1, MPI_INT, comm_root, comm );
    if (err != MPI_SUCCESS)
      MPI_Abort(comm,err);
#endif
    
    // If the read is successfull then update the get buffer pointers:
    if (nread > 0) {
#if defined( STK_HAS_MPI )
      // Broadcast the read buffer to all processors:
      err = MPI_Bcast(buf, nread, MPI_BYTE, comm_root, comm);
      if (err != MPI_SUCCESS)
	MPI_Abort(comm,err);
#endif
      // Set the get buffer:
      setg( comm_buffer, buf, buf + nread );

      // Return the next character from the file:
      comm_time += mpi_wall_time() - start_time ;

      return *buf ;
    }
  }

  // Failed: set the get buffer to NULL and return EOF
  setg(nullptr, nullptr, nullptr);

  comm_time += mpi_wall_time() - start_time ;

  return EOF;
}

/*--------------------------------------------------------------------*/
/* Overflow, a local call.
    Output complete lines of data on the root processor.
    Increase the buffer size on all other processors.
*/

int mpi_filebuf::overflow( int c )
{
  if ( nullptr != comm_buffer && comm_output ) { // open for write

    // Determine current offset and length:
    char * cur_buffer = comm_buffer ;
    size_t cur_offset = pptr()  - cur_buffer ;
    size_t cur_length = epptr() - cur_buffer ;

    assert( cur_offset <= cur_length /* detecting abuse by 'ostream' */ );

    if ( nullptr != comm_root_fp ) {
      if ( std::fwrite(cur_buffer,1,cur_offset,comm_root_fp) != cur_offset ) {
	return EOF ; // Write failed
      }
      cur_offset = 0 ;
    }
    else if ( cur_length <= cur_offset ) {
      // Not root processor, ran out of buffer space and
      // cannot write so increase the buffer size:
      cur_buffer = static_cast<char*>(std::realloc( cur_buffer , cur_length *= 2 ));
    }

    // If buffer is still good then reset the put-buffer

    if ( nullptr != cur_buffer ) {

      comm_buffer = cur_buffer ;

      setp( cur_buffer + cur_offset, cur_buffer + cur_length );

      if ( c != EOF ) {

	sputc(c);
	return c;
      }
      else {
	return 0;
      }
    }
  }
  return EOF ;
}

/*--------------------------------------------------------------------*/
/* Send output buffers to root processor and
   write them to the output file.
*/

#if defined( STK_HAS_MPI )
mpi_filebuf * mpi_filebuf::flush()
{
  const double start_time = mpi_wall_time();

  int result = -1 ; // Failure return value

  if ( nullptr != comm_buffer && comm_output ) { // Open for write

    int err = 0 ;

    result = 0 ;

    // Determine the local length:

    char * cur_buf = comm_buffer ;
    unsigned int cur_len = pptr() - cur_buf ;

    // Determine the global lengths

    char * recv_buf  = nullptr ;
    int  * recv_len  = nullptr ;
    int  * recv_disp = nullptr ;

    int nproc = 1 ;


//  if ( nullptr != comm_root_fp ) {

//  It should not be neccessary to allocate recv_len on non-root
//  nodes, but the MPI_Gatherv on Janus always accesses recv_len
//  even on non-root processors which causes a segmentaion
//  violation if recv_len is set to nullptr.

    if ( MPI_SUCCESS != ( err = MPI_Comm_size(comm,&nproc) ) )
      MPI_Abort( comm , err );
    recv_len = static_cast<int*>(std::malloc( sizeof(int) * nproc ));

    if ( nullptr == recv_len ) MPI_Abort( comm , MPI_ERR_UNKNOWN );

    for (int j = 0 ; j < nproc ; ++j )
      recv_len[j] = 0;
//  }

    // Gather buffer lengths on the root processor

    if ( MPI_SUCCESS != ( err =
	 MPI_Gather(&cur_len,1,MPI_INT,recv_len,1,MPI_INT,comm_root,comm)))
      MPI_Abort( comm , err );

    // Root processor must allocate enough buffer space:

    if ( nullptr != comm_root_fp ) {

      recv_len[ comm_root ] = 0 ; // Don't send to self

      if ( nullptr == ( recv_disp = static_cast<int*>(std::malloc( sizeof(int) * (nproc + 1) )) ) )
	result = -1 ;

      if ( 0 == result ) { // Allocation succeeded

	recv_disp[0] = 0 ;

	for (int i = 0 ; i < nproc ; ++i )
	  recv_disp[i+1] = recv_disp[i] + recv_len[i] ;

	if ( 0 < recv_disp[nproc] ) {
	  if ( nullptr == ( recv_buf = static_cast<char*>(std::malloc( recv_disp[nproc] ) ) ))
	    result = -1 ;
	}
	else {
	  result = 1 ; // No need to gather!
	}

	if ( -1 != result ) {

	  // Write the root processor's buffer

	  if ( 0 < cur_len ) {
	    if ( std::fwrite(cur_buf,1,cur_len,comm_root_fp) != cur_len )
	      result = -1 ; // Write failed

	    cur_len = 0 ; // Wrote this buffer
	  }
	}
      }
      std::fflush( comm_root_fp );
    }

    // Root process broadcasts that all is well with the allocation

    if ( MPI_SUCCESS != ( err = MPI_Bcast(&result,1,MPI_INT,comm_root,comm)))
      MPI_Abort( comm , err );

    if ( 0 == result ) { // All-is-well, need to gather and write

      // Gather the buffers to the root processor

      if ( MPI_SUCCESS != ( err =
	   MPI_Gatherv(cur_buf,  cur_len,             MPI_BYTE,
		       recv_buf, recv_len, recv_disp, MPI_BYTE,
		       comm_root, comm ) ) )
	MPI_Abort( comm , err );

       // Output the buffers, beginning with 'comm_root'

      if ( nullptr != comm_root_fp ) {

	for (int i = 1 ; i < nproc && 0 == result ; ++i ) {
	  const int j   = ( i + comm_root ) % nproc ;
	  const unsigned int len = recv_len[j] ;

	  if ( 0 < len )
	    if ( std::fwrite(recv_buf+recv_disp[j],1,len,comm_root_fp) != len )
	      result = -1 ; // Write failed
	}

	std::fflush( comm_root_fp );
      }

      // Broadcast that the write succeeded

      if ( MPI_SUCCESS != ( err = MPI_Bcast(&result,1,MPI_INT,comm_root,comm)))
	MPI_Abort( comm , err );
    }
    else if ( 1 == result ) {
      // Did not need to gather

      result = 0 ;
    }

    // Reset the output buffer

    setp( comm_buffer , epptr() );

    // Clean up allocated memory

    if ( nullptr != recv_buf  ) std::free( recv_buf );
    if ( nullptr != recv_len  ) std::free( recv_len );
    if ( nullptr != recv_disp ) std::free( recv_disp );
  }

  comm_time += mpi_wall_time() - start_time ;

  return -1 == result ? nullptr : this ;
}
#else
mpi_filebuf * mpi_filebuf::flush()
{
  sync();
  return this;
}
#endif
/*--------------------------------------------------------------------*/

int mpi_filebuf::sync()
{
  // The root processor will push to file, all others ignore

  if ( nullptr != comm_root_fp ) {

    // Determine the local length:

    char * cur_buf = comm_buffer ;
    int    cur_len = pptr() - cur_buf ;

    if ( 0 < cur_len ) std::fwrite(cur_buf,1,cur_len,comm_root_fp);

    std::fflush( comm_root_fp );

    setp( comm_buffer , epptr() );
  }

  return 0 ;
}

#if !defined(NOT_HAVE_STK_SEACASAPREPRO_LIB)
namespace stk
{
  void add_aprepro_defines(SEAMS::Aprepro &aprepro, const std::string &defines)
  {
    // See if any variables were defined on the command line...
    if (!defines.empty()) {
      // Defines are space/comma-separated pairs of the form 'var=value'
      // or "options" of the form '-W' or '--warning'
      // Split the string and then process each variable...
      std::vector<std::string> tokens = SEAMS::tokenize(defines, " ,\t");
      for (size_t i=0; i<tokens.size(); i++) {
	std::string token = tokens[i];
	if (token[0] == '-') {
	  aprepro.set_option(token);
	}
	else {
	  // It is an aprepro variable definition
	  std::vector<std::string> define = SEAMS::tokenize(token, "=");
	  if (define.size() == 2) {
	    // Determine whether the define is string type or double/int...
	    bool immutable = define[0][0] != '_';
	    std::stringstream ss(define[1]);
	    double d = 0;
	    ss >> d;
	    if (ss.fail() || !ss.eof()) {
	      // Not a valid number; treat as a string
	      aprepro.add_variable(define[0], define[1], immutable);
	    } else {
	      aprepro.add_variable(define[0], d, immutable);
	    }
	  } else {
	    std::cerr << "APREPRO: Invalid format for predefined variable: '" << token << "'\n"
		      << "         Required format is 'var=value' or 'var=\"value\"'\n";
	  }
	}
      }
    }
  }
}
#endif
