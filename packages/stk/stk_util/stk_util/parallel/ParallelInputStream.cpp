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

#include "stk_util/stk_config.h"        // for STK_HAS_MPI
#include <stk_util/parallel/ParallelInputStream.hpp>
#include <cstdio>                       // for NULL, EOF, fclose, fopen, etc
#include <cstring>                      // for memset
#include <iostream>                     // for cerr
#include <stdexcept>                    // for runtime_error
#include <string>                       // for string
#include "stk_util/parallel/Parallel.hpp"  // for ParallelMachine, etc

/*--------------------------------------------------------------------*/

namespace stk {
namespace {

#if defined( STK_HAS_MPI )

void broadcast( ParallelMachine comm , void * buf , int n )
{ MPI_Bcast( buf , n , MPI_BYTE , 0 , comm ); }

#else

void broadcast( ParallelMachine , void * , int ) {}

#endif

}
}

/*--------------------------------------------------------------------*/

namespace stk {
namespace {

//----------------------------------------------------------------------

class ParInBuf : public std::streambuf {
public:
  enum { BUFFER_LENGTH  = 0x010000 /* 64k bytes */ };
  enum { BUFFER_PUTBACK = 0x000010 /*  16 bytes */ };
  enum { MAX_READ       = BUFFER_LENGTH - BUFFER_PUTBACK };
  ParInBuf( ParallelMachine , const char * const );
  virtual ~ParInBuf();

protected:
   virtual int underflow(); // refill the input buffer
   virtual int overflow( int c = EOF ); // Not called
   virtual int sync(); // No-op
   virtual std::streambuf * setbuf( char * , std::streamsize ); // No-op

private:
  void close();

  ParallelMachine m_comm ;
  /// \todo REFACTOR: Investigate using C++ IO
  std::FILE     * m_root_fp ;
  char            m_buffer[ BUFFER_LENGTH ];
};

ParInBuf::ParInBuf( ParallelMachine comm , const char * const file_name )
  : m_comm( comm ), m_root_fp( NULL )
{
  std::memset(m_buffer, BUFFER_LENGTH, sizeof(char));
  int result = 1 ;

  if ( 0 == parallel_machine_rank( comm ) && NULL != file_name ) {
    result = NULL != ( m_root_fp = std::fopen( file_name , "r" ) );
  }

  broadcast( m_comm , & result , sizeof(int) );

  if ( ! result ) {
    std::string msg;
    msg.append("stk::ParallelInputStream( " );
    if ( 0 == parallel_machine_rank( comm ) && NULL != file_name ) {
      msg.append( file_name );
    }
    else {
      msg.append( "<NULL>" );
    }
    msg.append( " ) FAILED" );
    throw std::runtime_error(msg);
  }
}

void ParInBuf::close()
{
  if ( NULL != m_root_fp ) { std::fclose( m_root_fp ); m_root_fp = NULL ; }
  setg(NULL,NULL,NULL);
}

ParInBuf::~ParInBuf()
{
    try
    {
        close();
    }
    catch(...)
    {
        std::cerr << "Caught exception when trying to close file at " << __FILE__ ":" << __LINE__ << std::endl;
    }
}

int ParInBuf::underflow()
{
  char * const buf = m_buffer + BUFFER_PUTBACK ;
  int nread = 0 ;

  if ( gptr() == NULL || egptr() <= gptr() ) {
    if ( NULL != m_root_fp ) { nread = std::fread(buf,1,MAX_READ,m_root_fp); }
    broadcast( m_comm , & nread , sizeof(int) );
  }

  if ( 0 < nread ) {
    broadcast( m_comm , buf , nread );
    setg( m_buffer , buf , buf + nread );
  }
  else {
    close();
  }

  return 0 < nread ? *buf : EOF ;
}

namespace {

void throw_overflow()
{
  std::string msg ;
  msg.append("stk::ParallelInputStream::overflow CALL IS ERRONEOUS" );
  throw std::runtime_error(msg);
}

}

int ParInBuf::overflow( int )
{ throw_overflow(); return EOF ; }

int ParInBuf::sync()
{ return 0 ; }

std::streambuf * ParInBuf::setbuf( char * , std::streamsize )
{
  return this ;
}

//----------------------------------------------------------------------

} // namespace


ParallelInputStream::ParallelInputStream(
  ParallelMachine comm ,
  const char * const file_name )
  : std::istream( new ParInBuf( comm , file_name ) )
{}

ParallelInputStream::~ParallelInputStream()
{ delete rdbuf(); }

} // namespace stk


