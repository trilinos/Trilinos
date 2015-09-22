
/*
//@HEADER
// ***********************************************************************
// 
//        AztecOO: An Object-Oriented Aztec Linear Solver Package 
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER
*/

#include <stdio.h>
#include <stdarg.h>

#include <AztecOO_ConfigDefs.h>
#include <az_aztec.h>

/** az_ostream_out is a singleton class that holds a std::ostream pointer
  that can be used as the destination for Aztec's normal printf output.
  az_ostream_out is only used by the functions AZOO_printf_out and
  AZOO_flush_out below.
*/
class az_ostream_out {
 public:
  virtual ~az_ostream_out(){}

  std::ostream* get_ostream()
  { return strm_; }

  static az_ostream_out& get_instance()
  {
    static az_ostream_out az_ostrm_out_;
    return( az_ostrm_out_ );
  }

  void set_ostream(std::ostream& ostrm)
  {
    strm_ = &ostrm;
  }

  void clear_ostream()
  {
    strm_ = 0;
  }

 private:
  az_ostream_out() : strm_(0) {}

  std::ostream* strm_;
};

/** az_ostream_err is a singleton class that holds a std::ostream pointer
  that can be used as the destination for Aztec's normal stderr output.
  az_ostream_out is only used by the function AZOO_printf_err below.
  AZOO_flush_out below.
*/
class az_ostream_err {
 public:
  virtual ~az_ostream_err(){}

  std::ostream* get_ostream()
  { return strm_; }

  static az_ostream_err& get_instance()
  {
    static az_ostream_err az_ostrm_err_;
    return( az_ostrm_err_ );
  }

  void set_ostream(std::ostream& ostrm)
  {
    strm_ = &ostrm;
  }

  void clear_ostream()
  {
    strm_ = 0;
  }

 private:
  az_ostream_err() : strm_(0) {}

  std::ostream* strm_;
};


/** AZOO_printf_out gets called (through a call-back) when Aztec's C code
   calls the C function AZ_printf_out(), but only if the function
   AZOO_set_stream_out() below has been called. The purpose of this
   (AZOO_printf_out()) function is to send Aztec's printf output to a
   user-specified std::ostream.
*/
int AZOO_printf_out(const char* format, va_list aptr)
{
  az_ostream_out& azout = az_ostream_out::get_instance();
  std::ostream* ostrm = azout.get_ostream();

  if (ostrm == 0) {
    return( vprintf(format, aptr) );
  }

  static char buf[512];
  for(unsigned i=0; i<512; ++i) buf[i] = '\0';

  int ret = vsprintf(buf, format, aptr);

  *ostrm << buf;

  return( ret );
}

void AZOO_flush_out()
{
  az_ostream_out& azout = az_ostream_out::get_instance();
  std::ostream* ostrm = azout.get_ostream();

  if (ostrm == 0) {
    fflush(stdout);
  }
  else {
    *ostrm << std::flush;
  }
}

int AZOO_printf_err(const char* format, va_list aptr)
{
  az_ostream_err& azerr = az_ostream_err::get_instance();
  std::ostream* ostrm = azerr.get_ostream();

  if (ostrm == 0) {
    return( vfprintf(stderr, format, aptr) );
  }

  static char buf[512];
  for(unsigned i=0; i<512; ++i) buf[i] = '\0';

  int ret = vsprintf(buf, format, aptr);

  *ostrm << buf;

  return( ret );
}

/** The following function-pointers are declared/owned in the file
  az_printf.c
*/
extern "C" {
  extern int (*azoo_printf_out)(const char* format, va_list aptr);
  extern void (*azoo_flush_out)();
  extern int (*azoo_printf_err)(const char* format, va_list aptr);
}

/** AZOO_set_stream_out sets the specified ostrm as the target for the
  AZOO_printf_out function, and then sets the above function-pointer
  'azoo_printf_out' so that Aztec's printf output will be given to
  the std::ostream rather than to stdout. (azoo_flush_out is also set
  so that Aztec can flush the output stream.)
*/
void AZOO_set_stream_out(std::ostream& ostrm)
{
  az_ostream_out& azout = az_ostream_out::get_instance();
  azout.set_ostream(ostrm);

  azoo_printf_out = AZOO_printf_out;
  azoo_flush_out  = AZOO_flush_out;
}

/** AZOO_set_stream_err sets the specified ostrm as the target for the
  AZOO_printf_err function, and then sets the above function-pointer
  'azoo_printf_err' so that Aztec's stderr output will be given to
  the std::ostream rather than to stderr.
*/
void AZOO_set_stream_err(std::ostream& ostrm)
{
  az_ostream_err& azerr = az_ostream_err::get_instance();
  azerr.set_ostream(ostrm);

  azoo_printf_err = AZOO_printf_err;
}

/** AZOO_clear_streams sets the out and err streams to NULL
*/
void AZOO_clear_streams()
{
  az_ostream_out& azout = az_ostream_out::get_instance();
  azout.clear_ostream();

  az_ostream_err& azerr = az_ostream_err::get_instance();
  azerr.clear_ostream();
}

