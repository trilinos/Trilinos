/* ------------------------------------------------------------------ */
/* Copyright 2000 Sandia Corporation, Albuquerque, NM.                */
/* ------------------------------------------------------------------ */

#ifndef STK_UTIL_UTIL_null_streambuf_hpp
#define STK_UTIL_UTIL_null_streambuf_hpp

#include <iostream>
#include <cstdio>  /* Defines EOF */

//: Specialize the ANSI Standard C++ streambuf class
//: that throws away everything given to it without
//: generating an error.

class null_streambuf : public std::streambuf {
public:

  //: Constructor
  null_streambuf();

  //: Destructor
  virtual ~null_streambuf();

protected:

  //: Called when output buffer is filled
  virtual int overflow( int c = EOF );

  //: Sync is a no-op
  virtual int sync();

  //: Setbuf is a no-op
  virtual std::streambuf * setbuf( char * s , std::streamsize n );

private:

  null_streambuf( const null_streambuf & ); // Not allowed
  null_streambuf & operator = ( const null_streambuf & ); // Not allowed

  char buf[64]; // Throw away buffer
};

/*--------------------------------------------------------------------*/

#endif // STK_UTIL_UTIL_null_streambuf_hpp
