#ifndef MYZERO_HPP
#define MYZERO_HPP
#include <string.h>

inline void myzero( int *data, size_t num )
{
  memset( data, 0, num*sizeof(int) );
}

inline void myzero(short *data, size_t num )
{
  memset( data, 0, num*sizeof(short) );
}

inline void myzero(unsigned short *data, size_t num )
{
  memset( data, 0, num*sizeof(unsigned short) );
}

inline void myzero( float *data, size_t num )
{
  memset( data, 0, num*sizeof(float) );
}

inline void myzero( double *data, size_t num )
{
  memset( data, 0, num*sizeof(double) );
}

inline void myzero( signed char *data, size_t num )
{
  memset( data, 0, num*sizeof(signed char) );
}

inline void myzero( unsigned char *data, size_t num )
{
  memset( data, 0, num*sizeof(unsigned char) );
}
#endif // MYZERO_HPP
