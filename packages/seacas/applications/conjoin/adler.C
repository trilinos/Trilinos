// *********************************************************************
// Description: The following source code is borrowed heavily from the
// zlib project The zlib source is subject to the following copyright
// and license:
//
//   Copyright (C) 1995-2010 Jean-loup Gailly and Mark Adler
//
//   This software is provided 'as-is', without any express or implied
//   warranty.  In no event will the authors be held liable for any damages
//   arising from the use of this software.
//
//   Permission is granted to anyone to use this software for any purpose,
//   including commercial applications, and to alter it and redistribute it
//   freely, subject to the following restrictions:
//
//   1. The origin of this software must not be misrepresented; you must not
//      claim that you wrote the original software. If you use this software
//      in a product, an acknowledgment in the product documentation would be
//      appreciated but is not required.
//   2. Altered source versions must be plainly marked as such, and must not be
//      misrepresented as being the original software.
//   3. This notice may not be removed or altered from any source distribution.
//
//   Jean-loup Gailly        Mark Adler
//   jloup@gzip.org          madler@alumni.caltech.edu

#include "adler.h"

#define DO1(buf,i)  {s1 += buf[i]; s2 += s1;}
#define DO2(buf,i)  DO1(buf,i); DO1(buf,i+1);
#define DO4(buf,i)  DO2(buf,i); DO2(buf,i+2);
#define DO8(buf,i)  DO4(buf,i); DO4(buf,i+4);
#define DO16(buf)   DO8(buf,0); DO8(buf,8);


size_t adler(size_t adler, const void* vbuf, size_t len)
{
  const size_t BASE = 65521; /* largest prime smaller than 65536 */
  const size_t NMAX = 5552;  /* NMAX is the largest n such
				that 255n(n+1)/2 + (n+1)(BASE-1) <= 2^32-1 */

  const unsigned char *buf = (const unsigned char*)vbuf;

  size_t s1 = adler & 0xffff;
  size_t s2 = (adler >> 16) & 0xffff;
  int k = 0;

  if (buf == NULL || len == 0) return 1L;

  while (len > 0) {
    k = len < NMAX ? len : NMAX;
    len -= k;
    while (k >= 16) {
      DO16(buf);
      buf += 16;
      k -= 16;
    }
    if (k != 0) do {
      s1 += *buf++;
      s2 += s1;
    } while (--k);
    s1 %= BASE;
    s2 %= BASE;
  }
  return (s2 << 16) | s1;
}
