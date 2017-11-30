/*
* Copyright(C) 2005 National Technology & Engineering Solutions of
* Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
* NTESS, the U.S. Government retains certain rights in this software.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are
* met:
*
* *  Redistributions of source code must retain the above copyright
*    notice, this list of conditions and the following disclaimer.
*
* *  Redistributions in binary form must reproduce the above
*   copyright notice, this list of conditions and the following
*   disclaimer in the documentation and/or other materials provided
*   with the distribution.
*
* *  Neither the name of NTESS nor the names of its
*   contributors may be used to endorse or promote products derived
*  from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
* "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
* A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
* OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
* LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
* DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
* THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
* OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int count_digits(int value)
{
  char tmp[32];
  sprintf(tmp, "%d", value);
  return strlen(tmp);
}

/* Output a zero-padded sequence of digits from [0..limit).  For
   example zpad 10 will output 00, 01, 02, 03, ..., 09 each on a
   single line.
 */
int main(int argc, char **argv)
{
  int start = 0;
  if (argc < 2) {
    fprintf(stderr, "Usage: %s limit [start=0]\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  if (argc == 3) {
    start = strtol(argv[2], NULL, 10);
  }

  {
    int  i;
    char format[] = "%.0Xd\n";
    char digits[2];
    int  limit = strtol(argv[1], NULL, 10);

    /* Count number of digits needed to represent 'limit' */
    int width = count_digits(limit);
    sprintf(digits, "%d", width);

    /* Create an output format that will zero-pad to that width */
    format[3] = digits[0];
    for (i = start; i < limit; i++) {
      printf(format, i);
    }
  }
  return EXIT_SUCCESS;
}
