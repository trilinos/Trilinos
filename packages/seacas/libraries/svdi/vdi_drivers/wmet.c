/*
 * Copyright (C) 2009-2017 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *     * Neither the name of NTESS nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
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
 *
 */
/* $Id: wmet.c,v 1.23 2007/02/20 18:06:03 gdsjaar Exp $
 */

#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#if defined(CRAY)
#include <sys/types.h>
#endif
#include <stdlib.h>
#include <sys/file.h>
#include <unistd.h>

/*
 * SUN DEC/ULTRIX ALLIANT : C routines must have underscores
 * SGI CONVEX		  : C routines must have underscores
 *
 * VAX HP IBM/aix         : C routines do not have underscores
 *
 * CRAY/UNICOS            : C routines must be capitalized,
 *			      and no underscores
 */

#if defined(CRAY)
#define vdgnam WMETGN
#define cdrofs WMETOF
#define cdroff WMETFF
#define cdrcfs WMETCF
#define wmetbf WMETBF
#endif

#if !defined(CRAY) && !defined(ADDC_)
#define vdgnam wmetgn
#define cdrofs wmetof
#define cdroff wmetff
#define cdrcfs wmetcf
#endif

#if defined(ADDC_)
#define vdgnam wmetgn_
#define cdrofs wmetof_
#define cdroff wmetff_
#define cdrcfs wmetcf_
#define wmetbf wmetbf_
#endif

/* #define BUFFER_SIZE 132 */
#define BUFFER_SIZE 8192

extern char *getenv();

static char filename[100];       /* name of file */
static int  file_d = -1;         /* file descriptor - -1 if file not open */
static char buffer[BUFFER_SIZE]; /* for buffering data for output to file */
static int  buff_ptr = 0;        /* keep track of buffer space */

#if !defined(__convex__) && !defined(ultrix) && !defined(__osf__) && !defined(linux) &&            \
    !defined(__PARAGON__) && !defined(__hpux) && !defined(interix) && !defined(__CYGWIN__)
void vdgnam(name) char *name;
{
}
#endif

/*
Open file for sequential access.  FORTRAN unit number is ignored.
 */
void cdrofs(ifilcd) int *ifilcd; /* FORTRAN unit number ignored, provide for compatibility */
{
  char *fname;

  /* Determine name of file (if specified) in environment variable
     DUAL_FILENAME
     */
  if ((fname = getenv("DUAL_FILENAME")) == (char *)NULL)
    strcpy(filename, "dual_device.met");
  else
    strcpy(filename, fname);

  /* open the file  - if it doesn't exist, create it with mode 664 */
  /* -- open returns a file descriptor which is stored in the statelist */
  /* -- O_TRUNC ??? read/writing??? */
  if ((file_d = open(filename, (O_CREAT | O_TRUNC | O_RDWR), 0664)) == -1) {
    printf("cdrofs: error opening output file\n");
    exit(1);
  }
}

void cdroff(ifilcd, idum1, idum2, idum3) int *ifilcd, *idum1, *idum2, *idum3;
{
  cdrofs(ifilcd);
}

/*
  Close file sequential. FORTRAN unit number is ignored.
*/
void cdrcfs(ifilcd, eof) int *ifilcd, *eof;
{
  char buf[4];
  int  istat;

  /* if eof = 1 then write eof on file */
  if (*eof == 1) {
    *buf = EOF;
    if ((istat = write(file_d, buf, 4)) == -1) {
      printf("cdrcfs: error writing EOF\n");
      exit(1);
    }
  }
  /* close the file */
  close(file_d);
  file_d = -1;
}

/*
     Metafile buffering routine. (device dependent )
     Metafile driver passes OUTARY as an integer array, always with
     16 bits/word, right justified.  This buffering routine assumes
     16 bits/word and assumes 8 bit bytes.
*/
void     wmetbf(numwds, outary) int *numwds; /* number of words in OUTARY, 0 = buffer flush */
unsigned outary[];                           /* data to be buffered */
{
  static unsigned mask1 = ~(~0 << 8) << 8; /* mask off higher 8 bits */
  static unsigned mask2 = ~(~0 << 8);      /* mask off lower 8 bits */

  int i; /* loop variable */
  int istat;

  /* check for buffer flush and if there is something to flush */
  if (*numwds <= 0 && buff_ptr > 0) {

    /* write out the data as a byte stream. */
    if ((istat = write(file_d, buffer, buff_ptr)) == -1) {
      printf("wmetbf: write error\n");
    } /* end write buffer */

    /* reset buff pointer */
    buff_ptr = 0;

  } /* end if buffer flush */

  else /* pack outary into buffer */
  {
    if (buff_ptr + *numwds >= BUFFER_SIZE) {

      /* write out the data as a byte stream. */
      if ((istat = write(file_d, buffer, buff_ptr)) == -1) {
        printf("wmetbf: write error\n");
      } /* end write buffer */

      /* reset buff pointer */
      buff_ptr = 0;
    }
    for (i = 0; i < *numwds; i++) {
      buffer[buff_ptr++] = (outary[i] & mask1) >> 8;
      buffer[buff_ptr++] = (outary[i] & mask2);
    } /* end for i=... */
  }   /* end else */
} /* end wmetbf */
