/*
 * Copyright(C) 1999-2023 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
/* ifdef.h - ifdef file for cdr routines
 * This file is used to define the system dependent ways C is
 * called from FORTRAN.  Underscores are used by default.
 *
 * SUN DEC/ULTRIX ALLIANT : C routines must have underscores
 * SGI CONVEX             : C routines must have underscores
 *
 * VAX HP IBM/aix         : C routines do not have underscores
 *
 * CRAY/UNICOS            : C routines must be capitalized,
 *                            and no underscores
 *
 * This include file is used by all cdr routines
 */

#include <stdlib.h>
#include <unistd.h>

#if defined(ADDC_)
#endif
#if !defined(CRAY) && !defined(ADDC_) && !defined(COUGAR)
#define cdrcom_ cdrcom
#define cdrcm2_ cdrcm2
#define cdrunx_ cdrunx
#define vcjob_  vcjob
#define vconod_ vconod
#define cdr1ch_ cdr1ch
#define cdr1pk_ cdr1pk
#define cdr2cm_ cdr2cm
#define cdra2c_ cdra2c
#define cdrbit_ cdrbit
#define cdrbtr_ cdrbtr
#define cdrc2a_ cdrc2a
#define cdrc2h_ cdrc2h
#define cdrcfs_ cdrcfs
#define cdrcvt_ cdrcvt
#define cdrela_ cdrela
#define cdrgnm_ cdrgnm
#define cdri2c_ cdri2c
#define cdrinp_ cdrinp
#define cdrjob_ cdrjob
#define cdrlwr_ cdrlwr
#define cdrmon_ cdrmon
#define cdroab_ cdroab
#define cdrof3_ cdrof3
#define cdroff_ cdroff
#define cdrofs_ cdrofs
#define cdrono_ cdrono
#define cdrout_ cdrout
#define cdrpak_ cdrpak
#define cdrrfs_ cdrrfs
#define cdrrvt_ cdrrvt
#define cdrstd_ cdrstd
#define cdrtbk_ cdrtbk
#define cdrtim_ cdrtim
#define cdrtod_ cdrtod
#define cdrupk_ cdrupk
#define cdrupr_ cdrupr
#define cdrwfs_ cdrwfs
#endif
#if defined(CRA)
#define cdrcom_ CDRCOM
#define cdrcm2_ CDRCM2
#define cdrunx_ CDRUNX
#define vcjob_  VCJOB
#define vconod_ VCONOD
#define cdr1ch_ CDR1CH
#define cdr1pk_ CDR1PK
#define cdr2cm_ CDR2CM
#define cdra2c_ CDRA2C
#define cdrbit_ CDRBIT
#define cdrbtr_ CDRBTR
#define cdrc2a_ CDRC2A
#define cdrc2h_ CDRC2H
#define cdrcfs_ CDRCFS
#define cdrcvt_ CDRCVT
#define cdrela_ CDRELA
#define cdrgnm_ CDRGNM
#define cdri2c_ CDRI2C
#define cdrinp_ CDRINP
#define cdrjob_ CDRJOB
#define cdrlwr_ CDRLWR
#define cdrmon_ CDRMON
#define cdroab_ CDROAB
#define cdrof3_ CDROF3
#define cdroff_ CDROFF
#define cdrofs_ CDROFS
#define cdrono_ CDRONO
#define cdrout_ CDROUT
#define cdrpak_ CDRPAK
#define cdrrfs_ CDRRFS
#define cdrrvt_ CDRRVT
#define cdrstd_ CDRSTD
#define cdrtbk_ CDRTBK
#define cdrtim_ CDRTIM
#define cdrtod_ CDRTOD
#define cdrupk_ CDRUPK
#define cdrupr_ CDRUPR
#define cdrwfs_ CDRWFS
#endif
/* end ifdef.h */

/* Structures used to communicate with FORTRAN common blocks.
 * On all systems, except CRAY/UNICOS, structures must be declared
 * external.  On UNICOS, the structures are not external.
 */
#if defined(CRA)

/* cdrcm2_u.h */
struct cdr2
{
  char KGNAME[80];
} cdrcm2_;

/* cdrcom.h */
struct cdr
{
  int KWRTFL;
  int KRDFL;
  int KOUTFL;
  int KINFL;
  int KWRDSZ;
  int KBYTEL;
  int KCPW;
  int KBAUD;
  int KCOMTP;
} cdrcom_;

/* cdrunx.h */
struct unx
{
  int KUNTFD[1000];
} cdrunx_;

/* vcjob.h */
struct job
{
  int KIDSIZ;
  int KJOBID[4];
  int KUSRSZ;
  int KUSRID[4];
  int KSZROU;
  int KJROUT[4];
  int KSECUR;
  int KJTIME[4];
  int KJDATE[4];
  int MACHIN[3];
  int MACLEN;
} vcjob_;

/* for now, hard code this stuff in: */
#if defined(SECURITY_CODE)
#undef SECURITY_CODE
#endif
#define SECURITY_CODE 0

#if defined(CDR_MACHINE)
#undef CDR_MACHINE
#endif
#define CDR_MACHINE "UNICOS"

#if defined(ROUTE)
#undef ROUTE
#endif
#define ROUTE "BX 603"

/* vconod.h */
struct onode
{
  int ONODE;
} vconod_;

#else

/* cdrcm2.h  */
extern struct cdr2
{
  char KGNAME[80];
} cdrcm2_;

/* cdrcom.h */
extern struct cdr
{
  int KWRTFL;
  int KRDFL;
  int KOUTFL;
  int KINFL;
  int KWRDSZ;
  int KBYTEL;
  int KCPW;
  int KBAUD;
  int KCOMTP;
} cdrcom_;

/* cdrunx.h */
extern struct unx
{
  int KUNTFD[1000];
} cdrunx_;

/* vcjob.h */
extern struct job
{
  int KIDSIZ;
  int KJOBID[4];
  int KUSRSZ;
  int KUSRID[4];
  int KSZROU;
  int KJROUT[4];
  int KSECUR;
  int KJTIME[3];
  int KJDATE[3];
  int MACHIN[3];
  int MACLEN;
} vcjob_;

/* for now, hard code this stuff in: */
#if defined(SECURITY_CODE)
#undef SECURITY_CODE
#endif
#define SECURITY_CODE 0

#if defined(CDR_MACHINE)
#undef CDR_MACHINE
#endif
#define CDR_MACHINE "UNIX"

#if defined(ROUTE)
#undef ROUTE
#endif
#define ROUTE "BX 603"

/* vconod.h */
extern struct onode
{
  int ONODE;
} vconod_;

#endif

/* System includes */
#include <ctype.h>
#if !defined(__CYGWIN__)
#include <pwd.h>
#endif
#include <signal.h>
#include <stdio.h>
#include <string.h>

#if defined(CRA)
#include <fcntl.h>
#include <sys/types.h>
#endif
#include <sys/file.h>

#if defined(CRA) || defined(interix) || defined(linux)
#include <time.h>
#else
#include <sys/time.h>
#endif
#if defined(sun) && (defined(SYSV) || defined(SVR4))
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#endif

/*  *** CDR1CH ***

  Copy the Nth KBYTEL-width character of STRING to the
  lowest order character field of NUCHAR.

  Parameters:
    n - IN - character position, counted from the left, of
             the character to be extracted from STRING
    string - IN - word from which to copy the character
    nuchar - OUT - word in which to place the character

*/

void cdr1ch_(int *n, unsigned *string, unsigned *nuchar)
{
  if (*n > cdrcom_.KCPW) {
    return;
  }
  *nuchar = (*string >> (cdrcom_.KWRDSZ - *n * cdrcom_.KBYTEL) & ~(~0u << cdrcom_.KBYTEL));
}

/*  *** CDR1PK ***

  Copy the KBYTEL-width character of INCHAR (right-justified) to
  the Nth KBYTEL-width character position of BUFFER.

  Parameters:
    n - IN - character position, counted from the left, at which
             the KBYTEL-width character is copied
    inchar - IN - word containing one KBYTEL-width character, right
                  adjusted, zero-filled
    buffer - OUT - word into which the character is copied.

*/

void cdr1pk_(int *n, unsigned *inchar, unsigned *buffer)
{
  int      i;
  unsigned temp, mask;

  if (*n > cdrcom_.KCPW) {
    return;
  }
  temp    = *inchar << (i = cdrcom_.KWRDSZ - (*n * cdrcom_.KBYTEL));
  mask    = ~(~0u << cdrcom_.KBYTEL) << i;
  *buffer = ((mask & temp) | (~mask & *buffer));
}
/*  *** CDR2CM ***

  Bitwise complement, used for 2's complement.

  Parameters:
    in - IN - integer to complement
    iout - OUT - complemented integer

*/

void cdr2cm_(int *in, int *iout) { *iout = ~(*in) + 1; }

/*  *** CDRPAK ***

  Pack the IWIDTH lower order bits from NEXT into IBUF, starting at
  bit position IBUF(IWORD) + IBITLC.  Bits are counted left to right.
  IWORD and IBITLC are updated to reflect the new position in IBUF.

  Parameters:
    next - IN - word containing bits to pack into ibuf
    iwidth - IN - number of bits to get from NEXT
    iword  - IN/OUT -word in buffer to put the bits
    ibitlc - IN/OUT -number of bits currently used in IBUF(IWORD)
    ibuf - OUT - the array in which to put the bits

*/

void cdrpak_(unsigned *next, int *iwidth, int *iword, int *ibitlc, unsigned *ibuf)
{
  int      idx, i;
  unsigned temp, mask;

  /* check boundaries */
  if (*iwidth > cdrcom_.KWRDSZ) {
    return;
  }

  /* if bitlc is bigger than the word size, increment */
  if (*ibitlc >= cdrcom_.KWRDSZ) {
    *ibitlc = 0;
    (*iword)++;
  }

  /* "C" start arrays at 0 */
  idx = *iword - 1;

  if (*ibitlc + *iwidth > cdrcom_.KWRDSZ) {

    /* word overflow - break into parts */

    temp      = *next >> (i = *iwidth + *ibitlc - cdrcom_.KWRDSZ);
    mask      = ~0u << (cdrcom_.KWRDSZ - *ibitlc);
    ibuf[idx] = ((~mask & temp) | (mask & ibuf[idx]));
    *ibitlc   = i;
    (*iword)++;
    idx++;
    temp      = *next << (i = cdrcom_.KWRDSZ - *ibitlc);
    mask      = ~0u << i;
    ibuf[idx] = ((mask & temp) | (~mask & ibuf[idx]));
  }
  else {

    /* it fits all in one word */

    temp      = *next << (i = cdrcom_.KWRDSZ - *ibitlc - *iwidth);
    mask      = ~(~0u << *iwidth) << i;
    ibuf[idx] = ((mask & temp) | (~mask & ibuf[idx]));
    *ibitlc += *iwidth;
  }
}

/*  *** CDRA2C ***

  Convert from ascii to character

  Parameters:
    asci - IN - integer representing an ascii character
    charac - OUT - the character represented by asci

*/

void cdra2c_(int *asci, char *charac) { charac[0] = (char)*asci; }
/*  *** CDRBIT ***

  This is a dummy routine to satisfy externals. See CDRPAK

  Parameters:
    next - IN -
    iwidth - IN -
    iword - IN/OUT -
    ibitlc - IN/OUT -
    ibuf - OUT -

*/

void cdrbit_(unsigned *next, int *iwidth, int *iword, int *ibitlc, unsigned *ibuf)
{
  cdrpak_(next, iwidth, iword, ibitlc, ibuf);
}

/*  *** CDRUPK ***

  Copy IWIDTH bits from IBUF,starting at bit position IBUF(IWORD) +
  IBITLC, into low order of NEXT.  Bits are counted left to right.
  IWORD and IBITLC are updated to reflect the new position in IBUF.

  Parameters:
    ibuf - IN - the array containing the bits to copy into next
    iword  - IN/OUT - word in buffer from which to get the bits
    ibitlc - IN/OUT - number of bits already retrieved from IBUF(IWORD)
    iwidth - IN - number of bits to get from IBUF
    next - OUT - word to put the bits in, right justified, zero filled

*/

void cdrupk_(unsigned *ibuf, int *iword, int *ibitlc, int *iwidth, unsigned *next)
{
  int      idx, i;
  unsigned temp, mask;
#if defined(DEC) || defined(linux) || defined(interix)
  unsigned char *p_temp, ctemp;
#endif
  /* check boundaries */
  if (*iwidth > cdrcom_.KWRDSZ) {
    return;
  }

  /* if ibitlc is bigger than word size, increment */
  if (*ibitlc >= cdrcom_.KWRDSZ) {
    *ibitlc = 0;
    (*iword)++;
  }

  /* "C" arrays start at 0 */
  idx = *iword - 1;

  if (*ibitlc + *iwidth > cdrcom_.KWRDSZ) {

    /* crosses word boundary - break it up */

    temp = ibuf[idx];
#if defined(DEC) || defined(linux) || defined(interix)
    /* swap byte order around */
    p_temp        = (unsigned char *)(&temp);
    ctemp         = *p_temp;
    *p_temp       = *(p_temp + 3);
    *(p_temp + 3) = ctemp;
    ctemp         = *(p_temp + 1);
    *(p_temp + 1) = *(p_temp + 2);
    *(p_temp + 2) = ctemp;
#endif
    temp    = temp << (i = *iwidth + *ibitlc - cdrcom_.KWRDSZ);
    mask    = ~(~0u << (cdrcom_.KWRDSZ - *ibitlc)) << i;
    *next   = (mask & temp);
    *ibitlc = i;
    (*iword)++;
    idx++;
    temp = ibuf[idx];
#if defined(DEC) || defined(linux) || defined(interix)
    /* swap byte order around */
    p_temp        = (unsigned char *)(&temp);
    ctemp         = *p_temp;
    *p_temp       = *(p_temp + 3);
    *(p_temp + 3) = ctemp;
    ctemp         = *(p_temp + 1);
    *(p_temp + 1) = *(p_temp + 2);
    *(p_temp + 2) = ctemp;
#endif
    temp  = temp >> (cdrcom_.KWRDSZ - *ibitlc);
    mask  = ~(~0u << *ibitlc);
    *next = ((mask & temp) | *next);
  }
  else {

    /* doesn't cross word boundary */

    temp = ibuf[idx];
#if defined(DEC) || defined(linux) || defined(interix)
    /* swap byte order around */
    p_temp        = (unsigned char *)(&temp);
    ctemp         = *p_temp;
    *p_temp       = *(p_temp + 3);
    *(p_temp + 3) = ctemp;
    ctemp         = *(p_temp + 1);
    *(p_temp + 1) = *(p_temp + 2);
    *(p_temp + 2) = ctemp;
#endif
    temp  = temp >> (cdrcom_.KWRDSZ - *ibitlc - *iwidth);
    mask  = ~(~0u << *iwidth);
    *next = mask & temp;
    *ibitlc += *iwidth;
  }
}

/*  *** CDRBTR ***

  This is a dummy routine to satisfy externals. See CDRUPK

  Parameters:
    ibuf - IN -
    iword - IN/OUT -
    ibitlc - IN/OUT -
    iwidth - IN -
    next - OUT -

*/

void cdrbtr_(unsigned ibuf[], int *iword, int *ibitlc, int *iwidth, unsigned *next)
{
  cdrupk_(ibuf, iword, ibitlc, iwidth, next);
}
/*  *** CDRC2A ***

  Convert from character to ASCII

  Parameters:
    charac - IN - the character represented by asci
    asci - OUT - integer representing an ASCII character

*/

void cdrc2a_(char *charac, int *asci) { *asci = (int)charac[0]; }

/*  *** CDRC2H ***

  Convert character to hollerith. TITLE will contain at least
  LENGTH characters, and ITITLE will be dimensioned one larger
  in order to append the zero word that is required by Hollerith
  rules.

  Notes and Revisions:
    This routine just copies the input to the output and appends a
    zero word, the C/FORTRAN parameter passing handles the rest.

  Parameters:
    title - IN - character string of any length
    length - IN - length of title
    ititle - OUT - hollerith equivalent of title

*/

void cdrc2h_(unsigned *title, int *length, unsigned *ititle)
{
  int i, ilen;

  ilen = *length / cdrcom_.KCPW;
  for (i = 0; i < ilen; i++) {
    ititle[i] = title[i];
  }
  ititle[ilen + 1] = 0;
}

/*  *** CDRCFS ***

  Close file sequential.

  Notes and Revisions:
    This routine uses the array KUNTFD, in the FORTRAN common and
    the C external structure, cdrunx, to translate from FORTRAN
    unit number to UNIX/C file descriptors.

  Parameters:
    ifilcd - IN - FORTRAN unit number of the sequential file to close
    eof - IN - end-of-file request flag. If eof=1, an end of file
               marker is written on the end of the file.

*/

void cdrcfs_(int *ifilcd, int *eof)
{
  /* translate fortran unit number to file descriptor */
  int fd = cdrunx_.KUNTFD[*ifilcd];

  /* if eof = 1 then write eof on file */
  if (*eof == 1) {
    char buf[4];
    *buf = EOF;
    if (write(fd, buf, 4) == -1) {
      perror("CDRCFS error:");
    }
  }

  /* close the file */
  close(fd);
}
/*  *** CDRCVT ***

  Convert from internal character set to ASCII character set.

  Notes and Revisions:
    This routine is here to satisfy the call, since the internal
    character set is ASCII.

  Parameters:
    in - IN - character to be converted, in host internal format
    iout - OUT - character converted to ASCII

*/

void cdrcvt_(int *in, int *iout) { *iout = *in; }
/*  *** CDRELA ***

  Collect CPU time.  This is a dummy routine to satisfy calls from
  VDMONI. Monitoring isn't done on this system

  Parameters:
    icode - IN - 0 = initialize timer, 1 = save timer

*/
void cdrela_(int *icode) {}

/*  *** CDRGNM ***

  Name graphics output file or the post processor input file

  Notes and Revisions:
    This routine is called by VDGNAM (if the user called VDGNAM),
    and VBPKG.  Whoever gets here first, wins.

  Parameters:
    gname - IN - file name to call graphics output file
    glen - IN - length of file name
*/

void cdrgnm_(char gname[], int *glen)
{
  /* make sure name hasn't already been set */
  char blank = ' ';
  if (cdrcm2_.KGNAME[0] != blank) {
    return;
  }

  /* store it in the external (and common) variables */
  for (int i = 0; i < *glen; i++) {
    cdrcm2_.KGNAME[i] = gname[i];
  }

  /* append a backslash and an endofstring character */
  char bslash               = '\\';
  cdrcm2_.KGNAME[*glen]     = bslash;
  cdrcm2_.KGNAME[*glen + 1] = '\0';
}

/*  *** CDRI2C ***

  Convert positive integer to decimal character string equivalent

  Parameters:
    intp - IN - positive integer to be converted
    ndigit - IN - number of digits to be produced in string form
                 (pad left with zeros)
    string - OUT - character string of at least ndigit characters

*/

void cdri2c_(int *intp, int *ndigit, char string[])
{
  int inbr = *intp;

  for (int i = (*ndigit - 1); i >= 0; i--) {
    string[i] = inbr % 10 + '0';
    inbr      = inbr / 10;
  }
}
#define CDR_MAXLEN 512

/*  *** CDRINP ***

  Send the prompt, IPROMPT, to the terminal, and read ICOUNT
  characters from the terminal into buffer IBUFFER.

  Parameters:
    icount - IN - number of characters to be accepted from the
                  terminal
    ibuffer - OUT - buffer into which the characters are read
    iprompt - iN - prompt sequence to output to the terminal

*/

void cdrinp_(int *icount, int *ibuffer, int *iprompt)
{
  char buffer[CDR_MAXLEN];
  int  istat, actcnt = 0;

  /* check the input */
  int icnt = *icount;
  if (*icount > CDR_MAXLEN) {
    icnt = CDR_MAXLEN;
  }

  /* if prompt count is greater than 0, output the prompt */
  if (iprompt[0] > 0) {
    char prompt[CDR_MAXLEN];
    for (int i = 1; i <= iprompt[0]; i++) {
      prompt[i - 1] = (char)iprompt[i];
    }

    if ((istat = write(1, prompt, iprompt[0])) == -1) {
      perror("CDRINP error:");
    }
  }

  /* read from the terminal */
  if ((istat = read(1, buffer, icnt)) == -1) {
    perror("CDRINP error:");
  }
  else {
    /* strip off newline if there is one */
    if (buffer[istat - 1] == '\n') {
      actcnt = istat - 1;
    }
    else {
      actcnt = istat;
    }
  }

  /* load user's buffer  */
  for (int i = 0; i < actcnt; i++) {
    ibuffer[i] = (int)buffer[i];
  }

  /* zero fill left over buffer space */
  if (actcnt < icnt) {
    for (int i = actcnt; i < icnt; i++) {
      ibuffer[i] = 0;
    }

    /* if there wasn't a newline, flush input */
    if (buffer[istat - 1] != '\n') {
      read(1, buffer, CDR_MAXLEN);
    }
  }
}

/*  *** CDRJOB ***

Get job information and set the vars in the vcjob common block:
process id -  max 12 chars
length of process id
user name  - max 12 chars
length of user name
routing info - max 12 chars
length of routing
security classification
date  86-12-15
time  12:23:59
machine - max 12 chars
length of machine
*/

void cdrjob_(void) {}

/*  *** CDRLWR ***

Convert character string to all lower case

Notes and Revisions:
Only converts characters

Parameters:
in - IN - character string to be converted
iout - OUT - character string converted to all lower case
inlen - IN - length of character string in

*/

void cdrlwr_(char in[], char iout[], int *inlen)
{
  for (int i = 0; i <= *inlen - 1; i++) {
    if (isupper(in[i])) {
      iout[i] = tolower(in[i]);
    }
    else {
      iout[i] = in[i];
    }
  }
}
/*  *** CDRMON ***

Write out monitoring information.  This is a dummy routine to
satisfy calls from VDMONI. Monitoring isn't done on this system.

Parameters:
mdev - IN -
mpkg - IN -
mpage - IN -

*/

void cdrmon_(int *mdev, int *mpkg, int *mpages) {}

/*  *** CDROFS ***

Open a file for sequential access.

Notes and Revisions:
This routine uses the array KUNTFD, in the FORTRAN common and
the C external structure, cdrunx, to translate from FORTRAN
unit number to UNIX/C file descriptors.

Parameters:
ifilcd - IN - the FORTRAN unit number of the file to open

*/

void cdrofs_(int *ifilcd)
{
  /* check for graphics output file. file name is stored in KGNAME
     from external (common) cdrcm2. if file hasn't been named, or
     unit isn't KOUTFL, default to file{unit}          */

  char symbol[1024];
  char blank = ' ';
  if (*ifilcd == cdrcom_.KOUTFL && cdrcm2_.KGNAME[0] != blank) {
    int  i      = 0;
    int  j      = 0;
    char bslash = '\\';
    while (cdrcm2_.KGNAME[i] != bslash) {
      symbol[i++] = cdrcm2_.KGNAME[j++];
    }
    symbol[i] = '\0';
  }
  else {
    snprintf(symbol, 1024, "file%d", *ifilcd);
  }

  /* check the environment to see if a file name has been assigned */
  char *env = getenv(symbol);
  if (env != NULL && strlen(env) < 1024) {
    snprintf(symbol, 1024, "%s", env);
  }

  /* open the file  - if it doesn't exist, create it with mode 664 */
  int fd;
  if ((fd = open(symbol, (O_CREAT | O_RDWR), 0664)) == -1) {
    int  errnum = 722;
    int  errsev = 10;
    char err[50];
    snprintf(err, 50, "SVDI ERROR NUMBER %d SEVERITY CODE %d", errnum, errsev);
    perror(err);
  }
  else {

    /* associate fortran unit number with file descriptor */
    cdrunx_.KUNTFD[*ifilcd] = fd;
  }
}

/*  *** CDROAB ***

Open file for the Abekas

*/

void cdroab_(int *ifilcd, int *frame)
{
  char ic[5];
  int  len;

  /* can't call cdrgnm, because it only lets you name
     the file once, so set the name in the common block */

  len = 4;
  cdri2c_(frame, &len, ic);
  ic[len] = '\0';
  snprintf(cdrcm2_.KGNAME, 80, "%s.RGB", ic);

  cdrofs_(ifilcd);
}
/*  *** CDROF3 ***

Open a file for sequential access, fixed length, formatted

Notes and Revisions:
C file structures don't have the concept of record structure.
This routine calls CDROFS

Parameters:
ifilcd - IN - the FORTRAN unit number of the file to open
idum - IN - originally used for record length

*/

void cdrof3_(int *ifilcd, int *idum) { cdrofs_(ifilcd); }
/*  *** CDROFF ***

Open a file for sequential access, fixed length, formatted

Notes and Revisions:
C file structures don't have the concept of record structure.
This routine calls CDROFS

Parameters:
ifilcd - IN - the FORTRAN unit number of the file to open

*/

void cdroff_(int *ifilcd) { cdrofs_(ifilcd); }

/*  *** CDRONO ***

Set flag to designate this file is going to the output node

Parameters:
setono - IN - logical value

*/

void cdrono_(int *setono) { vconod_.ONODE = *setono; }
#define MAXBUFFER 512

/*  *** CDROUT ***

Transmit ICOUNT words from BUFFER to the terminal

Parameters:
ICOUNT - IN - number of characters to be sent to the terminal
BUFFER - IN - integer array containing the characters

*/

void cdrout_(int *icount, int *buffer)
{
  int count = *icount;
  int j     = 0;
  while (count > 0) {
    int c_count = (count > MAXBUFFER) ? MAXBUFFER : count;
    count       = count - c_count;

    char c_buffer[MAXBUFFER];
    for (int i = 0; i < c_count; i++) {
      c_buffer[i] = (char)buffer[j++];
    }

    /* write to standard output */
    if (write(1, c_buffer, c_count) == -1) {
      perror("CDROUT error");
    }
  }
}

/*  *** CDRRFS ***

Try to read LENGTH words from the file.  Return the actual number
of words read.

Notes and Revisions:
This routine uses the array KUNTFD, in the FORTRAN common and
the C external structure, cdrunx, to translate from FORTRAN
unit number to UNIX/C file descriptors.
read() will not read beyond EOF. read() returns the actual
number of bytes read, or -1 for fail.

Parameters:
ifilcd - IN - FORTRAN unit number of file to read
length - IN - number of words to try to read
buffer - OUT - array of size LENGTH into which the data is read
istat - OUT - status indicator: 0 = end of file;
>0 = actual number of words read;
<0 = error

*/

void cdrrfs_(int *ifilcd, int *length, char *buffer, int *istat)
{
  /* translate fortran unit number to file descriptor */
  int fd = cdrunx_.KUNTFD[*ifilcd];

  /* if the file hasn't been opened, open it */
  if (fd == 0) {
    cdrofs_(ifilcd);
    fd = cdrunx_.KUNTFD[*ifilcd];
  }

  /* read the data as a byte stream */
  if ((*istat = read(fd, buffer, (*length * cdrcom_.KCPW))) != -1) {
    *istat = *istat / cdrcom_.KCPW;
  }
}
/*  *** CDRRVT ***

Convert from ASCII character set to internal character set.

Notes and Revisions:
This routine is here to satisfy the call, since the internal
character set is ASCII.

Parameters:
in - IN - character to be converted in ASCII
iout - OUT - character converted to internal format

*/

void cdrrvt_(int *in, int *iout) { *in = *iout; }

/*  *** CDRSTD ***

Copy the left-most N characters from the input character string
to the output word.

Parameters:
n - IN - number of characters to copy from string
string - IN - input character string
nuchar - OUT - word in which to put the n characters

*/

void cdrstd_(int *n, unsigned *string, unsigned *nuchar)
{
  if (*n > cdrcom_.KCPW) {
    return;
  }
  *nuchar = (*string >> (cdrcom_.KWRDSZ - *n * cdrcom_.KBYTEL) & ~(~0u << (*n * cdrcom_.KBYTEL)));
}
/*  *** CDRTBK ***
    Generate a traceback.  This is a dummy routine to satisfy
    externals calls.
*/
void cdrtbk_(void) {}

/*  *** CDRTIM ***
    Pause TIM seconds. (dummy)
*/

void cdrtim_(float *tim) {}

void cdrtod_(void)
{
  snprintf((char *)vcjob_.KJTIME, sizeof(vcjob_.KJTIME), "%02d:%02d:%02d\\", 11, 11, 11);
  snprintf((char *)vcjob_.KJDATE, sizeof(vcjob_.KJDATE), "%02d-%02d-%02d\\", 5, 5, 5);
}

/*  *** CDRUPR ***

Convert character string to all upper case

Parameters:
in - IN - character string to be converted
iout - OUT - character string converted to all upper case
inlen - IN - length of character string in

*/

void cdrupr_(char in[], char iout[], int *inlen)
{
  for (int i = 0; i <= *inlen; i++) {
    if (islower(in[i])) {
      iout[i] = toupper(in[i]);
    }
  }
}

/*  *** CDRWFS ***

Write LENGTH words from BUFFER to the file.

Parameters:
ifilcd - IN - FORTRAN unit number of file to write
length - IN - number of words to be written
buffer - IN - array containing the words of data to be written
istat - OUT - status indicator. 0 = no error; 1 = error

*/

void cdrwfs_(int *ifilcd, int *length, char *buffer, int *istat)
{
  /* translate fortran unit number to file descriptor */
  int fd = cdrunx_.KUNTFD[*ifilcd];

  /* if the file hasn't been opened, open it */
  if (fd == 0) {
    cdrofs_(ifilcd);
    fd = cdrunx_.KUNTFD[*ifilcd];
  }

  /* write out the data as a byte stream. */
  if ((*istat = write(fd, buffer, (*length * cdrcom_.KCPW))) == -1) {
    *istat = 1;
  }
  else {
    *istat = 0;
  }
}
