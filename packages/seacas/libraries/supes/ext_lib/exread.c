/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
/*
 */

#include "fortranc.h"
#include "getline_int.h"
#include <stdio.h>
#include <string.h>
#include <unistd.h>

static int   my_getline(char *s, int len);
static char *copy_string(char *dest, char const *source, long int elements)
{
  char *d;
  for (d = dest; d + 1 < dest + elements && *source; d++, source++) {
    *d = *source;
  }
  *d = '\0';
  return d;
}

#if defined(ADDC_)
void exread_(char *prompt, char *input, FTNINT *iostat, long int PromptLength, long int InputLength)
#else
void exread(char *prompt, char *input, FTNINT *iostat, long int PromptLength, long int InputLength)
#endif

/*
************************************************************************

C     DESCRIPTION:
C     This routine prompts, reads, and echos from the standard input
C     device. For an interactive job, this would prompt for input from
C     the terminal and read (with echo) from the keyboard. For a batch
C     job, this would read from the main input file and echo to the
C     log file with the prompt string as a prefix. This routine should
C     assume the burden of assuring that the standard input and output
C     devices are properly opened.

C     FORMAL PARAMETERS:
C     PROMPT    CHARACTER       Prompt String
C     INPUT     CHARACTER       Input String
C     IOSTAT    INTEGER         I/O Status ( -1 = EOF, 0 = normal )

************************************************************************
*/
{
  static int debug = 0;
  if (debug == 1 || isatty(0) == 0 || isatty(1) == 0) {
    int icnt;

    (void)write(1, prompt, PromptLength);
    icnt = my_getline(input, InputLength);

    /* Next evaluate the error status. */
    /* For icnt < 0 indicate an error condition. */
    *iostat = (icnt >= 0) ? 0 : -1;
  }
  else {
    static char internal_prompt[128];
    char *      p = NULL;

    /* Fill line with blanks... */
    int   dlen = InputLength;
    char *ds   = input;
    while (dlen-- > 0) /* Blank out the entire string. */
      *ds++ = ' ';

    copy_string(internal_prompt, prompt, 128);
    internal_prompt[PromptLength - 1] = ' ';
    internal_prompt[PromptLength]     = '\0';

    p = getline_int(internal_prompt);

    if (p) {
      gl_histadd(p);
      int i = 0;
      /* Strip the trailing \n */
      p[strlen(p) - 1] = '\0';

      while (i < strlen(p) && i < InputLength) {
        input[i] = p[i];
        ++i;
      }
      *iostat = 0;
    }
    else {
      *iostat = -1;
    }
  }
}

static int my_getline(char *s, int len)
{
  char c;
  int  dlen;
  char nl = '\n';

  char *ds; /* A dummy argument used in nulling out s. */

  dlen = len;
  ds   = s;
  while (dlen-- > 0) /* Blank out the entire string. */
    *ds++ = ' ';

  dlen = len; /* Now take care of business. */
  for (; dlen > 0; dlen--) {
    if (read(0, &c, 1) != 1) {
      if (dlen == len) {
        (void)write(1, &nl, 1); /* We've encountered the End of File */
        return (-1);
      }
      else
        return (len - dlen + 1);
    }
    else if (c == '\n' || c == '\r')
      return (len - dlen + 1);

    *s++ = c;
  }
  /* If we get this far, we've read 'len' characters without hitting the
     end of the line.  Read until get end-of-line or end-of-file
  */
  while (read(0, &c, 1) == 1 && c != '\n' && c != '\r')
    ;
  return (len - dlen + 1);
}
