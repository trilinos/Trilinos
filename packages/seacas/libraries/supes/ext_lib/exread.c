/*
 * Copyright(C) 2008 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 * certain rights in this software
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 * * Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *           
 * * Redistributions in binary form must reproduce the above
 *   copyright notice, this list of conditions and the following
 *   disclaimer in the documentation and/or other materials provided
 *   with the distribution.
 *                         
 * * Neither the name of Sandia Corporation nor the names of its
 *   contributors may be used to endorse or promote products derived
 *   from this software without specific prior written permission.
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
/*
 * $Id: exread.c,v 1.20 2008/05/05 19:42:09 gdsjaar Exp $
 */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include "fortranc.h"
#include "getline_int.h"

static int my_getline(char *s, int len);

#if defined(ADDC_) 
void exread_(char *prompt, char *input, FTNINT *iostat,
	     long int PromptLength, long int InputLength )     
#else
void exread( char *prompt, char *input, FTNINT *iostat,
	     long int PromptLength, long int InputLength )
#endif
     
/*
************************************************************************
C
C     DESCRIPTION:
C     This routine prompts, reads, and echos from the standard input
C     device. For an interactive job, this would prompt for input from
C     the terminal and read (with echo) from the keyboard. For a batch
C     job, this would read from the main input file and echo to the
C     log file with the prompt string as a prefix. This routine should
C     assume the burden of assuring that the standard input and output
C     devices are properly openned.
C
C     FORMAL PARAMETERS:
C     PROMPT    CHARACTER       Prompt String
C     INPUT     CHARACTER       Input String
C     IOSTAT    INTEGER         I/O Status ( -1 = EOF, 0 = normal )
C
************************************************************************
*/
{
  static int debug = 0;
  if (debug == 1 || isatty(0) == 0 || isatty(1) == 0) {
    int icnt;
    
    write( 1, prompt, PromptLength );
    icnt = my_getline( input, InputLength );
    
    /* Next evaluate the error status. */
    /* For icnt <= 0 indicate an error condition. */
    *iostat = ( icnt > 0 ) ? 0 : -1;
  } else {
    static char internal_prompt[128];
    char *p = NULL;

    /* Fill line with blanks... */
    int dlen = InputLength;
    char *ds = input;
    while( dlen-- > 0 )		/* Blank out the entire string. */
      *ds++ = ' ';
    
    strncpy(internal_prompt, prompt, PromptLength);
    internal_prompt[PromptLength-1] = ' ';
    internal_prompt[PromptLength]   = '\0';

    p = getline_int(internal_prompt);
    gl_histadd(p);

    if (p) {
      int i = 0;
      /* Strip the trailing \n */
      p[strlen(p)-1] = '\0';
      
      while (i < strlen(p) && i < InputLength) {
	input[i] = p[i];
	++i;
      }
      *iostat = 0;
    } else {
      *iostat = -1;
    }
  }
}

static int my_getline(char *s, int len)
{
  char c;
  int dlen;
  char nl = '\n';

  char *ds;			/* A dummy argument used in nulling out s. */

  dlen = len;
  ds = s;
  while( dlen-- > 0 )		/* Blank out the entire string. */
    *ds++ = ' ';

  dlen = len;			/* Now take care of business. */
  for( ; dlen > 0; dlen-- ){
    if (read(0, &c, 1) != 1){
      if (dlen == len) {
	write(1,&nl,1);		/* We've encountered the End of File */
	return(-1);
      }
      else
	return(len - dlen + 1);
    }
    else
      if( c == '\n' || c == '\r' )
	return(len - dlen + 1);
   
    *s++ = c;
  }
  /* If we get this far, we've read 'len' characters without hitting the
     end of the line.  Read until get end-of-line or end-of-file
  */
  while (read(0, &c, 1) == 1 && c != '\n' && c != '\r');
  return(len - dlen + 1);
}
