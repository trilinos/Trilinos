/* 
 * @HEADER
 *
 * ***********************************************************************
 *
 *  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
 *                  Copyright 2012 Sandia Corporation
 *
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Government retains certain rights in this software.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the Corporation nor the names of the
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Questions? Contact Karen Devine	kddevin@sandia.gov
 *                    Erik Boman	egboman@sandia.gov
 *
 * ***********************************************************************
 *
 * @HEADER
 */


/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include "ch_input_const.h"
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#define LINE_LENGTH 200

static char line[LINE_LENGTH];	/* space to hold values */
static int offset = 0;		/* offset into line for next data */
static int break_pnt = LINE_LENGTH;	/* place in sequence to pause */
static int save_pnt;		/* place in sequence to save */

static void flush_line(ZOLTAN_FILE*);

double    read_val(
  ZOLTAN_FILE* infile,		/* file to read value from */
  int      *end_flag 		/* 0 => OK, 1 => EOL, -1 => EOF */
)
{
    double    val;		/* return value */
    char     *ptr;		/* ptr to next string to read */
    char     *ptr2;		/* ptr to next string to read */
    int       length;		/* length of line to read */
    int       length_left;	/* length of line still around */
    int       white_seen;	/* have I detected white space yet? */
    int       done;		/* checking for end of scan */
    int       i;		/* loop counter */

    *end_flag = 0;

    if (offset == 0 || offset >= break_pnt) {
	if (offset >= break_pnt) { /* Copy rest of line back to beginning. */
	    length_left = LINE_LENGTH - save_pnt - 1;
	    ptr2 = line;
	    ptr = &line[save_pnt];
	    for (i=length_left; i; i--) *ptr2++ = *ptr++;
	    length = save_pnt + 1;
	}
	else {
	    length = LINE_LENGTH;
	    length_left = 0;
	}

	/* Now read next line, or next segment of current one. */
	ptr2 = ZOLTAN_FILE_gets(&line[length_left], length, infile);

	if (ptr2 == (char *) NULL) {	/* We've hit end of file. */
	    *end_flag = -1;
	    return((double) 0.0);
	}

	if ((line[LINE_LENGTH - 2] != '\n') && (line[LINE_LENGTH - 2] != '\f')
	    && (strlen(line) == LINE_LENGTH - 1)){
	    /* Line too long.  Find last safe place in line. */
	    break_pnt = LINE_LENGTH - 1;
	    save_pnt = break_pnt;
	    white_seen = FALSE;
	    done = FALSE;
	    while (!done) {
		--break_pnt;
		if (line[break_pnt] != '\0') {
		    if (isspace((int)(line[break_pnt]))) {
			if (!white_seen) {
			    save_pnt = break_pnt + 1;
		            white_seen = TRUE;
			}
		    }
		    else if (white_seen) {
		        done= TRUE;
		    }
		}
	    }
	}
	else {
	    break_pnt = LINE_LENGTH;
	}

	offset = 0;
    }

    while (isspace((int)(line[offset])) && offset < LINE_LENGTH) offset++;
    if (line[offset] == '%' || line[offset] == '#') {
	*end_flag = 1;
	if (break_pnt < LINE_LENGTH) {
	    flush_line(infile);
	}
	return((double) 0.0);
    }

    ptr = &(line[offset]);
    val = strtod(ptr, &ptr2);

    if (ptr2 == ptr) {	/* End of input line. */
	offset = 0;
	*end_flag = 1;
	return((double) 0.0);
    }
    else {
	offset = (int) (ptr2 - line) / sizeof(char);
    }

    return(val);
}


int       read_int(
ZOLTAN_FILE *infile,		/* file to read value from */
int      *end_flag 		/* 0 => OK, 1 => EOL, -1 => EOF */
)
{
    int       val;		/* return value */
    char     *ptr;		/* ptr to next string to read */
    char     *ptr2;		/* ptr to next string to read */
    int       length;		/* length of line to read */
    int       length_left;	/* length of line still around */
    int       white_seen;	/* have I detected white space yet? */
    int       done;		/* checking for end of scan */
    int       i;		/* loop counter */

    *end_flag = 0;

    if (offset == 0 || offset >= break_pnt) {
	if (offset >= break_pnt) { /* Copy rest of line back to beginning. */
	    length_left = LINE_LENGTH - save_pnt - 1;
	    ptr2 = line;
	    ptr = &line[save_pnt];
	    for (i=length_left; i; i--) *ptr2++ = *ptr++;
	    length = save_pnt + 1;
	}
	else {
	    length = LINE_LENGTH;
	    length_left = 0;
	}

	/* Now read next line, or next segment of current one. */
	ptr2 = ZOLTAN_FILE_gets(&line[length_left], length, infile);

	if (ptr2 == (char *) NULL) {	/* We've hit end of file. */
	    *end_flag = -1;
	    return(0);
	}

	if ((line[LINE_LENGTH - 2] != '\n') && (line[LINE_LENGTH - 2] != '\f')
	    && (strlen(line) == LINE_LENGTH - 1)){
	    /* Line too long.  Find last safe place in line. */
	    break_pnt = LINE_LENGTH - 1;
	    save_pnt = break_pnt;
	    white_seen = FALSE;
	    done = FALSE;
	    while (!done) {
		--break_pnt;
		if (line[break_pnt] != '\0') {
		    if (isspace((int)(line[break_pnt]))) {
			if (!white_seen) {
			    save_pnt = break_pnt + 1;
		            white_seen = TRUE;
			}
		    }
		    else if (white_seen) {
		        done= TRUE;
		    }
		}
	    }
	}
	else {
	    break_pnt = LINE_LENGTH;
	}

	offset = 0;
    }

    while (isspace((int)(line[offset])) && offset < LINE_LENGTH) offset++;
    if (line[offset] == '%' || line[offset] == '#') {
	*end_flag = 1;
	if (break_pnt < LINE_LENGTH) {
	    flush_line(infile);
	}
	return(0);
    }

    ptr = &(line[offset]);
    val = (int) strtol(ptr, &ptr2, 10);

    if (ptr2 == ptr) {	/* End of input line. */
	offset = 0;
	*end_flag = 1;
	return(0);
    }
    else {
	offset = (int) (ptr2 - line) / sizeof(char);
    }

    return(val);
}


static void flush_line(
ZOLTAN_FILE *infile 		/* file to read value from */
)
{
    char      c;		/* character being read */

    c = ZOLTAN_FILE_getc(infile);
    while (c != '\n' && c != '\f')
	c = ZOLTAN_FILE_getc(infile);
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
