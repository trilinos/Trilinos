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

#include	<stdio.h>
#include        <string.h>
#include        "ch_input_const.h"
#include        "dr_externs.h"

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

static int input_assign_normal(ZOLTAN_FILE*, char *, int, short*),
           input_assign_inv(ZOLTAN_FILE*, char *, int, short*);

int chaco_input_assign(
ZOLTAN_FILE*     finassign,		/* input assignment file */
char     *inassignname,		/* name of input assignment file */
int       nvtxs,		/* number of vertices to output */
short    *assignment)		/* values to be printed */
{
    int       i;		/* return value */

    if (Chaco_In_Assign_Inv) {
	i = input_assign_inv(finassign, inassignname, nvtxs,
	    assignment);
    }
    else {
	i = input_assign_normal(finassign, inassignname, nvtxs,
	    assignment);
    }

    return(i);
}



static int input_assign_normal(
ZOLTAN_FILE*     finassign,		/* input assignment file */
char     *inassignname,		/* name of input assignment file */
int       nvtxs,		/* number of vertices to output */
short    *assignment)		/* values to be printed */
{
    const char *yo = "input_assign_normal";
    int       flag;		/* logical conditional */
    int       end_flag;		/* return flag from read_int() */
    int       i, j;		/* loop counter */

    DEBUG_TRACE_START(0, yo);

    /* Get the assignment vector one line at a time, checking as you go. */
    /* First read past any comments at top. */
    end_flag = 1;
    while (end_flag == 1) {
	assignment[0] = read_int(finassign, &end_flag);
    }

    if (assignment[0] < 0) {
	printf("ERROR: Entry %d in assignment file `%s' less than zero (%d)\n",
	       1, inassignname, assignment[0]);
	ZOLTAN_FILE_close(finassign);
        DEBUG_TRACE_END(0, yo);
	return (1);
    }

    if (end_flag == -1) {
	printf("ERROR: No values found in assignment file `%s'\n", inassignname);
	ZOLTAN_FILE_close(finassign);
        DEBUG_TRACE_END(0, yo);
	return (1);
    }

    flag = 0;
    if (assignment[0] > nvtxs)
	flag = assignment[1];
    for (i = 1; i < nvtxs; i++) {
	j = ZOLTAN_FILE_scanf(finassign, "%hd", &(assignment[i]));
	if (j != 1) {
	    printf("ERROR: Too few values in assignment file `%s'.\n", inassignname);
	    ZOLTAN_FILE_close(finassign);
            DEBUG_TRACE_END(0, yo);
	    return (1);
	}
	if (assignment[i] < 0) {
	    printf("ERROR: Entry %d in assignment file `%s' less than zero (%d)\n",
		   i+1, inassignname, assignment[i]);
	    ZOLTAN_FILE_close(finassign);
            DEBUG_TRACE_END(0, yo);
	    return (1);
	}
	if (assignment[i] > nvtxs) {	/* warn since probably an error */
	    if (assignment[i] > flag)
		flag = assignment[i];
	}
    }

    if (flag && Debug_Chaco_Input) {
	printf("WARNING: Possible error in assignment file `%s'\n", inassignname);
	printf("         More assignment sets (%d) than vertices (%d)\n", flag, nvtxs);
    }

    /* Check for spurious extra stuff in file. */
    flag = FALSE;
    end_flag = 0;
    while (!flag && end_flag != -1) {
	read_int(finassign, &end_flag);
	if (!end_flag)
	    flag = TRUE;
    }
    if (flag && Debug_Chaco_Input) {
	printf("WARNING: Possible error in assignment file `%s'\n", inassignname);
	printf("         Numerical data found after expected end of file\n");
    }

    ZOLTAN_FILE_close(finassign);
    DEBUG_TRACE_END(0, yo);
    return (0);
}


static int input_assign_inv(
ZOLTAN_FILE* finassign,		/* input assignment file */
char     *inassignname,		/* name of input assignment file */
int       nvtxs,		/* number of vertices to output */
short    *assignment)		/* values to be printed */
{
    const char *yo = "input_assign_inv";
    int       set;		/* set number being read */
    int       size;		/* number of vertices in set */
    int       total;		/* total number of vertices read */
    int       done;		/* have I hit end of file yet? */
    int       end_flag;		/* return flag from read_int() */
    int       i, j, k;		/* loop counter */

    DEBUG_TRACE_START(0, yo);

    /* Get the assignment vector one line at a time, checking as you go. */

    /* Initialize assignment to help error detection. */
    for (i = 0; i < nvtxs; i++) {
	assignment[i] = -1;
    }

    /* First read past any comments at top. */
    total = 0;
    set = 0;
    end_flag = 1;
    while (end_flag == 1) {
	size = read_int(finassign, &end_flag);
    }

    if (end_flag == -1) {
	printf("ERROR: In assignment file `%s'\n", inassignname);
	printf("       No values found\n");
	ZOLTAN_FILE_close(finassign);
        DEBUG_TRACE_END(0, yo);
	return (1);
    }

    if (size < 0) {
	printf("ERROR: In assignment file `%s'\n", inassignname);
	printf("       Size of set %d less than zero (%d)\n", set, size);
	ZOLTAN_FILE_close(finassign);
        DEBUG_TRACE_END(0, yo);
	return (1);
    }

    if (total + size > nvtxs) {
	printf("ERROR: In assignment file `%s'\n", inassignname);
	printf("       Total set sizes greater than nvtxs (%d)\n", nvtxs);
	ZOLTAN_FILE_close(finassign);
        DEBUG_TRACE_END(0, yo);
	return (1);
    }


    done = FALSE;
    while (!done && total < nvtxs) {
	for (i = 1; i <= size; i++) {
	    j = ZOLTAN_FILE_scanf(finassign, "%d", &k);
	    if (j != 1) {
	        printf("ERROR: Too few values in assignment file `%s'.\n",
		    inassignname);
	        ZOLTAN_FILE_close(finassign);
                DEBUG_TRACE_END(0, yo);
	        return (1);
	    }

	    if (k <= 0 || k > nvtxs) {
		printf("ERROR: In assignment file `%s'\n", inassignname);
		printf("       Entry %d of set %d invalid (%d)\n", 
		    total + i, set, k);
		ZOLTAN_FILE_close(finassign);
                DEBUG_TRACE_END(0, yo);
		return (1);
	    }

	    if ((int) assignment[k - 1] != -1) {
		printf("ERROR: In assignment file `%s'\n", inassignname);
		printf("       Vertex %d assigned to multiple sets\n", k);
		ZOLTAN_FILE_close(finassign);
                DEBUG_TRACE_END(0, yo);
		return (1);
	    }

	    assignment[k - 1] = (short) set;
	}

	total += size;
	j = ZOLTAN_FILE_scanf(finassign, "%d", &size);
	++set;
	if (j != 1) {
	    if (total != nvtxs) {
	        printf("ERROR: Too few values in assignment file `%s'.\n",
		    inassignname);
		ZOLTAN_FILE_close(finassign);
                DEBUG_TRACE_END(0, yo);
		return (1);
	    }
	    else {
		done = TRUE;
		size = 0;
	    }
	}


	if (size < 0) {
	    printf("ERROR: In assignment file `%s'\n", inassignname);
	    printf("       Size of set %d less than zero (%d)\n", set, size);
	    ZOLTAN_FILE_close(finassign);
            DEBUG_TRACE_END(0, yo);
	    return (1);
	}

	if (total + size > nvtxs) {
	    printf("ERROR: In assignment file `%s'\n", inassignname);
	    printf("       Total set sizes greater than nvtxs (%d)\n", nvtxs);
	    ZOLTAN_FILE_close(finassign);
            DEBUG_TRACE_END(0, yo);
	    return (1);
	}
    }

    ZOLTAN_FILE_close(finassign);
    DEBUG_TRACE_END(0, yo);
    return (0);
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
