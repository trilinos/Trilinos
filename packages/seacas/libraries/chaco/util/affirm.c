/*
 * Copyright (c) 2014, Sandia Corporation.
 * Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 * the U.S. Governement retains certain rights in this software.
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
 *     * Neither the name of Sandia Corporation nor the names of its
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
/* This software was developed by Bruce Hendrickson and Robert Leland   *
 * at Sandia National Laboratories under US Department of Energy        *
 * contract DE-AC04-76DP00789 and is copyrighted by Sandia Corporation. */

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "defs.h"

/* Record a return TRUE if answer is yes, FALSE if no. */
int       affirm(char *prompt)
{
    char      reply;		/* character typed in */
    int       done;		/* loop control */
    void      bail(char *msg, int status);

    if (prompt != NULL && (int) strlen(prompt) > 0) {
        printf("%s? ", prompt);
    }
    done = 0;
    while (!done) {
	reply = getchar();
	/* while (reply == ' ' || reply== '\n') reply= getchar(); */
	while (isspace(reply))
	    reply = getchar();

	if (reply == 'y' || reply == 'Y')
	    done = 1;
	else if (reply == 'n' || reply == 'N')
	    done = 2;
	else if (reply == 'q' || reply == 'Q')
	    done = 3;
	else if (reply == 'x' || reply == 'X')
	    done = 3;

	else {
	    printf("Valid responses begin with: y Y n N q Q x X\n");
	    if (prompt != NULL) printf("%s? ", prompt);
	    /* Flush rest of input line. */
	    while (reply != '\n')
		reply = getchar();
	}
    }
    if (done > 2)
	bail(NULL, 0);
    else if (done == 2)
	return (FALSE);
    return (TRUE);
}
