/*
 * Copyright(C) 2009 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 * certain rights in this software.
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
 */
/* $Id: help.h,v 1.3 2009/03/25 12:46:02 gdsjaar Exp $ */

/* Exit status returned by help() */
#define	H_FOUND		0	/* found the keyword */
#define	H_NOTFOUND	1	/* didn't find the keyword */
#define	H_ERROR		(-1)	/* didn't find the help file */

#ifndef	IO_SUCCESS	/* DECUS or VMS C will have defined these already */
#define	IO_SUCCESS	0
#endif
#ifndef	IO_ERROR
#define	IO_ERROR	1
#endif

extern void FreeHelp();		/* use this if you need memory */

typedef int BOOLEAN;

#define MAX_LINE_LEN	132
#define MAX_TOKENS	20

#define TRUE	1
#define FALSE	0
#define NO_CARET (-1)

#ifndef HELPFILE
#define HELPFILE "xhelp.gih"
#endif

#ifndef PROMPT
#define PROMPT	 "XHelp> "
#endif

enum DATA_TYPES {
	INT, CMPLX
};


#define END_OF_COMMAND (c_token >= num_tokens || equals(c_token,";"))

#ifndef MAIN
/* Some key global variables */
extern BOOLEAN screen_ok;
extern int	c_token, num_tokens;

extern char input_line[MAX_LINE_LEN+1];
#endif
#include <stdio.h>

struct cmplx {
	double real, imag;
};

struct value {
	enum DATA_TYPES type;
	union {
		int int_val;
		struct cmplx cmplx_val;
	} v;
};

struct lexical_unit {	/* produced by scanner */
	BOOLEAN is_token;	/* true if token, false if a value */ 
	struct value l_val;
	int start_index;	/* index of first char in token */
	int length;			/* length of token in chars */
};

