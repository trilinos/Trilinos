// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#include <stdio.h>

int RAD_Const_Warn_warned = 0, RAD_Const_Warn_warnlim = 10;

 typedef struct
ADvari_head {
	void *var;	/* address of the offending variable */
	struct ADvari_head *next;
	int gcgen;
	int opno;
	/* Double Val; */
	/* mutable Double aval; */
	} ADvari_head;

 int
RAD_Const_Warn(void *v)
{
	ADvari_head *V = (ADvari_head*)v;

	if (++RAD_Const_Warn_warned <= RAD_Const_Warn_warnlim)
		fprintf(stderr, "RAD_Const_Warn var #%lx had gcgen = %d, opno = %d\n",
			(unsigned long)V->var, V->gcgen, -1 - V->opno);
	return 1;
	}
