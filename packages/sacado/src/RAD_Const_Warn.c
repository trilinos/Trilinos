/*
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2007) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER
*/

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
