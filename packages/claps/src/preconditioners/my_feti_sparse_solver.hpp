//@HEADER
// ************************************************************************
// 
//         Claps: A Collection of Domain Decomposition Preconditioners
//                and Solvers
//         Copyright (2006) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Clark R. Dohrmann (crdohrm@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#ifndef MY_SPAK_H
#define MY_SPAK_H

#include "Claps_ConfigDefs.hpp"  // for definition of F77_FUNC

#define ORDMMD2_F77 F77_FUNC(ordmmd2,ORDMMD2)
#define SFINIT_F77 F77_FUNC(sfinit,SFINIT)
#define SYMFCT_F77 F77_FUNC(symfct,SYMFCT)
#define BFINIT_F77 F77_FUNC(bfinit,BFINIT)
#define BLKLDL_F77 F77_FUNC(blkldl,BLKLDL)
#define BLKSLV_F77 F77_FUNC(blkslv,BLKSLV)
#define BLKSLVN_F77 F77_FUNC(blkslvn,BLKSLVN)
#define BLKNS_F77 F77_FUNC(blkns,BLKNS)
extern "C" {
  
  void ORDMMD2_F77(int& n, int* , int* adj, int* invp, int* perm,
		   int& iwsize, int* , int& nofsub, int& iflag);
  
  void SFINIT_F77(int& n,   int& nnza,  int* ,  int* adj,   int* perm,
		  int* invp,int& maxsup,int& defblk,int* ,int& nnzl,
		  int& nsub,int& nsuper,int* ,int* snode, int& iwsize,
		  int* , int& iflag);
  
  void SYMFCT_F77(int& n, int& nnza, int* , int* adj, int* perm,
		  int* invp, int* , int& nsuper, int* ,
		  int* snode, int& nofsub, int* , int* ,
		  int* , int& iwsize, int* , int& iflag);
  
  void BFINIT_F77(int& nsuper, int* , int* snode, int* ,
		  int* , int& tmpsiz, int& rwsize);
  
  void BLKLDL_F77(int& nsuper, int* , int *pnode, int* ,
		  int* ,  int* ,   double *lnz, int& defblk,
		  int &asdef,  int& numZEM, int& lbdef, int *def,
		  double& tol,
		  int *iprow, int* ipcol, int& tmpsiz, double *temp,
		  int& iwsize, int* , int& rwsize, double *,
		  int& iflag);

  void BLKSLV_F77(int &nsuper, int* , int* , int *,
		  int* ,
		  double *lnx, int& defblk, int& numZEM, int& lbdef,
		  int* def, int* iprow, int* ipcol, int* perm, int* invp,
		  double *rhs, double *sol, double *temp);
  
  void BLKSLVN_F77(int &npanel, int* xpanel, int* xlindx, int *lindx,
		   int* xlnz,double *lnz, int& defblk, int& numZEM, int& lbdef,
		   int* def, int* iprow, int* ipcol, int* perm, int* invp,
		   int &lrhs, int &nrhs, double *rhs, int &lsoln,
		   double *sol, int &ltemp, double *temp);
  
  void BLKNS_F77(int &nsuper, int *, int *, int *,
		 int *, double *lnz, int &defblk, int &nrbm,
		 int &lbdef, int *def, int *ipcol, int *invp, double *ns,
		 int &numUncon, double *temp);
}
#endif
