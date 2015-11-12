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
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
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
