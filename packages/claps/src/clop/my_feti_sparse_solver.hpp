#ifndef MY_SPAK_H
#define MY_SPAK_H

#include "Claps_ConfigDefs.h"  // for definition of F77_FUNC

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
