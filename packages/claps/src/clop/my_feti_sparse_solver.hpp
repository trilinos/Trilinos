extern "C" {

  void ordmmd2_(int& n, int* , int* adj, int* invp, int* perm,
		int& iwsize, int* , int& nofsub, int& iflag);
  
  void sfinit_(int& n,   int& nnza,  int* ,  int* adj,   int* perm,
	       int* invp,int& maxsup,int& defblk,int* ,int& nnzl,
	       int& nsub,int& nsuper,int* ,int* snode, int& iwsize,
	       int* , int& iflag);
  
  void symfct_(int& n, int& nnza, int* , int* adj, int* perm,
	       int* invp, int* , int& nsuper, int* ,
	       int* snode, int& nofsub, int* , int* ,
	       int* , int& iwsize, int* , int& iflag);
  
  void bfinit_(int& nsuper, int* , int* snode, int* ,
	       int* , int& tmpsiz, int& rwsize);
  
  void blkldl_(int& nsuper, int* , int *pnode, int* ,
	       int* ,  int* ,   double *lnz, int& defblk,
	       int &asdef,  int& numZEM, int& lbdef, int *def,
	       double& tol,
	       int *iprow, int* ipcol, int& tmpsiz, double *temp,
	       int& iwsize, int* , int& rwsize, double *,
	       int& iflag);

  void blkslv_(int &nsuper, int* , int* , int *,
	       int* ,
	       double *lnx, int& defblk, int& numZEM, int& lbdef,
	       int* def, int* iprow, int* ipcol, int* perm, int* invp,
	       double *rhs, double *sol, double *temp);
  
  void blkslvn_(int &npanel, int* xpanel, int* xlindx, int *lindx,
		int* xlnz,double *lnz, int& defblk, int& numZEM, int& lbdef,
		int* def, int* iprow, int* ipcol, int* perm, int* invp,
		int &lrhs, int &nrhs, double *rhs, int &lsoln,
		double *sol, int &ltemp, double *temp);
  
  void blkns_(int &nsuper, int *, int *, int *,
	      int *, double *lnz, int &defblk, int &nrbm,
	      int &lbdef, int *def, int *ipcol, int *invp, double *ns,
	      int &numUncon, double *temp);
}
