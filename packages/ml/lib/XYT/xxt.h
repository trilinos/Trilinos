/*************************************xxt.h************************************
Module Name: xxt
Module Info: need xxt.{c,h} gs.{c,h} comm.{c,h} ivec.{c,h} error.{c,h} 
             bss_malloc.{c,h} bit_mask.{c,h} queue.{c,h} stack.{c,h}

author:  Henry M. Tufo III
e-mail:  hmt@asci.uchicago.edu
contact:
+--------------------------------+--------------------------------+
|MCS Division - Building 221     |Department of Computer Science  |
|Argonne National Laboratory     |Ryerson 152                     |
|9700 S. Cass Avenue             |The University of Chicago       |
|Argonne, IL  60439              |Chicago, IL  60637              |
|(630) 252-5354/5986 ph/fx       |(773) 702-6019/8487 ph/fx       |
+--------------------------------+--------------------------------+

Last Modification: 3.27.00
**************************************xxt.h***********************************/

/*************************************xxt.h************************************
File Description:
**************************************xxt.h***********************************/

/*************************************xxt.h************************************
Notes on Usage: 
**************************************xxt.h***********************************/


#ifndef _xxt_h
#define _xxt_h


typedef struct xxt_CDT *xxt_ADT;


/*************************************xxt.h************************************
Function: XXT_new()

Input :
Output:
Return: ADT ptr or NULL upon failure.
Description: This function allocates and returns an xxt handle
Usage: xxt_handle = xxt_new();
**************************************xxt.h***********************************/
extern xxt_ADT XXT_new(void);


/*************************************xxt.h************************************
Function: XXT_free()

Input : pointer to ADT.
Output:
Return:
Description: This function frees the storage associated with an xxt handle
Usage: XXT_free(xxt_handle);
**************************************xxt.h***********************************/
extern void XXT_free(xxt_ADT xxt_handle);


/*************************************xxt.h************************************
Function: XXT_factor

Input : ADT ptr,  and pointer to object
Output:
Return: 0 on failure, 1 on success
Description: This function sets the xxt solver 

xxt assumptions: given n rows of global coarse matrix (E_loc) where
   o global dofs N = sum_p(n), p=0,P-1 
   (i.e. row dist. with no dof replication)
   (5.21.00 will handle dif replication case)
   o m is the number of columns in E_loc (m>=n)
   o local2global holds global number of column i (i=0,...,m-1)
   o local2global holds global number of row    i (i=0,...,n-1)
   o mylocmatvec performs E_loc . x_loc where x_loc is an vector of
   length m in 1-1 correspondence with local2global
   (note that gs package takes care of communication).
   (note do not zero out upper m-n entries!)
   o mylocmatvec(void *grid_data, double *in, double *out)

ML beliefs/usage: move this to to ML_XXT_factor routine
   o my_ml holds address of ML struct associated w/E_loc, grid_data, grid_tag
   o grid_tag, grid_data, my_ml used in
      ML_Set_CSolve(my_ml, grid_tag, grid_data, ML_Do_CoarseDirect);
   o grid_data used in 
      A_matvec(grid_data,v,u);

Usage: 
**************************************xxt.h***********************************/
extern int XXT_factor(xxt_ADT xxt_handle,   /* prev. allocated xxt  handle */
                      int *local2global,    /* global column mapping       */
		      int n,                /* local num rows              */
		      int m,                /* local num cols              */
		      void *mylocmatvec,    /* b_loc=A_local.x_loc         */
		      void *grid_data       /* grid data for matvec        */
		      );


/*************************************xxt.h************************************
Function: XXT_solve

Input : ADT ptr, b (rhs)
Output: x (soln)
Return:
Description: This function performs x = E^-1.b
Usage: XXT_solve(xxt_handle, double *x, double *b)
**************************************xxt.h***********************************/
extern void XXT_solve(xxt_ADT xxt_handle, double *x, double *b);


/*************************************xxt.h************************************
Function: XXT_sp_1()

Input : pointer to ADT
Output: 
Return: 
Description: sets xxt parameter 1 in xxt_handle
Usage: implement later

ML_Tufo_Set_Parameters1
ML_Tufo_Set_Parameters2 ...

void XXT_sp_1(xxt_handle,parameter 1 value)
**************************************xxt.h***********************************/

#endif
