/*************************************xyt.h************************************
Module Name: xyt
Module Info: need xyt.{c,h} gs.{c,h} comm.{c,h} ivec.{c,h} error.{c,h} 
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
**************************************xyt.h***********************************/

/*************************************xyt.h************************************
File Description:
**************************************xyt.h***********************************/

/*************************************xyt.h************************************
Notes on Usage: 
**************************************xyt.h***********************************/


#ifndef _xyt_h
#define _xyt_h


typedef struct xyt_CDT *xyt_ADT;


/*************************************xyt.h************************************
Function: XYT_new()

Input :
Output:
Return: ADT ptr or NULL upon failure.
Description: This function allocates and returns an xyt handle
Usage: xyt_handle = xyt_new();
**************************************xyt.h***********************************/
extern xyt_ADT XYT_new(void);


/*************************************xyt.h************************************
Function: XYT_free()

Input : pointer to ADT.
Output:
Return:
Description: This function frees the storage associated with an xyt handle
Usage: XYT_free(xyt_handle);
**************************************xyt.h***********************************/
extern void XYT_free(xyt_ADT xyt_handle);


/*************************************xyt.h************************************
Function: XYT_factor

Input : ADT ptr,  and pointer to object
Output:
Return: 0 on failure, 1 on success
Description: This function sets the xyt solver 

xyt assumptions: given n rows of global coarse matrix (E_loc) where
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

ML beliefs/usage: move this to to ML_XYT_factor routine
   o my_ml holds address of ML struct associated w/E_loc, grid_data, grid_tag
   o grid_tag, grid_data, my_ml used in
      ML_Set_CSolve(my_ml, grid_tag, grid_data, ML_Do_CoarseDirect);
   o grid_data used in 
      A_matvec(grid_data,v,u);

Usage: 
**************************************xyt.h***********************************/
extern int XYT_factor(xyt_ADT xyt_handle,   /* prev. allocated xyt  handle */
                      int *local2global,    /* global column mapping       */
		      int n,                /* local num rows              */
		      int m,                /* local num cols              */
		      void *mylocmatvec,    /* b_loc=A_local.x_loc         */
		      void *grid_data       /* grid data for matvec        */
		      );


/*************************************xyt.h************************************
Function: XYT_solve

Input : ADT ptr, b (rhs)
Output: x (soln)
Return:
Description: This function performs x = E^-1.b
Usage: XYT_solve(xyt_handle, double *x, double *b)
**************************************xyt.h***********************************/
extern void XYT_solve(xyt_ADT xyt_handle, double *x, double *b);


/*************************************xyt.h************************************
Function: XYT_stats

Input :
Output:
Return:
Description:
**************************************xyt.h***********************************/
extern void xyt_stats(xyt_ADT xyt_handle);


/*************************************xyt.h************************************
Function: XYT_sp_1()

Input : pointer to ADT
Output: 
Return: 
Description: sets xyt parameter 1 in xyt_handle
Usage: implement later

ML_Tufo_Set_Parameters1
ML_Tufo_Set_Parameters2 ...

void XYT_sp_1(xyt_handle,parameter 1 value)
**************************************xyt.h***********************************/

#endif
