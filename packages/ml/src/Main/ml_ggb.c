
/*************************************************************************************       
  HAIM: GLOBAL EIGENVALUE CALCULATIONS FOR MULTIGRID           
  GLOBAL BASIS (GB) METHOD / GENERALIZED GLOBAL BASIS (GGB) METHOD USING ARPACK                      
**************************************************************************************/     
/*                                                            
   THE IDEA:                                                  
        We want to solve     (STS)*x = lam*x                      
        S - Smoothing iteration matrix
	T - Multilevel projector
		
        + We solve the problem by supplying ARPACK with the function that
	  produces Matrix-vector product 
	+ We use simple non-symmetric ARPACK mode                            
                                                              
	*/                                                            
#include "ml_ggb.h"
#include "ml_lapack.h"


/*****************************************************************************************
*****************************************************************************************/

void ML_ARPACK_GGB( struct ML_Eigenvalue_Struct *eigen_struct,ML *ml,
		     struct ML_CSR_MSRdata *mydata)
   
{

  int i ,  j, MatSize, level;
  
  /* Eigenvalue definitions */
  int      iparam[11];
  int      nev, ncv, info, mode, nconv;
  double   tol, tm, tmp_tol;
  char     bmat[2], which[3];
  ML_Operator *Amat;
  


  /* --------- Execution begins  ---------------- */

  printf("\n");
  ML_print_line("=", 80);
  printf("\tStarting Eigenvalue Calculation");
  tm = GetClock();

    /* Some matrix information */
    
  level         = ml->ML_finest_level;
  Amat          = &(ml->Amat[level]);
 
  /* Set parameters for ARPACK from the input file */

  iparam[MAX_ITRS] = eigen_struct->Max_Iter;
  nev              = eigen_struct->Num_Eigenvalues;
  ncv              = eigen_struct->Arnoldi;
  tol              = eigen_struct->Residual_Tol;
  printf("nev is %d and ncv is %d\n",nev,ncv);
    
     
  /* Set parameters for ARPACK: (2)then those that are fixed for MPSalsa */
  /* Setting ARPACK to: simple,nonsymetric, and finding Large Magnitude spectrum */ 
    
  which[0] = 'L';
  which[1] = 'M';
  which[2] = '\0';
  
  bmat[0] = 'I';
  bmat[1] = '\0';
  
  iparam[3] = 1;
  iparam[4] = 0;
  iparam[5] = 1;    /* We will check for convergence ourselves */
  iparam[7] = 0;
  iparam[8] = 0;
  iparam[9] = 0;
  iparam[10] = 0;
  
  iparam[SHIFTS] = 1;
  iparam[MODE] = 1;
    

  /* ARNOLDI driver. On output nconv is the number of converged eigenvalues */ 
  ML_ARPACK_driver(which, bmat, iparam, mode,
		     nev, ncv, tol, ml, mydata);
  
  printf("Time for eigenvalue computations is %g (sec.)\n",GetClock()-tm);
  
}


/********************************************************************/
/********************************************************************/
/********************************************************************/



void  ML_ARPACK_driver(char which[],
                              char bmat[], int iparam[], int mode, 
                              int nev, int ncv, double tol,  ML *ml, 
                              struct ML_CSR_MSRdata *mydata)
{

  int        i, j, kk, ldv, lworkl;
  int        nloc, nloc_max, nloc2, ido, flag, counter;
  int        count, nconv, ierr, info, comm;
  int        ipntr[14], m=0;
  int        one = 1, dummy1, dummy2, dummy3, dummy4;
  double     a1 , a2  , lamR, lamI;
  char       string[4];
  double     *v, *workl, *workd, *workev, *d, *resid;         /* Pointers used by ARPACK */  
  int        *select, rvec, *work;

  double     *vecx, *vecy, *rhs, *rhs1;                       /* Dummy Pointers */
  double     **eigvec;                                        /* Pointer to hold eigenvectors */ 
 
  FILE       *ifp;
  ML_Operator *Amat;

  /********************************  Begin ************************************/
  /*
    ifp   = fopen("CheckEigen1.m", "w");
  */

  /******************************************************
   * A standard eigenvalue problem is solved (BMAT = 'I').
   * NEV is the number of eigenvalues to be approximated.
   * NCV is the number of Krylov vectors kept. WHICH
   * determine what part of the spectrum is obtained.
   * The following conditions must be satisfied:
   *                  N <= MAXN,
   *                NEV <= MAXNEV,
   *             NEV + 2 <= NCV <= MAXNCV
   *
   ******************************************************/
  Amat = &(ml->Amat[ml->ML_finest_level]);
  nloc       = Amat->outvec_leng;
  ldv        = nloc;

  /******************************************************
   * The work array WORKL is used in P??AUPD as
   * workspace.  Its dimension LWORKL is set as
   * illustrated below.  The parameter TOL determines
   * the stopping criterion.  If TOL<=0, machine
   * precision is used.  The variable IDO is used for
   * reverse communication and is initially set to 0.
   * Setting INFO=0 indicates that a random vector is
   * generated in PNAUPD to start the Arnoldi iteration.
   ******************************************************/
  lworkl = 3*ncv*(ncv+2);
  ido    = 0;
  info   = 0;  
  
  /******************************************************
   * Use exact shifts with respect to the current Hessenberg
   * matrix (iparam[SHIFTS] = 1) where IPARAM[MAX_ITRS] is
   * the maximum number of Arnoldi iterations allowed.
   * Mode 1 of P??AUPD is used (IPARAM[MODE] = 1). For
   * details, see the documentation in P??AUPD.
   ******************************************************/

  nloc2  =  nloc; /* + data_org[AZ_N_external];  HAIM: nloc2=nloc */
  nloc_max = nloc2;


  select = (int    *) ML_allocate(ncv*sizeof(int));
  work   = (int    *) ML_allocate(ncv*sizeof(int));
  vecx   = (double *) ML_allocate(nloc2*sizeof(double));
  vecy   = (double *) ML_allocate(nloc2*sizeof(double));
  rhs    = (double *) ML_allocate(nloc2*sizeof(double));
  rhs1   = (double *) ML_allocate(nloc2*sizeof(double));
  d      = (double *) ML_allocate(4*ncv*sizeof(double));
  resid  = (double *) ML_allocate(nloc_max*sizeof(double));
  workd  = (double *) ML_allocate(3*nloc2*sizeof(double));
  workev = (double *) ML_allocate(3*ncv*sizeof(double));
  workl  = (double *) ML_allocate(lworkl*sizeof(double));
  v      = (double *) ML_allocate(ncv*ldv*sizeof(double));
 
  
  if (v == NULL) {
    fprintf(stderr,"Not enough space to allocate workl\n");
    exit(1);
  }

  
  /* Initialize dummy structures */
  for (j = 0; j < nloc2 ; j++ )  rhs[j] = rhs1[j] = vecx[j] = vecy[j] = 0.0;
  
  /* Initial guess is generated by ARPACK with info=0 */
  /* main loop */

  for (j = 0; j < 4*ncv ; j++ ) d[j] = 0.0; /* zero d[:] */
  flag = 1;
  count = 0;
  
  while ( flag == 1 ) {
    
    /*****************************************************
     * Repeatedly call P??AUPD and take actions indicated
     * by parameter IDO until either convergence is indicated
     * or the maximum iterations have been exceeded.
     *****************************************************/
    printf("nev is %d and ncv is %d\n",nev,ncv);
    dnaupd_(&ido, bmat, &nloc, which, &nev, &tol, resid,
	    &ncv, v, &ldv, iparam, ipntr, workd, workl,
	    &lworkl, &info );
    printf("after dnaupd\n");

    if ( (ido == -1) || (ido == 1) ) {
      printf("in if\n");
      count++;
      if (ml->comm->ML_mypid==0) printf("\t    Eigensolver iteration: %d\n",count);

      /***********************************************
       * matrix vector multiplication (using inverse)
       *   workd[ipntr[1]-1] <-- OP * workd[ipntr[0]-1]
       ***********************************************/
     
      /************************************/
      /* Compute (I-inv(P)*K)v = lam * v  */
      /************************************/
      
      /* First perform K*v = rhs */
      ML_Operator_Apply(Amat, Amat->invec_leng, &workd[ipntr[0]-1], 
			Amat->outvec_leng, rhs);

      /*******************************************************/
      /* Do ML cycle to find the iteration matrix eigenpairs */
      /* *****************************************************/
      
      /* Second perform inv(P)*rhs = rhs1 */
	ML_Solve_MGV(ml ,rhs , rhs1);
	
	/* Third perform  workd[ipntr[1]-1] <- v-rhs1 */ 
	for (kk = 0 ; kk < nloc2 ; kk++) 
	  workd[ipntr[1]+kk-1] = workd[ipntr[0]+kk-1] - rhs1[kk];
      
   }
      
      
    else flag = 0;
    printf("at the bottom\n");
  }


  if ( info < 0 ) {
      fprintf(stderr,"\nError with _naupd, info = %d\n",info);
      fprintf(stderr,"Check documentation in the end of file ml_ggb.c \n\n");
      nconv = 0;
      if ( info == -9999 ) {
        fprintf(stderr,"Size of Arnoldi factorization:%d\n",iparam[4]);
        fprintf(stderr,"Decrease ncv=%d to %d & rerun\n",ncv, iparam[4]);
      }
  }
  else {
    /***********************************************
     * No fatal errors occurred.
     * Post-Process using PSNEUPD.
     *
     * Extract computed eigenvalues.  Eigenvectors
     * may also be computed by setting rvec = 1.
     **********************************************/

    rvec = 1;
    ierr = 0;
    sprintf(string,"A");
    
    ml_dneupc__(&rvec, string, select, d, v, &ldv,
		  workev, bmat, &nloc, which, &nev,
		  &tol, resid, &ncv, iparam, ipntr, workd, workl,
		  &lworkl, &ierr);
    
    
    /*----------------------------------------------
      | The real part of the eigenvalue is returned   |
      | in the first column of the two dimensional    |
      | array D, and the imaginary part is returned   |
      | in the second column of D.  The corresponding |
      | eigenvectors are returned in the first NEV    |
      | columns of the two dimensional array V if     |
      | requested.  Otherwise, an orthogonal basis    |
      | for the invariant subspace corresponding to   |
      | the eigenvalues in D is returned in V.        |
      -----------------------------------------------*/
    
    if ( ierr !=  0) {
      /*-----------------------------------
        | Error condition:                   |
        | Check the documentation of PDNEUPC.|
        ------------------------------------*/
	fprintf(stderr,"\nError with _neupc, info = %d", ierr);
	fprintf(stderr,"\nCheck the documentation of _neupc.\n\n");
	exit(1);
      }

      else {
	nconv       =  iparam[4];    /* Number of eigenvalues to converge */

	/* Checking to see if real or complex case */
	  counter = 0;
	  for (j = 0; j < nconv ; j++ ) {                                    
	  
	    if (d[j+ncv] != 0.0) {
	      printf("\n");
	      ML_print_line("=", 80);
	      printf("\t\t GGB: We have complex eigenvectors \n");
	      ML_print_line("=", 80);
	      counter = 1;
	      break;
	    }
	  }
	  
	  if (counter == 0)  {
	    printf("\n");
	    ML_print_line("=", 80);
	    printf("\t\t GGB: We have real eigenvectors \n");
	    ML_print_line("=", 80);
	  }
	  
	  
	  eigvec = (double **) ML_allocate((nconv+1)*sizeof(double *));
	  
	  /* Dynamic memory allocation */                                                                           
	  for (kk = 0; kk < nconv+1 ; kk++ ) {
	    eigvec[kk] = (double *) ML_allocate(nloc2*sizeof(double));      
	    for (j = 0; j < nloc2; j++) eigvec[kk][j] = 0.;
	  }           
	
	  for (j = 0; j < nconv ; j++ ) {
	  
	  /*--------------------------
          | Compute the residual norm |
          |                           |
          |   ||  J*x - lambda*x ||   |
          |                           |
          | for the NCONV accurately  |
          | computed eigenvalues and  |
          | eigenvectors.  (iparam(5) |
          | indicates how many are    |
          | accurate to the requested |
          | tolerance)                |
	  ---------------------------*/
	  lamR = d[j];        /* Real Part of Eigenvalue */
	  lamI = d[j+ncv];    /* Imaginary Part of Eigenvalue */
	

	  if (lamI == 0.0){
	    /*-------------------
	      | Ritz value is real |
	      --------------------*/
	    /*
           *  d[j]     : j-th eigenvalue
           *
      	   *  v[j*ldv] : j-th eigenvector
	   */
	    
	    for (kk =0 ; kk < nloc2 ; kk++ ) {
	      rhs[kk]  = rhs1[kk] = 0.0; 
	      eigvec[j][kk] =  v[j*ldv+kk];
	    
	    }

	    /*   fprintf(ifp,"\t\t Q(%d,%d)= %20.13f; \n",kk+1,j+1,GLB->eigvecR[j][kk]); */
	    ML_Operator_Apply(Amat, Amat->invec_leng, &v[j*ldv],
			Amat->outvec_leng, rhs);

	    
	    ML_Solve_MGV(ml ,rhs , rhs1);
	    
	    for (kk =0 ; kk < nloc2 ; kk++ ) rhs1[kk] = (1-lamR)*v[j*ldv+kk] - rhs1[kk];
	    
	    
	    
	    a1 = sqrt(ddot_(&nloc2, rhs1, &one, rhs1, &one));
	    d[j+2*ncv] = a1; /*sqrt(lamR*lamR);*/
	    printf("a1 = %f \n",a1);
	    
	  }
	  
	  else{
	    /*----------------------
            | Ritz value is complex |
            -----------------------*/
          /*
           *  d[j]     : real part j-th eigenvalue
           *  d[j+ncv] : imag part j-th eigenvalue
           *
           *  v[j*ldv]     : real part j-th eigenvector
           *  v[(j+1)*ldv] : imag part j-th eigenvector
           */
	  
	  for (kk =0 ; kk < nloc2 ; kk++ ) {
	    rhs[kk] = rhs1[kk] = vecx[kk] = vecy[kk] = 0.0;
	    
	    eigvec[j][kk]   = v[j*ldv+kk]; 
	    eigvec[j+1][kk] = v[(j+1)*ldv+kk];
	
	    /*
	    fprintf(ifp,"\t\t Q(%d,%d)= %20.13f; \n",kk+1,j+1,GLB->eigvecR[j][kk]);
	    fprintf(ifp,"\t\t Q(%d,%d)= %20.13f; \n",kk+1,j+2,GLB->eigvecR[j+1][kk]);
	    */
	  }
	    ML_Operator_Apply(Amat, Amat->invec_leng, &v[j*ldv],
			      Amat->outvec_leng, vecx);
	    ML_Operator_Apply(Amat, Amat->invec_leng, &v[(j+1)*ldv],
			      Amat->outvec_leng, vecy);

	    
	    a1 = 1 - lamR;
	    a2 = 0.0;
	    
	    
	    
	    ML_Solve_MGV(ml ,vecx ,rhs);
	    ML_Solve_MGV(ml ,vecy ,rhs1);
	    
	    for (kk =0 ; kk < nloc ; kk++ ) {
	      rhs[kk]  = a1*v[j*ldv+kk]     + lamI*v[(j+1)*ldv+kk] - rhs[kk];
	      rhs1[kk] = a1*v[(j+1)*ldv+kk] - lamI*v[j*ldv+kk]     - rhs1[kk];
	    }
	      
	   
	    	  
	    
	    a2 = sqrt(ddot_(&nloc2, rhs1, &one, rhs1, &one) +
		      ddot_(&nloc2, rhs, &one, rhs, &one));
	    
	    
	    d[j+2*ncv] =   a2; //  /dlapy2_(&lamR,&lamI);
	    
	    d[j+1+2*ncv] = a2;
	    
	  
	    d[j+1]       =  d[j];
	    d[j+1+ncv]   = -d[j+ncv];
	    j = j + 1;
	  
	    
	  /* end of Rayleigh quotient stuff */
          /*-----------------------
            | Ritz value is complex. |
            | Residual of one Ritz   |
            | value of the conjugate |
            | pair is computed.      |
            ------------------------*/
	  }
      }
      

      /*   Display computed residuals   */

	  comm=0; dummy1 = 6; dummy2 = 3; dummy3 = ncv; dummy4 = -6;
	  ml_c_pdmout__(&comm, &dummy1, &nconv, &dummy2, d, &dummy3, &dummy4);

      }

    
    /* Haim: Convert Information into CSR Format */
    ML_GGB_2_CSR (eigvec, nconv, nloc2, mydata);
    
    
    /*  Print additional convergence information */
    
    if ( info == 1 ){
      printf("\nMaximum number of iterations reached.\n");
    }
    else if ( info == 3 ){
      printf("\nNo shifts could be applied during implicit\n");
      printf("Arnoldi update, try increasing NCV.\n\n");
    }
    else {
      printf("\nEigenvalue Calculation Summary\n\n");
      printf("The global size of the matrix = %d\n", nloc);
      printf("The number of Ritz values requested is %d\n", nev);
      printf("The number of Arnoldi vectors generated (NCV) is %d\n",ncv);
      printf("What portion of the spectrum: %s\n", which);
      printf("The number of converged Ritz values is %d\n",nconv);
      printf("Number of Implicit Arnoldi update iterations taken is %d\n",
              iparam[2]-1);
      printf("The number of OP*x is %d\n",iparam[8]);
      printf("The convergence criterion is %e\n", tol);

    }
  }
 

  ML_free((void *) select);
  ML_free((void *) work);
  ML_free((void *) vecx);
  ML_free((void *) vecy);
  ML_free((void *) rhs1);
  ML_free((void *) rhs);
  ML_free((void *) d);
  ML_free((void *) resid);
  ML_free((void *) workd);
  ML_free((void *) workev);
  ML_free((void *) workl);
  ML_free((void *) v);
  for (kk = 0; kk < nconv+1 ; kk++ ) {
    ML_free(eigvec[kk]);
  }           
  ML_free((void *) eigvec);

 


    /*
      fclose(ifp);  
    */

} /* end ML_ARPACK_driver */

/******************************************************************************/

void ML_GGB_2_CSR (double **eigvec, int nconv, int MatSize,
		   struct ML_CSR_MSRdata  *mydata)


{
  int          nrows,   ncolumns,  nnz;
  int          *rowptr, *columns;
  int          i, j , count;
  double       *values;
 
  FILE          *fp, *fp1, *eig;
  


  /******************************  begin *****************************/
  /*
    fp  = fopen("Rowptr.m","w");
    fp1 = fopen("Val_Col.m","w");
    eig = fopen("EIG_CHECK.m","w");
  */

  ncolumns = nconv;
  nrows    = MatSize;
  nnz      = nrows * ncolumns;
  
  /*****************************************************
    MEMORY ALLOCATION FOR CSR - SPARSE MATRIX STORAGE
   ******************************************************/
  
  rowptr = (int *) ML_allocate((nrows+1)*sizeof(int));
  columns = (int *) ML_allocate((nnz+1)*sizeof(int));
  values = (double *) ML_allocate((nnz+1)*sizeof(double));
  

  rowptr[0] = 0;
  count     = 0;




  for (i = 0; i < nrows; i++) {       /* Looping over the rows */
    
    rowptr[i+1] = rowptr[i] + nconv;

    
    for (j = 0; j < ncolumns; j++) {
    
      columns[count] = j;
      values[count]  = eigvec[j][i];
      count = count + 1;
    }
    
    
  }

  fprintf(stdout,"\n\t ********************************************************");
  fprintf(stdout,"\n\t   GGB:    Q PROLONGATION MATRIX (GLOBAL BASIS)         ");
  fprintf(stdout,"\n\t          ---------------------------------------    ");
  fprintf(stdout,"\n\t                 MATRIX SIZE:   %d * %d             ",nrows,ncolumns);
  fprintf(stdout,"\n\t                 NON ZEROS        = %d              ",nnz);
  fprintf(stdout,"\n\t ********************************************************\n");  
  
  
  
  /* DEBUGING */
  /*
  for (i = 0; i < nnz; i++)
    fprintf(fp1,"%f    %d \n",values[i],columns[i]);
  
  for (i = 0; i < nrows+1; i++) 
    fprintf(fp,"%d \n",rowptr[i]);
  
  for (i = 0; i < ncolumns; i++) {
    fprintf(eig,"EIG NUM = %d\n",i+1);
    for (j = 0; j < nrows; j++) {
      fprintf(eig,"\t %f\n ",eigvec[i][j]);
    }
  }
  */
 


  /*************************************************************************
     FEEDING THE NESECCARRY DATA FOR THE CSR SPARSE FORMAT INTO STRUCTURE
   **************************************************************************/
      
  mydata->Nrows   = nrows;
  mydata->Ncols   = ncolumns;
  mydata->Nnz     = nnz;
  mydata->rowptr  = rowptr;
  mydata->columns = columns;
  mydata->values  = values;


  
    
  /*
  fclose(fp);
  fclose(fp1);
  fclose(eig);
  */

} /* end function */

/*********************************************/
/*  ARPACK Documentation   if info < 0       */ 
/*********************************************/
/* Haim: taken from dnaupd.f file 

c  INFO    Integer.  (INPUT/OUTPUT)
c          If INFO .EQ. 0, a randomly initial residual vector is used.
c          If INFO .NE. 0, RESID contains the initial residual vector,
c                          possibly from a previous run.
c          Error flag on output.
c          =  0: Normal exit.
c          =  1: Maximum number of iterations taken.
c                All possible eigenvalues of OP has been found. IPARAM(5)  
c                returns the number of wanted converged Ritz values.
c          =  2: No longer an informational error. Deprecated starting
c                with release 2 of ARPACK.
c          =  3: No shifts could be applied during a cycle of the 
c                Implicitly restarted Arnoldi iteration. One possibility 
c                is to increase the size of NCV relative to NEV. 
c                See remark 4 below.
c          = -1: N must be positive.
c          = -2: NEV must be positive.
c          = -3: NCV-NEV >= 2 and less than or equal to N.
c          = -4: The maximum number of Arnoldi update iteration 
c                must be greater than zero.
c          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
c          = -6: BMAT must be one of 'I' or 'G'.
c          = -7: Length of private work array is not sufficient.
c          = -8: Error return from LAPACK eigenvalue calculation;
c          = -9: Starting vector is zero.
c          = -10: IPARAM(7) must be 1,2,3,4.
c          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatable.
c          = -12: IPARAM(1) must be equal to 0 or 1.
c          = -9999: Could not build an Arnoldi factorization.
c                   IPARAM(5) returns the size of the current Arnoldi
c                   factorization.

*/






