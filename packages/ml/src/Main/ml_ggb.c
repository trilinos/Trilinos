

/*************************************************************************************       
  HAIM: GLOBAL EIGENVALUE CALCULATIONS FOR MULTIGRID           
  THE GENERALIZED GLOBAL BASIS (GGB) METHOD USING ARPACK                      
***************************************************************************************/     
/*                                                            
   THE IDEA:                                                  
        We want to solve     (STS)*x = lam*x                      
        S - Smoothing iteration matrix
	T - Multilevel projector
		
        + We solve the problem by supplying ARPACK with the function that
	  produces Matrix-vector product 
	+ We use simple non-symmetric ARPACK mode                            
                                                              
	*/        


#ifdef HAVE_ML_PARPACK
#include <mpi.h>
#endif
                                             
       
#include "ml_ggb.h"



void ML_ARPACK_GGB( struct ML_Eigenvalue_Struct *eigen_struct,ML *ml,
		     struct ML_CSR_MSRdata *mydata, int Debug_Flag, 
		    int GGB_alp_flag) 
{

  int i ,  j, level;

  /* Eigenvalue definitions */
  int      iparam[11];
  int      nev, ncv, info, mode, nconv, Fattening;
  double   tol, tm, tmp_tol;
  char     bmat[2], which[3];
  ML_Operator *Amat;


  
  

  /* --------- Execution begins  ---------------- */
  
#if defined(HAVE_ML_ARPACK) || defined(HAVE_ML_PARPACK)
  
   
    if (ml->comm->ML_mypid==0) {
      printf("\n");
      ML_print_line("=", 80);
      printf("\t\tStarting Eigenvalue Calculation\n");
    }
    
  tm = GetClock();

  /* Some matrix information */
    
  level         = ml->ML_finest_level;
  Amat          = &(ml->Amat[level]);


  /* Set parameters for ARPACK from the input file */

  iparam[MAX_ITRS] = eigen_struct->Max_Iter;
  nev              = eigen_struct->Num_Eigenvalues;
  ncv              = eigen_struct->Arnoldi;
  tol              = eigen_struct->Residual_Tol;

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
		   nev, ncv, tol, ml, mydata,Fattening,
		   eigen_struct, Debug_Flag, GGB_alp_flag );
  
  if (ml->comm->ML_mypid==0)
    printf("Time for eigenvalue computations is %g (sec) on processor %d\n",GetClock()-tm,ml->comm->ML_mypid);
  



#else
  fprintf( stderr,
	   "ERROR: ML has not been configured with ARPACK support.\n"
	   "ERROR: Please reconfigure with the option `--with-ml_arpack'\n"
	   "ERROR: or `--with-ml_parpack'\n"
	   "ERROR: (file %s, line %d)\n",
	   __FILE__, __LINE__ );
  exit(-1);
#endif



  
}


/********************************************************************/
/********************************************************************/
/********************************************************************/



void  ML_ARPACK_driver(char which[],
		       char bmat[], int iparam[], int mode, 
		       int nev, int ncv, double tol,  ML *ml, 
		       struct ML_CSR_MSRdata *mydata, int Fattening,
		       struct ML_Eigenvalue_Struct *eigen_struct,
		       int Debug_Flag, int GGB_alp_flag)
{

  int        i, j, kk, ldv, lworkl;
  int        nloc,  nloc2, ido, flag, counter;
  int        count, nconv, ierr, info;
  int        ipntr[14], m=0;
  int        one = 1, dummy1, dummy2, dummy3, dummy4;
  double     a1 , a2  , lamR, lamI, time;
  char       string[4];
  double     *v, *workl, *workd, *workev, *d, *resid;         /* Pointers used by ARPACK */  
  int        *select, rvec, *work;
  int        proc_id;

  double     *vecx, *vecy, *rhs, *rhs1;                       /* Dummy Pointers */
  
  /* FILE       *ifp; */
  ML_Operator *Amat;

  int    comm;                                          /* MPI communicator */
  double tm;

  /********************************  Begin ************************************/

  /* Initialize some communications stuff for eig_driver */
  comm    = ml->comm->USR_comm; 
  proc_id = ml->comm->ML_mypid;

  /*  ifp   = fopen("CheckQmgs.m", "w");  */
  

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

  nloc2  =  nloc; 
  
  
  
  select = (int    *) ML_allocate(ncv*sizeof(int));
  work   = (int    *) ML_allocate(ncv*sizeof(int));
  rhs    = (double *) ML_allocate(nloc2*sizeof(double));
  rhs1   = (double *) ML_allocate(nloc2*sizeof(double));
  d      = (double *) ML_allocate(4*ncv*sizeof(double));
  resid  = (double *) ML_allocate(nloc2*sizeof(double));
  workd  = (double *) ML_allocate(3*nloc2*sizeof(double));
  workev = (double *) ML_allocate(3*ncv*sizeof(double));
  workl  = (double *) ML_allocate(lworkl*sizeof(double));
  v      = (double *) ML_allocate(ncv*ldv*sizeof(double));
  

  if (Debug_Flag > 2) {
    vecx   = (double *) ML_allocate(nloc2*sizeof(double));
    vecy   = (double *) ML_allocate(nloc2*sizeof(double));
  }
  
  if (v == NULL) {
    fprintf(stderr,"Not enough space to allocate workl\n");
    exit(-1);
  }
  
  
  /* Initialize dummy structures */
  for (j = 0; j < nloc2 ; j++ )  rhs[j] = rhs1[j] = resid[j] = 0.0;
  
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
    /*  tm = GetClock(); */
#ifdef HAVE_ML_PARPACK

    
    pdnaupd_(&comm, &ido, bmat, &nloc, which, &nev, &tol, resid,
	     &ncv, v, &ldv, iparam, ipntr, workd, workl,
	     &lworkl, &info );
    
#else    
#ifdef HAVE_ML_ARPACK
    
    dnaupd_(&ido, bmat, &nloc, which, &nev, &tol, resid,
	    &ncv, v, &ldv, iparam, ipntr, workd, workl,
	    &lworkl, &info );
#else
    fprintf(stderr, "ERROR with arpack/parpack\n");
    exit( EXIT_FAILURE );
#endif
#endif

    

    if ( (ido == -1) || (ido == 1) ) {
     
      count++;
      if (ml->comm->ML_mypid==0) printf("\t\t Eigensolver iteration: %d\n",count);
     
      /***********************************************
       * matrix vector multiplication (using inverse)
       *   workd[ipntr[1]-1] <-- OP * workd[ipntr[0]-1]
       ***********************************************/
     
      /************************************/
      /* Compute (I-inv(P)*K)v = lam * v  */
      /************************************/
      /*  time =  GetClock(); */
      /* First perform K*v = rhs */
      
      ML_Operator_Apply(Amat, Amat->invec_leng, &workd[ipntr[0]-1], 
			Amat->outvec_leng, rhs);


      /*******************************************************/
      /* Do ML cycle to find the iteration matrix eigenpairs */
      /* *****************************************************/
      

      /* Second perform inv(P)*rhs = rhs1 */
	ML_Solve_MGV(ml ,rhs , rhs1);
	
	/* printf("time for M^-1 * K * v = %g \n", GetClock()-time); */

	/* Third perform  workd[ipntr[1]-1] <- v-rhs1 */ 
	for (kk = 0 ; kk < nloc2 ; kk++) 
	  workd[ipntr[1]+kk-1] = workd[ipntr[0]+kk-1] - rhs1[kk];

	/*	
	  if (ml->comm->ML_mypid==0)
	  printf("Time for eigenvalue computations is %g (sec) on processor %d\n",GetClock()-tm,ml->comm->ML_mypid);
	*/
	
	
   }
      
      
    else flag = 0;
   
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
    
#if defined(HAVE_ML_ARPACK) || defined(HAVE_ML_PARPACK)
    /* C to Fortran wraper function */
    ml_pdneupc__(&comm,&rvec, string, select, d, v, &ldv,
		workev, bmat, &nloc, which, &nev,
		&tol, resid, &ncv, iparam, ipntr, workd, workl,
		&lworkl, &ierr, (ftnlen)1, (ftnlen)1, (ftnlen) 2);
#else
   fprintf( stderr, "ERROR with ARPACK or PARPACK\n");
   exit( EXIT_FAILURE );
#endif    
    
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
	      if (ml->comm->ML_mypid==0) {
		printf("\n");
		ML_print_line("=", 80);
		printf("\t\t GGB: We have complex eigenvectors \n");
		ML_print_line("=", 80);
	      }
	      counter = 1;
	      break;
	    }
	  }
	  
	  if (counter == 0)  {
	    if (ml->comm->ML_mypid==0) {
	      printf("\n");
	      ML_print_line("=", 80);
	      printf("\t\t GGB: We have real eigenvectors \n");
	      ML_print_line("=", 80);
	    }
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

	    /*
	    ML_Operator_Apply(Amat, Amat->invec_leng, &v[j*ldv],
			Amat->outvec_leng, rhs);
	    */
	    /* Need only to analyze the eigenvalues of the MG iteration matrix */
	    a1 = 0.0;
	    if (Debug_Flag > 2) {

	      ML_Operator_Apply(Amat, Amat->invec_leng, &v[j*ldv],
				Amat->outvec_leng, rhs);
	      
	      ML_Solve_MGV(ml ,rhs , rhs1);
	      
	      for (kk =0 ; kk < nloc2 ; kk++ ) 
		rhs1[kk]     = (1-lamR)*v[j*ldv+kk] - rhs1[kk];
	      

	      a1 = sqrt(ML_gdot(nloc, rhs1, rhs1 ,ml->comm)); 


	    }
	    
	    d[j+2*ncv] = a1; /*sqrt(lamR*lamR);*/
	   
	    
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
	    a2 = 0.0;
	    if (Debug_Flag > 2) {
	
	      /*
	      for (kk =0 ; kk < nloc2 ; kk++ ) {
		rhs[kk] = rhs1[kk] = vecx[kk] = vecy[kk] = 0.0;	      	     
	      }
	      */

	      /* The following mat vec we need always for MGGB */
	      /*
		ML_Operator_Apply(Amat, Amat->invec_leng, &v[j*ldv],
		Amat->outvec_leng, vecx);
	      */
	      
	      
	      /* Need only to analyze the eigenvalues of the MG iteration matrix */ 
	      a2 = 0.0;
	      
	      ML_Operator_Apply(Amat, Amat->invec_leng, &v[j*ldv],
				Amat->outvec_leng, vecx);
	      
	      ML_Operator_Apply(Amat, Amat->invec_leng, &v[(j+1)*ldv],
				Amat->outvec_leng, vecy);
	      
	      a1 = 1 - lamR;
	      
	      ML_Solve_MGV(ml ,vecx ,rhs);
	      ML_Solve_MGV(ml ,vecy ,rhs1);
	      
	      for (kk =0 ; kk < nloc ; kk++ ) {
		rhs[kk]  = a1*v[j*ldv+kk]     + lamI*v[(j+1)*ldv+kk] - rhs[kk];
		rhs1[kk] = a1*v[(j+1)*ldv+kk] - lamI*v[j*ldv+kk]     - rhs1[kk];
	      }
	      
	      
	      a2 = sqrt(ML_gdot(nloc2, rhs1, rhs1 ,ml->comm) +
			ML_gdot(nloc2, rhs, rhs ,ml->comm)); 
	      
	    }	    
	    
	    d[j+2*ncv] =   a2; 
	    
	    d[j+1+2*ncv] = a2;
	    
	    
	    d[j+1]       =  d[j];
	    d[j+1+ncv]   = -d[j+ncv];
	    j = j + 1;
	    
	    
	    /*-----------------------
	      | Ritz value is complex. |
	      | Residual of one Ritz   |
	      | value of the conjugate |
	      | pair is computed.      |
	      ------------------------*/
	  }
	    }
	    

      /*   Display computed residuals   */

	   dummy1 = 6; dummy2 = 3; dummy3 = ncv; dummy4 = -6;
	
	  /*	  cpdmout_(&comm, &dummy1, &nconv, &dummy2, d, &dummy3,&dummy4); */

     	  ml_pdmout__(&comm, &dummy1, &nconv, &dummy2, d, &dummy3, &dummy4);

      }


    /* Flag for the GGB alpha method */

    if (GGB_alp_flag == 1) {
      
      ML_GGBalp (v, nconv, nloc2, eigen_struct);
      
           
      ML_GGB2CSR (eigen_struct->Evec , eigen_struct->Pnconv, nloc2, proc_id, mydata, Debug_Flag);
      
      ML_free( v);
      ML_free( eigen_struct->Eval);



      eigen_struct->Nflag  = 2;  /* Flag to indicate the first Newton step */
      

    }
    else {
    /* Convert Information into CSR Format */

      ML_GGB2CSR (v, nconv, nloc2, proc_id, mydata, Debug_Flag);

      /* Keeping the eigenspace information to be reused with MGGB */
      eigen_struct->Evec   = v;
      eigen_struct->Eval   = d;
      eigen_struct->Pnconv = nconv; 
      eigen_struct->Nflag  = 1;
    }
    
    if (ml->comm->ML_mypid==0) {
#ifdef  HAVE_ML_PARPACK
      printf("\n\t\t Parallel arpack iterations\n"); 
#else
      printf("\n\t\t Serial arpack iterations\n");
#endif
    }

    /*  Print additional convergence information */
    if (ml->comm->ML_mypid==0) {
      
      if ( info == 1 ){
	printf("\nMaximum number of iterations reached.\n");
      }
      else if ( info == 3 ){
	printf("\nNo shifts could be applied during implicit\n");
	printf("Arnoldi update, try increasing NCV.\n\n");
      }
      else {
	 

	printf("\nEigenvalue Calculation Summary");
	printf("\n");    ML_print_line("=", 32);  
	printf("The number of processors is %d\n", ml->comm->ML_nprocs);
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
  }
 

  ML_free( select);
  ML_free( work);
  
  if (Debug_Flag > 2) {
    ML_free( vecx);
    ML_free( vecy);
  }
  ML_free( rhs1);
  ML_free( rhs);
  ML_free( resid);
  ML_free( workd);
  ML_free( workev);
  ML_free( workl);

    
  /*  fclose(ifp);   */


} /* end ML_ARPACK_driver */


/* Dense matrix-vector product */

extern int ML_OperatorGGB_Apply (double *densemat, int Nrows, int Ncols, double *din, double *dout, int Transpose)
{
  char     trans[2];
  double   alp=1.0, beta=0.0; 
  int      one =1, i;


  if (Transpose == 1)   trans[0] = 'T';
  if (Transpose == 0)   trans[0] = 'N';
 
  trans[1] = '\0';
  
  
   
  dgemv_(trans, &Nrows, &Ncols, &alp, densemat, &Nrows, din, &one, &beta, dout, &one, (ftnlen)1);


  return (1);
}






/******************************************************************************/
void ML_GGB2CSR (double *v, int nconv, int MatSize, int proc_id,
		   struct ML_CSR_MSRdata  *mydata, int Debug_Flag)

     /* Function to transfer dense columns of matrix to ML CSR format */
     /* The colums are input as one long vector */
     
{
  int          nrows,   ncolumns,  nnz;
  int          *rowptr, *columns;
  int          i, j , count;
  double       *values;
 
  FILE          *fp, *fp1, *eig;
  
  /******************************  begin *****************************/
  if (Debug_Flag == 10) {
    fp  = fopen("Rowptr.m","w");
    fp1 = fopen("Val_Col.m","w");
    eig = fopen("EIGvec.m","w");
  }

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

  
  for (i = 0; i < nrows; i++) {     
    
    rowptr[i+1] = rowptr[i] + nconv;

    
    for (j = 0; j < ncolumns; j++) {
    
      columns[count] = j;
      values[count]  = v[j*nrows+i];    
      count = count + 1;
    }
    
    
  }
  

  
  if (proc_id == 0) {
    fprintf(stdout,"\n\t *********************************************");
    fprintf(stdout,"\n\t       GGB PROLONGATION MATRIX (Processor 0)   ");
    fprintf(stdout,"\n\t      ---------------------------------------  ");
    fprintf(stdout,"\n\t            MATRIX SIZE:   %d * %d    ",nrows,ncolumns);
    fprintf(stdout,"\n\t            NON ZEROS        = %d     ",nnz);
    fprintf(stdout,"\n\t ***********************************************\n");  
  }
  
  

  /* DEBUGING */
  if (Debug_Flag == 10) {
    for (i = 0; i < nnz; i++)
      fprintf(fp1,"%20.13f    %d \n",values[i],columns[i]);
    
    for (i = 0; i < nrows+1; i++) 
      fprintf(fp,"%d \n",rowptr[i]);
    
    for (i = 0; i < ncolumns; i++) {
      fprintf(eig,"EIG NUM = %d\n",i+1);
      for (j = 0; j < nrows; j++) {
      	fprintf(eig,"\t %20.13f\n ",v[j*nrows+i]);
      }
    }
  }
 


  /*************************************************************************
     FEEDING THE NESECCARRY DATA FOR THE CSR SPARSE FORMAT INTO STRUCTURE
   **************************************************************************/
      
  mydata->Nrows   = nrows;
  mydata->Ncols   = ncolumns;
  mydata->Nnz     = nnz;
  mydata->rowptr  = rowptr;
  mydata->columns = columns;
  mydata->values  = values;


  
    
  if (Debug_Flag == 10) {
    fclose(fp);
    fclose(fp1);
    fclose(eig);
  }

} /* end function */



void  ML_GGBalp (double *NewVec, int nconv, int nloc2, struct ML_Eigenvalue_Struct
		   *eigen_struct)
{

  /**************************************************************/
/* Angle between two subspaces:                                 */
/* A  = Q1*R1                 A - GGB projection from step i    */ 
/* A1 = Q2*R2                 A1- GGB projection from step i+1  */  
/*                            QR factorization is based on      */
/*                            Househoulder reflectors           */
/*                                                              */
/* S = svd(Q1'*Q2)            SVD factorization of the          */ 
/*                            orthogonal basis                  */
/* cos(t) = min(S) = S(t,t)   Angle between the subspaces       */
/****************************************************************/

  double     *A , *current_vec;
  double     theta, eps= 5.0 ;

  int           m, n, i, j, k;
  int           nnew, nold, lwork;
  int           ind =1;
 


  /********************* begin ************************/
  nnew   =  nconv;
  ind    =  1;  

  nold   =  eigen_struct->Pnconv;
  m      =  nloc2;
  
  for (j=0; j<nnew ; j++) 
    {
      current_vec = (double *)  ML_allocate( ind*m* sizeof(double));
      k = 0;
      for (i=j*m; i<(j+1)*m; i++)  
	{
	  current_vec[k] = NewVec[i];
	  k = k+1;
	}


      /* Get the maximum principal angle and an orthogonal basis for the subspaces */
      theta = ML_subspace(m, eigen_struct->Evec, nold, current_vec, ind);
  
      theta = theta*57.2958;     /* Convert to degrees */
      printf("Angle between subspcaes is theta = %2.3f  \n",theta);
      
      
      /* Increase the space of the prolongation */
      if (theta > eps) {
	
	A = (double *)  ML_allocate( (nold+nnew)*m* sizeof(double)); 
	
	k = 0;
	for (i=0; i<(nold+ind)*m    ; i++) 
	  if (i<nold*m) A[i] = eigen_struct->Evec[i];
	  else          
	    {
	      A[i] =  current_vec[k];
	      k = k + 1;
	    }
	
	ML_free( eigen_struct->Evec);
	ML_free( current_vec);


	
	eigen_struct->Pnconv = nold+ind;
	eigen_struct->Evec   = A;    
      }
    }




}   


extern double  ML_subspace (int nrows, double *inp1, int ncols1, double *inp2, int ncols2)
{ 
  double     *tau, *work, *A, *S ,*U, *VT;
  double     *tau1, *work1, *A1, *B;
  double     theta;

  int         lda, lwork, info, ldv, ldu, ldvt;
  int         lwork1, info1, ldv1;


  int           m, n, i, j, k, one=1;
  int           nnew;
  char          jobu[2], jobvt[2];

  FILE          *fp2, *fp1, *fp;

  /*******************  begin  ***********************/

  if (ncols2 > ncols1) {
    printf("First matrix is assumed to be larger than the second. Change inputs and try again \n");
    exit(-1);
  } 
  
  /*
  fp2  = fopen("Amat.m","w");
  k = 0;
  for (i=0; i<ncols1; i++) 
    for (j=0; j<nrows; j++) {   
      fprintf(fp2,"A(%d,%d)=%20.13f;  \n",j+1,i+1,inp1[k]);
      k = k + 1;
    }    
  k = 0;
  for (i=0; i<ncols2; i++) 
    for (j=0; j<nrows; j++) {   
      fprintf(fp2,"A1(%d,%d)=%20.13f;  \n",j+1,i+1,inp2[k]);
      k = k + 1;
    }
  fclose(fp2);

  */

  m                     =  nrows; 
  lwork                 =  10*ncols1;
  lda                   =  m;
  lwork1                =  10*ncols2;


  tau    = (double *)  ML_allocate(ncols1* sizeof(double));
  tau1   = (double *)  ML_allocate(ncols2* sizeof(double));
  work   = (double *)  ML_allocate(lwork* sizeof(double));
  work1  = (double *)  ML_allocate(lwork1* sizeof(double));
  A1     = (double *)  ML_allocate(m*ncols2* sizeof(double));
  A      = (double *)  ML_allocate(m*ncols1* sizeof(double));


  /* Previous eigenvectors */
  for (i=0; i<m*ncols1; i++) A[i] = inp1[i];
  /* Current eigenvectors */
  for (i=0; i<m*ncols2; i++) A1[i] = inp2[i];


  /* Get QR factorization for both matrices */
  dgeqrf_(&m, &ncols1, A, &lda, tau, work, &lwork, &info);  
  dgeqrf_(&m, &ncols2, A1, &lda, tau1, work1, &lwork1, &info1);  
  

  if (info !=0 | info1 !=0)  {
    printf("Problem with QR factorization in ML_subspace function dgeqrf_\n");
    exit(-1);
  }
  

  /* Extract the R matrix from the output */
  /*
  fp1  = fopen("Rmat.m","w");
  for (i=0; i<ncols1; i++) 
    for (j=0; j<i+1; j++) 
      fprintf(fp1,"R(%d,%d)=%20.13f;  \n",j+1,i+1,A[m*i+j]);
  
  for (i=0; i<ncols2; i++) 
    for (j=0; j<i+1; j++) 
  fprintf(fp1,"R1(%d,%d)=%20.13f;  \n",j+1,i+1,A1[m*i+j]);
  
  fclose(fp1);
  */
   
  

  /* Now we obtain Q's  */
  dorgqr_(&m, &ncols1,&ncols1, A, &lda, tau, work, &lwork, &info);
  dorgqr_(&m, &ncols2, &ncols2, A1, &lda, tau1, work1, &lwork1, &info1);

  if (info !=0 | info1 !=0)  {
    printf("Problem with QR factorization in ML_subspace function dorgqr_\n");
    exit(-1);
  }
  /* Print the Q matrix of the first matrix */
  /*
  fp  = fopen("Qmat.m","w");
  k = 0;
  for (i=0; i<ncols1; i++) 
    for (j=0; j<m; j++) {   
      fprintf(fp,"Q(%d,%d)=%20.13f;  \n",j+1,i+1,A[k]);      
      k = k + 1;
    }
  k = 0;
  for (i=0; i<ncols2; i++) 
    for (j=0; j<m; j++) {   
      fprintf(fp,"Q1(%d,%d)=%20.13f;  \n",j+1,i+1,A1[k]);      
      k = k + 1;
    }
  fclose(fp);    
  */

  
  ML_free( tau);
  ML_free( work);
  ML_free( tau1);
  ML_free( work1);
  

  /* Perform  B = A1'*A, to compute the principal angle with SVD */
  B        = (double *)  ML_allocate(ncols1*ncols2* sizeof(double));
  
  k = 0;
  for (i=0; i<ncols2; i++) {
    for (j=0; j<ncols1; j++) {  
      B[k] =  ddot_(&m, &A[j*m], &one, &A1[i*m], &one);   
      k = k + 1;      
    }
  } 
  

  /* Compute SVD - we are interested only in the singular values */
  jobu[0]  = 'N'; 
  jobu[1]  = '\0';
  jobvt[0] = 'N';
  jobvt[1] = '\0';
  
  lda      =  ncols1;      
  lwork    =  10*ncols1;   
  ldu      =  ncols1;
  ldvt     =  ncols2;

  S        = (double *)  ML_allocate(ncols2* sizeof(double));
  work     = (double *)  ML_allocate(lwork* sizeof(double));

  /* SVD */
  dgesvd_(jobu, jobvt, &ncols1, &ncols2, B, &lda, S, U, &ldu, VT, &ldvt, 
	  work, &lwork, &info);

  if (info !=0 )  {
    printf("Problem with QR factorization in ML_subspace function dgesvd_\n");
    exit(-1);
  }

  /*  
  for (i=0; i<ncols2; i++) 
    printf("S[%d] = %3.30f , theta=%2.3f \n",i,S[i],acos(S[i])*180/3.14);
  */


  /* theta = acos(min(S)) */
  if (S[ncols2-1] > 1.0) theta = 0.0;
  else  theta = acos(S[ncols2-1]);

 


  
  ML_free( A);
  ML_free( A1);  
  ML_free( S);
  ML_free( B);
  ML_free( work);


  return theta;

}


/*********************************************************
   MGGB Variant - adaptive eigenspace computations
   based on the angle between vector and subspace (Q,R*q)
***********************************************************/
 
int  ML_MGGB_angle( struct ML_Eigenvalue_Struct *eigen_struct,ML *ml,
		     struct ML_CSR_MSRdata *mydata)
{
  int            Nrows, Ncols, i, level;  
  int            kk, ggb_flag = 0, count;
  double        *rhs, *rhs1;
  double         epsilon= 30.0;
  double         tm, theta, pi = 3.1415;


  double        *Rq, *vec, *dumm , *val;


  ML_Operator   *Amat;

  /*************************** begin *****************************/
  tm = GetClock();


  Ncols = mydata->Ncols;
  Nrows = mydata->Nrows;


  level         = ml->ML_finest_level;
  Amat          = &(ml->Amat[level]);

  /****************************************************************/
  
  vec    = eigen_struct->Evec;    /* information with all eigenvectors */
  val    = eigen_struct->Eval;    /* information with all eigenvalues */

  /* How many new vectors should we use to compute the angles  */
  if (Ncols >= 2) count = 2;
  else count = 1;
  
#ifdef HAVE_ML_PARPACK  
  
  /* see if the first eigenvector is imaginary or not  */
  if (count == 2 &&  val[eigen_struct->Arnoldi] != 0.0) count = 2;
  else count = 1;
  

  ggb_flag  = ML_Rayleigh(ml, Nrows, vec, count);  /* get flag based on Rayliegh quotient */
  
#else  /* get flag  based on angle between a subspace and a vector */
  
  /* How many new vectors should we compute the angle for */ 
  if (Ncols >= 2) count = 2;
  else count = 1;
  
  /* Rq will hold the approximated vectors: R*q */
  /* where R is the multigrid iteration operator */
  
  Rq      = (double *) ML_allocate(Nrows*count* sizeof(double));
  rhs    = (double *) ML_allocate(Nrows*sizeof(double));
  rhs1   = (double *) ML_allocate(Nrows*sizeof(double));
  dumm   = (double *) ML_allocate(Nrows*sizeof(double));
  

  for (i=0; i<count; i++) {
    
    for (kk =0 ; kk < Nrows ; kk++ ) 
      dumm[kk] = vec[i*Nrows+kk];
    
    
    ML_Operator_Apply(Amat, Amat->invec_leng, dumm, Amat->outvec_leng, rhs);
    ML_Solve_MGV(ml ,rhs , rhs1);
    
    for (kk =0 ; kk < Nrows ; kk++ ) 
      Rq[i*Nrows+kk] = dumm[kk]  - rhs1[kk];
    
  }
  

  theta = ML_subspace(Nrows, vec, Ncols, Rq, count);
  
  theta = theta*57.2958; 
  if (theta > epsilon) ggb_flag = 1;
  
  printf("\n");
  ML_print_line("=", 80);
  printf("Angle between subspcaes is theta = %2.3f\n",theta);
   
  ML_free( Rq);
  ML_free( rhs);
  ML_free( rhs1);
  ML_free( dumm);
  
#endif

  if (ml->comm->ML_mypid==0) {
    printf("Time for MGGB eigenspace angle measure is %g (sec.)\n",GetClock()-tm);
    if (ggb_flag == 1) printf("Recomputing eigenspace \n");
    else printf("Reusing previous eigenspace information \n");
    ML_print_line("=", 80);
    printf("\n");
  }
  
  return ggb_flag;

}

/*********************************************************************************
   Computes the norm of approximated (highest) eigenvector based on 
   Rayleigh quotient. Based on this returns flag whether to recompute new 
   eigenspace or not   

   || Rq - lam*q || < eps

    R    = I - (M^-1)A
    q    - computed eigenvector (q = qx + i*qy)
            q* R q 
    lam  = --------  is the rayliegh quotient
             q* q

**********************************************************************************/

extern int  ML_Rayleigh (ML *ml, int nrows, double *q, int count)
{ 
  int      i , j, level, flag=1;
  double   rayl,  theta = 0.0;
  double    *u, *v, *rhs, *rhs1;
  double   *compA, *compB, norm, eps ;
  ML_Operator  *Amat;
  FILE          *fp;

  /*******************  begin  ***********************/

  /* fp  = fopen("Comp_Norm.m","w");*/

  if (count == 2) {  /* imaginary eigenvector */

    /* allocate memory for working arrays */
    u      = (double *) ML_allocate(nrows*sizeof(double));
    v      = (double *) ML_allocate(nrows*sizeof(double));
    rhs    = (double *) ML_allocate(nrows*sizeof(double));
    rhs1   = (double *) ML_allocate(nrows*sizeof(double));
    
    level         = ml->ML_finest_level;
    Amat          = &(ml->Amat[level]);
    
    /* Get the product   u<-A*qx  and  v<-A*qy  */ 
    ML_Operator_Apply(Amat, Amat->invec_leng, &q[0],     Amat->outvec_leng, u);
    ML_Operator_Apply(Amat, Amat->invec_leng, &q[nrows], Amat->outvec_leng, v);
   
    /* Get the product rhs <- (M^-1)*u   and  rhs1 <- (M^-1)*v  */ 
    ML_Solve_MGV(ml ,u ,rhs);
    ML_Solve_MGV(ml ,v ,rhs1);
    
    /* Get the final rhs <- R*qx   and rhs1 <- R*qy          */
    for (i =0 ; i < nrows ; i++ ) {
      rhs[i]  = q[i] - rhs[i]; 
      rhs1[i] = q[i+nrows] - rhs1[i];
    }
    
    
    /* Get complex vector product */ 
    compA  =  ML_complex_gdot(nrows, &q[0], &q[nrows], rhs, rhs1, ml->comm);
    compB  =  ML_complex_gdot(nrows, &q[0], &q[nrows], &q[0], &q[nrows], ml->comm);
    
    /* since compB is real we can scale compA */
    compA[0] = compA[0]/compB[0];
    compA[1] = compA[1]/compB[0];
    
    /* Get the norm || R * q - lam * q ||_2  and put the real and imag parts in u and v*/
    for (i =0 ; i < nrows ; i++ ) {
      u[i]  =   rhs[i]  - (compA[0]*q[i] - compA[1]*q[i+nrows]);   /* real part */  
      v[i]  =   rhs1[i] - (compA[0]*q[i+nrows] + compA[1]*q[i]);   /* imag part */
    }
    compB  =  ML_complex_gdot(nrows, u, v, u, v, ml->comm);
    norm   =  sqrt(compB[0]);
    
    eps = sqrt(compA[0]*compA[0]+compA[1]*compA[1]);  /* absolute value of Rayleigh quotient */
    
    
    /* free complex arrays */
    ML_free(v);
    ML_free(rhs1);
    ML_free(compA);
    ML_free(compB);



  }
  else {    /* real eigenvectors */
        
    /* allocate memory for working arrays */
    u      = (double *) ML_allocate(nrows*sizeof(double));
    rhs    = (double *) ML_allocate(nrows*sizeof(double));
    
    level         = ml->ML_finest_level;
    Amat          = &(ml->Amat[level]);
    
   
    /* Get the product   u<-A*qx  */ 
    ML_Operator_Apply(Amat, Amat->invec_leng, &q[0], Amat->outvec_leng, u);

    /* Get the product rhs <- (M^-1)*u  */ 
    ML_Solve_MGV(ml ,u ,rhs);
    
    /* Get the final rhs <- R*qx    */
    for (i =0 ; i < nrows ; i++ )  rhs[i]  = q[i] - rhs[i]; 
          
    /* Get Rayliegh Quotient */            
    rayl  =  ML_gdot(nrows, &q[0], rhs, ml->comm);
   

    /* Get the norm || R * q - lam * q ||_2  and put the real and imag parts in u and v*/
    
    for (i =0 ; i < nrows ; i++ )  u[i]  =   rhs[i]  - rayl*q[i];   /* real part */  

    
    norm      =  sqrt(ML_gdot(nrows, u, u, ml->comm));
    eps       =  abs(rayl);

  }
    
    /* We use the following formula to determine the flag return 
       
        || R*q - lam*q ||
      -------------------- = tan(theta)    
           | lam |   
	   
	   theta is the angle between (R*q,q). For computations we set a critical angle to be 30 deg.
	   i.e., eps = tan(30)*|lam| = 0.5774 * |lam|
	   eps = tan(50)*|lam| = 1.1918 * |lam|
    */
    

    if (norm <= eps*1.1918) flag = 0;
    else flag =1;
    
    if (ml->comm->ML_mypid==0) {
      printf("\n");
      ML_print_line("=", 80);
      printf("Angle based on Rayliegh Quotient is %2.0f (deg.)\n", atan(norm/eps)*180/3.1415);
    }

     
    ML_free(u);
    ML_free(rhs);


    return flag; 
}

/***********************************************
     Dot product between to comlex vectors  
************************************************/
extern double *ML_complex_gdot(int leng, double *ureal, double *uimag, double *vreal, double *vimag, 
			      ML_Comm *comm) {

  /* input:
     leng   -  the length of the vectors
     ureal  -  pointer to array of the first vector real part 
     uimag  -  pointer to array of the first vector imaginary part
     vreal  -  pointer to array of the second vector real part 
     vimag  -  pointer to array of the second vector imaginary part
     comm   -  ml communication structure

     output:
     ptr[0] = real( u'*v);
     ptr[1] = imag( u'*v);
     pointer containing the real and imaginary part of the dot product 
	  
      
  */
  
  double    *ptr;
  /*--------------------------------------------------------------------------------------------*/
  
  ptr    = (double *) ML_allocate(2*sizeof(double));
  
  /* real part */
  ptr[0] = ML_gdot(leng, ureal, vreal ,comm) + ML_gdot(leng, uimag, vimag ,comm); 
  ptr[1] = ML_gdot(leng, ureal, vimag ,comm) - ML_gdot(leng, uimag, vreal ,comm);


  return(ptr);

}




extern double  ML_normc (double *real, double *imag, int leng)
{

  double    rl, im, dum1, dum2;
  double    norm;
  int       i;

  dum1 = dum2 = 0.0;

  for (i=0; i<leng; i++) {

    dum1 = dum1 + real[i]*real[i];
    dum2 = dum2 + imag[i]*imag[i];

  }

  rl = dum1 + dum2;  /* Real part */

  norm = sqrt(rl); 
  
  return norm;
}

/************************/
/* Free eigen structure */
/************************/
extern void ML_Eig_Destroy(void *data) 
{
  
  struct  ML_Eigenvalue_Struct *temp;

   temp = (struct ML_Eigenvalue_Struct *) data;
   if (temp != NULL) {
     if (temp->Evec != NULL) ML_free(temp->Evec);
     if (temp->Eval != NULL) ML_free(temp->Eval);
   
      ML_free(temp);
   }




}


/*********************************************/
/*  ARPACK Documentation   if info < 0       */ 
/*********************************************/

/* taken from dnaupd.f file 

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
