/******************************************************************************
* source for function
*
  double ML_Cycle_MG(ML_1Level *curr, double *sol, double *rhs,
                     int approx_all_zeros, ML_Comm *comm,
                     int res_norm_or_not, ML *ml)
*--------------------------------------------------------------------------- */

   int         i, lengc, lengf;
   double      *res,  *sol2, *rhs2, res_norm = 0., *normalscales;
   double      *rhss, *dtmp;
   ML_Operator *Amat, *Rmat;
   ML_Smoother *pre,  *post;
   ML_CSolve   *csolve;
#ifdef RAP_CHECK
   double    norm1, norm2;
#endif

#ifdef ML_ANALYSIS
   short   dummy;
   int     *cols, Nrows, j, allocated_space, ncols, lwork, info;
   double  *squareA, *vals, *eig, *work;
   char    instring[100], jobz, jobz2;
   FILE    *fp;
#endif

   Amat     = curr->Amat;
   Rmat     = curr->Rmat;
   pre      = curr->pre_smoother;
   post     = curr->post_smoother;
   csolve   = curr->csolve;
   lengf    = Amat->outvec_leng;

   /* ------------------------------------------------------------ */
   /* first do the normalization                                   */
   /* ------------------------------------------------------------ */

   rhss = (double *) ML_allocate( lengf * sizeof(double) );
   ML_DVector_GetDataPtr(curr->Amat_Normalization, &normalscales) ;
   for ( i = 0; i < lengf; i++ ) rhss[i] = rhs[i];

#ifdef ML_ANALYSIS
   if ( comm->ML_nprocs == 1 && curr->Pmat->to == NULL && lengf < 1000 )
   {
      fp = fopen("mlmatlab.m", "w");
      Nrows = lengf;
      fprintf(fp, "A = sparse(%d,%d);\n", Nrows, Nrows);
      allocated_space = 100;
      cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
      vals = (double *) ML_allocate(allocated_space*sizeof(double));
      for (i = 0; i < lengf; i++) {
         while(ML_Operator_Getrow(Amat,1,&i,allocated_space,cols,vals,&ncols)==0)
         {
            allocated_space = 2*allocated_space + 1;
            ML_free(vals); ML_free(cols);
            cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
            vals = (double *) ML_allocate(allocated_space*sizeof(double));
            if (vals == NULL) {
               printf("Not enough space to get matrix row. Row length of\n");
               printf("%d was not sufficient\n",(allocated_space-1)/2);
               exit(1);
            }
         }
         for (j = 0; j < ncols; j++)
            fprintf(fp, "A(%d,%d)=%25.16e;\n",i+1,cols[j]+1,vals[j]);
      }
      fprintf(fp, "[eigv,eig]=eig(full(A));\n");
      fprintf(fp, "rhs=zeros(%d,1);\n", Nrows);
      for (j = 0; j < Nrows; j++)
         fprintf(fp, "rhs(%d)=%25.16e;\n",j+1,rhss[j]);
      fprintf(fp, "res=eigv'*rhs;\n");
      fprintf(fp, "plot(res)\n");
      fclose(fp);
      printf("** BEFORE pre-smoothing -- \n");
      printf("** Now you can use matlab to call the file mlmatlab.m\n");
      printf("   and it should display the residual for each eigenvalue.\n");
      printf("Press y when you are done.");
      scanf("%s", instring);
   }
#endif

   /* ------------------------------------------------------------ */
   /* smoothing or coarse solve                                    */
   /* ------------------------------------------------------------ */
   if (Rmat->to == NULL) {    /* coarsest grid */
      if ( ML_CSolve_Check( csolve ) == 1 ) {
         ML_CSolve_Apply(csolve, lengf, sol, lengf, rhss);
      } else {
         ML_Smoother_Apply(pre, lengf, sol, lengf, rhss, approx_all_zeros);
         ML_Smoother_Apply(post, lengf, sol, lengf, rhss, ML_NONZERO);
      }
      if (res_norm_or_not == ML_COMPUTE_RES_NORM) {
         res = (double *) ML_allocate(lengf*sizeof(double));
         ML_Operator_Apply(Amat, lengf, sol, lengf, res);
         for ( i = 0; i < lengf; i++ ) res[i] = rhss[i] - res[i];
         res_norm = sqrt(ML_gdot(lengf, res, res, comm));
         ML_free(res);
      }
   }
   else {
      res = (double *) ML_allocate(lengf*sizeof(double));

      /* --------------------------------------------------------- */
      /* pre-smoothing and compute residual                        */
      /* --------------------------------------------------------- */
      ML_Smoother_Apply(pre, lengf, sol, lengf, rhss, approx_all_zeros);

      if ( ( approx_all_zeros != ML_ZERO ) ||
           ( pre->smoother->ML_id != ML_EMPTY ) )
      {
         ML_Operator_Apply(Amat, lengf, sol, lengf, res);
         for ( i = 0; i < lengf; i++ ) res[i] = rhss[i] - res[i];
      }
      else for ( i = 0; i < lengf; i++ ) res[i] = rhss[i];

#ifdef ML_ANALYSIS
      if ( comm->ML_nprocs == 1 && curr->Pmat->to == NULL && lengf < 1000 )
      {
         fp = fopen("mlmatlab.m", "w");
         Nrows = lengf;
         fprintf(fp, "A = sparse(%d,%d);\n", Nrows, Nrows);
         allocated_space = 100;
         cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
         vals = (double *) ML_allocate(allocated_space*sizeof(double));
         for (i = 0; i < lengf; i++) {
            while(ML_Operator_Getrow(Amat,1,&i,allocated_space,cols,vals,&ncols)
                  == 0)
            {
               allocated_space = 2*allocated_space + 1;
               ML_free(vals); ML_free(cols);
               cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
               vals = (double *) ML_allocate(allocated_space*sizeof(double));
               if (vals == NULL) {
                  printf("Not enough space to get matrix row. Row length of\n");
                  printf("%d was not sufficient\n",(allocated_space-1)/2);
                  exit(1);
               }
            }
            for (j = 0; j < ncols; j++)
               fprintf(fp, "A(%d,%d)=%25.16e;\n",i+1,cols[j]+1,vals[j]);
         }
         fprintf(fp, "[eigv,eig]=eig(full(A));\n");
         fprintf(fp, "rhs=zeros(%d,1);\n", Nrows);
         for (j = 0; j < Nrows; j++)
            fprintf(fp, "rhs(%d)=%25.16e;\n",j+1,res[j]);
         fprintf(fp, "res=eigv'*rhs;\n");
         fprintf(fp, "plot(res)\n");
         fclose(fp);
         printf("** AFTER  pre-smoothing -- \n");
         printf("** Now you can use matlab to call the file mlmatlab.m\n");
         printf("   and it should display the residual for each eigenvalue.\n");
         printf("Press y when you are done.");
         scanf("%s", instring);
     }
/*
         printf("Constructing eigensystem...\n");
         Nrows = lengf;
         squareA = ML_allocate( Nrows * Nrows * sizeof(double) );
         for ( i = 0; i < Nrows*Nrows; i++ ) squareA[i] = 0.0;
         allocated_space = 100;
         cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
         vals = (double *) ML_allocate(allocated_space*sizeof(double));
         for (i = 0; i < lengf; i++) {
            while(ML_Operator_Getrow(Amat,1,&i,allocated_space,cols,vals,&ncols)
                  == 0)
            {
               allocated_space = 2*allocated_space + 1;
               ML_free(vals); ML_free(cols);
               cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
               vals = (double *) ML_allocate(allocated_space*sizeof(double));
               if (vals == NULL) {
                  printf("Not enough space to get matrix row. Row length of\n");
                  printf("%d was not sufficient\n",(allocated_space-1)/2);
                  exit(1);
               }
            }
            for (j = 0; j < ncols; j++)
               squareA[cols[j]*Nrows+i] = vals[j];
         }
         ML_free(cols); ML_free(vals);
         eig = (double *) ML_allocate( Nrows * sizeof(double) );
         work = (double *) ML_allocate( 4 * Nrows * sizeof(double) );
         lwork = 4 * Nrows;
         jobz = 'V'; jobz2 = 'U';
         MLFORTRAN(dsyev)(&jobz,&jobz2,&Nrows,squareA,&Nrows,eig,work,&lwork,&info,
                dummy,dummy);
         printf("returning from dsyev ...\n");
         if ( info > 0 )
         {
            printf("No convergence in computing the eigenvalues.\n");
         } else if ( info < 0 )
         {
            printf("Eigenvalue computation: %d-th argument has error.\n",-info);
         } else {
            fp = fopen("mlmatlab.eig", "w");
            fprintf(fp, "%d\n", Nrows);
            for (i=0; i<Nrows; i++) fprintf(fp, "%25.16e\n",eig[i]);
            for (i=0; i<Nrows*Nrows; i++)
               fprintf(fp, "%25.16e\n",squareA[i]);
            fclose(fp);
         }
         ML_free(squareA);
         ML_free(eig);
         ML_free(work);
         fp = fopen("mlmatlab.m", "w");
         fprintf(fp, "res = [\n");
         for (i=0; i<Nrows; i++)
         {
            work[i] = 0.0;
            for (j=0; j<Nrows; j++)
               work[i] += ( squareA[i*Nrows+j] * rhss[j]);
            fprintf(fp, "    %25.16e \n", work[i]);
         }
         fclose(fp);
         printf("** BEFORE pre-smoothing -- \n");
         printf("** Now you can use matlab to call the file mlmatlab.m\n");
         printf("   and it should display the residual for each eigenvalue.\n");
         printf("Press y when you are done.");
         scanf("%s", instring);
*/
#endif

      if (res_norm_or_not == ML_COMPUTE_RES_NORM)
         res_norm = sqrt(ML_gdot(lengf, res, res, comm));

      lengc = Rmat->outvec_leng;

      rhs2 = (double *) ML_allocate(lengc*sizeof(double));
      sol2 = (double *) ML_allocate(lengc*sizeof(double));
      for ( i = 0; i < lengc; i++ ) sol2[i] = 0.0;

      /* --------------------------------------------------------- */
      /* normalization                                             */
      /* --------------------------------------------------------- */
      ML_DVector_GetDataPtr(curr->Amat_Normalization, &normalscales) ;
      if ( normalscales != NULL )
         for ( i = 0; i < lengf; i++ ) res[i] /= normalscales[i];

      /* ------------------------------------------------------------ */
      /* transform the data from equation to grid space, do grid      */
      /* transfer and then transfer back to equation space            */
      /* ------------------------------------------------------------ */
      if ( ML_Mapper_Check(curr->eqn2grid) == 1 )
      {
         dtmp = (double *) ML_allocate( lengf * sizeof( double ) );
         ML_Mapper_Apply(curr->eqn2grid, res, dtmp );
         for ( i = 0; i < lengf; i++ ) res[i] = dtmp[i];
         ML_free( dtmp );
      }
      ML_Operator_ApplyAndResetBdryPts(Rmat, lengf, res, lengc, rhs2);
      if ( ML_Mapper_Check(Rmat->to->grid2eqn) == 1 )
      {
         dtmp = (double *) ML_allocate( lengc * sizeof( double ) );
         ML_Mapper_Apply(Rmat->to->grid2eqn, rhs2, dtmp );
         for ( i = 0; i < lengc; i++ ) rhs2[i] = dtmp[i];
         ML_free( dtmp );
      }
      ML_DVector_GetDataPtr(Rmat->to->Amat_Normalization,&normalscales);
      if ( normalscales != NULL )
         for ( i = 0; i < lengc; i++ ) rhs2[i] = rhs2[i] * normalscales[i];

      /* --------------------------------------------------------- */
      /* process the next level and transfer back to this level    */
      /* --------------------------------------------------------- */
#ifdef ML_SINGLE_LEVEL_PROFILING
      i = curr->levelnum-1;
      switch(i) {
         case 9:
      ML_Cycle_MG_9( Rmat->to, sol2, rhs2, ML_ZERO,comm, ML_NO_RES_NORM, ml);
      if ( (ml->ML_scheme == ML_MGW) && (Rmat->to->Rmat->to != NULL))
	ML_Cycle_MG_9( Rmat->to, sol2, rhs2, ML_NONZERO,comm, ML_NO_RES_NORM,ml);
           break;
         case 8:
      ML_Cycle_MG_8( Rmat->to, sol2, rhs2, ML_ZERO,comm, ML_NO_RES_NORM, ml);
      if ( (ml->ML_scheme == ML_MGW) && (Rmat->to->Rmat->to != NULL))
	ML_Cycle_MG_8( Rmat->to, sol2, rhs2, ML_NONZERO,comm, ML_NO_RES_NORM,ml);
           break;
         case 7:
      ML_Cycle_MG_7( Rmat->to, sol2, rhs2, ML_ZERO,comm, ML_NO_RES_NORM, ml);
      if ( (ml->ML_scheme == ML_MGW) && (Rmat->to->Rmat->to != NULL))
	ML_Cycle_MG_7( Rmat->to, sol2, rhs2, ML_NONZERO,comm, ML_NO_RES_NORM,ml);
           break;
         case 6:
      ML_Cycle_MG_6( Rmat->to, sol2, rhs2, ML_ZERO,comm, ML_NO_RES_NORM, ml);
      if ( (ml->ML_scheme == ML_MGW) && (Rmat->to->Rmat->to != NULL))
	ML_Cycle_MG_6( Rmat->to, sol2, rhs2, ML_NONZERO,comm, ML_NO_RES_NORM,ml);
           break;
         case 5:
      ML_Cycle_MG_5( Rmat->to, sol2, rhs2, ML_ZERO,comm, ML_NO_RES_NORM, ml);
      if ( (ml->ML_scheme == ML_MGW) && (Rmat->to->Rmat->to != NULL))
	ML_Cycle_MG_5( Rmat->to, sol2, rhs2, ML_NONZERO,comm, ML_NO_RES_NORM,ml);
           break;
         case 4:
      ML_Cycle_MG_4( Rmat->to, sol2, rhs2, ML_ZERO,comm, ML_NO_RES_NORM, ml);
      if ( (ml->ML_scheme == ML_MGW) && (Rmat->to->Rmat->to != NULL))
	ML_Cycle_MG_4( Rmat->to, sol2, rhs2, ML_NONZERO,comm, ML_NO_RES_NORM,ml);
           break;
         case 3:
      ML_Cycle_MG_3( Rmat->to, sol2, rhs2, ML_ZERO,comm, ML_NO_RES_NORM, ml);
      if ( (ml->ML_scheme == ML_MGW) && (Rmat->to->Rmat->to != NULL))
	ML_Cycle_MG_3( Rmat->to, sol2, rhs2, ML_NONZERO,comm, ML_NO_RES_NORM,ml);
           break;
         case 2:
      ML_Cycle_MG_2( Rmat->to, sol2, rhs2, ML_ZERO,comm, ML_NO_RES_NORM, ml);
      if ( (ml->ML_scheme == ML_MGW) && (Rmat->to->Rmat->to != NULL))
	ML_Cycle_MG_2( Rmat->to, sol2, rhs2, ML_NONZERO,comm, ML_NO_RES_NORM,ml);
           break;
         case 1:
      ML_Cycle_MG_1( Rmat->to, sol2, rhs2, ML_ZERO,comm, ML_NO_RES_NORM, ml);
      if ( (ml->ML_scheme == ML_MGW) && (Rmat->to->Rmat->to != NULL))
	ML_Cycle_MG_1( Rmat->to, sol2, rhs2, ML_NONZERO,comm, ML_NO_RES_NORM,ml);
           break;
         case 0:
      ML_Cycle_MG_0( Rmat->to, sol2, rhs2, ML_ZERO,comm, ML_NO_RES_NORM, ml);
      if ( (ml->ML_scheme == ML_MGW) && (Rmat->to->Rmat->to != NULL))
	ML_Cycle_MG_0( Rmat->to, sol2, rhs2, ML_NONZERO,comm, ML_NO_RES_NORM,ml);
           break;
         default:
      ML_Cycle_MG( Rmat->to, sol2, rhs2, ML_ZERO,comm, ML_NO_RES_NORM, ml);
      if ( (ml->ML_scheme == ML_MGW) && (Rmat->to->Rmat->to != NULL))
	ML_Cycle_MG( Rmat->to, sol2, rhs2, ML_NONZERO,comm, ML_NO_RES_NORM,ml);
           break;
   }

#else
      ML_Cycle_MG( Rmat->to, sol2, rhs2, ML_ZERO,comm, ML_NO_RES_NORM, ml);
      if ( (ml->ML_scheme == ML_MGW) && (Rmat->to->Rmat->to != NULL))
	ML_Cycle_MG( Rmat->to, sol2, rhs2, ML_NONZERO,comm, ML_NO_RES_NORM,ml);
#endif

      /* ------------------------------------------------------------ */
      /* transform the data from equation to grid space, do grid      */
      /* transfer and then transfer back to equation space            */
      /* ------------------------------------------------------------ */
      if ( ML_Mapper_Check(Rmat->to->eqn2grid) == 1 )
      {
         dtmp = (double *) ML_allocate( lengc * sizeof( double ) );
         ML_Mapper_Apply(Rmat->to->eqn2grid, sol2, dtmp);
         for ( i = 0; i < lengc; i++ ) sol2[i] = dtmp[i];
         ML_free( dtmp );
      }
      ML_Operator_ApplyAndResetBdryPts(Rmat->to->Pmat,lengc,sol2,lengf,res);
      if ( ML_Mapper_Check(curr->grid2eqn) == 1 )
      {
         dtmp = (double *) ML_allocate( lengf * sizeof( double ) );
         ML_Mapper_Apply(curr->grid2eqn, res, dtmp);
         for ( i = 0; i < lengf; i++ ) res[i] = dtmp[i];
         ML_free( dtmp );
      }

      /* --------------------------------------------------------- */
      /* post-smoothing                                            */
      /* --------------------------------------------------------- */
      for ( i = 0; i < lengf; i++ ) sol[i] += res[i];
#if defined(RAP_CHECK) || defined(ANALYSIS)

   /* When using RAP, the restricted residual after the coarse grid */
   /* correction should be zero.                                    */

   ML_Operator_Apply(Amat, lengf, sol, lengf, res);
   for ( i = 0; i < lengf; i++ ) res[i] = rhs[i] - res[i];

#ifdef ML_ANALYSIS
      if ( comm->ML_nprocs == 1 && curr->Pmat->to == NULL && lengf < 1000 )
      {
         fp = fopen("mlmatlab.m", "w");
         Nrows = lengf;
         fprintf(fp, "A = sparse(%d,%d);\n", Nrows, Nrows);
         allocated_space = 100;
         cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
         vals = (double *) ML_allocate(allocated_space*sizeof(double));
         for (i = 0; i < lengf; i++) {
            while(ML_Operator_Getrow(Amat,1,&i,allocated_space,cols,vals,&ncols)
                  == 0)
            {
               allocated_space = 2*allocated_space + 1;
               ML_free(vals); ML_free(cols);
               cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
               vals = (double *) ML_allocate(allocated_space*sizeof(double));
               if (vals == NULL) {
                  printf("Not enough space to get matrix row. Row length of\n");
                  printf("%d was not sufficient\n",(allocated_space-1)/2);
                  exit(1);
               }
            }
            for (j = 0; j < ncols; j++)
               fprintf(fp, "A(%d,%d)=%25.16e;\n",i+1,cols[j]+1,vals[j]);
         }
         fprintf(fp, "[eigv,eig]=eig(full(A));\n");
         fprintf(fp, "rhs=zeros(%d,1);\n", Nrows);
         for (j = 0; j < ncols; j++)
            fprintf(fp, "rhs(%d)=%25.16e;\n",j+1,res[j]);
         fprintf(fp, "res=eigv'*rhs;\n");
         fprintf(fp, "plot(res)\n");
         fclose(fp);
         printf("** AFTER  coarse grid correction -- \n");
         printf("** Now you can use matlab to call the file mlmatlab.m\n");
         printf("   and it should display the residual for each eigenvalue.\n");
         printf("Press y when you are done.");
         scanf("%s", instring);
     }
#endif

   ML_DVector_GetDataPtr(Rmat->from->Amat_Normalization,&normalscales);

   if ( normalscales != NULL )
      for ( i = 0; i < lengf; i++ ) res[i] = res[i]/normalscales[i];

   ML_Operator_ApplyAndResetBdryPts(Rmat, lengf, res, lengc, rhs2);

   ML_DVector_GetDataPtr(Rmat->to->Amat_Normalization,&normalscales);
   if ( normalscales != NULL )
      for ( i = 0; i < lengc; i++ ) rhs2[i] = rhs2[i] * normalscales[i];

   norm1 = sqrt(ML_gdot(lengc, rhs2, rhs2, comm));
   norm2 = sqrt(ML_gdot(lengf, res, res, comm));
   if (comm->ML_mypid == 0) printf("|R r| = %e, |r| =  %e\n",norm1, norm2);

#endif
      ML_free(sol2);
      ML_free(rhs2);

      ML_Smoother_Apply(post, lengf, sol, lengf, rhss, ML_NONZERO);
#ifdef ML_ANALYSIS
      if ( comm->ML_nprocs == 1 && curr->Pmat->to == NULL && lengf < 1000 )
      {
         ML_Operator_Apply(Amat, lengf, sol, lengf, res);
         for ( i = 0; i < lengf; i++ ) res[i] = rhss[i] - res[i];
         fp = fopen("mlmatlab.m", "w");
         Nrows = lengf;
         fprintf(fp, "A = sparse(%d,%d);\n", Nrows, Nrows);
         allocated_space = 100;
         cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
         vals = (double *) ML_allocate(allocated_space*sizeof(double));
         for (i = 0; i < lengf; i++) {
            while(ML_Operator_Getrow(Amat,1,&i,allocated_space,cols,vals,&ncols)
                  == 0)
            {
               allocated_space = 2*allocated_space + 1;
               ML_free(vals); ML_free(cols);
               cols = (int    *) ML_allocate(allocated_space*sizeof(int   ));
               vals = (double *) ML_allocate(allocated_space*sizeof(double));
               if (vals == NULL) {
                  printf("Not enough space to get matrix row. Row length of\n");
                  printf("%d was not sufficient\n",(allocated_space-1)/2);
                  exit(1);
               }
            }
            for (j = 0; j < ncols; j++)
               fprintf(fp, "A(%d,%d)=%25.16e;\n",i+1,cols[j]+1,vals[j]);
         }
         fprintf(fp, "[eigv,eig]=eig(full(A));\n");
         fprintf(fp, "rhs=zeros(%d,1);\n", Nrows);
         for (j = 0; j < ncols; j++)
            fprintf(fp, "rhs(%d)=%25.16e;\n",j+1,res[j]);
         fprintf(fp, "res=eigv'*rhs;\n");
         fprintf(fp, "plot(res)\n");
         fclose(fp);
         printf("** AFTER  postsmoothing -- \n");
         printf("** Now you can use matlab to call the file mlmatlab.m\n");
         printf("   and it should display the residual for each eigenvalue.\n");
         printf("Press y when you are done.");
         scanf("%s", instring);
     }
#endif
      ML_free(res);
   }

   ML_free(rhss);
   return(res_norm);
