#include "ml_agg_reitzinger.h"

int  ML_Gen_MGHierarchy_UsingReitzinger(ML *ml_edges, ML* ml_nodes, 
				      int fine_level, int incr_or_decrease,
				      ML_Aggregate *ag, ML_Operator *Tmat,
				      ML_Operator *Tmat_trans,
				      ML_Operator ***Tmat_array,
				      ML_Operator ***Tmat_trans_array)
  {
    int coarsest_level, counter, Nghost, i, *Tcoarse_bindx = NULL;
    int *Tcoarse_rowptr, nz_ptr, row_length, j, *bindx = NULL;
    double *Tcoarse_val = NULL, *node2proc, *val = NULL;
    int allocated = 0, lower;
    ML_Operator *Kn_coarse, *Rn_coarse, *Tcoarse, *Pn_coarse;
    ML_Operator *Pe, *Tcoarse_trans;  
    struct ML_CSR_MSRdata *csr_data;
    int Nlevels_nodal;
	double *fido,*yyy, dtemp,*xxx;

 if (incr_or_decrease != ML_DECREASING) pr_error("Hiptmair: Only ML_DESCREASING is supported\n");

  /********************************************************************/
  /* Set up the operators corresponding to regular unsmoothed         */
  /* aggregation on the nodal matrix.                                 */
  /*------------------------------------------------------------------*/

  Nlevels_nodal = ML_Gen_MGHierarchy_UsingAggregation(ml_nodes, fine_level, 
                                            ML_DECREASING, ag);
  coarsest_level = fine_level - Nlevels_nodal + 1; 
  *Tmat_array = (ML_Operator **) ML_allocate(sizeof(ML_Operator *)*
					     ml_nodes->ML_num_levels);
  *Tmat_trans_array = (ML_Operator **) ML_allocate(sizeof(ML_Operator *)*
					     ml_nodes->ML_num_levels);
  (*Tmat_array)[fine_level] = Tmat;
  (*Tmat_trans_array)[fine_level] = Tmat_trans;

  printf("ready to build T on the coarse grids\n");

  /********************************************************************/
  /*                 Build T on the coarse grid.                      */
  /* I'm not completely sure that the space I allocate is sufficient. */
  /* I believe we made a quick calculation and decided that the number*/
  /* of nonzeros in Kn_coarse was an upper bound for the number of    */
  /* nonzeros in T_coarse.                                            */
  /*------------------------------------------------------------------*/

  counter = 0;

  /* node2proc = nodal to processor lookup table */
  
  Kn_coarse = &(ml_nodes->Amat[coarsest_level]);
  Rn_coarse = &(ml_nodes->Rmat[fine_level]);

  if (Kn_coarse->getrow->pre_comm != NULL)
   Nghost = Kn_coarse->getrow->pre_comm->total_rcv_length;
  else Nghost = 0;

   node2proc = (double *) ML_allocate(sizeof(double)*(Nghost+
						  Kn_coarse->invec_leng+1));
   for (i = 0; i < Kn_coarse->invec_leng+Nghost; i++)
      node2proc[i] = (double) Kn_coarse->comm->ML_mypid;

   if (Kn_coarse->getrow->pre_comm != NULL)
   ML_exchange_bdry(node2proc, Kn_coarse->getrow->pre_comm,
                    Kn_coarse->outvec_leng, Kn_coarse->comm, ML_OVERWRITE);


  printf("check allocation sizes for Tcoarse in ml_agg_reitzinger.c\n");

  Tcoarse_bindx = (int    *) ML_allocate(Kn_coarse->N_nonzeros*10*sizeof(int));
  Tcoarse_val = (double *) ML_allocate(Kn_coarse->N_nonzeros*10*sizeof(double));
  Tcoarse_rowptr= (int    *) ML_allocate(Kn_coarse->outvec_leng*10*sizeof(int));
  Tcoarse_rowptr[0] = 0;
  nz_ptr = 0;
  /* printf("Tcoarse\n"); */
  for (i = 0; i < Kn_coarse->outvec_leng; i++)
  {
     ML_get_matrix_row(Kn_coarse,1, &i,&allocated,&bindx,&val,&row_length, 0);
     ML_az_sort(bindx, row_length, NULL, NULL);

     for (j = 0; j < row_length; j++)
	 {
       if (node2proc[i] == node2proc[bindx[j]])
	   {
          if (bindx[j] > i)
		  {
             Tcoarse_bindx[nz_ptr]  =  bindx[j];
             Tcoarse_val[nz_ptr++]  =  1.;
             Tcoarse_bindx[nz_ptr]  =  i;
             Tcoarse_val[nz_ptr++]  = -1.;
             Tcoarse_rowptr[counter+1] = nz_ptr;
             counter++;
			 /*
			 printf("row %d:    %d (%e)    %d (%e)\n",
			        i,Tcoarse_bindx[nz_ptr-2],Tcoarse_val[nz_ptr-2],
					Tcoarse_bindx[nz_ptr-1],Tcoarse_val[nz_ptr-1]);
		     */
          }
	   }
       else if (node2proc[i] < node2proc[bindx[j]])
	      {
             Tcoarse_bindx[nz_ptr]  =  bindx[j];
             Tcoarse_val[nz_ptr++]  =  1.;
             Tcoarse_bindx[nz_ptr]  =  i;
             Tcoarse_val[nz_ptr++]  = -1.;
             Tcoarse_rowptr[counter+1] = nz_ptr;
             counter++;
	      }
     }
  }

  csr_data = (struct ML_CSR_MSRdata *) ML_allocate(sizeof(struct 
							  ML_CSR_MSRdata));
  csr_data->columns = Tcoarse_bindx;
  csr_data->values  = Tcoarse_val;
  csr_data->rowptr  = Tcoarse_rowptr;

  Tcoarse = ML_Operator_Create(ml_edges->comm);
  ML_Operator_Set_ApplyFuncData( Tcoarse, Kn_coarse->outvec_leng, counter, 
                                  ML_EMPTY, csr_data, counter, NULL, 0);
  ML_Operator_Set_Getrow(Tcoarse, ML_EXTERNAL, counter, CSR_getrows);
  ML_Operator_Set_ApplyFunc(Tcoarse, ML_INTERNAL, CSR_matvec);

  ML_CommInfoOP_Clone(&(Tcoarse->getrow->pre_comm),
                      ml_nodes->Amat[coarsest_level].getrow->pre_comm);
  (*Tmat_array)[coarsest_level] = Tcoarse;

  printf("T built on coarse grid\n");
  /*ML_Operator_Print(Tcoarse,"Tcoarse");*/
  fflush(stdout);
	

  Pn_coarse = &(ml_nodes->Pmat[coarsest_level]);
  Rn_coarse = &(ml_nodes->Rmat[fine_level]);

  /********************************************************************/
  /* Fix P and R so that they are not normalized. This is so that we  */
  /* can use a matrix triple product combined with post-processing to */
  /* generate Pe. The general idea is that the matrix                 */
  /*                        T_h P_n T_H^*                             */
  /* is almost Pe. If we make sure that P_n contains 1's and -1's, the*/
  /* matrix triple product will yield a matrix with +/- 1 and +/- 2's.*/
  /* If we remove all the 1's and divide the 2's by 2. we arrive at Pe*/
  /*------------------------------------------------------------------*/

/*
  Pn_coarse = &(ml_nodes->Pmat[coarsest_level]);
  csr_data = (struct ML_CSR_MSRdata *) Pn_coarse->data;
  for (i = 0; i < csr_data->rowptr[Pn_coarse->outvec_leng]; i++)
  {
    if (csr_data->values[i] < 0) csr_data->values[i] = -1.;
    else if (csr_data->values[i] > 0) csr_data->values[i] = 1.;
  }
*/

/* We want all +1 entries to avoid changing the sign of entries in Pe. */
  Pn_coarse = &(ml_nodes->Pmat[coarsest_level]);
  csr_data = (struct ML_CSR_MSRdata *) Pn_coarse->data;
  for (i = 0; i < csr_data->rowptr[Pn_coarse->outvec_leng]; i++)
    if (csr_data->values[i] != 0) csr_data->values[i] = 1.;

  Rn_coarse = &(ml_nodes->Rmat[fine_level]);
  csr_data = (struct ML_CSR_MSRdata *) Rn_coarse->data;
  for (i = 0; i < csr_data->rowptr[Rn_coarse->outvec_leng]; i++)
  {
    if (csr_data->values[i] < 0.) csr_data->values[i] = -1.;
    else if (csr_data->values[i] > 0) csr_data->values[i] = 1.;
  }

  /********************************************************************/
  /* Create Tcoarse_trans.                                            */
  /*------------------------------------------------------------------*/

  Tcoarse_trans = ML_Operator_Create(ml_edges->comm);
  /* ML_Operator_Print(Tcoarse,"Tcoarse"); */
  ML_Operator_Transpose_byrow(Tcoarse, Tcoarse_trans);
  (*Tmat_trans_array)[coarsest_level] = Tcoarse_trans;

  /********************************************************************/
  /* Here is some code that might work someday to generate Pe without */
  /* doing a matrix triple product.                                   */
  /*------------------------------------------------------------------*/
  /*
  csr_data = (struct ML_CSR_MSRdata *) Pn_coarse->data;
  for (i = 0; i < Tmat->outvec_leng; i++) {
    ML_get_matrix_row(Tmat, 1, &i, &allocated, &bindx, &val, &row_length, 0);
    if (row_length == 2) {
     agg1 = csr_data->columns[bindx[0]];
     agg2 = csr_data->columns[bindx[1]];
     printf("%d: agg1 and agg2 %d %d   | %d %d\n",i,agg1,agg2,bindx[0],bindx[1]);
    
     agg1 = aggr_num(bindx[0], PnMat);
     agg2 = aggr_num(bindx[1], PnMat);
     if (agg1 != agg2) {
        printf("Pe(%d,%d) = something;\n",i,thecolumn);

        To get the column suppose we store for each row in Kn_coarse the counter above.
        Thus, stored_counter[i] indicates the lowest coarse grid edge that is 
        assicated with coarse grid node i.
        To get the edge number we compute 
               min_agg = min(agg1,agg2), 
               max_agg = max(agg1,agg2);
        The column that we are looking for is given by 
          stored_counter[min_agg] + number of nonzeros in row min_agg of
		  Kn_coarse that are greater than min_agg and less than max_agg.


        'something' is either +1 or -1. To get the sign look at
        Tcoarse(thecolumn,) and Tmat(i,). My guess is that it is related to 
        whether agg1 is greater than agg2 (assuming that val[0] =1)
         }
  }
}
  */


  /********************************************************************/
  /* Matrix triple product for Pe where result is stuffed into MG     */
  /* hierarchy of the edge system.                                    */
  /*------------------------------------------------------------------*/

/*
  ML_Operator_Print( &(ml_nodes->Pmat[coarsest_level]),"Pn");
  exit(1);
*/
  ML_rap(Tmat, &(ml_nodes->Pmat[coarsest_level]), Tcoarse_trans, 
	 &(ml_edges->Pmat[coarsest_level]),ML_CSR_MATRIX);

  Pe = &(ml_edges->Pmat[coarsest_level]);
  csr_data = (struct ML_CSR_MSRdata *) Pe->data;

  /********************************************************************/
  /* weed out the 1's and divide the 2's by -2.                       */
  /* Note: my signs seem flipped from what we thought the triple      */
  /* matrix product would do. This could be a function of the Pnmat   */
  /* normalization.                                                   */
  /* Finally, convert the matrix to CSR and strip out the zeros.      */
  /* and make sure that it has the level transfer stuff to be in the  */
  /* MG grid hierarchy.                                               */
  /*------------------------------------------------------------------*/

#ifdef prepostcheck
   printf("Checking product Pe * e_i before post-processing\n");
   yyy = (double *) malloc( Pe->outvec_leng*sizeof(double));
   fido = (double *) malloc( Pe->invec_leng*sizeof(double));
/*
   printf("%d: Pe->invec_leng = %d\n",Pe->comm->ML_mypid,Pe->invec_leng);
   printf("%d: Tcoarse->outvec_leng = %d\n",Tcoarse->comm->ML_mypid,
          Tcoarse->outvec_leng);
   exit(1);
*/

   printf("Pe->invec_leng = %d\n",Pe->invec_leng);
   for (i=0; i< 137; i++)
   {
      for (j=0; j< Pe->invec_leng; j++) fido[j] = 0.;
	  if (Pe->comm->ML_nprocs == 1)
	     fido[i] = 1;
	  else
	  {
	     if ((Pe->comm->ML_mypid == 0) && (i < 120))
		    fido[i] = 1;
		 else if (Pe->comm->ML_mypid == 1 && i >= 120)
		    fido[i-120] = 1;
      }
	  /*
	  if (i==119)
      {
	     printf("e_119\n");
         for (j=0; j<Pe->invec_leng; j++)
		    printf("%d: e_119(%d) = %e\n", Pe->comm->ML_mypid,j,fido[j]);
	  }
	  if (i==120)
	  {
	     printf("e_120\n");
         for (j=0; j<Pe->invec_leng; j++)
		    printf("%d: e_120(%d) = %e\n", Pe->comm->ML_mypid,j,fido[j]);
      }
	  */

      ML_Operator_Apply(Pe, Pe->outvec_leng, fido,
                        Pe->outvec_leng,yyy);
      dtemp = ML_gdot(Pe->outvec_leng, yyy, yyy, Pe->comm);
      printf("norm(P(:,%d)) = %e\n",i,dtemp); 
      /*printf("%d: ||P_e * e_%d||^2 = %e\n",Pe->comm->ML_mypid,i,dtemp);*/
	  /*
	  if (i==50)
	  {
         for (j=0; j<Pe->invec_leng; j++)
		    printf("%d: yyy(%d) = %e\n", Pe->comm->ML_mypid,j,yyy[j]);
      }
      */

/*
      printf("%d: (%d) %e\n",ml_edges->comm->ML_mypid,i,dtemp);
*/

/*
      dtemp = sqrt(ML_gdot(Tcoarse->invec_leng, fido, fido, ml_edges->comm));
      printf("(%d): %e\n",i,dtemp); 
*/

   }

   exit(1);

#endif /* prepostcheck */

  for (j = 0; j < csr_data->rowptr[Pe->outvec_leng] ; j++) {
    if (csr_data->values[j] == 2) csr_data->values[j] = -1;
    else if (csr_data->values[j] == -2) csr_data->values[j] = 1;
    else if (csr_data->values[j] == -1) csr_data->values[j] = 0;
    else if (csr_data->values[j] ==  1) csr_data->values[j] = 0;
    else if (csr_data->values[j] != 0.0) printf("huh\n");
  }
 
  /*******************************************************************/
  /* weed out zeros in Pe.                                           */
  /*-----------------------------------------------------------------*/
 
  lower = csr_data->rowptr[0];
  nz_ptr = 0;
  for (i = 0; i < Pe->outvec_leng; i++) {
    for (j = lower; j < csr_data->rowptr[i+1]; j++) {
      if (csr_data->values[j] != 0.0) nz_ptr++;
    }
    lower = csr_data->rowptr[i+1];
    csr_data->rowptr[i+1] = nz_ptr;
  }
  nz_ptr = 0;
  for (i = 0; i < lower; i++) {
    if (csr_data->values[i] != 0.) {
      csr_data->values[nz_ptr] = csr_data->values[i];
      csr_data->columns[nz_ptr] = csr_data->columns[i];
      nz_ptr++;
    }
  }
 
  Pe->getrow->external = CSR_getrows;
  Pe->getrow->internal = NULL;
  Pe->getrow->ML_id    = ML_EXTERNAL;
  Pe->matvec->internal = CSR_matvec;
  Pe->matvec->external = NULL;
  Pe->matvec->ML_id = ML_INTERNAL;


/****************** Check the construction of Pe ***********************/

#ifdef postprocessscheck
   printf("Checking product Pe * e_i\n");
   yyy = (double *) malloc( Pe->outvec_leng*sizeof(double));
   fido = (double *) malloc( Pe->invec_leng*sizeof(double));
/*
   printf("%d: Pe->invec_leng = %d\n",Pe->comm->ML_mypid,Pe->invec_leng);
   printf("%d: Tcoarse->outvec_leng = %d\n",Tcoarse->comm->ML_mypid,
          Tcoarse->outvec_leng);
   exit(1);
*/

   printf("Pe->invec_leng = %d\n",Pe->invec_leng);
   for (i=0; i< 137; i++)
   {
      for (j=0; j< Pe->invec_leng; j++) fido[j] = 0.;
	  if (Pe->comm->ML_nprocs == 1)
	     fido[i] = 1;
	  else
	  {
	     if ((Pe->comm->ML_mypid == 0) && (i < 120))
		    fido[i] = 1;
		 else if (Pe->comm->ML_mypid == 1 && i >= 120)
		    fido[i-120] = 1;
      }
	  /*
	  if (i==84)
	  {
	     printf("e_120\n");
         for (j=0; j<Pe->invec_leng; j++)
		    printf("%d: e_84(%d) = %e\n", Pe->comm->ML_mypid,j,fido[j]);
      }
	  if (i==119)
      {
	     printf("e_119\n");
         for (j=0; j<Pe->invec_leng; j++)
		    printf("%d: e_119(%d) = %e\n", Pe->comm->ML_mypid,j,fido[j]);
	  }
	  if (i==120)
	  {
	     printf("e_120\n");
         for (j=0; j<Pe->invec_leng; j++)
		    printf("%d: e_120(%d) = %e\n", Pe->comm->ML_mypid,j,fido[j]);
      }
	  */

      ML_Operator_Apply(Pe, Pe->invec_leng, fido,
                        Pe->outvec_leng,yyy);
      dtemp = ML_gdot(Pe->outvec_leng, yyy, yyy, Pe->comm);
      printf("norm(P(:,%d)) = %e\n",i,dtemp); 
      /*printf("%d: norm(P(:,%d)) = %e\n",Pe->comm->ML_mypid,i,dtemp); */
	  /*
	  if (i==50)
	  {
         for (j=0; j<Pe->invec_leng; j++)
		    printf("%d: yyy(%d) = %e\n", Pe->comm->ML_mypid,j,yyy[j]);
      }
      */

/*
      printf("%d: (%d) %e\n",ml_edges->comm->ML_mypid,i,dtemp);
*/

/*
      dtemp = sqrt(ML_gdot(Tcoarse->invec_leng, fido, fido, ml_edges->comm));
      printf("(%d): %e\n",i,dtemp); 
*/

   }

   fflush(stdout);
  exit(1);
#endif /*postprocesscheck*/

  ML_Operator_Set_1Levels(&(ml_edges->Pmat[coarsest_level]),
			  &(ml_edges->SingleLevel[coarsest_level]), 
			  &(ml_edges->SingleLevel[fine_level]));
  ML_Gen_Restrictor_TransP(ml_edges, fine_level, coarsest_level);
  ML_Gen_AmatrixRAP(ml_edges, fine_level, coarsest_level);
  return(Nlevels_nodal);
  }

