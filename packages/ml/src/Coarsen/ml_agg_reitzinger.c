#include "ml_agg_reitzinger.h"
#define NEW_T_PE

extern int ML_Gen_SmoothPnodal(ML *ml,int level, int clevel, void *data,
			       double smoothP_damping_factor,
			       ML_Operator *SPn_mat);

int  ML_Gen_MGHierarchy_UsingReitzinger(ML *ml_edges, ML* ml_nodes, 
                    int fine_level, int incr_or_decrease,
                    ML_Aggregate *ag, ML_Operator *Tmat,
                    ML_Operator *Tmat_trans,
                    ML_Operator ***Tmat_array,
                    ML_Operator ***Tmat_trans_array,
                    int smooth_flag, double smooth_factor)
{
  int coarsest_level, counter, Nghost, i, *Tcoarse_bindx = NULL;
  int *Tcoarse_rowptr, nz_ptr, row_length, j, *bindx = NULL;
  double *Tcoarse_val = NULL, *node2proc, *val = NULL;
  int allocated = 0, lower;
  ML_Operator *Kn_coarse, *Rn_coarse, *Tcoarse, *Pn_coarse;
  ML_Operator *Pe, *Tcoarse_trans, *Tfine;
#ifdef LEASTSQ_SERIAL
  ML_Operator *SPn_mat;
#endif
  struct ML_CSR_MSRdata *csr_data;
  int Nlevels_nodal, grid_level;
  int created_ag_obj = 0;
  double *vec, *Tcoarse_vec, *Pn_vec, *Tfine_Pn_vec;
  int i1, old_nzptr, i3, *index;
  int *encoded_dir_node, *temp_bindx, Npos_dirichlet = 0, Nneg_dirichlet = 0;
  int Nnondirichlet, Nnz_finegrid, Nnz_allgrids;
  double *pos_coarse_dirichlet, *neg_coarse_dirichlet, *temp_val, d1, d2;
  double droptol;

  ML_CommInfoOP *getrow_comm; 
  int  N_input_vector;
  int  bail_flag;

  /*
  double *fido,*yyy, *vvv, dtemp;
  struct aztec_context *temp;
  */
  int nzctr;

  if (incr_or_decrease != ML_DECREASING)
    pr_error("Hiptmair: Only ML_DECREASING is supported\n");


  if (ag == NULL)
  {
     created_ag_obj = 1;
     ML_Aggregate_Create(&ag);
     ML_Aggregate_Set_CoarsenScheme_MIS(ag);
     ML_Aggregate_Set_DampingFactor(ag,0.0);
     ML_Aggregate_Set_Threshold( ag, 0.0);
  }
  Tfine = Tmat;
  /********************************************************************/
  /* Set up the operators corresponding to regular unsmoothed         */
  /* aggregation on the nodal matrix.                                 */
  /*------------------------------------------------------------------*/
  Nnz_finegrid = ml_edges->Amat[fine_level].N_nonzeros; 
  Nnz_allgrids = ml_edges->Amat[fine_level].N_nonzeros;

  Nlevels_nodal = ML_Gen_MGHierarchy_UsingAggregation(ml_nodes, fine_level, 
                                            ML_DECREASING, ag);

  coarsest_level = fine_level - Nlevels_nodal + 1; 
  i = ml_nodes->Amat[coarsest_level].invec_leng;
  ML_gsum_vec_int(&i,&j,1,ml_nodes->comm);
  if (i <= 1) {  /* no edges so cut out last level */
    Nlevels_nodal--;
    coarsest_level = fine_level - Nlevels_nodal + 1; 
  }

  *Tmat_array = (ML_Operator **) ML_allocate(sizeof(ML_Operator *)*
                         ml_nodes->ML_num_levels);
  *Tmat_trans_array = (ML_Operator **) ML_allocate(sizeof(ML_Operator *)*
                         ml_nodes->ML_num_levels);
  for (i = 0; i < ml_nodes->ML_num_levels; i++)
  {
     (*Tmat_array)[i] = NULL;
     (*Tmat_trans_array)[i] = NULL;
  }
  (*Tmat_array)[fine_level] = Tfine;
  (*Tmat_trans_array)[fine_level] = Tmat_trans;

  /********************************************************************/
  /*                 Build T on the coarse grid.                      */
  /*------------------------------------------------------------------*/

  /* main loop through grid hierarchy */
  for (grid_level = fine_level-1; grid_level >= coarsest_level ;grid_level--)
  {
     Kn_coarse = &(ml_nodes->Amat[grid_level]);
     Rn_coarse = &(ml_nodes->Rmat[grid_level+1]);
   
     if (Kn_coarse->getrow->pre_comm != NULL)
        Nghost = Kn_coarse->getrow->pre_comm->total_rcv_length;
     else Nghost = 0;
   
     /* node2proc = node to processor lookup table
        Main idea:
        
        1) Record the proc. id of all nodes (local + ghost) on the processor.
        2) Exchange boundary information.  Ghost nodes should now have the
           proc.  id of the "owning" processor.
   
        Tcoarse must be constructed in a way that processors do not duplicate
        work.  A processor owns a row if both nodes are local.  If one node
        is a ghost node, than the processor owns row i (edge i with endpoints
        k & j) if proc_id(k) < proc_id(j).  Otherwise, processor proc_id(j)
        will own row i.
     */
   
     node2proc = (double *) ML_allocate(sizeof(double)*(Nghost+
                            Kn_coarse->invec_leng+1));
     for (i = 0; i < Kn_coarse->invec_leng+Nghost; i++)
        node2proc[i] = (double) Kn_coarse->comm->ML_mypid;
   
     if (Kn_coarse->getrow->pre_comm != NULL)
     ML_exchange_bdry(node2proc, Kn_coarse->getrow->pre_comm,
                      Kn_coarse->outvec_leng, Kn_coarse->comm,
                      ML_OVERWRITE,NULL);
   
#ifdef DEBUG_T_BUILD
     printf("\n\n%d: Kn_coarse->N_nonzeros = %d "
            "Kn_coarse->invec_leng+Nghost = %d\n"
            "Kn_coarse->invec_leng = %d\n\n",
            Kn_coarse->comm->ML_mypid, Kn_coarse->N_nonzeros,
            Kn_coarse->invec_leng+Nghost, Kn_coarse->invec_leng);
#endif /* ifdef DEBUG_T_BUILD */

#ifdef NEW_T_PE
     /* let's figure out who corresponds to a Dirichlet point */

     allocated = 100;
     temp_bindx = (int    *)  ML_allocate( allocated*sizeof(int   ));
     temp_val   = (double *)  ML_allocate( allocated*sizeof(double));
     encoded_dir_node = (int *) ML_allocate((Tfine->outvec_leng+1)*sizeof(int));

     for (i = 0 ; i < Tfine->getrow->Nrows; i++) {
       ML_get_matrix_row(Tfine, 1, &i, &allocated, &temp_bindx, &temp_val,
			 &row_length, 0);
       if (row_length == 2) {
	 if (temp_val[1] == 0.) row_length--;
	 if (temp_val[0] == 0.) {
	   row_length--;
	   if (row_length != 0) {
	     temp_bindx[1] = temp_bindx[0];
	     temp_val[1] = temp_val[0];
	   }
	 }
       }
       if (row_length == 1) {
	 if      (temp_val[0] ==  1.) encoded_dir_node[i] = (1 + temp_bindx[0]);
	 else if (temp_val[0] == -1.) encoded_dir_node[i] = -(1 + temp_bindx[0]);
	 else printf("Warning uknown value T(%d,%d) = %e\n",
		     i,temp_bindx[0],temp_val[0]); 
       }
       else encoded_dir_node[i] = 0;
     }
     ML_free(temp_bindx);
     ML_free(temp_val);
     allocated = 0;


     vec = (double *) ML_allocate(sizeof(double)*(Rn_coarse->invec_leng
						     +1));

     pos_coarse_dirichlet = (double *) ML_allocate(sizeof(double)*(Rn_coarse->outvec_leng
						     +1));
     if (pos_coarse_dirichlet == NULL) {
        printf("\n\nML_Gen_MGHierarchy_UsingReitzinger: Not enough space"
               " allocated to check T.\n\n");
        exit(1);
     }
     for (i = 0; i < Rn_coarse->invec_leng; i++) vec[i] = 0.;
     for (i = 0; i < Tfine->outvec_leng; i++) {
       if(encoded_dir_node[i] > 0) vec[encoded_dir_node[i]-1] = 1.;
     }
     ML_Operator_Apply(Rn_coarse, Rn_coarse->invec_leng, vec,
		       Rn_coarse->outvec_leng,pos_coarse_dirichlet);
     Npos_dirichlet = 0;
     for (i = 0; i < Rn_coarse->outvec_leng; i++) {
       if (pos_coarse_dirichlet[i] != 0) {
	 Npos_dirichlet++;
       }
     }

     neg_coarse_dirichlet = (double *) ML_allocate(sizeof(double)*(Rn_coarse->outvec_leng
						     +1));
     if (neg_coarse_dirichlet == NULL) {
        printf("\n\nML_Gen_MGHierarchy_UsingReitzinger: Not enough space"
               " allocated to check T.\n\n");
        exit(1);
     }
     for (i = 0; i < Rn_coarse->invec_leng; i++) vec[i] = 0.;
     for (i = 0; i < Tfine->outvec_leng; i++) {
       if(encoded_dir_node[i] < 0) vec[-encoded_dir_node[i]-1] = 1.;
     }
     ML_Operator_Apply(Rn_coarse, Rn_coarse->invec_leng, vec,
		       Rn_coarse->outvec_leng,neg_coarse_dirichlet);
     ML_free(vec);
     Nneg_dirichlet = 0;
     for (i = 0; i < Rn_coarse->outvec_leng; i++) {
       if (neg_coarse_dirichlet[i] != 0) {
	 Nneg_dirichlet++;
       }
     }
#endif


     Tcoarse_bindx =(int *)
                    ML_allocate((2*Kn_coarse->N_nonzeros + Nneg_dirichlet +
				 Npos_dirichlet + 1)*sizeof(int) );
     Tcoarse_val = (double *)
                    ML_allocate((2*Kn_coarse->N_nonzeros + Nneg_dirichlet +
				 Npos_dirichlet + 1)*sizeof(double) );
     Tcoarse_rowptr= (int *)
                     ML_allocate((Kn_coarse->N_nonzeros + Nneg_dirichlet +
				  Npos_dirichlet + 1)*sizeof(int));
     Tcoarse_rowptr[0] = 0;
     counter = 0; nz_ptr = 0;
     nzctr = 0;
     for (i = 0; i < Kn_coarse->outvec_leng; i++)
     {
        ML_get_matrix_row(Kn_coarse,1, &i,&allocated,
                          &bindx,&val,&row_length, 0);
        ML_az_sort(bindx, row_length, NULL, NULL);
        nzctr += row_length;
   
        /* Step through unknowns bindx[j] connected to unknown i. */
        for (j = 0; j < row_length; j++)
        {
          {
             /* If nodes i and bindx[j] are owned by same processor ... */
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
                }
             }
             /* If node i is owned by a smaller processor than
                node bindx[j] ... */
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
     }

     if (nzctr > Kn_coarse->N_nonzeros && Kn_coarse->comm->ML_mypid == 0)
     {
        printf("\n\nML_Gen_MGHierarchy_UsingReitzinger: Not enough space"
               " allocated to build T.\n\n");
        exit(1);
     }
#ifdef NEW_T_PE
     for (i = 0; i < Rn_coarse->outvec_leng; i++) {
       if (pos_coarse_dirichlet[i] != 0) {
	 Tcoarse_bindx[nz_ptr] = i;
	 Tcoarse_val[nz_ptr++] = 1;
	 Tcoarse_rowptr[counter+1] = nz_ptr;
	 counter++;
       }
     }
     for (i = 0; i < Rn_coarse->outvec_leng; i++) {
       if (neg_coarse_dirichlet[i] != 0) {
	 Tcoarse_bindx[nz_ptr] = i;
	 Tcoarse_val[nz_ptr++] = -1;
	 Tcoarse_rowptr[counter+1] = nz_ptr;
	 counter++;
       }
     }
     ML_free(pos_coarse_dirichlet);
     ML_free(neg_coarse_dirichlet);
#endif

#ifdef DEBUG_T_BUILD
     else
        if (Kn_coarse->comm->ML_mypid == 0 && grid_level == 7)
        printf("%d (%d): ieqj = %d,\n"
               "Kn_coarse->N_nonzeros = %d, "
               "expected nnz (calc. from bindx) = %d,\n"
               "actual nnz(Tcoarse) = %d,\n"
               " invec_leng = %d,"
               " nghost = %d\n",
               Kn_coarse->comm->ML_mypid,
               grid_level,
               ieqj,
               Kn_coarse->N_nonzeros,
               nzctr,
               nz_ptr,
               Kn_coarse->invec_leng,
               Nghost);
     fflush(stdout);
#endif /* ifdef DEBUG_T_BUILD */

     ML_free(node2proc);
   
     csr_data = (struct ML_CSR_MSRdata *) ML_allocate(sizeof(struct 
                                 ML_CSR_MSRdata));
     csr_data->columns = Tcoarse_bindx;
     csr_data->values  = Tcoarse_val;
     csr_data->rowptr  = Tcoarse_rowptr;
   
     Tcoarse = ML_Operator_Create(ml_edges->comm);
     Tcoarse->data_destroy = ML_CSR_MSRdata_Destroy;
     ML_Operator_Set_ApplyFuncData( Tcoarse, Kn_coarse->outvec_leng, counter, 
                                     ML_EMPTY, csr_data, counter, NULL, 0);
     ML_Operator_Set_Getrow(Tcoarse, ML_EXTERNAL, counter, CSR_getrows);
     ML_Operator_Set_ApplyFunc(Tcoarse, ML_INTERNAL, CSR_matvec);
   
     ML_CommInfoOP_Clone(&(Tcoarse->getrow->pre_comm),
                         ml_nodes->Amat[grid_level].getrow->pre_comm);
     (*Tmat_array)[grid_level] = Tcoarse;
   
     Pn_coarse = &(ml_nodes->Pmat[grid_level]);
     Rn_coarse = &(ml_nodes->Rmat[grid_level+1]);

     /* Test every coarse grid edge to make sure that it is really */
     /* necessary. This test is needed because Kn has connections  */
     /* that do not correspond to edges in Ke. For example,        */
     /* consider the 2D Ke stencil:                                */
     /*                    1 ------ 2                              */
     /*                    |        |                              */
     /*                    |        |                              */
     /*                    3 ------ 4                              */
     /*                    |        |                              */
     /*                    |        |                              */
     /*                    5 ------ 6                              */
     /* To properly perform the aggregation, the graph of Kn       */
     /* will include the edges (2,3), (1,4), (3,6), (4,5)          */
     /* even though these edges do not exist. Since Kn is used     */
     /* to construct Tcoarse, this can lead to extraneous edges    */
     /* in Tcoarse. These extra edges can lead to Ke(i,i) = 0      */
     /* on the next coarser matrix. To avoid this, we attempt to   */
     /* eliminate these edges here.                                */

     /* Extraneous edges are detected by comparing                 */
     /*     T_H vec      vs.   T_h P_n vec                         */
     /* If an entry in T_H vec is not found in T_h P_n vec, it     */
     /* means that this edge can be removed from T_H. Note: right  */
     /* now we are doing an == with floating point numbers. This   */
     /* should work ... but it might break on some hardwares.      */

     vec = (double *) ML_allocate(sizeof(double)*(Tcoarse->invec_leng+1));
     Tcoarse_vec = (double *) ML_allocate(sizeof(double)*(Tcoarse->outvec_leng
							  + 1));
     if ((vec == NULL) || (Tcoarse_vec == NULL)) {
        printf("\n\nML_Gen_MGHierarchy_UsingReitzinger: Not enough space"
               " allocated to check T.\n\n");
        exit(1);
     }

     ML_random_vec(vec, Tcoarse->invec_leng, ml_edges->comm);
     ML_Operator_Apply(Tcoarse, Tcoarse->invec_leng, vec,
		       Tcoarse->outvec_leng,Tcoarse_vec);

     for (i1 = 0; i1 < Tcoarse->outvec_leng; i1++) 
       if (Tcoarse_vec[i1] < 0) Tcoarse_vec[i1] = -Tcoarse_vec[i1];

     Pn_coarse->matvec->internal = CSR_ones_matvec; /* turn of the scaling */
                                                    /* in Pn for the test. */
     Pn_vec = (double *) ML_allocate(sizeof(double)*(Pn_coarse->outvec_leng
						     +1));
     if (Pn_vec == NULL) {
        printf("\n\nML_Gen_MGHierarchy_UsingReitzinger: Not enough space"
               " allocated to check T.\n\n");
        exit(1);
     }


     ML_Operator_Apply(Pn_coarse, Pn_coarse->invec_leng, vec,
		       Pn_coarse->outvec_leng,Pn_vec);
     ML_free(vec);
     Tfine_Pn_vec = (double *) ML_allocate(sizeof(double)*(Tfine->outvec_leng
							   +1));
     if (Tfine_Pn_vec == NULL) {
        printf("\n\nML_Gen_MGHierarchy_UsingReitzinger: Not enough space"
               " allocated to check T.\n\n");
        exit(1);
     }
     ML_Operator_Apply(Tfine, Tfine->invec_leng, Pn_vec,
		       Tfine->outvec_leng,Tfine_Pn_vec);
     ML_free(Pn_vec);
     Pn_coarse->matvec->internal = CSR_matvec;
     for (i1 = 0; i1 < Tfine->outvec_leng; i1++) 
       if (Tfine_Pn_vec[i1] < 0) Tfine_Pn_vec[i1] = -Tfine_Pn_vec[i1];

     /* First sort the two vectors and then compare them.    */
     /* Mark Tcoarse_vec[i] with a '1.' if its corresponding */
     /* edge can be removed.                                 */

     index = (int *) ML_allocate((Tcoarse->outvec_leng+1)*sizeof(int));
     for (i1 = 0; i1 < Tcoarse->outvec_leng; i1++) index[i1] = i1;
                            /* use index to keep track of which edge */
                            /* is which in the sorted Tcoarse_vec    */
     ML_az_dsort2(Tcoarse_vec,Tcoarse->outvec_leng, index);
     ML_az_dsort(Tfine_Pn_vec,Tfine->outvec_leng);

     i3 = 0;
     for (i1 = 0; i1 < Tcoarse->outvec_leng; i1++) {
       while ( (i3 != Tfine->outvec_leng) && 
	       (Tcoarse_vec[i1] > Tfine_Pn_vec[i3])) {
	 i3++;
       }
       if ((i3 < Tfine->outvec_leng)&&(Tcoarse_vec[i1] == Tfine_Pn_vec[i3])) {
	 Tcoarse_vec[i1] = 0.;
       }
       else Tcoarse_vec[i1] = 1.;
     }
     ML_free(Tfine_Pn_vec);
     ML_az_sort(index, Tcoarse->outvec_leng, NULL, Tcoarse_vec);
     ML_free(index);

     csr_data = (struct ML_CSR_MSRdata *) Tcoarse->data;
     nz_ptr = 0;
     old_nzptr = 0;
     counter = 0;

     for (i1 = 0; i1 < Tcoarse->outvec_leng; i1++) {
       if (Tcoarse_vec[i1] == 1.) {
	 old_nzptr = csr_data->rowptr[i1+1];
       }
       else {
	 row_length = csr_data->rowptr[i1+1] - old_nzptr;
	 csr_data->rowptr[counter+1] = row_length + csr_data->rowptr[counter];
	 counter++;
	 for (i3 = 0; i3 < row_length; i3++) { 
	   csr_data->columns[nz_ptr] =  csr_data->columns[old_nzptr];
	   csr_data->values[nz_ptr++] =  csr_data->values[old_nzptr++];
	 }
       }
     }
     ML_free(Tcoarse_vec);
     Tcoarse->outvec_leng = counter;
     Tcoarse->getrow->Nrows = counter;

     csr_data->rowptr = (int *) realloc(csr_data->rowptr,
					sizeof(int)*(counter+1));
     csr_data->values = (double *) realloc(csr_data->values,
					   sizeof(double)*nz_ptr); 
     csr_data->columns = (int *) realloc(csr_data->columns,
					   sizeof(int)*nz_ptr); 
   
     /********************************************************************/
     /* Fix P and R so that they are not normalized. This is so that we  */
     /* can use a matrix triple product combined with post-processing to */
     /* generate Pe. The general idea is that the matrix                 */
     /*                        T_h P_n T_H^*                             */
     /* is almost Pe. If we make sure that P_n contains 1's and -1's, the*/
     /* matrix triple product will yield a matrix with +/- 1 and +/- 2's.*/
     /* If we remove all the 1's and divide the 2's by 2. we arrive at Pe*/
     /*------------------------------------------------------------------*/
   
     /* We want all +1 entries to avoid changing sign of entries in Pe. */
     Pn_coarse = &(ml_nodes->Pmat[grid_level]);
     csr_data = (struct ML_CSR_MSRdata *) Pn_coarse->data;
     for (i = 0; i < csr_data->rowptr[Pn_coarse->outvec_leng]; i++)
       if (csr_data->values[i] != 0) csr_data->values[i] = 1.;
   
     Rn_coarse = &(ml_nodes->Rmat[grid_level+1]);
     csr_data = (struct ML_CSR_MSRdata *) Rn_coarse->data;
   
     for (i = 0; i < csr_data->rowptr[Rn_coarse->outvec_leng]; i++)
       if (csr_data->values[i] != 0) csr_data->values[i] = 1.;
   
     /********************************************************************/
     /* Create Tcoarse_trans.                                            */
     /*------------------------------------------------------------------*/
   
     Tcoarse_trans = ML_Operator_Create(ml_edges->comm);

     /* Check that both dimensions of T are strictly greater than 0. 
        If not, clean up & break from main loop. */
     i = Tcoarse->outvec_leng;
     ML_gsum_vec_int(&i,&j,1,ml_nodes->comm);
     if (i==0)
     {
        if (Tcoarse->comm->ML_mypid == 0 && 5 < ML_Get_PrintLevel()) {
           printf("(%d) Bailing from AMG hierarchy build on level %d of "
                 " of levels %d down to %d because Tcoarse has zero rows....\n",
                 ml_edges->comm->ML_mypid,grid_level,fine_level,coarsest_level);
           fflush(stdout);
        }
        /*Don't delete Tcoarse because it is necessary for eigenvalue estimate
          in Hiptmair setup.*/
#ifdef NEW_T_PE
        ML_free(encoded_dir_node);
#endif
        /* Current level "grid_level" cannot be used b/c Tcoarse_trans
           would be used in Hiptmair smoother generation, but Tcoarse_trans
           hasn't been calculated. Hence no "+1".*/
        Nlevels_nodal = fine_level - grid_level;
        coarsest_level = grid_level + 1;
        if (Tcoarse->comm->ML_mypid == 0 && 5 < ML_Get_PrintLevel()) {
           printf("(%d) In ML_Gen_MGHierarchy_UsingReitzinger, "
               "Nlevels_nodal = %di, fine_level = %d, coarsest_level = %d\n",
              ml_nodes->comm->ML_mypid,Nlevels_nodal,fine_level,coarsest_level);
           fflush(stdout);
        }
        break; /* from main loop */
     }

     ML_Operator_Transpose_byrow(Tcoarse, Tcoarse_trans);
     (*Tmat_trans_array)[grid_level] = Tcoarse_trans;
   
     /********************************************************************/
     /* Here is some code that might work someday to generate Pe without */
     /* doing a matrix triple product.                                   */
     /*------------------------------------------------------------------*/
     /*
     csr_data = (struct ML_CSR_MSRdata *) Pn_coarse->data;
     for (i = 0; i < Tfine->outvec_leng; i++)
     {
       ML_get_matrix_row(Tfine, 1, &i, &allocated, &bindx, &val,
                         &row_length, 0);
       if (row_length == 2)
       {
        agg1 = csr_data->columns[bindx[0]];
        agg2 = csr_data->columns[bindx[1]];
        printf("%d: agg1 and agg2 %d %d   | %d %d\n",i,agg1,
               agg2,bindx[0],bindx[1]);
       
        agg1 = aggr_num(bindx[0], PnMat);
        agg2 = aggr_num(bindx[1], PnMat);
        if (agg1 != agg2) {
           printf("Pe(%d,%d) = something;\n",i,thecolumn);
   
           To get the column suppose we store for each row in Kn_coarse
           the counter above.
           Thus, stored_counter[i] indicates the lowest coarse grid edge that
           is associated with coarse grid node i.
           To get the edge number we compute 
                  min_agg = ML_min(agg1,agg2), 
                  max_agg = ML_max(agg1,agg2);
           The column that we are looking for is given by 
             stored_counter[min_agg] + number of nonzeros in row min_agg of
             Kn_coarse that are greater than min_agg and less than max_agg.
   
   
           'something' is either +1 or -1. To get the sign look at
           Tcoarse(thecolumn,) and Tfine(i,). My guess is that it is related to 
           whether agg1 is greater than agg2 (assuming that val[0] =1)
            }
        }
     }
     */
   
     /********************************************************************/
     /* Matrix triple product for Pe where result is stuffed into MG     */
     /* hierarchy of the edge system.                                    */
     /*------------------------------------------------------------------*/
   
     if (Tfine->invec_leng != ml_nodes->Pmat[grid_level].outvec_leng)
     {
        printf("In ML_Gen_MGHierarchy_UsingReitzinger: Tmat and Pnodal\n"
               "\tdimensions on grid level %d do not agree:\n"
               "\tTmat->invec_leng = %d, Pnodal->outvec_leng = %d\n",
               grid_level,Tfine->invec_leng,
               ml_nodes->Pmat[grid_level].outvec_leng);
        exit(1);
     }
     if (ml_nodes->Pmat[grid_level].invec_leng != Tcoarse_trans->outvec_leng)
     {
        printf("In ML_Gen_MGHierarchy_UsingReitzinger:"
               " Pnodal and Tmat_trans\n"
               "\tdimensions on grid level %d do not agree:\n"
               "\tPnodal->invec_leng = %d, Tcoarse_trans->outvec_leng = %d\n",
               grid_level, ml_nodes->Pmat[grid_level].outvec_leng,
               Tcoarse_trans->outvec_leng);
        exit(1);
     }

     /* We could try to patch Tcoarse_trans to avoid the error check failing
        in ML_rap.  Right now, we do the check that's done in ML_rap.  If
        the check fails, the number of AMG levels is reduced by one and we 
        exit the main loop. */
  
     bail_flag = 0;
     N_input_vector = Tcoarse_trans->invec_leng;
     getrow_comm = Tcoarse_trans->getrow->pre_comm;
     if ( getrow_comm != NULL)
     {
        for (i = 0; i < getrow_comm->N_neighbors; i++) {
           for (j = 0; j < getrow_comm->neighbors[i].N_send; j++) {
              if (getrow_comm->neighbors[i].send_list[j] >= N_input_vector) {
                 bail_flag = 1;
                /*fix code would go here*/
              }
           }
        }
     }
     /* If check has failed, clean up current level & break from main loop. */
     ML_gsum_vec_int(&bail_flag,&j,1,ml_nodes->comm);
     if (bail_flag)
     {
        if (Tcoarse->comm->ML_mypid == 0 && 5 < ML_Get_PrintLevel()) {
           printf("(%d) Bailing from AMG hierarchy build on level %d"
                  " of levels %d down to %d ........\n",
                  Tcoarse->comm->ML_mypid,grid_level,fine_level,coarsest_level);
           fflush(stdout);
        }
        ML_Operator_Destroy(Tcoarse_trans);
        (*Tmat_trans_array)[grid_level] = NULL;
#ifdef NEW_T_PE
        ML_free(encoded_dir_node);
#endif
        /* Current level "grid_level" cannot be used b/c Tcoarse_trans
           is used in Hiptmair smoother generation, but RAP w/ Tcoarse_trans
           will give an error (the reason for the bail). Hence no "+1".*/
        Nlevels_nodal = fine_level - grid_level;
        coarsest_level = grid_level + 1;
        if (Tcoarse->comm->ML_mypid == 0 && 5 < ML_Get_PrintLevel()) {
           printf("(%d) In ML_Gen_MGHierarchy_UsingReitzinger, "
               "Nlevels_nodal = %d fine_level = %d  coarsest_level = %d\n",
              ml_nodes->comm->ML_mypid,Nlevels_nodal,fine_level,coarsest_level);
           fflush(stdout);
        }
        break; /* from main loop */

     }

     ML_rap(Tfine, &(ml_nodes->Pmat[grid_level]), Tcoarse_trans, 
        &(ml_edges->Pmat[grid_level]),ML_CSR_MATRIX);

     Pe = &(ml_edges->Pmat[grid_level]);

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
   
     Nnondirichlet = Pe->invec_leng - Npos_dirichlet - Nneg_dirichlet;
     for (i = 0; i < Pe->outvec_leng ; i++) {
#ifdef NEW_T_PE
       if (encoded_dir_node[i] > 0) {
	 for (j = csr_data->rowptr[i]; j < csr_data->rowptr[i+1] ; j++) {
	   if ( (csr_data->columns[j] < Nnondirichlet) ||
		(csr_data->columns[j] >= Nnondirichlet+Npos_dirichlet)) {
	     csr_data->values[j] = 0.;
	   }
	 }
       }
       if (encoded_dir_node[i] < 0) {
	 for (j = csr_data->rowptr[i]; j < csr_data->rowptr[i+1] ; j++) {
	   if (csr_data->columns[j] < Nnondirichlet+Npos_dirichlet) {
	     csr_data->values[j] = 0.;
	   }
	 }
       }
       if (encoded_dir_node[i] == 0) {
#endif
	 for (j = csr_data->rowptr[i]; j < csr_data->rowptr[i+1] ; j++) {
	   if (csr_data->values[j] == 2) csr_data->values[j] = 1;
	   else if (csr_data->values[j] == -2) csr_data->values[j] = -1;
	   else if (csr_data->values[j] == -1) csr_data->values[j] = 0;
	   else if (csr_data->values[j] ==  1) csr_data->values[j] = 0;
	   else if (csr_data->values[j] != 0.0)
	     {
	       printf("ML_Gen_MGHierarchy_UsingReitzinger:"
		      " Error in building Pe.   Found entry %e, expecting"
		      " either +/-1 or -2.\n",csr_data->values[j]);
	       fflush(stdout);
	     }
	 }
#ifdef NEW_T_PE
       }
#endif
     }
#ifdef NEW_T_PE
     ML_free(encoded_dir_node);
#endif
    
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

#ifdef LEASTSQ_SERIAL
     SPn_mat = ML_Operator_Create(Pn_coarse->comm);
     ML_Gen_SmoothPnodal(ml_nodes,grid_level+1, grid_level, 
			 &(ml_nodes->Amat[grid_level+1]),1.5,
			 SPn_mat);
     Pn_coarse->N_nonzeros = 10*Pn_coarse->outvec_leng;
     ml_leastsq_edge_interp(Pn_coarse, SPn_mat, Tfine, Tcoarse,Pe,
			    Pn_coarse->N_nonzeros*3);
     Tcoarse_trans = ML_Operator_Create(ml_edges->comm);
     ML_Operator_Transpose_byrow(Tcoarse, Tcoarse_trans);
     (*Tmat_trans_array)[grid_level] = Tcoarse_trans;
     Pn_coarse = SPn_mat;
#endif

/****************** Check the construction of Pe ***********************/
     vec = (double *) ML_allocate(sizeof(double)*(Pn_coarse->invec_leng+1+
						  Pe->outvec_leng));
     Pn_vec = (double *) ML_allocate(sizeof(double)*(Pn_coarse->outvec_leng+
						     Tcoarse->outvec_leng+1));
     Tfine_Pn_vec = (double *) ML_allocate(sizeof(double)*(Tfine->outvec_leng+1));
     ML_random_vec(vec, Pn_coarse->invec_leng, ml_edges->comm);

     ML_Operator_Apply(Pn_coarse, Pn_coarse->invec_leng, vec,
		       Pn_coarse->outvec_leng,Pn_vec);
     ML_Operator_Apply(Tfine, Tfine->invec_leng, Pn_vec,
		       Tfine->outvec_leng,Tfine_Pn_vec);
     ML_Operator_Apply(Tcoarse, Tcoarse->invec_leng, vec,
		       Tcoarse->outvec_leng,Pn_vec);
     ML_Operator_Apply(Pe, Pe->invec_leng, Pn_vec,
		       Pe->outvec_leng,vec);
     ML_free(Pn_vec);
     d1 = ML_gdot(Pe->outvec_leng, vec,vec, Pe->comm);
     d2 = ML_gdot(Pe->outvec_leng, Tfine_Pn_vec,Tfine_Pn_vec, Pe->comm);
     if (ml_edges->comm->ML_mypid == 0 && 8 < ML_Get_PrintLevel())  {
       if ( fabs(d1 - d2) > 1.0e-3)  
	 printf("ML_agg_reitzinger: Pe TH != Th Pn %e %e\n",d1,d2);
     }
     ML_free(vec); ML_free(Tfine_Pn_vec);

#ifdef LEASTSQ_SERIAL
     ML_Operator_Destroy(SPn_mat);
#endif


   
#ifdef postprocessscheck
      printf("Checking product Pe * e_i\n");
      yyy = (double *) ML_allocate( Pe->outvec_leng*sizeof(double));
      fido = (double *) ML_allocate( Pe->invec_leng*sizeof(double));
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
            /* This is for debugging a 2 process run.  The indices will
               be split between 2 processes;  here, split occurs at 120. */
            if ((Pe->comm->ML_mypid == 0) && (i < 120))
               fido[i] = 1;
            else if (Pe->comm->ML_mypid == 1 && i >= 120)
               fido[i-120] = 1;
         }
   
         ML_Operator_Apply(Pe, Pe->invec_leng, fido,
                           Pe->outvec_leng,yyy);
         dtemp = ML_gdot(Pe->outvec_leng, yyy, yyy, Pe->comm);
         printf("norm(P(:,%d)) = %e\n",i,dtemp); 
      }
   
     fflush(stdout);
     exit(1);
#endif /*postprocesscheck*/


      /***************************
      * Smooth edge prolongator. *
      ***************************/
     if (smooth_flag == ML_YES)
     {
        if (Tmat->comm->ML_mypid == 0)
        if (Tfine->comm->ML_mypid==0 && 3 < ML_Get_PrintLevel())
           printf("Smoothing edge prolongator...\n");
        ML_Aggregate_Set_Flag_SmoothExistingTentativeP(ag, ML_YES);
        /* default for smoothing factor is 4.0/3.0 */
        if (smooth_factor == ML_DDEFAULT)
           ML_Aggregate_Set_DampingFactor(ag, 4.0/3.0);
        else
           ML_Aggregate_Set_DampingFactor(ag, smooth_factor);
        ML_AGG_Gen_Prolongator(ml_edges,grid_level+1,grid_level,
                               (void *) &(ml_edges->Amat[grid_level+1]), ag);
        
        /* Weed out small values in Pe. */
        droptol = 1e-12;
        if (Tfine->comm->ML_mypid==0 && ag->print_flag < ML_Get_PrintLevel()) {
           printf("Dropping Pe entries with absolute value smaller than %e\n",
                  droptol);
        }
        Pe = &(ml_edges->Pmat[grid_level]);
        csr_data = (struct ML_CSR_MSRdata *) Pe->data;
        lower = csr_data->rowptr[0];
        nz_ptr = 0;
        for (i = 0; i < Pe->outvec_leng; i++) {
          for (j = lower; j < csr_data->rowptr[i+1]; j++) {
            if (ML_abs(csr_data->values[j]) > droptol) nz_ptr++;
          }
          lower = csr_data->rowptr[i+1];
          csr_data->rowptr[i+1] = nz_ptr;
        }
        nz_ptr = 0;
        for (i = 0; i < lower; i++) {
          if (ML_abs(csr_data->values[i]) > droptol) {
            csr_data->values[nz_ptr] = csr_data->values[i];
            csr_data->columns[nz_ptr] = csr_data->columns[i];
            nz_ptr++;
          }
        }

        Pe->N_nonzeros = nz_ptr;
     }

    if (ag->print_flag < ML_Get_PrintLevel()) {
        Pe = &(ml_edges->Pmat[grid_level]);
        nz_ptr = ML_Comm_GsumInt(ml_edges->comm, Pe->N_nonzeros);
        if (Tfine->comm->ML_mypid==0)
           printf("Pe: Total nonzeros = %d (Nrows = %d)\n", nz_ptr,
                  Pe->outvec_leng);
     }

     ML_Operator_Set_1Levels(&(ml_edges->Pmat[grid_level]),
                 &(ml_edges->SingleLevel[grid_level]), 
                 &(ml_edges->SingleLevel[grid_level+1]));
     ML_Gen_Restrictor_TransP(ml_edges, grid_level+1, grid_level);
     ML_Gen_AmatrixRAP(ml_edges, grid_level+1, grid_level);
     Nnz_allgrids += ml_edges->Amat[grid_level].N_nonzeros;

     if (ag->print_flag < ML_Get_PrintLevel()) {
        Pe = &(ml_edges->Amat[grid_level]);
        nz_ptr = ML_Comm_GsumInt(ml_edges->comm, Pe->N_nonzeros);
        if (Tfine->comm->ML_mypid==0)
           printf("Ke: Total nonzeros = %d (Nrows = %d)\n", nz_ptr,
                  Pe->outvec_leng);
     }

     Tfine = Tcoarse;
  } /* Main FOR loop: for grid_level = fine_level-1 ... */

  if (ag->print_flag < ML_Get_PrintLevel())
  {
    ML_gsum_vec_int(&Nnz_allgrids,&j,1,ml_nodes->comm);
    ML_gsum_vec_int(&Nnz_finegrid,&j,1,ml_nodes->comm);

    if (Tfine->comm->ML_mypid==0 )
    {
      if (Nnz_finegrid == 0) 
         printf("Number of nonzeros on finest grid not given!"
               " Complexity not computed!\n");
      else
         printf("Multilevel complexity is %e\n",
               ((double) Nnz_allgrids)/((double) Nnz_finegrid));
    }
  }

  ML_free(bindx);
  ML_free(val);

  if (created_ag_obj == 1) ML_Aggregate_Destroy(&ag);
  return(Nlevels_nodal);
}

/******************************************************************************/

int ML_MGHierarchy_ReitzingerDestroy(int finest_level, int coarsest_level,
                     ML_Operator ***Tmat_array, ML_Operator ***Tmat_trans_array)
{
    int level;

    if (*Tmat_array != NULL)
    {
       for (level = finest_level; level >= coarsest_level; level--)
       {
          ML_Operator_Destroy((*Tmat_array)[level]);
          (*Tmat_array)[level] = NULL;
       }
       ML_free(*Tmat_array);
       *Tmat_array = NULL;
    }

    if (*Tmat_trans_array != NULL)
    {
       for (level = finest_level; level >= coarsest_level; level--)
       {
          ML_Operator_Destroy((*Tmat_trans_array)[level]);
          (*Tmat_trans_array)[level] = NULL;
       }
       ML_free(*Tmat_trans_array);
       *Tmat_trans_array = NULL;
    }
    return 0;
}
/* ************************************************************************* */
/* generate smooth prolongator                                               */
/* ------------------------------------------------------------------------- */

int ML_Gen_SmoothPnodal(ML *ml,int level, int clevel, void *data,
			   double smoothP_damping_factor,
			   ML_Operator *SPn_mat)

{
  int         Ncoarse, Nfine;
   double      max_eigen = -1.;
   ML_Operator *Amat, *Pmatrix = NULL, *AGGsmoother = NULL;
   struct      ML_AGG_Matrix_Context widget;
   ML_Krylov   *kdata;

   if ( smoothP_damping_factor == 0.0 ) return 0;

   Amat     = (ML_Operator *) data;
   Nfine    = Amat->outvec_leng;
   Pmatrix =  &(ml->Pmat[clevel]);
   Ncoarse = Pmatrix->invec_leng;

   kdata = ML_Krylov_Create( ml->comm );
   ML_Krylov_Set_PrintFreq( kdata, 0 );
   ML_Krylov_Set_ComputeEigenvalues( kdata );
   ML_Krylov_Set_Amatrix(kdata, Amat);
   ML_Krylov_Solve(kdata, Nfine, NULL, NULL);
   max_eigen = ML_Krylov_Get_MaxEigenvalue(kdata);
   max_eigen = 2.;
   Amat->lambda_max = max_eigen; 
   Amat->lambda_min = kdata->ML_eigen_min; 
   ML_Krylov_Destroy( &kdata );
   if ( max_eigen <= 0.0 ) {
     printf("Gen_Prolongator warning : max eigen <= 0.0 \n");
     max_eigen = 1.0;
   }
   widget.near_bdry = NULL;
   widget.omega  = smoothP_damping_factor / max_eigen;

   widget.drop_tol = 0.;
   widget.Amat   = &(ml->Amat[level]);
   AGGsmoother = ML_Operator_Create(ml->comm);
   ML_Operator_Set_ApplyFuncData(AGGsmoother, widget.Amat->invec_leng,
				 widget.Amat->outvec_leng, ML_EXTERNAL,&widget,
				 widget.Amat->matvec->Nrows, NULL, 0);
   ML_Operator_Set_Getrow(AGGsmoother, ML_EXTERNAL,
                          widget.Amat->getrow->Nrows, 
                          ML_AGG_JacobiSmoother_Getrows);
   ML_CommInfoOP_Clone(&(AGGsmoother->getrow->pre_comm),
		       widget.Amat->getrow->pre_comm);

   ML_2matmult(AGGsmoother, Pmatrix, SPn_mat);
   ML_Operator_Destroy(AGGsmoother);

   return 0;
}
 
struct linked_list {
  struct linked_list *next;
  int duplicate_row;
};

extern int ml_comp_Pe_entries(int coef_cols[], double coef_values[], 
			      int coef_count, int leftagg, 
			      struct linked_list **Trecorder,
			      struct ML_CSR_MSRdata *Tfine,
			      int *Trowcount, int *Tnzcount,
			      struct ML_CSR_MSRdata *Tcoarse,
			      int *Pnzcount, int Pe_columns[], 
			      double Pe_values[]);

extern int ml_record_entry(struct linked_list **Trecorder,int lower,int upper, 
			int therow);
extern int ml_dup_entry(int node1, int node2, struct linked_list **Trecorder,
	int Tcols[], int Trowptr[], int *lower, int *upper, 
	int *duplicate_row);

extern int ml_clean_Trecorder(struct linked_list ***Trecorder ,int N);

extern int ml_leastsq_edge_interp(ML_Operator *Pn_mat, ML_Operator *SPn_mat, 
			   ML_Operator *Tfine_mat, ML_Operator *Tcoarse_mat, 
			   ML_Operator *Pe_mat, int);

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int ml_clean_Trecorder(struct linked_list ***Trecorder ,int N)
{
  struct linked_list *current, *tmp;
  int i;

  for (i = 0; i < N; i++) {
    current = (*Trecorder)[i];
    while (current != NULL) {
      tmp = current->next;
      ML_free(current);
      current = tmp;
    }
  }
  ML_free(*Trecorder);
  *Trecorder = NULL;
  return 1;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int ml_record_entry(struct linked_list **Trecorder, int lower, int upper, 
		 int therow)
{
  /* Record the edge (lower, upper) in the linked list table pointed to */
  /* by Trecorder. The matrix row index associated with this edge is    */
  /* recorded as well. It is assumed that this edge is not already in   */
  /* the table and that lower < upper.                                  */
  /* This routine is used in conjunction with 'ml_dup_entry' to see if     */
  /* edges have already been added to the table.                        */

  struct linked_list *prev_head;

  prev_head = Trecorder[lower];
  Trecorder[lower] = (struct linked_list *) ML_allocate(sizeof(struct linked_list));
  Trecorder[lower]->next = prev_head;
  Trecorder[lower]->duplicate_row  = therow;
  return 1;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int ml_dup_entry(int node1, int node2, struct linked_list **Trecorder,
	int Tcols[], int Trowptr[], int *lower, int *upper, int *duplicate_row)
{
  /* This routine is used to check if a pair (node1,node2) is already  */
  /* in the table pointed to by Trecorder (this table is normally      */
  /* constructed via ml_record_entry). On return,                         */
  /*     duplicate_row    Set to -1 if the pair is not in the table.   */
  /*                      Set to indicate matrix row associated with   */
  /*                      pair if (node1,node2) is in the table.       */

  struct linked_list *current;
  int tmp, therow, flag = 0;
  *lower = node1;
  *upper = node2;
  if (*lower > node2) {
    *lower = node2;
    *upper = node1;
  }
  if (*lower == -1) {
    *lower = *upper;
    *upper = -1;
  }
  *duplicate_row = -1;
  current = Trecorder[*lower];
  while (( current != NULL) && (*duplicate_row == -1)) {
    therow = current->duplicate_row;
    tmp = Trowptr[therow];
    if (Tcols[tmp] == *upper) { flag = 1; }
    else {
      if (tmp+1 < Trowptr[therow+1]) {
	if (Tcols[tmp+1]== *upper)
	  { flag = 1; }
      }
      else if (*upper == -1) { flag = 1; }
    }
    if (flag == 1) {
      *duplicate_row = current->duplicate_row;
    }
    current = current->next;
  }
  return 1;
}




int ml_leastsq_edge_interp(ML_Operator *Pn_mat, ML_Operator *SPn_mat, 
			   ML_Operator *Tfine_mat, ML_Operator *Tcoarse_mat, 
			   ML_Operator *Pe_mat, int Tcoarse_nz_bound)
{
  /* Given SPn, Tfine and the nodal aggregates (given by Pn), compute */
  /* Tcoarse and Pe such that the following relation holds:           */
  /*              Pe Tcoarse = Tfine SPn                              */
  /* Each row of Pe is computed by solving a small least squares whose*/
  /* solution is known analytically. This algorithm will be described */
  /* in more detail in an upcoming paper for copper mtn (3/02).       */

  int i, j, max_nz_per_row;
  struct linked_list **Trecorder;
  struct ML_CSR_MSRdata *Pn, *SPn, *Tfine, *Tcoarse, *Pe;
  double thesign;
  int left,leftagg, right, rightagg;
  int coef_count, coef_cols[100], jcoef_count, k;
  double alpha0, beta0, alpha_m, beta_m, coef_vals[100];
  int Trowcount = 0, Tnzcount = 0, Pnzcount = 0;
  int *SPn_columns, *SPn_rowptr, *Tfine_columns, *Pn_columns, *Pn_rowptr;
  int *Tfine_rowptr;
  double *SPn_values, *Tfine_values;

  /* pull out a bunch of pointers */

  Pn            = Pn_mat->data;
  Pn_columns    = Pn->columns;
  Pn_rowptr     = Pn->rowptr;
  Tfine         = Tfine_mat->data;
  Tfine_columns = Tfine->columns;
  Tfine_values  = Tfine->values;
  Tfine_rowptr  = Tfine->rowptr;
  SPn           = SPn_mat->data;
  SPn_columns   = SPn->columns;
  SPn_values    = SPn->values;
  SPn_rowptr    = SPn->rowptr;

  ML_Operator_Init(Tcoarse_mat,Pn_mat->comm);
  Tcoarse = (struct ML_CSR_MSRdata *) ML_allocate(sizeof(struct ML_CSR_MSRdata));
  Tcoarse_mat->data = Tcoarse;
  ML_Operator_Init(Pe_mat,Pn_mat->comm);
  Pe = (struct ML_CSR_MSRdata *) ML_allocate(sizeof(struct ML_CSR_MSRdata));
  Pe_mat->data = Pe;

  /* allocate space for Tcoarse */

  Tcoarse->columns = (int    *) ML_allocate(Tcoarse_nz_bound*sizeof(int));
  Tcoarse->values  = (double *) ML_allocate(Tcoarse_nz_bound*sizeof(double));
  Tcoarse->rowptr  = (int    *) ML_allocate(Tcoarse_nz_bound*sizeof(int));

  /* build header information for a table that will be used to record */
  /* coarse edges that have already been created (i.e. there already  */
  /* exists a row in Tcoarse).                                        */

  Trecorder = (struct linked_list **) ML_allocate((Pn_mat->outvec_leng+1)*
					     sizeof(struct linked_list *));
  for (i = 0; i < Pn_mat->outvec_leng; i++) Trecorder[i] = NULL;
  Tcoarse->rowptr[0] = 0;

  /* lets try to bound the amount of space needed in Pe */

  coef_count = 0;
  for (i=0; i < Tfine_mat->outvec_leng; i++) {
    if (Tfine_rowptr[i+1] - Tfine_rowptr[i] == 1) {
      left  = Tfine_columns[Tfine_rowptr[i]];
      coef_count += (SPn_rowptr[left+1] - SPn_rowptr[left]);
    }
    else if (Tfine_rowptr[i+1] - Tfine_rowptr[i] == 2) {
      right = Tfine_columns[Tfine_rowptr[i]];
      left  = Tfine_columns[Tfine_rowptr[i]+1];
      coef_count += (SPn_rowptr[left+1] - SPn_rowptr[left]);
      coef_count += (SPn_rowptr[right+1] - SPn_rowptr[right]);
    }
  }
  coef_count = 2*coef_count + 10;

  Pe->columns = (int    *) ML_allocate(coef_count*sizeof(int));
  Pe->values  = (double *) ML_allocate(coef_count*sizeof(double));
  Pe->rowptr  = (int    *) ML_allocate((Tfine_mat->outvec_leng+1)*sizeof(int));
  Pe_mat->outvec_leng = Tfine_mat->outvec_leng;
  Pe->rowptr[0] = 0;

  coef_count = 0;
  for (i=0; i < Tfine_mat->outvec_leng; i++) {
    /* special case when Tfine(i,:) has only one entry */
    if (Tfine_rowptr[i+1] - Tfine_rowptr[i] == 1) {
      if (Tfine_values[Tfine_rowptr[i]] == 1.) thesign = 1.;
      else thesign = -1.;
      left  = Tfine_columns[Tfine_rowptr[i]];
      coef_count = 0;
      for (j = SPn_rowptr[left]; j < SPn_rowptr[left+1]; j++) {
	coef_cols[coef_count  ] = SPn_columns[j];
	coef_vals[coef_count++] = thesign*SPn_values[j];
      }
      ml_comp_Pe_entries(coef_cols, coef_vals, coef_count, -1, Trecorder,
			 Tfine, &Trowcount, &Tnzcount, Tcoarse, 
			 &Pnzcount,Pe->columns,Pe->values);
   
    }

    /* normal case when Tfine(i,:) has only two entries */

    else if (Tfine_rowptr[i+1] - Tfine_rowptr[i] == 2) {

      /* determine the two fine grid points and the  */
      /* corresponding aggregates where they reside. */

      j = Tfine_rowptr[i];
      if (Tfine_values[j] == 1.) {
	right = Tfine_columns[j];
	left  = Tfine_columns[j+1];
      }
      else {
	left  = Tfine_columns[j];
	right = Tfine_columns[j+1];
      }
      leftagg = Pn_columns[Pn_rowptr[left]];
      rightagg = Pn_columns[Pn_rowptr[right]];

      coef_count = 0;

      if (leftagg == rightagg) {
	/* case 1 */

	/* copy alpha into coef */
	for (j = SPn_rowptr[left]; j < SPn_rowptr[left+1]; j++) {
	  if (SPn_columns[j] != leftagg) {
	    coef_cols[coef_count  ] = SPn_columns[j];
	    coef_vals[coef_count++] = -SPn_values[j];
	  }
	  else alpha0 = SPn_values[j];
	}
	jcoef_count = coef_count;
	/* copy beta into coef */
	for (j = SPn_rowptr[right]; j < SPn_rowptr[right+1]; j++) {
	  if (SPn_columns[j] != rightagg) {
	    for (k = 0; k < jcoef_count; k++) {
	      if (SPn_columns[j] == coef_cols[k]) break;
	    }
	    if (k != jcoef_count) {
	      coef_vals[k] += SPn_values[j];
	    }
	    else {
	      coef_cols[coef_count  ] = SPn_columns[j];
	      coef_vals[coef_count++] = SPn_values[j];
	    }
	  }
	  else beta0 = SPn_values[j];
	}

	ml_comp_Pe_entries(coef_cols, coef_vals, coef_count, leftagg, 
			   Trecorder, Tfine, &Trowcount, &Tnzcount, Tcoarse,
			   &Pnzcount,Pe->columns,Pe->values);
	coef_count = 0;
      }
      else {
	alpha_m = 0.; alpha0 = 0.;
	beta_m = 0.; beta0 = 0.;
	/* case 2 */

	/* copy alpha into coef */

	for (j = SPn_rowptr[left]; j < SPn_rowptr[left+1]; j++) {
	  if (SPn_columns[j] == leftagg)
	    alpha_m = SPn_values[j];
	  else if (SPn_columns[j] == rightagg)
	    alpha0 = SPn_values[j];
	  else {
	    coef_cols[coef_count  ] = SPn_columns[j];
	    coef_vals[coef_count++] = -SPn_values[j];
	  }
	}
	ml_comp_Pe_entries(coef_cols, coef_vals, coef_count, leftagg,
			   Trecorder,Tfine, &Trowcount, &Tnzcount, Tcoarse, 
			   &Pnzcount,Pe->columns,Pe->values);
	coef_count = 0;

      /* copy beta into coef */
	for (j = SPn_rowptr[right]; j < SPn_rowptr[right+1]; j++) {
	  if (SPn_columns[j] == leftagg)
	    beta_m = SPn_values[j];
	  else if (SPn_columns[j] == rightagg)
	    beta0 = SPn_values[j];
	  else {
	    coef_cols[coef_count  ] = SPn_columns[j];
	    coef_vals[coef_count++] = SPn_values[j];
	  }
	}
	ml_comp_Pe_entries(coef_cols, coef_vals, coef_count, rightagg,
			   Trecorder,Tfine, &Trowcount, &Tnzcount, Tcoarse, 
			   &Pnzcount,Pe->columns,Pe->values);
	coef_count = 1;
	coef_cols[0] = leftagg;
	coef_vals[0] = beta_m+alpha0-1.;
	ml_comp_Pe_entries(coef_cols, coef_vals, coef_count, rightagg,
			   Trecorder,Tfine, &Trowcount, &Tnzcount, Tcoarse,
			   &Pnzcount,Pe->columns,Pe->values);
      }


    }
    Pe->rowptr[i+1] = Pnzcount;
    if (Pe->rowptr[i+1] - Pe->rowptr[i] > max_nz_per_row)
      max_nz_per_row =  Pe->rowptr[i+1] - Pe->rowptr[i];
  }
  Pe_mat->outvec_leng = Tfine_mat->outvec_leng;

  Tcoarse->columns = (int    *) realloc(Tcoarse->columns, sizeof(int)*
					(Tcoarse->rowptr[Trowcount] +1));
  Tcoarse->values  = (double *) realloc(Tcoarse->values, sizeof(double)*
					(Tcoarse->rowptr[Trowcount] +1));
  Tcoarse->rowptr  = (int    *) realloc(Tcoarse->rowptr, sizeof(int)*
					(Trowcount+1));
  Pe->columns      = (int    *) realloc(Pe->columns, sizeof(int)*
					(Pnzcount +1));
  Pe->values       = (double *) realloc(Pe->values, sizeof(double)*
					(Pnzcount +1));


   ML_Operator_Set_ApplyFuncData(Pe_mat,Trowcount, Tfine_mat->outvec_leng,
				 ML_EMPTY,Pe,Tfine_mat->outvec_leng,NULL,0);
   ML_Operator_Set_ApplyFunc (Pe_mat, ML_INTERNAL, CSR_matvec);
   ML_Operator_Set_Getrow(Pe_mat, ML_EXTERNAL, Tfine_mat->outvec_leng, CSR_getrows);
   Pe_mat->getrow->Nrows = Tfine_mat->outvec_leng;
   Pe_mat->max_nz_per_row = max_nz_per_row;
   Pe_mat->N_nonzeros     = Pnzcount;

   ML_Operator_Set_ApplyFuncData(Tcoarse_mat,Pn_mat->invec_leng, Trowcount,
				 ML_EMPTY,Tcoarse,Trowcount,NULL,0);
   ML_Operator_Set_ApplyFunc (Tcoarse_mat, ML_INTERNAL, CSR_matvec);
   ML_Operator_Set_Getrow(Tcoarse_mat, ML_EXTERNAL, Trowcount, CSR_getrows);
   Tcoarse_mat->getrow->Nrows = Trowcount;
   Tcoarse_mat->max_nz_per_row = 2;
   Tcoarse_mat->N_nonzeros     = Tcoarse->rowptr[Trowcount];
   /*
  ML_Operator_Print(Pe_mat,"Pe");
  ML_Operator_Print(Tcoarse_mat,"Tcoarse");
   */

  ml_clean_Trecorder(&Trecorder,Pn_mat->outvec_leng);

  return 1;

}

int ml_comp_Pe_entries(int coef_cols[], double coef_values[], int coef_count,
		       int leftagg, struct linked_list **Trecorder,
		       struct ML_CSR_MSRdata *Tfine,
		       int *Trowcount, int *Tnzcount,
		       struct ML_CSR_MSRdata *Tcoarse, 
		       int *Pnzcount, int Pe_columns[], double Pe_values[]) {

  /* Add entries in Pe using interpolation coefficients coef_values[j]   */
  /* corresponding to the coarse grid edge between coarse grid node      */
  /* 'leftagg' and the coarse grid node 'coef_cols[j]. In the process,'  */
  /* add edges to Tcoarse (leftagg,coef_cols[j]) if they are not already */
  /* present.                                                            */

  int j, k, edge, lower, upper;
  int duplicate_row, *Tcoarse_columns, *Tcoarse_rowptr;
  double *Tcoarse_values;

  Tcoarse_columns = Tcoarse->columns;
  Tcoarse_values  = Tcoarse->values;
  Tcoarse_rowptr  = Tcoarse->rowptr;

  for (j = 0; j < coef_count; j++) {
    ml_dup_entry(leftagg, coef_cols[j], Trecorder,
	      Tcoarse_columns, Tcoarse_rowptr,&lower, &upper, 
	      &duplicate_row);
    if (duplicate_row == -1) {

      /* record the edge in Trecorder so that we can find */
      /* it later given two node entries.                */

      ml_record_entry(Trecorder, lower, upper, *Trowcount);
      edge = *Trowcount;
      Tcoarse_columns[*Tnzcount     ] = coef_cols[j];
      Tcoarse_values[(*Tnzcount)++ ] = 1;
      if (leftagg != -1) {
	Tcoarse_columns[*Tnzcount     ] = leftagg;
	Tcoarse_values[(*Tnzcount)++ ] = -1;
      }
      (*Trowcount)++;
      Tcoarse_rowptr[*Trowcount] = *Tnzcount;
    }
    else {
      edge = duplicate_row;
      k = Tcoarse_rowptr[duplicate_row];
      if (Tcoarse_columns[k] == leftagg) k++;
      if ( Tcoarse_values[k] == -1.) 
	coef_values[j] = -coef_values[j];
    }
    Pe_columns[ (*Pnzcount) ] = edge;
    Pe_values[ (*Pnzcount)++] = coef_values[j];
  }
  return 1;

}

