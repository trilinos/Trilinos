#include "ml_agg_reitzinger.h"
#include "ml_vampir.h"
#define ML_NEW_T_PE
#ifdef GREG
#define ML_NEW_T_PE
#endif

#if defined(ML_NEW_ENRICH)
double checkit(ML_Operator *A, double *v);
#endif

int  ML_Gen_MGHierarchy_UsingReitzinger(ML *ml_edges, ML* ml_nodes, 
                    int fine_level, int incr_or_decrease,
                    ML_Aggregate *ag, ML_Operator *Tmat,
                    ML_Operator *Tmat_trans,
                    ML_Operator ***Tmat_array,
                    ML_Operator ***Tmat_trans_array,
                    int smooth_flag, double smooth_factor)
{
  int coarsest_level, counter, Nghost, i, *Tcoarse_bindx = NULL;
  int *Tcoarse_rowptr, nz_ptr, row_length, j, *temp_bindx = NULL;
  double *Tcoarse_val = NULL, *temp_val = NULL;
  int allocated = 0, lower;
  double *pid_coarse_node, *pid_fine_node;
  int    *pid_fine_edge;
  ML_Operator *Kn_coarse, *Rn_coarse, *Tcoarse, *Pn_coarse;
  ML_Operator *Pe, *Tcoarse_trans, *Tfine;
  char filename[80];
  /* FILE *fid1; */
  FILE *fid2;
#ifdef LEASTSQ_SERIAL
  ML_Operator *SPn_mat;
#endif
#if defined(ML_ENRICH) || defined(ML_NEW_ENRICH)
  ML_Operator *TTtransPe;
  ML_Operator *newPe;
  double dtemp;
  double beta = 2./3.;
#endif
  struct ML_CSR_MSRdata *csr_data;
  int Nlevels_nodal, grid_level;
  int created_ag_obj = 0;
  double *vec, *Tcoarse_vec, *Pn_vec, *Tfine_Pn_vec;
  int i1, old_nzptr, i3, *index;
  int Npos_dirichlet = 0, Nneg_dirichlet = 0;
  int Nnondirichlet, Nnz_finegrid, Nnz_allgrids;
  double d1;
  double droptol;
  int *num_msgs, *tmp_msgs, *Nsnd_info, *Nrcv_info, type, partner;
  int *tmp1, *tmp2, Nlocal, Nneighbors_rcv, Nneighbors_snd, new_length;
  double *new_Tfine_Pn_vec;
  USR_REQ         *request;
  int t1,t2,t3;
#ifdef ML_NEW_T_PE
  int    *encoded_dir_node, Pn_ghost = 0;
  double *pos_coarse_dirichlet, *neg_coarse_dirichlet, *edge_type;
#endif

  int *glob_fine_edge_nums, *glob_fine_node_nums;
  int *glob_coarse_edge_nums, *glob_coarse_node_nums;
  ML_Operator *Ke;

  ML_CommInfoOP *getrow_comm; 
  int  N_input_vector;
  int  bail_flag;

  int max_matrix_size;
  int nrows, ncols;

  int nzctr;

  double t0 = GetClock();

#if defined(ML_NEW_ENRICH)
  double *v;
  double *w;
  double *randvec=NULL;
  double *Pv;
  double *WPv;
  ML_Operator *W, *V, *Tfine_trans, *tmpmat;
  double numer, denom;
  double alpha;
  double *diag;
#endif

  if (incr_or_decrease != ML_DECREASING)
    pr_error("Hiptmair: Only ML_DECREASING is supported\n");

#ifdef ML_VAMPIR
  ML_Vampir_Setup();
#endif


#ifdef ML_VAMPIR
  VT_begin(ml_vt_aggregating_state);
#endif

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
  ml_edges->ML_finest_level = fine_level;

  Nlevels_nodal = ML_Gen_MGHierarchy_UsingAggregation(ml_nodes, fine_level, 
                                            ML_DECREASING, ag);

  coarsest_level = fine_level - Nlevels_nodal + 1; 
  i = ml_nodes->Amat[coarsest_level].invec_leng;
  ML_gsum_scalar_int(&i,&j,ml_nodes->comm);
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

  if (ag->print_flag < ML_Get_PrintLevel()) {
    Pe = &(ml_edges->Amat[fine_level]);
    nz_ptr = ML_Comm_GsumInt(ml_edges->comm, Pe->N_nonzeros);
    i = Pe->outvec_leng;
    ML_gsum_scalar_int(&i,&j,ml_edges->comm);
    if (ml_edges->comm->ML_mypid==0)
      printf("(level %d) Ke: Global nonzeros = %d, global rows = %d\n",
             fine_level, nz_ptr,i);
  }

  /* this is buggy -- don't use it right now 
  fid2 = fopen("ML_write_matrix_now","r");
  if (fid2 != NULL) {
    i = fscanf(fid2,"%d",&max_matrix_size);
    // if empty file, print everything
    fclose(fid2);
    if (i==EOF || i==0) {
       sprintf(filename,"T_%d",fine_level);
       ML_Operator_Print_UsingGlobalOrdering(Tfine,filename);
       sprintf(filename,"Kn_%d",fine_level);
       ML_Operator_Print_UsingGlobalOrdering(&(ml_nodes->Amat[fine_level]),filename);
    }
    else {
      ML_Operator_GetGlobalDimensions(Tfine,&nrows,&ncols);
      if (nrows < i && ncols < i) {
        sprintf(filename,"T_%d",fine_level);
        ML_Operator_Print_UsingGlobalOrdering(Tfine,filename);
      }
      ML_Operator_GetGlobalDimensions(&(ml_nodes->Amat[fine_level]),&nrows,&ncols);
      if (nrows < i && ncols < i) {
        sprintf(filename,"Kn_%d",fine_level);
        ML_Operator_Print_UsingGlobalOrdering(&(ml_nodes->Amat[fine_level]),filename);
      }
    }
  }
  */




#ifdef ML_VAMPIR
  VT_end(ml_vt_aggregating_state);
#endif

  /********************************************************************/
  /*                 Build T on the coarse grid.                      */
  /*------------------------------------------------------------------*/

  /* main loop through grid hierarchy */
  for (grid_level = fine_level-1; grid_level >= coarsest_level ;grid_level--)
  {

#ifdef ML_VAMPIR
  VT_begin(ml_vt_building_coarse_T_state);
#endif
     Kn_coarse = &(ml_nodes->Amat[grid_level]);
     Rn_coarse = &(ml_nodes->Rmat[grid_level+1]);
     if (Kn_coarse->getrow->pre_comm != NULL)
        Nghost = Kn_coarse->getrow->pre_comm->total_rcv_length;
     else Nghost = 0;
   
     /* pid_coarse_node = node to processor lookup table
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
   
     pid_coarse_node = (double *) ML_allocate(sizeof(double)*(Nghost+
                            Kn_coarse->invec_leng+1));
     for (i = 0; i < Kn_coarse->invec_leng+Nghost; i++)
        pid_coarse_node[i] = (double) Kn_coarse->comm->ML_mypid + 1;
   
     if (Kn_coarse->getrow->pre_comm != NULL)
     ML_exchange_bdry(pid_coarse_node, Kn_coarse->getrow->pre_comm,
                      Kn_coarse->outvec_leng, Kn_coarse->comm,
                      ML_OVERWRITE,NULL);
   
#ifdef DEBUG_T_BUILD
     printf("\n\n%d: Kn_coarse->N_nonzeros = %d ",
            Kn_coarse->comm->ML_mypid, Kn_coarse->N_nonzeros);
     printf("Kn_coarse->invec_leng+Nghost = %d\n",
            Kn_coarse->invec_leng+Nghost);
     printf("Kn_coarse->invec_leng = %d\n\n",Kn_coarse->invec_leng);
#endif /* ifdef DEBUG_T_BUILD */

#ifdef ML_NEW_T_PE
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
	         temp_bindx[0] = temp_bindx[1];
	         temp_val[0] = temp_val[1];
	       }
	     }
       }
       else if ( (row_length == 1) && (temp_val[0] == 0.)) row_length--;

       if (row_length == 1) {
	     if      (temp_val[0] == 1.) encoded_dir_node[i] = (1 + temp_bindx[0]);
	     else if (temp_val[0] ==-1.) encoded_dir_node[i] = -(1 + temp_bindx[0]);
	     else printf("Warning uknown value T(%d,%d) = %e\n",
		     i,temp_bindx[0],temp_val[0]); 
       }
       else encoded_dir_node[i] = 0;
     }

     if (ml_nodes->Pmat[grid_level].getrow->pre_comm != NULL) {
       ML_CommInfoOP_Compute_TotalRcvLength(ml_nodes->Pmat[grid_level].getrow->pre_comm);
       Pn_ghost = ml_nodes->Pmat[grid_level].getrow->pre_comm->total_rcv_length;
     }


     vec = (double *) ML_allocate(sizeof(double)*(Rn_coarse->invec_leng +
						  Pn_ghost+1));

     pos_coarse_dirichlet = (double *) ML_allocate(sizeof(double)*(Rn_coarse->outvec_leng
						     +1));
     if (pos_coarse_dirichlet == NULL) {
        printf("\n\nML_Gen_MGHierarchy_UsingReitzinger: Not enough space allocated to check T.\n\n");
        exit(1);
     }
     for (i = 0; i < Rn_coarse->invec_leng + Pn_ghost; i++) vec[i] = 0.;
     for (i = 0; i < Tfine->outvec_leng; i++) {
       if(encoded_dir_node[i] > 0) vec[encoded_dir_node[i]-1] = 1.;
     }
     ML_transposed_exchange_bdry(vec,ml_nodes->Pmat[grid_level].getrow->pre_comm, 
				 -1, Rn_coarse->comm, ML_ADD);

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
        printf("\n\nML_Gen_MGHierarchy_UsingReitzinger: Not enough space allocated to check T.\n\n");
        exit(1);
     }
     for (i = 0; i < Rn_coarse->invec_leng + Pn_ghost; i++) vec[i] = 0.;
     for (i = 0; i < Tfine->outvec_leng; i++) {
       if(encoded_dir_node[i] < 0) vec[-encoded_dir_node[i]-1] = 1.;
     }
     ML_transposed_exchange_bdry(vec,ml_nodes->Pmat[grid_level].getrow->pre_comm, 
				 -1, Rn_coarse->comm, ML_ADD);

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
                          &temp_bindx,&temp_val,&row_length, 0);
        ML_az_sort(temp_bindx, row_length, NULL, NULL);
        nzctr += row_length;
   
        /* Step through unknowns bindx[j] connected to unknown i. */
        for (j = 0; j < row_length; j++)
        {
          {
             /* If nodes i and bindx[j] are owned by same processor ... */
             if (pid_coarse_node[i] == pid_coarse_node[temp_bindx[j]])
             {
                if (temp_bindx[j] > i)
                {
                   Tcoarse_bindx[nz_ptr]  =  temp_bindx[j];
                   Tcoarse_val[nz_ptr++]  =  1.;
                   Tcoarse_bindx[nz_ptr]  =  i;
                   Tcoarse_val[nz_ptr++]  = -1.;
                   Tcoarse_rowptr[counter+1] = nz_ptr;
                   counter++;
                }
             }
             /* If node i is owned by a smaller processor than
                node bindx[j] ... */
             else if (pid_coarse_node[i] < pid_coarse_node[temp_bindx[j]])
             {
                Tcoarse_bindx[nz_ptr]  =  temp_bindx[j];
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
        printf("\n\nML_Gen_MGHierarchy_UsingReitzinger: Not enough space allocated to build T.\n\n");
        exit(1);
     }
     Nnondirichlet = counter;
#ifdef ML_NEW_T_PE
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
       if (Kn_coarse->comm->ML_mypid == 0 && grid_level == 7) {
	  printf("%d (%d): ieqj = %d,\n", Kn_coarse->comm->ML_mypid,
		 grid_level, ieqj);
	  printf("Kn_coarse->N_nonzeros = %d, ", Kn_coarse->N_nonzeros);
	  printf("expected nnz (calc. from bindx) = %d,\n",nzctr);
	  printf("actual nnz(Tcoarse) = %d,\n"nz_ptr);
	  printf(" invec_leng = %d, nghost = %d\n",
		 Kn_coarse->invec_leng,  Nghost);
	  fflush(stdout);
       }
#endif /* ifdef DEBUG_T_BUILD */

   
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
        printf("\n\nML_Gen_MGHierarchy_UsingReitzinger: Not enough space allocated to check T.\n\n");
        exit(1);
     }

     ML_random_vec(vec, Tcoarse->invec_leng, ml_edges->comm);
     ML_Operator_Apply(Tcoarse, Tcoarse->invec_leng, vec,
		       Tcoarse->outvec_leng,Tcoarse_vec);

     for (i1 = 0; i1 < Tcoarse->outvec_leng; i1++) 
       if (Tcoarse_vec[i1] < 0) Tcoarse_vec[i1] = -Tcoarse_vec[i1];

     Pn_coarse->matvec->internal = CSR_ones_matvec; /* turn off the scaling */
                                                    /* in Pn for the test. */
     Pn_vec = (double *) ML_allocate(sizeof(double)*(Pn_coarse->outvec_leng
						     +1));
     if (Pn_vec == NULL) {
        printf("\n\nML_Gen_MGHierarchy_UsingReitzinger: Not enough space allocated to check T.\n\n");
        exit(1);
     }


     ML_Operator_Apply(Pn_coarse, Pn_coarse->invec_leng, vec,
		       Pn_coarse->outvec_leng,Pn_vec);

     ML_free(vec);
     Tfine_Pn_vec = (double *) ML_allocate(sizeof(double)*(Tfine->outvec_leng
							   +1));
     if (Tfine_Pn_vec == NULL) {
        printf("\n\nML_Gen_MGHierarchy_UsingReitzinger: Not enough space allocated to check T.\n\n");
        exit(1);
     }
     ML_Operator_Apply(Tfine, Tfine->invec_leng, Pn_vec,
		       Tfine->outvec_leng,Tfine_Pn_vec);
     ML_free(Pn_vec);
     Pn_coarse->matvec->internal = CSR_matvec;
     for (i1 = 0; i1 < Tfine->outvec_leng; i1++) 
       if (Tfine_Pn_vec[i1] < 0) Tfine_Pn_vec[i1] = -Tfine_Pn_vec[i1];

     /* This code is to fix a bug that can occur in parallel   */
     /* Basically, the problem is that the fine grid edge that */
     /* we are trying to match with a coarse grid edge may     */
     /* reside on another processor. The way we will get around*/
     /* this is to figure out what processor will need each    */
     /* fine grid node and ship it to that processor. This is  */
     /* done as follows:                                       */
     /*   1) fill coarse nodes with proc id                    */
     /*   2) interpolate coarse nodes (no scaling in interp)   */
     /*   3) Use Tfine's pre_comm is communicate the proc ids  */
     /*   4) Go through each row of Tfine and record the lower */
     /*      node proc id for each edge. This is the processor */
     /*      that will need the information when looking for   */
     /*      bogus edges.                                      */
     /*   5) Send this information to the relevant processor   */
     /*      so that he can append it to his list.             */

     Pn_coarse = &(ml_nodes->Pmat[grid_level]);
     i = 0;
     if (Tfine->getrow->pre_comm != NULL) {
       ML_CommInfoOP_Compute_TotalRcvLength(Tfine->getrow->pre_comm);
       i = Tfine->getrow->pre_comm->total_rcv_length;
     }
     pid_fine_node = (double *) ML_allocate(sizeof(double)*(1+i+
					      Tfine->invec_leng));
     if (pid_fine_node == NULL) {
        printf("\n\nML_Gen_MGHierarchy_UsingReitzinger: Not enough space allocated to ship fine edges in T.\n\n");
        exit(1);
     }

     Pn_coarse->matvec->internal = CSR_ones_matvec; /* turn off the scaling */
                                                    /* in Pn for the test. */
     ML_Operator_Apply(Pn_coarse, Pn_coarse->invec_leng, pid_coarse_node,
		       Pn_coarse->outvec_leng, pid_fine_node);
     Pn_coarse->matvec->internal = CSR_matvec;      /* turn on the scaling */
     ML_free(pid_coarse_node);

     if (Tfine->getrow->pre_comm != NULL) {
       ML_exchange_bdry(pid_fine_node, Tfine->getrow->pre_comm,
			Tfine->invec_leng, Tfine->comm,
			ML_OVERWRITE,NULL);
     }
     pid_fine_edge = (int *) ML_allocate(sizeof(int)*(1+
					      Tfine->outvec_leng));
     if (pid_fine_edge == NULL) {
        printf("\n\nML_Gen_MGHierarchy_UsingReitzinger: Not enough space allocated to ship fine edges in T.\n\n");
        exit(1);
     }
     for (i = 0 ; i < Tfine->getrow->Nrows; i++) {
       pid_fine_edge[i] = Tfine->comm->ML_nprocs + 10;
       ML_get_matrix_row(Tfine, 1, &i, &allocated, &temp_bindx, &temp_val,
			 &row_length, 0);
       if ((row_length > 0) && (pid_fine_node[temp_bindx[0]] != 0.)) 
	 pid_fine_edge[i] = (int) pid_fine_node[temp_bindx[0]];
       if ((row_length > 1) && (pid_fine_node[temp_bindx[1]] != 0.)) {
	 if ( pid_fine_node[temp_bindx[1]] < pid_fine_node[temp_bindx[0]])
	   pid_fine_edge[i] = (int) pid_fine_node[temp_bindx[1]];
       }
       pid_fine_edge[i]--;
     }
     ML_free(pid_fine_node);
     ML_free(temp_bindx);
     ML_free(temp_val);
     allocated = 0;

     num_msgs = (int *) ML_allocate(sizeof(int)*(Tfine->comm->ML_nprocs+1));
     tmp_msgs = (int *) ML_allocate(sizeof(int)*(Tfine->comm->ML_nprocs+1));
     for (i = 0; i < Tfine->comm->ML_nprocs; i++) num_msgs[i] = 0;
     for (i = 0; i < Tfine->comm->ML_nprocs; i++) tmp_msgs[i] = 0;
     Nneighbors_snd = 0;
     for (i = 0 ; i < Tfine->getrow->Nrows; i++) {
       if ( (pid_fine_edge[i] < Tfine->comm->ML_nprocs) &&
	    (pid_fine_edge[i] != Tfine->comm->ML_mypid)) {
	 tmp_msgs[pid_fine_edge[i]]++;
	 if (num_msgs[pid_fine_edge[i]] == 0) {
	   num_msgs[pid_fine_edge[i]] = 1;
	   Nneighbors_snd++;
	 }
       }
     }
     /* Set up a send list and a receive list. In each list a pair of    */
     /* numbers is stored (neighbor, length) indicating the neighbor and */
     /* and the message length. For example, the ith neighbor and length */
     /* would be in Nsnd_info[2*i] and Nsnd_info[2*i+1] respectively.    */

     Nsnd_info = (int *) ML_allocate(sizeof(int)*(2*Nneighbors_snd+1));
     Nneighbors_snd = 0;
     for (i = 0; i < Tfine->comm->ML_nprocs; i++) {
       if (tmp_msgs[i] > 0) {
	 Nsnd_info[Nneighbors_snd++] = i;
	 Nsnd_info[Nneighbors_snd++] = tmp_msgs[i];
       }
     }
     Nneighbors_snd /= 2;

     /* Figure out how many neighbors we will recieve from */     
     ML_gsum_vec_int(&num_msgs,&tmp_msgs,Tfine->comm->ML_nprocs,Tfine->comm);
     Nneighbors_rcv = num_msgs[Tfine->comm->ML_mypid];
     ML_free(tmp_msgs);
     ML_free(num_msgs);

     /* Build a receive list similar to the send list */

     Nrcv_info = (int *) ML_allocate(sizeof(int)*(2*Nneighbors_rcv+1));
     if ( Nneighbors_rcv > 0 )  {
       request = (USR_REQ  *)  ML_allocate(Nneighbors_rcv*sizeof(USR_REQ ));
     } else request = NULL;


     type = 1235; 
     /* post receives */
     for (i = 0 ; i < Nneighbors_rcv ; i++) {
       partner = -1;
       Tfine->comm->USR_irecvbytes((void *) &(Nrcv_info[2*i]), 2*sizeof(int),
			    &partner, &type,Tfine->comm->USR_comm,request+i);
     }  
     /* send messages */
     for (i = 0 ; i < Nneighbors_snd ; i++) {
       partner = Nsnd_info[2*i]; 
       Nsnd_info[2*i] = Tfine->comm->ML_mypid;
       Tfine->comm->USR_sendbytes((void *) &(Nsnd_info[2*i]), 2*sizeof(int), 
			   partner, type, Tfine->comm->USR_comm);
       Nsnd_info[2*i] = partner;
     }
     /* wait */
     for (i = 0; i <  Nneighbors_rcv ; i++) {
       partner = -1;
       Tfine->comm->USR_cheapwaitbytes((void *) &(Nrcv_info[2*i]), 2*sizeof(int), 
			   &partner, &type, Tfine->comm->USR_comm, request+i);
     }

     /* sort `Nrcv_info' by proc id's */

     tmp1 = (int *) ML_allocate(sizeof(int)*(Nneighbors_rcv+1));
     tmp2 = (int *) ML_allocate(sizeof(int)*(Nneighbors_rcv+1));
     for (i = 0; i < Nneighbors_rcv; i++) {
       tmp1[i] = Nrcv_info[2*i];
       tmp2[i] = Nrcv_info[2*i+1];
     }
     ML_az_sort(tmp1, Nneighbors_rcv, tmp2, NULL);
     for (i = 0; i < Nneighbors_rcv; i++) {
       Nrcv_info[2*i]   = tmp1[i];
       Nrcv_info[2*i+1] = tmp2[i];
     }
     ML_free(tmp2); ML_free(tmp1);

     /* Reorder Tfine_Pn_vec so that all values needed locally */
     /* come first, followed by all values to be sent to the   */
     /* lowest numbered processor to highest numbered processor*/

     for (i = 0 ; i < Tfine->getrow->Nrows; i++) {
       if (pid_fine_edge[i] > Tfine->comm->ML_nprocs) 
	 pid_fine_edge[i] = -1;
      if (pid_fine_edge[i] == Tfine->comm->ML_mypid) 
	 pid_fine_edge[i] = -1;
     }
     ML_az_sort(pid_fine_edge, Tfine->outvec_leng, NULL, Tfine_Pn_vec);

     /* Organize a new Tfine_Pn_vec which will contain the locally */
     /* needed components as well as those received from other processors */

     Nlocal = Tfine->outvec_leng;
     for (i = 0; i <  Tfine->outvec_leng; i++) {
       if ( pid_fine_edge[i] != -1) break;
     }
     if (i < Tfine->outvec_leng ) Nlocal = i;
     ML_free(pid_fine_edge);

     i1 = 0;
     for (i = 0; i <  Nneighbors_rcv ; i++) 
       i1 += Nrcv_info[2*i+1];

     new_Tfine_Pn_vec = (double *) ML_allocate(sizeof(double)*(Nlocal+1+i1));
     for (i = 0; i < Nlocal; i++) new_Tfine_Pn_vec[i] = Tfine_Pn_vec[i];

     
     /* pack up data and communicate data */

     type = 1275; 
     /* post receives */
     i1 = Nlocal;
     for (i = 0 ; i < Nneighbors_rcv ; i++) {
       Tfine->comm->USR_irecvbytes((void *) &(new_Tfine_Pn_vec[i1]), 
				   Nrcv_info[2*i+1]*sizeof(double),
				   &(Nrcv_info[2*i]), &type,
				   Tfine->comm->USR_comm, request+i);
       i1 += Nrcv_info[2*i+1];
     }
     i1 = Nlocal;
     for (i = 0 ; i < Nneighbors_snd ; i++) {
       Tfine->comm->USR_sendbytes((void *) &(Tfine_Pn_vec[i1]), 
				  Nsnd_info[2*i+1]*sizeof(double), 
				  Nsnd_info[2*i],type,Tfine->comm->USR_comm);
       i1 += Nsnd_info[2*i+1];
     }

     i1 = Nlocal;
     for (i = 0; i <  Nneighbors_rcv ; i++) {
       Tfine->comm->USR_cheapwaitbytes((void *) &(new_Tfine_Pn_vec[i1]), 
				  Nrcv_info[2*i+1]*sizeof(double),
				  &(Nrcv_info[2*i]), &type,
				  Tfine->comm->USR_comm, request+i);
       i1 += Nrcv_info[2*i+1];
     }
     new_length = i1;
     ML_free(Nrcv_info);
     ML_free(Nsnd_info);
     ML_free(request);
     ML_free(Tfine_Pn_vec);
     Tfine_Pn_vec = new_Tfine_Pn_vec;


     /* Now sort the two vectors and then compare them.    */
     /* Mark Tcoarse_vec[i] with a '1.' if its corresponding */
     /* edge can be removed.                                 */

     index = (int *) ML_allocate((Tcoarse->outvec_leng+1)*sizeof(int));
     for (i1 = 0; i1 < Tcoarse->outvec_leng; i1++) index[i1] = i1;
                            /* use index to keep track of which edge */
                            /* is which in the sorted Tcoarse_vec    */
     ML_az_dsort2(Tcoarse_vec,Tcoarse->outvec_leng, index);
     ML_az_dsort(Tfine_Pn_vec,new_length);

     i3 = 0;
     for (i1 = 0; i1 < Tcoarse->outvec_leng; i1++) {
       while ( (i3 != new_length) && 
	       (Tcoarse_vec[i1] > Tfine_Pn_vec[i3])) {
	 i3++;
       }
       if ((i3 < new_length)&&(Tcoarse_vec[i1] == Tfine_Pn_vec[i3])) {
	 Tcoarse_vec[i1] = 0.;
       }
       else Tcoarse_vec[i1] = 1.;
     }
     ML_free(Tfine_Pn_vec);
     ML_az_sort(index, Tcoarse->outvec_leng, NULL, Tcoarse_vec);
     ML_free(index);

     t1 = Nnondirichlet;
     t2 = Npos_dirichlet;
     t3 = Nneg_dirichlet;
     for (i1 = 0; i1 < Tcoarse->outvec_leng; i1++) {
       if (Tcoarse_vec[i1] != 0.) {
	 if (i1 < Nnondirichlet) t1--; 
	 else if (i1 < Nnondirichlet + Npos_dirichlet) 
	   t2--;
	 else if (i1 < Nnondirichlet + Npos_dirichlet + Nneg_dirichlet) 
	   t3--;
       }
     }
     Nnondirichlet = t1;
     Npos_dirichlet = t2;
     Nneg_dirichlet = t3;
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
   
     /* Check that both dimensions of T are strictly greater than 0. 
        If not, clean up & break from main loop. */
     i = Tcoarse->outvec_leng;
     ML_gsum_scalar_int(&i,&j,ml_nodes->comm);
     if (i==0)
     {
        if (Tcoarse->comm->ML_mypid == 0 && 5 < ML_Get_PrintLevel()) {
           printf("(%d) Bailing from AMG hierarchy build on level %d of levels %d down to %d because Tcoarse has zero rows....\n",
                 ml_edges->comm->ML_mypid,grid_level,fine_level,coarsest_level);
           fflush(stdout);
        }
        /*Don't delete Tcoarse because it is necessary for eigenvalue estimate
          in Hiptmair setup.*/
#ifdef ML_NEW_T_PE
        ML_free(encoded_dir_node);
#endif
        /* Current level "grid_level" cannot be used b/c Tcoarse_trans
           would be used in Hiptmair smoother generation, but Tcoarse_trans
           hasn't been calculated. Hence no "+1".*/
        Nlevels_nodal = fine_level - grid_level;
        coarsest_level = grid_level + 1;
        if (Tcoarse->comm->ML_mypid == 0 && 5 < ML_Get_PrintLevel()) {
           printf("(%d) In ML_Gen_MGHierarchy_UsingReitzinger, Nlevels_nodal = %di, fine_level = %d, coarsest_level = %d\n",
              ml_nodes->comm->ML_mypid,Nlevels_nodal,fine_level,coarsest_level);
           fflush(stdout);
        }
        break; /* from main loop */
     }

     /********************************************************************/
     /* Create Tcoarse_trans.                                            */
     /*------------------------------------------------------------------*/

     Tcoarse_trans = ML_Operator_Create(ml_edges->comm);
     ML_Operator_Transpose_byrow(Tcoarse, Tcoarse_trans);
     (*Tmat_trans_array)[grid_level] = Tcoarse_trans;

#ifdef ML_VAMPIR
  VT_end(ml_vt_building_coarse_T_state);
#endif

   
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

#ifdef ML_VAMPIR
  VT_begin(ml_vt_build_Pe_state);
#endif
   
     if (Tfine->invec_leng != ml_nodes->Pmat[grid_level].outvec_leng)
     {
       printf("In ML_Gen_MGHierarchy_UsingReitzinger: Tmat and Pnodal\n");
       printf("\tdimensions on grid level %d do not agree:\n", grid_level);
       printf("\tTmat->invec_leng = %d, Pnodal->outvec_leng = %d\n",
               Tfine->invec_leng, ml_nodes->Pmat[grid_level].outvec_leng);
        exit(1);
     }
     if (ml_nodes->Pmat[grid_level].invec_leng != Tcoarse_trans->outvec_leng)
     {
       printf("In ML_Gen_MGHierarchy_UsingReitzinger:");
       printf(" Pnodal and Tmat_trans\n");
       printf("\tdimensions on grid level %d do not agree:\n", grid_level);
       printf("\tPnodal->invec_leng = %d, Tcoarse_trans->outvec_leng = %d\n",
	      ml_nodes->Pmat[grid_level].outvec_leng,
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
     ML_gsum_scalar_int(&bail_flag,&j,ml_nodes->comm);
     if (bail_flag)
     {
        if (Tcoarse->comm->ML_mypid == 0 && 5 < ML_Get_PrintLevel()) {
           printf("(%d) Bailing from AMG hierarchy build on level %d of levels %d down to %d ........\n",
                  Tcoarse->comm->ML_mypid,grid_level,fine_level,coarsest_level);
           fflush(stdout);
        }
        ML_Operator_Destroy(&Tcoarse_trans);
        (*Tmat_trans_array)[grid_level] = NULL;
#ifdef ML_NEW_T_PE
        ML_free(encoded_dir_node);
#endif
        /* Current level "grid_level" cannot be used b/c Tcoarse_trans
           is used in Hiptmair smoother generation, but RAP w/ Tcoarse_trans
           will give an error (the reason for the bail). Hence no "+1".*/
        Nlevels_nodal = fine_level - grid_level;
        coarsest_level = grid_level + 1;
        if (Tcoarse->comm->ML_mypid == 0 && 5 < ML_Get_PrintLevel()) {
           printf("(%d) In ML_Gen_MGHierarchy_UsingReitzinger, Nlevels_nodal = %d fine_level = %d  coarsest_level = %d\n",
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
   
#ifdef ML_NEW_T_PE
     i = 0;
     if (Pe->getrow->pre_comm != NULL) {
       ML_CommInfoOP_Compute_TotalRcvLength(Pe->getrow->pre_comm);
       i = Pe->getrow->pre_comm->total_rcv_length;
     }
     edge_type = (double *) ML_allocate(sizeof(double)*(1+i+Pe->invec_leng));
     if (edge_type == NULL) {
        printf("\n\nML_Gen_MGHierarchy_UsingReitzinger: Not enough space allocated to weed out Pe.\n\n");
        exit(1);
     }
     Nnondirichlet = Pe->invec_leng - Npos_dirichlet - Nneg_dirichlet;

     for (i = 0; i < Nnondirichlet; i++) edge_type[i] = 0.;
     for (i = 0; i < Npos_dirichlet; i++) edge_type[i+Nnondirichlet] = 1.;
     for (i = 0; i < Nneg_dirichlet; i++) edge_type[i+Nnondirichlet+
                                                    Npos_dirichlet] = -1.;
     if (Pe->getrow->pre_comm != NULL) {
       ML_exchange_bdry(edge_type, Pe->getrow->pre_comm,
			Pe->invec_leng, Pe->comm,
			ML_OVERWRITE,NULL);
     }
#endif
     for (i = 0; i < Pe->outvec_leng ; i++) {
#ifdef ML_NEW_T_PE
       if (encoded_dir_node[i] > 0) {
	 for (j = csr_data->rowptr[i]; j < csr_data->rowptr[i+1] ; j++) {
           if (edge_type[csr_data->columns[j]] != 1.0 ) {
	     csr_data->values[j] = 0.;
	   }
	 }
       }
       if (encoded_dir_node[i] < 0) {
	 for (j = csr_data->rowptr[i]; j < csr_data->rowptr[i+1] ; j++) {
           if (edge_type[csr_data->columns[j]] != -1.0 ) {
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
	       printf("ML_Gen_MGHierarchy_UsingReitzinger:");
	       printf(" Error in building Pe.   Found entry %e, expecting either +/-1 or -2.\n",csr_data->values[j]);
	       fflush(stdout);
	     }
	 }
#ifdef ML_NEW_T_PE
       }
#endif
     }
#ifdef ML_NEW_T_PE
     ML_free(encoded_dir_node);
     ML_free(edge_type);
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
     if (ml_nodes->comm->ML_mypid == 0)
        { printf("\n\nDoing LS prolongator\n\n");fflush(stdout);}
     SPn_mat = ML_Operator_Create(Pn_coarse->comm);
     ML_Gen_SmoothPnodal(ml_nodes,grid_level+1, grid_level, 
			 &(ml_nodes->Amat[grid_level+1]),1.5,
			 SPn_mat);

     /*sprintf(filename,"Spn%d",grid_level-1);
     ML_Operator_Print(SPn_mat,filename);*/

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
     d1 = sqrt(ML_gdot(Pn_coarse->invec_leng, vec, vec, ml_edges->comm));
     for (i = 0; i < Pn_coarse->invec_leng; i++) 
       vec[i] /= d1;

     /*sprintf(filename,"randvec%d",grid_level-1);
     fid = fopen(filename,"w");
     for (i=0; i<Pn_coarse->invec_leng; i++)
        fprintf(fid,"%d    %20.15e\n",i,vec[i]);
     fclose(fid);*/

     ML_Operator_Apply(Pn_coarse, Pn_coarse->invec_leng, vec,
		       Pn_coarse->outvec_leng,Pn_vec);
     ML_Operator_Apply(Tfine, Tfine->invec_leng, Pn_vec,
		       Tfine->outvec_leng,Tfine_Pn_vec);
     ML_Operator_Apply(Tcoarse, Tcoarse->invec_leng, vec,
		       Tcoarse->outvec_leng,Pn_vec);
     ML_Operator_Apply(Pe, Pe->invec_leng, Pn_vec,
		       Pe->outvec_leng,vec);
     ML_free(Pn_vec);

     /* this is buggy -- don't use it right now
     fid2 = fopen("ML_write_matrix_now","r");
     if (fid2 != NULL) {
       i = fscanf(fid2,"%d",&max_matrix_size);
       // if empty file, print everything 
       fclose(fid2);
       if (i==EOF || i==0) {
         sprintf(filename,"Kn_%d",grid_level);
         ML_Operator_Print_UsingGlobalOrdering(Kn_coarse,filename);
         sprintf(filename,"Pn_%d",grid_level);
         ML_Operator_Print_UsingGlobalOrdering(Pn_coarse,filename);
         sprintf(filename,"T_%d",grid_level);
         ML_Operator_Print_UsingGlobalOrdering(Tcoarse,filename);
         sprintf(filename,"Pe_%d",grid_level);
         ML_Operator_Print_UsingGlobalOrdering(Pe,filename);
       }
       else {
         ML_Operator_GetGlobalDimensions(Kn_coarse,&nrows,&ncols);
         if (nrows < i && ncols < i) {
           sprintf(filename,"Kn_%d",grid_level);
           ML_Operator_Print_UsingGlobalOrdering(Kn_coarse,filename);
         }
         ML_Operator_GetGlobalDimensions(Tcoarse,&nrows,&ncols);
         if (nrows < i && ncols < i) {
           sprintf(filename,"T_%d",grid_level);
           ML_Operator_Print_UsingGlobalOrdering(Tcoarse,filename);
         }
         ML_Operator_GetGlobalDimensions(Pn_coarse,&nrows,&ncols);
         if (nrows < i && ncols < i) {
           sprintf(filename,"Pn_%d",grid_level);
           ML_Operator_Print_UsingGlobalOrdering(Pn_coarse,filename);
         }
         ML_Operator_GetGlobalDimensions(Pe,&nrows,&ncols);
         if (nrows < i && ncols < i) {
           sprintf(filename,"Pe_%d",grid_level);
           ML_Operator_Print_UsingGlobalOrdering(Pe,filename);
         }
       }
     }
     */

     for (i = 0; i < Pe->outvec_leng ; i++) 
       vec[i] = vec[i] - Tfine_Pn_vec[i];

     if (8 < ML_Get_PrintLevel())  {
       d1 = sqrt(ML_gdot(Tfine->outvec_leng, Tfine_Pn_vec,Tfine_Pn_vec, Pe->comm));
     if (ml_edges->comm->ML_mypid == 0 )
       printf("\n\nML_agg_reitzinger:  ||Th Pn v|| = %15.10e\n\n",d1);
     }
     d1 = sqrt(ML_gdot(Pe->outvec_leng, vec,vec, Pe->comm));

     if (8 < ML_Get_PrintLevel())  {
       if (ml_edges->comm->ML_mypid == 0 )
         printf("ML_agg_reitzinger: %e\n",d1);
       if ( fabs(d1) > 1.0e-4)
       {
         if (ml_edges->comm->ML_mypid == 0 )
           printf("ML_agg_reitzinger: Pe TH != Th Pn %e (level %d)\n",
                  d1,grid_level);
         if (14 < ML_Get_PrintLevel() )  {

            Ke = ml_edges->Amat+grid_level+1;
            ML_build_global_numbering(Ke, Ke->comm, &glob_fine_edge_nums);
            Ke = ml_nodes->Amat+grid_level+1;
            ML_build_global_numbering(Ke,Ke->comm, &glob_fine_node_nums);
            Ke = ml_nodes->Amat+grid_level;
            ML_build_global_numbering(Ke,Ke->comm, &glob_coarse_node_nums);
            /* ml_edges->Amat+grid_level doesn't exist yet.. */
            ML_build_global_numbering(Tcoarse,Tcoarse->comm,
                                      &glob_coarse_edge_nums);

            ML_Operator_Print_UsingGlobalOrdering(Pn_coarse,"Pn_debug",
                              glob_fine_node_nums,glob_coarse_node_nums);
            ML_Operator_Print_UsingGlobalOrdering(Tfine,"Tfine_debug",
                              glob_fine_edge_nums,glob_fine_node_nums);
            ML_Operator_Print_UsingGlobalOrdering(Tcoarse,"Tcoarse_debug",
                              glob_coarse_edge_nums,glob_coarse_node_nums);
            ML_Operator_Print_UsingGlobalOrdering(Pe,"Pe_debug",
                              glob_fine_edge_nums,glob_coarse_edge_nums);

            ML_free(glob_fine_edge_nums);
            ML_free(glob_fine_node_nums);
            ML_free(glob_coarse_node_nums);
            ML_free(glob_coarse_edge_nums);
         }
         for (i = 0; i < Pe->outvec_leng; i++) {
           /* change this tolerance if you actually want */
           /* to see the individual components.          */
           if (fabs(vec[i]) > 1.0) 
             fprintf(stderr,"%d: ===> %d is %20.13e vs %20.13e\n",
                Pe->comm->ML_mypid,i,vec[i] + Tfine_Pn_vec[i],Tfine_Pn_vec[i]);
         }
       }
     }

     ML_free(vec); ML_free(Tfine_Pn_vec);

#ifdef LEASTSQ_SERIAL
     ML_Operator_Destroy(&SPn_mat);
#endif

#ifdef ML_VAMPIR
  VT_end(ml_vt_build_Pe_state);
#endif
   
      /***************************
      * Smooth edge prolongator. *
      ***************************/

#ifdef ML_VAMPIR
  VT_begin(ml_vt_smooth_Pe_state);
#endif
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
        ML_AGG_Gen_Prolongator(ml_edges,grid_level+1,grid_level,(void*) ag);
        
        /* Weed out small values in Pe. */
        droptol = 1e-4;
#ifdef GREG
        droptol = 1e-24;
#endif
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
            if (ML_abs(csr_data->values[j]) > droptol) {
              csr_data->values[nz_ptr] = csr_data->values[j];
              csr_data->columns[nz_ptr] = csr_data->columns[j];
              nz_ptr++;
            }
          }
          lower = csr_data->rowptr[i+1];
          csr_data->rowptr[i+1] = nz_ptr;
        }

        Pe->N_nonzeros = nz_ptr;
     }

#ifdef ML_VAMPIR
  VT_end(ml_vt_smooth_Pe_state);
#endif
#if defined(ML_ENRICH) || defined(ML_NEW_ENRICH)

#ifdef ML_VAMPIR
  VT_begin(ml_vt_enrich_Pe_state);
#endif

     /* ********************************************************************* */
     /* Generate enriched smoothed prolongator                                */
     /* The idea is that we want to smooth the enriched prolongator:          */
     /*  =  (I + alpha inv(diag(S))*S) * (I + beta inv(diag(T*T'))*(T*T'))    */
     /*  =  (I + alpha inv(diag(S))*S + beta inv(diag(T*T'))*(T*T'))          */
     /*  =  (I + alpha inv(diag(S))*S + T*T'*beta/2)                          */
     /* --------------------------------------------------------------------- */
#if !defined(ML_NEW_ENRICH)
     if (beta != 0.0) {
       if (ml_nodes->comm->ML_mypid == 0) {
         printf("\n\nDoing old enriched prolongator");
         printf("    *******  beta = %e *******\n\n",beta);
         fflush(stdout);
       }

       TTtransPe = ML_Operator_Create(Pn_coarse->comm);
       ML_rap((*Tmat_array)[grid_level+1], (*Tmat_trans_array)[grid_level+1], 
          Pe,TTtransPe,ML_CSR_MATRIX);

       /* scale TTtransPe by -beta/2 */

       dtemp = -beta/2;
       csr_data = (struct ML_CSR_MSRdata *) TTtransPe->data;
       for (i = 0; i < csr_data->rowptr[TTtransPe->outvec_leng]; i++)
         csr_data->values[i] *= dtemp;

       newPe = ML_Operator_Create(ml_edges->comm);
       ML_Operator_Add(Pe, TTtransPe, newPe, ML_CSR_MATRIX,1.);

       ML_Operator_Destroy(&TTtransPe);

       ML_Operator_Move2HierarchyAndDestroy_fragile(newPe,
                            &(ml_edges->Pmat[grid_level]));
     }
     else
       if (ml_nodes->comm->ML_mypid == 0 && 9 < ML_Get_PrintLevel()) {
         printf("\n\nnot compiled with enriched prolongator option\n\n");fflush(stdout);}
#else
       if (ml_nodes->comm->ML_mypid == 0) {
         printf("\n\nDoing new enriched prolongator\n\n");
         fflush(stdout);
       }
     v = (double *) ML_allocate( (Pe->invec_leng+Pe->outvec_leng) * sizeof(double));
     for (i=0; i < Pe->invec_leng; i++) v[i] = 1.0;
     /*
     w = (double *) ML_allocate( (Pe->invec_leng+Pe->outvec_leng) * sizeof(double));
     for (i=0; i < Pe->outvec_leng; i++) w[i] = 1.0;
     */
     /*
     randvec = (double *) ML_allocate( (Pe->invec_leng+Pe->outvec_leng) * sizeof(double));
     ML_random_vec(randvec, Pe->outvec_leng, Pe->comm);
     dtemp = sqrt(ML_gdot(Pe->outvec_leng, randvec,randvec, Pe->comm));
     for (i=0; i<Pe->outvec_leng; i++) randvec[i] = randvec[i] / dtemp; 
     */
     /*
     fid1 = fopen("randvec.txt","w");
     for (i=0; i<Pe->outvec_leng; i++) {
        randvec[i] = randvec[i] / dtemp; 
        fprintf(fid1,"%20.16e\n",randvec[i]);
     }
     fclose(fid1);
     */

     V = &(ml_edges->Amat[grid_level+1]);
     /* W = T * (D_n' \ (T' * M)) */
     tmpmat = ML_Operator_Create(Pe->comm);
     Tfine_trans = (*Tmat_trans_array)[grid_level+1];
     /* T'*M */
     ML_2matmult((*Tmat_trans_array)[grid_level+1],
                 &(ml_edges->Amat[grid_level+1]),
                 tmpmat, ML_CSR_MATRIX );
     /*
     ML_Operator_Print(&(ml_edges->Amat[grid_level+1]),"Amat");
     ML_Operator_Print((*Tmat_trans_array)[grid_level+1],"Ttmat");
     */
     /*
     printf("||T'*M*w||^2 = %10.6e\n",checkit(tmpmat,randvec));
     */
     /* Dn' \ (T'*M) */
     ML_Operator_Get_Diag(&(ml_nodes->Amat[grid_level+1]),
                          ml_nodes->Amat[grid_level+1].outvec_leng, &diag);
     /*
     fid1 = fopen("diag","w");
     for (i=0; i<ml_nodes->Amat[grid_level+1].outvec_leng; i++) {
        fprintf(fid1,"%20.16e\n",diag[i]);
     }
     fclose(fid1);
     ML_Operator_Print(ml_nodes->Amat+grid_level+1,"Anodal");
     */
     csr_data = (struct ML_CSR_MSRdata *) tmpmat->data;
     for (i = 0; i < tmpmat->outvec_leng; i++)
       for (j= csr_data->rowptr[i]; j < csr_data->rowptr[i+1]; j++)
         csr_data->values[j] /= diag[i];
     /*
     printf("||D_n' \\ T'*M*w||^2 = %10.6e\n",checkit(tmpmat,randvec));
     */
     W = ML_Operator_Create(Pn_coarse->comm);
     /* T * (Dn' \ (T'*M)) */
     ML_2matmult(Tfine, tmpmat, W, ML_CSR_MATRIX );
     ML_Operator_Destroy(&tmpmat);

     /* Pe * v */
     Pv = (double *) ML_allocate( Pe->outvec_leng*sizeof(double));
     WPv = (double *) ML_allocate( W->outvec_leng*sizeof(double));
     ML_Operator_Apply(Pe, Pe->invec_leng, v, Pe->outvec_leng, Pv);

     /*
     printf("||Pv|| = %lf\n",checkit(Pe,v));
     */

     /* W * Pe * v */
     ML_Operator_Apply(W, W->invec_leng, Pv, W->outvec_leng, WPv);
     /* WPv' * WPv */
     denom = ML_gdot(W->outvec_leng, WPv,WPv, W->comm);

     /* V * Pv */
     ML_Operator_Apply(V, V->invec_leng, Pv, V->outvec_leng, v);
     alpha = ag->smoothP_damping_factor / ml_edges->Amat[grid_level+1].lambda_max;
     for (i=0; i<V->outvec_leng; i++)
       v[i] = alpha*v[i] - Pv[i] + 1.0;     /* using 1.0 b/c constant vector
                                               on fine grid */
     numer = ML_gdot(W->outvec_leng, WPv,v, W->comm);
     beta = numer/denom;
     printf("\n\nbeta = %lf   (numer=%lf, denom = %lf)\n\n",beta,numer,denom);

     /* (T * Dn' \ T'*M) * Pe */
     TTtransPe = ML_Operator_Create(Pn_coarse->comm);
     ML_2matmult(W, Pe, TTtransPe, ML_CSR_MATRIX );
     /* scale new term by beta */
     csr_data = (struct ML_CSR_MSRdata *) TTtransPe->data;
     for (i = 0; i < csr_data->rowptr[TTtransPe->outvec_leng]; i++)
       csr_data->values[i] *= beta;

     newPe = ML_Operator_Create(ml_edges->comm);
     ML_Operator_Add(Pe, TTtransPe, newPe, ML_CSR_MATRIX,1.);

     ML_Operator_Destroy(&TTtransPe);
     ML_Operator_Destroy(&W);
     ML_free(v); /*ML_free(w);*/ ML_free(Pv); ML_free(WPv);
     if (randvec != NULL) ML_free(randvec);

     ML_Operator_Move2HierarchyAndDestroy_fragile(newPe,
                            &(ml_edges->Pmat[grid_level]));
#endif /*ML_NEW_ENRICH*/

#ifdef ML_VAMPIR
  VT_end(ml_vt_enrich_Pe_state);
#endif
#else
     if (ml_nodes->comm->ML_mypid == 0 && 9 < ML_Get_PrintLevel()) {
       printf("\n\nnot compiled with enriched prolongator option\n\n");fflush(stdout);}
#endif /*if defined(ML_ENRICH) ... */

#ifdef ML_VAMPIR
  VT_begin(ml_vt_make_coarse_A_with_rap_state);
#endif

    if (ag->print_flag < ML_Get_PrintLevel()) {
        Pe = &(ml_edges->Pmat[grid_level]);
        nz_ptr = ML_Comm_GsumInt(ml_edges->comm, Pe->N_nonzeros);
        j = Pe->outvec_leng;
        i = ML_Comm_GsumInt(ml_edges->comm, j);
        if (Tfine->comm->ML_mypid==0)
           printf("(level %d) Pe: Global nonzeros = %d, global rows = %d\n", grid_level,nz_ptr, i);
     }

     ML_Operator_Set_1Levels(&(ml_edges->Pmat[grid_level]),
                 &(ml_edges->SingleLevel[grid_level]), 
                 &(ml_edges->SingleLevel[grid_level+1]));
     ML_Gen_Restrictor_TransP(ml_edges, grid_level+1, grid_level);
     ML_Gen_AmatrixRAP(ml_edges, grid_level+1, grid_level);
     /*
     if (grid_level == fine_level-1)
       ML_Operator_Print(ml_edges->Amat+grid_level,"Amat_before");
     */
#if defined(ML_NEW_ENRICH)
     /* Weed out small values in A (MSR format). */
     droptol = 1e-9;
     tmpmat = ml_edges->Amat+grid_level;
     if (tmpmat->comm->ML_mypid==0 && ag->print_flag < ML_Get_PrintLevel())
        printf("Dropping Ke(i,j) if |Ke(i,j)| <  %e\n", droptol);
     csr_data = (struct ML_CSR_MSRdata *) tmpmat->data;
     lower = csr_data->columns[0];
     nz_ptr = tmpmat->outvec_leng+1;
     for (i = 0; i < tmpmat->outvec_leng; i++) {
       for (j = lower; j < csr_data->columns[i+1]; j++) {
         if (ML_abs(csr_data->values[j]) > droptol) {
           csr_data->values[nz_ptr] = csr_data->values[j];
           csr_data->columns[nz_ptr] = csr_data->columns[j];
           nz_ptr++;
         }
       }
       lower = csr_data->columns[i+1];
       csr_data->columns[i+1] = nz_ptr;
     }
     tmpmat->N_nonzeros = nz_ptr;
#endif /*if defined(ML_ENRICH) && ... */
     Nnz_allgrids += ml_edges->Amat[grid_level].N_nonzeros;
     /*
     if (grid_level == fine_level-1)
       ML_Operator_Print(ml_edges->Amat+grid_level,"Amat_after");
     */

     if (ag->print_flag < ML_Get_PrintLevel()) {
        Pe = &(ml_edges->Amat[grid_level]);
        nz_ptr = ML_Comm_GsumInt(ml_edges->comm, Pe->N_nonzeros);
        i = Pe->outvec_leng;
        ML_gsum_scalar_int(&i,&j,ml_nodes->comm);
        if (Tfine->comm->ML_mypid==0)
           printf("(level %d) Ke: Global nonzeros = %d, global rows = %d\n", grid_level, nz_ptr,i);
     }

     Tfine = Tcoarse;
#ifdef ML_VAMPIR
  VT_end(ml_vt_make_coarse_A_with_rap_state);
#endif
  } /* Main FOR loop: for grid_level = fine_level-1 ... */

  ml_edges->ML_coarsest_level = coarsest_level;

#ifdef ML_VAMPIR
  VT_begin(ml_vt_reitzinger_cleanup_state);
#endif

  if (ML_Get_PrintLevel() > 0)
  {
    ML_gsum_scalar_int(&Nnz_allgrids,&j,ml_nodes->comm);
    ML_gsum_scalar_int(&Nnz_finegrid,&j,ml_nodes->comm);

    if (Tfine->comm->ML_mypid==0 )
    {
      if (Nnz_finegrid <= 0) 
         printf("Number of nonzeros on finest grid not given! Complexity not computed!\n");
      else
         printf("Multilevel complexity is %e\n",
               ((double) Nnz_allgrids)/((double) Nnz_finegrid));
    }
  }

  if (created_ag_obj == 1) ML_Aggregate_Destroy(&ag);
#ifdef ML_VAMPIR
  VT_end(ml_vt_reitzinger_cleanup_state);
#endif

  t0 = GetClock() - t0;
  t0 = ML_gsum_double(t0, ml_nodes->comm);
  t0 = t0/((double) ml_nodes->comm->ML_nprocs);
  if (ML_Get_PrintLevel() > 0)
    if (Tfine->comm->ML_mypid==0)
      printf("AMG setup time \t= %e seconds\n",t0);

  return(Nlevels_nodal);
}

/******************************************************************************/

int ML_MGHierarchy_ReitzingerDestroy(int finest_level,
                     ML_Operator ***Tmat_array, ML_Operator ***Tmat_trans_array)
{
    int level;

    if (*Tmat_array != NULL)
    {
       for (level = finest_level; level >= 0; level--)
       {
          ML_Operator_Destroy((*Tmat_array)+level);
          (*Tmat_array)[level] = NULL;
       }
       ML_free(*Tmat_array);
       *Tmat_array = NULL;
    }

    if (*Tmat_trans_array != NULL)
    {
       for (level = finest_level; level >= 0; level--)
       {
          ML_Operator_Destroy((*Tmat_trans_array)+level);
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
   int         Nfine;
   double      max_eigen = -1.;
   ML_Operator *Amat, *Pmatrix = NULL, *AGGsmoother = NULL;
   struct      ML_AGG_Matrix_Context widget;
   ML_Krylov   *kdata;

   if ( smoothP_damping_factor == 0.0 ) return 0;

   Amat     = (ML_Operator *) data;
   Nfine    = Amat->outvec_leng;
   Pmatrix =  &(ml->Pmat[clevel]);


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

   ML_2matmult(AGGsmoother, Pmatrix, SPn_mat, ML_CSR_MATRIX );
   ML_Operator_Destroy(&AGGsmoother);

   return 0;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int ml_clean_Trecorder(struct ml_linked_list ***Trecorder ,int N)
{
  struct ml_linked_list *current, *tmp;
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

int ml_record_entry(struct ml_linked_list **Trecorder, int lower, int therow)
{
  /* Record the edge (lower, upper) in the linked list table pointed to */
  /* by Trecorder. The matrix row index associated with this edge is    */
  /* recorded as well. It is assumed that this edge is not already in   */
  /* the table and that lower < upper.                                  */
  /* This routine is used in conjunction with 'ml_dup_entry' to see if     */
  /* edges have already been added to the table.                        */

  struct ml_linked_list *prev_head;

  prev_head = Trecorder[lower];
  Trecorder[lower] = (struct ml_linked_list *) ML_allocate(sizeof(struct ml_linked_list));
  Trecorder[lower]->next = prev_head;
  Trecorder[lower]->duplicate_row  = therow;
  return 1;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int ml_dup_entry(int node1, int node2, struct ml_linked_list **Trecorder,
	int Tcols[], int Trowptr[], int *lower, int *upper, int *duplicate_row)
{
  /* This routine is used to check if a pair (node1,node2) is already  */
  /* in the table pointed to by Trecorder (this table is normally      */
  /* constructed via ml_record_entry). On return,                         */
  /*     duplicate_row    Set to -1 if the pair is not in the table.   */
  /*                      Set to indicate matrix row associated with   */
  /*                      pair if (node1,node2) is in the table.       */

  struct ml_linked_list *current;
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

  int i, j, max_nz_per_row = 1;
  struct ml_linked_list **Trecorder;
  struct ML_CSR_MSRdata *Pn, *SPn, *Tfine, *Tcoarse, *Pe;
  double thesign;
  int left,leftagg, right, rightagg;
  int coef_count, coef_cols[100], jcoef_count, k;
  double alpha0, beta_m, coef_vals[100];
  int Trowcount = 0, Tnzcount = 0, Pnzcount = 0;
  int *SPn_columns, *SPn_rowptr, *Tfine_columns, *Pn_columns, *Pn_rowptr;
  int *Tfine_rowptr;
  double *SPn_values, *Tfine_values;
  int jjj;
  int Trownnz;
  int tmpcol;

  /* pull out a bunch of pointers */

#ifdef ML_CPP
  Pn            = (ML_CSR_MSRdata *)Pn_mat->data;
#else
  Pn            = Pn_mat->data;
#endif
  Pn_columns    = Pn->columns;
  Pn_rowptr     = Pn->rowptr;
#ifdef ML_CPP
  Tfine         = (ML_CSR_MSRdata *)Tfine_mat->data;
#else
  Tfine         = Tfine_mat->data;
#endif
  Tfine_columns = Tfine->columns;
  Tfine_values  = Tfine->values;
  Tfine_rowptr  = Tfine->rowptr;
#ifdef ML_CPP
  SPn           = (ML_CSR_MSRdata *)SPn_mat->data;
#else
  SPn           = SPn_mat->data;
#endif
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

  Trecorder = (struct ml_linked_list **) ML_allocate((Pn_mat->outvec_leng+1)*
					     sizeof(struct ml_linked_list *));
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
    Trownnz = 0;
    for (jjj = Tfine_rowptr[i]; jjj < Tfine_rowptr[i+1]; jjj++)
       if (Tfine_values[jjj] != 0) Trownnz++;

    /* special case when Tfine(i,:) has only one entry */
    if (Trownnz == 1) {
      if (Tfine_values[Tfine_rowptr[i]] == 0.) 
      {
          Tfine_values[Tfine_rowptr[i]] = Tfine_values[Tfine_rowptr[i+1]];
          Tfine_values[Tfine_rowptr[i+1]] = 0.;
          tmpcol = Tfine_columns[Tfine_rowptr[i]];
          Tfine_columns[Tfine_rowptr[i]] = Tfine_columns[Tfine_rowptr[i+1]];
          Tfine_columns[Tfine_rowptr[i+1]] = tmpcol;
      }
      if (Tfine_values[Tfine_rowptr[i]] == 1.) thesign = 1.;
      else thesign = -1.;
      left  = Tfine_columns[Tfine_rowptr[i]];
      coef_count = 0;
      for (j = SPn_rowptr[left]; j < SPn_rowptr[left+1]; j++) {
         coef_cols[coef_count  ] = SPn_columns[j];
         coef_vals[coef_count++] = thesign*SPn_values[j];
      }
      ml_comp_Pe_entries(coef_cols, coef_vals, coef_count, -1, Trecorder,
             &Trowcount, &Tnzcount, Tcoarse, 
             &Pnzcount,Pe->columns,Pe->values);
    }

    /* normal case when Tfine(i,:) has only two entries */
    else if (Trownnz == 2) {

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
        /* else beta0 = SPn_values[j]; */
      }

      ml_comp_Pe_entries(coef_cols, coef_vals, coef_count, leftagg, 
			   Trecorder, &Trowcount, &Tnzcount, Tcoarse,
			   &Pnzcount,Pe->columns,Pe->values);
	coef_count = 0;
      }
      else {
	/* alpha_m = 0.; */ alpha0 = 0.;
	beta_m = 0.; /* beta0 = 0.; */
	/* case 2 */

	/* copy alpha into coef */

	for (j = SPn_rowptr[left]; j < SPn_rowptr[left+1]; j++) {
	  /*
	  if (SPn_columns[j] == leftagg)
	    alpha_m = SPn_values[j];
	  else if (SPn_columns[j] == rightagg)
	  */
	  if (SPn_columns[j] != leftagg) {
	    if (SPn_columns[j] == rightagg)
	      alpha0 = SPn_values[j];
	    else {
	      coef_cols[coef_count  ] = SPn_columns[j];
	      coef_vals[coef_count++] = -SPn_values[j];
	    }
	  }
	}
	ml_comp_Pe_entries(coef_cols, coef_vals, coef_count, leftagg,
			   Trecorder, &Trowcount, &Tnzcount, Tcoarse, 
			   &Pnzcount,Pe->columns,Pe->values);
	coef_count = 0;

      /* copy beta into coef */
	for (j = SPn_rowptr[right]; j < SPn_rowptr[right+1]; j++) {
	  if (SPn_columns[j] == leftagg)
	    beta_m = SPn_values[j];
	  /*
	  else if (SPn_columns[j] == rightagg)
	    beta0 = SPn_values[j];
	  else {
	  */
	  else if (SPn_columns[j] != rightagg) {
	    coef_cols[coef_count  ] = SPn_columns[j];
	    coef_vals[coef_count++] = SPn_values[j];
	  }
	}
	ml_comp_Pe_entries(coef_cols, coef_vals, coef_count, rightagg,
			   Trecorder, &Trowcount, &Tnzcount, Tcoarse, 
			   &Pnzcount,Pe->columns,Pe->values);
	coef_count = 1;
	coef_cols[0] = leftagg;
	coef_vals[0] = beta_m+alpha0-1.;
	ml_comp_Pe_entries(coef_cols, coef_vals, coef_count, rightagg,
			   Trecorder, &Trowcount, &Tnzcount, Tcoarse,
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

  ml_clean_Trecorder(&Trecorder,Pn_mat->outvec_leng);

  return 1;

}

int ml_comp_Pe_entries(int coef_cols[], double coef_values[], int coef_count,
		       int leftagg, struct ml_linked_list **Trecorder,
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

      ml_record_entry(Trecorder, lower, *Trowcount);
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

#include "ml_operator_blockmat.h"
/****************************************************************************/
/* Generate the new grid hierarchy for 2x2 block matrices from the original */
/* hierarchy stored in 'ml_edges'.                                          */
/****************************************************************************/

int ML_Gen_Hierarchy_ComplexMaxwell(ML *ml_edges, ML **new_ml , 
				    ML_Operator *originalM)
{

   int mesh_level, old_mesh_level, levels;
   ML_Operator *original, *blockmat, *newM, *lastM;
   ML  *block_ml;

   mesh_level = ml_edges->ML_finest_level;

   /* create a new empty hierarchy with the same number of levels */
   /* as in ml_edges.                                             */
   ML_Create(&block_ml,ml_edges->ML_num_levels);
   *new_ml = block_ml;


   /* Make A on the fine level into a 2x2 block matrix */

   levels = 1;
   original = &(ml_edges->Amat[mesh_level]);
   blockmat = &(block_ml->Amat[mesh_level]);
   ML_Operator_Gen_blockmat(blockmat, original , originalM );
   ML_Operator_blockmat_set_M_mat_destroy(blockmat,ML_NO);

   lastM = originalM;
   while( ml_edges->SingleLevel[mesh_level].Rmat->to != NULL) {
     levels++;
     old_mesh_level = mesh_level;
     mesh_level = ml_edges->SingleLevel[mesh_level].Rmat->to->levelnum;

     /* Make 2x2 block diagonal P */

     original = &(ml_edges->Pmat[mesh_level]);
     blockmat = &(block_ml->Pmat[mesh_level]);
     ML_Operator_Gen_blockmat(blockmat, original , NULL );
     /* This stuff sets the 'to' and 'from' field in P */
     /* which indicates from what level we interpolate */
     /* and to what level the interpolation goes.      */
     ML_Operator_Set_1Levels(blockmat, &(block_ml->SingleLevel[mesh_level]), 
			     &(block_ml->SingleLevel[old_mesh_level]));

     /* Make 2x2 block diagonal R */

     original = &(ml_edges->Rmat[old_mesh_level]);
     blockmat = &(block_ml->Rmat[old_mesh_level]);
     ML_Operator_Gen_blockmat(blockmat, original, NULL );
     /* This stuff sets the 'to' and 'from' field in P */
     /* which indicates from what level we interpolate */
     /* and to what level the interpolation goes.      */
     ML_Operator_Set_1Levels(blockmat,
                             &(block_ml->SingleLevel[old_mesh_level]), 
			     &(block_ml->SingleLevel[mesh_level]));
				  
     /* Make 2x2 block diagonal A */

     original = &(ml_edges->Amat[mesh_level]);
     blockmat = &(block_ml->Amat[mesh_level]);
     /*  newM = ML_Operator_Create(ml_edges->comm);
	 ML_rap(&(ml_edges->Rmat[old_mesh_level]), original, 
	 &(ml_edges->Pmat[mesh_level]), newM, ML_CSR_MATRIX);
     */
     newM = ML_Operator_Create(ml_edges->comm);
     ML_rap(&(ml_edges->Rmat[old_mesh_level]), lastM, 
            &(ml_edges->Pmat[mesh_level]), newM, ML_CSR_MATRIX);
     lastM = newM;
 
     /* comment these two out if you want to do rap */
     ML_Operator_Gen_blockmat(blockmat, original, newM);


#ifdef GREG
     blockmat->lambda_max = 3.286003165408382;
     blockmat->lambda_max_img = 3.148601063739487e-02;
     blockmat->lambda_max = 5.5*2.49357;
     blockmat->lambda_max_img = 5.5*.0142284;
#endif

     /* RAP works too */
     /* ML_Gen_AmatrixRAP(block_ml, old_mesh_level, mesh_level); */

     /* rst: I'm not sure ... but something like the following */
     /* should work for T if needed?                           */
     /*
       original = &(Tmat_array[mesh_level]);
       blockmat = &(blk_Tmat_array[mesh_level]);
       ML_Operator_Gen_blockmat(blockmat, original);
       original = &(Tmat_trans_array[mesh_level]);
       blockmat = &(blk_Tmat_trans_array[mesh_level]);
       ML_Operator_Gen_blockmat(blockmat, original);
     */

   }

   return levels;

}

#if defined(ML_NEW_ENRICH)
double checkit(ML_Operator *A, double *v)
{
  int i;
  double *y, result;

  y = (double *) ML_allocate(A->outvec_leng * sizeof(double));
  ML_Operator_Apply(A,A->invec_leng,v,A->outvec_leng,y);
  /* result = sqrt(ML_gdot(A->outvec_leng, y,y, A->comm)); */
  result = ML_gdot(A->outvec_leng, y,y, A->comm);
  ML_free(y);
  return result;
}
#endif
