#include "ml_operator_blockmat.h"
/*****************************************************************/
/* ML block matvec function that corresponds to an equivalent */
/* real form                                                     */
/*                    |  Ke    -M  |  | x |       | b |          */
/*                    |   M    Ke  |  | y |   =   | c |          */
/* to solve the problem                                          */
/*                    (Ke + i M ) (x + i y) = b + i c            */
/* where Ke and M are ML Operators.                               */ 
/*****************************************************************/

int ML_Operator_blockmat_matvec(ML_Operator *data, int inlen, double invec[],
		      int outlen, double outvec[])
{
  ML_Operator *mat;
  struct ML_Operator_blockmat_data *ML_Operator_blockmat_data;
  double *z;
  int i;

  mat = (ML_Operator *) data;
  ML_Operator_blockmat_data = (struct ML_Operator_blockmat_data *) mat->data;
  /* workspace vector */
  z = (double *) ML_allocate (outlen * sizeof(double) );

  /* multiply by (1,1) block (stiffness)*/
  ML_Operator_blockmat_data->Ke_matvec((ML_Operator *) ML_Operator_blockmat_data->Ke_matvec_data, inlen/2,
				     invec, outlen/2, outvec);

  /* multiply by (1,2) block (-mass)*/
  if (ML_Operator_blockmat_data->M_matvec != NULL) {
    ML_Operator_blockmat_data->M_matvec((ML_Operator *) ML_Operator_blockmat_data->M_matvec_data, inlen/2,
				     &(invec[inlen/2]), outlen/2, z);
    for (i=0; i< outlen/2; i++) outvec[i] -= z[i];
  }
  /***********/
  /* multiply by (2,2) block (stiffness)*/
  ML_Operator_blockmat_data->Ke_matvec((ML_Operator *) ML_Operator_blockmat_data->Ke_matvec_data, inlen/2,
				     &(invec[inlen/2]), outlen/2, 
				     &(outvec[outlen/2]));
  /* multiply by (2,1) block (mass)*/
  if (ML_Operator_blockmat_data->M_matvec != NULL) {
    ML_Operator_blockmat_data->M_matvec((ML_Operator *) ML_Operator_blockmat_data->M_matvec_data, inlen/2,
				     invec, outlen/2, z);
    for (i=0; i < outlen/2; i++) outvec[i+outlen/2] += z[i];
  }

  ML_free(z);
  return 1;
}

/*****************************************************************/
/* ML block communcation that corresponds to the communication   */
/* necessary to perform a matrix-vector multiply of an equivalent*/
/* real form matrix                                              */
/*                    |  Ke    -M  |  | x |       | b |          */
/*                    |   M    Ke  |  | y |   =   | c |          */
/*                                                               */
/*****************************************************************/

int ML_Operator_blockmat_comm( double *x, void *data)
{
  struct ML_Operator_blockmat_data *ML_Operator_blockmat_data;
  int N_Ke, Nghost, i;
  double *temp;
  ML_Operator *mat;

  mat = (ML_Operator *) data;
  ML_Operator_blockmat_data = (struct ML_Operator_blockmat_data *) mat->data;
  if (ML_Operator_blockmat_data->Ke_comm_data == NULL) return(0);

  N_Ke  = ML_Operator_blockmat_data->N_Ke;
  Nghost = ML_Operator_blockmat_data->Nghost;

  temp = (double *) ML_allocate(sizeof(double)*(N_Ke+Nghost));

  /* fill a temp vector call subblock communicate routine and restore */
  /* data back into larger vector.                                    */

  for (i = 0; i < N_Ke ; i++) temp[i] = x[i];
  for (i = 0; i < Nghost; i++) temp[i+N_Ke] = x[2*N_Ke+2*i];

#ifdef ML_CPP
  ML_exchange_bdry(temp, (ML_CommInfoOP *)ML_Operator_blockmat_data->Ke_comm_data, N_Ke,mat->comm,
#else
  ML_exchange_bdry(temp, ML_Operator_blockmat_data->Ke_comm_data, N_Ke,mat->comm,
#endif
                   ML_OVERWRITE,NULL);

  for (i = 0; i < N_Ke ; i++) x[i] = temp[i];
  for (i = 0; i < Nghost; i++) x[2*N_Ke+2*i] = temp[i+N_Ke];

  /* fill a temp vector call subblock communicate routine and restore */
  /* data back into larger vector.                                    */

  for (i = 0; i < N_Ke ; i++) temp[i] = x[i+N_Ke];
  for (i = 0; i < Nghost; i++) temp[i+N_Ke] = x[2*N_Ke+2*i + 1];

#ifdef ML_CPP
  ML_exchange_bdry(temp, (ML_CommInfoOP *)ML_Operator_blockmat_data->Ke_comm_data, N_Ke,mat->comm,
#else
  ML_exchange_bdry(temp, ML_Operator_blockmat_data->Ke_comm_data, N_Ke,mat->comm,
#endif
                   ML_OVERWRITE,NULL);

  for (i = 0; i < N_Ke ; i++) x[i+N_Ke] = temp[i];
  for (i = 0; i < Nghost; i++) x[2*N_Ke+2*i + 1] = temp[i+N_Ke];

  ML_free(temp);
  return(0);
}

/*****************************************************************/
/* ML getrow for the equivalent real form matrix                 */
/*                    |  Ke    -M  |                             */
/*                    |   M    Ke  |                             */
/*                                                               */
/* Note: WE ASSUME THAT KE AND M HAVE 'IDENTICAL' GHOST NODES.   */
/*****************************************************************/

int ML_Operator_blockmat_getrow(ML_Operator *data, int N_requested,
			  int requested_rows[], int allocated,
			  int columns[], double values[], int row_lengths[])
{
  ML_Operator *mat;
  struct ML_Operator_blockmat_data *ML_Operator_blockmat_data;
  int newrow, status = 1, i;
  int *workcol;
  double *workval;
  int work_lengths[1], invec_sub;

  mat = (ML_Operator *) data;
  ML_Operator_blockmat_data = (struct ML_Operator_blockmat_data *) mat->data;

  workcol = ML_Operator_blockmat_data->cols;
  workval = ML_Operator_blockmat_data->vals;

  work_lengths[0] = 0;
  row_lengths[0] = 0;

  if (N_requested != 1) return(1);
  invec_sub = mat->invec_leng/2;

  /* (1,1) and (1,2) blocks */
  if (requested_rows[0] < mat->outvec_leng/2)
  {
    /* (1,1) block */
    status = ML_Operator_blockmat_data->Ke_getrow((ML_Operator *) ML_Operator_blockmat_data->Ke_getrow_data,
				      N_requested, requested_rows, allocated, 
				      columns, values, row_lengths);

    if (status == 0) return(status);

    for (i=0; i< row_lengths[0]; i++) {
      if (columns[i] >= invec_sub) columns[i] *= 2;
    }


    if (ML_Operator_blockmat_data->M_getrow != NULL) {

      workval = &(values [ row_lengths[0] ]);
      workcol = &(columns[ row_lengths[0] ]);

      /* (1,2) block */
      status = ML_Operator_blockmat_data->M_getrow((ML_Operator *) ML_Operator_blockmat_data->M_getrow_data,
				      N_requested, requested_rows, allocated -
				      row_lengths[0], 
				      workcol, workval, work_lengths);

      if (status == 0) return(status);

      for (i=0; i< work_lengths[0]; i++) {
	workval[i] = -workval[i];     /* josh says -M */
      }
      for (i=0; i< work_lengths[0]; i++) {
	if (workcol[i] >= invec_sub) workcol[i] += workcol[i] + 1;
	else workcol[i] += invec_sub;
      }
    }
  }
  /* (2,1) and (2,2) blocks */
  else {
    /* shift requested row to (1,1) block */
    newrow = requested_rows[0] - invec_sub;

    /* (2,1) block */
    if (ML_Operator_blockmat_data->M_getrow != NULL) {
      status = ML_Operator_blockmat_data->M_getrow((ML_Operator *) 
		                      ML_Operator_blockmat_data->M_getrow_data,
				      N_requested, &newrow, allocated, 
				      columns, values, row_lengths);

      if (status == 0) return(status);

      for (i=0; i< row_lengths[0]; i++) {
	if (columns[i] >= invec_sub) columns[i] *= 2;
      }

    }
    workval = &(values [ row_lengths[0] ]);
    workcol = &(columns[ row_lengths[0] ]);

    /* (2,2) block */
    status = ML_Operator_blockmat_data->Ke_getrow((ML_Operator *) ML_Operator_blockmat_data->Ke_getrow_data,
					N_requested, &newrow, allocated-
					row_lengths[0], 
					workcol, workval, work_lengths);

    if (status == 0) return(status);

    /* post process data so that columns correspond to (2,2) block */

    for (i=0; i< work_lengths[0]; i++) {
      if (workcol[i] >= invec_sub) workcol[i] += workcol[i] + 1;
      else workcol[i] += invec_sub;
    }

  }
  /*make sure columns is long enough for the concatenation*/

  if (row_lengths[0] + work_lengths[0] > allocated) {
    return(0);
  }

  for (i=0; i<work_lengths[0]; i++) {
    columns[i + row_lengths[0]] = workcol[i];
    values[i + row_lengths[0]] = workval[i];
  }
  row_lengths[0] += work_lengths[0]; 

  return(1);
}

/*****************************************************************/
/* Convert 'Ke/M' into a 2x2 block matrix of the following form: */
/*                         ( Ke   -M  )                          */
/*      blockmat =         ( M    Ke  )                          */
/*****************************************************************/

int  ML_Operator_Gen_blockmat(ML_Operator *blockmat, ML_Operator *Ke,
                          ML_Operator *M)
{

  struct ML_Operator_blockmat_data *ML_Operator_blockmat_data;
  int scale_fact = 2; /* b/c we're building a 2x2 system */
  int  Nghost, Nneighbors, *neighbors, i;

  /* fill data structure which ML_Operator_blockmat_matvec, ML_Operator_blockmat_comm,   */
  /* and ML_Operator_blockmat_getrow will use. Shove new data structure into */
  /* into blockmat along with matvec/getrow/comm function ptrs.  */

  ML_Operator_Clean(blockmat);
  ML_Operator_Init(blockmat,Ke->comm);

  blockmat->max_nz_per_row = 30; /* we shouldn't need to set this */

  ML_Operator_blockmat_data = (struct ML_Operator_blockmat_data *) ML_allocate(
				       sizeof(struct ML_Operator_blockmat_data));

  ML_Operator_blockmat_data->Ke_diag = NULL;
  ML_Operator_blockmat_data->M_diag  = NULL;
  ML_Operator_blockmat_data->M_mat   = NULL;
  ML_Operator_blockmat_data->destroy_M_mat = ML_YES;

  ML_Operator_blockmat_data->N_Ke = Ke->invec_leng;
  Nghost = 0;
  if (Ke->getrow->pre_comm != NULL) {
    Nneighbors = ML_CommInfoOP_Get_Nneighbors(Ke->getrow->pre_comm);
    neighbors  = ML_CommInfoOP_Get_neighbors(Ke->getrow->pre_comm);
    for (i = 0; i < Nneighbors; i++) {
      Nghost += ML_CommInfoOP_Get_Nrcvlist(Ke->getrow->pre_comm,neighbors[i]);
    }
    ML_free(neighbors);
  }
  ML_Operator_blockmat_data->Nghost = Nghost;


  ML_Operator_blockmat_data->cols = (int    *) ML_allocate(100 * sizeof(int));
  ML_Operator_blockmat_data->vals = (double *) ML_allocate(100 * sizeof(double));

  /* setup matvec for diagonal part */

  ML_Operator_blockmat_data->Ke_matvec = Ke->matvec->func_ptr;
  ML_Operator_blockmat_data->Ke_matvec_data = Ke;

  /* setup matvec for offdiagonal part */

  ML_Operator_blockmat_data->M_matvec = NULL;
  ML_Operator_blockmat_data->M_matvec_data = NULL;
  if (M != NULL) {
    ML_Operator_blockmat_data->M_matvec = M->matvec->func_ptr;
    ML_Operator_blockmat_data->M_matvec_data = M;
  } 
  ML_Operator_Set_ApplyFuncData(blockmat, scale_fact*Ke->invec_leng, 
				scale_fact*Ke->outvec_leng, 
				ML_Operator_blockmat_data,
				scale_fact*Ke->outvec_leng,
				ML_Operator_blockmat_matvec,0);

  /* set getrow for diagonal block */

  ML_Operator_blockmat_data->Ke_getrow = Ke->getrow->func_ptr;
  ML_Operator_blockmat_data->Ke_getrow_data = Ke;

  ML_Operator_blockmat_data->M_getrow = NULL;
  ML_Operator_blockmat_data->M_getrow_data = NULL;

  /* set getrow for offdiagonal block */
  if (M != NULL) {
    ML_Operator_blockmat_data->M_getrow = M->getrow->func_ptr;
    ML_Operator_blockmat_data->M_getrow_data = M;
  }

  ML_Operator_Set_Getrow(blockmat, scale_fact*
			 Ke->outvec_leng, ML_Operator_blockmat_getrow);

  if (Ke->getrow->pre_comm != NULL) {
    ML_Operator_blockmat_data->Ke_comm_data = Ke->getrow->pre_comm;
    ML_CommInfoOP_Generate( &(blockmat->getrow->pre_comm), 
			    ML_Operator_blockmat_comm, blockmat, 
			    Ke->comm, 2*Ke->invec_leng, 2*Nghost);

  }
  else blockmat->getrow->pre_comm = ML_CommInfoOP_Create(); 

  blockmat->data_destroy = ML_Operator_blockmatdata_Destroy;

  if (M != NULL) {
    ML_Operator_blockmat_data->M_mat =  M;
  }
  if (Ke != NULL) {
    ML_Operator_blockmat_data->Ke_mat =  Ke;
  }

  return 1;
}
void  ML_Operator_blockmatdata_Destroy(void *data)
{
  struct ML_Operator_blockmat_data *temp;

  temp  = (struct ML_Operator_blockmat_data *) data;
  if (temp != NULL) {
    if (temp->cols != NULL) ML_free(temp->cols);
    if (temp->vals != NULL) ML_free(temp->vals);
    if (temp->Ke_diag != NULL) ML_free(temp->Ke_diag);
    if (temp->M_diag != NULL ) ML_free(temp->M_diag);

    if ( (temp->M_mat != NULL) && (temp->destroy_M_mat == ML_YES) ) {
      ML_Operator_Destroy(&(temp->M_mat));
    }
    ML_free(temp);
  }
}
int ML_Operator_blockmat_set_M_mat_destroy(ML_Operator *blockmat,
					   int yes_or_no)
{
  struct ML_Operator_blockmat_data *temp;

  temp  = (struct ML_Operator_blockmat_data *) blockmat->data;
  temp->destroy_M_mat = yes_or_no;
  return 1;
}
