#include "ml_config.h"
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS)
#include "ml_RefMaxwell_Utils.h"
#include "ml_epetra.h"
#include "ml_epetra_utils.h"
#include "ml_mat_formats.h"
#include <fstream>

// ================================================ ====== ==== ==== == =
int ML_Epetra::CSR_getrow_ones(ML_Operator *data, int N_requested_rows, int requested_rows[],
   int allocated_space, int columns[], double values[], int row_lengths[])
{
   /* register (C++11 deprecates this storage class specifier) */ int    *bindx, j;
   int     *rowptr,  row, itemp;
   struct ML_CSR_MSRdata *input_matrix;
   ML_Operator *mat_in;

   row            = *requested_rows;
   mat_in = (ML_Operator *) data;
   input_matrix = (struct ML_CSR_MSRdata *) ML_Get_MyGetrowData(mat_in);
   rowptr = input_matrix->rowptr;
   itemp = rowptr[row];
   *row_lengths = rowptr[row+1] - itemp;


   if (*row_lengths > allocated_space) {
    ML_avoid_unused_param( (void *) &N_requested_rows);
    return(0);
  }

   bindx  = &(input_matrix->columns[itemp]);
   for (j = 0 ; j < *row_lengths; j++) {
      *columns++ = *bindx++;
   }
   for (j = 0 ; j < *row_lengths; j++) {
      *values++  = 1;
   }
   return(1);
}

// ================================================ ====== ==== ==== == =
void ML_Epetra::IVOUT(const Epetra_IntVector & A, const char *of){
  std::ofstream os(of);
  int i;
  int NumProc=A.Map().Comm().NumProc();
  int MyPID  =A.Map().Comm().MyPID();

  for (int iproc=0; iproc < NumProc; iproc++) {
    if (MyPID==iproc) {
      int MyLength = A.MyLength();
      for (i=0; i<MyLength; i++) {
          os.width(20);
          os << A[i]<<std::endl;
      }
      os << std::flush;
    }/*end if*/
    // Do a few global ops to give I/O a chance to complete
    A.Map().Comm().Barrier();
    A.Map().Comm().Barrier();
    A.Map().Comm().Barrier();
  }/*end for*/
}/*end IVOUT*/

// ================================================ ====== ==== ==== == =
// Copied from ml_agg_genP.c
static int ML_Aux_Getrow(ML_Operator *data, int N_requested_rows, int requested_rows[],
                         int allocated_space, int columns[], double values[],
                         int row_lengths[])
{
  int ierr;
  int i, j, count, mod;
  int BlockCol, BlockRow, RowMod;
  double DiagValue = 0.0;
  int DiagID;
  int* Filter;

  ierr = (*(data->aux_data->aux_func_ptr))(data, N_requested_rows, requested_rows,
                                      allocated_space, columns, values, row_lengths);
  if (ierr == 0)
    return(0);

  if (N_requested_rows != 1) {
    fprintf(stderr, "ML_Aux_Getrow() works only if N_requested_rows == 1\n"
            "(file %s, line %d)\n",
            __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }

  /* new part */
  mod       = data->invec_leng / data->aux_data->filter_size;
  BlockRow  = requested_rows[0] / mod;
  RowMod    = requested_rows[0] % mod;
  Filter    = data->aux_data->filter[BlockRow];
  count     = 0;
  DiagID    = -1;
  DiagValue = 0.0;

  for (i = 0 ; i < row_lengths[0] ; ++i)
  {
    BlockCol = columns[i] / mod;

    if (BlockCol ==  BlockRow)
    {
      columns[count] = columns[i];
      values[count]  = values[i];
      if (requested_rows[0] == columns[i]) DiagID = count;
      ++count;
      goto after;
    }

    /* different block col, same equation */
    for (j = 0 ; j < Filter[0] ; ++j)
    {
      /* look for elements to discard */
      if (Filter[j + 1] == BlockCol)
      {
        /* We modify the matrix diagonal so that the constant vector */
        /* is still in the null space (assuming that the constant    */
        /* was in the original matrix null space). For block pde     */
        /* systems we need to check that both the column and the row */
        /* correspond to the same DOF within the node.               */
        /* Note: This code will not preserve the null space if it is */
        /* not given by piecewise constants.                         */

        if (columns[i] % mod == RowMod) {
          DiagValue += values[i];
        }
        goto after;
      }
    }
    columns[count] = columns[i];
    values[count]  = values[i];
    ++count;
after:
    continue;
  }

  if (DiagID == -1)
  {
    fprintf(stderr, "Diagonal not defined for row %d\n", requested_rows[0]);
  }

  values[DiagID] += DiagValue;
  row_lengths[0] = count;

  return(ierr);
}


// ================================================ ====== ==== ==== == =
int RefMaxwell_SetupCoordinates(ML_Operator* A, Teuchos::ParameterList &List_, double *&coordx, double *&coordy, double *&coordz)
// Coppied From int ML_Epetra::MultiLevelPreconditioner::SetupCoordinates()
{
  double* in_x_coord = 0;
  double* in_y_coord = 0;
  double* in_z_coord = 0;
  int NumPDEEqns_ =1;  // Always use 1 because A is a nodal matrix.

  // For node coordinates
  in_x_coord = List_.get("x-coordinates", (double *)0);
  in_y_coord = List_.get("y-coordinates", (double *)0);
  in_z_coord = List_.get("z-coordinates", (double *)0);

  if (!(in_x_coord == 0 && in_y_coord == 0 && in_z_coord == 0))
    {
      ML_Operator* AAA = A;

      int n = AAA->invec_leng, Nghost = 0;

      if (AAA->getrow->pre_comm)
        {
          if (AAA->getrow->pre_comm->total_rcv_length <= 0)
            ML_CommInfoOP_Compute_TotalRcvLength(AAA->getrow->pre_comm);
          Nghost = AAA->getrow->pre_comm->total_rcv_length;
        }

      std::vector<double> tmp(Nghost + n);
      for (int i = 0 ; i < Nghost + n ; ++i)
        tmp[i] = 0.0;

      n /= NumPDEEqns_;
      Nghost /= NumPDEEqns_;

      if (in_x_coord)
        {
          double* x_coord = (double *) ML_allocate(sizeof(double) * (Nghost+n));

          for (int i = 0 ; i < n ; ++i)
            tmp[i * NumPDEEqns_] = in_x_coord[i];

          ML_exchange_bdry(&tmp[0],AAA->getrow->pre_comm, NumPDEEqns_ * n,
                           AAA->comm, ML_OVERWRITE,NULL);

          for (int i = 0 ; i < n + Nghost ; ++i)
            x_coord[i] = tmp[i * NumPDEEqns_];

          coordx = x_coord;
        }

      if (in_y_coord)
        {
          double* y_coord = (double *) ML_allocate(sizeof(double) * (Nghost+n));

          for (int i = 0 ; i < n ; ++i)
            tmp[i * NumPDEEqns_] = in_y_coord[i];

          ML_exchange_bdry(&tmp[0],AAA->getrow->pre_comm, NumPDEEqns_ * n,
                           AAA->comm, ML_OVERWRITE,NULL);

          for (int i = 0 ; i < n + Nghost ; ++i)
            y_coord[i] = tmp[i * NumPDEEqns_];

          coordy = y_coord;
        }

      if (in_z_coord)
        {
          double* z_coord = (double *) ML_allocate(sizeof(double) * (Nghost+n));

          for (int i = 0 ; i < n ; ++i)
            tmp[i * NumPDEEqns_] = in_z_coord[i];

          ML_exchange_bdry(&tmp[0],AAA->getrow->pre_comm, NumPDEEqns_ * n,
                           AAA->comm, ML_OVERWRITE,NULL);

          for (int i = 0 ; i < n + Nghost ; ++i)
            z_coord[i] = tmp[i * NumPDEEqns_];

          coordz = z_coord;
        }

    } // if (!(in_x_coord == 0 && in_y_coord == 0 && in_z_coord == 0))

  return(0);
 }

// ================================================ ====== ==== ==== == =
// Copied from ml_agg_genP.c
static void ML_Init_Aux(ML_Operator* A, Teuchos::ParameterList &List) {
  int i, j, n, count, num_PDEs, BlockRow, BlockCol;
  double threshold;
  int* columns;
  double* values;
  int allocated, entries = 0;
  int N_dimensions;
  int DiagID;
  double DiagValue;
  int** filter;
  double dist;
  double *LaplacianDiag;
  int     Nghost;


  // Boundary exchange the coords
  double *x_coord=0, *y_coord=0, *z_coord=0;
  RefMaxwell_SetupCoordinates(A,List,x_coord,y_coord,z_coord);
  int dim=(x_coord!=0) + (y_coord!=0) + (z_coord!=0);

  /* Sanity Checks */
  if(dim == 0 || ((!x_coord && (y_coord || z_coord)) || (x_coord && !y_coord && z_coord))){
    std::cerr<<"Error: Coordinates not defined.  This is necessary for aux aggregation (found "<<dim<<" coordinates).\n";
    exit(-1);
  }

  num_PDEs = A->num_PDEs;
  N_dimensions = dim;
  threshold = A->aux_data->threshold;

  ML_Operator_AmalgamateAndDropWeak(A, num_PDEs, 0.0);
  n = A->invec_leng;
  Nghost = ML_CommInfoOP_Compute_TotalRcvLength(A->getrow->pre_comm);

  LaplacianDiag = (double *) ML_allocate((A->getrow->Nrows+Nghost+1)*
                                         sizeof(double));

  filter = (int**) ML_allocate(sizeof(int*) * n);

  allocated = 128;
  columns = (int *)    ML_allocate(allocated * sizeof(int));
  values  = (double *) ML_allocate(allocated * sizeof(double));

  for (i = 0 ; i < n ; ++i) {
    BlockRow = i;
    DiagID = -1;
    DiagValue = 0.0;

    ML_get_matrix_row(A,1,&i,&allocated,&columns,&values, &entries,0);

    for (j = 0; j < entries; j++) {
      BlockCol = columns[j];
      if (BlockRow != BlockCol) {
        dist = 0.0;

        switch (N_dimensions) {
        case 3:
          dist += (z_coord[BlockRow] - z_coord[BlockCol]) * (z_coord[BlockRow] - z_coord[BlockCol]);
        case 2:
          dist += (y_coord[BlockRow] - y_coord[BlockCol]) * (y_coord[BlockRow] - y_coord[BlockCol]);
        case 1:
          dist += (x_coord[BlockRow] - x_coord[BlockCol]) * (x_coord[BlockRow] - x_coord[BlockCol]);
        }

        if (dist == 0.0) {
          printf("node %d = %e ", i, x_coord[BlockRow]);
          if (N_dimensions > 1) printf(" %e ", y_coord[BlockRow]);
          if (N_dimensions > 2) printf(" %e ", z_coord[BlockRow]);
          printf("\n");
          printf("node %d = %e ", j, x_coord[BlockCol]);
          if (N_dimensions > 1) printf(" %e ", y_coord[BlockCol]);
          if (N_dimensions > 2) printf(" %e ", z_coord[BlockCol]);
          printf("\n");
          printf("Operator has inlen = %d and outlen = %d\n",
                 A->invec_leng, A->outvec_leng);
        }

        dist = 1.0 / dist;
        DiagValue += dist;
      }
      else if (columns[j] == i) {
        DiagID = j;
      }
    }

    if (DiagID == -1) {
      fprintf(stderr, "ERROR: matrix has no diagonal!\n"
              "ERROR: (file %s, line %d)\n",
              __FILE__, __LINE__);
      exit(EXIT_FAILURE);
    }
    LaplacianDiag[BlockRow] = DiagValue;
  }
  if ( A->getrow->pre_comm != NULL )
     ML_exchange_bdry(LaplacianDiag,A->getrow->pre_comm,A->getrow->Nrows,
                      A->comm, ML_OVERWRITE,NULL);


  for (i = 0 ; i < n ; ++i) {
    BlockRow = i;

    ML_get_matrix_row(A,1,&i,&allocated,&columns,&values, &entries,0);

    for (j = 0; j < entries; j++) {
      BlockCol = columns[j];
      if (BlockRow != BlockCol) {
        dist = 0.0;
        switch (N_dimensions) {
        case 3:
          dist += (z_coord[BlockRow] - z_coord[BlockCol]) * (z_coord[BlockRow] - z_coord[BlockCol]);
        case 2:
          dist += (y_coord[BlockRow] - y_coord[BlockCol]) * (y_coord[BlockRow] - y_coord[BlockCol]);
        case 1:
          dist += (x_coord[BlockRow] - x_coord[BlockCol]) * (x_coord[BlockRow] - x_coord[BlockCol]);
        }

        dist = 1.0 / dist;
        values[j] = dist;
      }
    }

    count = 0;
    for (j = 0 ; j < entries ; ++j) {
      if (  (i != columns[j]) &&
            (values[j]*values[j] <
       LaplacianDiag[BlockRow]*LaplacianDiag[columns[j]]*threshold*threshold)){
        columns[count++] = columns[j];
      }
    }

    /* insert the rows */
    filter[BlockRow] = (int*) ML_allocate(sizeof(int) * (count + 1));
    filter[BlockRow][0] = count;

    for (j = 0 ; j < count ; ++j)
      filter[BlockRow][j + 1] = columns[j];

  }

  ML_free(columns);
  ML_free(values);

  ML_free(LaplacianDiag);

  ML_Operator_UnAmalgamateAndDropWeak(A, num_PDEs, 0.0);

  A->aux_data->aux_func_ptr  = A->getrow->func_ptr;
  A->getrow->func_ptr = ML_Aux_Getrow;
  A->aux_data->filter = filter;
  A->aux_data->filter_size = n;

  // Cleanup
  ML_free(x_coord);
  ML_free(y_coord);
  ML_free(z_coord);

}

// ================================================ ====== ==== ==== == =
// Copied from ml_agg_genP.c
static void ML_Finalize_Aux(ML_Operator *A)
{
  int i;
  A->getrow->func_ptr = A->aux_data->aux_func_ptr;
  A->aux_data->aux_func_ptr = 0;
  for (i = 0 ; i < A->aux_data->filter_size ; ++i)
    ML_free((A->aux_data->filter[i]));
  ML_free(A->aux_data->filter);
}


// ================================================ ====== ==== ==== == =
int ML_Epetra::RefMaxwell_Aggregate_Nodes(const Epetra_CrsMatrix & A, Teuchos::ParameterList & List, ML_Comm * ml_comm, std::string PrintMsg,
                                          ML_Aggregate_Struct *& MLAggr,ML_Operator *&P, int &NumAggregates){

  /* Output level */
  bool verbose, very_verbose;
  int OutputLevel = List.get("ML output", -47);
  if(OutputLevel == -47) OutputLevel = List.get("output", 1);
  if(OutputLevel>=15) very_verbose=verbose=true;
  if(OutputLevel > 5) {very_verbose=false;verbose=true;}
  else very_verbose=verbose=false;

  /* Wrap A in a ML_Operator */
  ML_Operator* A_ML = ML_Operator_Create(ml_comm);
  ML_Operator_WrapEpetraCrsMatrix(const_cast<Epetra_CrsMatrix*>(&A),A_ML);

  /* Pull Teuchos Options */
  std::string CoarsenType  = List.get("aggregation: type", "Uncoupled");
  double Threshold         = List.get("aggregation: threshold", 0.0);
  double RowSum_Threshold  = List.get("aggregation: rowsum threshold", -1.0);
  int    NodesPerAggr = List.get("aggregation: nodes per aggregate",
                                  ML_Aggregate_Get_OptimalNumberOfNodesPerAggregate());
  bool UseAux         = List.get("aggregation: aux: enable",false);
  double AuxThreshold = List.get("aggregation: aux: threshold",0.0);
  int  MaxAuxLevels   = List.get("aggregation: aux: max levels",10);


  ML_Aggregate_Create(&MLAggr);
  ML_Aggregate_Set_MaxLevels(MLAggr, 2);
  ML_Aggregate_Set_StartLevel(MLAggr, 0);
  ML_Aggregate_Set_Threshold(MLAggr, Threshold);
  ML_Aggregate_Set_RowSum_Threshold(MLAggr, RowSum_Threshold);
  ML_Aggregate_Set_MaxCoarseSize(MLAggr,1);
  MLAggr->cur_level = 0;
  ML_Aggregate_Set_Reuse(MLAggr);
  MLAggr->keep_agg_information = 1;
  P = ML_Operator_Create(ml_comm);

  /* Process Teuchos Options */
  if (CoarsenType == "Uncoupled")
    ML_Aggregate_Set_CoarsenScheme_Uncoupled(MLAggr);
  else if (CoarsenType == "Uncoupled-MIS"){
    ML_Aggregate_Set_CoarsenScheme_UncoupledMIS(MLAggr);
  }
  else if (CoarsenType == "METIS"){
    ML_Aggregate_Set_CoarsenScheme_METIS(MLAggr);
    ML_Aggregate_Set_NodesPerAggr(0, MLAggr, 0, NodesPerAggr);
  }/*end if*/
  else {
    if(!A.Comm().MyPID()) printf("%s Unsupported (1,1) block aggregation type(%s), resetting to uncoupled-mis\n",PrintMsg.c_str(),CoarsenType.c_str());
    ML_Aggregate_Set_CoarsenScheme_UncoupledMIS(MLAggr);
  }

  /* Setup Aux Data */
  if(UseAux) {
    A_ML->aux_data->enable=1;
    A_ML->aux_data->threshold=AuxThreshold;
    A_ML->aux_data->max_level=MaxAuxLevels;
    ML_Init_Aux(A_ML,List);
    if(verbose && !A.Comm().MyPID()) {
      printf("%s Using auxiliary matrix\n",PrintMsg.c_str());
      printf("%s aux threshold = %e\n",PrintMsg.c_str(),A_ML->aux_data->threshold);
    }
  }

  /* Aggregate Nodes */
  int printlevel=ML_Get_PrintLevel();
  if(verbose) ML_Set_PrintLevel(10);
  NumAggregates = ML_Aggregate_Coarsen(MLAggr,A_ML, &P, ml_comm);
  if(verbose) ML_Set_PrintLevel(printlevel);

  if (NumAggregates == 0){
    std::cerr << "Found 0 aggregates, perhaps the problem is too small." << std::endl;
    ML_CHK_ERR(-2);
  }/*end if*/
  else if(very_verbose) printf("[%d] %s %d aggregates created invec_leng=%d\n",A.Comm().MyPID(),PrintMsg.c_str(),NumAggregates,P->invec_leng);

  if(verbose){
    int globalAggs=0;
    A.Comm().SumAll(&NumAggregates,&globalAggs,1);
    if(!A.Comm().MyPID()) {
      printf("%s Aggregation threshold = %e\n",PrintMsg.c_str(),Threshold);
      printf("%s Global aggregates     = %d\n",PrintMsg.c_str(),globalAggs);

    }
  }

  /* Cleanup */
  ML_qr_fix_Destroy();
  if(UseAux) ML_Finalize_Aux(A_ML);
  ML_Operator_Destroy(&A_ML);

  return 0;
}

// ================================================ ====== ==== ==== == =
// Build the edge nullspace
Epetra_MultiVector* ML_Epetra::Build_Edge_Nullspace(const Epetra_CrsMatrix & D0_Clean_Matrix,  const Teuchos::ArrayRCP<int> BCedges, Teuchos::ParameterList & List_,bool verbose_)
{
  Epetra_MultiVector *nullspace;
  double ** d_coords;
  const Epetra_Map & NodeMap = D0_Clean_Matrix.DomainMap();
  const Epetra_Map & EdgeMap = D0_Clean_Matrix.RangeMap();
  const Epetra_Comm & Comm   = D0_Clean_Matrix.Comm();

  /* Check the List - Do we have a nullspace pre-provided? */
  std::string nulltype=List_.get("null space: type","default vectors");
  double* nullvecs=List_.get("null space: vectors",(double*)0);
  int dim=List_.get("null space: dimension",0);
  if (nulltype=="pre-computed" && nullvecs && (dim==2 || dim==3)){
    /* Build a multivector out of it */
    if(verbose_ && !Comm.MyPID()) printf("Using pre-computed nullspace\n");
    int Ne=EdgeMap.NumMyElements();
    d_coords=new double*[dim];
    d_coords[0]=nullvecs;
    d_coords[1]=&nullvecs[Ne];
    if(dim==3) d_coords[2]=&nullvecs[2*Ne];
    nullspace=new Epetra_MultiVector(View,EdgeMap,d_coords,dim);
  }
  else{
    if(verbose_ && !Comm.MyPID()) printf("Building nullspace from scratch\n");
    /* Pull the (nodal) coordinates from Teuchos */
    double * xcoord=List_.get("x-coordinates",(double*)0);
    double * ycoord=List_.get("y-coordinates",(double*)0);
    double * zcoord=List_.get("z-coordinates",(double*)0);
    dim=(xcoord!=0) + (ycoord!=0) + (zcoord!=0);

    /* Sanity Checks */
    TEUCHOS_TEST_FOR_EXCEPT_MSG(
      dim == 0 || ((!xcoord && (ycoord || zcoord)) || (xcoord && !ycoord && zcoord)),
      "Error: Coordinates not defined and no nullspace is provided.  One of these are *necessary* for the EdgeMatrixFreePreconditioner (found "<<dim<<" coordinates).\n"
      );

#ifdef OLD
    /* Normalize */
    double d1 = sqrt(ML_gdot(NodeMap.NumMyElements(), xcoord, xcoord, ml_comm_));
    for (int i = 0; i < NodeMap.NumMyElements(); i++) xcoord[i] /= d1;
    d1 = sqrt(ML_gdot(NodeMap.NumMyElements(), ycoord, ycoord, ml_comm_));
    for (int i = 0; i < NodeMap.NumMyElements(); i++) ycoord[i] /= d1;
    if (dim==3) {
      d1 = sqrt(ML_gdot(NodeMap.NumMyElements(), zcoord, zcoord, ml_comm_));
      for (int i = 0; i < NodeMap.NumMyElements(); i++) zcoord[i] /= d1;
    }
#endif

    /* Build the MultiVector */
    d_coords=new double* [dim];
    d_coords[0]=xcoord; d_coords[1]=ycoord;
    if(dim==3) d_coords[2]=zcoord;
    Epetra_MultiVector n_coords(View,NodeMap,d_coords,dim);

    // CMS: We're removing the normalization here.

    /* Build the Nullspace */
    nullspace=new Epetra_MultiVector(EdgeMap,dim,true);
    D0_Clean_Matrix.Multiply(false,n_coords,*nullspace);
  }

  /* Nuke the BC edges */
  for(int j=0;j<dim;j++)
    for(int i=0;i<BCedges.size();i++)
      (*nullspace)[j][BCedges[i]]=0;

  /* Cleanup */
  delete [] d_coords ;
  return nullspace;
}/*end BuildNullspace*/



#endif
