#include "ml_config.h"
#include "ml_EdgeMatrixFreePreconditioner.h"
#include "ml_MultiLevelPreconditioner.h"
#include "ml_epetra.h"
#include "ml_epetra_utils.h"
#include "ml_mat_formats.h"
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS)

//mucho hax
void cms_residual_check(const char * tag, const Epetra_Operator * op,const Epetra_MultiVector& rhs, const Epetra_MultiVector& lhs);
double cms_compute_residual(const Epetra_Operator * op,const Epetra_MultiVector& rhs, const Epetra_MultiVector& lhs);


#define NO_OUTPUT

#ifdef NO_OUTPUT
#define MVOUT2(x,y,z) ;
#define MVOUT(x,y) ;
#define IVOUT(x,y) ;
#define Epetra_CrsMatrix_Print(x,y) ;
#define ML_Matrix_Print(w,x,y,z) ;
#else
extern void Epetra_CrsMatrix_Print(const Epetra_CrsMatrix& A, char* of);
extern void MVOUT (const Epetra_MultiVector & A, char * of);
extern void IVOUT(const Epetra_IntVector & A, char * of);
extern void MVOUT2(const Epetra_MultiVector & A,char* pref,int idx);
extern void ML_Matrix_Print(ML_Operator *ML,const Epetra_Comm &Comm,const Epetra_Map &Map, char *fname);
#endif

/* Turn this on the use the heavyweight wraps rather than the lightweight ones.
The lighter wraps are a bit more fragile and require that the Epetra_CrsMatrix
stick around through the life of the wrapped matrix */
//#define USE_HEAVYWEIGHT_WRAPS


// ================================================ ====== ==== ==== == =
/* This function does a "view" getrow in an ML_Operator.  This is intended to be
used in a ML_CSR_Matrix to Epetra_CrsMatrix (view) translator.  Inlined for
speed. */
inline void CSR_getrow_view(ML_Operator *M, int row, int *ncols,int **cols, double**vals){
  struct ML_CSR_MSRdata* M_= (struct ML_CSR_MSRdata*)ML_Get_MyGetrowData(M);
  *ncols=M_->rowptr[row+1]-M_->rowptr[row];
  *cols=&M_->columns[M_->rowptr[row]];
  *vals=&M_->values[M_->rowptr[row]];  
}/*end CSR_getrow_view*/


int CSR_getrow_ones(ML_Operator *data, int N_requested_rows, int requested_rows[],
   int allocated_space, int columns[], double values[], int row_lengths[])
{
   register int    *bindx, j;
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
ML_Epetra::EdgeMatrixFreePreconditioner::EdgeMatrixFreePreconditioner(const Epetra_Operator_With_MatMat & Operator, const Epetra_Vector& Diagonal, const Epetra_CrsMatrix & D0_Matrix,const Epetra_CrsMatrix & D0_Clean_Matrix,const Epetra_CrsMatrix &TMT_Matrix,const int* BCedges, const int numBCedges, const Teuchos::ParameterList &List,const bool ComputePrec):
  ML_Preconditioner(),Operator_(&Operator),D0_Matrix_(&D0_Matrix),D0_Clean_Matrix_(&D0_Clean_Matrix),TMT_Matrix_(&TMT_Matrix),BCedges_(BCedges),numBCedges_(numBCedges),Prolongator_(0),InvDiagonal_(0),CoarseMatrix(0),CoarsePC(0),Smoother_(0),verbose_(false)
{
  /* Set the Epetra Goodies */
  Comm_ = &(Operator_->Comm());
  if(verbose_ && !Comm_->MyPID()) printf("EMFP: Constructor Commencing\n");

  EdgeDomainMap_ = &(Operator_->OperatorDomainMap());
  EdgeRangeMap_ = &(Operator_->OperatorRangeMap());
  NodeDomainMap_ = &(TMT_Matrix_->OperatorDomainMap());
  NodeRangeMap_ = &(TMT_Matrix_->OperatorRangeMap());
  
  List_=List;
  Label_=strdup("ML edge matrix-free preconditioner");
  InvDiagonal_ = new Epetra_Vector(Diagonal);  
  if(ComputePrec) ML_CHK_ERRV(ComputePreconditioner());
}/*end constructor*/

// ================================================ ====== ==== ==== == =   
ML_Epetra::EdgeMatrixFreePreconditioner::~EdgeMatrixFreePreconditioner(){
  DestroyPreconditioner();
}/*end destructor*/


// ================================================ ====== ==== ==== == = 
// Computes the preconditioner
int ML_Epetra::EdgeMatrixFreePreconditioner::ComputePreconditioner(const bool CheckFiltering)
{
  Teuchos::ParameterList dummy, ListCoarse;
  ListCoarse=List_.get("edge matrix free: coarse",dummy);
  
  /* ML Communicator */
  ML_Comm_Create(&ml_comm_);

  /* Parameter List Options */
  int OutputLevel = List_.get("ML output", -47);
  int SmootherSweeps = List_.get("smoother: sweeps (level 0)", 0);
  MaxLevels = List_.get("max levels",10); 
  
  if (OutputLevel == -47) OutputLevel = List_.get("output", 1);
  if(OutputLevel > 5) verbose_=true;
  else verbose_=false;

  num_cycles  = List_.get("cycle applications",1);
  ML_Set_PrintLevel(OutputLevel);

  /* Sanity Checking*/
  int OperatorDomainPoints =  OperatorDomainMap().NumGlobalPoints();
  int OperatorRangePoints =  OperatorRangeMap().NumGlobalPoints();
  if (OperatorDomainPoints != OperatorRangePoints)
    ML_CHK_ERR(-1); // only square matrices

  /* Invert non-zeros on the diagonal */
  if(SmootherSweeps){
    for (int i = 0; i < InvDiagonal_->MyLength(); ++i)
      if ((*InvDiagonal_)[i] != 0.0)
        (*InvDiagonal_)[i] = 1.0 / (*InvDiagonal_)[i];   
    double nrm;
    InvDiagonal_->Norm2(&nrm);
    if(verbose_ && !Comm_->MyPID()) printf("Inverse Diagonal Norm = %6.4e\n",nrm);
  }/*end if*/
  
  /* Do the eigenvalue estimation for Chebyshev */
  if(SmootherSweeps) {
    if(verbose_ && !Comm_->MyPID()) printf("EMFP: Doing Smoother Setup\n");
    ML_CHK_ERR(SetupSmoother());
  }/*end if*/


  if(MaxLevels > 1) {  
    /* Build the Nullspace */
    if(verbose_ && !Comm_->MyPID()) printf("EMFP: Building Nullspace\n");
    Epetra_MultiVector *nullspace=BuildNullspace();
    if(!nullspace) ML_CHK_ERR(-1);
    
    MVOUT(*nullspace,"nullspace.dat");
    
    /* Build the prolongator */
    if(verbose_ && !Comm_->MyPID()) printf("EMFP: Building Prolongator\n");
    ML_CHK_ERR(BuildProlongator(*nullspace));
    
    /* DEBUG: Output matrices */
    Epetra_CrsMatrix_Print(*Prolongator_,"prolongator.dat");
    
    /* Form the coarse matrix */
    if(verbose_ && !Comm_->MyPID()) printf("EMFP: Building Coarse Matrix\n");
    ML_CHK_ERR(FormCoarseMatrix());

    /* DEBUG: Output matrices */
    Epetra_CrsMatrix_Print(*CoarseMatrix,"coarsemat.dat");
    
    /* Setup Preconditioner on Coarse Matrix */
    if(verbose_ && !Comm_->MyPID()) printf("EMFP: Building Coarse Precond\n");
    CoarsePC = new MultiLevelPreconditioner(*CoarseMatrix,ListCoarse);
    if(!CoarsePC) ML_CHK_ERR(-2);
  
    /* Clean Up */
    delete nullspace;
  }/*end if*/
    
  return 0;
}/*end ComputePreconditioner*/


// ================================================ ====== ==== ==== == = 
// Setup the Smoother
int ML_Epetra::EdgeMatrixFreePreconditioner::SetupSmoother()
{
  /* Variables */
  double lambda_min = 0.0;
  double lambda_max = 0.0;
  Teuchos::ParameterList IFPACKList;  

  /* Parameter-list Options */
  int PolynomialDegree = List_.get("smoother: degree", 3);
  PolynomialDegree = List_.get("smoother: sweeps", 3);// override if need be  
  int MaximumIterations = List_.get("eigen-analysis: max iters", 10);
  string EigenType_ = List_.get("eigen-analysis: type", "cg");
  double boost = List_.get("eigen-analysis: boost for lambda max", 1.0);


  MVOUT(*InvDiagonal_,"inv_diagonal.dat");
  
  /* Do the eigenvalue estimation*/
  if (EigenType_ == "power-method"){   
    ML_CHK_ERR(Ifpack_Chebyshev::PowerMethod(*Operator_,*InvDiagonal_,MaximumIterations,lambda_max));
    lambda_min=lambda_max/30.0;
  }/*end if*/
  else if(EigenType_ == "cg"){    
    ML_CHK_ERR(Ifpack_Chebyshev::CG(*Operator_,*InvDiagonal_,MaximumIterations,lambda_min,lambda_max));
  }/*end else if*/
  else
    ML_CHK_ERR(-1); // not recognized

  double alpha = List_.get("chebyshev: alpha",30.0001);
  lambda_min=lambda_max / alpha;

#ifndef NO_OUTPUT
  FILE *f=fopen("lambda_max.dat","w");
  fprintf(f,"%22.16e\n",lambda_max);
  fclose(f);
#endif
  
  /* Setup the Smoother's List*/
  IFPACKList.set("chebyshev: min eigenvalue", lambda_min);
  IFPACKList.set("chebyshev: max eigenvalue", boost * lambda_max);
  IFPACKList.set("chebyshev: ratio eigenvalue",alpha);
  IFPACKList.set("chebyshev: operator inv diagonal", InvDiagonal_);
  IFPACKList.set("chebyshev: degree", PolynomialDegree);
  IFPACKList.set("chebyshev: zero starting solution",false);

  if(verbose_ && !Comm_->MyPID()) printf("Chebyshev Smoother: lmin/lmax %6.4e/%6.4e\n",lambda_min,lambda_max);//DEBUG

  //NTS: Need to create two of these lists, one to use in the first iteration
  // (with zero starting solution set to true) and another to use when it's set
  // to false.
  
  /* Build the Smoother */
  Smoother_ = new Ifpack_Chebyshev(Operator_);
  if (Smoother_ == 0) ML_CHK_ERR(-1); 
  ML_CHK_ERR(Smoother_->SetParameters(IFPACKList));
  ML_CHK_ERR(Smoother_->Initialize());
  ML_CHK_ERR(Smoother_->Compute());

  return 0;
}/*end SetupSmoother */



// ================================================ ====== ==== ==== == = 
// Build the edge nullspace
Epetra_MultiVector * ML_Epetra::EdgeMatrixFreePreconditioner::BuildNullspace()
{
  /* Pull the coordinates from Teuchos */
  double * xcoord=List_.get("x-coordinates",(double*)0);
  double * ycoord=List_.get("y-coordinates",(double*)0);
  double * zcoord=List_.get("z-coordinates",(double*)0);
  dim=(xcoord!=0) + (ycoord!=0) + (zcoord!=0);
  
  /* Sanity Checks */
  if(dim == 0 || (!xcoord && (ycoord || zcoord) || (xcoord && !ycoord && zcoord))){
    cerr<<"Error: Coordinates not defined.  This is *necessary* for the EdgeMatrixFreePreconditioner.\n";
    return 0;
  }/*end if*/
  if(verbose_ && !Comm_->MyPID()) printf("BuildNullspace: Pulling %d vectors\n",dim);
  
  /* Build the MultiVector */
  double ** d_coords=new double* [3];
  d_coords[0]=xcoord; d_coords[1]=ycoord; d_coords[2]=zcoord;
  Epetra_MultiVector e_coords(View,*NodeDomainMap_,d_coords,dim);

  MVOUT(e_coords,"coords.dat");

  
  /* Build the Nullspace */
  Epetra_MultiVector * nullspace=new Epetra_MultiVector(*EdgeDomainMap_,dim,false);  
  D0_Clean_Matrix_->Multiply(false,e_coords,*nullspace);  

  /* Nuke the BC edges */
  for(int j=0;j<dim;j++)
    for(int i=0;i<numBCedges_;i++)
      (*nullspace)[j][BCedges_[i]]=0;
  
  /* Cleanup */
  delete [] d_coords ;
  return nullspace;
}/*end BuildNullspace*/




// ================================================ ====== ==== ==== == = 
//! Build the edge-to-vector-node prolongator described in Bochev, Hu, Siefert and Tuminaro (2006).
int ML_Epetra::EdgeMatrixFreePreconditioner::BuildProlongator(const Epetra_MultiVector & nullspace)
{

  /* Wrap TMT_Matrix in a ML_Operator */
  ML_Operator* TMT_ML = ML_Operator_Create(ml_comm_);
  ML_Operator_WrapEpetraCrsMatrix((Epetra_CrsMatrix*)TMT_Matrix_,TMT_ML);

  /* Pull Teuchos Options */
  string CoarsenType = List_.get("aggregation: type", "Uncoupled");
  double Threshold   = List_.get("aggregation: threshold", 0.0);  
  int    NodesPerAggr = List_.get("aggregation: nodes per aggregate", 
                                  ML_Aggregate_Get_OptimalNumberOfNodesPerAggregate());

  /* Setup the Aggregation */
  if(verbose_ && !Comm_->MyPID()) printf("EMFP: Building aggregates\n");
  ML_Aggregate_Struct * MLAggr;
  ML_Aggregate_Create(&MLAggr);
  ML_Aggregate_Set_MaxLevels(MLAggr, 2);
  ML_Aggregate_Set_StartLevel(MLAggr, 0);
  ML_Aggregate_Set_Threshold(MLAggr, Threshold);
  MLAggr->cur_level = 0;
  ML_Aggregate_Set_Reuse(MLAggr);
  MLAggr->keep_agg_information = 1;  
  ML_Operator *P = ML_Operator_Create(ml_comm_);
  
  /* Process Teuchos Options */
  if (CoarsenType == "Uncoupled")  MLAggr->coarsen_scheme = ML_AGGR_UNCOUPLED;
  else if (CoarsenType == "METIS"){
    MLAggr->coarsen_scheme = ML_AGGR_METIS;
    ML_Aggregate_Set_NodesPerAggr(0, MLAggr, 0, NodesPerAggr);
  }/*end if*/
  else if(verbose_ && !Comm_->MyPID()) {printf("aggregation: type = %d\n",CoarsenType.c_str()); ML_CHK_ERR(-1);}

  /* Aggregate Nodes */
  int NumAggregates = ML_Aggregate_Coarsen(MLAggr, TMT_ML, &P, ml_comm_);
  if (NumAggregates == 0){
    cerr << "Found 0 aggregates, perhaps the problem is too small." << endl;
    ML_CHK_ERR(-2);
  }/*end if*/
  else if(verbose_) printf("[%d] EMFP: %d aggregates created invec_leng=%d\n",Comm_->MyPID(),NumAggregates,P->invec_leng);

  
  if(verbose_) printf("[%d] Num Aggregates = %d\n",Comm_->MyPID(),NumAggregates);
  if(P==0) {fprintf(stderr,"ERROR: No tentative prolongator found\n");ML_CHK_ERR(-5);}
  
#ifndef NO_OUTPUT
  /* DEBUG: Dump aggregates */ 
  Epetra_IntVector AGG(View,*NodeDomainMap_,MLAggr->aggr_info[0]);
  IVOUT(AGG,"agg.dat");  
  /* DEBUG: Testing Stuff */
  ML_Operator_Print(P,"p.dat");
  Epetra_CrsMatrix_Print(*D0_Matrix_,"test.dat");
  Epetra_CrsMatrix_Print(*D0_Clean_Matrix_,"test_clean.dat");
#endif
  
  /* Create wrapper to do abs(T) */
  // NTS: Assume D0 has already been reindexed by now.
  if(verbose_ && !Comm_->MyPID()) printf("EMFP: abs(T) prewrap\n");
  ML_Operator* AbsD0_ML = ML_Operator_Create(ml_comm_);
  ML_CHK_ERR(ML_Operator_WrapEpetraCrsMatrix((Epetra_CrsMatrix*)D0_Matrix_,AbsD0_ML,verbose_));  
  
#ifdef USE_HEAVYWEIGHT_WRAPS
  ML_Operator_Set_Getrow(AbsD0_ML,AbsD0_ML->outvec_leng,ML_Epetra_CrsMatrix_get_one_row);
#else
  ML_Operator_Set_Getrow(AbsD0_ML,AbsD0_ML->outvec_leng,CSR_getrow_ones);
#endif
  
  /* Form abs(T) * P_n */
  if(verbose_ && !Comm_->MyPID()) printf("EMFP: Building abs(T) * P_n\n");
  ML_Operator* AbsD0P = ML_Operator_Create(ml_comm_);   
  ML_2matmult(AbsD0_ML,P,AbsD0P, ML_CSR_MATRIX);
  
  /* Wrap P_n into Epetra-land */
  if(verbose_ && !Comm_->MyPID()) printf("EMFP: Wrapping to PSparse\n");
  Epetra_CrsMatrix *Psparse;
  Epetra_CrsMatrix_Wrap_ML_Operator(AbsD0P,*Comm_,*EdgeRangeMap_,&Psparse,Copy,0);
  
  /* Nuke the rows in Psparse */
  Apply_BCsToMatrixRows(BCedges_,numBCedges_,*Psparse);
  
  /* DEBUG: output*/
#ifndef NO_OUTPUT
  Epetra_CrsMatrix_Print(*Psparse,"psparse.dat");
  MVOUT(nullspace,"nullspace.dat");
#endif
  
  /* Build the DomainMap of the new operator*/
  const Epetra_Map & FineColMap = Psparse->ColMap();
  CoarseMap_=new Epetra_Map(-1,NumAggregates*dim,0,*Comm_);
  
  /* Allocate the Prolongator_ */
  if(verbose_ && !Comm_->MyPID()) printf("EMFP: Building Prolongator\n");
  Prolongator_=new Epetra_CrsMatrix(Copy,*EdgeRangeMap_,0);
  int ne1, *idx1, *idx2;
  idx2=new int [dim*AbsD0P->max_nz_per_row];
  double *vals1, *vals2;
  vals2=new double[dim*AbsD0P->max_nz_per_row];
  int nonzeros;
  
  for(int i=0;i<Prolongator_->NumMyRows();i++){
    Psparse->ExtractMyRowView(i,ne1,vals1,idx1);
    nonzeros=0;
    for(int j=0;j<ne1;j++) nonzeros+=vals1[j]>0;

    for(int j=0;j<ne1;j++){
      for(int k=0;k<dim;k++) {
        idx2[j*dim+k]=FineColMap.GID(idx1[j])*dim+k;
        //FIX: This works only because there's an implicit linear mapping which
        //we're exploiting.
        if(idx2[j*dim+k]==-1) printf("[%d] ERROR: idx1[j]=%d / idx1[j]*dim+k=%d does not have a GID!\n",Comm_->MyPID(),idx1[j],idx1[j]*dim+k);
        if(vals1[j]==0 ) vals2[j*dim+k]=0;
        else vals2[j*dim+k]= nullspace[k][i] / nonzeros;  
      }/*end for*/
    }/*end for*/
    Prolongator_->InsertGlobalValues(EdgeRangeMap_->GID(i),dim*ne1,vals2,idx2);
  }/*end for*/
  
  
  /* FillComplete / OptimizeStorage for Prolongator*/
  if(verbose_ && !Comm_->MyPID()) printf("EMFP: Optimizing Prolongator\n");
  Prolongator_->FillComplete(*CoarseMap_,*EdgeRangeMap_);
  Prolongator_->OptimizeStorage();
  
  /* Post-wrapping to convert to ML indexing */
#ifdef HAVE_ML_EPETRAEXT
  Prolongator_ = dynamic_cast<Epetra_CrsMatrix*>(ModifyEpetraMatrixColMap(*Prolongator_,ProlongatorColMapTrans_,"Prolongator",verbose_));
#endif
  
  /* Cleanup */
  if(verbose_ && !Comm_->MyPID()) printf("EMFP: BuildProlongator Cleanup\n");  
  ML_Aggregate_Destroy(&MLAggr);
  ML_Operator_Destroy(&TMT_ML);
  ML_Operator_Destroy(&P);
  ML_Operator_Destroy(&AbsD0_ML);
  ML_Operator_Destroy(&AbsD0P);

  delete Psparse;
  delete [] idx2;
  delete [] vals2;
  return 0;
}/*end BuildProlongator_*/
  



// ================================================ ====== ==== ==== == = 
// Forms the coarse matrix, given the prolongator
int  ML_Epetra::EdgeMatrixFreePreconditioner::FormCoarseMatrix()
{
  ML_Operator *R= ML_Operator_Create(ml_comm_);
  ML_Operator *P= ML_Operator_Create(ml_comm_);
  ML_Operator *Temp_ML;
  CoarseMat_ML = ML_Operator_Create(ml_comm_);
  CoarseMat_ML->data_destroy=free;
  
  /* Build ML_Operator version of Prolongator_, Restriction Operator */
  if(verbose_ && !Comm_->MyPID()) printf("EMFP: Prolongator Prewrap\n");
  ML_CHK_ERR(ML_Operator_WrapEpetraCrsMatrix(Prolongator_,P,verbose_));
  P->num_rigid=P->num_PDEs=3;
  
  if(verbose_ && !Comm_->MyPID()) printf("EMFP: Prolongator Transpose\n");
  //NTS: ML_CHK_ERR won't work on this: it returns 1
  ML_Operator_Transpose_byrow(P, R);

  /* Do the A*P */
  if(verbose_ && !Comm_->MyPID()) printf("EMFP: AP\n");
  ML_CHK_ERR(Operator_->MatrixMatrix_Multiply(*Prolongator_,ml_comm_,&Temp_ML));  
  
  /* DEBUG: output*/
#ifndef NO_OUTPUT
  ML_Matrix_Print(Temp_ML,*Comm_,*EdgeRangeMap_,"coarse_temp.dat");
#endif

  /* Do R * AP */
  if(verbose_ && !Comm_->MyPID()) printf("EMFP: RAP\n");
  R->num_rigid=R->num_PDEs=3;
  //  ML_2matmult(R, Temp_ML,CoarseMat_ML,ML_CSR_MATRIX);
  ML_2matmult_block(R, Temp_ML,CoarseMat_ML,ML_CSR_MATRIX);
     
  Epetra_CrsMatrix_Wrap_ML_Operator(CoarseMat_ML,*Comm_,*CoarseMap_,&CoarseMatrix); 
  
#ifndef NO_OUTPUT
  ML_Matrix_Print(CoarseMat_ML,*Comm_,*EdgeRangeMap_,"coarsemat.dat");  
#endif
 
  /* Fix OAZ issues - is this needed?*/
  //  Remove_Zeroed_Rows(*CoarseMatrix);

#ifndef NO_OUTPUT
  Epetra_CrsMatrix_Print(*CoarseMatrix,"coarsemat0.dat");
#endif
  
  /* Cleanup */
  ML_Operator_Destroy(&P);
  ML_Operator_Destroy(&R);
  ML_Operator_Destroy(&Temp_ML);
  return 0;
}/*end FormCoarseMatrix*/

// ================================================ ====== ==== ==== == = 
// Print the individual operators in the multigrid hierarchy.
void ML_Epetra::EdgeMatrixFreePreconditioner::Print(const char *whichHierarchy)
{
  if(CoarseMatrix) CoarseMatrix->Print(cout);
}/*end Print*/
 

// ================================================ ====== ==== ==== == = 
// Destroys all structures allocated in \c ComputePreconditioner() if the preconditioner has been computed.
int ML_Epetra::EdgeMatrixFreePreconditioner::DestroyPreconditioner(){
  if (ml_comm_) { ML_Comm_Destroy(&ml_comm_); ml_comm_ = 0; }// will need this
  if (Prolongator_) {delete Prolongator_; Prolongator_=0;}
  if (InvDiagonal_) {delete InvDiagonal_; InvDiagonal_=0;}    
  if (CoarseMatrix) {delete CoarseMatrix; CoarseMatrix=0;}
  if (CoarseMat_ML) {ML_Operator_Destroy(&CoarseMat_ML);CoarseMat_ML=0;}
  if (CoarsePC) {delete CoarsePC; CoarsePC=0;}
  if (CoarseMap_) {delete CoarseMap_; CoarseMap_=0;}
  if (Smoother_) {delete Smoother_; Smoother_=0;}
  return 0;
}/*end DestroyPreconditioner*/



// ================================================ ====== ==== ==== == = 
//! Apply the preconditioner to an Epetra_MultiVector X, puts the result in Y
int ML_Epetra::EdgeMatrixFreePreconditioner::ApplyInverse(const Epetra_MultiVector& B, Epetra_MultiVector& X) const{
  static int iteration=0;  //HAQ
  double re1=0,rn1=0,re2=0,re3=0;


  
  /* Sanity Checks */
  int NumVectors=B.NumVectors();
  if (!B.Map().SameAs(*EdgeDomainMap_)) ML_CHK_ERR(-1);
  if (NumVectors != X.NumVectors()) ML_CHK_ERR(-1);

  Epetra_MultiVector r_edge(*EdgeDomainMap_,NumVectors,false);
  Epetra_MultiVector e_edge(*EdgeDomainMap_,NumVectors,false);
  Epetra_MultiVector e_node(*CoarseMap_,NumVectors,false);
  Epetra_MultiVector r_node(*CoarseMap_,NumVectors,false);

  for(int i=0;i<num_cycles;i++){    
    /* Pre-smoothing */
    MVOUT2(X,"xinit11",iteration);
   
    double re0=cms_compute_residual(Operator_,B,X);
    if(Smoother_) ML_CHK_ERR(Smoother_->ApplyInverse(B,X));
    if(verbose_) re1=cms_compute_residual(Operator_,B,X);
    
    MVOUT2(X,"sm11-1",iteration);

    if(MaxLevels > 1){

      /* Calculate Residual (r_e = b - (S+M+Addon) * x) */
      ML_CHK_ERR(Operator_->Apply(X,r_edge));
      ML_CHK_ERR(r_edge.Update(1.0,B,-1.0));
      
      MVOUT2(r_edge,"re11",iteration);
      
      /* Xfer to coarse grid (r_n = P' * r_e) */
      ML_CHK_ERR(Prolongator_->Multiply(true,r_edge,r_node));
      
      MVOUT2(r_node,"rn11",iteration);
      
      /* AMG on coarse grid  (e_n = (CoarseMatrix)^{-1} r_n) */
      ML_CHK_ERR(CoarsePC->ApplyInverse(r_node,e_node));
      
      MVOUT2(e_node,"en11",iteration);
      
      if(verbose_) rn1=cms_compute_residual(CoarseMatrix,r_node,e_node);      
      
      /* Xfer back to fine grid (e_e = P * e_n) */
      ML_CHK_ERR(Prolongator_->Multiply(false,e_node,e_edge));
      
      MVOUT2(e_edge,"ee11",iteration);
      
      /* Add in correction (x = x + e_e) */
      ML_CHK_ERR(X.Update(1.0,e_edge,1.0));
      
      MVOUT2(X,"xup11",iteration);
      
      if(verbose_) re2=cms_compute_residual(Operator_,B,X);
    }/*end if*/
    
    /* Post-Smoothing*/
    if(Smoother_) ML_CHK_ERR(Smoother_->ApplyInverse(B,X));

    if(verbose_) re3=cms_compute_residual(Operator_,B,X);
    
    if(verbose_ && !Comm_->MyPID())
      printf("11 Resid Reduction: %6.4e / %6.4e / %6.4e / %6.4e\n",re1/re0,re2/re1,re3/re2,re3/re0);

  }/*end for*/

  /* Debug - Dump Vectors */
  MVOUT2(B,"b11",iteration);
  MVOUT2(X,"x11",iteration);
  iteration++;//HAQ
  
  return 0;
}/*end ApplyInverse*/



#endif
