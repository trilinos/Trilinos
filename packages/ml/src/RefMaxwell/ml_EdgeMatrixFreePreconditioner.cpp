/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */
#include "ml_config.h"
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_EPETRAEXT)
#include "ml_EdgeMatrixFreePreconditioner.h"
#include "ml_MultiLevelPreconditioner.h"
#include "ml_epetra.h"
#include "ml_epetra_utils.h"
#include "ml_mat_formats.h"
#include "ml_RefMaxwell_11_Operator.h"
using namespace std;

//mucho hax
void cms_residual_check(const char * tag, const Epetra_Operator * op,const Epetra_MultiVector& rhs, const Epetra_MultiVector& lhs);
double cms_compute_residual(const Epetra_Operator * op,const Epetra_MultiVector& rhs, const Epetra_MultiVector& lhs);

#define ABS(x)((x)>0?(x):-(x))

#define NO_OUTPUT
extern void Epetra_CrsMatrix_Print(const Epetra_CrsMatrix& A, const char* of);
extern void ML_Matrix_Print(ML_Operator *ML,const Epetra_Comm &Comm,const Epetra_Map &Map, const char *fname);
extern void MVOUT (const Epetra_MultiVector & A, const char * of);
extern void IVOUT(const Epetra_IntVector & A, const char * of);
extern void MVOUT2(const Epetra_MultiVector & A,const char* pref,int idx);

//#define ENABLE_FAST_PTAP // This has a bug.  Leave it off for now -CMS
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_BlockMapOut.h"
#include "Teuchos_XMLObject.hpp"
#include "Teuchos_XMLParameterListWriter.hpp"
/*#define KDD_DEBUG*/

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
  ML_Preconditioner(),Operator_(&Operator),D0_Matrix_(&D0_Matrix),D0_Clean_Matrix_(&D0_Clean_Matrix),TMT_Matrix_(&TMT_Matrix),BCedges_(BCedges),numBCedges_(numBCedges),Prolongator_(0),InvDiagonal_(0),CoarseMatrix(0),CoarsePC(0),
#ifdef HAVE_ML_IFPACK
Smoother_(0),
#endif
verbose_(false),very_verbose_(false),print_hierarchy(false)
{
  /* Set the Epetra Goodies */
  Comm_ = &(Operator_->Comm());

  EdgeDomainMap_ = &(Operator_->OperatorDomainMap());
  EdgeRangeMap_ = &(Operator_->OperatorRangeMap());
  NodeDomainMap_ = &(TMT_Matrix_->OperatorDomainMap());
  NodeRangeMap_ = &(TMT_Matrix_->OperatorRangeMap());
  
  List_=List;
  Label_=new char[80];
  strcpy(Label_,"ML edge matrix-free preconditioner");
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
  if(OutputLevel == -47) OutputLevel = List_.get("output", 1);
  if(OutputLevel>=15) very_verbose_=verbose_=true;
  if(OutputLevel > 5) {very_verbose_=false;verbose_=true;}
  else very_verbose_=verbose_=false;  
  
  int SmootherSweeps = List_.get("smoother: sweeps (level 0)", 0);
  MaxLevels = List_.get("max levels",10); 
  print_hierarchy= List_.get("print hierarchy",false);  
  
  num_cycles  = List_.get("cycle applications",1);
  //  ML_Set_PrintLevel(OutputLevel);

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


  if(MaxLevels > 0) {  
    /* Build the Nullspace */
    if(verbose_ && !Comm_->MyPID()) printf("EMFP: Building Nullspace\n");
    Epetra_MultiVector *nullspace=BuildNullspace();
    if(!nullspace) ML_CHK_ERR(-1);
    if(print_hierarchy) MVOUT(*nullspace,"nullspace.dat");
    
    /* Build the prolongator */
    if(verbose_ && !Comm_->MyPID()) printf("EMFP: Building Prolongator\n");
    ML_CHK_ERR(BuildProlongator(*nullspace));
    
    /* DEBUG: Output matrices */
    if(print_hierarchy) Epetra_CrsMatrix_Print(*Prolongator_,"prolongator.dat");
    
    /* Form the coarse matrix */
    if(verbose_ && !Comm_->MyPID()) printf("EMFP: Building Coarse Matrix\n");
    ML_CHK_ERR(FormCoarseMatrix());

    /* DEBUG: Output matrices */
    if(print_hierarchy) Epetra_CrsMatrix_Print(*CoarseMatrix,"coarsemat.dat");
    
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
#ifdef HAVE_ML_IFPACK
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


  if(print_hierarchy) MVOUT(*InvDiagonal_,"inv_diagonal.dat");
  
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
#else
  if(!Comm_->MyPID())
    printf("ERROR: RefMaxwell must be compiled with --enable-ml-ifpack for this mode to work\n");
#endif
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

  /* Normalize */
  double d1 = sqrt(ML_gdot(NodeDomainMap_->NumMyElements(), xcoord, xcoord, ml_comm_));
  for (int i = 0; i < NodeDomainMap_->NumMyElements(); i++) xcoord[i] /= d1;
  d1 = sqrt(ML_gdot(NodeDomainMap_->NumMyElements(), ycoord, ycoord, ml_comm_));
  for (int i = 0; i < NodeDomainMap_->NumMyElements(); i++) ycoord[i] /= d1;
  if (dim==3) {
    d1 = sqrt(ML_gdot(NodeDomainMap_->NumMyElements(), zcoord, zcoord, ml_comm_));
    for (int i = 0; i < NodeDomainMap_->NumMyElements(); i++) zcoord[i] /= d1;
  }

  /* Build the MultiVector */
  double ** d_coords=new double* [dim];
  d_coords[0]=xcoord; d_coords[1]=ycoord;
  if(dim==3) d_coords[2]=zcoord;
  Epetra_MultiVector e_coords(View,*NodeDomainMap_,d_coords,dim);

  if(print_hierarchy) MVOUT(e_coords,"coords.dat");
  
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
    if(!Comm_->MyPID()) printf("RefMaxwell: Unsupported (1,1) block aggregation type(%s), resetting to uncoupled-mis\n",CoarsenType.c_str());
    ML_Aggregate_Set_CoarsenScheme_UncoupledMIS(MLAggr);
  }

  /* Aggregate Nodes */
  int NumAggregates = ML_Aggregate_Coarsen(MLAggr, TMT_ML, &P, ml_comm_);
  if (NumAggregates == 0){
    cerr << "Found 0 aggregates, perhaps the problem is too small." << endl;
    ML_CHK_ERR(-2);
  }/*end if*/
  else if(very_verbose_) printf("[%d] EMFP: %d aggregates created invec_leng=%d\n",Comm_->MyPID(),NumAggregates,P->invec_leng);

  
  if(very_verbose_) printf("[%d] Num Aggregates = %d\n",Comm_->MyPID(),NumAggregates);
  if(P==0) {fprintf(stderr,"ERROR: No tentative prolongator found\n");ML_CHK_ERR(-5);}
  
#ifndef NO_OUTPUT
  /* DEBUG: Dump aggregates, prolongator */ 
  Epetra_IntVector AGG(View,*NodeDomainMap_,MLAggr->aggr_info[0]);
  IVOUT(AGG,"agg.dat");  
  ML_Operator_Print(P,"p.dat");
#endif
  
  /* Create wrapper to do abs(T) */
  // NTS: Assume D0 has already been reindexed by now.
  if(verbose_ && !Comm_->MyPID()) printf("EMFP: abs(T) prewrap\n");
  ML_Operator* AbsD0_ML = ML_Operator_Create(ml_comm_);
  ML_CHK_ERR(ML_Operator_WrapEpetraCrsMatrix((Epetra_CrsMatrix*)D0_Matrix_,AbsD0_ML,verbose_));    
  ML_Operator_Set_Getrow(AbsD0_ML,AbsD0_ML->outvec_leng,CSR_getrow_ones);

  
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
    for(int j=0;j<ne1;j++) nonzeros+=ABS(vals1[j])>0;

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

  /* EXPERIMENTAL: Normalize Prolongator Columns */
  bool normalize_prolongator=List_.get("refmaxwell: normalize prolongator",false);
  if(normalize_prolongator){
    Epetra_Vector n_vector(*CoarseMap_,false);
    Prolongator_->InvColSums(n_vector);
    Prolongator_->RightScale(n_vector);
  }/*end if*/
  
  /* Post-wrapping to convert to ML indexing */
#ifdef HAVE_ML_EPETRAEXT
  Prolongator_ = dynamic_cast<Epetra_CrsMatrix*>(ModifyEpetraMatrixColMap(*Prolongator_,ProlongatorColMapTrans_,"Prolongator",(verbose_&&!Comm_->MyPID())));
#endif
  
  /* Cleanup */
  if(verbose_ && !Comm_->MyPID()) printf("EMFP: BuildProlongator Cleanup\n");  
  ML_qr_fix_Destroy();
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
  CoarseMat_ML = ML_Operator_Create(ml_comm_);
  CoarseMat_ML->data_destroy=free;
  ML_Operator *Temp_ML=0;
#ifndef ENABLE_FAST_PTAP  
  ML_Operator *R= ML_Operator_Create(ml_comm_);
  ML_Operator *P= ML_Operator_Create(ml_comm_);

  /* Build ML_Operator version of Prolongator_, Restriction Operator */
  if(verbose_ && !Comm_->MyPID()) printf("EMFP: Prolongator Prewrap\n");
  ML_CHK_ERR(ML_Operator_WrapEpetraCrsMatrix(Prolongator_,P,verbose_));
  P->num_rigid=P->num_PDEs=dim;
  
  if(verbose_ && !Comm_->MyPID()) printf("EMFP: Prolongator Transpose\n");
  //NTS: ML_CHK_ERR won't work on this: it returns 1
  ML_Operator_Transpose_byrow(P, R);
#else
  if(verbose_ && !Comm_->MyPID()) printf("EMFP: Running FAST_PTAP\n");
#endif

  /* EXPERIMENTAL: Disable the addon */
  bool disable_addon=List_.get("refmaxwell: disable addon",false);
  const ML_RefMaxwell_11_Operator *Op11 = dynamic_cast<const ML_RefMaxwell_11_Operator*>(Operator_);
#ifndef ENABLE_FAST_PTAP
  if(disable_addon && Op11){
    if(verbose_ && !Comm_->MyPID()) printf("EMFP: AP (*without* addon)\n");
    ML_Operator *SM_ML = ML_Operator_Create(ml_comm_);
    Temp_ML = ML_Operator_Create(ml_comm_);
    ML_Operator_WrapEpetraCrsMatrix((Epetra_CrsMatrix*)&(((ML_RefMaxwell_11_Operator *)Op11)->SM_Matrix()),SM_ML,verbose_);
    ML_2matmult(SM_ML,P,Temp_ML,ML_CSR_MATRIX);
    ML_Operator_Destroy(&SM_ML);
  }
  else{
#endif
#ifdef ENABLE_FAST_PTAP
    Op11->PtAP(*Prolongator_,ml_comm_,&CoarseMat_ML);
#else
    /* Do the A*P */
    if(verbose_ && !Comm_->MyPID()) printf("EMFP: AP\n");
    ML_CHK_ERR(Operator_->MatrixMatrix_Multiply(*Prolongator_,ml_comm_,&Temp_ML));  
  }
#endif

#ifndef ENABLE_FAST_PTAP
  /* DEBUG: output*/
#ifndef NO_OUTPUT
  ML_Matrix_Print(Temp_ML,*Comm_,*EdgeRangeMap_,"coarse_temp.dat");
#endif


  /* Do R * AP */
  if(verbose_ && !Comm_->MyPID()) printf("EMFP: RAP\n");
  R->num_rigid=R->num_PDEs=dim;
  ML_2matmult_block(R, Temp_ML,CoarseMat_ML,ML_CSR_MATRIX);
#endif

  
  /* Wrap to Epetra-land */
  //  Epetra_CrsMatrix_Wrap_ML_Operator(CoarseMat_ML,*Comm_,*CoarseMap_,&CoarseMatrix); 
  int nnz=100;
  double time;
  ML_Operator2EpetraCrsMatrix(CoarseMat_ML,CoarseMatrix,nnz,true,time,0,verbose_);
  // NTS: This is a hack to get around the sticking ones on the diagonal issue;
  
#ifndef NO_OUTPUT
  ML_Matrix_Print(CoarseMat_ML,*Comm_,*EdgeRangeMap_,"coarsemat.dat");  
  Epetra_CrsMatrix_Print(*CoarseMatrix,"coarsemat0.dat");
#endif
  
#ifdef KDD_DEBUG
  EpetraExt::BlockMapToMatrixMarketFile("coarsemat0_map.dat",CoarseMatrix->RowMap());
  EpetraExt::RowMatrixToMatrixMarketFile("coarsemat0.dat",*CoarseMatrix);

  ofstream opl("coarsemat0_params.dat");
  Teuchos::ParameterList dummy, ListCoarse;
  ListCoarse=List_.get("edge matrix free: coarse",dummy);
  Teuchos::XMLParameterListWriter XMLR;
  opl << XMLR.toXML(ListCoarse);
#endif



  /* Cleanup */
#ifndef ENABLE_FAST_PTAP
  ML_Operator_Destroy(&P);
  ML_Operator_Destroy(&R);
#else
  if(Temp_ML)
#endif
    ML_Operator_Destroy(&Temp_ML);

  ML_Operator_Destroy(&CoarseMat_ML);CoarseMat_ML=0;//HAX  
  return 0;
}/*end FormCoarseMatrix*/

// ================================================ ====== ==== ==== == = 
// Print the individual operators in the multigrid hierarchy.
void ML_Epetra::EdgeMatrixFreePreconditioner::Print(const char *whichHierarchy)
{
  ofstream ofs("Pmat.edge.m");
  if(Prolongator_) Prolongator_->Print(ofs);
  if(CoarsePC) CoarsePC->Print();
}/*end Print*/
 

// ================================================ ====== ==== ==== == = 
// Destroys all structures allocated in \c ComputePreconditioner() if the preconditioner has been computed.
int ML_Epetra::EdgeMatrixFreePreconditioner::DestroyPreconditioner(){
  if (ml_comm_) { ML_Comm_Destroy(&ml_comm_); ml_comm_ = 0; }// will need this
  if (Prolongator_) {delete Prolongator_; Prolongator_=0;}
  if (InvDiagonal_) {delete InvDiagonal_; InvDiagonal_=0;}    
  if (CoarsePC) {delete CoarsePC; CoarsePC=0;}
  if (CoarseMatrix) {delete CoarseMatrix; CoarseMatrix=0;}
  if (CoarseMat_ML) {ML_Operator_Destroy(&CoarseMat_ML);CoarseMat_ML=0;}
  if (CoarseMap_) {delete CoarseMap_; CoarseMap_=0;}
#ifdef HAVE_ML_IFPACK  
  if (Smoother_){delete Smoother_; Smoother_=0;}
#endif
  return 0;
}/*end DestroyPreconditioner*/

// ================================================ ====== ==== ==== == = 
//! Apply the preconditioner to an Epetra_MultiVector X, puts the result in Y
int ML_Epetra::EdgeMatrixFreePreconditioner::ApplyInverse(const Epetra_MultiVector& B, Epetra_MultiVector& X) const{
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
#ifdef HAVE_ML_IFPACK
    if(Smoother_) ML_CHK_ERR(Smoother_->ApplyInverse(B,X));
#endif
    
    if(MaxLevels > 0){
      if(i != 0
#ifdef HAVE_ML_IFPACK
         || Smoother_
#endif
         ){
        /* Calculate Residual (r_e = b - (S+M+Addon) * x) */
        ML_CHK_ERR(Operator_->Apply(X,r_edge));
        ML_CHK_ERR(r_edge.Update(1.0,B,-1.0));
        
        /* Xfer to coarse grid (r_n = P' * r_e) */
        ML_CHK_ERR(Prolongator_->Multiply(true,r_edge,r_node));
      }
      else{
        /* Xfer to coarse grid (r_n = P' * r_e) */
        ML_CHK_ERR(Prolongator_->Multiply(true,B,r_node));
      }
        
      /* AMG on coarse grid  (e_n = (CoarseMatrix)^{-1} r_n) */
      ML_CHK_ERR(CoarsePC->ApplyInverse(r_node,e_node));
      
      /* Xfer back to fine grid (e_e = P * e_n) */
      ML_CHK_ERR(Prolongator_->Multiply(false,e_node,e_edge));
      
      /* Add in correction (x = x + e_e)        */
      ML_CHK_ERR(X.Update(1.0,e_edge,1.0));
    }/*end if*/
    
    /* Post-Smoothing */
#ifdef HAVE_ML_IFPACK
    if(Smoother_) ML_CHK_ERR(Smoother_->ApplyInverse(B,X));
#endif
    
  }/*end for*/
  return 0;
}/*end ApplyInverse*/



#endif
