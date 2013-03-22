/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */
#include "ml_config.h"
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_EPETRAEXT)
#include "ml_FaceMatrixFreePreconditioner.h"
#include "ml_MultiLevelPreconditioner.h"
#include "ml_epetra.h"
#include "ml_epetra_utils.h"
#include "ml_mat_formats.h"
#include "ml_RefMaxwell_11_Operator.h"
#include "ml_RefMaxwell_Utils.h"
#include "ml_ifpack_epetra_wrap.h"

#define ABS(x)((x)>0?(x):-(x))

#define NO_OUTPUT

#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"
#include "EpetraExt_VectorOut.h"
// ================================================ ====== ==== ==== == =
inline void cross_product(const double *a,const double *b,double *c){
  c[0] = a[1]*b[2]-a[2]*b[1];
  c[1] = a[2]*b[0]-a[0]*b[2];
  c[2] = a[0]*b[1]-a[1]*b[0];
}


// ================================================ ====== ==== ==== == = 
ML_Epetra::FaceMatrixFreePreconditioner::FaceMatrixFreePreconditioner(Teuchos::RCP<const Epetra_Operator> Operator, 
								      Teuchos::RCP<const Epetra_Vector> Diagonal,
								      Teuchos::RCP<const Epetra_CrsMatrix> FaceNode_Matrix,
								      Teuchos::RCP<const Epetra_CrsMatrix> TMT_Matrix,
								      Teuchos::ArrayRCP<int> BCfaces,
								      const Teuchos::ParameterList &List,const bool ComputePrec):
  ML_Preconditioner(),Prolongator_(0),InvDiagonal_(0),CoarseMatrix(0),CoarsePC(0),
#ifdef HAVE_ML_IFPACK
Smoother_(0),
#endif
verbose_(false),very_verbose_(false),print_hierarchy(false)
{

  Operator_=Operator;
  Diagonal_=Diagonal;
  FaceNode_Matrix_=FaceNode_Matrix;
  TMT_Matrix_=TMT_Matrix;
  BCfaces_=BCfaces;
  
  /* Set the Epetra Goodies */
  Comm_ = &(Operator_->Comm());

  FaceDomainMap_ = &(Operator_->OperatorDomainMap());
  FaceRangeMap_ = &(Operator_->OperatorRangeMap());
  NodeDomainMap_ = &(TMT_Matrix_->OperatorDomainMap());
  NodeRangeMap_ = &(TMT_Matrix_->OperatorRangeMap());
  
  List_=List;
  Label_=new char[80];
  strcpy(Label_,"ML face matrix-free preconditioner");  

  /* Parameter List Options */
  int OutputLevel = List_.get("ML output", -47);
  if(OutputLevel == -47) OutputLevel = List_.get("output", 1);
  if(OutputLevel>=15) very_verbose_=verbose_=true;
  if(OutputLevel > 5) {very_verbose_=false;verbose_=true;}
  else very_verbose_=verbose_=false;  


  // Pull diagonal if needed
  const Epetra_CrsMatrix *Op11crs = dynamic_cast<const Epetra_CrsMatrix*>(&*Operator_);
  int SmootherSweeps = List_.get("smoother: sweeps", 0);
  if(SmootherSweeps){
    if(Diagonal==Teuchos::null && Op11crs){    
      if(verbose_ && !Comm_->MyPID()) printf("FMFP: Extracting matrix diagonal\n");
      InvDiagonal_ = new Epetra_Vector(Op11crs->RowMap());  
      Op11crs->ExtractDiagonalCopy(*InvDiagonal_);
    }	
    else {
      InvDiagonal_ = new Epetra_Vector(*Diagonal);  
    }
  }

  if(ComputePrec) ML_CHK_ERRV(ComputePreconditioner());
}/*end constructor*/

// ================================================ ====== ==== ==== == =   
ML_Epetra::FaceMatrixFreePreconditioner::~FaceMatrixFreePreconditioner(){
  DestroyPreconditioner();
}/*end destructor*/


// ================================================ ====== ==== ==== == = 
// Computes the preconditioner
int ML_Epetra::FaceMatrixFreePreconditioner::ComputePreconditioner(const bool CheckFiltering)
{
  Teuchos::ParameterList & ListCoarse=List_.sublist("face matrix free: coarse");
  
  /* ML Communicator */
  ML_Comm_Create(&ml_comm_);

  /* Parameter List Options */
  int SmootherSweeps = List_.get("smoother: sweeps", 0);
  MaxLevels = List_.get("max levels",10); 
  print_hierarchy= List_.get("print hierarchy",false);    
  num_cycles  = List_.get("cycle applications",1);

  /* Sanity Checking*/
  int OperatorDomainPoints =  OperatorDomainMap().NumGlobalPoints();
  int OperatorRangePoints =  OperatorRangeMap().NumGlobalPoints();
  if (OperatorDomainPoints != OperatorRangePoints)
    ML_CHK_ERR(-1); // only square matrices

 /* Output Header */
  if(verbose_ && !Comm_->MyPID()) {
    printf("------------------------------------------------------------------------------\n");
    printf("***\n");
    printf("*** ML_Epetra::FaceMatrixFreePreconditioner [%s]\n",Label());
    printf("***\n");
  }

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
  if(SmootherSweeps) 
    ML_CHK_ERR(SetupSmoother());


  if(MaxLevels > 0) {  
    /* Build the prolongator */
    ML_CHK_ERR(BuildProlongator());
    
    /* DEBUG: Output matrices */
    if(print_hierarchy) EpetraExt::RowMatrixToMatlabFile("prolongator.dat",*Prolongator_);
    
    /* Form the coarse matrix */
    ML_CHK_ERR(FormCoarseMatrix());

    /* DEBUG: Output matrices */
    if(print_hierarchy) EpetraExt::RowMatrixToMatlabFile("coarsemat.dat",*CoarseMatrix);
    
    /* Setup Preconditioner on Coarse Matrix */
    CoarsePC = new MultiLevelPreconditioner(*CoarseMatrix,ListCoarse);
    if(!CoarsePC) ML_CHK_ERR(-2);
  
    /* Clean Up */
    //delete nullspace;
  }/*end if*/
    
  return 0;
}/*end ComputePreconditioner*/


// ================================================ ====== ==== ==== == = 
// Setup the Smoother
int ML_Epetra::FaceMatrixFreePreconditioner::SetupSmoother()
{

#ifdef HAVE_ML_IFPACK
  Smoother_=ML_Epetra::ML_Gen_Smoother_Ifpack_Epetra(&*Operator_,InvDiagonal_,List_,"FMFP Smoother (level 0): ",verbose_);
#else
  if(!Comm_->MyPID())
    printf("ERROR: FMFP must be compiled with --enable-ml-ifpack for this mode to work\n");
#endif
  return 0;
}/*end SetupSmoother */

// ================================================ ====== ==== ==== == = 
//! Build pi operator described by Bochev, Siefert, Tuminaro, Xu and Zhu (2007).
int ML_Epetra::FaceMatrixFreePreconditioner::BuildNullspace(Epetra_MultiVector *& nullspace){
  int Nf=FaceRangeMap_->NumMyElements();

  /* Pull the coordinates from Teuchos */
  double * xcoord=List_.get("x-coordinates",(double*)0);
  double * ycoord=List_.get("y-coordinates",(double*)0);
  double * zcoord=List_.get("z-coordinates",(double*)0);
  dim=(xcoord!=0) + (ycoord!=0) + (zcoord!=0);    

  // Sanity check
  if(dim!=3){if(!Comm_->MyPID()) printf("ERROR: FaceMatrixFreePreconditioner only works in 3D"); ML_CHK_ERR(-1);}

  // Build the (unimported) coordinate multivector
  double **d_coords=new double* [dim];
  d_coords[0]=xcoord; d_coords[1]=ycoord;
  if(dim==3) d_coords[2]=zcoord;  
  Epetra_MultiVector n_coords_domain(View,*NodeDomainMap_,d_coords,dim);    
  Epetra_MultiVector *n_coords;

  // Import coordinate info
  if(FaceNode_Matrix_->Importer()){
    n_coords=new Epetra_MultiVector(FaceNode_Matrix_->ColMap(),dim);    
    n_coords->PutScalar(0.0);
    n_coords->Import(n_coords_domain,*FaceNode_Matrix_->Importer(),Add);
  }
  else n_coords=&n_coords_domain;

  // Sanity HAQ - Only works on Hexes
  if(FaceNode_Matrix_->GlobalMaxNumEntries()!=4)
    {if(!Comm_->MyPID()) printf("ERROR: FaceMatrixFreePreconditioner only works on Hexes"); ML_CHK_ERR(-2);}
    
  // Allocate vector
  nullspace=new Epetra_MultiVector(*FaceDomainMap_,3);

  // Fill matrix - NTS this will *NOT* do periodic BCs correctly.
  double *a=new double[dim];
  double *b=new double[dim];
  double *c=new double[dim];
  for(int i=0;i<Nf;i++){
    int Ni,*indices;
    double *values;
    FaceNode_Matrix_->ExtractMyRowView(i,Ni,values,indices);
    if(Ni != 4){printf("ERROR: Face %d has only %d nodes\n",i,Ni); ML_CHK_ERR(-1);}

    a[0] = (*n_coords)[0][indices[1]] - (*n_coords)[0][indices[0]];
    a[1] = (*n_coords)[1][indices[1]] - (*n_coords)[1][indices[0]];
    a[2] = (*n_coords)[2][indices[1]] - (*n_coords)[2][indices[0]];

    b[0] = (*n_coords)[0][indices[2]] - (*n_coords)[0][indices[0]];
    b[1] = (*n_coords)[1][indices[2]] - (*n_coords)[1][indices[0]];
    b[2] = (*n_coords)[2][indices[2]] - (*n_coords)[2][indices[0]];

    cross_product(a,b,c);

    // HAQ - Hardwiring for hexes
    // HAQ - Absolute value, presuming all hexes are actually pointed the same way.  This is a HAQ!!!
    (*nullspace)[0][i]=ABS(c[0])/6.0;
    (*nullspace)[1][i]=ABS(c[1])/6.0;
    (*nullspace)[2][i]=ABS(c[2])/6.0;
  }

  /* Cleanup */
  if(FaceNode_Matrix_->Importer()) delete n_coords;
  delete [] a; delete [] b; delete [] c;
  delete [] d_coords;
  return 0;
}

// ================================================ ====== ==== ==== == = 
//! Build the face-to-node prolongator described by Bochev, Siefert, Tuminaro, Xu and Zhu (2007).
int ML_Epetra::FaceMatrixFreePreconditioner::PBuildSparsity(ML_Operator *P, Epetra_CrsMatrix *&Psparse){

  /* Create wrapper to do abs(T) */
  // NTS: Assume D0 has already been reindexed by now.
  ML_Operator* AbsFN_ML = ML_Operator_Create(ml_comm_);
  ML_CHK_ERR(ML_Operator_WrapEpetraCrsMatrix(const_cast<Epetra_CrsMatrix*>(&*FaceNode_Matrix_),AbsFN_ML,verbose_));    
  ML_Operator_Set_Getrow(AbsFN_ML,AbsFN_ML->outvec_leng,CSR_getrow_ones);
  
  /* Form abs(T) * P_n */
  ML_Operator* AbsFNP = ML_Operator_Create(ml_comm_);   
  ML_2matmult(AbsFN_ML,P,AbsFNP, ML_CSR_MATRIX);
  
  /* Wrap P_n into Epetra-land */
  Epetra_CrsMatrix_Wrap_ML_Operator(AbsFNP,*Comm_,*FaceRangeMap_,&Psparse,Copy,0);
  
  /* Nuke the rows in Psparse */
  if(BCfaces_.size()>0) Apply_BCsToMatrixRows(BCfaces_.get(),BCfaces_.size(),*Psparse);

  // Cleanup
  ML_Operator_Destroy(&AbsFN_ML);
  ML_Operator_Destroy(&AbsFNP);
  
  return 0;
}

// ================================================ ====== ==== ==== == = 
//! Build the face-to-node prolongator described by Bochev, Siefert, Tuminaro, Xu and Zhu (2007).
int ML_Epetra::FaceMatrixFreePreconditioner::BuildProlongator()
{

  /* Wrap TMT_Matrix in a ML_Operator */
  ML_Operator* TMT_ML = ML_Operator_Create(ml_comm_);
  ML_Operator_WrapEpetraCrsMatrix(const_cast<Epetra_CrsMatrix*>(&*TMT_Matrix_),TMT_ML);

  /* Nodal Aggregation */
  ML_Aggregate_Struct *MLAggr=0;
  ML_Operator *P=0;
  int NumAggregates;
  int rv=ML_Epetra::RefMaxwell_Aggregate_Nodes(*TMT_Matrix_,List_,ml_comm_,std::string("FMFP (level 0) :"),MLAggr,P,NumAggregates);
  if(rv || !P) {if(!Comm_->MyPID()) printf("ERROR: Building nodal P\n");ML_CHK_ERR(-1);}

  /* Build 1-unknown sparsity of prolongator */
  Epetra_CrsMatrix *Psparse=0;
  PBuildSparsity(P,Psparse);
  if(!Psparse) {if(!Comm_->MyPID()) printf("ERROR: Building Psparse\n");ML_CHK_ERR(-2);}

  /* Build the "nullspace" */
  Epetra_MultiVector *nullspace;
  BuildNullspace(nullspace);
  if(!nullspace) {if(!Comm_->MyPID()) printf("ERROR: Building Nullspace\n");ML_CHK_ERR(-3);}

    
  /* Build the DomainMap of the new operator*/
  const Epetra_Map & FineColMap = Psparse->ColMap();
  CoarseMap_=new Epetra_Map(-1,NumAggregates*dim,0,*Comm_);

  /* Allocate the Prolongator_ */ 
  int max_nz_per_row=Psparse->MaxNumEntries();
  Prolongator_=new Epetra_CrsMatrix(Copy,*FaceRangeMap_,0);
  int ne1, *idx1, *idx2;
  idx2=new int [dim*max_nz_per_row];
  double *vals1, *vals2;
  vals2=new double[dim*max_nz_per_row];
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
	else vals2[j*dim+k]=(*nullspace)[k][i] / nonzeros;  
      }/*end for*/
    }/*end for*/
    Prolongator_->InsertGlobalValues(FaceRangeMap_->GID(i),dim*ne1,vals2,idx2);
  }/*end for*/
  
  
  /* FillComplete / OptimizeStorage for Prolongator*/
  Prolongator_->FillComplete(*CoarseMap_,*FaceRangeMap_);
  Prolongator_->OptimizeStorage();

#ifndef NO_OUTPUT
  /* DEBUG: Dump aggregates */ 
  Epetra_IntVector AGG(View,*NodeDomainMap_,MLAggr->aggr_info[0]);
  IVOUT(AGG,"agg.dat");  
  EpetraExt::RowMatrixToMatlabFile("fmfp_psparse.dat",*Psparse);
  EpetraExt::RowMatrixToMatlabFile("fmfp_prolongator.dat",*Prolongator_);
  EpetraExt::VectorToMatrixMarketFile("fmfp_null0.dat",*(*nullspace)(0));
  EpetraExt::VectorToMatrixMarketFile("fmfp_null1.dat",*(*nullspace)(1));
  EpetraExt::VectorToMatrixMarketFile("fmfp_null2.dat",*(*nullspace)(2));

#endif

  /* EXPERIMENTAL: Normalize Prolongator Columns */
  bool normalize_prolongator=List_.get("face matrix free: normalize prolongator",false);
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
  ML_qr_fix_Destroy();
  ML_Aggregate_Destroy(&MLAggr);
  ML_Operator_Destroy(&TMT_ML);
  ML_Operator_Destroy(&P);

  delete nullspace;
  delete Psparse;
  delete [] idx2;
  delete [] vals2;
  return 0;
}/*end BuildProlongator_*/
  



// ================================================ ====== ==== ==== == = 
// Forms the coarse matrix, given the prolongator
int  ML_Epetra::FaceMatrixFreePreconditioner::FormCoarseMatrix()
{
  CoarseMat_ML = ML_Operator_Create(ml_comm_);
  CoarseMat_ML->data_destroy=free;
  ML_Operator *Temp_ML=0;
  ML_Operator *R= ML_Operator_Create(ml_comm_);
  ML_Operator *P= ML_Operator_Create(ml_comm_);

  /* Build ML_Operator version of Prolongator_, Restriction Operator */
  ML_CHK_ERR(ML_Operator_WrapEpetraCrsMatrix(Prolongator_,P,verbose_));
  P->num_rigid=P->num_PDEs=dim;
  
  //NTS: ML_CHK_ERR won't work on this: it returns 1
  ML_Operator_Transpose_byrow(P, R);
  
  /* OPTION: Disable the addon */
  const Epetra_CrsMatrix *Op11crs = dynamic_cast<const Epetra_CrsMatrix*>(&*Operator_);
  const Epetra_Operator_With_MatMat *Op11mm = dynamic_cast<const Epetra_Operator_With_MatMat*>(&*Operator_);

  /* Do the A*P  with or without addon*/
  if(Op11crs){
    if(verbose_ && !Comm_->MyPID()) printf("FMFP: Running *without* addon\n");
    ML_Operator *SM_ML = ML_Operator_Create(ml_comm_);
    Temp_ML = ML_Operator_Create(ml_comm_);
    ML_Operator_WrapEpetraCrsMatrix((Epetra_CrsMatrix*)Op11crs,SM_ML,verbose_);
    ML_2matmult(SM_ML,P,Temp_ML,ML_CSR_MATRIX);
    ML_Operator_Destroy(&SM_ML);
  }
  else if(Op11mm){
    if(verbose_ && !Comm_->MyPID()) printf("FMFP: Running with addon\n");
    ML_CHK_ERR(Op11mm->MatrixMatrix_Multiply(*Prolongator_,ml_comm_,&Temp_ML));  
  }
  else{
    if(!Comm_->MyPID()) printf("ERROR: FMFP Illegal Operator\n");
    ML_CHK_ERR(-1);
  }

  /* Do R * AP */
  R->num_rigid=R->num_PDEs=dim;
  ML_2matmult_block(R, Temp_ML,CoarseMat_ML,ML_CSR_MATRIX);
  
  /* Wrap to Epetra-land */
  int nnz=100;
  double time;
  ML_Operator2EpetraCrsMatrix(CoarseMat_ML,CoarseMatrix,nnz,true,time,0,verbose_);
  // NTS: This is a hack to get around the sticking ones on the diagonal issue;
    
  /* Cleanup */
  ML_Operator_Destroy(&P);
  ML_Operator_Destroy(&R);
  ML_Operator_Destroy(&Temp_ML);
  ML_Operator_Destroy(&CoarseMat_ML);CoarseMat_ML=0;//HAX  
  return 0;
}/*end FormCoarseMatrix*/

// ================================================ ====== ==== ==== == = 
// Print the individual operators in the multigrid hierarchy.
void ML_Epetra::FaceMatrixFreePreconditioner::Print(int whichHierarchy)
{
  if(CoarsePC) CoarsePC->Print();
}/*end Print*/
 

// ================================================ ====== ==== ==== == = 
// Destroys all structures allocated in \c ComputePreconditioner() if the preconditioner has been computed.
int ML_Epetra::FaceMatrixFreePreconditioner::DestroyPreconditioner(){
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
int ML_Epetra::FaceMatrixFreePreconditioner::ApplyInverse(const Epetra_MultiVector& B_, Epetra_MultiVector& X) const{
  const Epetra_MultiVector *B;
  Epetra_MultiVector *Bcopy=0;

  /* Sanity Checks */
  int NumVectors=B_.NumVectors();
  if (!B_.Map().SameAs(*FaceDomainMap_)) ML_CHK_ERR(-1);
  if (NumVectors != X.NumVectors()) ML_CHK_ERR(-1);

  Epetra_MultiVector r_edge(*FaceDomainMap_,NumVectors,false);
  Epetra_MultiVector e_edge(*FaceDomainMap_,NumVectors,false);
  Epetra_MultiVector e_node(*CoarseMap_,NumVectors,false);
  Epetra_MultiVector r_node(*CoarseMap_,NumVectors,false);

  /* Deal with the B==X case */
  if (B_.Pointers()[0] == X.Pointers()[0]){
    Bcopy=new Epetra_MultiVector(B_);
    B=Bcopy;
    X.PutScalar(0.0);  
  }    
  else B=&B_;


  for(int i=0;i<num_cycles;i++){    
    /* Pre-smoothing */
#ifdef HAVE_ML_IFPACK
    if(Smoother_) ML_CHK_ERR(Smoother_->ApplyInverse(*B,X));
#endif
    
    if(MaxLevels > 0){
      if(i != 0
#ifdef HAVE_ML_IFPACK
         || Smoother_
#endif
         ){
        /* Calculate Residual (r_e = b - (S+M+Addon) * x) */
        ML_CHK_ERR(Operator_->Apply(X,r_edge));
        ML_CHK_ERR(r_edge.Update(1.0,*B,-1.0));
        
        /* Xfer to coarse grid (r_n = P' * r_e) */
        ML_CHK_ERR(Prolongator_->Multiply(true,r_edge,r_node));
      }
      else{
        /* Xfer to coarse grid (r_n = P' * r_e) */
        ML_CHK_ERR(Prolongator_->Multiply(true,*B,r_node));
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
    if(Smoother_) ML_CHK_ERR(Smoother_->ApplyInverse(*B,X));
#endif
    
  }/*end for*/

  /* Cleanup */
  if(Bcopy) delete Bcopy;
  
  return 0;
}/*end ApplyInverse*/



#endif
