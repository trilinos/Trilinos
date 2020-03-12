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
#include "ml_RefMaxwell_Utils.h"
#include "ml_ifpack_epetra_wrap.h"
#ifdef ML_MPI
#include "Epetra_MpiComm.h"
#endif

#define ABS(x)((x)>0?(x):-(x))

//#define ENABLE_FAST_PTAP // This has a bug.  Leave it off for now -CMS
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"


// ================================================ ====== ==== ==== == =
ML_Epetra::EdgeMatrixFreePreconditioner::EdgeMatrixFreePreconditioner(Teuchos::RCP<const Epetra_Operator> Operator,
								      Teuchos::RCP<const Epetra_Vector> Diagonal,
								      Teuchos::RCP<const Epetra_CrsMatrix> D0_Matrix,
								      Teuchos::RCP<const Epetra_CrsMatrix> D0_Clean_Matrix,
								      Teuchos::RCP<const Epetra_CrsMatrix> TMT_Matrix,
								      Teuchos::ArrayRCP<int> BCedges,
								      const Teuchos::ParameterList &List,const bool ComputePrec):
  ML_Preconditioner(),
  Prolongator_(0),InvDiagonal_(0),CoarseMatrix(0),CoarsePC(0),CoarseNullspace_(0),
#ifdef HAVE_ML_IFPACK
Smoother_(0),
#endif
CoarseMap_(0),
verbose_(false),very_verbose_(false),print_hierarchy(false)
{
  Operator_=Operator;
  D0_Matrix_=D0_Matrix;
  D0_Clean_Matrix_=D0_Clean_Matrix;
  TMT_Matrix_=TMT_Matrix;
  BCedges_=BCedges;

  /* Set the Epetra Goodies */
  Comm_ = &(Operator_->Comm());

  EdgeDomainMap_ = &(Operator_->OperatorDomainMap());
  EdgeRangeMap_ = &(Operator_->OperatorRangeMap());
  NodeDomainMap_ = &(TMT_Matrix_->OperatorDomainMap());
  NodeRangeMap_ = &(TMT_Matrix_->OperatorRangeMap());

  List_=List;
  Label_=new char[80];
  strcpy(Label_,"ML edge matrix-free preconditioner");


  // Pull diagonal if needed
  const Epetra_CrsMatrix *Op11crs = dynamic_cast<const Epetra_CrsMatrix*>(&*Operator_);
  if(Diagonal==Teuchos::null && Op11crs){
    InvDiagonal_ = new Epetra_Vector(Op11crs->RowMap());
    Op11crs->ExtractDiagonalCopy(*InvDiagonal_);
  }
  else
    InvDiagonal_ = new Epetra_Vector(*Diagonal);

  if(ComputePrec) ML_CHK_ERRV(ComputePreconditioner());
}/*end constructor*/

// ================================================ ====== ==== ==== == =
ML_Epetra::EdgeMatrixFreePreconditioner::~EdgeMatrixFreePreconditioner(){
  DestroyPreconditioner();
}/*end destructor*/


// ================================================ ====== ==== ==== == =
// Computes the preconditioner
int ML_Epetra::EdgeMatrixFreePreconditioner::ComputePreconditioner(const bool /* CheckFiltering */)
{
  Teuchos::ParameterList dummy, ListCoarse;

  /* ML Communicator */
  ML_Comm_Create(&ml_comm_);
#ifdef ML_MPI
  // Use the same communicator as A if we're using MPI.
  const Epetra_MpiComm * Mcomm=dynamic_cast<const Epetra_MpiComm*>(Comm_);
  if(Mcomm) ML_Comm_Set_UsrComm(ml_comm_,Mcomm->GetMpiComm());
#endif

  /* Parameter List Options */
  int OutputLevel = List_.get("ML output", -47);
  if(OutputLevel == -47) OutputLevel = List_.get("output", 1);
  if(OutputLevel>=15) very_verbose_=verbose_=true;
  if(OutputLevel > 5) {very_verbose_=false;verbose_=true;}
  else very_verbose_=verbose_=false;

  int SmootherSweeps = List_.get("smoother: sweeps", 0);
  MaxLevels = List_.get("max levels",10);
  print_hierarchy= List_.get("print hierarchy",false);
  num_cycles  = List_.get("cycle applications",1);



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
  if(SmootherSweeps)
    ML_CHK_ERR(SetupSmoother());


  if(MaxLevels > 0) {
    /* Build the Nullspace */
    Epetra_MultiVector *nullspace=ML_Epetra::Build_Edge_Nullspace(*D0_Clean_Matrix_,BCedges_,List_,verbose_);
    dim = nullspace->NumVectors();

    if(!nullspace) ML_CHK_ERR(-1);
    if(print_hierarchy) EpetraExt::MultiVectorToMatrixMarketFile("nullspace.dat",*nullspace,0,0,false);

    /* Build the prolongator */
    ML_CHK_ERR(BuildProlongator(*nullspace));
     
    /* DEBUG: Output matrices */
    if(print_hierarchy)
      EpetraExt::RowMatrixToMatlabFile("prolongator.dat",*Prolongator_);

    /* Form the coarse matrix */
    ML_CHK_ERR(FormCoarseMatrix());

    /* DEBUG: Output matrices */
    if(print_hierarchy) EpetraExt::RowMatrixToMatlabFile("coarsemat.dat",*CoarseMatrix);

    /* Setup Preconditioner on Coarse Matrix */
    ListCoarse=List_.get("edge matrix free: coarse",dummy);
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
  Smoother_=ML_Epetra::ML_Gen_Smoother_Ifpack_Epetra(&*Operator_,InvDiagonal_,List_,"EMFP Smoother (level 0): ",verbose_);
#else
  if(!Comm_->MyPID())
    printf("ERROR: EMFP must be compiled with --enable-ml-ifpack for this mode to work\n");
#endif

  return 0;
}/*end SetupSmoother */


// ================================================ ====== ==== ==== == =
//! Build the edge-to-vector-node prolongator described in Bochev, Hu, Siefert and Tuminaro (2006).
int ML_Epetra::EdgeMatrixFreePreconditioner::BuildProlongator(const Epetra_MultiVector & nullspace)
{
  /* Pull the (nodal) coordinates from Teuchos */
  double * xcoord=List_.get("x-coordinates",(double*)0);
  double * ycoord=List_.get("y-coordinates",(double*)0);
  double * zcoord=List_.get("z-coordinates",(double*)0);
  double * mcoord=List_.get("material coordinates",(double*)0);
  bool build_coarse_coords=true;
  if(!mcoord && dim!=(xcoord!=0) + (ycoord!=0) + (zcoord!=0) ) build_coarse_coords=false;
  
  /* Do the aggregation */
  ML_Aggregate_Struct * MLAggr=0;
  ML_Operator *P;
  int NumAggregates;
  int rv=ML_Epetra::RefMaxwell_Aggregate_Nodes(*TMT_Matrix_,List_,ml_comm_,std::string("EMFP (level 0) :"),
					       MLAggr,P,NumAggregates,CoarseCoord_);
  if(rv!=0) ML_CHK_ERR(-2);

  print_hierarchy= List_.get("print hierarchy",false);
  if (print_hierarchy) {
    /* Wrap P_n into Epetra-land */
    Epetra_CrsMatrix *P_epetra;
    Epetra_CrsMatrix_Wrap_ML_Operator(P,*Comm_,*NodeRangeMap_,&P_epetra,Copy,0);
    EpetraExt::RowMatrixToMatlabFile("P.dat",*P_epetra);
    EpetraExt::RowMatrixToMatrixMarketFile("P2.dat",*P_epetra, "P", "P", true);
  }

#if 0    
  /* Get aggregate information for nullspace */
  // NOTE: Should list the aggregate id for each row
  int * aggr_info;
  std::vector<int> dofs_per_agg(NumAggregates,0);
  if(MLAggr) {
    ML_Aggregate_Get_AggrMap(MLAggr,0,&aggr_info);
    for(int i=0; i<TMT_Matrix_->NumMyRows(); i++)
      dofs_per_agg[ aggr_info[i] ] ++;


    printf("DEBUG CMS: dofs_per_agg: ");
    for(int i=0; i<NumAggregates; i++)
      printf("%d ",dofs_per_agg[i]);
    printf("\n");
  }
#endif

  /* Create wrapper to do abs(T) */
  // NTS: Assume D0 has already been reindexed by now.
  ML_Operator* AbsD0_ML = ML_Operator_Create(ml_comm_);
  ML_CHK_ERR(ML_Operator_WrapEpetraCrsMatrix(const_cast<Epetra_CrsMatrix*>(&*D0_Matrix_),AbsD0_ML,verbose_));
  ML_Operator_Set_Getrow(AbsD0_ML,AbsD0_ML->outvec_leng,CSR_getrow_ones);

  /* Form abs(T) * P_n */
  ML_Operator* AbsD0P = ML_Operator_Create(ml_comm_);
  ML_2matmult(AbsD0_ML,P,AbsD0P, ML_CSR_MATRIX);

  /* Wrap P_n into Epetra-land */
  Epetra_CrsMatrix *Psparse;
  Epetra_CrsMatrix_Wrap_ML_Operator(AbsD0P,*Comm_,*EdgeRangeMap_,&Psparse,Copy,0);

  if (print_hierarchy) EpetraExt::RowMatrixToMatlabFile("P_intermediate.dat",*Psparse);

  /* Nuke the rows in Psparse */
  if(BCedges_.size()>0) Apply_BCsToMatrixRows(BCedges_.get(),BCedges_.size(),*Psparse);

  if (print_hierarchy) EpetraExt::RowMatrixToMatlabFile("P_intermediate2.dat",*Psparse);

  /* Build the DomainMap of the new operator*/
  const Epetra_Map & FineColMap = Psparse->ColMap();
  CoarseMap_=new Epetra_Map(-1,NumAggregates*dim,0,*Comm_);

  /* Allocate the Prolongator_ */
  Prolongator_=new Epetra_CrsMatrix(Copy,*EdgeRangeMap_,0);
  int ne1, *idx1, *idx2;
  idx2=new int [dim*AbsD0P->max_nz_per_row];
  double *vals1, *vals2;
  vals2=new double[dim*AbsD0P->max_nz_per_row];
  int nonzeros;

  /* EXPERIMENTAL: Normalize the aggregates */
  // NOTE: If we're normalizing aggregates we're going to exploit the fact that Psparse
  // is (a) not smoothed and (b) normalized.  This will enable us to get a plausible nullspace later on
  bool normalize_aggregates=MLAggr && List_.get("edge matrix free: normalize aggregates",false);
  if(verbose_ && !Comm_->MyPID() && normalize_aggregates) printf("EMFP: Normalizing aggregates\n");

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
        else vals2[j*dim+k]= normalize_aggregates ? (nullspace[k][i] * vals1[j] / 2.0) : (nullspace[k][i] / nonzeros);
      }/*end for*/
    }/*end for*/
    Prolongator_->InsertGlobalValues(EdgeRangeMap_->GID(i),dim*ne1,vals2,idx2);
  }/*end for*/

  /* FillComplete / OptimizeStorage for Prolongator*/
  Prolongator_->FillComplete(*CoarseMap_,*EdgeRangeMap_);
  Prolongator_->OptimizeStorage();

  /* Build the coarse nullspace (only works if we're normalizing aggregates) */
  bool build_coarse_nullspace= normalize_aggregates && List_.get("edge matrix free: explicit coarse nullspace",false);
  if(build_coarse_nullspace) {    
    if(verbose_ && !Comm_->MyPID()) printf("EMFP: Using explicit nullspace\n");
    int Ncoarse = Prolongator_->DomainMap().NumMyElements();
    CoarseNullspace_ = new double[Ncoarse*dim];
    memset(CoarseNullspace_,0,Ncoarse*dim*sizeof(double));

    // Use Psparse to migrate the constant and then split it up between dofs
    Epetra_Vector ones(Psparse->RangeMap(),false), output(Psparse->DomainMap(),true);
    Psparse->Multiply(true,ones,output);
    for(int i=0; i<Psparse->DomainMap().NumMyElements(); i++) {
      // Remember, the vector comes first
      for(int j=0; j<dim; j++)
        CoarseNullspace_[Ncoarse*j + i*dim + j] = 1.0;
    }    
    List_.sublist("edge matrix free: coarse").set("null space: dimension",dim);
    List_.sublist("edge matrix free: coarse").set("null space: vectors",CoarseNullspace_);
    List_.sublist("edge matrix free: coarse").set("null space: type","pre-computed");
    List_.sublist("edge matrix free: coarse").set("null space: add default vectors",false);
  }



  /* EXPERIMENTAL: Normalize Prolongator Columns */
  bool normalize_prolongator=List_.get("edge matrix free: normalize prolongator",false);
  if(normalize_prolongator){
    Epetra_Vector n_vector(*CoarseMap_,false);
    Prolongator_->InvColSums(n_vector);
    Prolongator_->RightScale(n_vector);
  }/*end if*/

  /* Post-wrapping to convert to ML indexing */
#ifdef HAVE_ML_EPETRAEXT
  Prolongator_ = dynamic_cast<Epetra_CrsMatrix*>(ModifyEpetraMatrixColMap(*Prolongator_,ProlongatorColMapTrans_,"Prolongator",(verbose_&&!Comm_->MyPID())));
#endif

  /* Build the Coarse Coordinates according to the colum map*/
  if(build_coarse_coords) { 
    /* Set coordinates on ListCoarse */
    Teuchos::ParameterList & ListCoarse=List_.sublist("edge matrix free: coarse");
    if(xcoord) ListCoarse.set("x-coordinates",CoarseCoord_.x.getRawPtr());
    if(ycoord) ListCoarse.set("y-coordinates",CoarseCoord_.y.getRawPtr());
    if(zcoord) ListCoarse.set("z-coordinates",CoarseCoord_.z.getRawPtr());
    if(mcoord) ListCoarse.set("material coordinates",CoarseCoord_.material.getRawPtr());
  }


  /* Cleanup */
  if(MLAggr) ML_Aggregate_Destroy(&MLAggr);
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
  ML_CHK_ERR(ML_Operator_WrapEpetraCrsMatrix(Prolongator_,P,verbose_));
  P->num_rigid=P->num_PDEs=dim;

  //NTS: ML_CHK_ERR won't work on this: it returns 1
  ML_Operator_Transpose_byrow(P, R);
#else
  if(verbose_ && !Comm_->MyPID()) printf("EMFP: Running FAST_PTAP\n");
#endif

  /* OPTION: Disable the addon */
  const Epetra_CrsMatrix *Op11crs = dynamic_cast<const Epetra_CrsMatrix*>(&*Operator_);
  const Epetra_Operator_With_MatMat *Op11mm = dynamic_cast<const Epetra_Operator_With_MatMat*>(&*Operator_);
#ifndef ENABLE_FAST_PTAP
  if(Op11crs){
    if(verbose_ && !Comm_->MyPID()) printf("EMFP: Running *without* addon\n");
    ML_Operator *SM_ML = ML_Operator_Create(ml_comm_);
    Temp_ML = ML_Operator_Create(ml_comm_);
    ML_Operator_WrapEpetraCrsMatrix((Epetra_CrsMatrix*)Op11crs,SM_ML,verbose_);
    ML_2matmult(SM_ML,P,Temp_ML,ML_CSR_MATRIX);
    ML_Operator_Destroy(&SM_ML);
  }
  else if(Op11mm){
#endif
#ifdef ENABLE_FAST_PTAP
    Op11mm->PtAP(*Prolongator_,ml_comm_,&CoarseMat_ML);
#else
    /* Do the A*P */
    if(verbose_ && !Comm_->MyPID()) printf("EMFP: Running with addon\n");
    ML_CHK_ERR(Op11mm->MatrixMatrix_Multiply(*Prolongator_,ml_comm_,&Temp_ML));
  }
  else{
    if(!Comm_->MyPID()) printf("EMFP: Error Invalid Operator\n");
    ML_CHK_ERR(-1);
  }
#endif

#ifndef ENABLE_FAST_PTAP
  /* Do R * AP */
  R->num_rigid=R->num_PDEs=dim;
  ML_2matmult_block(R, Temp_ML,CoarseMat_ML,ML_CSR_MATRIX);
#endif


  /* Wrap to Epetra-land */
  //  Epetra_CrsMatrix_Wrap_ML_Operator(CoarseMat_ML,*Comm_,*CoarseMap_,&CoarseMatrix);
  int nnz=100;
  double time;
  ML_Operator2EpetraCrsMatrix(CoarseMat_ML,CoarseMatrix,nnz,true,time,0,very_verbose_);
  // NTS: This is a hack to get around the sticking ones on the diagonal issue;

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
void ML_Epetra::EdgeMatrixFreePreconditioner::Print(int /* whichHierarchy */)
{
  /*ofstream ofs("Pmat.edge.m");
    if(Prolongator_) Prolongator_->Print(ofs);*/
  if(Prolongator_) EpetraExt::RowMatrixToMatlabFile("prolongator.dat",*Prolongator_);
  if(CoarseMatrix) EpetraExt::RowMatrixToMatlabFile("coarsemat.dat",*CoarseMatrix);


  if(CoarsePC) CoarsePC->Print();
}/*end Print*/

// ================================================ ====== ==== ==== == =
// Return operator complexity and #nonzeros in fine grid matrix.
void ML_Epetra::EdgeMatrixFreePreconditioner::Complexities(double &complexity, double &fineNnz){
  fineNnz=0.0;  complexity=1.0;

  if(!Operator_.is_null()) {
    const Epetra_RowMatrix * rm = dynamic_cast<const Epetra_RowMatrix*>(&*Operator_);
    if(!rm) return;

    fineNnz=rm->NumGlobalNonzeros();
    double coarse_oc=0.0, coarse_nnz=0.0;
    if(CoarsePC) CoarsePC->Complexities(coarse_oc,coarse_nnz);

    complexity = 1.0 + coarse_oc*coarse_nnz / fineNnz;
  }
}/*end Complexities*/

// ================================================ ====== ==== ==== == =
// Destroys all structures allocated in \c ComputePreconditioner() if the preconditioner has been computed.
int ML_Epetra::EdgeMatrixFreePreconditioner::DestroyPreconditioner(){
  if (ml_comm_) { ML_Comm_Destroy(&ml_comm_); ml_comm_ = 0; }// will need this
  if (Prolongator_) {delete Prolongator_; Prolongator_=0;}
  if (InvDiagonal_) {delete InvDiagonal_; InvDiagonal_=0;}
  if (CoarsePC) {delete CoarsePC; CoarsePC=0;}
  if (CoarseMatrix) {delete CoarseMatrix; CoarseMatrix=0;}
  if (CoarseNullspace_) {delete [] CoarseNullspace_; CoarseNullspace_=0;}
  if (CoarseMat_ML) {ML_Operator_Destroy(&CoarseMat_ML);CoarseMat_ML=0;}
  if (CoarseMap_) {delete CoarseMap_; CoarseMap_=0;}
#ifdef HAVE_ML_IFPACK
  if (Smoother_){delete Smoother_; Smoother_=0;}
#endif
  return 0;
}/*end DestroyPreconditioner*/

// ================================================ ====== ==== ==== == =
//! Apply the preconditioner to an Epetra_MultiVector X, puts the result in Y
int ML_Epetra::EdgeMatrixFreePreconditioner::ApplyInverse(const Epetra_MultiVector& B_, Epetra_MultiVector& X) const{
  const Epetra_MultiVector *B;
  Epetra_MultiVector *Bcopy=0;

  /* Sanity Checks */
  int NumVectors=B_.NumVectors();
  if (!B_.Map().SameAs(*EdgeDomainMap_)) ML_CHK_ERR(-1);
  if (NumVectors != X.NumVectors()) ML_CHK_ERR(-2);

  Epetra_MultiVector r_edge(*EdgeDomainMap_,NumVectors,false);
  Epetra_MultiVector e_edge(*EdgeDomainMap_,NumVectors,false);
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
