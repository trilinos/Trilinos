#include "ml_config.h"
#include "ml_LevelWrap.h"
#include "ml_include.h"

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_IFPACK)
#include "ml_utils.h"
#include "ml_epetra_utils.h"
#include "Ifpack.h"
#include "ml_ValidateParameters.h"
#include "ml_ifpack_epetra_wrap.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#endif

using Teuchos::rcp;
using Teuchos::RCP;
using std::string;

// Sets default parameters for aggregation-based 2-level domain decomposition preconditioners.
int ML_Epetra::SetDefaultsLevelWrap(Teuchos::ParameterList & inList,bool OverWrite){
  Teuchos::ParameterList ilist,&ilist_in=inList.sublist("smoother: ifpack list");
  Teuchos::ParameterList ListLW,sublist, &sublist_in=inList.sublist("levelwrap: coarse list");

  // sublist
  ML_Epetra::SetDefaults("SA",sublist);
  ML_Epetra::UpdateList(sublist,sublist_in,OverWrite);

  // ilist
  ilist.set("relaxation: type", "symmetric Gauss-Seidel");
  ilist.set("relaxation: zero starting solution",false);
  ilist.set("relaxation: sweeps",2);
  ML_Epetra::UpdateList(ilist,ilist_in,OverWrite);

  // overall list
  ML_Epetra::SetDefaults("SA",ListLW);
  ListLW.set("smoother: type","IFPACK");
  ListLW.set("smoother: ifpack type","point relaxation stand-alone");
  ListLW.set("levelwrap: coarse list",sublist_in);
  ListLW.set("smoother: ifpack list",ilist_in);

  ML_Epetra::UpdateList(ListLW,inList,OverWrite);
  return 0;
}

// Constructs a Level Wrap (using R=P^T)
ML_Epetra::LevelWrap::LevelWrap(Teuchos::RCP<Epetra_CrsMatrix> A0,
                                Teuchos::RCP<Epetra_CrsMatrix> P0,
                                const Teuchos::ParameterList& List,
                                const bool ComputePrec):
  use_pt_(true),
  user_A1_(false),
  use_mlsmoother_(false),   // ML smoothing is experimental
  ml_subproblem_(NULL),
  pre_or_post(ML_BOTH),
  verbose_(false) /* (unused) ,
  NumApplications_(0),
  ApplicationTime_(0.0),
  FirstApplication_(true),
  FirstApplicationTime_(0.0),
  NumConstructions_(0),
  ConstructionTime_(0.0)
                  */
{
  A0_=A0;
  P0_=P0;
  List_=List;
  if(ComputePrec) ComputePreconditioner();
}


// Constructs a Level Wrap (using different R)
ML_Epetra::LevelWrap::LevelWrap(Teuchos::RCP<Epetra_CrsMatrix> A0,
                                Teuchos::RCP<Epetra_CrsMatrix> P0,
                                Teuchos::RCP<Epetra_CrsMatrix> R0,
                                const Teuchos::ParameterList& List,
                                const bool ComputePrec):
  use_pt_(false),
  user_A1_(false),
  use_mlsmoother_(false),   // ML smoothing is experimental
  ml_subproblem_(NULL),
  pre_or_post(ML_BOTH),
  verbose_(false) /* (unused) ,
  NumApplications_(0),
  ApplicationTime_(0.0),
  FirstApplication_(true),
  FirstApplicationTime_(0.0),
  NumConstructions_(0),
  ConstructionTime_(0.0)
                  */
{
  A0_=A0;
  P0_=P0;
  R0_=R0;
  List_=List;
  if(ComputePrec) ComputePreconditioner();
}


// Destructor
ML_Epetra::LevelWrap::~LevelWrap(){
}


// Computes the preconditioner
int ML_Epetra::LevelWrap::ComputePreconditioner(const bool CheckFiltering){
#ifdef ML_TIMING
  double t_time,t_diff;
  StartTimer(&t_time);
#endif
  int output_level=List_.get("ML output",0);
  if(output_level>0) verbose_=true;

  // Sanity check
  if(!A0_->DomainMap().SameAs(P0_->RangeMap())) return -1;
  if(!use_pt_ && !A0_->RangeMap().SameAs(R0_->DomainMap())) return -2;

  //********************
  // Setup smoother
  //********************
  string PreOrPostSmoother = List_.get("smoother: pre or post","both");
  if(PreOrPostSmoother == "post") pre_or_post = ML_POSTSMOOTHER;
  else if(PreOrPostSmoother == "pre")  pre_or_post = ML_PRESMOOTHER;
  else if(PreOrPostSmoother == "both") pre_or_post = ML_BOTH;

  if (use_mlsmoother_) {   // ML smoothing is experimental
    ML_Create(&ml_subproblem_,1);
#ifdef HAVE_MPI
    const Epetra_MpiComm *epcomm = dynamic_cast<const Epetra_MpiComm*>(&(A0_->Comm()));
    // Get the MPI communicator, as it may not be MPI_COMM_W0RLD, and update the ML comm object
   if (epcomm) ML_Comm_Set_UsrComm(ml_subproblem_->comm,epcomm->Comm());
#endif
    int numMyRows = A0_->NumMyRows();
    int N_ghost   = A0_->NumMyCols() - numMyRows;
    if (N_ghost < 0) N_ghost = 0;
    const Epetra_Operator *Z = const_cast<Epetra_CrsMatrix*>(&*A0_);
    const Epetra_RowMatrix *Arow=dynamic_cast<const Epetra_RowMatrix*>(Z);
    ML_Init_Amatrix(ml_subproblem_,0,numMyRows, numMyRows,(void *) Arow);
    ml_subproblem_->Amat[0].type = ML_TYPE_ROW_MATRIX;
    ml_subproblem_->Amat[0].N_nonzeros = A0_->NumMyNonzeros();
    ML_Set_Amatrix_Getrow(ml_subproblem_, 0, ML_Epetra_RowMatrix_getrow,
                          ML_Epetra_comm_wrapper, numMyRows+N_ghost);
    ML_Set_Amatrix_Matvec(ml_subproblem_, 0, ML_Epetra_matvec);
    int Nparts =  List_.sublist("smoother: ifpack list").get("partitioner: local parts",-1);
    if (Nparts == -1) { printf("must supply partitioner: local parts\n"); return(-1); }

    ml_subproblem_->Amat[0].num_PDEs = numMyRows/Nparts;
    ML_Gen_Smoother_BlockGaussSeidel(ml_subproblem_, 0, ML_PRESMOOTHER,
                       1, 1.0 , ml_subproblem_->Amat[0].num_PDEs);
  }
  else Smoother_=rcp(ML_Gen_Smoother_Ifpack_Epetra(const_cast<Epetra_CrsMatrix*>(&*A0_),0,List_,"LevelWrap Smoother (level 0): ",verbose_));



  //********************
  // Setup A1
  //********************
  Epetra_CrsMatrix * A1=List_.get("user coarse matrix",(Epetra_CrsMatrix *)0);
  if(A1){
    // User A1
    user_A1_=true;
    A1_=rcp<Epetra_CrsMatrix>(A1,false);
  }
  else{
    // Do RAP
    if(use_pt_) ML_Epetra_PtAP(*A0_,*P0_,A1,false);
    else ML_Epetra_RAP(*A0_,*P0_,*R0_,A1,false);
    A1_=rcp<Epetra_CrsMatrix>(A1);
  }

  //********************
  // Setup Coarse Solver
  //********************
  Teuchos::ParameterList & CoarseList=List_.sublist("levelwrap: coarse list");
  A1prec_=rcp<MultiLevelPreconditioner>(new MultiLevelPreconditioner(*A1_,CoarseList,true));

#ifdef ML_TIMING
  StopTimer(&t_time,&t_diff);
  /* Output */
  ML_Comm *comm_;
  ML_Comm_Create(&comm_);
  int printl=ML_Get_PrintLevel();
  ML_Set_PrintLevel(output_level);
  ReportTimer(t_diff,"ML_LevelWrap::ComputePreconditioner",comm_);
  ML_Set_PrintLevel(printl);
  ML_Comm_Destroy(&comm_);
  NumConstructions_++;
  ConstructionTime_+=t_diff;
#endif
  return 0;
}

// ================================================ ====== ==== ==== == =
// Return operator complexity and #nonzeros in fine grid matrix.
void ML_Epetra::LevelWrap::Complexities(double &complexity, double &fineNnz){
  fineNnz= (!A0_.is_null()) ? A0_->NumGlobalNonzeros() : 0.0;
  complexity=1.0;

  if(!A1prec_.is_null()) {
    double coarse_oc=0.0, coarse_nnz=0.0;
    A1prec_->Complexities(coarse_oc,coarse_nnz);
    complexity = 1.0 + coarse_oc*coarse_nnz / fineNnz;
  }
}/*end Complexities */


// ================================================ ====== ==== ==== == =
// Apply the preconditioner w/ RHS B and get result X
int ML_Epetra::LevelWrap::ApplyInverse(const Epetra_MultiVector& B, Epetra_MultiVector& X_) const{
#ifdef ML_TIMING
  double t_time,t_diff;
  StartTimer(&t_time);
#endif

  // Sanity Checks
  if (!B.Map().SameAs(OperatorDomainMap())) return -1;
  if (!X_.Map().SameAs(OperatorRangeMap())) return -1;
  if (!X_.Map().SameAs(B.Map())) return -1;
  if (B.NumVectors() != X_.NumVectors()) return -1;

  // Build new work vector X
  Epetra_MultiVector X(X_.Map(),X_.NumVectors(),true);
  Epetra_MultiVector tmp0(X_.Map(),X_.NumVectors(),true);
  Epetra_MultiVector tmp1(P0_->DomainMap(),X_.NumVectors(),true);
  Epetra_MultiVector tmp2(P0_->DomainMap(),X_.NumVectors(),true);

  double *Bptr, *Xptr;
  int    BLDA, XLDA;

  // Pre Smoother
  if(pre_or_post==ML_BOTH || pre_or_post==ML_PRESMOOTHER){
    if (use_mlsmoother_) {   // ML smoothing is experimental
      B.ExtractView(&Bptr,&BLDA);
      X.ExtractView(&Xptr,&XLDA);
      ML_Smoother_Apply(&(ml_subproblem_->pre_smoother[0]),
                        X.MyLength(), Xptr,
                        X.MyLength(), Bptr, ML_ZERO);
    }
    else Smoother_->ApplyInverse(B,X);

    A0_->Apply(X,tmp0);
    tmp0.Update(1.0,B,-1.0);
  }
  else  tmp0 = B;

  // Form coarse residual
  if(use_pt_) P0_->Multiply(true,tmp0,tmp1);
  else R0_->Multiply(false,tmp0,tmp1);

  // Solve coarse problem
  A1prec_->ApplyInverse(tmp1,tmp2);

  // Update solution
  P0_->Multiply(false,tmp2,tmp0);
  X.Update(1.0,tmp0,1.0);

  // Post Smoother
  if(pre_or_post==ML_BOTH || pre_or_post==ML_POSTSMOOTHER){
    if (use_mlsmoother_) {   // ML smoothing is experimental
      Epetra_MultiVector tmp3(X_.Map(),X_.NumVectors(),true);

      A0_->Apply(X,tmp0);
      tmp0.Update(1.0,B,-1.0);
      tmp3.PutScalar(0.0);
      tmp0.ExtractView(&Bptr,&BLDA);
      tmp3.ExtractView(&Xptr,&XLDA);

      ML_Smoother_Apply(&(ml_subproblem_->pre_smoother[0]),X.MyLength(),Xptr,
                        X.MyLength(), Bptr, ML_ZERO);
      X.Update(1.0,tmp3,1.0);
    }
    else Smoother_->ApplyInverse(B,X);
  }

  // Copy to output
  X_=X;

#ifdef ML_TIMING
  StopTimer(&t_time,&t_diff);
  /* Output */
  ML_Comm *comm_;
  ML_Comm_Create(&comm_);
  ApplicationTime_+= t_diff;
  if(FirstApplication_){
    FirstApplication_=false;
    FirstApplicationTime_=ApplicationTime_;
  }/*end if*/
  ML_Comm_Destroy(&comm_);
#endif

  return 0;
}


// Print the individual operators in the multigrid hierarchy.
void ML_Epetra::LevelWrap::Print(int level){
  if(A1prec_!=Teuchos::null) A1prec_->Print(-1);
}

// Destroys all structures allocated in \c ComputePreconditioner() if the preconditioner has been computed.
int ML_Epetra::LevelWrap::DestroyPreconditioner(){
  int printl=ML_Get_PrintLevel();
  int output_level=List_.get("ML output",0);
  output_level=List_.get("output",output_level);
  ML_Set_PrintLevel(output_level);

  Smoother_=Teuchos::null;
  A0_=Teuchos::null;
  P0_=Teuchos::null;
  A1prec_=Teuchos::null;
  A1_=Teuchos::null;


#ifdef ML_TIMING
  ML_Comm *comm_;
  ML_Comm_Create(&comm_);
  ReportTimer(ConstructionTime_ ,   "ML_LevelWrap (construction  )",comm_);
  ReportTimer(FirstApplicationTime_,"ML_LevelWrap (1st iter time )",comm_);
  ReportTimer(ApplicationTime_ ,    "ML_LevelWrap (total itr cost)",comm_);
  ML_Set_PrintLevel(printl);
  ML_Comm_Destroy(&comm_);
#else
  ML_Set_PrintLevel(printl);
#endif

  if (ml_subproblem_ != NULL) ML_Destroy(&ml_subproblem_);
  return 0;
}

#endif
