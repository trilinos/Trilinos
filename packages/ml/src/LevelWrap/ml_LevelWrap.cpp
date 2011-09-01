#include "ml_LevelWrap.h"
#include "ml_include.h"

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_IFPACK)
#include "ml_utils.h"
#include "ml_epetra_utils.h"
#include "Ifpack.h"
#include "ml_ValidateParameters.h"

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
  pre_or_post(ML_BOTH),
  verbose_(false),
  NumApplications_(0),
  ApplicationTime_(0.0),
  FirstApplication_(true),
  FirstApplicationTime_(0.0),
  NumConstructions_(0),
  ConstructionTime_(0.0)
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
  pre_or_post(ML_BOTH),
  verbose_(false),
  NumApplications_(0),
  ApplicationTime_(0.0),
  FirstApplication_(true),
  FirstApplicationTime_(0.0),
  NumConstructions_(0),
  ConstructionTime_(0.0)
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

  Teuchos::ParameterList & IfpackList=List_.sublist("smoother: ifpack list");
  string MyIfpackType = List_.get("smoother: type", "Gauss-Seidel");
  if(MyIfpackType=="IFPACK") MyIfpackType = List_.get("smoother: ifpack type", MyIfpackType);
  int MyIfpackOverlap = List_.get("smoother: ifpack overlap", 0);

  Ifpack Factory;
  Ifpack_Preconditioner* Prec= Factory.Create(MyIfpackType,const_cast<Epetra_CrsMatrix*>(&*A0_), MyIfpackOverlap);
  Prec->SetParameters(IfpackList);
  Prec->Compute();
  Smoother_=rcp<Ifpack_Preconditioner>(Prec);
  Smoother_->Print(cout);

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
  StopTimer(&t_time_curr,&t_diff);
  /* Output */
  ML_Comm *comm_;
  ML_Comm_Create(&comm_);
  int printl=ML_Get_PrintLevel();
  ML_Set_PrintLevel(output_level);  
  ReportTimer(t_diff,"ML_LevelWrap::ComputePreconditioner",comm_);
  ML_Set_PrintLevel(printl);
  ML_Comm_Destroy(&comm_);
  NumConstructions_++;
  ConstructionTime_+=t_time_curr-t_time_start;
#endif  
  return 0;
}


// Apply the preconditioner w/ RHS B and get result X
int ML_Epetra::LevelWrap::ApplyInverse(const Epetra_MultiVector& B, Epetra_MultiVector& X_) const{
#ifdef ML_TIMING
  double t_time,t_diff;
  StartTimer(&t_time);
#endif
   
  // Sanity Checks
  if (!B.Map().SameAs(OperatorDomainMap())) return -1;
  if (B.NumVectors() != X_.NumVectors()) return -1;

  // Build new work vector X 
  Epetra_MultiVector X(X_.Map(),X_.NumVectors(),true);
  Epetra_MultiVector tmp0(X_.Map(),X_.NumVectors(),true);
  Epetra_MultiVector tmp1(P0_->RangeMap(),X_.NumVectors(),true);
  Epetra_MultiVector tmp2(P0_->RangeMap(),X_.NumVectors(),true);
  
  // Pre Smoother
  if(pre_or_post==ML_BOTH || pre_or_post==ML_PRESMOOTHER){
    Smoother_->ApplyInverse(B,X);
  }

  // Form coarse residual
  A0_->Apply(X,tmp0);
  tmp0.Update(1.0,B,-1.0); 
  if(use_pt_) P0_->Multiply(true,tmp0,tmp1);
  else R0_->Multiply(false,tmp0,tmp1);

  // Solve coarse problem
  A1prec_->ApplyInverse(tmp1,tmp2);

  // Update solution
  P0_->Multiply(false,tmp2,tmp0);
  X.Update(1.0,tmp0,1.0);

  // Post Smoother
  if(pre_or_post==ML_BOTH || pre_or_post==ML_PRESMOOTHER){
    Smoother_->ApplyInverse(B,X);
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
  return 0;
}

#endif
