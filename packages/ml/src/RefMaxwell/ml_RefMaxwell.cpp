/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */
#include "ml_config.h"
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_EPETRAEXT) 
#include <iostream>
#include <vector>
#include <string.h>
#include "ml_RefMaxwell.h"
#include "ml_epetra.h"
#include "ml_epetra_utils.h"
#include "ml_MultiLevelPreconditioner.h"
#include "ml_RefMaxwell_11_Operator.h"
#include "ml_RefMaxwell_Utils.h"
#include "ml_EdgeMatrixFreePreconditioner.h"
#include "ml_ValidateParameters.h"
#include "Teuchos_ArrayRCP.hpp"

#include "EpetraExt_RowMatrixOut.h"
#include "ml_ifpack_epetra_wrap.h"


using Teuchos::rcp;
using Teuchos::ArrayRCP;

#ifdef HAVE_ML_IFPACK
#include "Ifpack.h"
#endif


// ============================================================================
int* FindLocalDiricheltAndInflowRowsFromOnesAndZeros(const Epetra_CrsMatrix & Matrix, int &numBCRows){
  int *dirichletRows = new int[Matrix.NumMyRows()];
  numBCRows = 0;
  for (int i=0; i<Matrix.NumMyRows(); i++) {
    int numEntries, *cols;
    double *vals;
    int ierr = Matrix.ExtractMyRowView(i,numEntries,vals,cols);
    if (ierr == 0) {
      int nz=0;
      for (int j=0; j<numEntries; j++) if (vals[j] != 0.0) nz++;
      if (nz == 1) dirichletRows[numBCRows++] = i;      

      // EXPERIMENTAL: Treat Special Inflow Boundaries as Dirichlet Boundaries
      if(nz==2) dirichletRows[numBCRows++] = i;  
	
    }/*end if*/
  }/*end fpr*/
  return dirichletRows;
}/*end FindLocalDiricheltRowsFromOnesAndZeros*/


// ================================================ ====== ==== ==== == = 
ML_Epetra::RefMaxwellPreconditioner::RefMaxwellPreconditioner(const Epetra_CrsMatrix& SM_Matrix,      //S+M
                                                              const Epetra_CrsMatrix& D0_Clean_Matrix,//T or D0 w/ nothing zero'd
#ifdef ENABLE_MS_MATRIX
                                                              const Epetra_CrsMatrix& Ms_Matrix,      //M1(sigma)
#endif
                                                              const Epetra_CrsMatrix& M0inv_Matrix,   //M0^{-1}
                                                              const Epetra_CrsMatrix& M1_Matrix,      //M1(1)
                                                              const Teuchos::ParameterList& List,
                                                              const bool ComputePrec):
  ML_Preconditioner(),SM_Matrix_(&SM_Matrix),D0_Matrix_(0), D0_Clean_Matrix_(&D0_Clean_Matrix),
#ifdef ENABLE_MS_MATRIX
  Ms_Matrix_(&Ms_Matrix),
#endif
  M0inv_Matrix_(&M0inv_Matrix),M1_Matrix_((Epetra_CrsMatrix*)&M1_Matrix),TMT_Matrix_(0),TMT_Agg_Matrix_(0),
  BCrows(0),numBCrows(0),HasOnlyDirichletNodes(false),Operator11_(0),EdgePC(0),NodePC(0),PreEdgeSmoother(0),PostEdgeSmoother(0),
#ifdef HAVE_ML_IFPACK
  IfSmoother(0),
#endif
  aggregate_with_sigma(false),lump_m1(false),
  use_local_nodal_solver(false),LocalNodalSolver(0),NodesToLocalNodes(0),LocalNodalMatrix(0),
  verbose_(false),very_verbose_(false)
{
  /* Set the Epetra Goodies */
  Comm_ = &(SM_Matrix_->Comm());
  DomainMap_ = &(SM_Matrix_->OperatorDomainMap());
  RangeMap_ = &(SM_Matrix_->OperatorRangeMap());
  NodeMap_ = &(D0_Clean_Matrix_->OperatorDomainMap());

  Label_=new char [80];
  strcpy(Label_,"ML reformulated Maxwell preconditioner");
  List_=List;
  SetDefaultsRefMaxwell(List_,false);
  
#ifdef ML_TIMING
  /* Internal Timings */
  NumApplications_ = 0;
  ApplicationTime_ = 0.0;
  FirstApplication_ = true;
  FirstApplicationTime_ = 0.0;
  NumConstructions_ = 0;
  ConstructionTime_ = 0.0;
#endif

  if(ComputePrec) ML_CHK_ERRV(ComputePreconditioner());
}/*end constructor*/


// ================================================ ====== ==== ==== == = 
ML_Epetra::RefMaxwellPreconditioner::~RefMaxwellPreconditioner()
{
  if (IsComputePreconditionerOK_) 
    DestroyPreconditioner(); 
}/*end destructor*/


// ================================================ ====== ==== ==== == = 
// Print the individual operators in the multigrid hierarchy.
void ML_Epetra::RefMaxwellPreconditioner::Print(int whichHierarchy){
  if(IsComputePreconditionerOK_ && EdgePC && whichHierarchy==11) EdgePC->Print(-1);
  if(IsComputePreconditionerOK_ && NodePC && whichHierarchy==22) NodePC->Print(-1);  
}/*end Print*/


// ================================================ ====== ==== ==== == = 
// Computes the preconditioner
int ML_Epetra::RefMaxwellPreconditioner::ComputePreconditioner(const bool CheckFiltering)
{
#ifdef ML_TIMING
  double t_time_start, t_time_curr, t_diff[7];
  StartTimer(&t_time_start);
  t_time_curr=t_time_start;
#endif

  int output_level=List_.get("ML output",0);
  output_level=List_.get("output",output_level);
  
  Teuchos::ParameterList dummy;

  /* Validate List */
  Teuchos::ParameterList newList;
  ML_CreateSublists(List_,newList);
  List_ = newList;
  //TODO check for failure of validation, and print a helpful message,
  //     just like in MultiLevelPreconditioner
  ValidateRefMaxwellParameters(List_);
  
  /* Pull Solver Mode, verbosity, matrix output */
  mode=List_.get("refmaxwell: mode","additive");
  print_hierarchy= List_.get("print hierarchy",-2);  
  int vb_level=List_.get("ML output",0);
  if(vb_level >= 15) {very_verbose_=true;verbose_=true;}
  else if (vb_level >= 5) {very_verbose_=false;verbose_=true;}
  else very_verbose_=verbose_=false;
  aggregate_with_sigma= List_.get("refmaxwell: aggregate with sigma",false);  
  
  /* Nuke everything if we've done this already */
  if(IsComputePreconditionerOK_) DestroyPreconditioner();

  /* Output Header */
  if(verbose_ && !Comm_->MyPID()) {
    printf("------------------------------------------------------------------------------\n");
    printf("***\n");
    printf("*** ML_Epetra::RefMaxwellPreconditioner [%s]\n",Label());
    printf("***\n");
  }

  /* Find the Dirichlet Rows (using SM_Matrix_) and columns (using D0_Clean_Matrix_) */
  BCrows=FindLocalDiricheltAndInflowRowsFromOnesAndZeros(*SM_Matrix_,numBCrows);
  Epetra_IntVector * BCnodes=FindLocalDirichletColumnsFromRows(BCrows,numBCrows,*D0_Clean_Matrix_);   
  int Nn=BCnodes->MyLength();
  int numBCnodes=0;
  for(int i=0;i<Nn;i++){
    if((*BCnodes)[i]) numBCnodes++;
  }
  
  /* Sanity Check: We have at least some Dirichlet nodes */
  int HasInterior = numBCnodes != Nn;
  if(very_verbose_ && !Comm_->MyPID()) {
    printf("[%d] HasInterior = %d %d/%d\n",Comm_->MyPID(),HasInterior,Nn-numBCnodes,Nn);
    printf("Edge Matrix: Unknowns = %d Nonzeros = %d\n",SM_Matrix_->NumGlobalRows(),SM_Matrix_->NumGlobalNonzeros());
  }
  int globalInterior=HasInterior;
  Comm_->MaxAll(&HasInterior,&globalInterior,1);  
  if(!globalInterior){
    HasOnlyDirichletNodes=true;
    if(!Comm_->MyPID()) printf("WARNING: All nodes are Dirichlet nodes.\n");
  }/*end if*/

  /* Do the Nuking for D0_Matrix_ */ 
  D0_Matrix_ = new Epetra_CrsMatrix(*D0_Clean_Matrix_);
  Apply_BCsToMatrixRows(BCrows,numBCrows,*D0_Matrix_);
  Apply_BCsToMatrixColumns(*BCnodes,*D0_Matrix_);   
  D0_Matrix_->OptimizeStorage();

  /* EXPERIMENTAL - Lump M1, if requested */
  lump_m1 = List_.get("refmaxwell: lump m1",false);
  if(lump_m1){
    Epetra_Vector mvec1(*RangeMap_,false), mvec2(*RangeMap_,false);
    mvec1.PutScalar(1.0);
    M1_Matrix_->Multiply(false,mvec1,mvec2);
    M1_Matrix_ = new Epetra_CrsMatrix(Copy,*RangeMap_,1);
    for(int i=0;i<mvec2.MyLength();i++){
       int idx = RangeMap_->GID(i);     
       M1_Matrix_->InsertGlobalValues(idx,1,&(mvec2[i]),&idx);    
     }
    M1_Matrix_->FillComplete();
    M1_Matrix_->OptimizeStorage();
  }/*end if*/  


#ifdef ML_TIMING
  StopTimer(&t_time_curr,&(t_diff[0]));
#endif
  
  /* Build the TMT Matrix */
  if(!HasOnlyDirichletNodes){
    //    ML_Epetra_PtAP(*Ms_Matrix_,*D0_Matrix_,TMT_Matrix_,verbose_);
    ML_Epetra_PtAP(*SM_Matrix_,*D0_Matrix_,TMT_Matrix_,verbose_);
    Remove_Zeroed_Rows(*TMT_Matrix_,1e-10);
  }/*end if */
  
  /* Build the TMT-Agg Matrix */
  if(aggregate_with_sigma)
#ifdef ENABLE_MS_MATRIX
    ML_Epetra_PtAP(*Ms_Matrix_,*D0_Clean_Matrix_,TMT_Agg_Matrix_,verbose_);
#else
  ML_Epetra_PtAP(*SM_Matrix_,*D0_Clean_Matrix_,TMT_Agg_Matrix_,verbose_);
#endif
  else ML_Epetra_PtAP(*M1_Matrix_,*D0_Clean_Matrix_,TMT_Agg_Matrix_,verbose_);
  Remove_Zeroed_Rows(*TMT_Agg_Matrix_);
  
#ifdef ML_TIMING
  StopTimer(&t_time_curr,&(t_diff[1]));
#endif
  
  /* Boundary nuke the edge matrices */
#ifdef ENABLE_MS_MATRIX
  Apply_OAZToMatrix(BCrows,numBCrows,*Ms_Matrix_);
#endif
  Apply_OAZToMatrix(BCrows,numBCrows,*M1_Matrix_);    


  /* DEBUG: Output matrices */
  if(print_hierarchy == -1 || print_hierarchy == 0){
    if(verbose_ && !Comm_->MyPID()) printf("Dumping Matrices to Disk\n");
    EpetraExt::RowMatrixToMatlabFile("sm_matrix.dat",*SM_Matrix_);
#ifdef ENABLE_MS_MATRIX
    EpetraExt::RowMatrixToMatlabFile("ms_matrix.dat",*Ms_Matrix_);
#endif
    EpetraExt::RowMatrixToMatlabFile("ms1_nuked.dat",*M1_Matrix_);
    EpetraExt::RowMatrixToMatlabFile("ms0inv_nuked.dat",*M0inv_Matrix_);
    EpetraExt::RowMatrixToMatlabFile("d0_nuked.dat",*D0_Matrix_);
    EpetraExt::RowMatrixToMatlabFile("d0_clean.dat",*D0_Clean_Matrix_);
    if(TMT_Matrix_) EpetraExt::RowMatrixToMatlabFile("tmt_matrix.dat",*TMT_Matrix_);
    EpetraExt::RowMatrixToMatlabFile("tmt_agg_matrix.dat",*TMT_Agg_Matrix_);

  }

  /* Cleanup from the Boundary Conditions */
  delete BCnodes; 
  
#ifdef HAVE_ML_EPETRAEXT
  /* Fix the solver maps for ML / Epetra compatibility */
  SM_Matrix_ = dynamic_cast<Epetra_CrsMatrix*>(ModifyEpetraMatrixColMap(*SM_Matrix_,SM_Matrix_Trans_,"SM",(verbose_&&!Comm_->MyPID())));
  D0_Matrix_ = dynamic_cast<Epetra_CrsMatrix*>(ModifyEpetraMatrixColMap(*D0_Matrix_,D0_Matrix_Trans_,"D0",(verbose_&&!Comm_->MyPID())));
  D0_Clean_Matrix_ = dynamic_cast<Epetra_CrsMatrix*>(ModifyEpetraMatrixColMap(*D0_Clean_Matrix_,D0_Clean_Matrix_Trans_,"D0Clean",(verbose_&&!Comm_->MyPID())));
#ifdef ENABLE_MS_MATRIX
  if(Ms_Matrix_!= M1_Matrix_) Ms_Matrix_ = dynamic_cast<Epetra_CrsMatrix*>(ModifyEpetraMatrixColMap(*Ms_Matrix_,Ms_Matrix_Trans_,"Ms",(verbose_&&!Comm_->MyPID())));
#endif
  M1_Matrix_ = dynamic_cast<Epetra_CrsMatrix*>(ModifyEpetraMatrixColMap(*M1_Matrix_,M1_Matrix_Trans_,"M1",(verbose_&&!Comm_->MyPID())));
  M0inv_Matrix_ = dynamic_cast<Epetra_CrsMatrix*>(ModifyEpetraMatrixColMap(*M0inv_Matrix_,M0inv_Matrix_Trans_,"M0inv",(verbose_&&!Comm_->MyPID())));
  if(TMT_Matrix_) TMT_Matrix_ = dynamic_cast<Epetra_CrsMatrix*>(ModifyEpetraMatrixColMap(*TMT_Matrix_,TMT_Matrix_Trans_,"TMT",(verbose_&&!Comm_->MyPID()))); 
  TMT_Agg_Matrix_ = dynamic_cast<Epetra_CrsMatrix*>(ModifyEpetraMatrixColMap(*TMT_Agg_Matrix_,TMT_Agg_Matrix_Trans_,"TMTA",(verbose_&&!Comm_->MyPID())));
#endif

#ifdef ML_TIMING
  StopTimer(&t_time_curr,&(t_diff[2]));
#endif
  
  /* Build the (1,1) Block Operator, if needed */
  if(List_.get("refmaxwell: disable addon",true))
    Operator11_=rcp((Epetra_CrsMatrix*)SM_Matrix_,false);
  else
    Operator11_=rcp(new ML_RefMaxwell_11_Operator(*SM_Matrix_,*D0_Matrix_,*M0inv_Matrix_,*M1_Matrix_));

#ifdef ML_TIMING
  StopTimer(&t_time_curr,&(t_diff[3]));
#endif

  /* This is just placeholder code that should probably be removed */
  Diagonal_=new Epetra_Vector(*DomainMap_,false);
  Diagonal_->PutScalar(1.0);

#ifdef ML_TIMING
  StopTimer(&t_time_curr,&(t_diff[4]));
#endif

  // BC edges
  ArrayRCP<int> BCedges_arcp(BCrows,0,numBCrows,false);
  
  /* Build the (1,1) Block Preconditioner */ 
  string solver11=List_.get("refmaxwell: 11solver","edge matrix free");
  Teuchos::ParameterList List11=List_.get("refmaxwell: 11list",dummy);
  if (List11.name() == "ANONYMOUS") List11.setName("refmaxwell: 11list");
  if(solver11=="edge matrix free")
    //    EdgePC=new EdgeMatrixFreePreconditioner(*Operator11_,*Diagonal_,*D0_Matrix_,*D0_Clean_Matrix_,*TMT_Agg_Matrix_,BCrows,numBCrows,List11,true);
    EdgePC=new EdgeMatrixFreePreconditioner(Operator11_,rcp(Diagonal_,false),rcp(D0_Matrix_,false),rcp(D0_Clean_Matrix_,false),rcp(TMT_Agg_Matrix_,false),BCedges_arcp,List11,true);

  else {printf("RefMaxwellPreconditioner: ERROR - Illegal (1,1) block preconditioner\n");return -1;}
#ifdef ML_TIMING
  StopTimer(&t_time_curr,&(t_diff[5]));
#endif
  if(print_hierarchy) EdgePC->Print();
  /* Build the (2,2) Block Preconditioner */
  if(!HasOnlyDirichletNodes){
    string solver22=List_.get("refmaxwell: 22solver","multilevel");
    Teuchos::ParameterList List22=List_.get("refmaxwell: 22list",dummy);
    if (List22.name() == "ANONYMOUS") List11.setName("refmaxwell: 11list");
    SetDefaults("SA",List22,0,0,false);
    if(solver22=="multilevel") NodePC=new MultiLevelPreconditioner(*TMT_Matrix_,List22);
    else {printf("RefMaxwellPreconditioner: ERROR - Illegal (2,2) block preconditioner\n");return -1;}
    //NTS: Add Adaptive, MatrixFree
    if(print_hierarchy == -1) NodePC->Print(print_hierarchy);
  }/*end if*/
 
  /* Setup the Edge Smoother */
  //  if(mode=="additive")
  SetEdgeSmoother(List_);

#ifdef ML_TIMING
  StopTimer(&t_time_curr,&(t_diff[6]));
  /* Output */
  ML_Comm *comm_;
  ML_Comm_Create(&comm_);
  int printl=ML_Get_PrintLevel();
  ML_Set_PrintLevel(output_level);  
  ReportTimer(t_diff[0],"ML_RMP::ComputePreconditioner (BCs / nukes 1 )",comm_);
  ReportTimer(t_diff[1],"ML_RMP::ComputePreconditioner (TMT builds    )",comm_);
  ReportTimer(t_diff[2],"ML_RMP::ComputePreconditioner (remaps/nukes 2)",comm_);
  ReportTimer(t_diff[3],"ML_RMP::ComputePreconditioner (build 1,1 op  )",comm_);
  ReportTimer(t_diff[4],"ML_RMP::ComputePreconditioner (diagonal flip )",comm_);
  ReportTimer(t_diff[5],"ML_RMP::ComputePreconditioner (build 1,1 prec)",comm_);
  ReportTimer(t_diff[6],"ML_RMP::ComputePreconditioner (build 2,2 prec)",comm_);
  ReportTimer(t_time_curr-t_time_start,"ML_RMP::ComputePreconditioner (total         )",comm_);
  ML_Set_PrintLevel(printl);
  ML_Comm_Destroy(&comm_);

  NumConstructions_++;
  ConstructionTime_+=t_time_curr-t_time_start;
#endif  

  /* Local Nodal Stuff, if active */
  use_local_nodal_solver = List_.get("refmaxwell: enable local nodal solver",false);
  NodesToLocalNodes=List_.get("refmaxwell: global to local nodal transfer matrix",(Epetra_CrsMatrix*)0);
  if(!NodesToLocalNodes) use_local_nodal_solver=false;
  if(use_local_nodal_solver){
    if(verbose_ && !Comm_->MyPID()) printf("RefMaxwell: Local Nodal Solver ENABLED\n");
    // EXPERIMENTAL: Local nodal solver
    ML_Epetra::ML_Epetra_PtAP(*TMT_Matrix_,*NodesToLocalNodes,LocalNodalMatrix,false);
    Teuchos::ParameterList ListN=List_.get("refmaxwell: local nodal list",dummy);
    if (ListN.name() == "refmaxwell: local nodal list") ListN.setName("refmaxwell: 11list");
    SetDefaults("SA",ListN,0,0,false);
    LocalNodalSolver=new MultiLevelPreconditioner(*LocalNodalMatrix,ListN);
  }
  
  /* Output footer */
  if(verbose_ && !Comm_->MyPID()) {
    printf("------------------------------------------------------------------------------\n");
  }

  IsComputePreconditionerOK_=true;
  return 0;
}/*end ComputePreconditioner*/




// ================================================ ====== ==== ==== == = 
// Destroys all structures allocated in \c ComputePreconditioner() if the preconditioner has been computed.
int ML_Epetra::RefMaxwellPreconditioner::DestroyPreconditioner(){
  Operator11_=Teuchos::null;
  if(Diagonal_)  {delete Diagonal_;Diagonal_=0;}

  int printl=ML_Get_PrintLevel();
  int output_level=List_.get("ML output",0);
  output_level=List_.get("output",output_level);

  ML_Set_PrintLevel(output_level);
  if(EdgePC) {delete EdgePC; EdgePC=0;}
  if(NodePC) {delete NodePC; NodePC=0;}
  ML_Set_PrintLevel(printl);  
  if(D0_Matrix_) {delete D0_Matrix_; D0_Matrix_=0;}
  if(TMT_Matrix_) {delete TMT_Matrix_; TMT_Matrix_=0;}
  if(TMT_Agg_Matrix_) {delete TMT_Agg_Matrix_; TMT_Agg_Matrix_=0;}
  if(lump_m1) {delete M1_Matrix_; M1_Matrix_=0;}
  if(BCrows) {delete [] BCrows; BCrows=0;numBCrows=0;}
  ML_Set_PrintLevel(output_level);  
  if(PreEdgeSmoother)  {delete PreEdgeSmoother; PreEdgeSmoother=0;}
  if(PostEdgeSmoother)  {delete PostEdgeSmoother; PostEdgeSmoother=0;}
  if(LocalNodalSolver) {delete LocalNodalSolver; LocalNodalSolver=0;}
  if(LocalNodalMatrix) {delete LocalNodalMatrix; LocalNodalMatrix=0;} 
#ifdef HAVE_ML_IFPACK
  if(IfSmoother) {delete IfSmoother;IfSmoother=0;}
#endif

#ifdef ML_TIMING
  ML_Comm *comm_;
  ML_Comm_Create(&comm_);
  ReportTimer(ConstructionTime_ ,   "ML_RMP::ComputePreconditioner (construction  )",comm_);  
  ReportTimer(FirstApplicationTime_,"ML_RMP::ComputePreconditioner (1st iter time )",comm_);  
  ReportTimer(ApplicationTime_ ,    "ML_RMP::ComputePreconditioner (total itr cost)",comm_);
  ML_Set_PrintLevel(printl);  
  ML_Comm_Destroy(&comm_);
#else
  ML_Set_PrintLevel(printl);  
#endif

  return 0;
}/*end DestroyPreconditioner*/



// ================================================ ====== ==== ==== == = 
// Apply the preconditioner to an Epetra_MultiVector X, puts the result in Y
int ML_Epetra::RefMaxwellPreconditioner::ApplyInverse(const Epetra_MultiVector& B, Epetra_MultiVector& X_) const
{
  int rv;
  /* Sanity Checks */
  if (!B.Map().SameAs(*DomainMap_)) ML_CHK_ERR(-1);
  if (B.NumVectors() != X_.NumVectors()) ML_CHK_ERR(-1);

  /* Check for zero RHS */
  bool norm0=true;
  double *norm=new double[B.NumVectors()]; 
  B.Norm2(norm);
  for(int i=0;norm0==true && i<B.NumVectors();i++) norm0=norm0 && (norm[i]==0);
  delete [] norm;
  if(norm0) return 0;

  /* Build new work vector X */
  Epetra_MultiVector X(X_.Map(),X_.NumVectors());
  X.PutScalar(0);
  
  /* What mode to run in? */
  if(mode=="212") rv=ApplyInverse_Implicit_212(B,X);
  else if(mode=="additive") rv=ApplyInverse_Implicit_Additive(B,X);
  else if(mode=="121") rv=ApplyInverse_Implicit_121(B,X);
  else {fprintf(stderr,"%s","RefMaxwellPreconditioner ERROR: Invalid ApplyInverse mode set in Teuchos list");ML_CHK_ERR(-2);}
  ML_CHK_ERR(rv);

  /* Copy work vector to output */
  X_=X;
  
  /* Timer Stuff */
#ifdef ML_TIMING
  ML_Epetra::RefMaxwellPreconditioner* This = const_cast<ML_Epetra::RefMaxwellPreconditioner *>(this);
  if(FirstApplication_){
    This->FirstApplication_=false;
    This->FirstApplicationTime_=ApplicationTime_;
  }/*end if*/
  This->NumApplications_++;
#endif

  return 0;
}/*end ApplyInverse*/


// ================================================ ====== ==== ==== == = 
int ML_Epetra::RefMaxwellPreconditioner::SetEdgeSmoother(Teuchos::ParameterList &List){  

  string smoother=List.get("smoother: type","Chebyshev");
  int smoother_sweeps=List.get("smoother: sweeps",2);
  int output=List.get("ML output",0);

  if(smoother == "Hiptmair"){
    /* Smoother info output */
    int NumGlobalRows=SM_Matrix_->NumGlobalRows();
    int NumGlobalNNZ=SM_Matrix_->NumGlobalNonzeros();
    if(verbose_ && !Comm_->MyPID()) {
      int edge_sweeps=List.get("subsmoother: edge sweeps",2);
      int node_sweeps=List.get("subsmoother: node sweeps",2);  
      char msg[80];
      sprintf(msg,"%s","RefMaxwell (level 0) :");
      printf("%s # global rows = %d, # estim. global nnz = %d\n",msg,NumGlobalRows,NumGlobalNNZ);     
      printf("%s %s %d (e=%d/n=%d)\n",msg,smoother.c_str(),smoother_sweeps,edge_sweeps,node_sweeps);
    }

    /* Setup Teuchos Lists - Hiptmair */
    int edge_sweeps=List.get("subsmoother: edge sweeps",2);
    int node_sweeps=List.get("subsmoother: node sweeps",2);  

    Teuchos::ParameterList PreList;
    PreList.setName("refmaxwell: edge presmoother");
    PreList.set("coarse: type",smoother);
    PreList.set("coarse: sweeps",1);
    PreList.set("coarse: edge sweeps",smoother_sweeps*edge_sweeps);
    PreList.set("coarse: node sweeps",smoother_sweeps*node_sweeps);
    PreList.set("coarse: subsmoother type",List.get("subsmoother: type","symmetric Gauss-Seidel"));  
    PreList.set("coarse: Chebyshev alpha",List.get("subsmoother: Chebyshev alpha",30.0));
    PreList.set("coarse: MLS alpha",List.get("subsmoother: MLS alpha",30.0));
    PreList.set("coarse: damping factor",List.get("subsmoother: SGS damping factor",1.0));  
    PreList.set("smoother: Hiptmair efficient symmetric",true);
    PreList.set("PDE equations",1);
    PreList.set("max levels",1);
    
    Teuchos::ParameterList PostList(PreList);
    PostList.setName("refmaxwell: edge postsmoother");
    PreList.set("coarse: pre or post","pre");
    PreList.set("smoother: pre or post","pre");
    PreList.set("zero starting solution", true);
    PreList.set("ML label","(1,1) block fine pre-smoother");    
    PreList.set("ML output",output);
    PostList.set("coarse: pre or post","post");
    PostList.set("smoother: pre or post","post"); 
    PostList.set("zero starting solution", false);
    PostList.set("ML label","(1,1) block fine post-smoother");
    PostList.set("ML output",output);
    
    if(HasOnlyDirichletNodes){
      if(PreList.get("coarse: type","dummy") == "Hiptmair"){
        smoother=PreList.get("coarse: subsmoother type","symmetric Gauss-Seidel");
        PreList.set("coarse: type",smoother);
        PostList.set("coarse: type",smoother);
      }/*end if*/
#ifdef HAVE_ML_IFPACK
      IfSmoother=ML_Gen_Smoother_Ifpack_Epetra(const_cast<Epetra_CrsMatrix*>(&*SM_Matrix_),0,List_,"RefMaxwell (level 0): ",verbose_);
#else
      PreEdgeSmoother  = new MultiLevelPreconditioner(*SM_Matrix_,PreList);
      PostEdgeSmoother = new MultiLevelPreconditioner(*SM_Matrix_,PostList);    
#endif
    }/*end if*/
    else{
      PreEdgeSmoother  = new MultiLevelPreconditioner(*SM_Matrix_,*D0_Matrix_,*TMT_Matrix_,PreList,true,true);
      PostEdgeSmoother = new MultiLevelPreconditioner(*SM_Matrix_,*D0_Matrix_,*TMT_Matrix_,PostList,true,true);
    }/*end if*/
  }/*end if*/
#ifdef HAVE_ML_IFPACK
  else {
    IfSmoother=ML_Gen_Smoother_Ifpack_Epetra(const_cast<Epetra_CrsMatrix*>(&*SM_Matrix_),0,List_,"RefMaxwell (level 0): ",verbose_);
  }/*end else*/
#endif
         
  return 0;
}/*end SetEdgeSmoother*/

// ================================================ ====== ==== ==== == = 
//! Implicitly applies in the inverse in a 2-1-2 format
int ML_Epetra::RefMaxwellPreconditioner::ApplyInverse_Implicit_212(const Epetra_MultiVector& B, Epetra_MultiVector& X) const
{

#ifdef ML_TIMING
  double t_time,t_diff;
  StartTimer(&t_time);
#endif
   
  int NumVectors=B.NumVectors();

  /* Setup Temps */  
  Epetra_MultiVector node_sol1(*NodeMap_,NumVectors,true);
  Epetra_MultiVector node_sol2(*NodeMap_,NumVectors,false);
  Epetra_MultiVector node_rhs(*NodeMap_,NumVectors,false);
  Epetra_MultiVector edge_temp1(*DomainMap_,NumVectors,false);
  Epetra_MultiVector edge_rhs(*DomainMap_,NumVectors,false);
  Epetra_MultiVector edge_sol(*DomainMap_,NumVectors,true);


  /* Build Nodal RHS */
  ML_CHK_ERR(D0_Matrix_->Multiply(true,B,node_rhs));

  /* Precondition (2,2) Block */
  ML_CHK_ERR(NodePC->ApplyInverse(node_rhs,node_sol1));
  
  /* Build Residual */
  ML_CHK_ERR(D0_Matrix_->Multiply(false,node_sol1,edge_temp1));
  ML_CHK_ERR(edge_rhs.Update(1.0,B,-1.0));

  /* Precondition (1,1) Block */
  //  _CHK_ERR(PreEdgeSmoother->ApplyInverse(B,X));
  ML_CHK_ERR(EdgePC->ApplyInverse(edge_rhs,edge_sol));

  /* Build Nodal RHS */  
  ML_CHK_ERR(edge_temp1.Update(1.0,edge_rhs,-1.0));
  ML_CHK_ERR(D0_Matrix_->Multiply(true,edge_temp1,node_rhs));

  /* Precondition (2,2) Block */  
  ML_CHK_ERR(NodePC->ApplyInverse(node_rhs,node_sol2));
  
  /* Assemble solution (x = xe + T*(xn1 + xn2)) */
  ML_CHK_ERR(node_sol1.Update(1.0,node_sol2,1.0));
  
  ML_CHK_ERR(D0_Matrix_->Multiply(false,node_sol1,X));
  ML_CHK_ERR(X.Update(1.0,edge_sol,1.0));

#ifdef ML_TIMING
  StopTimer(&t_time,&t_diff);
  /* Output */
  ML_Comm *comm_;
  ML_Comm_Create(&comm_);
  this->ApplicationTime_+= t_diff;  
  ML_Comm_Destroy(&comm_);
#endif  

  
  return 0;
}/*end ApplyInverse_Implicit_212*/

// ================================================ ====== ==== ==== == = 
//! Implicitly applies in the inverse in an additive format
int  ML_Epetra::RefMaxwellPreconditioner::ApplyInverse_Implicit_Additive(const Epetra_MultiVector& B, Epetra_MultiVector& X) const
{
#ifdef ML_TIMING
  double t_time,t_diff;
  StartTimer(&t_time);
#endif

  int NumVectors=B.NumVectors();
  Epetra_MultiVector TempE1(X.Map(),NumVectors,false);
  Epetra_MultiVector TempE2(X.Map(),NumVectors,true);
  Epetra_MultiVector TempN1(*NodeMap_,NumVectors,false);
  Epetra_MultiVector TempN2(*NodeMap_,NumVectors,true);
  Epetra_MultiVector Resid(B.Map(),NumVectors);

  /* Pre-Smoothing */
#ifdef HAVE_ML_IFPACK
  if(IfSmoother) {ML_CHK_ERR(IfSmoother->ApplyInverse(B,X));}
  else 
#endif
  if(PreEdgeSmoother) ML_CHK_ERR(PreEdgeSmoother->ApplyInverse(B,X));

  /* Build Residual */
  ML_CHK_ERR(SM_Matrix_->Multiply(false,X,TempE1));
  ML_CHK_ERR(Resid.Update(-1.0,TempE1,1.0,B,0.0));  

  if(!HasOnlyDirichletNodes){
    ML_CHK_ERR(D0_Matrix_->Multiply(true,Resid,TempN1));
  }

  /* Precondition (1,1) block (additive)*/
  ML_CHK_ERR(EdgePC->ApplyInverse(Resid,TempE2));

  /* Precondition (2,2) block (additive)*/
  if(!HasOnlyDirichletNodes){
    ML_CHK_ERR(NodePC->ApplyInverse(TempN1,TempN2));             

    /* EXPERIMENTAL: Local Nodal Stuff, if active */
    if(use_local_nodal_solver){
      const Epetra_Map& LocalMap=LocalNodalMatrix->DomainMap();
      Epetra_MultiVector TempNL1(LocalMap,NumVectors,true);
      Epetra_MultiVector TempNL2(LocalMap,NumVectors,true);
      Epetra_MultiVector TempN3(*NodeMap_,NumVectors,true);
      
      NodesToLocalNodes->Multiply(true,TempN1,TempNL1);
      LocalNodalSolver->ApplyInverse(TempNL1,TempNL2);
      NodesToLocalNodes->Multiply(false,TempNL2,TempN3);
      TempN2.Update(1.0,TempN3,1.0);
    }/*end if*/

    D0_Matrix_->Multiply(false,TempN2,TempE1);
  }/*end if*/

  /* Update solution */
  if(HasOnlyDirichletNodes) X.Update(1.0,TempE2,1.0);
  else X.Update(1.0,TempE1,1.0,TempE2,1.0);

  /* Post-Smoothing */
#ifdef HAVE_ML_IFPACK
  if(IfSmoother) {ML_CHK_ERR(IfSmoother->ApplyInverse(B,X));}
  else 
#endif
    if(PostEdgeSmoother) ML_CHK_ERR(PostEdgeSmoother->ApplyInverse(B,X));

  
#ifdef ML_TIMING
  StopTimer(&t_time,&t_diff);
  /* Output */
  ML_Comm *comm_;
  ML_Comm_Create(&comm_);
  this->ApplicationTime_+= t_diff;  
  ML_Comm_Destroy(&comm_);
#endif  
  
  return 0;
}

// ================================================ ====== ==== ==== == = 
//! Implicitly applies in the inverse in an 1-2-1 format
int  ML_Epetra::RefMaxwellPreconditioner::ApplyInverse_Implicit_121(const Epetra_MultiVector& B, Epetra_MultiVector& X) const
{
#ifdef ML_TIMING
  double t_time,t_diff;
  StartTimer(&t_time);
#endif
  
  int NumVectors=B.NumVectors();
  Epetra_MultiVector TempE1(X.Map(),NumVectors,false);
  Epetra_MultiVector TempE2(X.Map(),NumVectors,true);
  Epetra_MultiVector TempN1(*NodeMap_,NumVectors,false);
  Epetra_MultiVector TempN2(*NodeMap_,NumVectors,true);
  Epetra_MultiVector Resid(B);
  

  /* Pre-Smoothing */
  ML_CHK_ERR(PreEdgeSmoother->ApplyInverse(B,X));
  
  /* Precondition (1,1) Block */
  ML_CHK_ERR(EdgePC->ApplyInverse(Resid,TempE2));
  ML_CHK_ERR(X.Update(1.0,TempE2,1.0));;  
 
  /* Build Residual */
  ML_CHK_ERR(SM_Matrix_->Multiply(false,X,TempE1));
  ML_CHK_ERR(Resid.Update(-1.0,TempE1,1.0,B,0.0));  
  if(!HasOnlyDirichletNodes){
    ML_CHK_ERR(D0_Matrix_->Multiply(true,Resid,TempN1));
  }
  
  /* Precondition (2,2) Block */
  if(!HasOnlyDirichletNodes){
    ML_CHK_ERR(NodePC->ApplyInverse(TempN1,TempN2));             
    D0_Matrix_->Multiply(false,TempN2,TempE1);
  }/*end if*/
  if(!HasOnlyDirichletNodes) X.Update(1.0,TempE1,1.0);  

  /* Build Residual */
  ML_CHK_ERR(SM_Matrix_->Multiply(false,X,TempE1));
  ML_CHK_ERR(Resid.Update(-1.0,TempE1,1.0,B,0.0));  
  
  /* Precondition (1,1) Block */
  TempE2.PutScalar(0.0);
  ML_CHK_ERR(EdgePC->ApplyInverse(Resid,TempE2));
  ML_CHK_ERR(X.Update(1.0,TempE2,1.0));;  

  /* Post-Smoothing */
  ML_CHK_ERR(PostEdgeSmoother->ApplyInverse(B,X));
  
#ifdef ML_TIMING
  StopTimer(&t_time,&t_diff);
  /* Output */
  ML_Comm *comm_;
  ML_Comm_Create(&comm_);
  this->ApplicationTime_+= t_diff;  
  ML_Comm_Destroy(&comm_);
#endif  
  
  return 0;
}

// ================================================ ====== ==== ==== == = 
int ML_Epetra::SetDefaultsRefMaxwell(Teuchos::ParameterList & inList,bool OverWrite)
{  
  /* Sublists */
  Teuchos::ParameterList ListRF,List11, List11c, List22, dummy;
  Teuchos::ParameterList & List11_=inList.sublist("refmaxwell: 11list");
  Teuchos::ParameterList & List22_=inList.sublist("refmaxwell: 22list");
  Teuchos::ParameterList & List11c_=List11_.sublist("edge matrix free: coarse");

  /* Build Teuchos List: (1,1) coarse */    
  ML_Epetra::SetDefaults("SA",List11c);
  List11c.set("cycle applications",1);
  List11c.set("smoother: type","Chebyshev");
  List11c.set("aggregation: threshold",.01);
  List11c.set("coarse: type","Amesos-KLU");  
  List11c.set("ML label","coarse (1,1) block");
  ML_Epetra::UpdateList(List11c,List11c_,OverWrite);
  
  /* Build Teuchos List: (1,1) */
  ML_Epetra::SetDefaults("SA",List11);
  List11.set("cycle applications",1);
  List11.set("aggregation: type","Uncoupled");
  List11.set("smoother: sweeps",0);
  List11.set("edge matrix free: coarse",List11c);
  List11.set("aggregation: threshold",.01);  
  ML_Epetra::UpdateList(List11,List11_,OverWrite);
  
  /* Build Teuchos List: (2,2) */  
  ML_Epetra::SetDefaults("SA",List22);  
  List22.set("cycle applications",1);
  List22.set("smoother: type","Chebyshev");
  List22.set("aggregation: type","Uncoupled");
  List22.set("aggregation: threshold",.01);
  List22.set("coarse: type","Amesos-KLU");
  List22.set("ML label","(2,2) block");

  // This line is commented out due to IFPACK issues
  //  List22.set("smoother: sweeps (level 0)",0);
  
  ML_Epetra::UpdateList(List22,List22_,OverWrite);    
  
  /* Build Teuchos List: Overall */  
  SetDefaults("maxwell",ListRF,0,0,false);
  ListRF.set("smoother: type","Chebyshev");
  ListRF.set("smoother: sweeps",2);  
  ListRF.set("refmaxwell: 11solver","edge matrix free");
  ListRF.set("refmaxwell: 11list",List11);
  ListRF.set("refmaxwell: 22solver","multilevel");
  ListRF.set("refmaxwell: 22list",List22);
  ListRF.set("refmaxwell: mode","additive");
  ListRF.set("default values","RefMaxwell");
  ListRF.set("zero starting solution",false);

  ML_Epetra::UpdateList(ListRF,inList,OverWrite);
  
  return 0;  
}/*end SetDefaultsRefMaxwell*/

#endif
