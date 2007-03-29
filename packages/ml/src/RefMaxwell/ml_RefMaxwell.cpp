#include "ml_config.h"
#include <string.h>
#include "ml_RefMaxwell.h"
#include "ml_epetra.h"
#include "ml_epetra_utils.h"
#include "ml_MultiLevelPreconditioner.h"
#include "ml_RefMaxwell_11_Operator.h"
#include "ml_EdgeMatrixFreePreconditioner.h"

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_EPETRAEXT) && defined(HAVE_ML_IFPACK)
#include "EpetraExt_MatrixMatrix.h" //haq
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_CrsMatrixIn.h"//haq
//#include "Ifpack_Hiptmair.h"

static int c_iteration=0;//DEBUG

//#define NO_OUTPUT


extern "C"{
  void cms_vec_dump3(double *x,int N,char* fname, int iteration);
  void cms_vec_dump2(double *x,int N,char* fname){    
#ifndef NO_OUTPUT
    cms_vec_dump3(x,N,fname,c_iteration);
#endif
  }
}

/* Stuff to avoid errors wrt rpc.cpp */
#ifdef NO_OUTPUT
#define MVOUT2(x,y,z) ;
#define MVOUT(x,y) ;
#define Epetra_CrsMatrix_Print(x,y) ;
#else
void cms_vec_dump3(double *x,int N,char* fname, int iteration){
  char fn[80];
  sprintf(fn,"%s.%d.dat",fname,iteration);
  FILE *f=fopen(fn,"w");
  for(int i=0;i<N;i++)
    fprintf(f,"%22.16e\n",x[i]);
  fclose(f);
}

void MVOUT(const Epetra_MultiVector & A, char *of){
  ofstream os(of);
  int i,j;
  int NumProc=A.Map().Comm().NumProc();
  int MyPID  =A.Map().Comm().MyPID();
  int NumVectors=A.NumVectors();
  
  for (int iproc=0; iproc < NumProc; iproc++) {
    if (MyPID==iproc) {
      int MyLength = A.MyLength();
      for (i=0; i<MyLength; i++) {        
	for (j = 0; j < NumVectors ; j++) {
          os.width(20);
          os.precision(16);
          os.setf(ios_base::scientific,ios_base::floatfield);
          os << A[j][i];
          os << "   ";
        }
        os<<endl;
      }/*end for*/
      os << flush;      
    }/*end if*/
    // Do a few global ops to give I/O a chance to complete
    A.Map().Comm().Barrier();
    A.Map().Comm().Barrier();
    A.Map().Comm().Barrier();
  }/*end for*/
}/*end MultiVectorToMatlabFile*/


void MVOUT2(const Epetra_MultiVector & A,char* pref,int idx){
  char c[80];
  sprintf(c,"%s.%d.dat",pref,idx);
  MVOUT(A,c);
}/* end MVOUT2*/

void Epetra_CrsMatrix_Print(const Epetra_CrsMatrix& A, char* of) {
  EpetraExt::RowMatrixToMatlabFile(of,A);      
}/*end Epetra_CrsMatrix_Print*/

void IVOUT(const Epetra_IntVector & A, char *of){
  ofstream os(of);
  int i;
  int NumProc=A.Map().Comm().NumProc();
  int MyPID  =A.Map().Comm().MyPID();
  
  for (int iproc=0; iproc < NumProc; iproc++) {
    if (MyPID==iproc) {
      int MyLength = A.MyLength();
      for (i=0; i<MyLength; i++) {        
          os.width(20);
          os << A[i]<<endl;
      }
      os << flush;      
    }/*end if*/
    // Do a few global ops to give I/O a chance to complete
    A.Map().Comm().Barrier();
    A.Map().Comm().Barrier();
    A.Map().Comm().Barrier();
  }/*end for*/
}/*end MultiVectorToMatlabFile*/

//#else
//extern void MVOUT (const Epetra_MultiVector & A, ostream & os);//HAQ
//extern void MVOUT2(const Epetra_MultiVector & A,char* pref,int idx);//HAQ
//extern void Epetra_CrsMatrix_Print(const Epetra_CrsMatrix& A, ostream& os);//HAQ
#endif



// ================================================ ====== ==== ==== == = 
ML_Epetra::RefMaxwellPreconditioner::RefMaxwellPreconditioner(const Epetra_CrsMatrix& SM_Matrix,      //S+M
                                                              const Epetra_CrsMatrix& D0_Clean_Matrix,//T or D0 w/ nothing zero'd
                                                              const Epetra_CrsMatrix& Ms_Matrix,      //M1(sigma)
                                                              const Epetra_CrsMatrix& M0inv_Matrix,   //M0^{-1}
                                                              const Epetra_CrsMatrix& M1_Matrix,      //M1(1)
                                                              //                                                              const Epetra_CrsMatrix& TMT_Matrix,     //T' M1(sigma) T
                                                              const Teuchos::ParameterList& List,
                                                              const bool ComputePrec):
  ML_Preconditioner(),SM_Matrix_(&SM_Matrix),D0_Matrix_(0), D0_Clean_Matrix_(&D0_Clean_Matrix),Ms_Matrix_(&Ms_Matrix),
  M0inv_Matrix_(&M0inv_Matrix),M1_Matrix_(&M1_Matrix),TMT_Matrix_(0),TMT_Agg_Matrix_(0),
  Diagonal_(0),Operator11_(0),BCrows(0),numBCrows(0),EdgePC(0),NodePC(0),HasOnlyDirichletNodes(false)
{
  /* Set the Epetra Goodies */
  Comm_ = &(SM_Matrix_->Comm());
  DomainMap_ = &(SM_Matrix_->OperatorDomainMap());
  RangeMap_ = &(SM_Matrix_->OperatorRangeMap());
  NodeMap_ = &(D0_Clean_Matrix_->OperatorDomainMap());

  Label_=strdup("ML reformulated Maxwell preconditioner");
  List_=List;
  
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
void ML_Epetra::RefMaxwellPreconditioner::Print(const char *whichHierarchy){
  if(IsComputePreconditionerOK_ && EdgePC && !strcmp(whichHierarchy,"11")) EdgePC->Print("main");
  if(IsComputePreconditionerOK_ && NodePC && !strcmp(whichHierarchy,"22")) NodePC->Print("main");  
}/*end Print*/


// ================================================ ====== ==== ==== == = 
// Computes the preconditioner
int ML_Epetra::RefMaxwellPreconditioner::ComputePreconditioner(const bool CheckFiltering)
{
  Teuchos::ParameterList dummy;

  /* Pull Solver Mode */
  mode=List_.get("refmaxwell: mode","212");
  
  /* Nuke everything if we've done this already */
  if(IsComputePreconditionerOK_) DestroyPreconditioner();

  /* Find the Dirichlet Rows (using SM_Matrix_) and columns (using D0_Clean_Matrix_) */
  BCrows=FindLocalDiricheltRowsFromOnesAndZeros(*SM_Matrix_,numBCrows);
  Epetra_IntVector * BCnodes=FindLocalDirichletColumnsFromRows(BCrows,numBCrows,*D0_Clean_Matrix_);
   
#ifdef STUFF_WE_PROBABLY_DONT_NEED
  /* Convert to Dirichlet node list */
  int Nn=BCnodes->MyLength();
  int *BCnodes_int= new int[Nn];
  int numBCnodes=0;
  for(int i=0;i<Nn;i++)
    if((*BCnodes)[i]) {
      int cid=D0_Clean_Matrix_->GCID(i);
      int lid=NodeMap_->LID(cid);
      if(cid==0) fprintf(stderr,"[%d] Error gcid %d map fails\n",Comm_->MyPID(),i);
      else if(lid==-1) fprintf(stderr,"[%d] Error nlid %d/%d map fails\n",Comm_->MyPID(),i,cid);      
      BCnodes_int[numBCnodes++]=lid;
    }
  for(int i=0;i<numBCnodes;i++){
    if(BCnodes_int[i] == -1) fprintf(stderr,"[%d] Error row %d/%d map fails\n",Comm_->MyPID(),i,BCnodes_int[i]);
  }
#else
  int Nn=BCnodes->MyLength();
  int numBCnodes=0;
  for(int i=0;i<Nn;i++) if((*BCnodes)[i]) numBCnodes++;
#endif


  /* Sanity Check: We have at least some Dirichlet nodes */
  /* NTS: This should get fixed later for robustness */
  int HasInterior = numBCnodes != Nn;
  printf("HasInterior = %d %d/%d\n",HasInterior,Nn-numBCnodes,Nn);
  int globalInterior=HasInterior;
#ifdef DOESNT_WORK_DONT_KNOW_WHY
  Comm_->MaxAll(&HasInterior,&globalInterior,1);  
  //  if(!globalInterior){
#endif
  if(!globalInterior){
    HasOnlyDirichletNodes=true;
    if(!Comm_->MyPID()) printf("WARNING: All nodes are Dirichlet nodes!  Refmax may crash horribly.  YMMV,VWP,USR18OO.\n");
  }/*end if*/

  D0_Matrix_ = new Epetra_CrsMatrix(*D0_Clean_Matrix_);
  /* Do the Nuking for D0_Matrix_ */
  Apply_BCsToMatrixRows(BCrows,numBCrows,*D0_Matrix_);
  Apply_BCsToMatrixColumns(*BCnodes,*D0_Matrix_);   

  /* Build the TMT Matrix */  
  /* NTS: When ALEGRA builds this matrix itself, we get rid of these lines */
  if(!HasOnlyDirichletNodes){
    Epetra_CrsMatrix temp1(Copy,*RangeMap_,1);
    TMT_Matrix_=new Epetra_CrsMatrix(Copy,*NodeMap_,1);
    EpetraExt::MatrixMatrix::Multiply(*Ms_Matrix_,false,*D0_Matrix_,false,temp1);
    EpetraExt::MatrixMatrix::Multiply(*D0_Matrix_,true,temp1,false,*TMT_Matrix_);
    Remove_Zeroed_Rows(*TMT_Matrix_);
    TMT_Matrix_->OptimizeStorage();
  }/*end if*/

  /* Build the TMT-Agg Matrix  (used for aggregating the (1,1) block*/
  /* NTS: When ALEGRA builds this matrix itself, we get rid of these lines */
  Epetra_CrsMatrix temp2(Copy,*RangeMap_,0);
  TMT_Agg_Matrix_=new Epetra_CrsMatrix(Copy,*NodeMap_,0);  
  EpetraExt::MatrixMatrix::Multiply(*M1_Matrix_,false,*D0_Clean_Matrix_,false,temp2);
  EpetraExt::MatrixMatrix::Multiply(*D0_Clean_Matrix_,true,temp2,false,*TMT_Agg_Matrix_);
  //EpetraExt::MatrixMatrix::Multiply(*M1_Matrix_,false,*D0_Matrix_,false,temp2);
  //  EpetraExt::MatrixMatrix::Multiply(*D0_Matrix_,true,temp2,false,*TMT_Agg_Matrix_);
  Remove_Zeroed_Rows(*TMT_Agg_Matrix_);


  TMT_Agg_Matrix_->OptimizeStorage();

  // NTS: Should I be building w/ D0_Clean for aggregation?  
  
#ifndef NO_OUTPUT
  if(TMT_Matrix_) Epetra_CrsMatrix_Print(*TMT_Matrix_,"tmt_matrix.dat");
  Epetra_CrsMatrix_Print(*TMT_Agg_Matrix_,"tmt_agg_matrix.dat");
#endif

  
  // NTS: SHOULD M0 Be nuked in here???

  /* Boundary nuke the edge matrices */
  Apply_OAZToMatrix(BCrows,numBCrows,*Ms_Matrix_);
  Apply_OAZToMatrix(BCrows,numBCrows,*M1_Matrix_);    

  /* DEBUG: Output matrices */
#ifndef NO_OUTPUT
  Epetra_CrsMatrix_Print(*SM_Matrix_,"sm_matrix.dat");
  Epetra_CrsMatrix_Print(*M1_Matrix_,"m1_nuked.dat");  
  Epetra_CrsMatrix_Print(*D0_Matrix_,"d9_nuked.dat");  
  Epetra_CrsMatrix_Print(*D0_Clean_Matrix_,"d0_clean.dat");  
#endif
  
  /* Cleanup from the Boundary Conditions */
  delete BCnodes; 
#ifdef STUFF_WE_PROBABLY_DONT_NEED
  delete [] BCnodes_int;
#endif
  
  /* Build the (1,1) Block Operator */
  Operator11_ = new ML_RefMaxwell_11_Operator(*SM_Matrix_,*D0_Matrix_,*M0inv_Matrix_,*M1_Matrix_);
  
  if(mode!="additive"){  
    /* Approximate the Diagonal of the (1,1) Block Operator */
    /*  M1L=spdiags(M1*ones(Ne,1),0,Ne,Ne);
        M1LT=M1L*T;
        DIAG11=diag(SM)+diag(M1LT*(M0L\M1LT')); */
    Epetra_Vector temp_edge1(*DomainMap_,false);
    Epetra_Vector temp_edge2(*DomainMap_,false);
    Diagonal_=new Epetra_Vector(*DomainMap_,false);
    temp_edge1.PutScalar(1.0);
    
    
    /* Build a lumped approximation of M1 */
    M1_Matrix_->Multiply(false,temp_edge1,temp_edge2);
    Epetra_CrsMatrix M1diag(Copy,*DomainMap_,1,true);
    for(int i=0;i<M1diag.NumMyRows();i++){
      int gid=DomainMap_->GID(i);
      M1diag.InsertGlobalValues(gid,1,&(temp_edge2[i]),&gid);
    }/*end for*/
    M1diag.FillComplete();
    M1diag.OptimizeStorage();
    
    //  ofstream ofs0("m1l.dat");
    //  MVOUT(temp_edge2,ofs0);
    
    
    /* Build lumped-M1 approximation of the whole matrix and pull the diagonal */
    Epetra_CrsMatrix temp_matrix1(Copy,*NodeMap_,0);
    Epetra_CrsMatrix temp_matrix2(Copy,*NodeMap_,0);
    Epetra_CrsMatrix temp_matrix3(Copy,*DomainMap_,0);
    EpetraExt::MatrixMatrix::Multiply(*D0_Matrix_,true,M1diag,false,temp_matrix1);
    EpetraExt::MatrixMatrix::Multiply(*M0inv_Matrix_,false,temp_matrix1,false,temp_matrix2);
    EpetraExt::MatrixMatrix::Multiply(temp_matrix1,true,temp_matrix2,false,temp_matrix3);
    temp_matrix3.ExtractDiagonalCopy(*Diagonal_);
    
#ifndef NO_OUTPUT
    MVOUT(*Diagonal_,"addon_diagonal.dat");
#endif
 
    SM_Matrix_->ExtractDiagonalCopy(temp_edge2);
    Diagonal_->Update(1.0,temp_edge2,1.0);
    
#ifndef NO_OUTPUT
    MVOUT(*Diagonal_,"diagonal.dat");
#endif
  
  /* Sanity Check the Diagonal */
    double min_val; Diagonal_->MinValue(&min_val);
    if(Comm_->MyPID()==0 && min_val <1e-16) {printf("ERROR: Minimum Estimated Diagonal <1e-16 (%6.4e)\n",min_val);return -2;}
  }/*end if*/
  else{
    // THIS IS A HACK!!!!!
    Diagonal_=new Epetra_Vector(*DomainMap_,false);
    Diagonal_->PutScalar(1.0);
  }

    
  /* Build the (2,2) Block Preconditioner */
  if(!HasOnlyDirichletNodes){
    ML_reseed_random_vec(8675309);//DEBUG
    string solver22=List_.get("refmaxwell: 22solver","multilevel");
    Teuchos::ParameterList List22=List_.get("refmaxwell: 22list",dummy);
    SetDefaultsSA(List22,0,0,false);
    if(solver22=="multilevel") NodePC=new MultiLevelPreconditioner(*TMT_Matrix_,List22);
    else {printf("RefMaxwellPreconditioner: ERROR - Illegal (2,2) block preconditioner\n");return -1;}
    //NTS: Add Adaptive, MatrixFree
  }/*end if*/

  const ML_Aggregate* nodal_aggregates=NULL;
  if(!HasOnlyDirichletNodes) nodal_aggregates=NodePC->GetML_Aggregate();
  // NTS: Delete me
    
  /* Build the (1,1) Block Preconditioner */ 
  ML_reseed_random_vec(8675309);//DEBUG 
  string solver11=List_.get("refmaxwell: 11solver","edge matrix free");
  Teuchos::ParameterList List11=List_.get("refmaxwell: 11list",dummy);
  //  if(solver11=="edge matrix free") EdgePC=new EdgeMatrixFreePreconditioner(*Operator11_,*Diagonal_,*D0_Matrix_,*D0_Clean_Matrix_,*TMT_Matrix_,List11,true);
  if(solver11=="edge matrix free")
    //    EdgePC=new EdgeMatrixFreePreconditioner(*Operator11_,*Diagonal_,*D0_Matrix_,*D0_Clean_Matrix_,*TMT_Matrix_,nodal_aggregates,BCrows,numBCrows,List11,true);
    EdgePC=new EdgeMatrixFreePreconditioner(*Operator11_,*Diagonal_,*D0_Matrix_,*D0_Clean_Matrix_,*TMT_Agg_Matrix_,nodal_aggregates,BCrows,numBCrows,List11,true);
    else {printf("RefMaxwellPreconditioner: ERROR - Illegal (1,1) block preconditioner\n");return -1;}


  /* Setup the Hiptmair smoother in additive mode */
  if(mode=="additive"){
    SetEdgeSmoother(List_);
  }/*end if*/
    
  IsComputePreconditionerOK_=true;
  return 0;
}/*end ComputePreconditioner*/




// ================================================ ====== ==== ==== == = 
// Destroys all structures allocated in \c ComputePreconditioner() if the preconditioner has been computed.
int ML_Epetra::RefMaxwellPreconditioner::DestroyPreconditioner(){
  if(Operator11_) {delete Operator11_;Operator11_=0;}
  if(Diagonal_)  {delete Diagonal_;Diagonal_=0;}
  if(EdgePC) {delete EdgePC; EdgePC=0;}
  if(NodePC) {delete NodePC; NodePC=0;}
  if(D0_Matrix_) {delete D0_Matrix_; D0_Matrix_=0;}
  if(TMT_Matrix_) {delete TMT_Matrix_; TMT_Matrix_=0;}
  if(TMT_Agg_Matrix_) {delete TMT_Agg_Matrix_; TMT_Agg_Matrix_=0;}
  if(BCrows) {delete [] BCrows; BCrows=0;numBCrows=0;}
  if(PreEdgeSmoother)  {delete PreEdgeSmoother; PreEdgeSmoother=0;}
  if(PostEdgeSmoother) {delete PostEdgeSmoother; PostEdgeSmoother=0;}
  return 0;
}/*end DestroyPreconditioner*/



// ================================================ ====== ==== ==== == = 
// Apply the preconditioner to an Epetra_MultiVector X, puts the result in Y
int ML_Epetra::RefMaxwellPreconditioner::ApplyInverse(const Epetra_MultiVector& B, Epetra_MultiVector& X) const
{
  int rv;
  /* Sanity Checks */
  if (!B.Map().SameAs(*DomainMap_)) ML_CHK_ERR(-1);
  if (B.NumVectors() != X.NumVectors()) ML_CHK_ERR(-1);

  /* Check for zero RHS */
  bool norm0=true;
  double *norm=new double[B.NumVectors()]; 
  B.Norm2(norm);
  for(int i=0;norm0==true && i<B.NumVectors();i++) norm0=norm0 && (norm[i]==0);
  delete [] norm;
  if(norm0) return 0;

  /* Build new work vector if X and B are the same */
  //  Epetra_MultiVector *X=&X_;
  //  if(&X_==&B)
  //    X=new Epetra_MultiVector(X_);
  //X=new Epetra_MultiVector(X_.Map(),X_.NumVectors(),true);

  
  /* What mode to run in? */
  if(mode=="212") rv=ApplyInverse_Implicit_212(B,X);
  else if(mode=="additive") rv=ApplyInverse_Implicit_Additive(B,X);
  else if(mode=="121") rv=ApplyInverse_Implicit_121(B,X);
  else {fprintf(stderr,"RefMaxwellPreconditioner ERROR: Invalid ApplyInverse mode set in Teuchos list");ML_CHK_ERR(-2);}
  ML_CHK_ERR(rv);

  /* Copy work vector to output if needed */
  //  if(&X_==&B) {
  //    X_=*X;
  //    delete X;
  //  }/*end if*/
  
  return 0;
}/*end ApplyInverse*/


// ================================================ ====== ==== ==== == = 
int ML_Epetra::RefMaxwellPreconditioner::SetEdgeSmoother(Teuchos::ParameterList &List1){  
  Teuchos::ParameterList dummy;
  Teuchos::ParameterList List = List1.get("refmaxwell: additive smoother",dummy);
  cout<<"************** SetEdgeSmoother **************"<<endl;

  /* Setup Teuchos Lists*/
  Teuchos::ParameterList PreList(List);
  PreList.set("PDE equations",1);
  PreList.set("max levels",1);
  PreList.set("coarse: pre or post","pre");
  PreList.set("smoother: pre or post","pre");
  PreList.set("zero starting solution", true);
  Teuchos::ParameterList PostList(List);
  PostList.set("PDE equations",1);
  PostList.set("max levels",1); 
  PostList.set("coarse: pre or post","post");
  PostList.set("smoother: pre or post","post");
  PostList.set("zero starting solution", false);// this should be FALSE!!!
  
  //NTS: TMT_Matrix_ isn't quite the same as the Node_Matrix_ we want.
  //NTS: This is a total kludge, since Hiptmair will build its own copies of
  // this anyway.

  if(HasOnlyDirichletNodes){
    if(List.get("coarse: type","dummy") == "Hiptmair"){
      string smoother;
      smoother=List.get("coarse: subsmoother type","symmetric Gauss-Seidel");
      PreList.set("coarse: type",smoother);
      PostList.set("coarse: type",smoother);
    }/*end if*/      
    PreEdgeSmoother  = new MultiLevelPreconditioner(*SM_Matrix_,PreList);
    PostEdgeSmoother = new MultiLevelPreconditioner(*SM_Matrix_,PostList);    
  }/*end if*/
  else{
    PreEdgeSmoother  = new MultiLevelPreconditioner(*SM_Matrix_,*D0_Matrix_,*TMT_Matrix_,PreList);
    PostEdgeSmoother = new MultiLevelPreconditioner(*SM_Matrix_,*D0_Matrix_,*TMT_Matrix_,PostList);
  }/*end if*/
    
  cout<<"*********************************************"<<endl;  
  return 0;
}/*end SetEdgeSmoother*/


// ================================================ ====== ==== ==== == = 

void cms_residual_check(const char * tag, const Epetra_Operator * op,const Epetra_MultiVector& rhs, const Epetra_MultiVector& lhs){
  int NumVectors=rhs.NumVectors();
  double *norm_old, *norm_new;
  norm_old=new double[NumVectors];
  norm_new=new double[NumVectors];

  Epetra_MultiVector temp(rhs);
  op->Apply(lhs,temp);
  temp.Update(1.0,rhs,-1.0);
  
  rhs.Norm2(norm_old);
  temp.Norm2(norm_new);
  if(op->Comm().MyPID()==0)
    for(int i=0;i<NumVectors;i++)
      printf("%s[%d]: Norm Reduction %6.4e [%6.4e]\n",tag,i,norm_new[i] / norm_old[i],norm_old[i]);  

  delete [] norm_old; delete [] norm_new;
}

double cms_compute_residual(const Epetra_Operator * op,const Epetra_MultiVector& rhs, const Epetra_MultiVector& lhs){
  int NumVectors=rhs.NumVectors();
  double *norm_old, *norm_new;
  norm_old=new double[NumVectors];
  norm_new=new double[NumVectors];

  Epetra_MultiVector temp(rhs);
  op->Apply(lhs,temp);
  temp.Update(1.0,rhs,-1.0);  
  temp.Norm2(norm_new);
  rhs.Norm2(norm_old);  
  double rv=norm_new[0] / norm_old[0];


  delete [] norm_old; delete [] norm_new;
  return rv;
}


// ================================================ ====== ==== ==== == = 
//! Implicitly applies in the inverse in a 2-1-2 format
int ML_Epetra::RefMaxwellPreconditioner::ApplyInverse_Implicit_212(const Epetra_MultiVector& B, Epetra_MultiVector& X) const
{
  int NumVectors=B.NumVectors();
  double r0,r1,r2,r3,r4;

  double norm;
  
  MVOUT2(B,"b",c_iteration);//DEBUG

  r0=cms_compute_residual(SM_Matrix_,B,X);//DEBUG  
  
  /* Setup Temps */  
  Epetra_MultiVector node_sol1(*NodeMap_,NumVectors,true);
  Epetra_MultiVector node_sol2(*NodeMap_,NumVectors,false);
  Epetra_MultiVector node_rhs(*NodeMap_,NumVectors,false);
  Epetra_MultiVector edge_temp1(*DomainMap_,NumVectors,false);
  Epetra_MultiVector edge_rhs(*DomainMap_,NumVectors,false);
  Epetra_MultiVector edge_sol(*DomainMap_,NumVectors,true);

  /* Build Nodal RHS  (nrhs = D0'*b) */
  ML_CHK_ERR(D0_Matrix_->Multiply(true,B,node_rhs));


  MVOUT2(node_rhs,"nrhs1",c_iteration);//DEBUG
  
  /* (2,2) Block Solve (xn1 = TMT^{-1} nrhs) */
  ML_reseed_random_vec(8675309);  
  ML_CHK_ERR(NodePC->ApplyInverse(node_rhs,node_sol1));
  MVOUT2(node_sol1,"xn1",c_iteration);//DEBUG

  r1=cms_compute_residual(TMT_Matrix_,node_rhs,node_sol1);//DEBUG
  
  /* Build Edge RHS  (erhs = b - Ms *D0 * xn1) */
  ML_CHK_ERR(D0_Matrix_->Multiply(false,node_sol1,edge_temp1));
  ML_CHK_ERR(Ms_Matrix_->Multiply(false,edge_temp1,edge_rhs));
  ML_CHK_ERR(edge_rhs.Update(1.0,B,-1.0));
  MVOUT2(edge_rhs,"erhs",c_iteration);//DEBUG
  
  /* (1,1) Block Solve (xe = (S+M+Addon)^{-1} erhs) */
  ML_CHK_ERR(EdgePC->ApplyInverse(edge_rhs,edge_sol));
  r2=cms_compute_residual(SM_Matrix_,edge_rhs,edge_sol);//DEBUG
  MVOUT2(edge_sol,"xe1",c_iteration);//DEBUG
  
  /* Build Nodal RHS  (nrhs = D0'* (erhs - Ms * xe)) */
  ML_CHK_ERR(Ms_Matrix_->Multiply(false,edge_sol,edge_temp1));
  ML_CHK_ERR(edge_temp1.Update(1.0,edge_rhs,-1.0));
  ML_CHK_ERR(D0_Matrix_->Multiply(true,edge_temp1,node_rhs));
  MVOUT2(node_rhs,"nrhs2",c_iteration);//DEBUG
  
  /* (2,2) Block Solve (xn2 = TMT^{-1} nrhs) */
  ML_CHK_ERR(NodePC->ApplyInverse(node_rhs,node_sol2));
  MVOUT2(node_sol2,"xn2",c_iteration);//DEBUG
  
  /* Assemble solution (x = xe + T*(xn1 + xn2)) */
  ML_CHK_ERR(node_sol1.Update(1.0,node_sol2,1.0));

  r3=cms_compute_residual(TMT_Matrix_,node_rhs,node_sol1);//DEBUG

  
  ML_CHK_ERR(D0_Matrix_->Multiply(false,node_sol1,X));
  ML_CHK_ERR(X.Update(1.0,edge_sol,1.0));

  MVOUT2(X,"x",c_iteration);//DEBUG
  
  r4=cms_compute_residual(SM_Matrix_,B,X);//DEBUG  
  if(Comm_->MyPID()==0)
    printf("Residual Norms: %22.16e / %22.16e / %22.16e / %22.16e\n",r1,r2,r3,r4/r0);
    //    printf("Residual Norms: %6.4e / %6.4e / %6.4e / %6.4e\n",r1,r2,r3,r4/r0);

  c_iteration++;//DEBUG
  
  return 0;
}/*end ApplyInverse_Implicit_212*/

// ================================================ ====== ==== ==== == = 
//! Implicitly applies in the inverse in an additive format
int  ML_Epetra::RefMaxwellPreconditioner::ApplyInverse_Implicit_Additive(const Epetra_MultiVector& B, Epetra_MultiVector& X_) const
{

  int NumVectors=B.NumVectors();
  Epetra_MultiVector X(X_);
  Epetra_MultiVector TempE1(X.Map(),NumVectors,false);
  Epetra_MultiVector TempE2(X.Map(),NumVectors,true);
  Epetra_MultiVector TempN1(*NodeMap_,NumVectors,false);
  Epetra_MultiVector TempN2(*NodeMap_,NumVectors,true);
  Epetra_MultiVector Resid(B);
  
  double r0,r1,r,r2,r3,r4,r5;
  r0=cms_compute_residual(SM_Matrix_,B,X);//DEBUG


  MVOUT2(X,"a-x0",c_iteration);//DEBUG 
  MVOUT2(B,"a-b1",c_iteration);//DEBUG 

  
  /* Pre-Smoothing */
  ML_CHK_ERR(PreEdgeSmoother->ApplyInverse(B,X));

  MVOUT2(X,"a-x1",c_iteration);//DEBUG 
  r1=cms_compute_residual(SM_Matrix_,B,X);//DEBUG
  
  /* Build Residual */
  ML_CHK_ERR(SM_Matrix_->Multiply(false,X,TempE1));
  ML_CHK_ERR(Resid.Update(-1.0,TempE1,1.0));  
  if(!HasOnlyDirichletNodes) ML_CHK_ERR(D0_Matrix_->Multiply(true,Resid,TempN1));

  MVOUT2(Resid,"a-r1",c_iteration);//DEBUG
  MVOUT2(TempN1,"a-nr1",c_iteration);//DEBUG

  
  /* Precondition (1,1) block (additive)*/
  ML_CHK_ERR(EdgePC->ApplyInverse(Resid,TempE2));

  MVOUT2(TempE2,"a-p11",c_iteration);//DEBUG  
  r2=cms_compute_residual(SM_Matrix_,Resid,TempE2);//DEBUG
  
  /* Precondition (2,2) block (additive)*/
  if(!HasOnlyDirichletNodes){
    ML_CHK_ERR(NodePC->ApplyInverse(TempN1,TempN2));             

    MVOUT2(TempN2,"a-p22",c_iteration);//DEBUG  
    r3=cms_compute_residual(TMT_Matrix_,TempN1,TempN2);//DEBUG
    D0_Matrix_->Multiply(false,TempN2,TempE1);
  }/*end if*/
    
  /* Update solution */
  if(HasOnlyDirichletNodes) X.Update(1.0,TempE2,1.0);
  else X.Update(1.0,TempE1,1.0,TempE2,1.0);
  MVOUT2(X,"a-x2",c_iteration);//DEBUG
  r4=cms_compute_residual(SM_Matrix_,B,X);//DEBUG
  
  /* Post-Smoothing */
  ML_CHK_ERR(PostEdgeSmoother->ApplyInverse(B,X));
  MVOUT2(X,"a-x3",c_iteration);//DEBUG    
  r5=cms_compute_residual(SM_Matrix_,B,X);//DEBUG  
  c_iteration++;
  
  if(Comm_->MyPID()==0)
    printf("Residual Norms: %22.16e / %22.16e / %22.16e / %22.16e / %22.16e\n",r1/r0,r2,r3,r4/r0,r5/r0);

  X_=X;
  
  return 0;
}

// ================================================ ====== ==== ==== == = 
//! Implicitly applies in the inverse in an 1-2-1 format
int  ML_Epetra::RefMaxwellPreconditioner::ApplyInverse_Implicit_121(const Epetra_MultiVector& B, Epetra_MultiVector& X) const
{
  return -1;
}
  

// ================================================ ====== ==== ==== == = 

#ifdef STUFF_NOT_WORKING_YET
int ML_Epetra::SetDefaultsRefMaxwell(ParameterList & inList, 
			  int * options, double * params, bool OverWrite)
{
  /* Build Teuchos List: (1,1) */
  int dim=3, smooth=3;
  Teuchos::ParameterList List11;
  ML_Epetra::SetDefaultsSA(List11);
  List11.set("cycle applications",1);
  List11.set("aggregation: type","Uncoupled");
  List11.set("PDE equations",dim);
  List11.set("smoother: type","ML symmetric Gauss-Seidel");
  List11.set("smoother: sweeps (level 0)",0);
  //        List11.set("smoother: type","Chebyshev");  
  List11.set("smoother: sweeps",smooth);
  List11.set("x-coordinates",coordx);
  List11.set("y-coordinates",coordy);
  if(dim==3) List11.set("z-coordinates",coordz);
  else List11.set("z-coordinates",(double*)0);
  List11.set("output",10);
  
  /* Build Teuchos List: (2,2) */  
  Teuchos::ParameterList List22;
  ML_Epetra::SetDefaultsSA(List22);  
  List22.set("cycle applications",1);
  List22.set("smoother: type","ML symmetric Gauss-Seidel");
  //List22.set("smoother: type","Chebyshev");  
  List22.set("smoother: sweeps (level 0)",0);
  List22.set("smoother: sweeps",smooth);
  List22.set("x-coordinates",coordx);
  List22.set("y-coordinates",coordy); 
  if(dim==3) List22.set("z-coordinates",coordz);
  else List22.set("z-coordinates",(double*)0); 
  List22.set("output",10);
  
  /* Build Teuchos List: Fine Smoother */
  Teuchos::ParameterList ListSM;
  ListSM.set("coarse: type","Hiptmair");
  ListSM.set("coarse: sweeps",1);
  ListSM.set("smoother: Hiptmair efficient symmetric",true);
  ListSM.set("coarse: subsmoother type","symmetric Gauss-Seidel");
  ListSM.set("coarse: edge sweeps",smooth);
  ListSM.set("coarse: node sweeps",smooth);  
  ListSM.set("zero starting solution",false);  
  
  
  /* Build Teuchos List: Overall */  
  Teuchos::ParameterList ListRF;
  ListRF.set("refmaxwell: 11solver","edge matrix free");
  ListRF.set("refmaxwell: 11list",List11);
  ListRF.set("refmaxwell: 22solver","multilevel");
  ListRF.set("refmaxwell: 22list",List22);
  ListRF.set("refmaxwell: mode","additive");
  ListRF.set("refmaxwell: additive smoother",ListSM);
  
  for(ParameterList::ConstIterator param=List.begin(); param!=List.end(); param++)
    if ( inList.isParameter(inList.name(param)) == false || OverWrite )
      inList.setEntry( inList.name(param) , inList.entry(param) );

  return 0;
  
} //ML_Epetra::SetDefaultsMaxwell()
#endif


#ifdef STUFF_THAT_DOESNT_WORK_YET_BECUASE_ML_IS_FULL_OF_VISCIOUS_HACKS
  /* Variables + Constants we'll need */
  char parameter[80];
  int level=0,NumLevels_=1;
  int LevelID_[1]={ML_BOTH};
  int Mynum_smoother_steps=1;

  /* Build ML Object and stuff in SM */
  ML_Create(&ml_,NumLevels_);
  ML_Init_Amatrix(ml_,LevelID_[0],SM_Matrix_->NumMyRows(),SM_Matrix_->NumMyRows(), (void *) SM_Matrix_);
  ml_->Amat[LevelID_[0]].type = ML_TYPE_CRS_MATRIX;
  ml_->Amat[LevelID_[0]].N_nonzeros = SM_Matrix_->NumGlobalNonzeros();
  
  /* Allocate stuff the smoother needs */  
  nodal_args_ = ML_Smoother_Arglist_Create(2);
  edge_args_  = ML_Smoother_Arglist_Create(2);
  
  /* Pull List Options (General) */
  int num_smoother_steps = List_.get("smoother: sweeps", 2);
  double omega = List_.get("smoother: damping factor",1.0);
  int pre_or_post = 0;
  string PreOrPostSmoother = List_.get("smoother: pre or post","both");
  string Smoother = List_.get("smoother: type","Hiptmair");

  /* Pull List Options (Chebyshev) */
  int ChebyshevPolyOrder = List_.get("smoother: MLS polynomial order",-97);
  if (ChebyshevPolyOrder == -97) 
     ChebyshevPolyOrder = List_.get("smoother: polynomial order",-97);
  double ChebyshevAlpha = List_.get("smoother: MLS alpha",-2.0);
  if (ChebyshevAlpha == -2.) ChebyshevAlpha = List_.get("smoother: Chebyshev alpha", -2.0);
  int SmootherLevels = (NumLevels_>1)?(NumLevels_-1):1;

  /* Pull List Options (Hiptmair) */
  string SubSmootherType;
  int nodal_its = 1, edge_its = 1;
  SubSmootherType = List_.get("subsmoother: type","symmetric Gauss-Seidel"); //HAQ
  nodal_its = List_.get("subsmoother: node sweeps", 2);
  edge_its = List_.get("subsmoother: edge sweeps", 2);

  /* This stuff is cut and pasted out of MultiLevelPreconditioner::SetSmoothers */
  if( Smoother == "Hiptmair" ) {
    double subsmOmega;
    if (Comm().NumProc() == 1) subsmOmega = 1.0;
    else                       subsmOmega = ML_DDEFAULT;
    subsmOmega = List_.get("subsmoother: damping factor",subsmOmega);
    sprintf(parameter,"subsmoother: type (level %d)", level);
    string MySubSmootherType = List_.get(parameter,SubSmootherType);

    int logical_level = LevelID_[level];
    void *edge_smoother = 0, *nodal_smoother = 0;
    double node_coarsening_rate=0.0;
    double edge_coarsening_rate=0.0;
      
    sprintf(parameter,"subsmoother: node sweeps (level %d)", level);
    int Mynodal_its = List_.get(parameter, nodal_its);
    sprintf(parameter,"subsmoother: edge sweeps (level %d)", level);
    int Myedge_its = List_.get(parameter, edge_its);
    
    // Only SGS is supported
    if (MySubSmootherType == "symmetric Gauss-Seidel") {
      sprintf(parameter,"subsmoother: damping factor (level %d)",
              logical_level);
      double MysubsmOmega = List_.get(parameter,subsmOmega);
      nodal_smoother=(void *) ML_Gen_Smoother_SymGaussSeidel;
      ML_Smoother_Arglist_Set(nodal_args_, 0, &Mynodal_its);
      ML_Smoother_Arglist_Set(nodal_args_, 1, &MysubsmOmega);
      edge_smoother=(void *) ML_Gen_Smoother_SymGaussSeidel;
      ML_Smoother_Arglist_Set(edge_args_, 0, &Myedge_its);
      ML_Smoother_Arglist_Set(edge_args_, 1, &MysubsmOmega);
    }
    else if (Comm().MyPID() == 0)
      cerr << "Only SGS is supported as a Hiptmair subsmoother ... not " << MySubSmootherType << endl;
    
    int hiptmair_type = (int)
      List_.get("smoother: Hiptmair efficient symmetric", true);
    
    ML_Gen_Smoother_Hiptmair(ml_, logical_level, ML_BOTH,
                             Mynum_smoother_steps, Tmat_array, Tmat_trans_array, NULL, 
                             MassMatrix_array,
                             edge_smoother, edge_args_, nodal_smoother, nodal_args_,
                             hiptmair_type);
    
  }/*end if*/
    else if (Comm().MyPID() == 0)
      cerr << "Hiptmair is the only supported smoother type, not " << Smoother << endl;

  /* Cleanup */
#endif

#endif
