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

#define BASE_IDX 0
#define NO_OUTPUT


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
  if(A.Comm().NumProc()==1){
    int MaxNumIndices = A.MaxNumEntries();
    int* Indices;// = new int[MaxNumIndices];
    double* Values; //= new double[MaxNumIndices];
    int NumIndices;
    int i, j,NumMyRows=A.NumMyRows();
    FILE *f=fopen(of,"w");  
    for (i=0; i<NumMyRows; i++) {
      int Row = A.GRID(i); // Get global row number
      //      A.ExtractGlobalRowCopy(Row, MaxNumIndices, NumIndices, Values, Indices);      
      A.ExtractMyRowView(i, NumIndices, Values, Indices);      
      for (j = 0; j < NumIndices ; j++)
        fprintf(f,"%8d %8d %22.16e\n",Row,A.GCID(Indices[j]),Values[j]);   
    }/*end for*/ 
    fclose(f);
    //  delete [] Indices;
    //  delete [] Values;
  }/*end if*/
  else
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


void ML_Matrix_Print(ML_Operator *ML,const Epetra_Comm &Comm,const Epetra_Map &Map, char *fname){
  ML_Operator_Print(ML,fname);
}
#endif



// ================================================ ====== ==== ==== == = 
ML_Epetra::RefMaxwellPreconditioner::RefMaxwellPreconditioner(const Epetra_CrsMatrix& SM_Matrix,      //S+M
                                                              const Epetra_CrsMatrix& D0_Clean_Matrix,//T or D0 w/ nothing zero'd
                                                              const Epetra_CrsMatrix& Ms_Matrix,      //M1(sigma)
                                                              const Epetra_CrsMatrix& M0inv_Matrix,   //M0^{-1}
                                                              const Epetra_CrsMatrix& M1_Matrix,      //M1(1)
                                                              const Teuchos::ParameterList& List,
                                                              const bool ComputePrec):
  ML_Preconditioner(),SM_Matrix_(&SM_Matrix),D0_Matrix_(0), D0_Clean_Matrix_(&D0_Clean_Matrix),Ms_Matrix_(&Ms_Matrix),
  M0inv_Matrix_(&M0inv_Matrix),M1_Matrix_(&M1_Matrix),TMT_Matrix_(0),TMT_Agg_Matrix_(0),
    BCrows(0),numBCrows(0),HasOnlyDirichletNodes(false),Operator11_(0),EdgePC(0),NodePC(0),PreEdgeSmoother(0),PostEdgeSmoother(0),verbose_(false)
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

  /* Pull Solver Mode, verbosity */
  mode=List_.get("refmaxwell: mode","212");
  int vb_level=List_.get("output",0);
  if(vb_level >= 5) verbose_=true;
  else verbose_=false;
  
  /* Nuke everything if we've done this already */
  if(IsComputePreconditionerOK_) DestroyPreconditioner();

  /* Find the Dirichlet Rows (using SM_Matrix_) and columns (using D0_Clean_Matrix_) */
  BCrows=FindLocalDiricheltRowsFromOnesAndZeros(*SM_Matrix_,numBCrows);
  Epetra_IntVector * BCnodes=FindLocalDirichletColumnsFromRows(BCrows,numBCrows,*D0_Clean_Matrix_);   
  int Nn=BCnodes->MyLength();
  int numBCnodes=0;
  for(int i=0;i<Nn;i++){
    if((*BCnodes)[i]) numBCnodes++;
  }

  
  /* Sanity Check: We have at least some Dirichlet nodes */
  /* NTS: This should get fixed later for robustness */
  int HasInterior = numBCnodes != Nn;
  if(verbose_) printf("[%d] HasInterior = %d %d/%d\n",Comm_->MyPID(),HasInterior,Nn-numBCnodes,Nn);
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
  
  /* Build the TMT Matrix */  
  /* NTS: When ALEGRA builds this matrix itself, we get rid of these lines */
  //#define BUILD_WITH_EPETRAEXT
#ifdef BUILD_WITH_EPETRAEXT
  if(!HasOnlyDirichletNodes){
    Epetra_CrsMatrix temp1(Copy,*RangeMap_,BASE_IDX);
    TMT_Matrix_=new Epetra_CrsMatrix(Copy,*NodeMap_,BASE_IDX);
    //    EpetraExt::MatrixMatrix::Multiply(*Ms_Matrix_,false,*D0_Matrix_,false,temp1);
    EpetraExt::MatrixMatrix::Multiply(*SM_Matrix_,false,*D0_Matrix_,false,temp1);
    EpetraExt::MatrixMatrix::Multiply(*D0_Matrix_,true,temp1,false,*TMT_Matrix_);
    TMT_Matrix_->OptimizeStorage();
    Remove_Zeroed_Rows(*TMT_Matrix_);
  }/*end if*/
#else
  /* Build the TMT Matrix */
  if(!HasOnlyDirichletNodes){
    //    ML_Epetra_PtAP(*Ms_Matrix_,*D0_Matrix_,TMT_Matrix_,verbose_);
    ML_Epetra_PtAP(*SM_Matrix_,*D0_Matrix_,TMT_Matrix_,verbose_);
    TMT_Matrix_->OptimizeStorage();  
    Remove_Zeroed_Rows(*TMT_Matrix_);
  }/*end if */
#endif        


#ifdef BUILD_WITH_EPETRAEXT  
  /* Build the TMT-Agg Matrix  (used for aggregating the (1,1) block*/
  /* NTS: When ALEGRA builds this matrix itself, we get rid of these lines */
  Epetra_CrsMatrix temp2(Copy,*RangeMap_,0);
  TMT_Agg_Matrix_=new Epetra_CrsMatrix(Copy,*NodeMap_,0);  
  EpetraExt::MatrixMatrix::Multiply(*M1_Matrix_,false,*D0_Clean_Matrix_,false,temp2);
  EpetraExt::MatrixMatrix::Multiply(*D0_Clean_Matrix_,true,temp2,false,*TMT_Agg_Matrix_);
  Remove_Zeroed_Rows(*TMT_Agg_Matrix_);
#else
  /* Build the TMT-Agg Matrix */
  ML_Epetra_PtAP(*M1_Matrix_,*D0_Clean_Matrix_,TMT_Agg_Matrix_,verbose_);
  TMT_Agg_Matrix_->OptimizeStorage();    
  Remove_Zeroed_Rows(*TMT_Agg_Matrix_);
#endif

#ifndef NO_OUTPUT
  if(TMT_Matrix_) Epetra_CrsMatrix_Print(*TMT_Matrix_,"tmt_matrix.dat");
  Epetra_CrsMatrix_Print(*TMT_Agg_Matrix_,"tmt_agg_matrix.dat");
#endif

  /* Boundary nuke the edge matrices */
  Apply_OAZToMatrix(BCrows,numBCrows,*Ms_Matrix_);
  Apply_OAZToMatrix(BCrows,numBCrows,*M1_Matrix_);    


  /* DEBUG: Output matrices */
#ifndef NO_OUTPUT
  Epetra_CrsMatrix_Print(*SM_Matrix_,"sm_matrix.dat");
  Epetra_CrsMatrix_Print(*Ms_Matrix_,"ms_matrix.dat");  
  Epetra_CrsMatrix_Print(*M1_Matrix_,"m1_nuked.dat");  
  Epetra_CrsMatrix_Print(*D0_Matrix_,"d0_nuked.dat");  
  Epetra_CrsMatrix_Print(*D0_Clean_Matrix_,"d0_clean.dat");  
#endif


  /* Cleanup from the Boundary Conditions */
  delete BCnodes; 
  
#ifdef HAVE_ML_EPETRAEXT
  /* Fix the solver maps for ML / Epetra compatibility */
  SM_Matrix_ = dynamic_cast<Epetra_CrsMatrix*>(ModifyEpetraMatrixColMap(*SM_Matrix_,SM_Matrix_Trans_,"SM",verbose_));
  D0_Matrix_ = dynamic_cast<Epetra_CrsMatrix*>(ModifyEpetraMatrixColMap(*D0_Matrix_,D0_Matrix_Trans_,"D0",verbose_));
  D0_Clean_Matrix_ = dynamic_cast<Epetra_CrsMatrix*>(ModifyEpetraMatrixColMap(*D0_Clean_Matrix_,D0_Clean_Matrix_Trans_,"D0Clean",verbose_));
  Ms_Matrix_ = dynamic_cast<Epetra_CrsMatrix*>(ModifyEpetraMatrixColMap(*Ms_Matrix_,Ms_Matrix_Trans_,"Ms",verbose_));
  M1_Matrix_ = dynamic_cast<Epetra_CrsMatrix*>(ModifyEpetraMatrixColMap(*M1_Matrix_,M1_Matrix_Trans_,"M1",verbose_));
  M0inv_Matrix_ = dynamic_cast<Epetra_CrsMatrix*>(ModifyEpetraMatrixColMap(*M0inv_Matrix_,M0inv_Matrix_Trans_,"M0inv",verbose_));
  TMT_Matrix_ = dynamic_cast<Epetra_CrsMatrix*>(ModifyEpetraMatrixColMap(*TMT_Matrix_,TMT_Matrix_Trans_,"TMT",verbose_));
  TMT_Agg_Matrix_ = dynamic_cast<Epetra_CrsMatrix*>(ModifyEpetraMatrixColMap(*TMT_Agg_Matrix_,TMT_Agg_Matrix_Trans_,"TMTA",verbose_));
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
    Epetra_CrsMatrix M1diag(Copy,*DomainMap_,BASE_IDX,true);
    for(int i=0;i<M1diag.NumMyRows();i++){
      int gid=DomainMap_->GID(i);
      M1diag.InsertGlobalValues(gid,1,&(temp_edge2[i]),&gid);
    }/*end for*/
    M1diag.FillComplete();
    M1diag.OptimizeStorage();        
    
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

  /* Build the (1,1) Block Preconditioner */ 
  ML_reseed_random_vec(8675309);//DEBUG 
  string solver11=List_.get("refmaxwell: 11solver","edge matrix free");
  Teuchos::ParameterList List11=List_.get("refmaxwell: 11list",dummy);
  if(solver11=="edge matrix free")
    EdgePC=new EdgeMatrixFreePreconditioner(*Operator11_,*Diagonal_,*D0_Matrix_,*D0_Clean_Matrix_,*TMT_Agg_Matrix_,BCrows,numBCrows,List11,true);
  else {printf("RefMaxwellPreconditioner: ERROR - Illegal (1,1) block preconditioner\n");return -1;}

  
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
  Epetra_MultiVector X(X_);
  X.PutScalar(0);
  
  /* What mode to run in? */
  if(mode=="212") rv=ApplyInverse_Implicit_212(B,X);
  else if(mode=="additive") rv=ApplyInverse_Implicit_Additive(B,X);
  else if(mode=="121") rv=ApplyInverse_Implicit_121(B,X);
  else {fprintf(stderr,"RefMaxwellPreconditioner ERROR: Invalid ApplyInverse mode set in Teuchos list");ML_CHK_ERR(-2);}
  ML_CHK_ERR(rv);

  /* Copy work vector to output */
  X_=X;
  
  return 0;
}/*end ApplyInverse*/


// ================================================ ====== ==== ==== == = 
int ML_Epetra::RefMaxwellPreconditioner::SetEdgeSmoother(Teuchos::ParameterList &List1){  
  Teuchos::ParameterList dummy;
  Teuchos::ParameterList List = List1.get("refmaxwell: additive smoother",dummy);

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
    PreEdgeSmoother  = new MultiLevelPreconditioner(*SM_Matrix_,*D0_Matrix_,*TMT_Matrix_,PreList,true,true);
    PostEdgeSmoother = new MultiLevelPreconditioner(*SM_Matrix_,*D0_Matrix_,*TMT_Matrix_,PostList,true,true);
  }/*end if*/
  
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
  double r0=1,r1=1,r2=1,r3=1,r4=1;

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
  if(verbose_ && Comm_->MyPID()==0)
    printf("Residual Norms: %22.16e / %22.16e / %22.16e / %22.16e\n",r1,r2,r3,r4/r0);

  c_iteration++;//DEBUG
  
  return 0;
}/*end ApplyInverse_Implicit_212*/

// ================================================ ====== ==== ==== == = 
//! Implicitly applies in the inverse in an additive format
int  ML_Epetra::RefMaxwellPreconditioner::ApplyInverse_Implicit_Additive(const Epetra_MultiVector& B, Epetra_MultiVector& X) const
{

  int NumVectors=B.NumVectors();
  Epetra_MultiVector TempE1(X.Map(),NumVectors,false);
  Epetra_MultiVector TempE2(X.Map(),NumVectors,true);
  Epetra_MultiVector TempN1(*NodeMap_,NumVectors,false);
  Epetra_MultiVector TempN2(*NodeMap_,NumVectors,true);
  Epetra_MultiVector Resid(B);
  
  double r0=1,r1=1,r2=1,r3=1,r4=1,r5=1;
  r0=cms_compute_residual(SM_Matrix_,B,X);//DEBUG

  MVOUT2(X,"a-x0",c_iteration);//DEBUG 
  MVOUT2(B,"a-b1",c_iteration);//DEBUG 
  
  /* Pre-Smoothing */
  ML_CHK_ERR(PreEdgeSmoother->ApplyInverse(B,X));

  MVOUT2(X,"a-x1",c_iteration);//DEBUG 
  if(verbose_) r1=cms_compute_residual(SM_Matrix_,B,X);//DEBUG
  
  /* Build Residual */
  ML_CHK_ERR(SM_Matrix_->Multiply(false,X,TempE1));
  ML_CHK_ERR(Resid.Update(-1.0,TempE1,1.0));  
  if(!HasOnlyDirichletNodes) ML_CHK_ERR(D0_Matrix_->Multiply(true,Resid,TempN1));

  MVOUT2(Resid,"a-r1",c_iteration);//DEBUG
  MVOUT2(TempN1,"a-nr1",c_iteration);//DEBUG
  
  /* Precondition (1,1) block (additive)*/
  ML_CHK_ERR(EdgePC->ApplyInverse(Resid,TempE2));

  MVOUT2(TempE2,"a-p11",c_iteration);//DEBUG  
  if(verbose_) r2=cms_compute_residual(SM_Matrix_,Resid,TempE2);//DEBUG
  
  /* Precondition (2,2) block (additive)*/
  if(!HasOnlyDirichletNodes){
    ML_CHK_ERR(NodePC->ApplyInverse(TempN1,TempN2));             

    MVOUT2(TempN2,"a-p22",c_iteration);//DEBUG  
    if(verbose_) r3=cms_compute_residual(TMT_Matrix_,TempN1,TempN2);//DEBUG
    D0_Matrix_->Multiply(false,TempN2,TempE1);
  }/*end if*/
    
  /* Update solution */
  if(HasOnlyDirichletNodes) X.Update(1.0,TempE2,1.0);
  else X.Update(1.0,TempE1,1.0,TempE2,1.0);
  MVOUT2(X,"a-x2",c_iteration);//DEBUG
  if(verbose_) r4=cms_compute_residual(SM_Matrix_,B,X);//DEBUG
  
  /* Post-Smoothing */
  ML_CHK_ERR(PostEdgeSmoother->ApplyInverse(B,X));
  MVOUT2(X,"a-x3",c_iteration);//DEBUG    
  if(verbose_) r5=cms_compute_residual(SM_Matrix_,B,X);//DEBUG  
  c_iteration++;
  
  if(verbose_ && Comm_->MyPID()==0)
    printf("Residual Norms: %22.16e / %22.16e / %22.16e / %22.16e / %22.16e\n",r1/r0,r2,r3,r4/r0,r5/r0);
  
  return 0;
}

// ================================================ ====== ==== ==== == = 
//! Implicitly applies in the inverse in an 1-2-1 format
int  ML_Epetra::RefMaxwellPreconditioner::ApplyInverse_Implicit_121(const Epetra_MultiVector& B, Epetra_MultiVector& X) const
{
  return -1;
}


int ML_Epetra::UpdateList(Teuchos::ParameterList &source, Teuchos::ParameterList &dest, bool OverWrite){
  for(Teuchos::ParameterList::ConstIterator param=source.begin(); param!=source.end(); param++)
    if ( dest.isParameter(source.name(param)) == false || OverWrite )
      dest.setEntry(source.name(param),source.entry(param));
  return 0;
}

// ================================================ ====== ==== ==== == = 
int ML_Epetra::SetDefaultsRefMaxwell(Teuchos::ParameterList & inList,bool OverWrite)
{  
  const int smooth=3;
  const int default_output=1;
  //  const int default_output=10;
  
  /* Sublists */
  Teuchos::ParameterList ListRF,List11, List11c, List22, ListSM, dummy;
  Teuchos::ParameterList & List11_=inList.sublist("refmaxwell: 11list");
  Teuchos::ParameterList & List22_=inList.sublist("refmaxwell: 22list");
  Teuchos::ParameterList & ListSM_=inList.sublist("refmaxwell: additive smoother"); 
  Teuchos::ParameterList & List11c_=List11_.sublist("edge matrix free: coarse");

  /* Build Teuchos List: (1,1) coarse */    
  ML_Epetra::SetDefaultsSA(List11c);
  List11c.set("cycle applications",1);
  List11c.set("smoother: type","Chebyshev");
  //List11c.set("smoother: type","symmetric Gauss-Seidel");
  List11c.set("smoother: sweeps",smooth);
  List11c.set("coarse: type","Chebyshev");
  //List11c.set("coarse: type","symmetric Gauss-Seidel");  
  List11c.set("output",default_output);
  ML_Epetra::UpdateList(List11c,List11c_,OverWrite);
  
  /* Build Teuchos List: (1,1) */
  ML_Epetra::SetDefaultsSA(List11);
  List11.set("cycle applications",1);
  List11.set("aggregation: type","Uncoupled");
  List11.set("smoother: sweeps",0);
  List11.set("output",default_output);
  List11.set("edge matrix free: coarse",List11c);
  ML_Epetra::UpdateList(List11,List11_,OverWrite);
  
  /* Build Teuchos List: (2,2) */  
  ML_Epetra::SetDefaultsSA(List22);  
  List22.set("cycle applications",1);
  List22.set("smoother: type","Chebyshev");
  //List22.set("smoother: type","symmetric Gauss-Seidel");  
  List22.set("aggregation: type","MIS");
  List22.set("smoother: sweeps (level 0)",0);
  List22.set("smoother: sweeps",smooth);
  List22.set("coarse: type","Chebyshev");
  //List22.set("coarse: type","symmetric Gauss-Seidel");    
  List22.set("output",default_output);
  ML_Epetra::UpdateList(List22,List22_,OverWrite);  
  
  /* Build Teuchos List: Fine Smoother */
  ListSM.set("coarse: type","Hiptmair");
  ListSM.set("coarse: sweeps",1);
  ListSM.set("smoother: Hiptmair efficient symmetric",true);
  ListSM.set("coarse: subsmoother type","Chebyshev");
  //ListSM.set("coarse: subsmoother type","symmetric Gauss-Seidel");
  ListSM.set("coarse: edge sweeps",smooth);
  ListSM.set("coarse: node sweeps",smooth);  
  ListSM.set("zero starting solution",false);
  ListSM.set("max levels",1);  
  ListSM.set("output",default_output);
  ML_Epetra::UpdateList(ListSM,ListSM_,OverWrite);
  
  /* Build Teuchos List: Overall */  
  ListRF.set("refmaxwell: 11solver","edge matrix free");
  ListRF.set("refmaxwell: 11list",List11);
  ListRF.set("refmaxwell: 22solver","multilevel");
  ListRF.set("refmaxwell: 22list",List22);
  ListRF.set("refmaxwell: mode","additive");
  ListRF.set("refmaxwell: additive smoother",ListSM);
  ListRF.set("default values","RefMaxwell");
  ListRF.set("output",default_output);
  ML_Epetra::UpdateList(ListRF,inList,OverWrite);
  
  return 0;  
}/*end SetDefaultsRefMaxwell*/

#endif
