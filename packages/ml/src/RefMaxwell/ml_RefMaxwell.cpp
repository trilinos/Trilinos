#include "ml_config.h"
#include <string.h>
#include "ml_RefMaxwell.h"
#include "ml_epetra.h"
#include "ml_MultiLevelPreconditioner.h"
#include "ml_RefMaxwell_11_Operator.h"
#include "ml_EdgeMatrixFreePreconditioner.h"
#include "EpetraExt_MatrixMatrix.h" //haq
extern void MVOUT (const Epetra_MultiVector & A, ostream & os);//HAQ
extern void MVOUT2(const Epetra_MultiVector & A,char* pref,int idx);//HAQ
extern void Epetra_CrsMatrix_Print(const Epetra_CrsMatrix& A, ostream& os);
#include "EpetraExt_CrsMatrixIn.h"//haq
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) 


// ================================================ ====== ==== ==== == = 
ML_Epetra::RefMaxwellPreconditioner::RefMaxwellPreconditioner(const Epetra_CrsMatrix& SM_Matrix,    //S+M
                                                              const Epetra_CrsMatrix& D0_Matrix,    //T or D0
                                                              const Epetra_CrsMatrix& Ms_Matrix,    //M1(sigma)
                                                              const Epetra_CrsMatrix& M0inv_Matrix, //M0^{-1}
                                                              const Epetra_CrsMatrix& M1_Matrix,    //M1(1)
                                                              const Epetra_CrsMatrix& TMT_Matrix,   //T' M1(sigma) T
                                                              const Teuchos::ParameterList& List,
                                                              const bool ComputePrec):
  ML_Preconditioner(),SM_Matrix_(&SM_Matrix),D0_Matrix_(&D0_Matrix),Ms_Matrix_(&Ms_Matrix),
  //M0inv_Matrix_(&M0inv_Matrix),M1_Matrix_(&M1_Matrix),Diagonal_(0),Operator11_(0)
  M0inv_Matrix_(&M0inv_Matrix),M1_Matrix_(&M1_Matrix),TMT_Matrix_(&TMT_Matrix),Diagonal_(0),Operator11_(0)
{

  // TOTAL HAQ
  //    EpetraExt::MatlabFileToCrsMatrix("tmt_matrix.dat" ,SM_Matrix_->Comm(),TMT_Matrix_);
  //    TMT_Matrix_->OptimizeStorage();
  //END HAQ

  /* Set the Epetra Goodies */
  Comm_ = &(SM_Matrix_->Comm());
  DomainMap_ = &(SM_Matrix_->OperatorDomainMap());
  RangeMap_ = &(SM_Matrix_->OperatorRangeMap());
  NodeMap_ = &(M0inv_Matrix_->OperatorDomainMap());
  
  Label_=strdup("ML reformulated Maxwell preconditioner"),
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
  
  /* Nuke everything if we've done this already */
  if(IsComputePreconditionerOK_) DestroyPreconditioner();

  /* Build the (1,1) Block Operator */
  Operator11_ = new ML_RefMaxwell_11_Operator(*SM_Matrix_,*D0_Matrix_,*M0inv_Matrix_,*M1_Matrix_);

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

  //  ofstream ofs1("addon_diagonal.dat");
  //  MVOUT(*Diagonal_,ofs1);


  SM_Matrix_->ExtractDiagonalCopy(temp_edge2);
  Diagonal_->Update(1.0,temp_edge2,1.0);
  
  //  ofstream ofs2("diagonal.dat");
  //  MVOUT(*Diagonal_,ofs2);
    
  /* Sanity Check the Diagonal */
  double min_val; Diagonal_->MinValue(&min_val);
  if(Comm_->MyPID()==0 && min_val <1e-16) {printf("ERROR: Minimum Estimated Diagonal <1e-16 (%6.4e)\n",min_val);return -2;}
  
  /* Build the (1,1) Block Preconditioner */
  string solver11=List_.get("refmaxwell: 11solver","edge matrix free");
  Teuchos::ParameterList List11=List_.get("refmaxwell: 11list",dummy);
  if(solver11=="edge matrix free") EdgePC=new EdgeMatrixFreePreconditioner(*Operator11_,*Diagonal_,*D0_Matrix_,*TMT_Matrix_,List11,true);
  else {printf("RefMaxwellPreconditioner: ERROR - Illegal (1,1) block preconditioner\n");return -1;}

  /* Build the (2,2) Block Preconditioner */
  string solver22=List_.get("refmaxwell: 22solver","multilevel");
  Teuchos::ParameterList List22=List_.get("refmaxwell: 22list",dummy);
  SetDefaultsSA(List22,0,0,false);
  List22.print(cout);//DEBUG


  ofstream ofs("tmt_matrix.dat");
  Epetra_CrsMatrix_Print(*TMT_Matrix_,ofs);//DEBUG

  
  ML_reseed_random_vec(8675309);//DEBUG
  if(solver22=="multilevel") NodePC=new MultiLevelPreconditioner(*TMT_Matrix_,List22);
  else {printf("RefMaxwellPreconditioner: ERROR - Illegal (2,2) block preconditioner\n");return -1;}
  //NTS: Add Adaptive, MatrixFree

  /* Pull Solver Mode */
  mode=List_.get("refmaxwell: mode","212");
  
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
}/*end DestroyPreconditioner*/



// ================================================ ====== ==== ==== == = 
// Apply the preconditioner to an Epetra_MultiVector X, puts the result in Y
int ML_Epetra::RefMaxwellPreconditioner::ApplyInverse(const Epetra_MultiVector& B, Epetra_MultiVector& X) const
{
  int rv;
  /* Sanity Checks */
  if (!B.Map().SameAs(*DomainMap_)) ML_CHK_ERR(-1);
  if (B.NumVectors() != X.NumVectors()) ML_CHK_ERR(-1);
  
  /* What mode to run in? */
  if(mode=="212") rv=ApplyInverse_Implicit_212(B,X);
  else if(mode=="additive") rv=ApplyInverse_Implicit_Additive(B,X);
  else if(mode=="121") rv=ApplyInverse_Implicit_121(B,X);
  else {fprintf(stderr,"RefMaxwellPreconditioner ERROR: Invalid ApplyInverse mode set in Teuchos list");ML_CHK_ERR(-2);}
  ML_CHK_ERR(rv);
  
  return 0;
}/*end ApplyInverse*/


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

  delete norm_old; delete_norm_new;
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


  delete norm_old; delete_norm_new;
  return rv;
}


// ================================================ ====== ==== ==== == = 
//! Implicitly applies in the inverse in a 2-1-2 format
int ML_Epetra::RefMaxwellPreconditioner::ApplyInverse_Implicit_212(const Epetra_MultiVector& B, Epetra_MultiVector& X) const
{
  int NumVectors=B.NumVectors();
  double r0,r1,r2,r3,r4;

  double norm;
  
  static int iteration=0;//DEBUG
  MVOUT2(B,"b",iteration);//DEBUG

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


  MVOUT2(node_rhs,"nrhs1",iteration);//DEBUG
  
  /* (2,2) Block Solve (xn1 = TMT^{-1} nrhs) */
  ML_reseed_random_vec(8675309);  
  ML_CHK_ERR(NodePC->ApplyInverse(node_rhs,node_sol1));
  MVOUT2(node_sol1,"xn1",iteration);//DEBUG

  r1=cms_compute_residual(TMT_Matrix_,node_rhs,node_sol1);//DEBUG
  
  /* Build Edge RHS  (erhs = b - Ms *D0 * xn1) */
  ML_CHK_ERR(D0_Matrix_->Multiply(false,node_sol1,edge_temp1));
  ML_CHK_ERR(Ms_Matrix_->Multiply(false,edge_temp1,edge_rhs));
  ML_CHK_ERR(edge_rhs.Update(1.0,B,-1.0));
  MVOUT2(edge_rhs,"erhs",iteration);//DEBUG
  
  /* (1,1) Block Solve (xe = (S+M+Addon)^{-1} erhs) */
  ML_CHK_ERR(EdgePC->ApplyInverse(edge_rhs,edge_sol));
  r2=cms_compute_residual(SM_Matrix_,edge_rhs,edge_sol);//DEBUG
  MVOUT2(edge_sol,"xe1",iteration);//DEBUG
  
  /* Build Nodal RHS  (nrhs = D0'* (erhs - Ms * xe)) */
  ML_CHK_ERR(Ms_Matrix_->Multiply(false,edge_sol,edge_temp1));
  ML_CHK_ERR(edge_temp1.Update(1.0,edge_rhs,-1.0));
  ML_CHK_ERR(D0_Matrix_->Multiply(true,edge_temp1,node_rhs));
  MVOUT2(node_rhs,"nrhs2",iteration);//DEBUG
  
  /* (2,2) Block Solve (xn2 = TMT^{-1} nrhs) */
  ML_CHK_ERR(NodePC->ApplyInverse(node_rhs,node_sol2));
  MVOUT2(node_sol2,"xn2",iteration);//DEBUG
  
  /* Assemble solution (x = xe + T*(xn1 + xn2)) */
  ML_CHK_ERR(node_sol1.Update(1.0,node_sol2,1.0));

  r3=cms_compute_residual(TMT_Matrix_,node_rhs,node_sol1);//DEBUG

  
  ML_CHK_ERR(D0_Matrix_->Multiply(false,node_sol1,X));
  ML_CHK_ERR(X.Update(1.0,edge_sol,1.0));

  MVOUT2(X,"x",iteration);//DEBUG
  
  r4=cms_compute_residual(SM_Matrix_,B,X);//DEBUG  
  if(Comm_->MyPID()==0)
    printf("Residual Norms: %6.4e / %6.4e / %6.4e / %6.4e\n",r1,r2,r3,r4/r0);


  iteration++;//DEBUG
  
  return 0;
}/*end ApplyInverse_Implicit_212*/

// ================================================ ====== ==== ==== == = 
//! Implicitly applies in the inverse in an additive format
int  ML_Epetra::RefMaxwellPreconditioner::ApplyInverse_Implicit_Additive(const Epetra_MultiVector& B, Epetra_MultiVector& X) const
{
  return -1;
}

// ================================================ ====== ==== ==== == = 
//! Implicitly applies in the inverse in an 1-2-1 format
int  ML_Epetra::RefMaxwellPreconditioner::ApplyInverse_Implicit_121(const Epetra_MultiVector& B, Epetra_MultiVector& X) const
{
  return -1;
}
  


#endif
