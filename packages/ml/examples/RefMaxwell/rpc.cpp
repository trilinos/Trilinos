#include <stdio.h>
#include <math.h>
#include "ml_config.h"
#include "ml_include.h"
//#include "Epetra_CrsMatrix.h"
//#include "Epetra_Map.h"
//#include "Epetra_Time.h"
#include "Epetra_Export.h"
#include "EpetraExt_MatrixMatrix.h"
//#include "Epetra_Multi_CrsMatrix.h"
#include "ml_RefMaxwell_11_Operator.h"
#include "ml_RefMaxwell.h"
#include "ml_EdgeMatrixFreePreconditioner.h"
//#include "ml_operator.h"
//#include "ml_rap.h"
#include "ml_epetra_utils.h"
//#include "ml_utils.h"
//#include "ml_xyt.h"

#include "Epetra_Util.h"


//#include "Trilinos_Util.h"
#include "Epetra_LocalMap.h"
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_Transpose_RowMatrix.h"
#include "EpetraExt_SolverMap_CrsMatrix.h"
#include "EpetraExt_MultiVectorIn.h"

/* Aztec Stuff*/
#include "Epetra_LinearProblem.h"
#include "Teuchos_ParameterList.hpp"
#include "AztecOO.h"


#ifdef USE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

using namespace ML_Epetra;

//#define DIAGS


//#define PAVEL
//#define FILE_OUTPUT

int    NumIters =50;


#ifdef PAVEL
int dim=2;
#else
int dim=3;
#endif


Epetra_RowMatrix* ModifyEpetraMatrixColMap(const Epetra_RowMatrix &A,
                                           EpetraExt::CrsMatrix_SolverMap &transform);

int MatlabFileToMultiVector(const char *filename, const Epetra_BlockMap & map, int N, Epetra_MultiVector * & A);
void MVOUT(const Epetra_MultiVector & A, ostream & os);

#ifdef USE_MPI
void matrix_read(Epetra_MpiComm &Comm);
void rpc_test(Epetra_MpiComm &Comm,
              const Epetra_CrsMatrix &S,
              Epetra_CrsMatrix &SM,
              const Epetra_CrsMatrix &Ms,
              const Epetra_CrsMatrix &M1,
              const Epetra_CrsMatrix &M0inv,
              const Epetra_CrsMatrix &D0,
              const Epetra_MultiVector &coords,
              const Epetra_Vector &x_exact,
              const Epetra_Vector &x0,
              const Epetra_Vector &b);
#else
void matrix_read(Epetra_SerialComm &Comm);
void rpc_test(Epetra_SerialComm &Comm,
              const Epetra_CrsMatrix &S,
              Epetra_CrsMatrix &SM,              
              const Epetra_CrsMatrix &Ms,
              const Epetra_CrsMatrix &M1,
              const Epetra_CrsMatrix &M0inv,
              const Epetra_CrsMatrix &D0,
              const Epetra_MultiVector &coords,
              const Epetra_Vector &x_exact,
              const Epetra_Vector &x0,
              const Epetra_Vector &b);
#endif

int main(int argc, char* argv[]){ 
  char *dir;
  int rv,slen;
  /* Initialize */
#ifdef USE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  /* Sanity */
  if(Comm.MyPID()==0){
   if(argc<2) exit(1);
   dir=argv[1];
   slen=strlen(dir);
  }
  else dir=new char[80];

   
  /* I can't believe that Epetra_Comm doesn't support moving strings */
  //NTS: This bombs with path names below 5 characters
#ifdef USE_MPI
  Comm.Broadcast(&slen,1,0);
  MPI_Bcast(dir,slen,MPI_CHAR,0,MPI_COMM_WORLD);
  rv=chdir(dir);
  if(rv) {printf("Directory does not exist\n");exit(1);}
#endif

  /* Readin  + Test*/
  matrix_read(Comm);


  /* Cleanup*/
  if(Comm.MyPID()!=0) delete dir;
  
  
#ifdef USE_MPI
  MPI_Finalize();
#endif

  return 0;
}/*end main*/


void Epetra_CrsMatrix_Print(const Epetra_CrsMatrix& A, ostream& os);

//   EdgeMatrixFreePreconditioner(const Epetra_Operator_With_MatMat & Operator, const Epetra_Vector& Diagonal,
//                                 const Epetra_CrsMatrix & D0_matrix,const Epetra_CrsMatrix & TMT_matrix,
//                                 const Teuchos::ParameterList &List,const bool ComputePrec = true);

void print_stats(const Epetra_CrsMatrix& A, char *label){
  printf("[%d] %s[global]: (size) %dx%d (R/D) %dx%d (R/C) %d/%d\n",A.Comm().MyPID(),label,
         A.NumGlobalRows(),A.NumGlobalCols(),
         A.RangeMap().NumGlobalElements(),A.DomainMap().NumGlobalElements(),
         A.RowMap().NumGlobalElements(),A.ColMap().NumGlobalElements());
  printf("[%d] %s[local]: (size) %dx%d (R/D) %dx%d (R/C) %d/%d\n",A.Comm().MyPID(),label,
         A.NumMyRows(),A.NumMyCols(),
         A.RangeMap().NumMyElements(),A.DomainMap().NumMyElements(),
         A.RowMap().NumMyElements(),A.ColMap().NumMyElements());
}/*end print_stats*/


/******************************************/
/******************************************/
/******************************************/
#ifdef USE_MPI
void matrix_read(Epetra_MpiComm &Comm){
#else
void matrix_read(Epetra_SerialComm &Comm){
#endif
  Epetra_CrsMatrix *SM,*SMe,*Se,*S,*Ms,*Mse, *D0,*D0e,*M0,*M1, *M1e;

  /* Read Matrices */
  EpetraExt::MatlabFileToCrsMatrix("S.dat" ,Comm,Se);
  EpetraExt::MatlabFileToCrsMatrix("M1.dat",Comm,M1e);
  EpetraExt::MatlabFileToCrsMatrix("M0.dat",Comm,M0);
  EpetraExt::MatlabFileToCrsMatrix("Tclean.dat",Comm,D0e);
#ifdef PAVEL
  EpetraExt::MatlabFileToCrsMatrix("Ms.dat" ,Comm,Mse);
#else
  EpetraExt::MatlabFileToCrsMatrix("SM.dat" ,Comm,SMe);
#endif


  /* Optimize Storage*/
  M1e->OptimizeStorage();
  D0e->OptimizeStorage();
#ifdef PAVEL
  Mse->OptimizeStorage();
#else
  SMe->OptimizeStorage();
#endif


  
 /* Fix up the column maps (since ML needs all the columns, and Epetra doesn't
    supply them).

    This code will need to migrate into the actual preconditioner wrapper, but
    putting inside the xfer function or the (1,1) block matvec doesn't make sense.
 */
  EpetraExt::CrsMatrix_SolverMap S_CMT,SM_CMT, D0_CMT, Ms_CMT, M0inv_CMT,M1_CMT;
  S =dynamic_cast<Epetra_CrsMatrix*>(ModifyEpetraMatrixColMap(*Se, S_CMT ));
  M1=dynamic_cast<Epetra_CrsMatrix*>(ModifyEpetraMatrixColMap(*M1e,M1_CMT));
  D0=dynamic_cast<Epetra_CrsMatrix*>(ModifyEpetraMatrixColMap(*D0e,D0_CMT));

#ifdef PAVEL
  Ms=dynamic_cast<Epetra_CrsMatrix*>(ModifyEpetraMatrixColMap(*Mse,Ms_CMT));
#else
  SM=dynamic_cast<Epetra_CrsMatrix*>(ModifyEpetraMatrixColMap(*SMe,SM_CMT));
#endif
  
  
  /* Build SM  (pavel matrices) or Ms (alegra)*/
#ifdef PAVEL
  SM=new Epetra_CrsMatrix(*S);
  EpetraExt::MatrixMatrix::Add(*Ms,false,1,*SM,1);
#else
  Ms=new Epetra_CrsMatrix(*SM);  
  EpetraExt::MatrixMatrix::Add(*S,false,-1,*Ms,1);  
#endif

  /* Optimize the storage*/
  S->OptimizeStorage();
  SM->OptimizeStorage();
  Ms->OptimizeStorage();
  M1->OptimizeStorage();
  M0->OptimizeStorage();
  D0->OptimizeStorage();


  /* Build LHS/RHS */
  Epetra_Map EdgeMap=SM->DomainMap();  
  Epetra_Vector rhs(EdgeMap,false);  
  Epetra_Vector lhs(EdgeMap,true);
  Epetra_Vector x_exact(EdgeMap,false);
  x_exact.PutScalar(1.0);
  SM->Multiply(false,x_exact,rhs);

  /* Build Lumped M0^{-1} */
  Epetra_Map NodeMap=M0->DomainMap();
  Epetra_Vector ones(NodeMap),diag(NodeMap),invdiag(NodeMap);
  ones.PutScalar(1.0);
  M0->Multiply(false,ones,diag);
  
  Epetra_CrsMatrix M0inve(Copy,NodeMap,1);
  invdiag.Reciprocal(diag);
  for(int i=0;i<M0inve.NumMyRows();i++){
    int gid=NodeMap.GID(i);
    M0inve.InsertGlobalValues(gid,1,&(invdiag[i]),&gid);
  }/*end for*/  
  M0inve.FillComplete();
  M0inve.OptimizeStorage();

  
  /* Remap this bad boy */
  Epetra_CrsMatrix * M0inv=dynamic_cast<Epetra_CrsMatrix*>(ModifyEpetraMatrixColMap(M0inve,M0inv_CMT));  
  M0inv->OptimizeStorage();

  /* Read in coordinates*/
  Epetra_MultiVector *coords;
  MatlabFileToMultiVector("coord_node.txt",NodeMap,dim,coords);
  
  /* Tests */
  rpc_test(Comm,*S,*SM,*Ms,*M1,*M0inv,*D0,*coords,x_exact,lhs,rhs);
  
  /* Cleanup */
  // The CrsMatrix_SolverMap cleans up the non-E matrices.  I'd call this
  // "convenient," but that would be a savage, savage, lie.
  delete M0; delete M1e;
  delete D0e;delete Se; 
  delete coords; 
#ifdef PAVEL
  delete SM;
  delete Mse;
#else
  delete Ms;
  delete SMe;
#endif


} 

/******************************************/
/******************************************/
/******************************************/
#ifdef USE_MPI
void rpc_test(Epetra_MpiComm &Comm,
#else
void rpc_test(Epetra_SerialComm &Comm,
#endif
              const Epetra_CrsMatrix &S,
              Epetra_CrsMatrix &SM,              
              const Epetra_CrsMatrix &Ms,
              const Epetra_CrsMatrix &M1,
              const Epetra_CrsMatrix &M0inv,
              const Epetra_CrsMatrix &D0,
              const Epetra_MultiVector &coords,
              const Epetra_Vector &x_exact,
              const Epetra_Vector &x0,
              const Epetra_Vector &b){
  
  /* Pull the Maps */
  const Epetra_Map &EdgeMap=SM.DomainMap();
  const Epetra_Map &NodeMap=M0inv.DomainMap();

  /* Build the TMT (Nodal Laplacian, L0) Matrix */
  printf("[%d] EPC: Building TMT/L0 Matrix\n",Comm.MyPID());
  Epetra_CrsMatrix m_temp(Copy,EdgeMap,0);
  Epetra_CrsMatrix L0(Copy,NodeMap,0);  
  //  print_stats(D0,"D0");
  EpetraExt::MatrixMatrix::Multiply(Ms,false,D0,false,m_temp);
  m_temp.FillComplete(NodeMap,EdgeMap);
  EpetraExt::MatrixMatrix::Multiply(D0,true,m_temp,false,L0);
  L0.FillComplete();
  L0.OptimizeStorage();

  //  ofstream ofs("l0.dat");
  //  Epetra_CrsMatrix_Print(L0,ofs);

  int smooth=3;
  
  printf("[%d] RPC: Building Teuchos Lists\n",Comm.MyPID());
  /* Build Teuchos List: (1,1) */  
  Teuchos::ParameterList List11;
  SetDefaultsSA(List11);
  List11.set("cycle applications",1);
  List11.set("aggregation: type","Uncoupled");
  List11.set("PDE equations",dim);
  List11.set("smoother: type","MLS");  
  List11.set("smoother: sweeps",smooth);
  List11.set("smoother: MLS polynomial order",smooth);
  List11.set("eigen-analysis: type", "power-method");
  List11.set("eigen-analysis: max iters",100);
  List11.set("chebyshev: alpha",30.0001);
  List11.set("x-coordinates",coords[0]);
  List11.set("y-coordinates",coords[1]);
  if(dim==3) List11.set("z-coordinates",coords[2]);
  else List11.set("z-coordinates",(double*)0);
  List11.set("output",10);
  
  /* Build Teuchos List: (2,2) */  
  Teuchos::ParameterList List22;  
  List22.set("cycle applications",1);
  List22.set("smoother: type","MLS");
  List22.set("smoother: sweeps",smooth);
  List22.set("smoother: MLS polynomial order",smooth);
  List22.set("smoother: MLS alpha",30.0001);
  List22.set("coarse: type","MLS");
  List22.set("coarse: MLS polynomial order",smooth);
  List22.set("eigen-analysis: type", "power-method");  
  List22.set("x-coordinates",coords[0]);
  List22.set("y-coordinates",coords[1]); 
  if(dim==3) List22.set("z-coordinates",coords[2]);
  else List22.set("z-coordinates",(double*)0); 
  List22.set("output",10);

  /* Build Teuchos List: Overall */  
  Teuchos::ParameterList ListRF;
  ListRF.set("refmaxwell: 11solver","edge matrix free");
  ListRF.set("refmaxwell: 11list",List11);
  ListRF.set("refmaxwell: 22solver","multilevel");
  ListRF.set("refmaxwell: 22list",List22);
  ListRF.set("refmaxwell: mode","212");
  
  /* Build the (1,1) Block preconditioner */
  printf("[%d] RPC: Building the RexMaxwellPreconditioner\n",Comm.MyPID());
  ML_reseed_random_vec(8675309);
  //  RefMaxwellPreconditioner PrecRF(SM,D0,Ms,M0inv,M1,L0,ListRF);
  RefMaxwellPreconditioner PrecRF(SM,D0,Ms,M0inv,M1,ListRF);
  
  /* Build Sample Vector */
  Epetra_Vector x0_(x0);
  //  Epetra_Vector x_exact(EdgeMap,true);
  //  Epetra_Vector x0(EdgeMap,true);
  //  Epetra_Vector rhs(EdgeMap,false);
  //  x_exact.PutScalar(1.0);
  //  SM.Apply(x_exact,rhs);

  /* Aztec Setup */
  Epetra_LinearProblem Problem(&SM, &x0_, (Epetra_MultiVector*)&b);
  AztecOO solver(Problem);
  solver.SetPrecOperator(&PrecRF);

  /* Get solver options from Teuchos list */
  double Tol      = 1e-10;
  string type     = "gmres";
  int    output   = 1;
  string conv     = "r0";

  /* Set solver options - Solver type*/
  if (type == "cg") solver.SetAztecOption(AZ_solver, AZ_cg);
  else if (type == "cg_condnum") solver.SetAztecOption(AZ_solver, AZ_cg_condnum);
  else if (type == "gmres") solver.SetAztecOption(AZ_solver, AZ_gmres);
  else if (type == "gmres_condnum") solver.SetAztecOption(AZ_solver, AZ_gmres_condnum);
  else if (type == "fixed point") solver.SetAztecOption(AZ_solver, AZ_fixed_pt);
  
  /* Set solver options - Convergence Criterion*/
  if(conv == "r0") solver.SetAztecOption(AZ_conv,AZ_r0);
  else if(conv == "rhs") solver.SetAztecOption(AZ_conv,AZ_rhs);
  else if(conv == "Anorm") solver.SetAztecOption(AZ_conv,AZ_Anorm);
  else if(conv == "noscaled") solver.SetAztecOption(AZ_conv,AZ_noscaled);
  else if(conv == "sol") solver.SetAztecOption(AZ_conv,AZ_sol);

  /* Set solver options - other */
  solver.SetAztecOption(AZ_output, output);

  /* Do the solve */
  solver.Iterate(NumIters, Tol);

  /* Check out the solution */
  double nxe,nd;
  x_exact.Norm2(&nxe);  
  Epetra_Vector diff(x_exact);
  diff.Update(1.0,x0_,-1.0);
  diff.Norm2(&nd);  
  if(Comm.MyPID()==0) printf("||sol-exact||/||exact||=%6.4e\n",nd/nxe);
  
}/*end epc_test*/
             
                 
/******************************************/
/******************************************/
/******************************************/
                 
void Epetra_CrsMatrix_Print(const Epetra_CrsMatrix& A, ostream& os) {
  int MyPID = A.RowMap().Comm().MyPID();
  int NumProc = A.RowMap().Comm().NumProc();
       
  for (int iproc=0; iproc < NumProc; iproc++) {
    if (MyPID==iproc) {
      int NumMyRows1 = A.NumMyRows();
      int MaxNumIndices = A.MaxNumEntries();
      int * Indices  = new int[MaxNumIndices];
      double * Values  = new double[MaxNumIndices];
      int NumIndices;
      int i, j;
      
      for (i=0; i<NumMyRows1; i++) {
	int Row = A.GRID(i); // Get global row number
	A.ExtractGlobalRowCopy(Row, MaxNumIndices, NumIndices, Values, Indices);
        
	for (j = 0; j < NumIndices ; j++) {   
	  os.width(10);
	  os <<  Row +1; os << "    ";	
	  os.width(10);
	  os <<  Indices[j]+1; os << "    ";
	  os.width(20);os.precision(16);os.setf(ios_base::scientific,ios_base::floatfield);
	  os <<  Values[j]; os << "    ";
	  os << endl;
	}/*end for*/
      }/*end for*/
			
      delete [] Indices;
      delete [] Values;
      
      os << flush;      
    }/*end if*/
    // Do a few global ops to give I/O a chance to complete
    A.RowMap().Comm().Barrier();
    A.RowMap().Comm().Barrier();
    A.RowMap().Comm().Barrier();
  }/*end for*/
	
  return;
}/*end Epetra_CrsMatrix_Print*/
 

             
Epetra_RowMatrix* ModifyEpetraMatrixColMap(const Epetra_RowMatrix &A,
                                           EpetraExt::CrsMatrix_SolverMap &transform){
    Epetra_RowMatrix *B;
    Epetra_CrsMatrix *Acrs;

    const Epetra_CrsMatrix *Atmp = dynamic_cast<const Epetra_CrsMatrix*>(&A);
    if (Atmp != 0) {
      Acrs = const_cast<Epetra_CrsMatrix*>(Atmp);
      B = &(transform(*Acrs));
    }
    else
      B = const_cast<Epetra_RowMatrix *>(&A);


    if (B != &A)
      printf("** Transforming column map of matrix\n");
    else
      printf("** Leaving column map of matrix unchanged\n");
    
    return B;
} //ModifyEpetraMatrixColMap()


/******************************************/
/******************************************/
/******************************************/
void MVOUT(const Epetra_MultiVector & A, ostream & os){

#ifdef FILE_OUTPUT
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
#endif
  
}/*end MultiVectorToMatlabFile*/


void MVOUT2(const Epetra_MultiVector & A,char* pref,int idx){
#ifdef FILE_OUTPUT
  char c[80];
  sprintf(c,"%s.%d.dat",pref,idx);
  ofstream ofs(c);
  MVOUT(A,ofs);
#endif
}/* end MVOUT2*/


              
/******************************************/
/******************************************/
/******************************************/
void IVOUT(const Epetra_IntVector & A, ostream & os){
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

              
/******************************************/
/******************************************/
/******************************************/
int MatlabFileToMultiVector( const char *filename, const Epetra_BlockMap & map, int N, Epetra_MultiVector * & A) {

  const int lineLength = 1025;
  char line[lineLength];

  FILE * handle = 0;

  handle = fopen(filename,"r");  // Open file
  if (handle == 0)
    EPETRA_CHK_ERR(-1); // file not found

  // Next, strip off header lines (which start with "%")
  //  do {
  //    if(fgets(line, lineLength, handle)==0) return(-1);
  //  } while (line[0] == '%');

  // Compute the offset for each processor for when it should start storing values
  int numMyPoints = map.NumMyPoints();
  int offset=0;
  map.Comm().ScanSum(&numMyPoints, &offset, 1); // ScanSum will compute offsets for us
  offset -= numMyPoints; // readjust for my PE
  if(map.Comm().NumProc() == 1) offset=0;//CMS
  
  // Now construct vector/multivector
  if (N==1)
    A = new Epetra_Vector(map);
  else
    A = new Epetra_MultiVector(map, N);

  double ** Ap = A->Pointers();

  // Now read in lines that we will discard
  for (int i=0; i<offset; i++)
    if(fgets(line, lineLength, handle)==0) return(-3);
  for (int i=0; i<numMyPoints; i++) {
    for (int j=0; j<N; j++) {
      double * v = Ap[j];    
      // Now read in each value and store to the local portion of the the  if the row is owned.
      double V;
      if(fscanf(handle, "%le ", &V)==0) return(-5);
      v[i] = V;
    }
  }

  if (fclose(handle)) return(-1);
  
  return(0);
}
