/*
 * Goal of this test:
 * - verify that the RefMaxwell Solver functions correctly
 *
 * This test will:
 * - Read matrices from disk
 * - Form a RefMaxwell preconditioner
 * - Test several combinations of options
 * - Verify the answer against a reference solution
 *
 * \date May 18, 2007
 *
 * \author Chris Siefert, 01414
 *
 */

#include <iostream>
#include <math.h>
#ifndef ICL
#ifndef WIN32
#include <unistd.h>
#endif
#endif
#include "ml_include.h"
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_EPETRAEXT)


#include "EpetraExt_MatrixMatrix.h"
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_SolverMap_CrsMatrix.h"
#include "EpetraExt_MultiVectorIn.h"
#include "ml_epetra_utils.h"

#include "ml_RefMaxwell.h"
#include "Epetra_LinearProblem.h"
#include "Teuchos_ParameterList.hpp"
#include "AztecOO.h"

#include "Teuchos_RCP.hpp"

using namespace ML_Epetra;

#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#define Epetra_ActiveComm Epetra_MpiComm
#else
#include "Epetra_SerialComm.h"
#define Epetra_ActiveComm Epetra_SerialComm
#endif

int MatlabFileToMultiVector(const char *filename, const Epetra_BlockMap & map, int N, Epetra_MultiVector * & A);
int dim=3;

/*******************************************************/
Teuchos::ParameterList Build_Teuchos_List(int LDA, double *coord_ptr,const char *str_tag, const char* str_val,const char *int_tag, int int_val,
					  const char *double_tag=0, double double_val=0,  const char *bool_tag=0, bool bool_val=0) {
  Teuchos::ParameterList List_Coarse, RMList;

  /* Pass in given options */
  List_Coarse.set("x-coordinates",&coord_ptr[0]);
  List_Coarse.set("y-coordinates",&coord_ptr[LDA]);
  List_Coarse.set("z-coordinates",&coord_ptr[2*LDA]);
  if(int_tag) List_Coarse.set(int_tag,int_val);
  if(str_tag) List_Coarse.set(str_tag,str_val);
  if(double_tag) List_Coarse.set(double_tag,double_val);
  if(bool_tag) List_Coarse.set(bool_tag,bool_val);
  List_Coarse.set("ML output",10);
  List_Coarse.set("coarse: type","Amesos-KLU");

  Teuchos::ParameterList List11(List_Coarse), List11c(List_Coarse), List22(List_Coarse);

  /* Set other necessary parameters */
  List11.set("smoother: sweeps",0);
  List11.set("PDE equations",3);
  List11c.set("PDE equations",3);

  /* Only needed because this problem is *really* small */
  //  List22.set("smoother: sweeps (level 0)",2);

  /* Setup basic list structure */
  List11.set("edge matrix free: coarse",List11c);
  RMList.setName("refmaxwell list");
  RMList.set("smoother: sweeps",2);
  RMList.set("refmaxwell: 11list",List11);
  RMList.set("refmaxwell: 22list",List22);
  RMList.set("ML output",10);
  if(int_tag) RMList.set(int_tag,int_val);
  if(str_tag) RMList.set(str_tag,str_val);
  SetDefaultsRefMaxwell(RMList,false);
  return RMList;
}


/*******************************************************/
void rpc_test_additive(Epetra_ActiveComm &Comm,
                       Teuchos::ParameterList & List,
                       Epetra_CrsMatrix &SM,
                       const Epetra_CrsMatrix &M1,
                       const Epetra_CrsMatrix &M0inv,
                       const Epetra_CrsMatrix &D0,
                       const Epetra_Vector &x_exact,
                       const Epetra_Vector &x0,
                       const Epetra_Vector &b,
		       bool run_gmres){

  // Generate Fake Ms
  Epetra_CrsMatrix *Ms=0;
  EpetraExt::MatrixMatrix::Add(M1,false,1.0,M1,false,1.0,Ms);
  Ms->FillComplete();

  RefMaxwellPreconditioner PrecRF(SM,D0,*Ms,M0inv,M1,List);
  Epetra_Vector x0_(x0);


  Epetra_LinearProblem Problem(&SM, &x0_, (Epetra_MultiVector*)&b);
  AztecOO solver(Problem);
  solver.SetPrecOperator(&PrecRF);

  if(run_gmres) solver.SetAztecOption(AZ_solver, AZ_gmres);
  else solver.SetAztecOption(AZ_solver, AZ_cg);
  solver.SetAztecOption(AZ_conv,AZ_r0);
  solver.SetAztecOption(AZ_output,1);
  solver.Iterate(100,1e-9);

  /* Check out the solution */
  double nxe,nd;
  x_exact.Norm2(&nxe);
  Epetra_Vector diff(x_exact);
  diff.Update(1.0,x0_,-1.0);
  diff.Norm2(&nd);
  if(Comm.MyPID()==0) printf("||sol-exact||/||exact||=%6.4e\n",nd/nxe);
  if(nd/nxe > 1e-8) exit(1);


  delete Ms;
}


/*******************************************************/
void rpc_test_additive_newconstructor(Epetra_ActiveComm &Comm,
                       Teuchos::ParameterList & List,
                       Epetra_CrsMatrix &SM,
                       const Epetra_CrsMatrix &M1,
                       const Epetra_CrsMatrix &M0inv,
                       const Epetra_CrsMatrix &D0,
                       const Epetra_Vector &x_exact,
                       const Epetra_Vector &x0,
                       const Epetra_Vector &b,
		       bool run_gmres){

  using Teuchos::rcp;
  // Generate Fake Ms
  Epetra_CrsMatrix *Ms=0;
  EpetraExt::MatrixMatrix::Add(M1,false,1.0,M1,false,1.0,Ms);
  Ms->FillComplete();

  // New List
  Teuchos::ParameterList NewList = List;
  NewList.set("M0inv",rcp(&M0inv,false));
  NewList.set("M1",rcp(&M1,false));
  NewList.set("D0",rcp(&D0,false));
  NewList.set("Ms",rcp((const Epetra_CrsMatrix*)Ms,false));

  RefMaxwellPreconditioner PrecRF(SM,NewList);
  Epetra_Vector x0_(x0);


  Epetra_LinearProblem Problem(&SM, &x0_, (Epetra_MultiVector*)&b);
  AztecOO solver(Problem);
  solver.SetPrecOperator(&PrecRF);

  if(run_gmres) solver.SetAztecOption(AZ_solver, AZ_gmres);
  else solver.SetAztecOption(AZ_solver, AZ_cg);
  solver.SetAztecOption(AZ_conv,AZ_r0);
  solver.SetAztecOption(AZ_output,1);
  solver.Iterate(100,1e-9);

  /* Check out the solution */
  double nxe,nd;
  x_exact.Norm2(&nxe);
  Epetra_Vector diff(x_exact);
  diff.Update(1.0,x0_,-1.0);
  diff.Norm2(&nd);
  if(Comm.MyPID()==0) printf("||sol-exact||/||exact||=%6.4e\n",nd/nxe);
  if(nd/nxe > 1e-8) exit(1);


  delete Ms;
}

/*******************************************************/
int rpc_test_additive_repeat(Epetra_ActiveComm &Comm,
                       Teuchos::ParameterList & List,
                       Epetra_CrsMatrix &SM,
                       const Epetra_CrsMatrix &M1,
                       const Epetra_CrsMatrix &M0inv,
                       const Epetra_CrsMatrix &D0,
                       const Epetra_Vector &x_exact,
                       const Epetra_Vector &x0,
                       const Epetra_Vector &b,
		       bool run_gmres){
  // Generate Fake Ms
  Epetra_CrsMatrix *Ms=0;
  EpetraExt::MatrixMatrix::Add(M1,false,1.0,M1,false,1.0,Ms);
  Ms->FillComplete();

  // First solve
  ML_set_random_seed(24601);
  RefMaxwellPreconditioner PrecRF(SM,D0,*Ms,M0inv,M1,List);
  Epetra_Vector x0_(x0);

  Epetra_LinearProblem Problem(&SM, &x0_, (Epetra_MultiVector*)&b);
  AztecOO solver(Problem);
  solver.SetPrecOperator(&PrecRF);

  if(run_gmres) solver.SetAztecOption(AZ_solver, AZ_gmres);
  else solver.SetAztecOption(AZ_solver, AZ_cg);
  solver.SetAztecOption(AZ_conv,AZ_r0);
  solver.SetAztecOption(AZ_output,1);
  solver.Iterate(100,1e-9);

  // Second solve
  ML_set_random_seed(24601);
  RefMaxwellPreconditioner PrecRF2(SM,D0,*Ms,M0inv,M1,List);
  Epetra_Vector x0_2(x0);

  Epetra_LinearProblem Problem2(&SM, &x0_2, (Epetra_MultiVector*)&b);
  AztecOO solver2(Problem2);
  solver2.SetPrecOperator(&PrecRF2);

  if(run_gmres) solver2.SetAztecOption(AZ_solver, AZ_gmres);
  else solver2.SetAztecOption(AZ_solver, AZ_cg);
  solver2.SetAztecOption(AZ_conv,AZ_r0);
  solver2.SetAztecOption(AZ_output,1);
  solver2.Iterate(100,1e-9);

  delete Ms;

  if (solver.NumIters() != solver2.NumIters()) {
    if(!x0.Comm().MyPID()) printf("ERROR: RNG control test failed (iters = %d vs. %d)!\n",solver.NumIters(),solver2.NumIters());
    return 1;
  }
  else return 0;
}





/*******************************************************/
bool matrix_read(Epetra_ActiveComm &Comm){
  Epetra_CrsMatrix *SM=0, *S, *Se,*D0,*D0e,*M0,*M1, *M1e;

  /* Read Matrices */
  EpetraExt::MatlabFileToCrsMatrix("S.dat" ,Comm,Se);
  EpetraExt::MatlabFileToCrsMatrix("M1.dat",Comm,M1e);
  EpetraExt::MatlabFileToCrsMatrix("M0.dat",Comm,M0);
  EpetraExt::MatlabFileToCrsMatrix("Tclean.dat",Comm,D0e);

  /* Fix Column Maps */
  EpetraExt::CrsMatrix_SolverMap S_CMT, D0_CMT, M0inv_CMT,M1_CMT;
  S =dynamic_cast<Epetra_CrsMatrix*>(ModifyEpetraMatrixColMap(*Se,S_CMT ));
  M1=dynamic_cast<Epetra_CrsMatrix*>(ModifyEpetraMatrixColMap(*M1e,M1_CMT));
  D0=dynamic_cast<Epetra_CrsMatrix*>(ModifyEpetraMatrixColMap(*D0e,D0_CMT));


  /* Optimize Storage */
  M1->FillComplete();
  M0->FillComplete();
  D0->FillComplete(D0e->DomainMap(),D0e->RangeMap());
  S->FillComplete();

  /* Build SM */
  EpetraExt::MatrixMatrix::Add(*S,false,1.0,*M1,false,1.0,SM);
  SM->FillComplete();

  /* Build RHS */
  Epetra_Map EdgeMap=SM->DomainMap();
  Epetra_Vector rhs(EdgeMap,false);
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

  /* Remap this bad boy */
  Epetra_CrsMatrix * M0inv=dynamic_cast<Epetra_CrsMatrix*>(ModifyEpetraMatrixColMap(M0inve,M0inv_CMT));
  M0inv->FillComplete();

  /* Read in nodal coordinates */
  int N;
  double *coord_ptr;
  Epetra_MultiVector *coords=0;
  MatlabFileToMultiVector("coord_node.dat",NodeMap,dim,coords);
  coords->ExtractView(&coord_ptr,&N);

  coords->Print(std::cout);

  /* Build the edge coordinates */
  Epetra_MultiVector n_coords_ghost(D0->ColMap(),dim);  
  if(D0->Importer()) n_coords_ghost.Import(*coords,*D0->Importer(),Add);
  else n_coords_ghost = *coords;

  Epetra_MultiVector e_coords(D0->RowMap(),dim);
  for(int i=0; i<D0->NumMyRows(); i++) { 
    int NumEntries;
    int * Indices;
    double * Values;      
    D0->ExtractMyRowView(i,NumEntries,Values,Indices);
    for(int j=0; j<NumEntries; j++)
      for (int k=0; k<dim; k++)
	e_coords[k][i] += n_coords_ghost[k][Indices[j]] / NumEntries;
  }     

  /* Build Lists */
  Teuchos::ParameterList List_2level = Build_Teuchos_List(N,coord_ptr,"coarse: type","Amesos-KLU","max levels",1);
  Teuchos::ParameterList List_SGS    = Build_Teuchos_List(N,coord_ptr,"smoother: type","symmetric Gauss-Seidel",0,1);
  Teuchos::ParameterList List_Cheby  = Build_Teuchos_List(N,coord_ptr,"smoother: type","Chebyshev",0,1);
  Teuchos::ParameterList List_SORa   = Build_Teuchos_List(N,coord_ptr,"smoother: type","Chebyshev",0,1);
  List_SORa.set("smoother: type","IFPACK");
  List_SORa.set("smoother: ifpack type","SORa");

  Teuchos::ParameterList List_Aux    = Build_Teuchos_List(N,coord_ptr,"smoother: type","Chebyshev",0,1,
							  "aggregation: aux: threshold",0.01, "aggregation: aux: enable",true);

  Teuchos::ParameterList List_LineSGS = Build_Teuchos_List(N,coord_ptr,"smoother: type","line Gauss-Seidel",
							   "coarse: max size",10,"smoother: line detection threshold",0.1);
  List_LineSGS.set("smoother: type","symmetric block Gauss-Seidel");
  List_LineSGS.set("smoother: line detection threshold",0.1);
  if(dim >= 1) List_LineSGS.set("x-coordinates",&e_coords[0][0]);
  if(dim >= 2) List_LineSGS.set("y-coordinates",&e_coords[1][0]);
  if(dim == 3) List_LineSGS.set("z-coordinates",&e_coords[2][0]);


  /* Do Tests */
  Epetra_Vector lhs(EdgeMap,true);
  if(!Comm.MyPID()) printf("*** Test 1 ***\n");
  rpc_test_additive(Comm,List_2level,*SM,*M1,*M0inv,*D0,x_exact,lhs,rhs,false);
  lhs.PutScalar(0.0);
  if(!Comm.MyPID()) printf("*** Test 2 ***\n");
  rpc_test_additive(Comm,List_SGS,*SM,*M1,*M0inv,*D0,x_exact,lhs,rhs,false);
  lhs.PutScalar(0.0);
  if(!Comm.MyPID()) printf("*** Test 3 ***\n");
  rpc_test_additive(Comm,List_Cheby,*SM,*M1,*M0inv,*D0,x_exact,lhs,rhs,false);
  lhs.PutScalar(0.0);
  if(!Comm.MyPID()) printf("*** Test 4 ***\n");
  rpc_test_additive(Comm,List_SORa,*SM,*M1,*M0inv,*D0,x_exact,lhs,rhs,true);
  lhs.PutScalar(0.0);
  if(!Comm.MyPID()) printf("*** Test 5 ***\n");
  rpc_test_additive(Comm,List_Aux,*SM,*M1,*M0inv,*D0,x_exact,lhs,rhs,false);

  // This only works in serial, for reasons which are not entirely obvious
  if(!Comm.MyPID()) printf("*** Test 6 ***\n");
  if(Comm.NumProc()==1) {
    lhs.PutScalar(0.0);
    rpc_test_additive(Comm,List_LineSGS,*SM,*M1,*M0inv,*D0,x_exact,lhs,rhs,false);
  }

  /* Check RNG control */
  int status1, status2 = 0;
  lhs.PutScalar(0.0);
  if(!Comm.MyPID()) printf("*** Test 7 ***\n");
  status1 = rpc_test_additive_repeat(Comm,List_SORa,*SM,*M1,*M0inv,*D0,x_exact,lhs,rhs,true);
  lhs.PutScalar(0.0);
  if(!Comm.MyPID()) printf("*** Test 8 ***\n");
  status2 = rpc_test_additive_repeat(Comm,List_Cheby,*SM,*M1,*M0inv,*D0,x_exact,lhs,rhs,true);

  /* Test w/ ParameterList-based constructor */
  lhs.PutScalar(0.0);
  if(!Comm.MyPID()) printf("*** Test 9 ***\n");
  rpc_test_additive_newconstructor(Comm,List_2level,*SM,*M1,*M0inv,*D0,x_exact,lhs,rhs,false);
  lhs.PutScalar(0.0);
  if(!Comm.MyPID()) printf("*** Test 10 ***\n");
  rpc_test_additive_newconstructor(Comm,List_SGS,*SM,*M1,*M0inv,*D0,x_exact,lhs,rhs,false);
  lhs.PutScalar(0.0);
  if(!Comm.MyPID()) printf("*** Test 11 ***\n");
  rpc_test_additive_newconstructor(Comm,List_Cheby,*SM,*M1,*M0inv,*D0,x_exact,lhs,rhs,false);
  lhs.PutScalar(0.0);
  if(!Comm.MyPID()) printf("*** Test 12 ***\n");
  rpc_test_additive_newconstructor(Comm,List_SORa,*SM,*M1,*M0inv,*D0,x_exact,lhs,rhs,true);
  lhs.PutScalar(0.0);
  if(!Comm.MyPID()) printf("*** Test 13 ***\n");
  rpc_test_additive_newconstructor(Comm,List_Aux,*SM,*M1,*M0inv,*D0,x_exact,lhs,rhs,false);

  delete M0; delete M1e;
  delete D0e;delete Se;
  delete SM;
  delete coords;

  return (status1||status2);
}




/*******************************************************/
int main(int argc, char* argv[]){
  /* Initialize */
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  int status = matrix_read(Comm);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  if(status==0)
    return EXIT_SUCCESS;
  else
    return EXIT_FAILURE;
}/*end main*/



/*******************************************************/
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
  for (int j=0; j<N; j++) {
    double * v = Ap[j];
    for (int i=0; i<numMyPoints; i++) {
      // Now read in each value and store to the local portion of the the  if the row is owned.
      double V;
      if(fscanf(handle, "%le ", &V)==0) return(-5);
      v[i] = V;
    }
  }

  if (fclose(handle)) return(-1);

  return(0);
}




#else

#ifdef HAVE_MPI
#include "mpi.h"
#endif

int main(int argc, char *argv[])
{

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  puts("This test requires:");
  puts("--enable-epetra");
  puts("--enable-aztecoo");
  puts("--enable-epetraext");
  puts("--enable-ifpack");
  puts("--enable-teuchos");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  // not to break tests
  return(EXIT_SUCCESS);
}

#endif
