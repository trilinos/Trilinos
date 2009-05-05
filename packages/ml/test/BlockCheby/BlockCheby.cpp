#ifndef HAVE_CONFIG_H
#define HAVE_CONFIG_H
#endif

#include "ml_config.h"

#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_EPETRAEXT) && defined(HAVE_ML_AZTECOO)

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Time.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_VbrMatrix.h"
#include "Teuchos_ParameterList.hpp"
#include "ml_MultiLevelPreconditioner.h"
#include "AztecOO.h"
#include "ml_utils.h"
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_VectorIn.h"
#include "EpetraExt_MultiVectorIn.h"

using namespace Teuchos;

void PrintLine() 
{
  cout << endl;
  for( int i=0 ; i<80 ; ++i )
    cout << "=";
  cout << endl;
  cout << endl;
  
  return;
}

#define MAX(x,y) ((x)>(y)?(x):(y))


int TestMultiLevelPreconditioner(char ProblemType[],
				 Teuchos::ParameterList & MLList,
				 Epetra_LinearProblem & Problem,double & TotalErrorResidual,
				 double & TotalErrorExactSol,bool cg=false)
{
  
  Epetra_MultiVector* lhs = Problem.GetLHS();
  Epetra_MultiVector* rhs = Problem.GetRHS();
  Epetra_RowMatrix* A = Problem.GetMatrix();
  
  // ======================================== //
  // create a rhs corresponding to lhs or 1's //
  // ======================================== //

  ML_set_random_seed(987654);

  lhs->PutScalar(1.0);
    
  A->Multiply(false,*lhs,*rhs);
  lhs->PutScalar(0.0);
  
  Epetra_Time Time(A->Comm());
  
  // =================== //
  // call ML and AztecOO //
  // =================== //
  
  AztecOO solver(Problem);
  ML_Epetra::MultiLevelPreconditioner * MLPrec = new ML_Epetra::MultiLevelPreconditioner(*A, MLList, true);
  
  // tell AztecOO to use this preconditioner, then solve
  solver.SetPrecOperator(MLPrec);
  solver.SetAztecOption(AZ_solver, AZ_cg);
  solver.SetAztecOption(AZ_output, 10);
  
  solver.Iterate(1000, 1e-100);
  
  delete MLPrec;
  
  // ==================================================== //
  // compute difference between exact solution and ML one //
  // ==================================================== //
  
  double d = 0.0, d_tot = 0.0;
  
  for( int i=0 ; i<lhs->Map().NumMyElements() ; ++i )
    d += ((*lhs)[0][i] - 1.0) * ((*lhs)[0][i] - 1.0);
  
  A->Comm().SumAll(&d,&d_tot,1);
  
  // ================== //
  // compute ||Ax - b|| //
  // ================== //
  
  double Norm;
  Epetra_Vector Ax(rhs->Map());
  A->Multiply(false, *lhs, Ax);
  Ax.Update(1.0, *rhs, -1.0);
  Ax.Norm2(&Norm);
  
  string msg = ProblemType;
  
  if (A->Comm().MyPID() == 0) {
    cout << msg << "......Using " << A->Comm().NumProc() << " processes" << endl;
    cout << msg << "......||A x - b||_2 = " << Norm << endl;
    cout << msg << "......||x_exact - x||_2 = " << sqrt(d_tot) << endl;
    cout << msg << "......Total Time = " << Time.ElapsedTime() << endl;
  }
  
  TotalErrorExactSol += sqrt(d_tot);
  TotalErrorResidual += Norm;
  
  return( solver.NumIters() );
  
}

int main(int argc, char *argv[]) {

#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif
  
  // initialize the random number generator
  int ml_one = 1;
  ML_srandom1(&ml_one);
  // ===================== //
  // create linear problem //
  // ===================== //
  Epetra_CrsMatrix *BadMatrix;
  int * i_blockids;
  int numblocks;

  EpetraExt::MatlabFileToCrsMatrix("samplemat.dat",Comm,BadMatrix);
  const Epetra_Map *Map=&BadMatrix->RowMap();

  // Read in the block GLOBAL ids
  Epetra_Vector* d_blockids;
  int rv=EpetraExt::MatrixMarketFileToVector("blockids.dat",*Map,d_blockids);
  i_blockids=new int[d_blockids->MyLength()];
  numblocks=-1;
  for(int i=0;i<d_blockids->MyLength();i++){
    i_blockids[i]=(int)(*d_blockids)[i]-1;
    numblocks=MAX(numblocks,i_blockids[i]);
  }
  numblocks++;

  BadMatrix->FillComplete();  BadMatrix->OptimizeStorage();
  int N=BadMatrix->RowMatrixRowMap().NumMyElements();
  
  // Create the trivial blockID list
  int * trivial_blockids=new int[N];
  for(int i=0;i<N;i++)
    trivial_blockids[i]=i;    
  
  Epetra_Vector LHS(*Map); 
  Epetra_Vector RHS(*Map);
  Epetra_LinearProblem BadProblem(BadMatrix, &LHS, &RHS);

  Teuchos::ParameterList MLList;
  double TotalErrorResidual = 0.0, TotalErrorExactSol = 0.0;
  char mystring[80];

#ifdef OLD
  // ====================== //
  // ML Cheby 
  // ====================== //
  if (Comm.MyPID() == 0) PrintLine();
  ML_Epetra::SetDefaults("SA",MLList);
  MLList.set("smoother: type","Chebyshev");
  MLList.set("coarse: type","Amesos-KLU");    
  MLList.set("max levels",2);
  MLList.set("ML output",10);
  MLList.set("smoother: polynomial order",2);  
  strcpy(mystring,"Cheby");
  TestMultiLevelPreconditioner(mystring, MLList, BadProblem,
                               TotalErrorResidual, TotalErrorExactSol);

  // These tests only work in serial due to how the block id's are written.
  if(Comm.NumProc()==1){
    // ====================== //
    // ML Block Cheby (Trivial)
    // ====================== //
    if (Comm.MyPID() == 0) PrintLine(); 
    ML_Epetra::SetDefaults("SA",MLList); 
    MLList.set("smoother: type","Block Chebyshev");
    MLList.set("smoother: Block Chebyshev number of blocks",N);
    MLList.set("smoother: Block Chebyshev block list",trivial_blockids);
    MLList.set("coarse: type","Amesos-KLU");  
    MLList.set("max levels",2);
    MLList.set("ML output",10);  
    MLList.set("smoother: polynomial order",2);
    strcpy(mystring,"ML Block Cheby (Trivial)");
    TestMultiLevelPreconditioner(mystring, MLList, BadProblem,
                                 TotalErrorResidual, TotalErrorExactSol);
  
  
    // ====================== //
    // ML Block Cheby (Smart)
    // ====================== //  
    if (Comm.MyPID() == 0) PrintLine();
    ML_Epetra::SetDefaults("SA",MLList);  
    MLList.set("smoother: type","Block Chebyshev");
    MLList.set("smoother: Block Chebyshev number of blocks",numblocks);
    MLList.set("smoother: Block Chebyshev block list",i_blockids);    
    MLList.set("coarse: type","Amesos-KLU");  
    MLList.set("max levels",2);
    MLList.set("ML output",10);  
    MLList.set("smoother: polynomial order",2);
    strcpy(mystring,"ML Block Cheby (Smart)");
    TestMultiLevelPreconditioner(mystring, MLList, BadProblem,
                                 TotalErrorResidual, TotalErrorExactSol);
  }
#endif
  
#if defined(HAVE_ML_IFPACK)
#ifdef OLD
  // ====================== //
  // IFPACK Cheby 
  // ====================== //
  if (Comm.MyPID() == 0) PrintLine();
  ML_Epetra::SetDefaults("SA",MLList);
  MLList.set("smoother: type","IFPACK-Chebyshev");
  MLList.set("coarse: type","Amesos-KLU");    
  MLList.set("max levels",2);
  MLList.set("ML output",10);
  MLList.set("smoother: polynomial order",2);  
  strcpy(mystring,"IFPACK Cheby");
  TestMultiLevelPreconditioner(mystring, MLList, BadProblem,
                               TotalErrorResidual, TotalErrorExactSol);
#endif

  // ====================== //
  // IFPACK Block Cheby (Trivial)
  // ====================== //  
  int NumBlocks=Map->NumMyElements();
  int *BlockStarts=new int[NumBlocks+1];
  int *Blockids=new int [NumBlocks];
#ifdef OLD
  for(int i=0;i<NumBlocks;i++){
    BlockStarts[i]=i;
    Blockids[i]=Map->GID(i);
  }
  BlockStarts[NumBlocks]=NumBlocks;
  if (Comm.MyPID() == 0) PrintLine();
  ML_Epetra::SetDefaults("SA",MLList);
  MLList.set("smoother: type","IFPACK-Block Chebyshev");
  MLList.set("smoother: Block Chebyshev number of blocks",NumBlocks);
  MLList.set("smoother: Block Chebyshev block starts",&BlockStarts[0]);
  MLList.set("smoother: Block Chebyshev block list",&Blockids[0]);    
  MLList.set("coarse: type","Amesos-KLU");  
  MLList.set("max levels",2);
  MLList.set("ML output",10);  
  MLList.set("smoother: polynomial order",2);
  strcpy(mystring,"IFPACK Block Cheby (Trivial)");
  TestMultiLevelPreconditioner(mystring, MLList, BadProblem,
                               TotalErrorResidual, TotalErrorExactSol);
#endif
  delete [] BlockStarts; delete [] Blockids;

  
  // ====================== //
  // IFPACK Block Cheby (Smart)
  // ====================== //
  // Figure out how many blocks we actually have and build a map...
  Epetra_Map* IfpackMap;
  int g_NumBlocks=-1,g_MaxSize=-1;
  if(Comm.MyPID() == 0){
    const int lineLength = 1025;
    char line[lineLength];

    FILE * f=fopen("localids_in_blocks.dat","r");
    assert(f!=0);

    // Next, strip off header lines (which start with "%")
    do {
      if(fgets(line, lineLength,f)==0) return(-4);
    } while (line[0] == '%');

    // Grab the number we actually care about
    sscanf(line, "%d %d", &g_NumBlocks, &g_MaxSize);
    fclose(f);
  }
  Comm.Broadcast(&g_NumBlocks,1,0);
  Comm.Broadcast(&g_MaxSize,1,0);
  Epetra_Map BlockMap(g_NumBlocks,0,Comm);   
  Epetra_MultiVector *blockids_disk=0;
  rv=EpetraExt::MatrixMarketFileToMultiVector("localids_in_blocks.dat",BlockMap,blockids_disk);

  // Put all the block info into the right place
  NumBlocks=BlockMap.NumMyElements();
  BlockStarts=new int[NumBlocks+1];
  Blockids= new int[g_MaxSize*NumBlocks];
  // NTS: Blockids_ is overallocated because I don't want to write a counting loop
  int i,cidx;
  for(i=0,cidx=0;i<NumBlocks;i++){
    BlockStarts[i]=cidx;
    Blockids[cidx]=(int)(*blockids_disk)[0][i];cidx++;
    if((*blockids_disk)[1][i] > 1e-2){
      Blockids[cidx]=(int)(*blockids_disk)[1][i];cidx++;
    }    
  }
  BlockStarts[NumBlocks]=cidx;
  
  
  if (Comm.MyPID() == 0) PrintLine();
  ML_Epetra::SetDefaults("SA",MLList);
  MLList.set("smoother: type","IFPACK-Block Chebyshev");
  MLList.set("smoother: Block Chebyshev number of blocks",NumBlocks);
  MLList.set("smoother: Block Chebyshev block starts",&BlockStarts[0]);
  MLList.set("smoother: Block Chebyshev block list",&Blockids[0]);    
  MLList.set("coarse: type","Amesos-KLU");  
  MLList.set("max levels",2);
  MLList.set("ML output",10);  
  MLList.set("smoother: polynomial order",2);
  strcpy(mystring,"IFPACK Block Cheby (Smart)");
  TestMultiLevelPreconditioner(mystring, MLList, BadProblem,
                               TotalErrorResidual, TotalErrorExactSol);
  
  delete blockids_disk; delete [] BlockStarts; delete [] Blockids;
#endif
  


  
  // ===================== //
  // print out total error //
  // ===================== //

  if (Comm.MyPID() == 0) {
    cout << endl;
    cout << "......Total error for residual        = " << TotalErrorResidual << endl;
    cout << "......Total error for exact solution  = " << TotalErrorExactSol << endl;
    cout << endl;
  }

  delete [] i_blockids;
  delete [] trivial_blockids;
  delete BadMatrix;
  delete d_blockids;
  if (TotalErrorResidual > 1e-5) {
    cerr << "Error: `BlockCheby.exe' failed!" << endl;
    exit(EXIT_FAILURE);
  }

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  if (Comm.MyPID() == 0)
    cerr << "`BlockCheby.exe' passed!" << endl;

  return (EXIT_SUCCESS);
}

#else

#include <stdlib.h>
#include <stdio.h>
#ifdef HAVE_MPI
#include "mpi.h"
#endif

int main(int argc, char *argv[])
{
  // still need to deal with MPI, some architecture don't like
  // an exit(0) without MPI_Finalize()
#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  puts("Please configure ML with --enable-epetra --enable-teuchos --enable-galeri --enable-aztecoo");

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);
}

#endif /* #if defined(ML_WITH_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_GALERI) */
