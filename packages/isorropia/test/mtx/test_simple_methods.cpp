
#include <cstring>
#include <cstdio>
#include <iostream>
#include <fstream>

#include <Epetra_ConfigDefs.h>
#include <Epetra_Time.h>
#include <Epetra_Import.h>
#include <Epetra_Map.h>
#include <Epetra_LocalMap.h>
#include <Epetra_CrsGraph.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#include <Epetra_MultiVector.h>
#include <EpetraExt_MatrixMatrix.h>
#include <EpetraExt_BlockMapIn.h>
#include <EpetraExt_CrsMatrixIn.h>
#include <EpetraExt_RowMatrixOut.h>


#define SUCCESS 0
#define FAIL 1

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#ifdef EPETRA_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include <Teuchos_ParameterList.hpp> 

#include <Isorropia_EpetraPartitioner.hpp>
#include <ispatest_utils.hpp>
#include <ispatest_epetra_utils.hpp>

///////////////////////////////////////
// Answers
///////////////////////////////////////
int blockAnswer2[100] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                         1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                         1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
int blockAnswer4[100] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                         1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                         2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
                         3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3};
int cyclicAnswer2[100] = {0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,
                          1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,
                          0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,
                          1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1};
int cyclicAnswer4[100] = {0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,
                          1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,
                          2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,
                          3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3};
///////////////////////////////////////



int getMatrixDim(char *mmfilename, Epetra_Comm &comm,
                 int *N, int *M,  int *nnz);

int read_matrixmarket_linear(const char *, Epetra_Comm &, Epetra_CrsMatrix *&);

bool testPartitioningMVInput(Teuchos::RCP<const Epetra_MultiVector> mv, const char *method, int worldsize, int myrank);
bool testPartitioningMatrixInput(Teuchos::RCP<const Epetra_RowMatrix> mat, const char *method, int worldsize, int myrank);
bool testPartitioningMapInput(Teuchos::RCP<const Epetra_BlockMap> map, const char *method, int worldsize, int myrank);
bool checkPartition1(const Isorropia::Epetra::Partitioner &partitioner);
bool checkPartition2(const Isorropia::Epetra::Partitioner &partitioner, const Epetra_BlockMap &map, const char* method);
bool checkPartition4(const Isorropia::Epetra::Partitioner &partitioner, const Epetra_BlockMap &map, const char* method);


////////////////////////////////////////////////////////////////////////////////
// Testing simple methods with Multivector
////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv) 
{

#ifdef EPETRA_MPI
  MPI_Init(&argc,&argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  Epetra_CrsMatrix *A;

  int worldsize = Comm.NumProc();
  int myrank = Comm.MyPID();
  bool success = true;

  /////////////////////////////////////////////////////////////
  /// Read and distribute matrices
  /////////////////////////////////////////////////////////////
  read_matrixmarket_linear("arrow100.mtx",Comm,A);
  /////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////
  // Test partitioning with matrix input
  /////////////////////////////////////////////////////////////
  Teuchos::RCP<const Epetra_RowMatrix> rowmatrix = Teuchos::rcp(A);

  success = success && testPartitioningMatrixInput(rowmatrix,"BLOCK",worldsize,myrank);
  success = testPartitioningMatrixInput(rowmatrix,"CYCLIC",worldsize,myrank) && success; // test lhs of && so it is run
  /////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////
  // Test partitioning with multivector input
  /////////////////////////////////////////////////////////////
  Teuchos::RCP<const Epetra_MultiVector> mv = Teuchos::rcp(new Epetra_MultiVector(A->RowMap(),1));

  success = testPartitioningMVInput(mv,"BLOCK",worldsize,myrank) && success;
  success = testPartitioningMVInput(mv,"CYCLIC",worldsize,myrank) && success;
  /////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////
  // Test partitioning with BlockMap input
  /////////////////////////////////////////////////////////////
  Teuchos::RCP<const Epetra_BlockMap> map = Teuchos::rcp(&(A->RowMap()),false);

  success = testPartitioningMapInput(map,"BLOCK",worldsize,myrank) && success;
  success = testPartitioningMapInput(map,"CYCLIC",worldsize,myrank) && success;
  /////////////////////////////////////////////////////////////

#ifdef EPETRA_MPI
  MPI_Finalize();
#endif

  if(!success)
  {
    return FAIL;
  }

  return(SUCCESS);
}
////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
bool testPartitioningMVInput(Teuchos::RCP<const Epetra_MultiVector> mv, const char *method, int worldsize, int myrank)
{
  bool success=true;

  /////////////////////////////////////////////////////////////
  // Partition
  /////////////////////////////////////////////////////////////
  Teuchos::ParameterList paramlist;
  paramlist.set("IMBALANCE TOL","1.03");

  if(strcmp(method,"BLOCK")==0)
  {
    paramlist.set("PARTITIONING_METHOD", "BLOCK");
  }
  else if(strcmp(method,"CYCLIC")==0)
  {
    paramlist.set("PARTITIONING_METHOD", "CYCLIC");
  }
  else 
  {
    if(myrank==0)
    {
      std::cout   << "Error in testPartitionMVInput() with method argument " << method << "." << std::endl
                  << "      Method must be 'BLOCK' or 'CYCLIC'" << std::endl;
    }
    return(false);
  }

  Isorropia::Epetra::Partitioner partitioner(mv, paramlist, false);
  partitioner.partition();
  /////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////
  /// Check partition
  /////////////////////////////////////////////////////////////
  if (worldsize==1)
  {
    success = checkPartition1(partitioner);
  }
  else if (worldsize==2)
  {
    success = checkPartition2(partitioner,mv->Map(),method);
  }
  else if (worldsize==4)
  {
    success = checkPartition4(partitioner,mv->Map(),method);
  }
  else 
  {
    if(myrank==0)
    {
      std::cout << "WARNING for testPartitionMVInput() with number of processors = " << worldsize << std::endl
                << "        Number of processes not supported .... (by default) returning success" << std::endl;
    }
  }

  if(myrank==0)
  {
    if (success == false)
    {
      std::cout   << "Failure in testPartitionMVInput(), method argument " << method 
                  << ", numProcs = " << worldsize << std::endl;
    }
    else
    {
      std::cout   << "Success in testPartitionMVInput(), method argument " << method 
                  << ", numProcs = " << worldsize << std::endl;
    }
  }
  /////////////////////////////////////////////////////////////

  return success;
}
////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
bool testPartitioningMatrixInput(Teuchos::RCP<const Epetra_RowMatrix> mat, const char *method, int worldsize, int myrank)
{
  bool success=true;

  /////////////////////////////////////////////////////////////
  // Partition
  /////////////////////////////////////////////////////////////
  Teuchos::ParameterList paramlist;
  paramlist.set("IMBALANCE TOL","1.03");

  if(strcmp(method,"BLOCK")==0)
  {
    paramlist.set("PARTITIONING_METHOD", "BLOCK");
  }
  else if(strcmp(method,"CYCLIC")==0)
  {
    paramlist.set("PARTITIONING_METHOD", "CYCLIC");
  }
  else 
  {
    if(myrank==0) 
    {
      std::cout   << "Error in testPartitionMatrixInput() with method argument " << method << "." << std::endl
                  << "      Method must be 'BLOCK' or 'CYCLIC'" << std::endl;
    }
    return(false);
  }

  Isorropia::Epetra::Partitioner partitioner(mat, paramlist, false);
  partitioner.partition();
  /////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////
  /// Check partition
  /////////////////////////////////////////////////////////////
  if (worldsize==1)
  {
    success = checkPartition1(partitioner);
  }
  else if (worldsize==2)
  {
    success = checkPartition2(partitioner,mat->Map(),method);
  }
  else if (worldsize==4)
  {
    success = checkPartition4(partitioner,mat->Map(),method);
  }
  else 
  {
    if(myrank==0)
    {
      std::cout << "WARNING for testPartitionMatrixInput() with number of processors = " << worldsize << std::endl
                << "        Number of processes not supported .... (by default) returning success" << std::endl;
    }
  }

  if(myrank==0)
  {
    if (success == false)
    {
      std::cout   << "Failure in testPartitionMatrixInput(), method argument " << method 
                  << ", numProcs = " << worldsize << std::endl;
    }
    else
    {
      std::cout   << "Success in testPartitionMatrixInput(), method argument " << method 
                  << ", numProcs = " << worldsize << std::endl;
    }
  }
  /////////////////////////////////////////////////////////////

  return success;
}
////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
bool testPartitioningMapInput(Teuchos::RCP<const Epetra_BlockMap> map, const char *method, int worldsize, int myrank)
{
  bool success=true;

  /////////////////////////////////////////////////////////////
  // Partition
  /////////////////////////////////////////////////////////////
  Teuchos::ParameterList paramlist;
  paramlist.set("IMBALANCE TOL","1.03");

  if(strcmp(method,"BLOCK")==0)
  {
    paramlist.set("PARTITIONING_METHOD", "BLOCK");
  }
  else if(strcmp(method,"CYCLIC")==0)
  {
    paramlist.set("PARTITIONING_METHOD", "CYCLIC");
  }
  else 
  {
    if(myrank==0)
    {
      std::cout   << "Error in testPartitionMapInput() with method argument " << method << "." << std::endl
                  << "      Method must be 'BLOCK' or 'CYCLIC'" << std::endl;
    }
    return(false);
  }

  Isorropia::Epetra::Partitioner partitioner(map, paramlist, false);
  partitioner.partition();
  /////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////
  /// Check partition
  /////////////////////////////////////////////////////////////
  if (worldsize==1)
  {
    success = checkPartition1(partitioner);
  }
  else if (worldsize==2)
  {
    success = checkPartition2(partitioner,*map, method);
  }
  else if (worldsize==4)
  {
    success = checkPartition4(partitioner,*map, method);
  }
  else
  {
    if(myrank==0)
    {
      std::cout << "WARNING for testPartitionMapInput() with number of processors = " << worldsize << std::endl
                << "        Number of processes not supported .... (by default) returning success" << std::endl;
    }
  }

  if(myrank==0)
  {
    if (success == false)
    {
      std::cout   << "Failure in testPartitionMapInput(), method argument " << method 
                  << ", numProcs = " << worldsize << std::endl;
    }
    else
    {
      std::cout   << "Success in testPartitionMapInput(), method argument " << method 
                  << ", numProcs = " << worldsize << std::endl;
    }
  }
  /////////////////////////////////////////////////////////////

  return success;
}
////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
bool checkPartition1(const Isorropia::Epetra::Partitioner &partitioner)
{
  int numRows = 100;

  for(int i=0;i<numRows;i++)
  {
    if(partitioner[i] !=0) // For worldsize=1, everything should be assigned to part 0
    {
      return false;
    }
  }
  return true;
}
////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
bool checkPartition2(const Isorropia::Epetra::Partitioner &partitioner, 
                     const Epetra_BlockMap &map, const char *method)
{
  int success = 1;
  int numRows = map.NumMyElements();

  for(int i=0;i<numRows;i++)
  {
    if ( strcmp(method,"BLOCK")==0 && partitioner[i] != blockAnswer2[map.GID(i)])
    {      
      success = 0;
      break;
    }
    else if ( strcmp(method,"CYCLIC")==0 && partitioner[i] != cyclicAnswer2[map.GID(i)])
    {      
      success = 0;
      break;
    }
  }

#ifdef EPETRA_MPI
  int returnFlag;
  MPI_Allreduce(&success,&returnFlag,1,MPI_INT,MPI_LAND,MPI_COMM_WORLD);

  if(returnFlag==0)
  {
    return false;
  }
#endif

  return true;
}
////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
bool checkPartition4(const Isorropia::Epetra::Partitioner &partitioner, 
                     const Epetra_BlockMap &map, const char *method)
{
  int success = 1;
  int numRows = map.NumMyElements();

  for(int i=0;i<numRows;i++)
  {
    if ( strcmp(method,"BLOCK")==0 && partitioner[i] != blockAnswer4[map.GID(i)])
    {      
      success = 0;
      break;
    }
    else if ( strcmp(method,"CYCLIC")==0 && partitioner[i] != cyclicAnswer4[map.GID(i)])
    {      
      success = 0;
      break;
    }
  }

#ifdef EPETRA_MPI
  int returnFlag;
  MPI_Allreduce(&success,&returnFlag,1,MPI_INT,MPI_LAND,MPI_COMM_WORLD);

  if(returnFlag==0)
  {
    return false;
  }
#endif

  return true;
}
////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
int read_matrixmarket_linear(const char *infilename, Epetra_Comm &comm,
			     Epetra_CrsMatrix *&A)
{

  char mmfilename[256];
  strcpy(mmfilename,infilename);

  Epetra_Map *rowMap;

  // Construct a Map that puts approximately the same number of
  // equations on each processor.  Row-based decomposition.
  // Peek into MatrixMarket file to get Matrix Dimensions
  int nrows, ncols, nnz;
  if (!getMatrixDim(mmfilename, comm, &nrows, &ncols, &nnz))
  {
    return 0;
  }

  rowMap = new Epetra_Map(nrows, 0, comm);

  // Read MatrixMarket file into matrix A. 1D row partitioning
  //rowMap, colMap, rangeMap, domainMap
  EpetraExt::MatrixMarketFileToCrsMatrix(mmfilename, *rowMap, A, false, false);

  return 0;
}
////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////
int getMatrixDim(char *mmfilename, Epetra_Comm &comm,
		 int *N,   // Output:  Number of rows.
		 int *M,   // Output:  Number of columns.
		 int *nnz  // Output:  Number of nonzeros.
		 )
{
  //  Peek into the MatrixMarket file to get the matrix dimensions.
  //  Proc 0 does the peeking and broadcasts to other procs.
  int tmp[3]={0,0,0};
  int aok = 1;

  if (comm.MyPID() == 0) {

    std::ifstream mm(mmfilename);

    if (!mm) {
      cout << "getMatrixDim:  Cannot open MatrixMarket " << mmfilename << endl;
      *N = *M = -1;
      aok = 0;
    }
    else {
      // Skip comments in MatrixMarket header.
      char c = mm.peek();
      char line[100];
      while (c == '%') {
        mm.getline(line, 100);
        c = mm.peek();
        continue;
      }
      // Get N, M and nnz from MatrixMarket file.
      mm >> tmp[0] >> tmp[1] >> tmp[2];
      mm.close();
    }
  }
  comm.Broadcast(tmp, 3, 0);
  *N = tmp[0];
  *M = tmp[1];
  *nnz = tmp[2];
  if (*N < 1 || *M < 1 || *N != *M) {
    if (comm.MyPID()==0) {
      cout << "getMatrixDim:  Invalid input in matrix " << mmfilename << endl;
      cout << "getMatrixDim:  N=" << *N << " M=" << *M
           << " nnz=" << *nnz << endl;
    }
    aok = 0;
  }

  return aok;
}
////////////////////////////////////////////////////////////////////////////




