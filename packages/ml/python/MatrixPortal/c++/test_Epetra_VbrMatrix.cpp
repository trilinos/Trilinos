// Construct a VBR matrix of type
//
//     |  *   *            |
//     |      *    *       | 
// A = |           ... ... |
//     |                *  |
//
// (`*' being a nonzero block). The block size of block-row
// is will be i+1.

#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
#include "Epetra_VbrMatrix.h"
#include "Epetra_IntSerialDenseVector.h"

int main(int argc, char *argv[]) 
{
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  // set global dimension to 5, could be any number
  int NumGlobalElements = 5;
  
  // create a linear map
  Epetra_Map Map(NumGlobalElements,0,Comm);
  
  // local number of rows
  int NumMyElements = Map.NumMyElements();
  // get update list
  int* MyGlobalElements = Map.MyGlobalElements( );

  // dimension of each block
  Epetra_IntSerialDenseVector ElementSizeList(NumMyElements);

  // now construct a funky matrix. Diagonal block of block row i will have
  // dimension i+1 (don't run this code with too many nodes...). The
  // dimension of each block row is recordered in ElementSizeList.
  // Here ElementSizeList is declared as Epetra_IntSerialDenseVector,
  // but an int array is fine as well. 
  // max_blk keeps trace of the max block dimension
  
  int max_blk = 0;
  
  for (int i = 0; i < NumMyElements; ++i) 
  {
    ElementSizeList[i] = 1 + MyGlobalElements[i];
    if (ElementSizeList[i] > max_blk) max_blk =  ElementSizeList[i];
  }
  
  // create a block map based on the already declared point map
  // (used to determine NumMyElements and MyGlobalElements).
  // The same point map can be used for more block maps,
  // just change the input value of ElementSizeList
  Epetra_BlockMap BlockMap(NumGlobalElements,NumMyElements,
			   MyGlobalElements, 
			   ElementSizeList.Values(),0,Comm);

  // create a VBR matrix based on BlockMap
  Epetra_VbrMatrix A(Copy, BlockMap, 2);

  int MaxBlockSize = max_blk * max_blk * 100;
  
  int Indices[2];
  double* Values; Values = new double[MaxBlockSize];

  // cycle over all the local rows. 
  
  for (int i = 0; i < NumMyElements; ++i) 
  {
    // get GID of local row
    int GlobalNode = MyGlobalElements[i];
    // all lines but the last one will have to nonzero block-elements
    Indices[0] = GlobalNode;
    int NumEntries = 1;
    if (GlobalNode != NumGlobalElements - 1) 
    {
      Indices[1] = GlobalNode + 1;
      NumEntries++;
    }

    // with VBR matrices, we have to insert one block at time.
    // This required two more instructions, one to start this
    // process (BeginInsertGlobalValues), and another one to
    // commit the end of submissions (EndSubmitEntries).
    
    A.BeginInsertGlobalValues(GlobalNode, NumEntries, Indices);
    // insert diagonal
    int BlockRows = ElementSizeList[i];
    for (int k = 0 ; k < BlockRows * BlockRows ; ++k)
      Values[k] = 1.0 * i;
    A.SubmitBlockEntry(Values,BlockRows,BlockRows,BlockRows);

    // insert off diagonal if any
    if (GlobalNode != NumGlobalElements-1) 
    {
      int BlockCols = BlockRows + 1;
      for (int k = 0 ; k < BlockRows * BlockCols ; ++k)
	Values[k] = 1.0 * i;
      A.SubmitBlockEntry(Values,BlockRows,BlockRows,BlockCols);
    }
    
    A.EndSubmitEntries();
  }

  A.FillComplete();

  cout << A;

  delete[] Values;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);
}
