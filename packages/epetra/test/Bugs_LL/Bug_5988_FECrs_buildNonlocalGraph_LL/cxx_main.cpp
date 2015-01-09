#include <Epetra_FECrsGraph.h>
#include <Epetra_FECrsMatrix.h>
#include <Epetra_MpiComm.h>


int main(int argc, char**argv)
{
  MPI_Init(&argc,&argv);
  int rank; // My process ID

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  Epetra_MpiComm Comm( MPI_COMM_WORLD );

  int NumMyEquations = 2;
  Epetra_Map Map(-1LL, NumMyEquations, 0, Comm);

  Epetra_FECrsGraph graph(Copy, Map, 0, false, true);
  for (int i=0; i<NumMyEquations; ++i)
    {
      long long entry = i+rank*2;
      graph.InsertGlobalIndices(i+rank*2, 1, &entry);
    }
  {
    long long row = 0;
    long long col = 1;
    if (rank == 1)
      graph.InsertGlobalIndices(1, &row, 1, &col);
  }
  graph.GlobalAssemble();
  graph.FillComplete();
  graph.OptimizeStorage();

  Epetra_FECrsMatrix matrix(Copy, graph);
  {
    long long row = 0;
    long long col = 1;
    double value = 1;
    if (rank == 1)
      matrix.SumIntoGlobalValues(1, &row, 1, &col, &value);
  }
  matrix.GlobalAssemble();
  matrix.FillComplete();
  double norm = matrix.NormFrobenius();
  if (rank == 0) {
    std::cout << "Frobenius norm (should be 1 with 2 or more processors): "
              << norm << std::endl;
    if(norm != 1)
        std::cout << "tests FAILED" << std::endl;
  }
  MPI_Finalize();
}

