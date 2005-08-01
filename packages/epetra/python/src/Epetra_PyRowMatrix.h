#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Import.h"

class PyRowMatrix: public Epetra_RowMatrix 
{
  public:
    PyRowMatrix(const Epetra_Comm& Comm) :
      Comm_(Comm)
    {}

    virtual ~PyRowMatrix() {}

    virtual int MaxNumEntries() const = 0;

    virtual int ExtractDiagonalCopy(Epetra_Vector & Diagonal) const = 0;

    virtual int Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const  = 0;

    virtual int Solve(bool Upper, bool Trans, bool UnitDiagonal, const Epetra_MultiVector& X, 
                      Epetra_MultiVector& Y) const = 0;

    virtual int InvRowSums(Epetra_Vector& x) const = 0;

    virtual int LeftScale(const Epetra_Vector& x)  = 0;

    virtual int InvColSums(Epetra_Vector& x) const = 0;

    virtual int RightScale(const Epetra_Vector& x) = 0;

    virtual bool Filled() const  = 0;

    virtual double NormInf() const = 0;

    virtual double NormOne() const = 0;

    virtual int NumGlobalNonzeros() const = 0;

    virtual int NumGlobalRows() const = 0;

    virtual int NumGlobalCols() const = 0;

    virtual int NumGlobalDiagonals() const = 0;

    virtual int NumMyNonzeros() const = 0;

    virtual int NumMyRows() const = 0;

    virtual int NumMyCols() const = 0;

    virtual int NumMyDiagonals() const = 0;

    virtual bool LowerTriangular() const = 0;

    virtual bool UpperTriangular() const = 0;

    virtual const Epetra_Map & RowMatrixRowMap() const = 0;

    virtual const Epetra_Map & RowMatrixColMap() const = 0;

    virtual const Epetra_Import * RowMatrixImporter() const = 0;

    virtual int NumMyRowEntries(int MyRow, int & NumEntries) const = 0;

    virtual int ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, double *Values, int * Indices) const  = 0;

    virtual const Epetra_Map & Map() const = 0;

    virtual const Epetra_Map & OperatorRangeMap() const = 0;

    virtual const Epetra_Map & OperatorDomainMap() const = 0;

    virtual const char * Label() const
    {
      return("PyRowMatrix");
    }

    virtual bool UseTranspose() const = 0;

    virtual bool HasNormInf() const = 0;

    virtual const Epetra_Comm& Comm() const
    {
      return(Comm_);
    }

    int SetUseTranspose(bool UseTranspose) = 0;

    int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const = 0;

    int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const = 0;

  private:
    const Epetra_Comm& Comm_;
}; 

void SetDouble(double* array, int pos, double val)
{
  array[pos] = val;
}

double GetDouble(double* array, int pos)
{
  return(array[pos]);
}

void SetInt(int* array, int pos, int val)
{
  array[pos] = val;
}

int GetInt(int* array, int pos)
{
  return(array[pos]);
}

int* NewInt(int size)
{
  return new int[size];
}

double* NewDouble(int size)
{
  return new double[size];
}

void DeleteInt(int* array)
{
  delete[] array;
}

void DeleteDouble(double* array)
{
  delete[] array;
}
