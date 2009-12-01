#ifndef TPICRSMATRIX_H
#define TPICRSMATRIX_H

#include <Epetra_BasicRowMatrix.h>
#include <Epetra_SerialDenseVector.h>
#include <Epetra_IntSerialDenseVector.h>
#include <TPI.h>

class TPICrsMatrix : public virtual Epetra_BasicRowMatrix {
  public:

    //! Constructor
    TPICrsMatrix(const Epetra_Map& RowMap, const Epetra_Map& ColMap, const int* NumEntriesPerRow);

    //! Destructor
    virtual ~TPICrsMatrix();

    //! Insert entries into the matrix using local indices
    int InsertMyValues(int LocalRow, int NumEntries, const double *values, const int *indices);

    //! Indicate that we are finished inserting entries.
    int FillComplete(const Epetra_Map& DomainMap, const Epetra_Map& RangeMap);

    //! @name Implementing Epetra_BasicRowMatrix
    //@{ 

    //! Query the matrix for the number of entries on a particular row.
    int NumMyRowEntries(int MyRow, int & NumEntries) const;

    //! Extract a copy of a particular row.
    int ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, double *Values, int * Indices) const;

    //! Extract a view of a particular row.
    int ExtractMyEntryView(int CurEntry, double * & Value, int & RowIndex, int & ColIndex);

    //! Extract a const view of a particular row.
    int ExtractMyEntryView(int CurEntry, double const * & Value, int & RowIndex, int & ColIndex) const;

    //@}

    //! @name Over-ridden from Epetra_BasicRowMatrix
    //@{ 

    //! Transpose not currently supported by this class; will return -1 and have no effect.
    int SetUseTranspose(bool use_transpose);

    //! Multiply() method over-ridden to add support for TPI. Error if \c TransA is \c true.
    int Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    //! Multiply() method over-ridden to add support for TPI. Error if \c TransA is \c true.
    int Multiply(bool TransA, double alpha, const Epetra_MultiVector& X, double beta, Epetra_MultiVector& Y) const;

    //@}

  private:
    int RNNZ(int LocalRow) const;
    int NumMyRows_, NumMyEntries_;
    Epetra_IntSerialDenseVector IndexOffset_;
    Epetra_IntSerialDenseVector Indices_;
    Epetra_SerialDenseVector Values_;

    struct tpi_crs_matrix {
      int      nRow ;
      const int    * A_pc ;
      const int    * A_ia ;
      const double * A_a ;
      int      numVectors;
      double alpha, beta;
      const double ** x ;
      double ** y ;
    };

    static void tpi_work_crs_matrix_apply( TPI_Work * work );

    static void tpi_crs_matrix_apply(
        const int      nRow ,
        const int    * A_pc ,
        const int    * A_ia ,
        const double * A_a ,
        const int      numVectors, 
        double alpha, 
        const double ** x ,
        double beta,
              double ** y );
};

#endif
