#ifndef _Epetra_ESI_Operator_h_
#define _Epetra_ESI_Operator_h_

#include "Epetra_Object.h"
#include "Epetra_RowMatrix.h"

//
//Epetra_ESI_Vector.h includes everything "above" it,
//including Epetra_ESI_Object.h, Epetra_Comm.h, Epetra_Map.h,
//Epetra_ESI_Vector.h, Epetra_Vector.h, and the ESI headers.
//
#include "Epetra_ESI_Vector.h"


namespace epetra_esi {

/**
This class implements these ESI interfaces:
<ul>
<li> esi::Object
<li> esi::Operator
<li> esi::OperatorTranspose
</ul>
Note that although this class is templated, it may only be instantiated on
the type-pair double,int.
*/

template<class Scalar, class Ordinal>
class Operator : public virtual esi::Object,
                 public virtual epetra_esi::Object,
                 public virtual esi::Operator<Scalar, Ordinal>,
                 public virtual esi::OperatorTranspose<Scalar, Ordinal>,
                 public virtual Epetra_RowMatrix
{
 public:
  /** Constructor. */
  Operator(esi::Operator<Scalar, Ordinal>& esi_op);

  /** Destructor. */
  virtual ~Operator() { delete petra_row_map_; delete petra_import_map_;
                                  delete petra_comm_; }


  //
  //esi::Operator functions.
  //

  /** Function for signalling that the initialization and data-loading is
    complete.
   @return error-code 0 if successful.
  */
  virtual esi::ErrorCode setup( void )
    { setupHasBeenCalled_ = true; return(esi_op_.setup()); }


  /** Function for applying this operator to a vector. (This is a matvec.) If
   the setup() function has not yet been called, it will be called internally,
   automatically at this point.
   @param x Input.
   @param y On exit, the contents of this will be the result of the matvec.
   @return error-code 0 if successful.
  */
  virtual esi::ErrorCode apply(esi::Vector<Scalar, Ordinal>& x, esi::Vector<Scalar, Ordinal>& y)
    { int err = 0;
      if (!setupHasBeenCalled_) err = setup();
      if (err) return(-1);
      err = esi_op_.apply(x, y);
      return(err);
    }


  //
  //esi::OperatorTranspose
  //

  /** Function for applying the transpose of this operator to a vector.
    (Transpose matvec.)
   @param x Input.
   @param y On exit, the contents of this will be the result of the matvec.
   @return error-code 0 if successful.
  */
  virtual esi::ErrorCode applyTranspose(esi::Vector<Scalar, Ordinal>& x, esi::Vector<Scalar, Ordinal>& y)
    { if (!haveTranspose_) return(-1);
      int err = 0;
      if (!setupHasBeenCalled_) err = setup();
      if (err) return(-1);
      return( esi_op_trans_->applyTranspose(x, y) ); }


  //
  //Epetra_RowMatrix functions
  //

  bool Filled() const { return( true ); }

  int ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, double *Values, int * Indices) const
    { if (!haveRowReadAccess_) return(-1);
      int GlobalRow = MyRow + localOffset_;
//cout << "ExtractMyRowCopy MyRow " << MyRow<< ", GlobalRow " << GlobalRow<<endl;
      return( esi_row_read_->copyOutRow(GlobalRow, Values, Indices, Length,
                                        NumEntries) );
    }

  int ExtractDiagonalCopy(Epetra_Vector & Diagonal) const
    { //cout << "ExtractDiagonalCopy" << endl;
     return(-1); }

  int Multiply(bool TransA, const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
    {
       const Epetra_Vector* px = dynamic_cast<const Epetra_Vector*>(&X);
       Epetra_Vector* py = dynamic_cast<Epetra_Vector*>(&Y);
       if (px == NULL || py == NULL) {
         cout << "epetra_esi::Operator::Multiply ERROR, input multi-vectors must"
           << " be castable to type Epetra_Vector." << endl;
         return(-1);
       }

       //ok, we have Epetra_Vectors, let's create epetra_esi::Vector to use
       //with the esi::Operator's apply function.
       epetra_esi::Vector<double,int>* pesi_x =
             new epetra_esi::Vector<double,int>(*px);
       epetra_esi::Vector<double,int>* pesi_y = 
             new epetra_esi::Vector<double,int>(*py);
       if (pesi_x == NULL || pesi_y == NULL) {
         cout << "epetra_esi::Operator::Multiply ERROR, failed to create epetra_esi::Vector." << endl;
         return(-1);
       }

       int err = -1;
       if (!TransA) {
         err = esi_op_.apply(*pesi_x, *pesi_y);
       }
       else {
         if (haveTranspose_) {
           err = esi_op_trans_->applyTranspose(*pesi_x, *pesi_y);
         }
       }

       delete pesi_x;
       delete pesi_y;

       return(err);
    }

  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
    {
       cout << "epetra_esi::Operator::Apply not implemented." << endl;
       return(-1);
    }

  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
    {
       cout << "epetra_esi::Operator::ApplyInverse not implemented." << endl;
       return(-1);
    }

  char * Label() const
    {
      static char label[] = "epetra_esi::Operator";
      return(label);
    }

  bool UseTranspose() const {return(useTranspose_);}

  int SetUseTranspose(bool UseTranspose)
    { useTranspose_ = UseTranspose; return(0); }

  int Solve(bool Upper, bool Trans, bool UnitDiagonal,
            const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
    { //cout << "Solve" << endl;
      return(-1); }

  int InvRowSums(Epetra_Vector& x) const
    { //cout << "InvRowSums" << endl; 
      return(-1); }

  int LeftScale(const Epetra_Vector& x)
    { //cout << "LeftScale" << endl; 
      return(-1); }

  int InvColSums(Epetra_Vector& x) const
    { //cout << "InvColSums" << endl; 
      return(-1); }

  int RightScale(const Epetra_Vector& x)
    { //cout << "RightScale" << endl; 
      return(-1); }

  bool HasNormInf() const { return(false); }

  double NormInf() const
    { //cout << "NormInf" << endl; 
      return(100.0); } // !!!! Garbage return value !!!!

  double NormOne() const
    { //cout << "NormOne" << endl; 
      return(-1); }

  int NumGlobalNonzeros() const
    { //cout << "NumGlobalNonzeros" << endl; 
      return(-1); }

  int NumGlobalRows() const
    { //cout << "NumGlobalRows" << endl; 
      return(numGlobalRows_); }

  int NumGlobalCols() const
    { //cout << "NumGlobalCols" << endl; 
      return(numGlobalCols_); }

  int NumGlobalDiagonals() const
    { //cout << "NumGlobalDiagonals" << endl; 
      return(-1); }

  int NumMyNonzeros() const
    { //cout << "NumMyNonzeros" << endl; 
      return(-1); }

  int NumMyRowEntries(int MyRow, int & NumEntries) const
    { //cout << "NumMyRowEntries" << endl; 
      return(-1); }

  int NumMyRows() const
    { 
      //AztecOO assumes ExtractMyRowCopy will work if NumMyCols and NumMyRows
      //both work. So here we only return a valid value if haveRowReadAccess_
      //is true.
      if (!haveRowReadAccess_) 
      return(-1);
      //cout << "NumMyRows " << numLocalRows_ << endl; 
      return(numLocalRows_);
    }

  int NumMyCols() const
    {
      //AztecOO assumes ExtractMyRowCopy will work if NumMyCols and NumMyRows
      //both work. So here we only return a valid value if haveRowReadAccess_
      //is true.
      if (!haveRowReadAccess_) return(-1);
      //cout << "NumMyCols " << numLocalCols_ << endl; 
      return(numLocalCols_);
    }

  int NumMyDiagonals() const
    { //cout << "NumMyDiagonals" << endl; 
      return(-1); }

  bool LowerTriangular() const
    { //cout << "LowerTriangular" << endl; 
      return(false); }

  bool UpperTriangular() const
    { //cout << "UpperTriangular" << endl; 
      return(false); }

  const Epetra_Comm & Comm() const
    { //cout << "Comm" << endl; 
      return(*petra_comm_); }

  const Epetra_BlockMap & DomainMap() const
    { //cout << "DomainMap" << endl;
      return(*petra_domain_map_); }

  const Epetra_BlockMap & RangeMap() const
    { //cout << "RangeMap" << endl;
      return(*petra_range_map_); }

  const Epetra_BlockMap & BlockRowMap() const
    { //cout << "RowMap" << endl; 
      return(*petra_row_map_); }

  const Epetra_BlockMap & BlockImportMap() const
    { //cout << "ImportMap" << endl; 
      return(*petra_import_map_); }

  const Epetra_Import * Importer() const
    { //cout << "Importer" << endl; 
      return(petra_import_); }

 private:
  esi::Operator<Scalar, Ordinal>& esi_op_;
  esi::OperatorTranspose<Scalar, Ordinal>* esi_op_trans_;
  esi::MatrixData<Ordinal>* esi_matrix_data_;
  esi::MatrixRowReadAccess<Scalar,Ordinal>* esi_row_read_;
  bool setupHasBeenCalled_;
  bool haveTranspose_;
  bool haveMatrixData_;
  bool haveRowReadAccess_;
  bool useTranspose_;
  Epetra_Comm* petra_comm_;
  Epetra_Map* petra_domain_map_, *petra_range_map_;
  Epetra_Map* petra_row_map_, *petra_import_map_;
  Epetra_Import* petra_import_;
  Ordinal numGlobalRows_, numGlobalCols_, numLocalRows_, numLocalCols_;
  Ordinal localOffset_;
};

}; //namespace epetra_esi

#ifdef EPETRA_ESI_INCLUDE_IMPLEMENTATION
#include "Epetra_ESI_Operator.cpp"
#endif

#endif

