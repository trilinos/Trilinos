#include "ml_include.h"

#if defined(HAVE_ML_EPETRA)

#include "Epetra_Map.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Import.h"
#include "ml_RowMatrix.h"

//==============================================================================
ML_Epetra::RowMatrix::RowMatrix(ML_Operator* Op,
				Epetra_Comm& Comm) :
  Op_(0),
  Comm_(Comm),
  NumMyRows_(-1),
  NumGlobalRows_(-1),
  NumMyCols_(-1),
  NumGlobalCols_(-1),
  RowMap_(0),
  ColMap_(0),
  NumMyNonzeros_(0),
  Allocated_(128),
  NumGlobalNonzeros_(0),
  MaxNumEntries_(0),
  NormInf_(-1.0),
  NumMyDiagonals_(0),
  NumGlobalDiagonals_(0)
{
 
  Op_ = Op;

  NumMyRows_ = Op->outvec_leng;
  NumMyCols_ = Op->invec_leng;
 
  Comm_.SumAll(&NumMyRows_,&NumGlobalRows_,1);
  Comm_.SumAll(&NumMyCols_,&NumGlobalCols_,1);

  // right now only square matrices
  if (NumMyRows_ != NumMyCols_)
    ML_CHK_ERRV(-1); // FIXME

  if (NumGlobalRows_ != NumGlobalCols_)
    ML_CHK_ERRV(-2); // FIXME

  // create a row map based on linear distribution
  // I need this map because often codes use the RowMatrixRowMap.
  // Also, I need to check that the map of the input vector
  // and of the output vector are consistent with what I have here
  RowMap_ = new Epetra_Map(-1,NumMyRows_,0,Comm_);

  if (NumGlobalRows_ != RowMap_->NumMyElements())
    ML_CHK_ERRV(-3); // something went wrong

  // FIXME: use or delete ColMap!


  // need to perform a getrow() on all lines to setup some
  // parameters. This should not be too expensive

  Diagonal_.resize(NumMyRows());
  NumMyRowEntries_.resize(NumMyRows());

  Indices_.resize(Allocated_);
  Values_.resize(Allocated_);

  int ncnt; // ma che cavolo di nome e`???
  int MaxMyNumEntries;
  double MyNormInf = -1.0;

  for (int i = 0; i < NumMyCols() ; i++) {

    int ierr = ML_Operator_Getrow(Op_,1,&i,Allocated_,
				  &Indices_[0],&Values_[0],&ncnt);

    if (ierr == 0) {
      do {
	Allocated_ *= 2;
	Indices_.resize(Allocated_);
	Values_.resize(Allocated_);
	ierr = ML_Operator_Getrow(Op_,1,&i,Allocated_,
				  &Indices_[0],&Values_[0],&ncnt);
      } while (ierr == 0);
    }

    NumMyRowEntries_[i] = ncnt;
    NumMyNonzeros_ += ncnt;

    if (ncnt > MaxMyNumEntries) 
      MaxMyNumEntries = ncnt;
    
    // extract diagonal element and compute NormInf
    Diagonal_[i] = 0.0;
    double RowNormInf = 0.0;
    for (int j = 0 ; j < ncnt ; ++j) {
      MyNormInf += Values_[j];
      if (Indices_[j] == i) {
	Diagonal_[i] = Values_[j];
	++NumMyDiagonals_;
      }
    }

    if (RowNormInf > MyNormInf)
      MyNormInf = RowNormInf;

  } // for each col

  // fix a couple of global integers

  Comm_.SumAll(&NumMyNonzeros_,&NumGlobalNonzeros_,1);
  Comm_.SumAll(&NumMyDiagonals_,&NumGlobalDiagonals_,1);
  Comm_.MaxAll(&MaxMyNumEntries,&MaxNumEntries_,1);
  Comm_.MaxAll(&MyNormInf,&NormInf_,1);

  return;
}
      
//==============================================================================
ML_Epetra::RowMatrix::~RowMatrix()
{
  if (RowMap_) 
    delete RowMap_;

  return;

}

//==============================================================================
int ML_Epetra::RowMatrix::NumMyRowEntries(int MyRow, int & NumEntries)
{
  
  if (MyRow < 0 || MyRow >= NumMyRows())
    ML_CHK_ERR(-1); // out of range

  NumEntries = NumMyRowEntries_[MyRow];

  return(0);

}

//==============================================================================
int ML_Epetra::RowMatrix::MaxNumEntries()
{
  return(MaxNumEntries_);
}

//==============================================================================
int ML_Epetra::RowMatrix::
ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, 
		 double *Values, int * Indices)
{

  if (MyRow < 0 || MyRow >= NumMyRows())
    ML_CHK_ERR(-1); // not a local row

  if (NumMyRowEntries_[MyRow] > Length)
    ML_CHK_ERR(-2); // need more space


  int ierr = ML_Operator_Getrow(Op_,1,&MyRow,Length,
				Indices,Values,&NumEntries);

  ML_CHK_ERR(ierr); // if ierr is zero, the ML_Operator has changed

  if (NumEntries != NumMyRowEntries_[MyRow])
    ML_CHK_ERR(-4); // something went wrong

  return(0);

}

//==============================================================================
int ML_Epetra::RowMatrix::ExtractDiagonalCopy(Epetra_Vector & Diagonal)
{
  if (Diagonal.Map().SameAs(*RowMap_) != true)
    ML_CHK_ERR(-1); // incompatible maps

  for (int i = 0 ; i < NumMyRows() ; ++i)
    Diagonal[i] = Diagonal_[i];

  return(0);
}


//==============================================================================
// FIXME: some problems with ghost nodes here ??
int ML_Epetra::RowMatrix::Multiply(bool TransA, 
				   const Epetra_MultiVector& X, 
				   Epetra_MultiVector& Y)
{

  if (TransA == true)
    ML_CHK_ERR(-1); // not supported right now, easy fix

  if (X.Map().SameAs(*RowMap_) != true)
    ML_CHK_ERR(-2); // incompatible maps
  
  if (Y.Map().SameAs(*RowMap_) != true)
    ML_CHK_ERR(-3); // incompatible maps

  int ierr;

  ierr = ML_Operator_Apply(Op_,X.MyLength(),X.Values(),
			   Y.MyLength(), Y.Values());

  ML_RETURN(ierr);
}


//==============================================================================
double ML_Epetra::RowMatrix::NormInf()
{

 
  return(NormInf_);
}


//==============================================================================
int ML_Epetra::RowMatrix::NumGlobalNonzeros()
{
  return(NumGlobalNonzeros_);
}

//==============================================================================
int ML_Epetra::RowMatrix::NumGlobalRows()
{
  return(NumGlobalRows_);
}


//==============================================================================
int ML_Epetra::RowMatrix::NumGlobalCols()
{
  return(NumGlobalCols_);
}

//==============================================================================
int ML_Epetra::RowMatrix::NumGlobalDiagonals()
{
  return(NumGlobalDiagonals_);
}

//==============================================================================
int ML_Epetra::RowMatrix::NumMyNonzeros()
{
  return(NumMyNonzeros_);
}

//==============================================================================
int ML_Epetra::RowMatrix::NumMyRows() 
{
  return(NumMyRows_);
}

//==============================================================================
int ML_Epetra::RowMatrix::NumMyCols() {

  return(NumMyCols_);
}

//==============================================================================
int ML_Epetra::RowMatrix::NumMyDiagonals()
{
  return(NumMyDiagonals_);
}

//==============================================================================
bool ML_Epetra::RowMatrix::LowerTriangular()
{
  // FIXME ???
  return(false);
}

//==============================================================================
bool ML_Epetra::RowMatrix::UpperTriangular()
{
  // FIXME ???

  return(false);

}

//==============================================================================
const Epetra_Map& ML_Epetra::RowMatrix::RowMatrixRowMap(){

  return(*RowMap_);
}

//==============================================================================
const Epetra_Map& ML_Epetra::RowMatrix::RowMatrixColMap()
{
  ML_EXIT(-1); // still to fix
  return(*ColMap_);
}

//==============================================================================
const Epetra_Import* ML_Epetra::RowMatrix::RowMatrixImporter()
{

  return(Importer_);
}

#endif /* HAVE_ML_EPETRA */

