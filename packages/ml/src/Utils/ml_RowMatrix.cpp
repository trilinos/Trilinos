/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */
#include "ml_include.h"

#if defined(HAVE_ML_EPETRA)

#include "Epetra_Map.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_Import.h"
#include "Epetra_MultiVector.h"
#include "ml_RowMatrix.h"
#include "ml_lapack.h"

//==============================================================================
ML_Epetra::RowMatrix::RowMatrix(ML_Operator* Op,
                                const Epetra_Comm* UserComm,
                                const bool cheap, 
                                const USR_COMM comm) :
  Op_(0),
  FreeCommObject_(false),
  NumMyRows_(-1),
  NumGlobalRows_(-1),
  NumMyCols_(-1),
  NumGlobalCols_(-1),
  RangeMap_(0),
  ColMap_(0),
  MaxNumEntries_(0),
  Allocated_(128),
  NormInf_(-1.0),
  NumMyNonzeros_(0),
  NumGlobalNonzeros_(0),
  NumMyDiagonals_(0),
  NumGlobalDiagonals_(0),
  Importer_(0)
{
 
  if (UserComm) {
    // simply stick it to the pointer if the user passed Comm
    Comm_ = UserComm;
  }
  else {
    // otherwise we create a default communicator object,
    // and track down the need of freeing it in the dtor
#ifdef HAVE_MPI
    Comm_ = new Epetra_MpiComm(comm);
#else
    Comm_ = new Epetra_SerialComm;
#endif
    FreeCommObject_ = true;
  }

  Op_ = Op;

  Label_ = new char[80];
  sprintf(Label_,"%s","ML_Epetra::RowMatrix");

  NumMyRows_ = Op->outvec_leng; 
  NumMyCols_ = Op->invec_leng;  // this is fixed at the end of this method 
 
  // create a row map based on linear distribution
  // I need this map because often codes use the RowMatrixRowMap.
  // Also, I need to check that the map of the input vector
  // and of the output vector are consistent with what I have here
  RangeMap_ = new Epetra_Map(-1,NumMyRows_,0,Comm());
//  if (NumMyCols_ == NumMyRows_) // FIXME: not necessarily true for global values
//    DomainMap_ = RangeMap_;
//  else
      DomainMap_ = new Epetra_Map(-1,NumMyCols_,0,Comm());
  NumGlobalRows_ = RangeMap_->NumGlobalElements();
  NumGlobalCols_ = DomainMap_->NumGlobalElements();

  // need to perform a getrow() on all lines to setup some
  // parameters. Note that NormInf is no longer computed.

  int MaxMyNumEntries;

  if (cheap) {

    NumMyNonzeros_ = Op_->N_nonzeros;
    MaxMyNumEntries = Op_->max_nz_per_row;
    Allocated_ = MaxMyNumEntries;
    Indices_.resize(Allocated_);
    Values_.resize(Allocated_);
    NumMyDiagonals_ = NumMyRows_;
    Diagonal_.resize(0);

  }
  else {
    NumMyRowEntries_.resize(NumMyRows());

    Diagonal_.resize(NumMyRows()); 
    Indices_.resize(Allocated_);
    Values_.resize(Allocated_);

    int ncnt; // ma che cavolo di nome e`???
    MaxMyNumEntries = 0;
    double MyNormInf = -1.0;

    // Note: I do not compute norm inf any more, it was quite
    // expensive.

    for (int i = 0; i < NumMyRows() ; ++i) {

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

      double tmp = 0.0;

      // extract diagonal element
      // NOTE: the diagonal value is set to zero if not present
      Diagonal_[i] = 0.0;
      for (int j = 0 ; j < ncnt ; ++j) {
        if (Indices_[j] == i) {
          Diagonal_[i] = Values_[j];
          ++NumMyDiagonals_;
        }
        tmp += fabs(Values_[j]);
      }

      if (tmp > MyNormInf) MyNormInf = tmp;

    } // for each row

    Comm().MaxAll(&MyNormInf, &NormInf_, 1);
  }

  // fix a couple of global integers

//TODO: CJ
  long long NumMyNonzeros_tmp = NumMyNonzeros_;
  Comm().SumAll(&NumMyNonzeros_tmp,&NumGlobalNonzeros_,1);
  long long NumMyDiagonals_tmp = NumMyDiagonals_;
  Comm().SumAll(&NumMyDiagonals_tmp,&NumGlobalDiagonals_,1);
  Comm().MaxAll(&MaxMyNumEntries,&MaxNumEntries_,1);

  // build a list of global indices for columns

  if (Comm().NumProc() == 1) {
    ColMap_ = DomainMap_; 
  }
  else {

    int Nghost;
    int Ncols;
    int Ncols_offset;

    if (Op->getrow->pre_comm == NULL)
      Nghost = 0;
    else {
      if (Op->getrow->pre_comm->total_rcv_length <= 0)
        ML_CommInfoOP_Compute_TotalRcvLength(Op->getrow->pre_comm);
      Nghost = Op->getrow->pre_comm->total_rcv_length;
    }

    Ncols = Op->invec_leng;

    Comm().ScanSum(&Ncols,&Ncols_offset,1); 
    Ncols_offset -= Ncols;

    std::vector<double> global_col_id; 
    global_col_id.resize(Ncols + Nghost + 1);

    std::vector<int> global_col_id_as_int; 
    global_col_id_as_int.resize(Ncols + Nghost + 1);

    for (int i = 0 ; i < Ncols ; ++i) {
      global_col_id[i] = (double) (Ncols_offset + i);
      global_col_id_as_int[i] = Ncols_offset + i;
    }

    for (int i = 0 ; i < Nghost; i++) 
      global_col_id[i + Ncols] = -1;

    ML_exchange_bdry(&global_col_id[0],Op->getrow->pre_comm,
                     Op->invec_leng,Op->comm,ML_OVERWRITE,NULL);

    for (int j = 0; j < Ncols + Nghost; ++j) {
      global_col_id_as_int[j] = (int) global_col_id[j];
    }

    // create the column map

    ColMap_ = new Epetra_Map(-1,Ncols + Nghost,
                             &global_col_id_as_int[0],0,Comm());
    NumMyCols_ = ColMap_->NumMyElements();
  }

  return;
}
      
//==============================================================================
ML_Epetra::RowMatrix::~RowMatrix()
{
  if (ColMap_)
    if (ColMap_ != DomainMap_) {
      delete ColMap_;
      ColMap_ = 0;
    }

  // the one that is always allocated is RangeMap_;
  if (DomainMap_) 
    if (DomainMap_ != RangeMap_) {
      delete DomainMap_;
      DomainMap_ = 0;
    }

  if (RangeMap_) {
    delete RangeMap_;
    RangeMap_ = 0;
  }

  if (Label_)
    delete [] Label_;

  if (FreeCommObject_)
    delete Comm_;

  if (Importer_)
    delete Importer_;
  
  return;
}

//==============================================================================
int ML_Epetra::RowMatrix::NumMyRowEntries(int MyRow, int & NumEntries) const
{
  if (NumMyRowEntries_.size() == 0)
    ML_CHK_ERR(-2); // prec was built as "cheap"
    
  if (MyRow < 0 || MyRow >= NumMyRows())
    ML_CHK_ERR(-1); // out of range

  NumEntries = NumMyRowEntries_[MyRow];

  return(0);
}

//==============================================================================
int ML_Epetra::RowMatrix::MaxNumEntries() const
{
  return(MaxNumEntries_);
}

//==============================================================================
int ML_Epetra::RowMatrix::
ExtractMyRowCopy(int MyRow, int Length, int & NumEntries, 
		 double *Values, int * Indices) const
{

  if (MyRow < 0 || MyRow >= NumMyRows())
    ML_CHK_ERR(-1); // not a local row

#if 0
  if (NumMyRowEntries_[MyRow] > Length) {
    cerr << MyRow << " " << NumMyRowEntries_[MyRow] << " " << Length << endl;
    ML_CHK_ERR(-2); // need more space
  }
#endif

  int ierr = ML_Operator_Getrow(Op_,1,&MyRow,Length,
				Indices,Values,&NumEntries);

  if (ierr < 0)
    ML_CHK_ERR(ierr);
  
#if 0
  if (NumEntries != NumMyRowEntries_[MyRow])
    ML_CHK_ERR(-4); // something went wrong
#endif

  return(0);

}

//==============================================================================
int ML_Epetra::RowMatrix::ExtractDiagonalCopy(Epetra_Vector & Diagonal) const
{
  if (Diagonal.Map().SameAs(*RangeMap_) != true)
    ML_CHK_ERR(-1); // incompatible maps

  if (Diagonal_.size() == 0) {
    int size = NumMyRows();
    int incr = 1;
    double* ptr;
    ML_Operator_Get_Diag(Op_, NumMyRows(), &ptr);

    DCOPY_F77(&size, ptr, &incr, Diagonal.Values(), &incr);
  }
  else {
    for (int i = 0 ; i < NumMyRows() ; ++i)
      Diagonal[i] = Diagonal_[i];
  }

  return(0);
}


//==============================================================================
// FIXME: some problems with ghost nodes here ??
int ML_Epetra::RowMatrix::Multiply(bool TransA, 
				   const Epetra_MultiVector& X, 
				   Epetra_MultiVector& Y) const
{

  if (TransA == true)
    ML_CHK_ERR(-1); // not supported right now, easy fix

  if (X.Map().SameAs(*DomainMap_) != true)
    ML_CHK_ERR(-2); // incompatible maps
  
  if (Y.Map().SameAs(*RangeMap_) != true)
    ML_CHK_ERR(-3); // incompatible maps

  int ierr;

  for (int i = 0 ; i < X.NumVectors() ; ++i) {
    ierr = ML_Operator_Apply(Op_,X.MyLength(),X[i],
			     Y.MyLength(), Y[i]);
    ML_CHK_ERR(ierr);
  }

  return(0);
}


//==============================================================================
double ML_Epetra::RowMatrix::NormInf() const
{
  return(NormInf_);
}

//TODO: CJ correct the int returning functions.
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
//==============================================================================
int ML_Epetra::RowMatrix::NumGlobalNonzeros() const
{
  return(NumGlobalNonzeros_);
}

//==============================================================================
int ML_Epetra::RowMatrix::NumGlobalRows() const
{
  return(NumGlobalRows_);
}

//==============================================================================
int ML_Epetra::RowMatrix::NumGlobalCols() const
{
  return(NumGlobalCols_);
}

//==============================================================================
int ML_Epetra::RowMatrix::NumGlobalDiagonals() const
{
  return(NumGlobalDiagonals_);
}
#endif

//==============================================================================
long long ML_Epetra::RowMatrix::NumGlobalNonzeros64() const
{
  return(NumGlobalNonzeros_);
}

//==============================================================================
long long ML_Epetra::RowMatrix::NumGlobalRows64() const
{
  return(NumGlobalRows_);
}

//==============================================================================
long long ML_Epetra::RowMatrix::NumGlobalCols64() const
{
  return(NumGlobalCols_);
}

//==============================================================================
long long ML_Epetra::RowMatrix::NumGlobalDiagonals64() const
{
  return(NumGlobalDiagonals_);
}

//==============================================================================
int ML_Epetra::RowMatrix::NumMyNonzeros() const
{
  return(NumMyNonzeros_);
}

//==============================================================================
int ML_Epetra::RowMatrix::NumMyRows() const
{
  return(NumMyRows_);
}

//==============================================================================
int ML_Epetra::RowMatrix::NumMyCols() const
{
  return(NumMyCols_);
}

//==============================================================================
int ML_Epetra::RowMatrix::NumMyDiagonals() const
{
  return(NumMyDiagonals_);
}

//==============================================================================
bool ML_Epetra::RowMatrix::LowerTriangular() const
{
  // FIXME ???
  return(false);
}

//==============================================================================
bool ML_Epetra::RowMatrix::UpperTriangular() const
{
  // FIXME ???

  return(false);

}

//==============================================================================
const Epetra_Map& ML_Epetra::RowMatrix::RowMatrixRowMap() const
{
  return(*RangeMap_);
}

//==============================================================================
const Epetra_Map& ML_Epetra::RowMatrix::RowMatrixColMap() const
{
  return(*ColMap_);
}

//==============================================================================
const Epetra_Import* ML_Epetra::RowMatrix::RowMatrixImporter() const
{
  if (!Importer_)
    Importer_ = new Epetra_Import(RowMatrixColMap(),RowMatrixRowMap());
  return(Importer_);
}

//==============================================================================
int ML_Epetra::RowMatrix::Print() const
{
   ML_Operator_Print_UsingGlobalOrdering(Op_,Label_,0,0);
   return(0);
}

#endif /* HAVE_ML_EPETRA */

