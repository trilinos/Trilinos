#ifndef IFPACK_BLOCKGAUSSSEIDEL_H
#define IFPACK_BLOCKGAUSSSEIDEL_H

#include "Ifpack_ConfigDefs.h"
#ifdef HAVE_IFPACK_TEUCHOS
#include "Ifpack_BlockPreconditioner.h" 
#include "Ifpack_Utils.h"

//! Ifpack_BlockGaussSeidel: a class to define block-GaussSeidel preconditioners preconditioners of Epetra_RowMatrix's.

/*!
  The Ifpack_BlockGaussSeidel class enables the construction of block-GaussSeidel 
  preconditioners.
*/

template<class T>
class Ifpack_BlockGaussSeidel : public Ifpack_BlockPreconditioner<T> {

public:

  //! Constructor.
  Ifpack_BlockGaussSeidel(Epetra_RowMatrix* Matrix) :
    Ifpack_BlockPreconditioner<T>(Matrix)
  {
  };
  
  //! Destructor.
  virtual ~Ifpack_BlockGaussSeidel()
  {};

  //! Applies the block GaussSeidel preconditioner to X, returns the result in Y.
  /*! 
    \param In
    X - A Epetra_MultiVector of dimension NumVectors to be preconditioned.
    \param Out
    Y -A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.

    \warning I do consider that X and Y are two distinct objects.
    This means that this method will not work with AztecOO; the
    class must be used through Ifpack_AdditiveSchwarz only.
    */
  int ApplyInverse(const Epetra_MultiVector& X, 
		   Epetra_MultiVector& Y) const
  {

    if (IsComputed() == false)
      IFPACK_CHK_ERR(-4); // need to compute the prec first

    if (OverlapLevel_ != 0)
      IFPACK_CHK_ERR(-1); // only zero overlap is supported

    if (X.NumVectors() != Y.NumVectors())
      IFPACK_CHK_ERR(-1); // not valid

    if (ZeroStartingSolution_)
      Y.PutScalar(0.0);

    // ------------ //
    // single sweep //
    // ------------ //

    if (NumSweeps() == 1 && ZeroStartingSolution_
	&& (PrintFrequency() != 0)) {
      IFPACK_CHK_ERR(ApplyBGS(X,Y));
      return(0);
    }

    // --------------------- //
    // general case (solver) //
    // --------------------- //

    if (PrintFrequency())
      Ifpack_PrintResidual(Label(),Matrix(),Y,X);

    // starting solution
    Epetra_MultiVector AX(X);
    Epetra_MultiVector Ynew(Y);

    for (int j = 0; j < NumSweeps() ; j++) {

      // compute the residual, I can skip first iteration
      // if starting solution is zero
      if (j || !ZeroStartingSolution_) {
	IFPACK_CHK_ERR(Apply(Y,AX));
	AX.Update(1.0,X,-1.0);
      }

      // apply the block diagonal of A and update the residual
      ApplyBGS(AX,Ynew);

      Y.Update(DampingFactor(),Ynew,1.0);

      if (PrintFrequency() && (j != 0) && (j % PrintFrequency() == 0))
	Ifpack_PrintResidual(j,Matrix(),Y,X);

    }

    if (PrintFrequency())
      Ifpack_PrintResidual(NumSweeps(),Matrix(),Y,X);

    return(0);

  }

  //! Applies one sweep of block Jacobi.
  virtual int ApplyBGS(const Epetra_MultiVector& X, 
		       Epetra_MultiVector& Y) const
  {

    // cycle over all local subdomains

    int Length = Matrix().MaxNumEntries();
    vector<int> Indices;
    vector<double> Values;
    Indices.resize(Length);
    Values.resize(Length);

    int NumMyRows = Matrix().NumMyRows();

    Y.PutScalar(0.0);

    for (int i = 0 ; i < NumLocalBlocks() ; ++i) {

      int LID, GID;

      // update from previous block

      for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
	LID = Containers_[i]->ID(j);

	int NumEntries;
	IFPACK_CHK_ERR(Matrix().ExtractMyRowCopy(LID, Length,NumEntries,
						 &Values[0], &Indices[0]));

	for (int k = 0 ; k < NumEntries ; ++k) {
	  int col = Indices[k];
	  if (col >= NumMyRows ) 
	    continue;

	  if ((*Partitioner_)(col) < i) {
	    for (int kk = 0 ; kk < Y.NumVectors() ; ++kk) {
	      X[kk][LID] -= Values[k] * Y[kk][col];
	    }
	  }
	}
      }

      // solve with this block

      for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
	LID = Containers_[i]->ID(j);
	for (int k = 0 ; k < Y.NumVectors() ; ++k) {
	  Containers_[i]->RHS(j,k) = X[k][LID];
	}
      }

      IFPACK_CHK_ERR(Containers_[i]->ApplyInverse());

      for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
	LID = Containers_[i]->ID(j);
	for (int k = 0 ; k < Y.NumVectors() ; ++k) {
	  Y[k][LID] += Containers_[i]->LHS(j,k);
	}
      }
    }

    return(0);
  }

private:

//! Sets the label of \c this object.
virtual int SetLabel()
{
  Label_ = "Ifpack_BlockGaussSeidel, # blocks = "
    + Ifpack_toString(NumLocalBlocks());
  return(0);
}

};

#endif // HAVE_IFPACK_TEUCHOS
#endif // class Ifpack_BlockGaussSeidel
