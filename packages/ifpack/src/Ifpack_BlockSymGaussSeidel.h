#ifndef IFPACK_BLOCKSYMGAUSSSEIDEL_H
#define IFPACK_BLOCKSYMGAUSSSEIDEL_H

#include "Ifpack_ConfigDefs.h"
#ifdef HAVE_IFPACK_TEUCHOS
#include "Ifpack_BlockPreconditioner.h" 
#include "Ifpack_Utils.h"

//! Ifpack_BlockSymGaussSeidel: a class to define block-GaussSeidel preconditioners preconditioners of Epetra_RowMatrix's.

/*!
  The Ifpack_BlockSymGaussSeidel class enables the construction of block-GaussSeidel 
  preconditioners.
*/

template<class T>
class Ifpack_BlockSymGaussSeidel : public Ifpack_BlockPreconditioner<T> {

public:

  //! Constructor.
  Ifpack_BlockSymGaussSeidel(Epetra_RowMatrix* Matrix) :
    Ifpack_BlockPreconditioner<T>(Matrix)
  {
    KeepNonFactoredMatrix_ = true;
  };
  
  //! Destructor.
  virtual ~Ifpack_BlockSymGaussSeidel()
  {};

  //! Applies the block GaussSeidel preconditioner to X, returns the result in Y.
  /*! 
    \param In
    X - A Epetra_MultiVector of dimension NumVectors to be preconditioned.
    \param Out
    Y -A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.

    \warning In order to work with AztecOO, any implementation of this method 
    must support the case where X and Y are the same object.
    */
  int ApplyInverse(const Epetra_MultiVector& X, 
		   Epetra_MultiVector& Y) const
  {

    if (OverlapLevel_ != 0)
      IFPACK_CHK_ERR(-1); // required zero-overlap mode

    if (IsComputed() == false)
      IFPACK_CHK_ERR(-4); // need to compute the prec first

    if (X.NumVectors() != Y.NumVectors())
      IFPACK_CHK_ERR(-1); // not valid

    // need an additional vector for AztecOO preconditioning
    // (as X and Y both point to the same memory space)
    // FIXME: is overlap between blocks is zero, I can save
    // a vector
    Epetra_MultiVector Xtmp(X);

    if (ZeroStartingSolution_)
      Y.PutScalar(0.0);

    // ------------ //
    // single sweep //
    // ------------ //

    if (NumSweeps() == 1 && ZeroStartingSolution_
	&& (PrintFrequency() != 0)) {
      IFPACK_CHK_ERR(ApplyBSGS(Xtmp,Y));
      return(0);
    }

    // --------------------- //
    // general case (solver) //
    // --------------------- //

    if (PrintFrequency())
      Ifpack_PrintResidual(Label(),Matrix(),Y,Xtmp);

    // starting solution
    Epetra_MultiVector AX(Xtmp);
    Epetra_MultiVector Ynew(Y);

    for (int j = 0; j < NumSweeps() ; j++) {

      // compute the residual, I can skip first iteration
      // if starting solution is zero
      if (j || !ZeroStartingSolution_) {
	IFPACK_CHK_ERR(Apply(Y,AX));
	AX.Update(1.0,Xtmp,-1.0);
      }

      // apply the block diagonal of A and update the residual
      ApplyBSGS(AX,Ynew);

      Y.Update(DampingFactor(),Ynew,1.0);

      if (PrintFrequency() && (j != 0) && (j % PrintFrequency() == 0))
	Ifpack_PrintResidual(j,Matrix(),Y,Xtmp);

    }

    if (PrintFrequency())
      Ifpack_PrintResidual(NumSweeps(),Matrix(),Y,Xtmp);

    return(0);

  }

  //! Applies one sweep of block Jacobi.
  virtual int ApplyBSGS(Epetra_MultiVector& X, 
			Epetra_MultiVector& Y) const
  {

    // cycle over all local subdomains

    int Length = Matrix().MaxNumEntries();
    vector<int> Indices;
    vector<double> Values;
    Indices.resize(Length);
    Values.resize(Length);

    int NumMyRows = Matrix().NumMyRows();

    // ================== //
    // apply (D - E)^{-1} //
    // ================== //

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
	  Y[k][LID] = Containers_[i]->LHS(j,k);
	}
      }
    }

    // ============ //
    // apply D to Y //
    // ============ //

    for (int i = 0 ; i < NumLocalBlocks() ; ++i) {

      int LID, GID;

      for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
	LID = Containers_[i]->ID(j);
	for (int k = 0 ; k < Y.NumVectors() ; ++k) {
	  Containers_[i]->RHS(j,k) = Y[k][LID];
	}
      }

      IFPACK_CHK_ERR(Containers_[i]->Apply());

      for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
	LID = Containers_[i]->ID(j);
	for (int k = 0 ; k < Y.NumVectors() ; ++k) {
	  Y[k][LID] = Containers_[i]->LHS(j,k);
	}
      }
    }
 
    X = Y;

    // ================== //
    // apply (D - F)^{-1} //
    // ================== //

    for (int i = NumLocalBlocks() - 1; i >=0 ; --i) {

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

	  if ((*Partitioner_)(col) > i) {
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
	  Y[k][LID] = Containers_[i]->LHS(j,k);
	}
      }
    }
    return(0);
  }


private:

//! Sets the label of \c this object.
virtual int SetLabel()
{
  Label_ = "Ifpack_BlockSymGaussSeidel, # blocks = "
    + Ifpack_toString(NumLocalBlocks());
  return(0);
}

};

#endif // HAVE_IFPACK_TEUCHOS
#endif // class IFPACK_BLOCKSYMgAUSSSEIDEL
