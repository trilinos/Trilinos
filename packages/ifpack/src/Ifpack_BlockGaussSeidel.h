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

  Ifpack_BlockGaussSeidel(Epetra_RowMatrix* Matrix) :
    Ifpack_BlockPreconditioner<T>(Matrix)
  {
  };

  virtual ~Ifpack_BlockGaussSeidel()
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

    if (IsComputed() == false)
      IFPACK_CHK_ERR(-4); // need to compute the prec first

    if (X.NumVectors() != Y.NumVectors())
      IFPACK_CHK_ERR(-1); // not valid

    if (NumSweeps() == 1) {

      // simple case, I apply the preconditioner to the input     //
      // vector. No matrix-vector product, nor other allocations  //
      Y = X;
      IFPACK_CHK_ERR(ApplyBGS(Y));

    }
    else {

      // generic case with more than one sweep. This requires room for two
      // additional vectors, and a matrix-vector product for each sweeps. The
      // starting solution is the vector to be preconditioned.
      // copy rhs 
      Epetra_MultiVector Xtmp(X);

      // starting solution
      Epetra_MultiVector AX(Y);

      for (int j = 0; j < NumSweeps() ; j++) {

	// compute the residual
	IFPACK_CHK_ERR(Apply(Y,AX));

	AX.Update(1.0,Xtmp,-1.0);

	// apply the lower block triangular part of A
	ApplyBGS(AX);

	// update the residual
	Y.Update(DampingFactor(), AX, 1.0);

      }
    }
    return(0);
  }

private:

  //!Applies one sweep of Gauss-Seidel to Y, overwrites results on Y.
  int ApplyBGS(Epetra_MultiVector& Y) const
  {
    // cycle over all local subdomains

    int Length = Matrix()->MaxNumEntries();
    vector<int> Indices;
    vector<double> Values;
    Indices.resize(Length);
    Values.resize(Length);

    for (int i = 0 ; i < NumLocalBlocks() ; ++i) {

      int LID, GID;

      // update from previous block

      for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
	LID = Containers_[i]->ID(j);

	int NumEntries;
	IFPACK_CHK_ERR(Matrix()->ExtractMyRowCopy(LID, Length,NumEntries,
						    &Values[0], &Indices[0]));

	for (int k = 0 ; k < NumEntries ; ++k) {
	  int col = Indices[k];
	  if (col >= Matrix()->NumMyRows() ) 
	    continue;

	  if ((*Partitioner_)(col) < i) {
	    for (int kk = 0 ; kk < Y.NumVectors() ; ++kk) {
	      Y[kk][LID] = Y[kk][LID] - Values[k] * Y[kk][col];
	    }
	  }
	}
      }

      // solve with this block

      for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
	LID = Containers_[i]->ID(j);
	for (int k = 0 ; k < Y.NumVectors() ; ++k) {
	  Containers_[i]->RHS(j,k) = Y[k][LID];
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

  //! Sets the label of \c this object.
  virtual int SetLabel()
  {
    Label_ = "Amesos_BlockGaussSeidel, # blocks = "
      + Ifpack_toString(NumLocalBlocks());
   }

};

#endif // HAVE_IFPACK_TEUCHOS
#endif // class Ifpack_BlockGaussSeidel
