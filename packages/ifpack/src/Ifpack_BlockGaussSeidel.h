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

    // need this for AztecOO
    Epetra_MultiVector Xtmp(X);
    Epetra_MultiVector Xtmp2(X);

    if (ZeroStartingSolution_)
      Y.PutScalar(0.0);

    for (int j = 0; j < NumSweeps() ; j++) {

      // cycle over all local subdomains

      int Length = Matrix().MaxNumEntries();
      vector<int> Indices;
      vector<double> Values;
      Indices.resize(Length);
      Values.resize(Length);

      int NumMyRows = Matrix().NumMyRows();

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
		Xtmp2[kk][LID] = Xtmp[kk][LID] - Values[k] * Y[kk][col];
	      }
	    }
	  }
	}

	// solve with this block

	for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
	  LID = Containers_[i]->ID(j);
	  for (int k = 0 ; k < Y.NumVectors() ; ++k) {
	    Containers_[i]->RHS(j,k) = Xtmp2[k][LID];
	  }
	}

	IFPACK_CHK_ERR(Containers_[i]->ApplyInverse());

	for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
	  LID = Containers_[i]->ID(j);
	  for (int k = 0 ; k < Y.NumVectors() ; ++k) {
	    Y[k][LID] = Y[k][LID] + DampingFactor() * Containers_[i]->LHS(j,k);
	  }
	}
      }

    }
    return(0);
  }

private:

//! Sets the label of \c this object.
virtual int SetLabel()
{
  Label_ = "Amesos_BlockGaussSeidel, # blocks = "
    + Ifpack_toString(NumLocalBlocks());
}

};

#endif // HAVE_IFPACK_TEUCHOS
#endif // class Ifpack_BlockGaussSeidel
