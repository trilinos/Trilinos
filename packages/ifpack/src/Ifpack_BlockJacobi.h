#ifndef IFPACK_BLOCKJACOBI_H
#define IFPACK_BLOCKJACOBI_H

#include "Ifpack_ConfigDefs.h"
#ifdef HAVE_IFPACK_TEUCHOS
#include "Ifpack_BlockPreconditioner.h" 
#include "Ifpack_Utils.h"

//! Ifpack_BlockJacobi: a class to define block-Jacobi preconditioners preconditioners of Epetra_RowMatrix's.

/*!
  The Ifpack_BlockJacobi class enables the construction of block-Jacobi 
  preconditioners.
*/

template<class T>
class Ifpack_BlockJacobi : public Ifpack_BlockPreconditioner<T> {

public:

  //! Constructor.
  Ifpack_BlockJacobi(Epetra_RowMatrix* Matrix) :
    Ifpack_BlockPreconditioner<T>(Matrix)
  {
  };
  
  //! Destructor.
  virtual ~Ifpack_BlockJacobi()
  {};

  //! Applies the block Jacobi preconditioner to X, returns the result in Y.
  /*! 
    \param 
    X - (In) A Epetra_MultiVector of dimension NumVectors to be preconditioned.
    \param
    Y -(Out) A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.

    \warning I do suppose that X and Y are two different objects.
    This means that this class cannot be directly used with AztecOO,
    and it should be used through Ifpack_AdditiveSchwarz only.
    */
  virtual int ApplyInverse(const Epetra_MultiVector& X, 
			   Epetra_MultiVector& Y) const
  {

    if (IsComputed() == false)
      IFPACK_CHK_ERR(-4); // need to compute the prec first

    if (X.NumVectors() != Y.NumVectors())
      IFPACK_CHK_ERR(-1); // not valid

    if (ZeroStartingSolution_)
      Y.PutScalar(0.0);

    // ------------ //
    // single sweep //
    // ------------ //

    if (NumSweeps() == 1 && ZeroStartingSolution_
	&& (PrintFrequency() == 0)) {
      IFPACK_CHK_ERR(ApplyBJ(X,Y));
      return(0);
    }

    // --------------------- //
    // general case (solver) //
    // --------------------- //

    if (PrintFrequency())
      Ifpack_PrintResidual(Label(),Matrix(),Y,X);

    // starting solution
    Epetra_MultiVector AX(Y);

    for (int j = 0; j < NumSweeps() ; j++) {

      // compute the residual
      IFPACK_CHK_ERR(Apply(Y,AX));

      AX.Update(1.0,X,-1.0);

      // apply the block diagonal of A and update
      // the residual
      ApplyBJ(AX,Y);

      if (PrintFrequency() && (j != 0) && (j % PrintFrequency() == 0))
	Ifpack_PrintResidual(j,Matrix(),Y,X);

    }

    if (PrintFrequency())
      Ifpack_PrintResidual(NumSweeps(),Matrix(),Y,X);

    return(0);
  }

  //! Applies one sweep of block Jacobi.
  virtual int ApplyBJ(const Epetra_MultiVector& X, 
		      Epetra_MultiVector& Y) const
  {
    // cycle over all local subdomains
    for (int i = 0 ; i < NumLocalBlocks() ; ++i) {

      int LID, GID;

      // extract RHS from X
      for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
	LID = Containers_[i]->ID(j);
	for (int k = 0 ; k < X.NumVectors() ; ++k) {
	  Containers_[i]->RHS(j,k) = X[k][LID];
	}
      }

      // apply the inverse of each block
      IFPACK_CHK_ERR(Containers_[i]->ApplyInverse());

      // copy back into solution vector Y
      for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
	LID = Containers_[i]->ID(j);
	for (int k = 0 ; k < Y.NumVectors() ; ++k) {
	  Y[k][LID] = Y[k][LID] + DampingFactor() * (*W_)[LID] * Containers_[i]->LHS(j,k);
	}
      }
    }
    return(0);
  }

  //! Sets label.
  virtual int SetLabel()
  {
    Label_ = "Ifpack_BlockJacobi, # blocks = "
      + Ifpack_toString(NumLocalBlocks());
    return(0);
  }

};

#endif // HAVE_IFPACK_TEUCHOS
#endif // class Ifpack_BlockJacobi
