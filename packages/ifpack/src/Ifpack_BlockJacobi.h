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

  Ifpack_BlockJacobi(Epetra_RowMatrix* Matrix) :
    Ifpack_BlockPreconditioner<T>(Matrix)
  {
  };

  virtual ~Ifpack_BlockJacobi()
  {};

  //! Applies the block Jacobi preconditioner to X, returns the result in Y.
  /*! 
    \param In
    X - A Epetra_MultiVector of dimension NumVectors to be preconditioned.
    \param Out
    Y -A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.

    \warning In order to work with AztecOO, any implementation of this method 
    must support the case where X and Y are the same object.
    */
  virtual int ApplyInverse(const Epetra_MultiVector& X, 
			   Epetra_MultiVector& Y) const
  {

    // cycle over all local subdomains
    for (int i = 0 ; i < NumLocalBlocks() ; ++i) {

      int LID, GID;

      for (int j = 0 ; j < Partitioner_->NumRowsInPart(i) ; ++j) {
	LID = Containers_[i]->ID(j);
	for (int k = 0 ; k < X.NumVectors() ; ++k) {
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

  virtual int SetLabel()
  {
    Label_ = "Amesos_BlockJacobi, # blocks = "
      + Ifpack_toString(NumLocalBlocks());
   }

};

#endif // HAVE_IFPACK_TEUCHOS
#endif // class Ifpack_BlockJacobi
