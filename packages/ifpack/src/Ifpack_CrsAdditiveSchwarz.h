#ifndef IFPACK_CRSADDITIVESCHWARZ_H
#define IFPACK_CRSADDITIVESCHWARZ_H

#include "Ifpack_Preconditioner.h"
#include "Ifpack_Utils.h"
#include "Ifpack_Partitioner.h"
#include "Ifpack_METISPartitioner.h"
#include "Ifpack_AdditiveSchwarz.h"
#include "Ifpack_Graph_Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"

//! Ifpack_CrsAdditiveSchwarz: a class to define overlapping Schwarz preconditioner for Epetra_CrsMatrix's.

/*!
  Class Ifpack_CrsAdditiveSchwarz enables the creation of overlapping
  Schwarz preconditioner, with arbitrary overlap, for Epetra_CrsMatrix's.

  \date Sep-04
*/
template<class T>
class Ifpack_CrsAdditiveSchwarz : public Ifpack_AdditiveSchwarz<T> {

public:

  //! Creates an Ifpack_CrsAdditiveSchwarz object.
  /*! Creates an instance of Ifpack_CrsAdditiveSchwarz class.
   * \param In
   * Matrix - matrix to be preconditioned
   *
   * \param In
   * OverlapLevel - level of overlap (any positive number or zero).
   */
  Ifpack_CrsAdditiveSchwarz(Epetra_CrsMatrix* Matrix, int OverlapLevel);

  //! Destructor.
  ~Ifpack_CrsAdditiveSchwarz();

  //! Computes the preconditioners.
  /*! Computes the preconditioners: 
   *
   * \param In
   * List - list specifying the parameters for Jacobi. See above.
   *    
   * \return
   * 0 if successful, 1 if a zero element has been found on the
   * diagonal. In this latter case, the preconditioner can still be
   * applied, as the inverse of zero elements is replaced by 1.0.
   */
  virtual int Compute();

private:
    
  //! Sets up the preconditioner
  int SetUp();

  //! Pointer to the matrix to be preconditioned.
  Epetra_CrsMatrix* CrsMatrix_;
  //! Pointer to the graph of \c Matrix_.
  Ifpack_Graph* Graph_;
  //! Pointer to the graph of \c OverlappingMatrix_.
  Ifpack_Graph* OverlappingGraph_;
  //! Pointer to the partitioner.
  Ifpack_Partitioner* Partitioner_;
  //! Level of overlap among the processors.
  int OverlapLevel_;

};

//==============================================================================
template<class T>
Ifpack_CrsAdditiveSchwarz<T>::
Ifpack_CrsAdditiveSchwarz(Epetra_CrsMatrix* CrsMatrix,
			  const int OverlapLevel) :
  Ifpack_AdditiveSchwarz<T>(CrsMatrix),
  CrsMatrix_(CrsMatrix),
  Graph_(0),
  OverlappingGraph_(0),
  Partitioner_(0),
  OverlapLevel_(OverlapLevel)
{
  if (CrsMatrix_->Comm().NumProc() == 1)
    OverlapLevel_ = 0;

  if (OverlapLevel_)
    IsOverlapping_ = true;
}

//==============================================================================
template<class T>
Ifpack_CrsAdditiveSchwarz<T>::~Ifpack_CrsAdditiveSchwarz()
{
  if (OverlappingMatrix_)
    delete OverlappingMatrix_;
  
  if (OverlappingGraph_) 
    delete OverlappingGraph_;

  if (Graph_)
    delete Graph_;

  if (Partitioner_)
    delete Partitioner_;

}

//==============================================================================
template<typename T>
int Ifpack_CrsAdditiveSchwarz<T>::Compute()
{

  IFPACK_CHK_ERR(SetUp());
  IFPACK_CHK_ERR(Ifpack_AdditiveSchwarz<T>::SetUp());

  if (Inverse_ == 0)
    IFPACK_CHK_ERR(-1);

  if (LocalizedMatrix_ == 0)
    IFPACK_CHK_ERR(-1);

  // FIXME: add overlap
  Label_ = "Ifpack Additive Schwarz";

  IFPACK_CHK_ERR(Inverse_->SetParameters(List_));

  IFPACK_CHK_ERR(Inverse_->Compute());

  IsComputed_ = true;

  return(0);
}
//==============================================================================
template<class T>
int Ifpack_CrsAdditiveSchwarz<T>::SetUp()
{

  if (Matrix_ == 0) 
    IFPACK_CHK_ERR(-1);

  // FIXME: Graph is not always usefull (only if overlap = 0:
  Graph_ = new Ifpack_Graph_Epetra_CrsGraph(&(CrsMatrix_->Graph()));
  
  if (OverlapLevel_ > 0) {

    Epetra_CrsMatrix* CrsOverlappingMatrix;
    CrsOverlappingMatrix = Ifpack_CreateOverlappingCrsMatrix(CrsMatrix_,
							     OverlapLevel_);
    if (CrsOverlappingMatrix == 0)
      IFPACK_CHK_ERR(-1);

    OverlappingGraph_ = 
      new Ifpack_Graph_Epetra_CrsGraph(&(CrsOverlappingMatrix->Graph()));
    if (OverlappingGraph_ == 0)
      IFPACK_CHK_ERR(-1);

    // FIXME: more general??
    Partitioner_ = new Ifpack_METISPartitioner(OverlappingGraph_);
    if (Partitioner_ == 0)
      IFPACK_CHK_ERR(-1);

    OverlappingMatrix_ = CrsOverlappingMatrix;
  }
  else
    Partitioner_ = new Ifpack_METISPartitioner(Graph_);
  
  // might need to set if OverlapLevel_ has been changed
  List_.set("overlap level", OverlapLevel_);

  IFPACK_CHK_ERR(Partitioner_->SetParameters(List_));
  IFPACK_CHK_ERR(Partitioner_->Compute());
  List_.set("partitioner object", Partitioner_);

  return(0);
}
  
#endif
