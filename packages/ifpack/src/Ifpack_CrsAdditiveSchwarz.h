#ifndef IFPACK_CRSADDITIVESCHWARZ_H
#define IFPACK_CRSADDITIVESCHWARZ_H

#include "Ifpack_Preconditioner.h"
#include "Ifpack_Utils.h"
#include "Ifpack_Partitioner.h"
#include "Ifpack_LinearPartitioner.h"
#include "Ifpack_GreedyPartitioner.h"
#include "Ifpack_METISPartitioner.h"
#include "Ifpack_EquationPartitioner.h"
#include "Ifpack_AdditiveSchwarz.h"
#include "Ifpack_Graph_Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"

//! Ifpack_CrsAdditiveSchwarz: a class to define overlapping Schwarz preconditioner for Epetra_CrsMatrix's.

/*!
Class Ifpack_CrsAdditiveSchwarz enables the creation of overlapping
Schwarz preconditioner, with arbitrary overlap, for Epetra_CrsMatrix's.

The only difference between this class and Ifpack_AdditiveSchwarz is that
the user can pass just one matrix (the one to be preconditioned), plus
and integer defining the amount of overlap. By convention, the
minimal-overlap case is indicated by 0-overlap.
As the input matrix is an Epetra_CrsMatrix, \c Import() method can
be used to create the overlapping matrix. This matrix is automatically deleted
by the destructor. 

\note No overlapping matrix is constructed if the overlap level is 0.

This class is templated with the local solver type. Please refer to
the documentation of Ifpack_AdditiveSchwarz for examples of use.

\date Sep-04
*/
template<class T>
class Ifpack_CrsAdditiveSchwarz : public Ifpack_AdditiveSchwarz<T> {

public:

  //! Creates an Ifpack_CrsAdditiveSchwarz object.
  /*! Creates an instance of Ifpack_CrsAdditiveSchwarz class.
   * \param 
   * Matrix - (In) matrix to be preconditioned
   *
   * \param 
   * OverlapLevel - (In) level of overlap (any positive number or zero).
   */
  Ifpack_CrsAdditiveSchwarz(Epetra_CrsMatrix* Matrix, int OverlapLevel);

  //! Destructor.
  ~Ifpack_CrsAdditiveSchwarz();

  //! Computes the preconditioners.
  /*! Computes the preconditioners: 
   *
   * \return
   * 0 if successful, 1 otherwise.
   */
  virtual int Compute();

private:
    
  //! Sets up the preconditioner
  int Setup();

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
  Label_ = "Ifpack_CrsAdditiveSchwarz";

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

  IFPACK_CHK_ERR(Setup());
  IFPACK_CHK_ERR(Ifpack_AdditiveSchwarz<T>::Setup());

  if (Inverse_ == 0)
    IFPACK_CHK_ERR(-1);

  if (LocalizedMatrix_ == 0)
    IFPACK_CHK_ERR(-1);

  IFPACK_CHK_ERR(Inverse_->SetParameters(List_));
  IFPACK_CHK_ERR(Inverse_->Compute());

  Label_ = "Ifpack_AdditiveSchwarz, ov = "
    + Ifpack_toString(OverlapLevel_)
    + " (" + Inverse_->Label() + ")";

  IsComputed_ = true;

  return(0);
}
//==============================================================================
template<class T>
int Ifpack_CrsAdditiveSchwarz<T>::Setup()
{

  if (Matrix_ == 0) 
    IFPACK_CHK_ERR(-1);

  string Type = List_.get("partitioner: type", "metis");
  
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


    if (Type == "metis") 
      Partitioner_ = new Ifpack_METISPartitioner(OverlappingGraph_);
    else if (Type == "linear")
      Partitioner_ = new Ifpack_LinearPartitioner(OverlappingGraph_);
    else if (Type == "greedy")
      Partitioner_ = new Ifpack_GreedyPartitioner(OverlappingGraph_);
    else if (Type == "equation")
      Partitioner_ = new Ifpack_EquationPartitioner(OverlappingGraph_);
    else
      IFPACK_CHK_ERR(-1);

    if (Partitioner_ == 0)
      IFPACK_CHK_ERR(-1);

    OverlappingMatrix_ = CrsOverlappingMatrix;
  }
  else {

    Graph_ = new Ifpack_Graph_Epetra_CrsGraph(&(CrsMatrix_->Graph()));

    if (Type == "metis") 
      Partitioner_ = new Ifpack_METISPartitioner(Graph_);
    else if (Type == "linear")
      Partitioner_ = new Ifpack_LinearPartitioner(Graph_);
    else if (Type == "greedy")
      Partitioner_ = new Ifpack_GreedyPartitioner(Graph_);
    else if (Type == "equation")
      Partitioner_ = new Ifpack_EquationPartitioner(Graph_);
    else
      IFPACK_CHK_ERR(-1);

    if (Partitioner_ == 0)
      IFPACK_CHK_ERR(-1);
  }
  
  IFPACK_CHK_ERR(Partitioner_->SetParameters(List_));
  IFPACK_CHK_ERR(Partitioner_->Compute());
  List_.set("partitioner: object", Partitioner_);

  return(0);
}
  
#endif
