#ifndef ML_INVERSEOPERATOR_H
#define ML_INVERSEOPERATOR_H
#include "Epetra_Vector.h"
#include "ml_include.h"
#include "ml_RowMatrix.h"
#include "Teuchos_RefCountPtr.hpp"
#include "Ifpack_Preconditioner.h"
#include "Ifpack_PointRelaxation.h"
#include "Ifpack_ICT.h"
#include "Ifpack_Amesos.h"
#include "MLAPI_Space.h"
#include "MLAPI_DoubleVector.h"

using namespace std;

namespace MLAPI {

/*!
 * \class InverseOperator
 *
 * \brief InverseOperator: basic class to define smoother and coarse solvers.
 *
 *
 * \author Marzio Sala, SNL 9214
 *
 * \date Last updated on 07-Jan-05
 */

class InverseOperator {

public:

  //! Basic constructor.
  InverseOperator()
  {}
  
  //! Copy constructor.
  InverseOperator(const InverseOperator& RHS)
  {
    Op_ = RHS.GetOperator();
    RowMatrix_ = RHS.RowMatrix();
    Data_ = RHS.Data();
  }

  //! Destructor.
  ~InverseOperator()
  {}

  //! Operator =.
  InverseOperator& operator=(const InverseOperator& RHS)
  {
    if (this == &RHS)
      return(*this);

    Op_ = RHS.GetOperator();
    RowMatrix_ = RHS.RowMatrix();
    Data_ = RHS.Data();
    return(*this);
  }

  //! Reshapes the objecct by setting the Operator and the specified type.
  /*! Reshapes the object by preparing \c this object to apply the inverse
   * of the input operator \c Op, as specified by \c Type.
   * 
   * \param Op (In) - Operator of which \c this object defines the inverse
   *
   * \param Type (In) - String containing the method that will be used
   *              to define the action of the inverse.
   *
   * \param List (In) - List containing all parameters that may define
   *              the action of the inverse.
   *
   * At the moment, \c Type can be one of the following:
   * - \c "SGS": symmetric Gauss-Seidel (through IFPACK)
   * - \c "Amesos": direct solver (through Amesos)
   */
  void Reshape(const Operator& Op, const string Type,
               Teuchos::ParameterList& List)
  {
    Op_ = Op;

    RowMatrix_ = Teuchos::rcp(new ML_Epetra::RowMatrix(Op.GetOperator(),
                                                       &GetEpetraComm()));
    Ifpack_Preconditioner* Prec;
    if (Type == "SGS") 
      Prec = new Ifpack_PointRelaxation(RowMatrix_.get());
    else if (Type == "ILU")
      Prec = new Ifpack_ICT(RowMatrix_.get());
    else if (Type == "Amesos") 
      Prec = new Ifpack_Amesos(RowMatrix_.get());
    else
      throw("preconditioner type not recognized");

    Data_ = Teuchos::rcp(Prec);

    List.set("relaxation: zero starting solution", false);
    Data_->SetParameters(List);
    Data_->Initialize();
    Data_->Compute();
  }

  //! Returns a reference to the Operator of which \c this object defines the inverse.
  const Operator& GetOperator() const
  {
    return(Op_);
  }

  //! Applies \c this object to vector \c lhs, returns values in \c rhs.
  int ApplyInverse(const DoubleVector& lhs, DoubleVector& rhs) const
  {
    Epetra_Vector elhs(View,RowMatrix_->OperatorDomainMap(),
                       (double*)&(lhs(0)));
    Epetra_Vector erhs(View,RowMatrix_->OperatorRangeMap(),
                       (double*)&(rhs(0)));

    Data_->ApplyInverse(elhs,erhs);
    return(0);
  }

  //! Returns a reference to the range space of \c this object.
  const Space& RangeSpace() const {
    return(Op_.RangeSpace());
  }

  //! Returns a reference to the domain space of \c this object.
  const Space& DomainSpace() const {
    return(Op_.DomainSpace());
  }

  //! Returns pointer of the internally stored ML_Epetra::RowMatrix object.
  const Teuchos::RefCountPtr<ML_Epetra::RowMatrix> RowMatrix() const
  {
    return(RowMatrix_);
  }

  //! Returns a pointer to the internally stored IFPACK preconditioner.
  const Teuchos::RefCountPtr<Ifpack_Preconditioner> Data() const
  {
    return(Data_);
  }

private:
  //! Operator of which \c this object define the inverse.
  Operator Op_;
  //! Wrapper for IFPACK
  Teuchos::RefCountPtr<ML_Epetra::RowMatrix> RowMatrix_;
  //! IFPACK preconditioner.
  Teuchos::RefCountPtr<Ifpack_Preconditioner> Data_;

}; // InverseOperator

} // namespace MLAPI
#endif // ML_INVERSEOPERATOR_H
