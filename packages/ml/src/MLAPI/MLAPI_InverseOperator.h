#ifndef ML_INVERSEOPERATOR_H
#define ML_INVERSEOPERATOR_H
#include "Epetra_Vector.h"
#include "ml_include.h"
#include "ml_RowMatrix.h"
#include "Teuchos_RefCountPtr.hpp"
#include "Ifpack_Preconditioner.h"
#include "Ifpack_PointRelaxation.h"
#include "Ifpack_IC.h"
#include "Ifpack_ICT.h"
#include "Ifpack_ILU.h"
#include "Ifpack_ILUT.h"
#include "Ifpack_AdditiveSchwarz.h"
#include "Ifpack_Amesos.h"
#include "MLAPI_Space.h"
#include "MLAPI_DoubleVector.h"
#include "MLAPI_DoubleVector_Utils.h"
#include "MLAPI_Workspace.h"
#include "MLAPI_DataBase.h"

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

class InverseOperator : public BaseObject {

public:
  //@{ Constructors and destructors.

  //! Basic constructor.
  InverseOperator()
  {}
  
  //! Copy constructor.
  InverseOperator(const InverseOperator& RHS)
  {
    Op_ = RHS.GetOperator();
    RowMatrix_ = RHS.RowMatrix();
    Data_ = RHS.GetData();
    SetLabel(RHS.GetLabel());
  }

  //! Destructor.
  ~InverseOperator()
  {}

  // @}
  // @{ Overloaded operators

  //! Operator =.
  InverseOperator& operator=(const InverseOperator& RHS)
  {
    if (this == &RHS)
      return(*this);

    Op_ = RHS.GetOperator();
    RowMatrix_ = RHS.RowMatrix();
    Data_ = RHS.GetData();

    SetLabel(RHS.GetLabel());
    return(*this);
  }

  // @}
  // @{ Reshaping methods

  //! Reshapes the object by setting the Operator and the specified type.
  /*! Reshapes the object by preparing \c this object to apply the inverse
   * of the input operator \c Op, as specified by \c Type.
   * 
   * \param Op (In) - Operator of which \c this object defines the inverse
   *
   * \param Type (In) - String containing the method that will be used
   *              to define the action of the inverse.
   *
   *                                      a double value.
   */
  void Reshape(const Operator& Op, const SmootherDataBase& DB)
  {
    Op_ = Op;

    RowMatrix_ = Teuchos::rcp(new ML_Epetra::RowMatrix(Op.GetData(),
                                                       &GetEpetra_Comm()));

    string Type = DB.GetType();

    // FIXME: to add overlap and level-of-fill
    int NumSweeps   = DB.GetSweeps();
    double Damping  = DB.GetDamping();
    int LOF_ilu     = DB.GetILUFill();
    double LOF_ict  = DB.GetICTFill();
    double LOF_ilut = DB.GetILUTFill();

    Teuchos::ParameterList IFPACKList;
    IFPACKList.set("relaxation: sweeps", NumSweeps);
    IFPACKList.set("relaxation: damping factor", Damping);
    IFPACKList.set("fact: level-of-fill", LOF_ilu);
    IFPACKList.set("fact: ict level-of-fill", LOF_ict);
    IFPACKList.set("fact: ilut level-of-fill", LOF_ilut);
    IFPACKList.set("relaxation: zero starting solution", false);
    
    bool verbose = (MyPID() == 0 && GetPrintLevel() > 5);

    Ifpack_Preconditioner* Prec;

    if (verbose) {
      cout << "Build smoother `" << Type << "'" << endl;
//      cout << "# global rows = " << RowMatrix_.get()->NumGlobalRows() << endl;
    }

    if (Type == "Jacobi") {
      if (verbose) {
        cout << "Damping factor = " << Damping 
             << ", sweeps = " << NumSweeps << endl;
        cout << endl;
      }
      IFPACKList.set("relaxation: type", "Jacobi");
      Prec = new Ifpack_PointRelaxation(RowMatrix_.get());
    }
    else if (Type == "Gauss-Seidel") {
      if (verbose) {
        cout << "Damping factor = " << Damping 
             << ", sweeps = " << NumSweeps << endl;
        cout << endl;
      }
      IFPACKList.set("relaxation: type", "Gauss-Seidel");
      Prec = new Ifpack_PointRelaxation(RowMatrix_.get());
    }
    else if (Type == "symmetric Gauss-Seidel") {
      if (verbose) {
        cout << "Damping factor = " << Damping 
             << ", sweeps = " << NumSweeps << endl;
        cout << endl;
      }
      IFPACKList.set("relaxation: type", "symmetric Gauss-Seidel");
      Prec = new Ifpack_PointRelaxation(RowMatrix_.get());
    }
    else if (Type == "ILU") {
      if (verbose) {
        cout << "ILU factorization, ov = 0, no reordering, LOF = "
             << LOF_ilu << endl;
        cout << endl;
      }
      Prec = new Ifpack_ILU(RowMatrix_.get());
    }
    else if (Type == "ILUT") {
      if (verbose) {
        cout << "ILUT factorization, ov = 0, no reordering, LOF = "
             << LOF_ilu << endl;
        cout << endl;
      }
      Prec = new Ifpack_ILUT(RowMatrix_.get());
    }
    else if (Type == "IC") {
      if (verbose) {
        cout << "IC factorization, ov = 0, no reordering, LOF = "
             << LOF_ilu << endl;
        cout << endl;
      }
      Prec = new Ifpack_IC(RowMatrix_.get());
    }
    else if (Type == "ICT") {
      if (verbose) {
        cout << "ICT factorization, ov = 0, no reordering, LOF = "
             << LOF_ilu << endl;
        cout << endl;
      }
      Prec = new Ifpack_ICT(RowMatrix_.get());
    }
    else if (Type == "LU") {
      if (verbose) {
        cout << "LU factorization, ov = 0, local solver = KLU" << endl;
        cout << endl;
      }
      Prec = new Ifpack_AdditiveSchwarz<Ifpack_Amesos>(RowMatrix_.get());
    }
    else if (Type == "Amesos" || Type == "Amesos-KLU")  {
      if (verbose) {
        cout << "Amesos-KLU direct solver" << endl;
        cout << endl;
      }
      Prec = new Ifpack_Amesos(RowMatrix_.get());
    }
    else
      ML_THROW("Requested type (" + Type + ") not recognized", -1);

    Data_ = Teuchos::rcp(Prec);

    Data_->SetParameters(IFPACKList);
    Data_->Initialize();
    Data_->Compute();
  }

  void Reshape(const Operator& Op, const CoarseSolverDataBase& DB)
  {
    Op_ = Op;

    RowMatrix_ = Teuchos::rcp(new ML_Epetra::RowMatrix(Op.GetData(),
                                                       &GetEpetra_Comm()));

    string Type = DB.GetType();

    bool verbose = (MyPID() == 0 && GetPrintLevel() > 5);

    Ifpack_Preconditioner* Prec;

    if (verbose) {
      cout << "Build coarse solver `" << Type << "'" << endl;
//      cout << "# global rows = " << RowMatrix_.get()->NumGlobalRows() << endl;
    }

    if (Type == "Amesos" || Type == "Amesos-KLU")  {
      if (verbose) {
        cout << "Amesos-KLU direct solver" << endl;
        cout << endl;
      }
      Prec = new Ifpack_Amesos(RowMatrix_.get());
    }
    else
      ML_THROW("Requested type (" + Type + ") not recognized", -1);

    Data_ = Teuchos::rcp(Prec);

    Data_->Initialize();
    Data_->Compute();
  }

  // @}
  // @{ Query methods
  
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
  const Teuchos::RefCountPtr<Ifpack_Preconditioner> GetData() const
  {
    return(Data_);
  }

  //! Returns a reference to the Operator of which \c this object defines the inverse.
  const Operator& GetOperator() const
  {
    return(Op_);
  }

  //! Prints out basic information about \c this object.
  ostream& Print(std::ostream& os, const bool verbose = true) const
  {

    // FIXME: to be completed in some way???
    if (MyPID() == 0) {
      os << "InverseOperator `" << GetLabel() << "'" << endl;
    }

    return(os);

  }

  // @}
  // @{ Mathematical methods
  
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

  //! Applies the operator to LHS, returns the results.
  DoubleVector operator()(const DoubleVector& LHS)
  {
    DoubleVector RHS(LHS.VectorSpace());
    RHS = 0.0;
    ApplyInverse(LHS,RHS);
    return(RHS);
  }

  //! Applies the operator to LHS using RHS as initial solution, returns the results.
  DoubleVector operator()(const DoubleVector& LHS,
                          const DoubleVector& RHS)
  {
    DoubleVector RHS2 = Duplicate(RHS);
    ApplyInverse(LHS,RHS2);
    return(RHS2);
  }

private:

  // @}
  // @{ Internal data

  //! Operator of which \c this object define the inverse.
  Operator Op_;
  //! Wrapper for IFPACK
  Teuchos::RefCountPtr<ML_Epetra::RowMatrix> RowMatrix_;
  //! IFPACK preconditioner.
  Teuchos::RefCountPtr<Ifpack_Preconditioner> Data_;

  // @}
  
}; // InverseOperator

} // namespace MLAPI
#endif // ML_INVERSEOPERATOR_H
