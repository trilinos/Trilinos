#ifndef ML_INVERSEOPERATOR_H
#define ML_INVERSEOPERATOR_H

#include "ml_common.h"
#include "ml_include.h"
#include "ml_RowMatrix.h"
#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_Vector.h"
#include "Ifpack_Preconditioner.h"
#include "Ifpack_PointRelaxation.h"
#include "Ifpack_IC.h"
#include "Ifpack_ICT.h"
#include "Ifpack_ILU.h"
#include "Ifpack_ILUT.h"
#include "Ifpack_AdditiveSchwarz.h"
#include "Ifpack_Amesos.h"
#include "MLAPI_Space.h"
#include "MLAPI_BaseOperator.h"
#include "MLAPI_CompObject.h"
#include "MLAPI_TimeObject.h"
#include "MLAPI_MultiVector.h"
#include "MLAPI_Workspace.h"

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
 * \date Last updated on Feb-05.
 */

class InverseOperator : public BaseOperator, public CompObject, public TimeObject {

public:

  // @{ \name Constructors and destructors.

  //! Empty constructor.
  InverseOperator()
  {}
  
  //! Constructor for a given Operator and type, and default parameters.
  InverseOperator(const Operator& Op, const string Type)
  {
    Reshape(Op, Type);
  }

  //! Constructor for a given Operator, type and parameters.
  InverseOperator(const Operator& Op, const string Type,
                  Teuchos::ParameterList& List)
  {
    Reshape(Op, Type, List);
  }

  //! Copy constructor.
  InverseOperator(const InverseOperator& RHS)
  {
    Op_           = RHS.GetOperator();
    RCPRowMatrix_ = RHS.RCPRowMatrix();
    RCPData_      = RHS.GetRCPData();
    SetLabel(RHS.GetLabel());
  }

  //! Destructor.
  ~InverseOperator()
  {}

  // @}
  // @{ \name Overloaded operators

  //! Operator =.
  InverseOperator& operator=(const InverseOperator& RHS)
  {
    if (this == &RHS)
      return(*this);

    Op_           = RHS.GetOperator();
    RCPRowMatrix_ = RHS.RCPRowMatrix();
    RCPData_      = RHS.GetRCPData();

    SetLabel(RHS.GetLabel());
    return(*this);
  }

  // @}
  // @{ Reshaping methods

  //! Reshapes the object with default values.
  void Reshape(const Operator& Op, const string Type)
  {
    Teuchos::ParameterList List;
    Reshape(Op, Type, List);
  }

  //! Reshapes the object by setting the Operator and the specified type.
  void Reshape(const Operator& Op, const string Type,
               Teuchos::ParameterList& List)
  {
    ResetTimer();

    Op_ = Op;

    RCPRowMatrix_ = Teuchos::rcp(new ML_Epetra::RowMatrix(Op.GetML_Operator(),
                                                       &GetEpetra_Comm()));

    // FIXME: to add overlap and level-of-fill
    int NumSweeps   = List.get("smoother: sweeps", 1);
    double Damping  = List.get("smoother: damping factor", 0.67); 
    int LOF_ilu     = List.get("smoother: ilu fill", 0);
    double LOF_ict  = List.get("smoother: ilut fill", 1.0);
    double LOF_ilut = List.get("smoother: ict fill", 1.0);

    Teuchos::ParameterList IFPACKList;
    IFPACKList.set("relaxation: sweeps", NumSweeps);
    IFPACKList.set("relaxation: damping factor", Damping);
    IFPACKList.set("fact: level-of-fill", LOF_ilu);
    IFPACKList.set("fact: ict level-of-fill", LOF_ict);
    IFPACKList.set("fact: ilut level-of-fill", LOF_ilut);
    IFPACKList.set("relaxation: zero starting solution", false);
    
    bool verbose = false; //(GetMyPID() == 0 && GetPrintLevel() > 5);

    Ifpack_Preconditioner* Prec;

    if (Type == "Jacobi") {
      if (verbose) {
        cout << "Damping factor = " << Damping 
             << ", sweeps = " << NumSweeps << endl;
        cout << endl;
      }
      IFPACKList.set("relaxation: type", "Jacobi");
      Prec = new Ifpack_PointRelaxation(RowMatrix());
    }
    else if (Type == "Gauss-Seidel") {
      if (verbose) {
        cout << "Damping factor = " << Damping 
             << ", sweeps = " << NumSweeps << endl;
        cout << endl;
      }
      IFPACKList.set("relaxation: type", "Gauss-Seidel");
      Prec = new Ifpack_PointRelaxation(RowMatrix());
    }
    else if (Type == "symmetric Gauss-Seidel") {
      if (verbose) {
        cout << "Damping factor = " << Damping 
             << ", sweeps = " << NumSweeps << endl;
        cout << endl;
      }
      IFPACKList.set("relaxation: type", "symmetric Gauss-Seidel");
      Prec = new Ifpack_PointRelaxation(RowMatrix());
    }
    else if (Type == "ILU") {
      if (verbose) {
        cout << "ILU factorization, ov = 0, no reordering, LOF = "
             << LOF_ilu << endl;
        cout << endl;
      }
      Prec = new Ifpack_ILU(RowMatrix());
    }
    else if (Type == "ILUT") {
      if (verbose) {
        cout << "ILUT factorization, ov = 0, no reordering, LOF = "
             << LOF_ilu << endl;
        cout << endl;
      }
      Prec = new Ifpack_ILUT(RowMatrix());
    }
    else if (Type == "IC") {
      if (verbose) {
        cout << "IC factorization, ov = 0, no reordering, LOF = "
             << LOF_ilu << endl;
        cout << endl;
      }
      Prec = new Ifpack_IC(RowMatrix());
    }
    else if (Type == "ICT") {
      if (verbose) {
        cout << "ICT factorization, ov = 0, no reordering, LOF = "
             << LOF_ilu << endl;
        cout << endl;
      }
      Prec = new Ifpack_ICT(RowMatrix());
    }
    else if (Type == "LU") {
      if (verbose) {
        cout << "LU factorization, ov = 0, local solver = KLU" << endl;
        cout << endl;
      }
      Prec = new Ifpack_AdditiveSchwarz<Ifpack_Amesos>(RowMatrix());
    }
    else if (Type == "Amesos" || Type == "Amesos-KLU")  {
      if (verbose) {
        cout << "Amesos-KLU direct solver" << endl;
        cout << endl;
      }
      Prec = new Ifpack_Amesos(RowMatrix());
    }
    else
      ML_THROW("Requested type (" + Type + ") not recognized", -1);

    RCPData_ = Teuchos::rcp(Prec);

    RCPData_->SetParameters(IFPACKList);
    RCPData_->Initialize();
    RCPData_->Compute();

    UpdateFlops(RCPData_->InitializeFlops());
    UpdateFlops(RCPData_->ComputeFlops());
    UpdateTime();

  }

  // @}
  // @{ Get and Set methods.
  
  //! Returns a reference to the range space of \c this object.
  const Space GetOperatorRangeSpace() const {
    return(Op_.GetRangeSpace());
  }

  //! Returns a reference to the domain space of \c this object.
  const Space GetOperatorDomainSpace() const {
    return(Op_.GetDomainSpace());
  }

  //! Returns a reference to the range space of \c this object.
  const Space GetRangeSpace() const {
    return(Op_.GetRangeSpace());
  }

  //! Returns a reference to the domain space of \c this object.
  const Space GetDomainSpace() const {
    return(Op_.GetDomainSpace());
  }

  //! Returns pointer of the internally stored ML_Epetra::RowMatrix object.
  const Teuchos::RefCountPtr<ML_Epetra::RowMatrix> RCPRowMatrix() const
  {
    return(RCPRowMatrix_);
  }

  //! Returns pointer of the internally stored ML_Epetra::RowMatrix object.
  ML_Epetra::RowMatrix* RowMatrix() const
  {
    return(RCPRowMatrix_.get());
  }

  //! Returns a reference to the Operator of which \c this object defines the inverse.
  const Operator& GetOperator() const
  {
    return(Op_);
  }

  //! Returns a pointer to the internally stored IFPACK preconditioner.
  Teuchos::RefCountPtr<Ifpack_Preconditioner>& GetRCPData() 
  {
    return(RCPData_);
  }

  //! Returns a pointer to the internally stored IFPACK preconditioner.
  const Teuchos::RefCountPtr<Ifpack_Preconditioner>& GetRCPData() const
  {
    return(RCPData_);
  }

  // @}
  // @{ Mathematical methods
  
  //! Applies \c this object to vector \c lhs, returns values in \c rhs.
  int Apply(const MultiVector& x, MultiVector& y) const
  {
    ResetTimer();

    int x_nv = x.GetNumVectors();
    int y_nv = y.GetNumVectors();
    double FL = RCPData_->ComputeFlops();

    if (x_nv != y_nv)
      ML_THROW("Number of vectors of x and y differ (" +
               GetString(x_nv) + " vs. " + GetString(x_nv), -1);
                
    for (int v = 0 ; v < x_nv ; ++v) {

      Epetra_Vector x_Epetra(View,RowMatrix()->OperatorDomainMap(),
                             (double*)&(x(0,v)));
      Epetra_Vector y_Epetra(View,RowMatrix()->OperatorRangeMap(),
                             (double*)&(y(0,v)));

      RCPData_->ApplyInverse(x_Epetra,y_Epetra);
    }

    UpdateFlops(RCPData_->ComputeFlops() - FL);
    UpdateTime();

    return(0);
  }

  //! Applies the operator to LHS, returns the results.
  MultiVector operator()(const MultiVector& LHS)
  {
    MultiVector RHS(LHS.GetVectorSpace());
    RHS = 0.0;
    Apply(LHS,RHS);
    return(RHS);
  }

  //! Applies the operator to LHS using RHS as initial solution, returns the results.
  MultiVector operator()(const MultiVector& LHS,
                         const MultiVector& RHS)
  {
    MultiVector RHS2 = Duplicate(RHS);
    Apply(LHS,RHS2);
    return(RHS2);
  }

  // @}
  // @{ \name Miscellaneous methods

  //! Prints out basic information about \c this object.
  ostream& Print(std::ostream& os, const bool verbose = true) const
  {

    if (GetMyPID() == 0) {
      os << "***MLAPI::InverseOperator" << endl;
      os << "Label             = " << GetLabel() << endl;
      os << "Number of rows    = " << GetRangeSpace().GetNumGlobalElements() << endl;
      os << "Number of columns = " << GetRangeSpace().GetNumGlobalElements() << endl;
      os << "Flop count        = " << GetFlops() << endl;
      os << "Cumulative time   = " << GetTime() << endl;
      os << "MFlops rate       = " << 1.0e-6 * GetFlops() / GetTime() << endl;
      os << endl;
    }

    return(os);

  }

  // @}
  
private:

  //! Operator of which \c this object define the inverse.
  Operator Op_;
  //! Wrapper for IFPACK
  Teuchos::RefCountPtr<ML_Epetra::RowMatrix> RCPRowMatrix_;
  //! IFPACK preconditioner.
  Teuchos::RefCountPtr<Ifpack_Preconditioner> RCPData_;

  // @}
  
}; // InverseOperator

} // namespace MLAPI

#endif // ML_INVERSEOPERATOR_H
