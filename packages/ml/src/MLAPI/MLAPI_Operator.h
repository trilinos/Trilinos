#ifndef ML_OPERATOR_H
#define ML_OPERATOR_H
#include "ml_include.h"
#include <iostream>
#include "ml_operator.h"
#include "ml_epetra.h"
#include "ml_amesos.h"
#include "ml_epetra_utils.h"
#include "ml_amesos_wrap.h"
#ifdef HAVE_ML_ANASAZI
#include "ml_anasazi.h"
#endif
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "MLAPI_Space.h"
#include "MLAPI_DoubleVector.h"
#include "MLAPI_Workspace.h"
#include "MLAPI_Operator_Box.h"
#include "Epetra_MpiComm.h"
#include "Epetra_CrsMatrix.h"

using namespace std;

namespace Teuchos {
  class ParameterList;
}

namespace MLAPI {

#ifdef HAVE_ML_EPETRAz
const int MatrixType = ML_EpetraCRS_MATRIX;
#else
const int MatrixType = ML_CSR_MATRIX;
#endif

/*!
 * \class Operator
 *
 * \brief Operator: basic class to define operators within MLAPI.
 *
 * \author Marzio Sala, SNL 9214
 *
 * \date Last updated on 07-Jan-05
 */

class Operator : public BaseObject {

public:
  //@{ Constructors and destructors.

  //! Default constructor.
  Operator() 
  {
    OperatorBox_ = Teuchos::null;
  }

  //! Constructor with given already computed ML_Operator pointer.
  Operator(const Space& DomainSpace, const Space& RangeSpace,
           ML_Operator* Op, bool Ownership = true)
  {
    OperatorBox_ = Teuchos::rcp(new ML_Operator_Box(Op,Ownership));
    RangeSpace_ = RangeSpace;
    DomainSpace_ = DomainSpace;
    BuildColumnSpace();
  }

  //! Constructor with given already FillComplete()'d object.
  Operator(const Space& DomainSpace, const Space& RangeSpace,
           Epetra_RowMatrix& Matrix)
  {
    RangeSpace_ = RangeSpace;
    DomainSpace_ = DomainSpace;

    ML_Operator* Op = ML_Operator_Create(MLAPI::GetML_Comm());
    OperatorBox_ = Teuchos::rcp(new ML_Operator_Box(Op,true));
    Epetra2MLMatrix(&Matrix, OperatorBox_->GetData());
    BuildColumnSpace();

  }

  //! Constructor with given already FillComplete()'d object.
  Operator(const Space& DomainSpace, const Space& RangeSpace,
           Epetra_RowMatrix* Matrix)
  {
    RangeSpace_ = RangeSpace;
    DomainSpace_ = DomainSpace;

    RowMatrix_ = Teuchos::rcp(Matrix);

    ML_Operator* Op = ML_Operator_Create(MLAPI::GetML_Comm());
    OperatorBox_ = Teuchos::rcp(new ML_Operator_Box(Op,true));
    Epetra2MLMatrix(RowMatrix_.get(), OperatorBox_->GetData());
    BuildColumnSpace();

  }

  //! Copy constructor.
  Operator(const Operator& RHS) 
  {
    DomainSpace_ = RHS.DomainSpace();
    RangeSpace_  = RHS.RangeSpace();
    ColumnSpace_ = RHS.ColumnSpace();
    OperatorBox_ = RHS.OperatorBox();
    RowMatrix_   = RHS.RowMatrix();
    
    SetLabel(RHS.GetLabel());
  }

  // Destructor.
  ~Operator()
  {
    Destroy();
  }

  // @}
  // @{ Overloaded operators

  //! Makes \c this object equivalent to \c RHS.
  Operator& operator=(const Operator& RHS) 
  {
    Destroy();

    DomainSpace_ = RHS.DomainSpace();
    RangeSpace_  = RHS.RangeSpace();
    ColumnSpace_ = RHS.ColumnSpace();
    OperatorBox_ = RHS.OperatorBox();
    RowMatrix_   = RHS.RowMatrix();
    
    SetLabel(RHS.GetLabel());
    return(*this);
  }

  //! Sets the label of \c this object.
  Operator& operator=(const string& Label)
  {
    SetLabel(Label);
    return(*this);
  }

  // @}
  // @{ Query methods
  
  //! Returns a reference to the internally stored domain space.
  const Space& DomainSpace() const {
    return(DomainSpace_);
  }

  //! Returns a reference to the internally stored range space.
  const Space& RangeSpace() const {
    return(RangeSpace_);
  }

  //! Returns a reference to the internally stored column space.
  const Space& ColumnSpace() const 
  {
    return(ColumnSpace_);
  }

  //! Returns a pointer to the internally stored ML_Operator struct.
  ML_Operator* GetData() const
  {
    return(OperatorBox_->GetData());
  }
  
  // @}
  // @{ Mathematical methods.
  
  //! Applies \c this operator to LHS, returns the result in \c RHS.
  int Apply(const DoubleVector& LHS, DoubleVector& RHS) const
  {
    assert (DomainSpace() == LHS.VectorSpace());
    assert (RangeSpace() == RHS.VectorSpace());
    assert (GetData() != 0);

    int DomainSize = DomainSpace().NumMyElements();
    int RangeSize  = RangeSpace().NumMyElements();
    
    int (*func)(ML_Operator*,int,double*,int,double*) = GetData()->matvec->func_ptr;

    (*func)(GetData(),DomainSize,(double*)&LHS(0),
            RangeSize,(double*)&RHS(0));
    return(0);
  }

  void ComputeEigenValues(const string Type, double* Er, double* Ei, 
                          double* V) const
  {
    if (Type == "LAPACK")
    {
      int ierr;
      ierr = ML_Operator_Eigensolver_Dense(GetData(), Er, Ei, V);
    }
    else {
      cerr << "ERROR: In CompareEigenValues()" << endl;
      cerr << "ERROR: (file " << __FILE__ << ", line " << __LINE__ << ")" << endl;
      cerr << "ERROR: Requested type (" << Type
           << ") not supported" << endl;
      throw(-1);
    }

    return;
  }
  
  //! Returns the RefCountPtr of OperatorBox_.
  const Teuchos::RefCountPtr<ML_Operator_Box>& OperatorBox() const
  {
    return(OperatorBox_);
  }

  //! Returns the RefCountPtr of OperatorBox_.
  const Teuchos::RefCountPtr<Epetra_RowMatrix>& RowMatrix() const
  {
    return(RowMatrix_);
  }

  //! Prints basic information about \c this object.
  ostream& Print(std::ostream& os, const bool verbose = true) const
  {
    int    *bindx;
    double *val;
    int    allocated, row_length;
    ML_Operator* matrix = GetData();

    if (matrix->getrow == NULL) 
      throw("getrow not set");

    if (MyPID() == 0) {
      os << endl;
      os << "*** MLAPI::Operator ***" << endl;
      os << "Label             = " << GetLabel() << endl;
      os << "Number of rows    = " << RangeSpace().NumGlobalElements() << endl;
      os << "Number of columns = " << DomainSpace().NumGlobalElements() << endl;
      os << endl;
    }

    if (!verbose) 
      return(os);

    allocated = 100;
    bindx = (int    *)  ML_allocate(allocated*sizeof(int   ));
    val   = (double *)  ML_allocate(allocated*sizeof(double));

    for (int iproc = 0 ; iproc < NumProc() ; ++iproc) {

      if (MyPID() == 0) {
        os.width(10);
        os << "ProcID";
        os.width(20);
        os << "Global Row";
        os.width(20);
        os << "Global Col";
        os.width(20);
        os << "Value" << endl;
        os << endl;
      }

      if (MyPID() == iproc) {

        for (int i = 0 ; i < matrix->getrow->Nrows; i++) {
          ML_get_matrix_row(matrix, 1, &i, &allocated, &bindx, &val,
                            &row_length, 0);
          for  (int j = 0; j < row_length; j++) {
            int GlobalRow = DomainSpace()(i);
            int GlobalCol = ColumnSpace()(bindx[j]);
            os.width(10);
            os << iproc;
            os.width(20);
            os << GlobalRow;
            os.width(20);
            os << GlobalCol;
            os.width(20);
            os << val[j] << endl;
          }
        }
      }
      Barrier();
    }

    ML_free(val);
    ML_free(bindx);
    return (os);
  }

private:
  
  // @}
  // @{ Internally used methods.
  
  //! Build the column space, by computing the GID of all local columns.
  void BuildColumnSpace()
  {

    vector<double> dtemp;
    vector<int> GlobalElements;

    int Nrows = GetData()->getrow->Nrows;
    int Nghosts;
    if (GetData()->getrow->pre_comm == NULL) Nghosts = 0;
    else {
      if (GetData()->getrow->pre_comm->total_rcv_length <= 0)
        ML_CommInfoOP_Compute_TotalRcvLength(GetData()->getrow->pre_comm);
      Nghosts = GetData()->getrow->pre_comm->total_rcv_length;
    }

    dtemp.resize(Nrows + Nghosts);

    for (int i = 0 ; i < Nrows ; ++i) 
      dtemp[i] = 1.0 * DomainSpace()(i);
    for (int i = 0 ; i < Nghosts; ++i) 
      dtemp[i + Nrows] = -1;

    ML_exchange_bdry(&dtemp[0],GetData()->getrow->pre_comm,
                     GetData()->outvec_leng,
                     GetML_Comm(), ML_OVERWRITE,NULL);

    GlobalElements.resize(Nrows + Nghosts);

    for (int i = 0 ; i < Nrows + Nghosts ; ++i)
      GlobalElements[i] = (int)dtemp[i];

    ColumnSpace_.Reshape(-1, Nrows + Nghosts, &GlobalElements[0]);

    return;
  }

  //! Destroys all internal data.
  void Destroy() { }

  // @}
  // @{ Internal data

  //! Domain space.
  Space DomainSpace_;
  //! Range space.
  Space RangeSpace_;
  //! Column space.
  Space ColumnSpace_;
  //! Container for the underlying ML_Operator pointer.
  Teuchos::RefCountPtr<ML_Operator_Box> OperatorBox_;
  //! Container for the underlying Epetra_RowMatrix pointer
  Teuchos::RefCountPtr<Epetra_RowMatrix> RowMatrix_;

  // @}

}; // Operator

} // namespace MLAPI
#endif // ML_OPERATOR_H

