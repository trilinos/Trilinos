#ifndef ML_OPERATOR_H
#define ML_OPERATOR_H
#include "ml_include.h"
#include <iostream>
#include "ml_operator.h"
#include "ml_epetra.h"
#include "ml_amesos.h"
#include "ml_epetra_utils.h"
#include "ml_amesos_wrap.h"
#include "ml_RowMatrix.h"
#ifdef HAVE_ML_ANASAZI
#include "ml_anasazi.h"
#endif
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "MLAPI_Space.h"
#include "MLAPI_MultiVector.h"
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
    RCPOperatorBox_ = Teuchos::null;
  }

  //! Constructor with given already computed ML_Operator pointer.
  Operator(const Space& DomainSpace, const Space& RangeSpace,
           ML_Operator* Op, bool Ownership = true)
  {
    Reshape(DomainSpace, RangeSpace, Op, Ownership);
  }

  //! Constructor with given already FillComplete()'d object.
  Operator(const Space& DomainSpace, const Space& RangeSpace,
           Epetra_RowMatrix* Matrix, bool Ownership = true)
  {
    Reshape(DomainSpace, RangeSpace, Matrix, Ownership);
  }

  //! Copy constructor.
  Operator(const Operator& RHS) 
  {
    DomainSpace_    = RHS.GetDomainSpace();
    RangeSpace_     = RHS.GetRangeSpace();
    ColumnSpace_    = RHS.GetColumnSpace();
    RCPOperatorBox_ = RHS.GetRCPOperatorBox();
    RCPRowMatrix_   = RHS.GetRCPRowMatrix();
    
    SetLabel(RHS.GetLabel());
  }

  //! Destructor.
  ~Operator()
  {
    Destroy();
  }

  // @{
  // @} Reshape methods

  //! Reshape with given already computed ML_Operator pointer.
  void Reshape(const Space& DomainSpace, const Space& RangeSpace,
               ML_Operator* Op, bool Ownership = true)
  {
    RangeSpace_ = RangeSpace;
    DomainSpace_ = DomainSpace;

    RCPOperatorBox_ = Teuchos::rcp(new ML_Operator_Box(Op,Ownership));
    if (RangeSpace == DomainSpace) // FIXME !!!!!!!
      RCPRowMatrix_ = Teuchos::rcp(new ML_Epetra::RowMatrix(Op,&(GetEpetra_Comm())),Ownership);
    else
      ML_THROW("FIX BUG IN ml_RowMatrix (for rect matrices)", -1);

    BuildColumnSpace();
  }

  //! Reshape with given already FillComplete()'d object.
  void Reshape(const Space& DomainSpace, const Space& RangeSpace,
               Epetra_RowMatrix* Matrix, bool Ownership = true)
  {
    RangeSpace_ = RangeSpace;
    DomainSpace_ = DomainSpace;


    ML_Operator* Op = ML_Operator_Create(MLAPI::GetML_Comm());
    RCPOperatorBox_ = Teuchos::rcp(new ML_Operator_Box(Op,true));

    RCPRowMatrix_ = Teuchos::rcp(Matrix,Ownership);
    Epetra2MLMatrix(RCPRowMatrix_.get(), GetML_Operator());

    BuildColumnSpace();

  }

  // @}
  // @{ Overloaded operators

  //! Makes \c this object equivalent to \c RHS.
  Operator& operator=(const Operator& RHS) 
  {
    Destroy();

    DomainSpace_ = RHS.GetDomainSpace();
    RangeSpace_  = RHS.GetRangeSpace();
    ColumnSpace_ = RHS.GetColumnSpace();
    RCPOperatorBox_ = RHS.GetRCPOperatorBox();
    RCPRowMatrix_   = RHS.GetRCPRowMatrix();
    
    SetLabel(RHS.GetLabel());
    return(*this);
  }

  //! Sets the label of \c this object.
  inline Operator& operator=(const string& Label)
  {
    SetLabel(Label);
    return(*this);
  }

  // @}
  // @{ Query methods
  
  //! Returns a reference to the internally stored domain space.
  inline const Space& GetDomainSpace() const {
    return(DomainSpace_);
  }

  //! Returns a reference to the internally stored range space.
  inline const Space& GetRangeSpace() const {
    return(RangeSpace_);
  }

  //! Returns a reference to the internally stored column space.
  inline const Space& GetColumnSpace() const 
  {
    return(ColumnSpace_);
  }

  // @}
  // @{ Mathematical methods.
  
  //! Applies \c this operator to LHS, returns the result in \c RHS.
  int Apply(const MultiVector& LHS, MultiVector& RHS) const
  {
    assert (GetDomainSpace() == LHS.GetVectorSpace());
    assert (GetRangeSpace() == RHS.GetVectorSpace());
    assert (GetML_Operator() != 0);

    int DomainSize = GetDomainSpace().GetNumMyElements();
    int RangeSize  = GetRangeSpace().GetNumMyElements();
    
    int (*func)(ML_Operator*,int,double*,int,double*) = GetML_Operator()->matvec->func_ptr;

    (*func)(GetML_Operator(),DomainSize,(double*)&LHS(0),
            RangeSize,(double*)&RHS(0));
    return(0);
  }

  //! Returns the RefCountPtr of OperatorBox_.
  inline const Teuchos::RefCountPtr<ML_Operator_Box>& GetRCPOperatorBox() const
  {
    return(RCPOperatorBox_);
  }

  //! Returns the RefCountPtr of RowMatrix_
  inline const Teuchos::RefCountPtr<Epetra_RowMatrix> GetRCPRowMatrix() const
  {
    return(RCPRowMatrix_);
  }

  //! Returns the RefCountPtr of OperatorBox_.
  inline const Epetra_RowMatrix* GetRowMatrix() const
  {
    return(RCPRowMatrix_.get());
  }
  
  //! Returns the RefCountPtr of OperatorBox_.
  inline ML_Operator* GetML_Operator() const
  {
    return(GetRCPOperatorBox()->GetData());
  }

  //! Prints basic information about \c this object.
  ostream& Print(std::ostream& os, const bool verbose = true) const
  {
    int    *bindx;
    double *val;
    int    allocated, row_length;
    ML_Operator* matrix = GetML_Operator();

    if (matrix->getrow == NULL) 
      throw("getrow not set");

    if (GetMyPID() == 0) {
      os << endl;
      os << "*** MLAPI::Operator ***" << endl;
      os << "Label             = " << GetLabel() << endl;
      os << "Number of rows    = " << GetRangeSpace().GetNumGlobalElements() << endl;
      os << "Number of columns = " << GetDomainSpace().GetNumGlobalElements() << endl;
      os << endl;
    }

    if (!verbose) 
      return(os);

    allocated = 100;
    bindx = (int    *)  ML_allocate(allocated*sizeof(int   ));
    val   = (double *)  ML_allocate(allocated*sizeof(double));

    for (int iproc = 0 ; iproc < GetNumProcs() ; ++iproc) {

      if (GetMyPID() == 0) {
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

      if (GetMyPID() == iproc) {

        for (int i = 0 ; i < matrix->getrow->Nrows; i++) {
          ML_get_matrix_row(matrix, 1, &i, &allocated, &bindx, &val,
                            &row_length, 0);
          for  (int j = 0; j < row_length; j++) {
            int GlobalRow = GetDomainSpace()(i);
            int GlobalCol = GetColumnSpace()(bindx[j]);
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

    int Nrows = GetML_Operator()->getrow->Nrows;
    int Nghosts;
    if (GetML_Operator()->getrow->pre_comm == NULL) Nghosts = 0;
    else {
      if (GetML_Operator()->getrow->pre_comm->total_rcv_length <= 0)
        ML_CommInfoOP_Compute_TotalRcvLength(GetML_Operator()->getrow->pre_comm);
      Nghosts = GetML_Operator()->getrow->pre_comm->total_rcv_length;
    }

    dtemp.resize(Nrows + Nghosts);

    for (int i = 0 ; i < Nrows ; ++i) 
      dtemp[i] = 1.0 * GetDomainSpace()(i);
    for (int i = 0 ; i < Nghosts; ++i) 
      dtemp[i + Nrows] = -1;

    ML_exchange_bdry(&dtemp[0],GetML_Operator()->getrow->pre_comm,
                     GetML_Operator()->outvec_leng,
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
  Teuchos::RefCountPtr<ML_Operator_Box> RCPOperatorBox_;
  //! Container for the underlying Epetra_RowMatrix pointer
  Teuchos::RefCountPtr<Epetra_RowMatrix> RCPRowMatrix_;

  // @}

}; // Operator

} // namespace MLAPI
#endif // ML_OPERATOR_H

