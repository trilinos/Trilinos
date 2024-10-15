#ifndef ML_LBOPERATOR_H
#define ML_LBOPERATOR_H

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

/*!
\file MLAPI_LoadBalanceOperator.h

\brief wraps an MLAPI operator with zero rows on some processors.

\author Michael Gee, TU Munich.

\date Last updated on July-10.
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
//#include "ml_common.h"
//#include <iostream>
//#include "ml_operator.h"
//#include "ml_epetra_utils.h"
//#include "ml_RowMatrix.h"
//#include "Teuchos_RefCountPtr.hpp"
//#include "MLAPI_Error.h"
//#include "MLAPI_Space.h"
//#include "MLAPI_MultiVector.h"
#include "MLAPI_Operator.h"
//#include "MLAPI_CompObject.h"
//#include "MLAPI_TimeObject.h"
//#include "MLAPI_Workspace.h"
//#include "MLAPI_Operator_Box.h"

namespace MLAPI {

/*!
 * \class Operator
 *
 * \brief Operator: basic class to define operators within MLAPI.
 *
 * \author Michael Gee, TU Munich.
 *
 * \date Last updated on July-10.
 */

class LoadBalanceOperator : public Operator {

public:
  //@{ \name Constructors and destructors.

  //! Default constructor.
  LoadBalanceOperator()
  {
    participate_ = true;
    mlcomm_ = NULL;
    RCPOperatorBox_ = Teuchos::null;
  }

  //! Constructor with given already computed ML_Operator pointer.
  LoadBalanceOperator(const Space& DomainSpace, const Space& RangeSpace,
           ML_Operator* Op, bool Ownership = true,
           Teuchos::RefCountPtr<ML_Operator_Box> AuxOp = Teuchos::null)
  {
    participate_ = true;
    mlcomm_ = NULL;
    Reshape(DomainSpace, RangeSpace, Op, Ownership, AuxOp);
  }

  //! Constructor with given already FillComplete()'d object.
  LoadBalanceOperator(const Space& DomainSpace, const Space& RangeSpace,
           Epetra_RowMatrix* Matrix, bool Ownership = true,
           Teuchos::RefCountPtr<ML_Operator_Box> AuxOp = Teuchos::null)
  {
    if (Matrix) participate_ = true;
    else        participate_ = false;
    mlcomm_ = NULL;
    Reshape(DomainSpace, RangeSpace, Matrix, Ownership, AuxOp);
  }

  //! Copy constructor.
  LoadBalanceOperator(const LoadBalanceOperator& RHS)
  {
    participate_       = RHS.participate_;
    if (RHS.mlcomm_ && participate_)
    {
      ML_Comm_Create(&mlcomm_);
      ML_Comm_Set_UsrComm(mlcomm_,RHS.mlcomm_->USR_comm);
    }
    DomainSpace_       = RHS.GetDomainSpace();
    RangeSpace_        = RHS.GetRangeSpace();
    ColumnSpace_       = RHS.GetColumnSpace();
    RCPOperatorBox_    = RHS.GetRCPOperatorBox();
    RCPAuxOperatorBox_ = RHS.GetRCPAuxOperatorBox();
    RCPRowMatrix_      = RHS.GetRCPRowMatrix();

    SetLabel(RHS.GetLabel());
  }


  //! Destructor.
  ~LoadBalanceOperator()
  {
    Destroy();
  }

  // @}
  // @{ \name Reshape methods

  //! Resets \c this object.
  void Reshape()
  {
    Destroy();
  }

  //! Reshape with given already computed ML_Operator pointer.
  void Reshape(const Space& DomainSpace, const Space& RangeSpace,
               ML_Operator* Op, bool Ownership = true,
               Teuchos::RefCountPtr<ML_Operator_Box> AuxOp = Teuchos::null)
  {
    StackPush();

    RangeSpace_        = RangeSpace;
    DomainSpace_       = DomainSpace;

    RCPOperatorBox_    = Teuchos::rcp(new ML_Operator_Box(Op,Ownership));
    RCPRowMatrix_      = Teuchos::rcp(new ML_Epetra::RowMatrix(Op,&(GetEpetra_Comm()),
                                                               false));
    RCPAuxOperatorBox_ = AuxOp;

    StackPop();
  }

  //! Reshape with given already FillComplete()'d object.
  void Reshape(const Space& DomainSpace, const Space& RangeSpace,
               Epetra_RowMatrix* Matrix, bool Ownership = true,
               Teuchos::RCP<ML_Operator_Box> AuxOp = Teuchos::null)
  {
    StackPush();
    if (!Matrix) participate_ = false;
    else         participate_ = true;

    // these are outer spaces, all participate in their comm
    RangeSpace_ = RangeSpace;
    DomainSpace_ = DomainSpace;


#ifdef ML_MPI
    if (participate_)
    {
      // create an ML_Comm that still contains MPI_COMM_WORLD
      ML_Comm_Create(&mlcomm_);
      // have to use the comm from Matrix for the subset of procs here
      MPI_Comm mpicomm = dynamic_cast<const Epetra_MpiComm&>(Matrix->Comm()).GetMpiComm();
      ML_Comm_Set_UsrComm(mlcomm_,mpicomm);
      ML_Operator* Op = ML_Operator_Create(mlcomm_);

      RCPOperatorBox_ = Teuchos::rcp(new ML_Operator_Box(Op,true));
      RCPAuxOperatorBox_ = AuxOp;
      RCPRowMatrix_ = Teuchos::rcp(Matrix,Ownership);
      ML_Operator_WrapEpetraMatrix(RCPRowMatrix_.get(), GetML_Operator());
    }
#else
    ML_Operator* Op = ML_Operator_Create(MLAPI::GetML_Comm());
    RCPOperatorBox_ = Teuchos::rcp(new ML_Operator_Box(Op,true));
    RCPAuxOperatorBox_ = AuxOp;

    RCPRowMatrix_ = Teuchos::rcp(Matrix,Ownership);
    ML_Operator_WrapEpetraMatrix(RCPRowMatrix_.get(), GetML_Operator());
#endif

    StackPop();
  }

  // @}
  // @{ \name Overloaded operators

  //! Makes \c this object equivalent to \c RHS.
  LoadBalanceOperator& operator=(const LoadBalanceOperator& RHS)
  {
    StackPush();

    Destroy();

    participate_       = RHS.participate_;
    if (RHS.mlcomm_ && participate_)
    {
      ML_Comm_Create(&mlcomm_);
      ML_Comm_Set_UsrComm(mlcomm_,RHS.mlcomm_->USR_comm);
    }
    DomainSpace_    = RHS.GetDomainSpace();
    RangeSpace_     = RHS.GetRangeSpace();
    ColumnSpace_    = RHS.GetColumnSpace();
    RCPOperatorBox_ = RHS.GetRCPOperatorBox();
    RCPRowMatrix_   = RHS.GetRCPRowMatrix();

    SetLabel(RHS.GetLabel());

    StackPop();

    return(*this);
  }

  //! Sets the label of \c this object.
  inline Operator& operator=(const std::string& Label)
  {
    SetLabel(Label);
    return(*this);
  }

  // @}
  // @{ \name Get and Set methods

  //! Returns a bool indicating whether this proc participates in the operator application
  virtual inline bool GetParticipation() const {
    return(participate_);
  }

  //! Returns a reference to the internally stored domain space.
  const Space GetOperatorDomainSpace() const {
    return(DomainSpace_);
  }

  //! Returns a reference to the internally stored range space.
  const Space GetOperatorRangeSpace() const {
    return(RangeSpace_);
  }

  //! Returns a reference to the internally stored domain space.
  inline const Space GetDomainSpace() const {
    return(DomainSpace_);
  }

  //! Returns a reference to the internally stored range space.
  inline const Space GetRangeSpace() const {
    return(RangeSpace_);
  }

  //! Returns a reference to the internally stored column space.
  inline const Space GetColumnSpace() const
  {
    return(ColumnSpace_);
  }

  //! Returns the number of global rows.
  inline int GetNumGlobalRows() const
  {
    return(GetRangeSpace().GetNumGlobalElements());
  }

  //! Returns the number of local rows.
  inline int GetNumMyRows() const
  {
    return(GetRangeSpace().GetNumMyElements());
  }

  //! Returns the number of global columns.
  inline int GetNumGlobalCols() const
  {
    if (participate_) return(GetRowMatrix()->NumGlobalCols());
    else              return(0);
  }

  //! Returns the number of local columns.
  inline int GetNumMyCols() const
  {
    if (participate_) return(GetRowMatrix()->NumMyCols());
    else              return(0);
  }

  //! Returns the global number of nonzeros.
  inline int GetNumGlobalNonzeros() const
  {
    if (participate_) return(GetRowMatrix()->NumGlobalNonzeros());
    else              return(0);
  }

  //! Returns the local number of nonzeros.
  inline int GetNumMyNonzeros() const
  {
    if (participate_) return(GetRowMatrix()->NumMyNonzeros());
    else              return(0);
  }

  //! Returns the RefCountPtr of OperatorBox_.
  inline const Epetra_RowMatrix* GetRowMatrix() const
  {
    if (participate_) return(RCPRowMatrix_.get());
    else              return NULL;
  }

  //! Returns the RefCountPtr of OperatorBox_.
  inline ML_Operator* GetML_Operator() const
  {
    if (participate_) return(GetRCPOperatorBox()->GetData());
    else              return NULL;
  }

  //! Returns the RefCountPtr of OperatorBox_.
  inline const Teuchos::RCP<ML_Operator_Box>& GetRCPOperatorBox() const
  {
    return(RCPOperatorBox_);
  }

  //! Returns the RefCountPtr of AuxOperatorBox_.
  inline const Teuchos::RCP<ML_Operator_Box>& GetRCPAuxOperatorBox() const
  {
    return(RCPAuxOperatorBox_);
  }
  //! Returns the RefCountPtr of RowMatrix_
  inline const Teuchos::RCP<Epetra_RowMatrix>& GetRCPRowMatrix() const
  {
    return(RCPRowMatrix_);
  }

  //! Returns the global ID of local row ID \c LRID.
  int GetGRID(const int LRID) const
  {
#ifdef MLAPI_CHECK
    if (LRID < 0 || LRID >= GetNumMyRows())
      ML_THROW("LRID in invalid", -1);
#endif
    return(GetRangeSpace()(LRID));
  }

  //! Returns the global ID of local column ID \c LCID.
  int GetGCID(const int LCID) const
  {
#ifdef MLAPI_CHECK
    if (LCID < 0 || LCID >= GetRowMatrix()->NumMyCols())
      ML_THROW("LRID in invalid", -1);
#endif
    return(GetRowMatrix()->RowMatrixColMap().GID(LCID));
  }

  // @}
  // @{ \name Mathematical methods.

  //! Applies \c this operator to LHS, returns the result in \c RHS.
  int Apply(const MultiVector& X, MultiVector& Y) const
  {
    ResetTimer();
    StackPush();

    if (GetDomainSpace() != X.GetVectorSpace())
      ML_THROW("Domain spaces differ", -1);
    if (GetRangeSpace() != Y.GetVectorSpace())
      ML_THROW("Range spaces differ", -1);
    if (X.GetNumVectors() != Y.GetNumVectors())
      ML_THROW("Number of vectors differ", -1);

    if (participate_)
    {

      if (GetML_Operator() == 0)
        ML_THROW("Operator not set", -1);

      int (*func)(ML_Operator*,int,double*,int,double*) =
        GetML_Operator()->matvec->func_ptr;

      for (int v = 0 ; v < X.GetNumVectors() ; ++v) {
        double* x_ptr = (double*)&X(0) + v * X.GetMyLength();
        double* y_ptr = (double*)&Y(0) + v * Y.GetMyLength();
        (*func)(GetML_Operator(),X.GetMyLength(),x_ptr,
                Y.GetMyLength(), y_ptr);
      }
    }

    StackPop();

    UpdateFlops(2.0 * GetNumGlobalNonzeros());
    UpdateTime();

    return(0);
  }

  // @}
  // @{ \name Miscellaneous methods

  //! Prints basic information about \c this object.
  std::ostream& Print(std::ostream& os, const bool verbose = true) const
  {
    ML_THROW("Print(...) not implemented in MLAPI::LoadBalanceOperator",-1);
    if (GetRCPOperatorBox().get() == 0) {
      if (GetMyPID() == 0) {
        os << std::endl;
        os << "*** MLAPI::Operator ***" << std::endl;
        os << "Label  = " << GetLabel() << std::endl;
        os << "Status = empty" << std::endl;
        os << std::endl;
      }
      return(os);
    }

    StackPush();

    int    *bindx;
    double *val;
    int    allocated, row_length;
    ML_Operator* matrix = GetML_Operator();

    if (matrix->getrow == NULL)
      ML_THROW("getrow not set", -1);

    if (GetMyPID() == 0) {
      os << std::endl;
      os << "*** MLAPI::Operator ***" << std::endl;
      os << "Label             = " << GetLabel() << std::endl;
      os << "Number of rows    = " << GetRangeSpace().GetNumGlobalElements() << std::endl;
      os << "Number of columns = " << GetDomainSpace().GetNumGlobalElements() << std::endl;
      os << "Flop count        = " << GetFlops() << std::endl;
      os << "Cumulative time   = " << GetTime() << std::endl;
      if (GetTime() != 0.0)
        os << "MFlops rate       = " << 1.0e-6 * GetFlops() / GetTime() << std::endl;
      else
        os << "MFlops rate       = 0.0" << std::endl;
      os << std::endl;
    }

    if (!verbose)
      return(os);

    allocated = 100;
    bindx = (int    *)  ML_allocate(allocated*sizeof(int   ));
    val   = (double *)  ML_allocate(allocated*sizeof(double));

    if (GetMyPID() == 0) {
      os.width(10);
      os << "ProcID";
      os.width(20);
      os << "Global Row";
      os.width(20);
      os << "Global Col";
      os.width(20);
      os << "Value" << std::endl;
      os << std::endl;
    }

    for (int iproc = 0 ; iproc < GetNumProcs() ; ++iproc) {

      if (GetMyPID() == iproc) {

        for (int i = 0 ; i < matrix->getrow->Nrows; i++) {
          ML_get_matrix_row(matrix, 1, &i, &allocated, &bindx, &val,
                            &row_length, 0);
          for  (int j = 0; j < row_length; j++) {
            int GlobalRow = GetRangeSpace()(i);
            //int GlobalCol = GetColumnSpace()(bindx[j]);
            int GlobalCol = GetRowMatrix()->RowMatrixColMap().GID(bindx[j]);
            os.width(10);
            os << iproc;
            os.width(20);
            os << GlobalRow;
            os.width(20);
            os << GlobalCol;
            os.width(20);
            os << val[j] << std::endl;
          }
        }
      }
      Barrier();
    }

    if (GetMyPID() == 0)
      os << std::endl;

    Barrier();

    ML_free(val);
    ML_free(bindx);

    StackPop();

    return (os);
  }

  //! Build the column space, by computing the GID of all local columns.
  void BuildColumnSpace()
  {
    StackPush();

    if (participate_)
    {
      if (mlcomm_->ML_nprocs == 1) {
        ColumnSpace_ = DomainSpace_;
        return;
      }

      std::vector<double> dtemp;
      std::vector<int> GlobalElements;

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
    }

    StackPop();

    return;
  }

  // @}

private:

  //! Destroys all internal data and resets \c this object.
  void Destroy()
  {
    participate_ = true;
    if (mlcomm_) ML_Comm_Destroy(&mlcomm_);
    RangeSpace_.Reshape();
    DomainSpace_.Reshape();
    RCPOperatorBox_    = Teuchos::null;
    RCPRowMatrix_      = Teuchos::null;
    RCPAuxOperatorBox_ = Teuchos::null;
  }

  //! indicate whether I participate in this operator or not
  bool participate_;
  ML_Comm* mlcomm_;
  //! Domain space.
  Space DomainSpace_;
  //! Range space.
  Space RangeSpace_;
  //! Column space.
  Space ColumnSpace_;
  //! Container for the underlying ML_Operator pointer.
  Teuchos::RCP<ML_Operator_Box> RCPOperatorBox_;
  //! Container for the underlying ML_Operator pointer.
  Teuchos::RCP<ML_Operator_Box> RCPAuxOperatorBox_;
  //! Container for the underlying Epetra_RowMatrix pointer
  Teuchos::RCP<Epetra_RowMatrix> RCPRowMatrix_;

}; // Operator

} // namespace MLAPI

#endif // ML_OPERATOR_H
