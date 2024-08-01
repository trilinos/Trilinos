// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
/*!
 * \file ml_nox_preconditioner.H
 *
 * \class ML_Nox_Preconditioner
 *
 * \brief ML nonlinear preconditioner and solver
 *
 * \date Last update do Doxygen: 31-Mar-05
 *
 */
// ML-headers
#include "ml_common.h"
#include "TrilinosCouplings_config.h"

#if defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_AMESOS) && defined(HAVE_ML_EPETRAEXT)

#ifndef ML_NOX_CONSTRAINEDMULTILEVELOPERATOR_H
#define ML_NOX_CONSTRAINEDMULTILEVELOPERATOR_H

#include <vector>

// epetra-headers
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Object.h"
#include "Epetra_Comm.h"
#include "Epetra_SerialComm.h"
#ifdef ML_MPI
#include "Epetra_MpiComm.h"
#endif
#include "Epetra_SerialComm.h"

// NOX-headers
#include "NOX_Epetra_Interface.H"

// ML-headers
#include "ml_common.h"
#include "ml_include.h"
#include "ml_epetra_utils.h"
#include "ml_epetra.h"
#include "ml_epetra_operator.h"
#include "ml_MultiLevelOperator.h"
#include "ml_agg_VBMETIS.h"
// utils
#include "ml_nox_preconditioner_utils.H"
// this is the wrapper for the Epetra_CrsGraph used by ml in matrix-free method
// it is derived from Epetra_RowMatrix so the ml - Epetra_RowMatrix - wrapper will take it
#include "ml_Epetra_wrap_CrsGraph_as_RowMatrix.H"

// interface to user code on fine level
#include "ml_nox_fineinterface.H"

// interface to user code on coarse level
#include "ml_nox_coarseinterface.H"



// ---------- Forward Declarations ----------

//! ML_NOX: default namespace for all ML/NOX nonlinear methods
namespace ML_NOX
{
class  ML_Nox_ConstrainedMultiLevelOperator : public virtual Epetra_Operator
{
public:

  

  explicit ML_Nox_ConstrainedMultiLevelOperator(
                           ML_Epetra::MultiLevelOperator* ml_operator,
                           ML_NOX::Nox_CoarseProblem_Interface& coarseinterface);

  
  virtual ~ML_Nox_ConstrainedMultiLevelOperator();


  //! Returns the label of this class.
  /*! 
   Returns the label of this class.
   */
  const char* Label() const {return(&(label_[0]));}

  //! get Comm of this class
  /*! 
   Derived from Epetra_Operator, returns ref to the Epetra_Comm of this class 
   */
  const Epetra_Comm & Comm() const {return(comm_);}

  //! Get fine level OperatorDomainMap.
  /*! 
   Derived from Epetra_Operator, get fine level OperatorDomainMap.
   */
  const Epetra_Map & OperatorDomainMap() const {return ml_operator_->OperatorDomainMap();}

  //! Get fine level OperatorRangeMap.
  /*! 
   Derived from Epetra_Operator, get fine level OperatorRangeMap.
   */
  const Epetra_Map & OperatorRangeMap() const {return ml_operator_->OperatorRangeMap();}
  
  //! ApplyInverse the preconditioner
  /*! 
   ApplyInverse the preconditioner. Method is derived from Epetra_Operator. 
   \param X   (In) : Epetra_MultiVector matching the fine level map of this
                     preconditioner
   \param Y (Out)  : Epetra_MultiVector containing the result on output
   */
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const ;

  //! Apply the preconditioner, not implemented.
  /*! 
   Derived from Epetra_Operator, not implemented. 
   */
  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const 
  { cout << "ML_Nox_ConstrainedMultiLevelOperator::Apply does not make sense for an ML-Multilevel Operator\n"
         << "Use ApplyInverse\n";throw -1;return(-1);}
    
  //! SetUseTranspose, not implemented.
  /*! 
   Derived from Epetra_Operator, not implemented. 
   */
  int SetUseTranspose(bool UseTranspose) { usetranspose_=UseTranspose; return 0;}

  //! NormInf, not implemented.
  /*! 
   Derived from Epetra_Operator, not implemented. 
   */
  double NormInf() const {cout << "ML_Nox_ConstrainedMultiLevelOperator::NormInf not impl.\n";
                          throw -1;return(-1.0);}


  //! UseTranspose, not implemented.
  /*! 
   Derived from Epetra_Operator, not implemented. 
   */
  bool UseTranspose() const {return usetranspose_;}

  //! HasNormInf, not implemented.
  /*! 
   Derived from Epetra_Operator, always returns false. 
   */
  bool HasNormInf() const {return false;}

private:
  // don't want copy-ctor and = operator
  ML_Nox_ConstrainedMultiLevelOperator(ML_Nox_ConstrainedMultiLevelOperator& old);
  ML_Nox_ConstrainedMultiLevelOperator operator = (const ML_Nox_ConstrainedMultiLevelOperator& old);
  
  const char*             label_;             // the label of this class
  bool                                  usetranspose_;
  const Epetra_Comm&                    comm_;
  ML_NOX::Nox_CoarseProblem_Interface&  coarseinterface_;
  ML_Epetra::MultiLevelOperator*        ml_operator_;

}; // class  ML_Nox_ConstrainedMultiLevelOperator : public virtual Epetra_Operator
}  // namespace ML_NOX


#endif // ML_NOX_CONSTRAINEDMULTILEVELOPERATOR_H
#endif // defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) 
