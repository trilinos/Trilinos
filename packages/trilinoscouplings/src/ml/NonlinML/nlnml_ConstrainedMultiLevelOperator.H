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
 * \file nlnml_ConstrainedMultiLevelOperator.H
 *
 * \class NLNML_ConstrainedMultiLevelOperator
 *
 * \brief a wrapper for ML's MultiLevelOperator that additionally enforces constraints
 *
 * \date Last update do Doxygen: 31-Mar-05
 *
 */
#include "TrilinosCouplings_config.h"

#if defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_AMESOS) && defined(HAVE_ML_EPETRAEXT)

#ifndef NLNML_CONSTRAINEDMULTILEVELOPERATOR_H
#define NLNML_CONSTRAINEDMULTILEVELOPERATOR_H

// ML-headers
#include "ml_common.h"
#include "ml_include.h"
#include "ml_epetra.h"
#include "ml_epetra_operator.h"

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

// Teuchos
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"

#include "nlnml_coarselevelnoxinterface.H"


using namespace Teuchos;
namespace NLNML
{

class NLNML_ConstrainedMultiLevelOperator : public virtual Epetra_Operator
{
public:

  explicit NLNML_ConstrainedMultiLevelOperator(
              RefCountPtr<ML_Epetra::MultiLevelOperator>        ml_operator,
              RefCountPtr<NLNML::NLNML_CoarseLevelNoxInterface> coarseinterface,
              bool applyconstraints);

  virtual ~NLNML_ConstrainedMultiLevelOperator() {}

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
  { std::cout << "NLNML_ConstrainedMultiLevelOperator::Apply does not make sense for an ML-Multilevel Operator\n"
              << "Use ApplyInverse\n";throw -1;return(-1);}

  //! SetUseTranspose, not implemented.
  /*!
   Derived from Epetra_Operator, not implemented.
   */
  int SetUseTranspose(bool UseTranspose) { return 0;}

  //! NormInf, not implemented.
  /*!
   Derived from Epetra_Operator, not implemented.
   */
  double NormInf() const {std::cout << "NLNML_ConstrainedMultiLevelOperator::NormInf not impl.\n";
                          return(-1.0);}

  //! UseTranspose, not implemented.
  /*!
   Derived from Epetra_Operator, not implemented.
   */
  bool UseTranspose() const {return false;}

  //! HasNormInf, not implemented.
  /*!
   Derived from Epetra_Operator, always returns false.
   */
  bool HasNormInf() const {return false;}

private:

  explicit NLNML_ConstrainedMultiLevelOperator(NLNML::NLNML_ConstrainedMultiLevelOperator& old);
  NLNML_ConstrainedMultiLevelOperator operator = (const NLNML_ConstrainedMultiLevelOperator& old);

  const char*                                       label_;
  const Epetra_Comm&                                comm_;
  RefCountPtr<NLNML::NLNML_CoarseLevelNoxInterface> coarseinterface_;
  RefCountPtr<ML_Epetra::MultiLevelOperator>        ml_operator_;
  bool                                              applyconstraints_;

}; // class NLNML_ConstrainedMultiLevelOperator
} // namespace NLNML

#endif
#endif
