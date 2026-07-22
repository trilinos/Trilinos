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
 * \file nlnml_prepostoperator.H
 *
 * \class NLNML_PrePostOperator
 *
 * \brief class NLNML_PrePostOperator used by the nonlinear preconditioner
 *
 * \date Last update do Doxygen: 31-Mar-05
 *
 */
#include "ml_common.h"
#include "ml_include.h"
#include "TrilinosCouplings_config.h"

#if defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_AMESOS) && defined(HAVE_ML_EPETRAEXT)
#ifndef NLNML_PREPOSTOPERATOR_H
#define NLNML_PREPOSTOPERATOR_H

#include <iostream>

// Epetra headers
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_Comm.h"
#include "Epetra_SerialComm.h"
#ifdef ML_MPI
#include "Epetra_MpiComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"

// NOX headers
#include "NOX_Abstract_PrePostOperator.H"
#include "NOX_Epetra_Group.H"
#include "NOX_Solver_Generic.H"

// own headers
#include "nlnml_finelevelnoxinterface.H"
#include "nlnml_coarselevelnoxinterface.H"

// forward declarations
class Epetra_CrsMatrix;

using namespace Teuchos;
namespace NLNML
{


class NLNML_PrePostOperator : public virtual NOX::Abstract::PrePostOperator
{
public:

  NLNML_PrePostOperator(
              RefCountPtr<NLNML::NLNML_CoarseLevelNoxInterface> coarseinterface,
              RefCountPtr<NLNML::NLNML_FineLevelNoxInterface> finterface);

  NLNML_PrePostOperator(NLNML_PrePostOperator& old);
  NLNML_PrePostOperator(const NLNML::NLNML_PrePostOperator& old);

  virtual ~NLNML_PrePostOperator();

  void runPreSolve(const NOX::Solver::Generic& solver);

  const std::string& getType() const { return type_; }

private:

  NLNML_PrePostOperator operator = (const NLNML_PrePostOperator& old);

protected:

  std::string                                       type_;
  RefCountPtr<NLNML::NLNML_FineLevelNoxInterface>   fineinterface_;
  RefCountPtr<NLNML::NLNML_CoarseLevelNoxInterface> coarseinterface_;

}; // class NLNML_PREPOSTOPERATOR_H
} // namespace NLNML
#endif
#endif
