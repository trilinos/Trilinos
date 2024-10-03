// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Piro_Epetra_PerformAnalysis.hpp"

#include "Epetra_Vector.h"
#include "Epetra_Map.h"

#include "Piro_PerformAnalysis.hpp"

#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraModelEvaluator.hpp"

int Piro::Epetra::PerformAnalysis(
    EpetraExt::ModelEvaluator &piroModel,
    Teuchos::ParameterList &analysisParams,
    Teuchos::RCP<Epetra_Vector> &p)
{
  const Teuchos::RCP<const EpetraExt::ModelEvaluator> epetraModel = Teuchos::rcpFromRef(piroModel);
  const Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactory = Teuchos::null;
  Thyra::EpetraModelEvaluator thyraModel(epetraModel, lowsFactory);

  Teuchos::RCP<Thyra::VectorBase<double> > p_thyra;
  const int status = ::Piro::PerformAnalysis(thyraModel, analysisParams, p_thyra);

  if (Teuchos::nonnull(p_thyra)) {
    const int l = 0;
    const Teuchos::RCP<const Epetra_Map> p_map = epetraModel->get_p_map(l);
    p = Thyra::get_Epetra_Vector(*p_map, p_thyra);
  } else {
    p = Teuchos::null;
  }

  return status;
}
