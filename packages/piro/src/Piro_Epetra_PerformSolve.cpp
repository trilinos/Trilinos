// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Piro_Epetra_PerformSolve.hpp"

#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Map.h"

#include "Piro_PerformSolve.hpp"

#include "Thyra_EpetraModelEvaluator.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"

#include <iterator>

namespace Piro {

namespace Epetra {

namespace Detail {

template <typename VectorType, typename MultiVectorType>
void PerformSolveImpl(
    const EpetraExt::ModelEvaluator &piroSolver,
    Teuchos::ParameterList &solveParams,
    Teuchos::Array<Teuchos::RCP<VectorType> > &responses,
    Teuchos::Array<Teuchos::Array<Teuchos::RCP<MultiVectorType> > > &sensitivities)
{
  const Teuchos::RCP<const EpetraExt::ModelEvaluator> epetraSolver = Teuchos::rcpFromRef(piroSolver);
  const Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactory = Teuchos::null;
  const Thyra::EpetraModelEvaluator thyraSolver(epetraSolver, lowsFactory);

  typedef Teuchos::Array<Teuchos::RCP<Thyra::VectorBase<double> > > ThyraResponseArray;
  ThyraResponseArray thyraResponses;

  typedef Teuchos::Array<Teuchos::Array<Teuchos::RCP<Thyra::MultiVectorBase<double> > > > ThyraSensitivityArray;
  ThyraSensitivityArray thyraSensitivities;

  PerformSolveBase(thyraSolver, solveParams, thyraResponses, thyraSensitivities);

  responses.clear();
  responses.reserve(thyraResponses.size());
  for (ThyraResponseArray::const_iterator it_begin = thyraResponses.begin(),
      it_end = thyraResponses.end(),
      it = it_begin;
      it != it_end;
      ++it) {
    const int j = std::distance(it_begin, it);
    const Epetra_Map g_map = *piroSolver.get_g_map(j);
    const Teuchos::RCP<Thyra::VectorBase<double> > g_thyra = *it;
    const Teuchos::RCP<VectorType> g =
      Teuchos::nonnull(g_thyra) ? Thyra::get_Epetra_Vector(g_map, g_thyra) : Teuchos::null;
    responses.push_back(g);
  }

  sensitivities.clear();
  sensitivities.resize(thyraSensitivities.size());
  for (ThyraSensitivityArray::const_iterator it_begin = thyraSensitivities.begin(),
      it_end = thyraSensitivities.end(),
      it = it_begin;
      it != it_end;
      ++it) {
    const int j = std::distance(it_begin, it);
    const Epetra_Map g_map = *piroSolver.get_g_map(j);
    for (ThyraSensitivityArray::value_type::const_iterator jt = it->begin(), jt_end = it->end(); jt != jt_end; ++jt) {
      const Teuchos::RCP<Thyra::MultiVectorBase<double> > dgdp_thyra = *jt;
      const Teuchos::RCP<MultiVectorType> dgdp =
        Teuchos::nonnull(dgdp_thyra) ? Thyra::get_Epetra_MultiVector(g_map, dgdp_thyra) : Teuchos::null;
      sensitivities[j].push_back(dgdp);
    }
  }
}

} // namespace Detail

void PerformSolve(
    const EpetraExt::ModelEvaluator &piroSolver,
    Teuchos::ParameterList &solveParams,
    Teuchos::Array<Teuchos::RCP<Epetra_Vector> > &responses,
    Teuchos::Array<Teuchos::Array<Teuchos::RCP<Epetra_MultiVector> > > &sensitivities)
{
  Detail::PerformSolveImpl(piroSolver, solveParams, responses, sensitivities);
}

void PerformSolve(
    const EpetraExt::ModelEvaluator &piroSolver,
    Teuchos::ParameterList &solveParams,
    Teuchos::Array<Teuchos::RCP<const Epetra_Vector> > &responses,
    Teuchos::Array<Teuchos::Array<Teuchos::RCP<const Epetra_MultiVector> > > &sensitivities)
{
  Detail::PerformSolveImpl(piroSolver, solveParams, responses, sensitivities);
}

} // namespace Epetra

} // namespace Piro
