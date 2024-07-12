// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PIRO_PERFORMSOLVE_HPP
#define PIRO_PERFORMSOLVE_HPP

#include "Thyra_ResponseOnlyModelEvaluatorBase.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_VectorBase.hpp"

#include "Piro_SolutionObserverBase.hpp"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"

namespace Piro {

//! \name Top-level Thyra solve drivers
//@{
//! \brief Evaluates the solved model and returns the first response.
//! \details .
//! \ingroup Piro_Thyra_solve_driver_grp
template <typename Scalar>
void PerformSolve(
    const Thyra::ResponseOnlyModelEvaluatorBase<Scalar> &piroModel,
    Teuchos::RCP<Thyra::VectorBase<Scalar> > &response);

//! \brief Evaluates the solved model and returns the specified response.
//! \details Returns the specified (first by default) response.
//! \ingroup Piro_Thyra_solve_driver_grp
template <typename Scalar>
void PerformSolve(
    const Thyra::ResponseOnlyModelEvaluatorBase<Scalar> &piroModel,
    Teuchos::ParameterList &solveParams,
    Teuchos::RCP<Thyra::VectorBase<Scalar> > &response);

//! \brief Evaluates the solved model and returns the specified response and sensitivity.
//! \details Returns the specified (first by default) response and optionally the corresponding sensitivity with respect to the first parameter.
//! \ingroup Piro_Thyra_solve_driver_grp
template <typename Scalar>
void PerformSolve(
    const Thyra::ResponseOnlyModelEvaluatorBase<Scalar> &piroModel,
    Teuchos::ParameterList &solveParams,
    Teuchos::RCP<Thyra::VectorBase<Scalar> > &response,
    Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > &sensitivity);

//! \brief Evaluates the solved model and returns specified responses and sensitivities.
//! \details Returns the requested responses and optionally the corresponding sensitivities with respect to all parameters.
//!          This version accepts pointers to non-<tt>const</tt> objects.
//! \ingroup Piro_Thyra_solve_driver_grp
template <typename Scalar>
void PerformSolve(
    const Thyra::ResponseOnlyModelEvaluatorBase<Scalar> &piroModel,
    Teuchos::ParameterList &solveParams,
    Teuchos::Array<Teuchos::RCP<Thyra::VectorBase<Scalar> > > &responses,
    Teuchos::Array<Teuchos::Array<Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > > > &sensitivities);

//! \brief Evaluates the solved model and returns specified responses and sensitivities.
//! \details Returns the requested responses and optionally the corresponding sensitivities with respect to all parameters.
//!          This version accepts pointers to <tt>const</tt>-qualified objects.
//! \ingroup Piro_Thyra_solve_driver_grp
template <typename Scalar>
void PerformSolve(
    const Thyra::ResponseOnlyModelEvaluatorBase<Scalar> &piroModel,
    Teuchos::ParameterList &solveParams,
    Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > > &responses,
    Teuchos::Array<Teuchos::Array<Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > > > &sensitivities);

//! \brief Evaluates the solved model and returns specified responses and sensitivities.
//! \details Returns the requested responses and optionally the corresponding sensitivities with respect to all parameters.
//!          This version accepts pointers to non-<tt>const</tt> objects.
//! \ingroup Piro_Thyra_solve_driver_grp
template <typename Scalar>
void PerformSolve(
    const Thyra::ResponseOnlyModelEvaluatorBase<Scalar> &piroModel,
    Teuchos::ParameterList &solveParams,
    Teuchos::Array<Teuchos::RCP<Thyra::VectorBase<Scalar> > > &responses,
    Teuchos::Array<Teuchos::Array<Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > > > &sensitivities,
    Teuchos::Array<Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > > &directions,
    Teuchos::Array<Teuchos::Array<Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > > > &reducedHessian);

//! \brief Evaluates the solved model and returns specified responses and sensitivities.
//! \details Returns the requested responses and optionally the corresponding sensitivities with respect to all parameters.
//!          This version accepts pointers to <tt>const</tt>-qualified objects.
//! \ingroup Piro_Thyra_solve_driver_grp
template <typename Scalar>
void PerformSolve(
    const Thyra::ResponseOnlyModelEvaluatorBase<Scalar> &piroModel,
    Teuchos::ParameterList &solveParams,
    Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > > &responses,
    Teuchos::Array<Teuchos::Array<Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > > > &sensitivities,
    Teuchos::Array<Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > > &directions,
    Teuchos::Array<Teuchos::Array<Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > > > &reducedHessian);
//@}

//! \name Other Thyra solve drivers
//! \brief The drivers do not statically check that the model is of the response-only variety.
//@{
//! \brief Evaluates the solved model and returns the first response.
//! \details Returns the first (i.e. with index 0) reponse.
//! \ingroup Piro_Thyra_solve_driver_grp
template <typename Scalar>
void PerformSolveBase(
    const Thyra::ModelEvaluator<Scalar> &piroModel,
    Teuchos::RCP<Thyra::VectorBase<Scalar> > &response);

//! \brief Evaluates the solved model and returns the specified response.
//! \details Returns the specified (first by default) response.
//! \ingroup Piro_Thyra_solve_driver_grp
template <typename Scalar>
void PerformSolveBase(
    const Thyra::ModelEvaluator<Scalar> &piroModel,
    Teuchos::ParameterList &solveParams,
    Teuchos::RCP<Thyra::VectorBase<Scalar> > &response);

//! \brief Evaluates the solved model and returns the specified response and sensitivity.
//! \details Returns the specified (first by default) response and optionally the corresponding sensitivity with respect to the first parameter.
//! \ingroup Piro_Thyra_solve_driver_grp
template <typename Scalar>
void PerformSolveBase(
    const Thyra::ModelEvaluator<Scalar> &piroModel,
    Teuchos::ParameterList &solveParams,
    Teuchos::RCP<Thyra::VectorBase<Scalar> > &response,
    Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > &sensitivity);

//! \brief Evaluates the solved model and returns specified responses and sensitivities.
//! \details Returns the requested responses and optionally the corresponding sensitivities with respect to all parameters.
//!          This version accepts pointers to non-<tt>const</tt> objects.
//! \ingroup Piro_Thyra_solve_driver_grp
template <typename Scalar>
void PerformSolveBase(
    const Thyra::ModelEvaluator<Scalar> &piroModel,
    Teuchos::ParameterList &solveParams,
    Teuchos::Array<Teuchos::RCP<Thyra::VectorBase<Scalar> > > &responses,
    Teuchos::Array<Teuchos::Array<Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > > > &sensitivities);

//! \brief Evaluates the solved model and returns specified responses and sensitivities.
//! \details Returns the requested responses and optionally the corresponding sensitivities with respect to all parameters.
//!          This version accepts pointers to <tt>const</tt>-qualified objects.
//! \ingroup Piro_Thyra_solve_driver_grp
template <typename Scalar>
void PerformSolveBase(
    const Thyra::ModelEvaluator<Scalar> &piroModel,
    Teuchos::ParameterList &solveParams,
    Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > > &responses,
    Teuchos::Array<Teuchos::Array<Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > > > &sensitivities);

//! \brief Evaluates the solved model and returns specified responses and sensitivities.
//! \details Returns the requested responses and optionally the corresponding sensitivities with respect to all parameters.
//!          This version accepts pointers to <tt>const</tt>-qualified objects.
//! \ingroup Piro_Thyra_solve_driver_grp
template <typename Scalar>
void PerformSolveBase(
    const Thyra::ModelEvaluator<Scalar> &piroModel,
    Teuchos::ParameterList &solveParams,
    Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > > &responses,
    Teuchos::Array<Teuchos::Array<Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > > > &sensitivities,
    Teuchos::RCP<SolutionObserverBase<Scalar, const Thyra::VectorBase<Scalar> > > observer);

//! \brief Evaluates the solved model and returns specified responses and sensitivities.
//! \details Returns the requested responses and optionally the corresponding sensitivities with respect to all parameters.
//!          This version accepts pointers to non-<tt>const</tt> objects.
//! \ingroup Piro_Thyra_solve_driver_grp
template <typename Scalar>
void PerformSolveBase(
    const Thyra::ModelEvaluator<Scalar> &piroModel,
    Teuchos::ParameterList &solveParams,
    Teuchos::Array<Teuchos::RCP<Thyra::VectorBase<Scalar> > > &responses,
    Teuchos::Array<Teuchos::Array<Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > > > &sensitivities,
    Teuchos::Array<Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > > &directions,
    Teuchos::Array<Teuchos::Array<Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > > > &reducedHessian);

//! \brief Evaluates the solved model and returns specified responses and sensitivities.
//! \details Returns the requested responses and optionally the corresponding sensitivities with respect to all parameters.
//!          This version accepts pointers to <tt>const</tt>-qualified objects.
//! \ingroup Piro_Thyra_solve_driver_grp
template <typename Scalar>
void PerformSolveBase(
    const Thyra::ModelEvaluator<Scalar> &piroModel,
    Teuchos::ParameterList &solveParams,
    Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > > &responses,
    Teuchos::Array<Teuchos::Array<Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > > > &sensitivities,
    Teuchos::Array<Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > > &directions,
    Teuchos::Array<Teuchos::Array<Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > > > &reducedHessian);

//! \brief Evaluates the solved model and returns specified responses and sensitivities.
//! \details Returns the requested responses and optionally the corresponding sensitivities with respect to all parameters.
//!          This version accepts pointers to <tt>const</tt>-qualified objects.
//! \ingroup Piro_Thyra_solve_driver_grp
template <typename Scalar>
void PerformSolveBase(
    const Thyra::ModelEvaluator<Scalar> &piroModel,
    Teuchos::ParameterList &solveParams,
    Teuchos::Array<Teuchos::RCP<const Thyra::VectorBase<Scalar> > > &responses,
    Teuchos::Array<Teuchos::Array<Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > > > &sensitivities,
    Teuchos::Array<Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > > &directions,
    Teuchos::Array<Teuchos::Array<Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > > > &reducedHessian,
    Teuchos::RCP<SolutionObserverBase<Scalar, const Thyra::VectorBase<Scalar> > > observer);
//@}

} // namespace Piro

#include "Piro_PerformSolve_Def.hpp"

#endif /* PIRO_PERFORMSOLVE_HPP */
