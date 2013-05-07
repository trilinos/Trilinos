// @HEADER
// ************************************************************************
//
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER

#ifndef PIRO_PERFORMSOLVE_HPP
#define PIRO_PERFORMSOLVE_HPP

//! \file Piro_PerformSolve.hpp
//! \brief Drivers for evaluating the responses and sensitivities of a solved model.

#include "Thyra_ResponseOnlyModelEvaluatorBase.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_VectorBase.hpp"

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
//@}

} // namespace Piro

#include "Piro_PerformSolve_Def.hpp"

#endif /* PIRO_PERFORMSOLVE_HPP */
