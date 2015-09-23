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

#ifndef PIRO_LOCASOLVER_HPP
#define PIRO_LOCASOLVER_HPP

#include "Piro_SteadyStateSolver.hpp"

#include "Piro_ObserverBase.hpp"

#include "LOCA.H"
#include "LOCA_Thyra.H"
#include "LOCA_Thyra_SaveDataStrategy.H"

#include "Teuchos_ParameterList.hpp"

namespace Piro {

/** \brief Thyra-based Model Evaluator for LOCA solves
 *  \ingroup Piro_Thyra_solver_grp
 * */
template <typename Scalar>
class LOCASolver : public SteadyStateSolver<Scalar> {
public:
  /** \name Constructor/Destructor */
  //@{
  /** \brief Constructs a LOCASolver instance given a model and optionally a data saving strategy . */
  LOCASolver(
      const Teuchos::RCP<Teuchos::ParameterList> &piroParams,
      const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &model,
      const Teuchos::RCP<LOCA::Thyra::SaveDataStrategy> &saveDataStrategy);

  ~LOCASolver();
  //@}

private:
  /** \name Overridden from Thyra::ModelEvaluatorDefaultBase . */
  //@{
  /** \brief . */
  void evalModelImpl(
      const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
      const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs) const;
  //@}

  Teuchos::RCP<Teuchos::ParameterList> piroParams_;
  Teuchos::RCP<LOCA::Thyra::SaveDataStrategy> saveDataStrategy_;

  Teuchos::RCP<LOCA::GlobalData> globalData_;
  mutable LOCA::ParameterVector paramVector_;
  Teuchos::RCP<LOCA::Thyra::Group> group_;
  Teuchos::RCP<LOCA::StatusTest::Abstract> locaStatusTests_;
  Teuchos::RCP<NOX::StatusTest::Generic> noxStatusTests_;

  Teuchos::RCP<LOCA::Stepper> stepper_;
  mutable bool first_;
};


template <typename Scalar>
Teuchos::RCP<LOCASolver<Scalar> >
observedLocaSolver(
    const Teuchos::RCP<Teuchos::ParameterList> &piroParams,
    const Teuchos::RCP<Thyra::ModelEvaluator<Scalar> > &model,
    const Teuchos::RCP<Piro::ObserverBase<Scalar> > &observer);

} // namespace Piro

#endif /* PIRO_LOCASOLVER_HPP */
