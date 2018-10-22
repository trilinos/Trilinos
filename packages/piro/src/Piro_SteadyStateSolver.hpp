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

#ifndef PIRO_STEADYSTATESOLVER_HPP
#define PIRO_STEADYSTATESOLVER_HPP

#include "Thyra_ResponseOnlyModelEvaluatorBase.hpp"

namespace Piro {

/** \brief Thyra-based abstract Model Evaluator for steady-states solves
 *  \ingroup Piro_Thyra_solver_grp
 * */
template <typename Scalar>
class SteadyStateSolver
    : public Thyra::ResponseOnlyModelEvaluatorBase<Scalar>
{
  public:

  /** \name Constructors/initializers */
  //@{
  /** \brief . */
  explicit SteadyStateSolver(const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > &model);

  /** \brief . */
  SteadyStateSolver(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > &model,
      int numParameters);
  //@}

  /** \name Overridden from Thyra::ModelEvaluatorBase . */
  //@{
  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;
  /** \brief . */
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_p_space(int l) const;
  /** \brief . */
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_g_space(int j) const;
  //@}

  /** \name Overridden from Thyra::ResponseOnlyModelEvaluatorBase . */
  //@{
  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;
  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<Scalar> getLowerBounds() const;
  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<Scalar> getUpperBounds() const;
  //@}

  protected:
  /** \name Service methods for subclasses. */
  //@{
  /** \brief . */
  const Thyra::ModelEvaluator<Scalar> &getModel() const;

  /** \brief . */
  int num_p() const;
  /** \brief . */
  int num_g() const;

  /** \brief . */
  void evalConvergedModel(
      const Thyra::ModelEvaluatorBase::InArgs<Scalar>& modelInArgs,
      const Thyra::ModelEvaluatorBase::OutArgs<Scalar>& outArgs) const;
  //@}

  private:
  /** \name Overridden from Thyra::ModelEvaluatorDefaultBase . */
  //@{
  /** \brief . */
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;

  /** \brief . */
  Teuchos::RCP<Thyra::LinearOpBase<Scalar> > create_DgDp_op_impl(int j, int l) const;
  //@}

  /** \name Internal implemention methods. */
  //@{
  /** \brief Implementation of createInArgs . */
  Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgsImpl() const;
  //@}

  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > model_;

  int num_p_;
  int num_g_;
};

}

#include "Piro_SteadyStateSolver_Def.hpp"
#endif /*PIRO_STEADYSTATESOLVER_HPP*/
