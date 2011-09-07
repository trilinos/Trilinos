// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_DEFAULT_FINITE_DIFFERENCE_MODEL_EVALUATOR_DECL_HPP
#define THYRA_DEFAULT_FINITE_DIFFERENCE_MODEL_EVALUATOR_DECL_HPP

#include "Thyra_ModelEvaluatorDelegatorBase.hpp"
#include "Thyra_DirectionalFiniteDiffCalculator.hpp"
#include "Teuchos_StandardMemberCompositionMacros.hpp"
#include "Teuchos_StandardCompositionMacros.hpp"


namespace Thyra {


/** \brief This class wraps any ModelEvaluator object and computes certain
 * derivatives using finite differences.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup Thyra_Nonlin_ME_support_grp
 */
template<class Scalar>
class DefaultFiniteDifferenceModelEvaluator
  : virtual public ModelEvaluatorDelegatorBase<Scalar>
{
public:

  /** \name Constructors/initializers/accessors/utilities. */
  //@{

  /** \brief Utility object that computes directional finite differences */
  STANDARD_COMPOSITION_MEMBERS(
    DirectionalFiniteDiffCalculator<Scalar>, direcFiniteDiffCalculator );

  /** \brief . */
  DefaultFiniteDifferenceModelEvaluator();

  /** \brief . */
  void initialize(
    const RCP<ModelEvaluator<Scalar> > &thyraModel
    ,const RCP<DirectionalFiniteDiffCalculator<Scalar> > &direcFiniteDiffCalculator
    );

  //@}

  /** \name Public functions overridden from Teuchos::Describable. */
  //@{

  /** \brief . */
  std::string description() const;

  //@}

private:

  /** \name Private functions overridden from ModelEvaulatorDefaultBase. */
  //@{

  /** \brief . */
  ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;
  /** \brief . */
  void evalModelImpl(
    const ModelEvaluatorBase::InArgs<Scalar> &inArgs,
    const ModelEvaluatorBase::OutArgs<Scalar> &outArgs
    ) const;

  //@}
 
};


/** \brief Nonmember constructor.
 *
 * \relates DefaultFiniteDifferenceModelEvaluator
 */
template<class Scalar>
RCP<DefaultFiniteDifferenceModelEvaluator<Scalar> >
defaultFiniteDifferenceModelEvaluator()
{
  return Teuchos::rcp(new DefaultFiniteDifferenceModelEvaluator<Scalar>());
}


/** \brief Nonmember constructor.
 *
 * \relates DefaultFiniteDifferenceModelEvaluator
 */
template<class Scalar>
RCP<DefaultFiniteDifferenceModelEvaluator<Scalar> >
defaultFiniteDifferenceModelEvaluator(
  const RCP<ModelEvaluator<Scalar> > &thyraModel,
  const RCP<DirectionalFiniteDiffCalculator<Scalar> > &direcFiniteDiffCalculator
  )
{
  RCP<DefaultFiniteDifferenceModelEvaluator<Scalar> > fdModel =
    defaultFiniteDifferenceModelEvaluator<Scalar>();
  fdModel->initialize(thyraModel, direcFiniteDiffCalculator);
  return fdModel;
}


} // namespace Thyra


#endif // THYRA_DEFAULT_FINITE_DIFFERENCE_MODEL_EVALUATOR_DECL_HPP
