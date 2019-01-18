/*
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
// Questions? Contact Roscoe A. Bartlett (bartlettra@ornl.gov) 
// 
// ***********************************************************************
// @HEADER
*/


#ifndef DUMMY_TEST_MODEL_EVALUATOR_DECL_HPP
#define DUMMY_TEST_MODEL_EVALUATOR_DECL_HPP


#include "Thyra_ModelEvaluatorDefaultBase.hpp"


namespace Thyra {

// Mock Extended InArgs and OutArgs Objects. In practice, these
// objects will be defined by a solver that needs to extend the
// InArgs/OutArgs for specialized data without cluttering the core
// model evaluator interface.
template<class Scalar>
struct MockExtendedInArgs
{
  Teuchos::RCP<Thyra::VectorBase<Scalar> > a;
};

template<class Scalar>
struct MockExtendedOutArgs
{
  Teuchos::RCP<Thyra::VectorBase<Scalar> > b;
};


template<class Scalar> class DummyTestModelEvaluator;


/** \brief Nonmember constuctor.
 *
 * \relates DummyTestModelEvaluator
 */
template<class Scalar>
RCP<DummyTestModelEvaluator<Scalar> >
dummyTestModelEvaluator(
  const Ordinal x_size = 2,
  const ArrayView<const Ordinal> &p_sizes = Teuchos::null,
  const ArrayView<const Ordinal> &g_sizes = Teuchos::null,
  const bool supports_x_dot = false,
  const bool supports_x_dot_dot = false,
  const bool supports_extended_inargs = true,
  const bool supports_extended_outargs = true
  );


/** \brief Test helper ModelEvaluator.
 *
 * This class is used to help unit test the various ModelEvaluator support
 * software, that is it.
 */
template<class Scalar>
class DummyTestModelEvaluator : public ModelEvaluatorDefaultBase<Scalar>
{
public:

  /** \name Initializers/Accessors */
  //@{

  /** \brief . */
  DummyTestModelEvaluator(
    const Ordinal x_size,
    const ArrayView<const Ordinal> &p_sizes,
    const ArrayView<const Ordinal> &g_sizes,
    const bool supports_x_dot = false,
    const bool supports_x_dot_dot = false,
    const bool supports_extended_inargs = true,
    const bool supports_extended_outargs = true
    );

  //@}

  /** \name Public functions overridden from ModelEvaulator. */
  //@{

  /** \brief . */
  RCP<const VectorSpaceBase<Scalar> > get_x_space() const;
  /** \brief . */
  RCP<const VectorSpaceBase<Scalar> > get_p_space(int l) const;
  /** \brief . */
  RCP<const Teuchos::Array<std::string> > get_p_names(int l) const;
  /** \brief . */
  RCP<const VectorSpaceBase<Scalar> > get_f_space() const;
  /** \brief . */
  RCP<const VectorSpaceBase<Scalar> > get_g_space(int j) const;
  /** \brief . */
  Teuchos::ArrayView<const std::string> get_g_names(int j) const;
  /** \brief . */
  ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;
  /** \brief . */
  ModelEvaluatorBase::InArgs<Scalar> getLowerBounds() const;
  /** \brief . */
  ModelEvaluatorBase::InArgs<Scalar> getUpperBounds() const;
  /** \brief . */
  RCP<LinearOpBase<Scalar> > create_W_op() const;
  /** \brief . */
  RCP<PreconditionerBase<Scalar> > create_W_prec() const;
  /** \brief . */
  RCP<const LinearOpWithSolveFactoryBase<Scalar> > get_W_factory() const;
  /** \brief . */
  ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;
  /** \brief . */
  void reportFinalPoint(
    const ModelEvaluatorBase::InArgs<Scalar> &finalPoint,
    const bool wasSolved
    );

  //@}

  // For unit testing
  void change_p_size_incorrectly(const Ordinal new_size);
  void change_p_size_correctly(const Ordinal new_size);

private: // functions

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

private: // data members

  RCP<const VectorSpaceBase<Scalar> > x_space_;
  Array<RCP<const VectorSpaceBase<Scalar> > > p_space_;
  RCP<const VectorSpaceBase<Scalar> > f_space_;
  Array<RCP<const VectorSpaceBase<Scalar> > > g_space_;
  Array<std::string> g_names_;
  RCP<const LinearOpWithSolveFactoryBase<Scalar> > W_factory_;
  ModelEvaluatorBase::InArgs<Scalar> nominalValues_;
  RCP<VectorBase<Scalar> > x0_;
  ModelEvaluatorBase::InArgs<Scalar> prototypeInArgs_;
  ModelEvaluatorBase::OutArgs<Scalar> prototypeOutArgs_;

};


} // namespace Thyra


#endif // DUMMY_TEST_MODEL_EVALUATOR_DECL_HPP
