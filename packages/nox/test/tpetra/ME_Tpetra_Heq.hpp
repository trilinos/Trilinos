// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#ifndef NOX_TPETRA_ME_HEQ_DECL_HPP
#define NOX_TPETRA_ME_HEQ_DECL_HPP

#include "Thyra_StateFuncModelEvaluatorBase.hpp"
#include "Tpetra_CrsMatrix.hpp"

template<class Scalar, class LO, class GO, class Node>
class EvaluatorTpetraHeq;

template<class Scalar, class LO, class GO, class Node>
class HeqJacobianOperator;

/** \brief Nonmember constuctor.
 *
 * \relates EvaluatorTpetraHeq
 */
template<class Scalar, class LO, class GO, class Node>
Teuchos::RCP<EvaluatorTpetraHeq<Scalar, LO, GO, Node> >
evaluatorTpetraHeq(const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                   const Tpetra::global_size_t numGlobalElements,
                   const Scalar omega);


/** \brief Chandrasehkar H-equation
 *
 * The equation modeled is:

 \verbatim

   F(x)_i = (1 - \frac{\omega}{2N}\sum_{j=1}^N \frac{\mu_ix_j}{\mu_i + \mu_j})^{-1} - x_i = 0, 0 <= i <= N-1

   where \omega \in [0,1] is a given parameter and:
      \mu_i = (i + 0.5)/N

 \endverbatim

 * The Jacobian for this system is dense, so this creates a matrix-free Jacobian operator class.
 * The preconditioner this computes is the inverse diagonal of the Jacobian.
 */
template<class Scalar, class LO, class GO, class Node>
class EvaluatorTpetraHeq
  : public Thyra::StateFuncModelEvaluatorBase<Scalar>
{
public:

  // Public typedefs
  typedef Scalar scalar_type;
  typedef LO local_ordinal_type;
  typedef GO global_ordinal_type;
  typedef Node node_type;
  typedef HeqJacobianOperator<Scalar, LO, GO, Node> jac_op;
  typedef Tpetra::Map<LO, GO, Node> tpetra_map;
  typedef Tpetra::Vector<Scalar, LO, GO, Node> tpetra_vec;
  typedef Tpetra::Operator<Scalar, LO, GO, Node> tpetra_op;
  typedef Tpetra::CrsMatrix<Scalar, LO, GO, Node> tpetra_mat;
  typedef Thyra::VectorSpaceBase<Scalar> thyra_vec_space;
  typedef Thyra::VectorBase<Scalar> thyra_vec;
  typedef Thyra::LinearOpBase<Scalar> thyra_op;
  typedef Thyra::PreconditionerBase<Scalar> thyra_prec;

  // Constructor
  EvaluatorTpetraHeq(const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                     const Tpetra::global_size_t numGlobalElements,
                     const Scalar omega);

  /** \name Initializers/Accessors */
  //@{

  /** \brief . */
  void setShowGetInvalidArgs(bool showGetInvalidArg);

  void set_W_factory(const Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> >& W_factory);

  //@}

  /** \name Public functions overridden from ModelEvaulator. */
  //@{

  /** \brief . */
  Teuchos::RCP<const thyra_vec_space> get_x_space() const;
  /** \brief . */
  Teuchos::RCP<const thyra_vec_space> get_f_space() const;
  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;
  /** \brief . */
  Teuchos::RCP<thyra_op> create_W_op() const;
  /** \brief . */
  Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> > get_W_factory() const;
  /** \brief . */
  Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;
  /** \brief . */
  Teuchos::RCP<thyra_prec> create_W_prec() const;
  //@}

private:

  /** \name Private functions overridden from ModelEvaulatorDefaultBase. */
  //@{

  /** \brief . */
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;
  /** \brief . */
  void evalModelImpl(
    const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
    const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs
    ) const;

  //@}

private: // data members

  const Teuchos::RCP<const Teuchos::Comm<int> >  comm_;
  const Tpetra::global_size_t numGlobalElements_;
  std::vector<std::size_t> procNumElements_;
  std::vector<GO> procMinGIDs_;

  Teuchos::RCP<const thyra_vec_space> xSpace_;
  Teuchos::RCP<const tpetra_map>   xMap_;

  Teuchos::RCP<const thyra_vec_space> fSpace_;
  Teuchos::RCP<const tpetra_map>   fMap_;

  Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> > W_factory_;

  const Scalar omega_;
  Teuchos::RCP<tpetra_vec> integralOp_;

  Thyra::ModelEvaluatorBase::InArgs<Scalar> nominalValues_;
  Teuchos::RCP<thyra_vec> x0_;
  bool showGetInvalidArg_;
  Thyra::ModelEvaluatorBase::InArgs<Scalar> prototypeInArgs_;
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> prototypeOutArgs_;

  mutable Teuchos::RCP<Teuchos::Time> residTimer_;
  mutable Teuchos::RCP<Teuchos::Time> intOpTimer_;
};

template<class Scalar, class LO, class GO, class Node>
class HeqJacobianOperator
  : public Tpetra::Operator<Scalar, LO, GO, Node>
{
public:

  // Public typedefs
  typedef Scalar scalar_type;
  typedef LO local_ordinal_type;
  typedef GO global_ordinal_type;
  typedef Node node_type;
  typedef Tpetra::Map<LO, GO, Node> tpetra_map;
  typedef Tpetra::Vector<Scalar, LO, GO, Node> tpetra_vec;
  typedef Tpetra::MultiVector<Scalar, LO, GO, Node> tpetra_mv;

  // Constructor
  HeqJacobianOperator(const Teuchos::RCP<const tpetra_map>& map,
                      const std::vector<std::size_t>& procNumElements,
                      const std::vector<GO>& procMinGIDs);

  void initialize(const Scalar& omega,
                  const Teuchos::RCP<tpetra_vec>& integralOp);

  void unitialize();

  virtual Teuchos::RCP<const tpetra_map> getDomainMap () const;

  virtual Teuchos::RCP<const tpetra_map> getRangeMap () const;

  virtual void
  apply (const tpetra_mv& X,
         tpetra_mv& Y,
         Teuchos::ETransp mode = Teuchos::NO_TRANS,
         Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
         Scalar beta = Teuchos::ScalarTraits<Scalar>::zero()) const;

private:

  Teuchos::RCP<const tpetra_map> map_;
  std::vector<std::size_t> procNumElements_;
  std::vector<GO> procMinGIDs_;
  Scalar omega_;
  Teuchos::RCP<tpetra_vec> integralOp_;
  Teuchos::RCP<tpetra_vec> integralOpX_;

};

//==================================================================
#include "ME_Tpetra_Heq_def.hpp"
//==================================================================

#endif
