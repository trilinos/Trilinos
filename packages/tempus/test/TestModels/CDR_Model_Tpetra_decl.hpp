//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef TEMPUS_CDR_MODEL_TPETRA_DECL_HPP
#define TEMPUS_CDR_MODEL_TPETRA_DECL_HPP

#include "Thyra_StateFuncModelEvaluatorBase.hpp"
#include "Thyra_TpetraThyraWrappers_decl.hpp"

#include <Tpetra_CrsGraph_decl.hpp>
#include <Tpetra_Import_decl.hpp>
#include <Tpetra_Map_decl.hpp>
#include <Tpetra_Vector_decl.hpp>

namespace Tempus_Test {

/** \brief 1D CGFEM model for convection/diffusion/reaction
 *
 * The equation modeled is:

 \verbatim

dT     dT   d^2(T)
-- + a -- + ------ - K * T**2 = 0
dt     dz    dz^2

   subject to:
      T = 1.0 @ z = z_min

 \endverbatim

 * The Matrix <tt>W = d(f)/d(x)</tt> is implemented as a
 * <tt>Thyra::MultiVectorBase</tt> object and the class
 * <tt>Thyra::DefaultSerialDenseLinearOpWithSolveFactory</tt> is used to
 * create the linear solver.
 */
template <typename SC, typename LO, typename GO, typename Node>
class CDR_Model_Tpetra : public ::Thyra::StateFuncModelEvaluatorBase<SC> {
 public:
  using tpetra_map    = Tpetra::Map<LO, GO, Node>;
  using tpetra_graph  = Tpetra::CrsGraph<LO, GO, Node>;
  using tpetra_matrix = Tpetra::CrsMatrix<SC, LO, GO, Node>;
  using tpetra_vec    = Tpetra::Vector<SC, LO, GO, Node>;
  using tpetra_extract =
      ::Thyra::TpetraOperatorVectorExtraction<SC, LO, GO, Node>;

  CDR_Model_Tpetra(const Teuchos::RCP<const Teuchos::Comm<int>> &comm,
                   const GO numGlobalElements, const SC zMin, const SC zMax,
                   const SC a,   // convection
                   const SC k);  // source

  /** \name Initializers/Accessors */
  //@{

  void set_x0(const Teuchos::ArrayView<const SC> &x0);

  void setShowGetInvalidArgs(bool showGetInvalidArg);

  void set_W_factory(
      const Teuchos::RCP<const ::Thyra::LinearOpWithSolveFactoryBase<SC>>
          &W_factory);

  //@}

  /** \name Public functions overridden from ModelEvaluator. */
  //@{

  Teuchos::RCP<const ::Thyra::VectorSpaceBase<SC>> get_x_space() const;
  Teuchos::RCP<const ::Thyra::VectorSpaceBase<SC>> get_f_space() const;
  ::Thyra::ModelEvaluatorBase::InArgs<SC> getNominalValues() const;
  Teuchos::RCP<Thyra::LinearOpWithSolveBase<double>> create_W() const;
  Teuchos::RCP<::Thyra::LinearOpBase<SC>> create_W_op() const;
  Teuchos::RCP<const ::Thyra::LinearOpWithSolveFactoryBase<SC>> get_W_factory()
      const;
  ::Thyra::ModelEvaluatorBase::InArgs<SC> createInArgs() const;
  Teuchos::RCP<::Thyra::PreconditionerBase<SC>> create_W_prec() const;
  //@}

 private:
  /** Allocates and returns the Jacobian matrix graph */
  virtual Teuchos::RCP<const Tpetra::CrsGraph<LO, GO, Node>> createGraph();

  /** \name Private functions overridden from ModelEvaluatorDefaultBase. */
  //@{

  ::Thyra::ModelEvaluatorBase::OutArgs<SC> createOutArgsImpl() const;
  void evalModelImpl(
      const ::Thyra::ModelEvaluatorBase::InArgs<SC> &inArgs,
      const ::Thyra::ModelEvaluatorBase::OutArgs<SC> &outArgs) const;

  //@}

 private:  // data members
  const Teuchos::RCP<const Teuchos::Comm<int>> comm_;
  const int numGlobalElements_;
  const SC zMin_;
  const SC zMax_;
  const SC a_;
  const SC k_;

  Teuchos::RCP<const ::Thyra::VectorSpaceBase<SC>> xSpace_;
  Teuchos::RCP<const tpetra_map> xOwnedMap_;
  Teuchos::RCP<const tpetra_map> xGhostedMap_;
  Teuchos::RCP<const Tpetra::Import<LO, GO, Node>> importer_;

  Teuchos::RCP<const ::Thyra::VectorSpaceBase<SC>> fSpace_;
  Teuchos::RCP<const tpetra_map> fOwnedMap_;

  Teuchos::RCP<const tpetra_graph> wGraph_;

  Teuchos::RCP<const ::Thyra::LinearOpWithSolveFactoryBase<SC>> wFactory_;

  Teuchos::RCP<tpetra_vec> nodeCoordinates_;
  Teuchos::RCP<tpetra_vec> GhostedNodeCoordinates_;

  mutable Teuchos::RCP<tpetra_vec> uPtr_;
  mutable Teuchos::RCP<tpetra_vec> uDotPtr_;
  mutable Teuchos::RCP<tpetra_vec> xPtr_;

  mutable Teuchos::RCP<tpetra_vec> jDiag_;

  ::Thyra::ModelEvaluatorBase::InArgs<SC> nominalValues_;
  Teuchos::RCP<::Thyra::VectorBase<SC>> x0_;
  Teuchos::Array<SC> p_;
  bool showGetInvalidArg_;
  ::Thyra::ModelEvaluatorBase::InArgs<SC> prototypeInArgs_;
  ::Thyra::ModelEvaluatorBase::OutArgs<SC> prototypeOutArgs_;
};

}  // namespace Tempus_Test

#endif  // TEMPUS_CDR_MODEL_TPETRA_DECL_HPP
