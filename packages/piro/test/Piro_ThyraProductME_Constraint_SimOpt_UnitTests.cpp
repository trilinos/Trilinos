// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"

#include "Piro_ThyraProductME_Constraint_SimOpt.hpp"
#include "MockModelEval_A_Tpetra.hpp"

#include "ROL_ThyraVector.hpp"
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_DefaultModelEvaluatorWithSolveFactory.hpp"
#include "Thyra_ModelEvaluatorDelegatorBase.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_RowMatrixTransposer.hpp"

namespace {

using Teuchos::RCP;
using Teuchos::rcp;

/* Adjoint model that, as the adjoint models of Albany and MockModelEval_A_Tpetra, evaluates the Jacobian of the
 * underlying (forward) model in the W operator and then replaces the matrix of W with its transpose.
 * The transpose is stored with a reversed column map, so an operator that holds a transpose cannot be filled again by
 * the underlying model, which assumes the column map of its Jacobian graph. Hence the W operator must be created
 * anew for each evaluation. */
class TransposingAdjointModel : public Thyra::ModelEvaluatorDelegatorBase<double> {
public:
  explicit TransposingAdjointModel(const RCP<Thyra::ModelEvaluator<double>>& model) :
    Thyra::ModelEvaluatorDelegatorBase<double>(model) {}

private:
  void evalModelImpl(const Thyra::ModelEvaluatorBase::InArgs<double>& inArgs,
                     const Thyra::ModelEvaluatorBase::OutArgs<double>& outArgs) const override
  {
    const RCP<const Thyra::ModelEvaluator<double>> model = this->getUnderlyingModel();
    Thyra::ModelEvaluatorBase::InArgs<double> modelInArgs = model->createInArgs();
    modelInArgs.setArgs(inArgs, true);
    Thyra::ModelEvaluatorBase::OutArgs<double> modelOutArgs = model->createOutArgs();
    modelOutArgs.setArgs(outArgs, true);
    model->evalModel(modelInArgs, modelOutArgs);

    const RCP<Thyra::LinearOpBase<double>> W = outArgs.get_W_op();
    if (Teuchos::nonnull(W)) {
      const RCP<Tpetra_CrsMatrix> W_crs =
          Teuchos::rcp_dynamic_cast<Tpetra_CrsMatrix>(ConverterT::getTpetraOperator(W), true);
      Tpetra::RowMatrixTransposer<double,LO,GO> transposer(W_crs);
      const RCP<Tpetra_CrsMatrix> T = transposer.createTranspose();

      // Copy of the transpose with the column map in reversed order
      const auto colGIDs = T->getColMap()->getLocalElementList();
      Teuchos::Array<GO> reversedColGIDs(colGIDs.size());
      for (int i = 0; i < static_cast<int>(colGIDs.size()); ++i)
        reversedColGIDs[i] = colGIDs[colGIDs.size()-1-i];
      const RCP<const Tpetra_Map> reversedColMap = rcp(new Tpetra_Map(
          Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid(), reversedColGIDs(), T->getColMap()->getIndexBase(),
          T->getColMap()->getComm()));

      const RCP<Tpetra_CrsMatrix> T_reversed = rcp(new Tpetra_CrsMatrix(T->getRowMap(), reversedColMap, T->getGlobalMaxNumRowEntries()));
      typename Tpetra_CrsMatrix::nonconst_global_inds_host_view_type indices("indices", T->getGlobalMaxNumRowEntries());
      typename Tpetra_CrsMatrix::nonconst_values_host_view_type values("values", T->getGlobalMaxNumRowEntries());
      for (std::size_t i = 0; i < T->getLocalNumRows(); ++i) {
        const GO row = T->getRowMap()->getGlobalElement(i);
        std::size_t numEntries = 0;
        T->getGlobalRowCopy(row, indices, values, numEntries);
        T_reversed->insertGlobalValues(row, numEntries, values.data(), indices.data());
      }
      T_reversed->fillComplete(T->getDomainMap(), T->getRangeMap());

      *W_crs = *T_reversed;
    }
  }
};

RCP<Thyra::LinearOpWithSolveFactoryBase<double>> createLinearSolverFactory() {
  Stratimikos::DefaultLinearSolverBuilder builder;
  const RCP<Teuchos::ParameterList> params = Teuchos::parameterList();
  params->set("Linear Solver Type", "Belos");
  params->set("Preconditioner Type", "None");
  Teuchos::ParameterList& belosParams = params->sublist("Linear Solver Types").sublist("Belos");
  belosParams.set("Solver Type", "Block GMRES");
  Teuchos::ParameterList& gmresParams = belosParams.sublist("Solver Types").sublist("Block GMRES");
  gmresParams.set("Convergence Tolerance", 1.0e-12);
  gmresParams.set("Maximum Iterations", 100);
  builder.setParameterList(params);
  return Thyra::createLinearSolveStrategy(builder);
}

} // namespace

// The inverse adjoint Jacobian computed with an adjoint model must match the one computed with the transposed
// Jacobian, also when the adjoint Jacobian is recomputed at different states and when the solvers are reused.
TEUCHOS_UNIT_TEST(Piro_ThyraProductME_Constraint_SimOpt, InverseAdjointJacobianWithAdjointModel)
{
  const auto comm = Tpetra::getDefaultComm();
  const RCP<Thyra::LinearOpWithSolveFactoryBase<double>> lowsFactory = createLinearSolverFactory();

  const RCP<Thyra::ModelEvaluator<double>> model = rcp(new Thyra::DefaultModelEvaluatorWithSolveFactory<double>(
      rcp(new MockModelEval_A_Tpetra(comm)), lowsFactory));
  const RCP<Thyra::ModelEvaluator<double>> adjointModel = rcp(new Thyra::DefaultModelEvaluatorWithSolveFactory<double>(
      rcp(new TransposingAdjointModel(rcp(new MockModelEval_A_Tpetra(comm)))), lowsFactory));

  Teuchos::ParameterList piroParams;
  Piro::ThyraProductME_Constraint_SimOpt<double> constrWithAdjointModel(model, adjointModel, piroParams, Teuchos::VERB_NONE);
  Piro::ThyraProductME_Constraint_SimOpt<double> constrWithTranspose(model, Teuchos::null, piroParams, Teuchos::VERB_NONE);

  ROL::ThyraVector<double> p(model->getNominalValues().get_p(0)->clone_v());
  ROL::ThyraVector<double> v(Thyra::createMember(model->get_f_space()));
  Thyra::seed_randomize<double>(42);
  Thyra::randomize<double>(-1.0, 1.0, v.getVector().ptr());

  for (int i = 0; i < 3; ++i) {
    // State with entries in [1, 3], for which the Jacobian is nonsingular
    ROL::ThyraVector<double> u(Thyra::createMember(model->get_x_space()));
    Thyra::randomize<double>(1.0, 3.0, u.getVector().ptr());
    constrWithAdjointModel.update(u, p, ROL::UpdateType::Temp);
    constrWithTranspose.update(u, p, ROL::UpdateType::Temp);

    // Twice at the same state, the second time reusing the solvers
    for (int j = 0; j < 2; ++j) {
      const auto iajv = u.clone(), iajv_ref = u.clone();
      double tol = 1.0e-12;
      constrWithAdjointModel.applyInverseAdjointJacobian_1(*iajv, v, u, p, tol);
      constrWithTranspose.applyInverseAdjointJacobian_1(*iajv_ref, v, u, p, tol);

      const double ref_norm = iajv_ref->norm();
      iajv->axpy(-1.0, *iajv_ref);
      out << "State " << i << ", solve " << j << ": |iajv - iajv_ref| = " << iajv->norm()
          << ", |iajv_ref| = " << ref_norm << std::endl;
      TEST_COMPARE(iajv->norm(), <=, 1.0e-8 * ref_norm);
    }
  }
}
