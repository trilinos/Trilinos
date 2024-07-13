// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_StackedTimer.hpp"

// NOX Objects
#include "NOX.H"
#include "NOX_Thyra.H"

// Trilinos Objects
#include "Teuchos_Comm.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_AbstractFactoryStd.hpp"

#include "Tpetra_Core.hpp"
#include "Tpetra_Vector.hpp"

#include "BelosTypes.hpp"
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "ME_Tpetra_1DFEM.hpp"

#include "NOX_Thyra_MatrixFreeJacobianOperator.hpp"
#include "NOX_MatrixFree_ModelEvaluatorDecorator.hpp"
#include "NOX_TpetraTypedefs.hpp"

TEUCHOS_UNIT_TEST(NOX_Tpetra_1DFEM, Responses_g4_p2)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  // Get default Tpetra template types
  using Scalar = NOX::Scalar;
  using LO = NOX::LocalOrdinal;
  using GO = NOX::GlobalOrdinal;
  using Node = NOX::NodeType;
  using converter = ::Thyra::TpetraOperatorVectorExtraction<Scalar,LO,GO,Node>;

  // Create the model evaluator object
  Scalar x00 = 0.0;
  Scalar x01 = 1.0;
  const Tpetra::global_size_t numGlobalElements = 100;
  Teuchos::RCP<EvaluatorTpetra1DFEM<Scalar,LO,GO,Node> > model =
    evaluatorTpetra1DFEM<Scalar,LO,GO,Node>(comm, numGlobalElements, x00, x01);

  auto g_thyra = ::Thyra::createMember(*model->get_g_space(4),"g");
  auto DfDp_thyra = Teuchos::rcp_dynamic_cast<::Thyra::MultiVectorBase<Scalar>>(model->create_DfDp_op(2),true);
  auto DgDx_thyra = Teuchos::rcp_dynamic_cast<::Thyra::MultiVectorBase<Scalar>>(model->create_DgDx_op(4),true);
  auto Dg4Dp2_thyra = Teuchos::rcp_dynamic_cast<::Thyra::MultiVectorBase<Scalar>>(model->create_DgDp_op(4,2),true);
  auto Dg4Dp4_thyra = Teuchos::rcp_dynamic_cast<::Thyra::MultiVectorBase<Scalar>>(model->create_DgDp_op(4,4),true);

  auto inArgs = model->createInArgs();
  auto x = ::Thyra::createMember(model->get_x_space(),"x");
  ::Thyra::assign(x.ptr(),1.0);
  inArgs.set_x(x);
  auto outArgs = model->createOutArgs();

  outArgs.set_g(4,::Thyra::ModelEvaluatorBase::Evaluation<::Thyra::VectorBase<Scalar>>(g_thyra));
  outArgs.set_DfDp(2,::Thyra::ModelEvaluatorBase::Derivative<Scalar>(DfDp_thyra));
  outArgs.set_DgDx(4,::Thyra::ModelEvaluatorBase::Derivative<Scalar>(DgDx_thyra,Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM));
  outArgs.set_DgDp(4,2,::Thyra::ModelEvaluatorBase::Derivative<Scalar>(Dg4Dp2_thyra));
  outArgs.set_DgDp(4,4,::Thyra::ModelEvaluatorBase::Derivative<Scalar>(Dg4Dp4_thyra));

  model->evalModel(inArgs,outArgs);

  TEST_EQUALITY(model->Np(),5);
  TEST_EQUALITY(model->Ng(),7);

  auto g = converter::getTpetraMultiVector(g_thyra);
  auto DfDp = converter::getTpetraMultiVector(DfDp_thyra);
  auto DgDx = converter::getTpetraMultiVector(DgDx_thyra);
  auto Dg4Dp2 = converter::getTpetraMultiVector(Dg4Dp2_thyra);
  auto Dg4Dp4 = converter::getTpetraMultiVector(Dg4Dp4_thyra);

  // g = T(Zmax) - 2.0
  // x=1.0
  auto g_host = g->getLocalViewHost(Tpetra::Access::ReadOnly);
  out << "g4 = " << g_host(0,0) << std::endl;
  auto tol = std::numeric_limits<Scalar>::epsilon()*100.0;
  TEST_FLOATING_EQUALITY(g_host(0,0),-1.0,tol);

  // DgDx: Right end node is 1, the rest of the vector is zero.
  out << "Dg4Dx:\n";
  DgDx->describe(out,Teuchos::VERB_EXTREME);
  {
    auto DgDx_host = DgDx->getLocalViewHost(Tpetra::Access::ReadOnly);
    if (comm->getRank() == (comm->getSize() - 1)) {
      TEST_FLOATING_EQUALITY(DgDx_host(DgDx_host.extent(0) - 1, 0), 1.0, tol);
    }
  }

  Teuchos::Array<NOX::TMultiVector::mag_type> norms(1);
  DgDx->norm2(norms);
  TEST_FLOATING_EQUALITY(norms[0],1.0,tol);

  // DfDp
  out << "DfDp2:\n";
  DfDp->describe(out,Teuchos::VERB_EXTREME);
  auto DfDp_host = DfDp->getLocalViewHost(Tpetra::Access::ReadOnly);
  auto rank = comm->getRank();
  auto size = comm->getSize();
  for (size_t i=0; i < DfDp_host.extent(0); ++i) {
    if ((rank == 0) && (i == 0)) { // left end
      TEST_FLOATING_EQUALITY(DfDp_host(0,0),0.0,tol);
    }
    else if (rank == (size - 1) && (i==(DfDp_host.extent(0)-1))) { // right end
      TEST_FLOATING_EQUALITY(DfDp_host(DfDp_host.extent(0)-1,0),0.005,tol);
    }
    else { // interior nodes
      out << "rank=" << rank << ", i=" << i <<std::endl;
      TEST_FLOATING_EQUALITY(DfDp_host(i,0),0.01,tol);
    }
  }

  out << "Dg4Dp2:\n";
  Dg4Dp2->describe(out,Teuchos::VERB_EXTREME);
  auto Dg4Dp2_host = Dg4Dp2->getLocalViewHost(Tpetra::Access::ReadOnly);
  TEST_FLOATING_EQUALITY(Dg4Dp2_host(0,0),0.0,tol);

  out << "Dg4Dp4:\n";
  Dg4Dp4->describe(out,Teuchos::VERB_EXTREME);
  auto Dg4Dp4_host = Dg4Dp4->getLocalViewHost(Tpetra::Access::ReadOnly);
  TEST_FLOATING_EQUALITY(Dg4Dp4_host(0,0),0.0,tol);
}

TEUCHOS_UNIT_TEST(NOX_Tpetra_1DFEM, Responses_g6_p4)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  // Get default Tpetra template types
  using Scalar = NOX::Scalar;
  using LO = NOX::LocalOrdinal;
  using GO = NOX::GlobalOrdinal;
  using Node = NOX::NodeType;
  using converter = ::Thyra::TpetraOperatorVectorExtraction<Scalar,LO,GO,Node>;

  // Create the model evaluator object
  Scalar x00 = 0.0;
  Scalar x01 = 1.0;
  const Tpetra::global_size_t numGlobalElements = 100;
  Teuchos::RCP<EvaluatorTpetra1DFEM<Scalar,LO,GO,Node> > model =
    evaluatorTpetra1DFEM<Scalar,LO,GO,Node>(comm, numGlobalElements, x00, x01);

  auto g_thyra = ::Thyra::createMember(*model->get_g_space(6),"g");
  auto DfDp4_thyra = Teuchos::rcp_dynamic_cast<::Thyra::MultiVectorBase<Scalar>>(model->create_DfDp_op(4),true);
  auto DgDx_thyra = Teuchos::rcp_dynamic_cast<::Thyra::MultiVectorBase<Scalar>>(model->create_DgDx_op(6),true);
  auto Dg6Dp2_thyra = Teuchos::rcp_dynamic_cast<::Thyra::MultiVectorBase<Scalar>>(model->create_DgDp_op(6,2),true);
  auto Dg6Dp4_thyra = Teuchos::rcp_dynamic_cast<::Thyra::MultiVectorBase<Scalar>>(model->create_DgDp_op(6,4),true);

  auto inArgs = model->createInArgs();
  auto x = ::Thyra::createMember(model->get_x_space(),"x");
  ::Thyra::assign(x.ptr(),1.0);
  inArgs.set_x(x);
  auto outArgs = model->createOutArgs();

  outArgs.set_g(6,::Thyra::ModelEvaluatorBase::Evaluation<::Thyra::VectorBase<Scalar>>(g_thyra));
  outArgs.set_DfDp(4,::Thyra::ModelEvaluatorBase::Derivative<Scalar>(DfDp4_thyra));
  outArgs.set_DgDx(6,::Thyra::ModelEvaluatorBase::Derivative<Scalar>(DgDx_thyra,Thyra::ModelEvaluatorBase::DERIV_MV_GRADIENT_FORM));
  outArgs.set_DgDp(6,2,::Thyra::ModelEvaluatorBase::Derivative<Scalar>(Dg6Dp2_thyra));
  outArgs.set_DgDp(6,4,::Thyra::ModelEvaluatorBase::Derivative<Scalar>(Dg6Dp4_thyra));

  model->evalModel(inArgs,outArgs);

  TEST_EQUALITY(model->Np(),5);
  TEST_EQUALITY(model->Ng(),7);

  auto g = converter::getTpetraMultiVector(g_thyra);
  auto DfDp4 = converter::getTpetraMultiVector(DfDp4_thyra);
  auto DgDx = converter::getTpetraMultiVector(DgDx_thyra);
  auto Dg6Dp2 = converter::getTpetraMultiVector(Dg6Dp2_thyra);
  auto Dg6Dp4 = converter::getTpetraMultiVector(Dg6Dp4_thyra);

  // g = 2.0 * T(Zmax) - T(zMin)
  // x=1.0
  auto g_host = g->getLocalViewHost(Tpetra::Access::ReadOnly);
  out << "g6 = " << g_host(0,0) << std::endl;
  auto tol = std::numeric_limits<Scalar>::epsilon()*100.0;
  TEST_FLOATING_EQUALITY(g_host(0,0),1.0,tol);

  // DgDx: Left end node is 2, right end node is -1, the rest of the vector is zero.
  out << "Dg6Dx:\n";
  DgDx->describe(out,Teuchos::VERB_EXTREME);

  {
    auto DgDx_host = DgDx->getLocalViewHost(Tpetra::Access::ReadOnly);
    if (comm->getRank() == 0) {
      TEST_FLOATING_EQUALITY(DgDx_host(0, 0), 2.0, tol);
    }
    else if (comm->getRank() == (comm->getSize() - 1)) {
      TEST_FLOATING_EQUALITY(DgDx_host(DgDx_host.extent(0) - 1, 0), -1.0, tol);
    }
  }

  Teuchos::Array<NOX::TMultiVector::mag_type> norms(1);
  DgDx->norm2(norms);
  TEST_FLOATING_EQUALITY(norms[0],std::sqrt(5.0),tol);

  // DfDp4
  out << "DfDp4:\n";
  DfDp4->describe(out,Teuchos::VERB_EXTREME);
  auto DfDp4_host = DfDp4->getLocalViewHost(Tpetra::Access::ReadOnly);
  auto rank = comm->getRank();
  for (size_t i=0; i < DfDp4_host.extent(0); ++i) {
    if ((rank == 0) && (i == 0)) { // left end
      TEST_FLOATING_EQUALITY(DfDp4_host(0,0),-1.0,tol);
    }
    else { // all other nodes are zero
      out << "rank=" << rank << ", i=" << i <<std::endl;
      TEST_FLOATING_EQUALITY(DfDp4_host(i,0),0.0,tol);
    }
  }

  out << "Dg6Dp2:\n";
  Dg6Dp2->describe(out,Teuchos::VERB_EXTREME);
  auto Dg6Dp2_host = Dg6Dp2->getLocalViewHost(Tpetra::Access::ReadOnly);
  TEST_FLOATING_EQUALITY(Dg6Dp2_host(0,0),0.0,tol);

  out << "Dg6Dp4:\n";
  Dg6Dp4->describe(out,Teuchos::VERB_EXTREME);
  auto Dg6Dp4_host = Dg6Dp4->getLocalViewHost(Tpetra::Access::ReadOnly);
  TEST_FLOATING_EQUALITY(Dg6Dp4_host(0,0),0.0,tol);
}

TEUCHOS_UNIT_TEST(NOX_Tpetra_1DFEM,get_index)
{
  // Get default Tpetra template types
  using Scalar = NOX::Scalar;
  using LO = NOX::LocalOrdinal;
  using GO = NOX::GlobalOrdinal;
  using Node = NOX::NodeType;

  // Create the model evaluator object
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = 100;
  Scalar x00 = 0.0;
  Scalar x01 = 1.0;
  Teuchos::RCP<EvaluatorTpetra1DFEM<Scalar,LO,GO,Node> > model =
    evaluatorTpetra1DFEM<Scalar,LO,GO,Node>(comm, numGlobalElements, x00, x01);

  TEST_EQUALITY(model->get_p_index("Dummy p(0)"),std::make_pair(0,0));
  TEST_EQUALITY(model->get_p_index("Dummy p(1)"),std::make_pair(1,0));
  TEST_EQUALITY(model->get_p_index("k"),std::make_pair(2,0));
  TEST_EQUALITY(model->get_p_index("Dummy p(3)"),std::make_pair(3,0));
  TEST_EQUALITY(model->get_p_index("T_left"),std::make_pair(4,0));

  TEST_EQUALITY(model->get_g_index("Dummy g(0)"),std::make_pair(0,0));
  TEST_EQUALITY(model->get_g_index("Dummy g(1)"),std::make_pair(1,0));
  TEST_EQUALITY(model->get_g_index("Dummy g(2)"),std::make_pair(2,0));
  TEST_EQUALITY(model->get_g_index("Dummy g(3)"),std::make_pair(3,0));
  TEST_EQUALITY(model->get_g_index("Constraint: T_right=2"),std::make_pair(4,0));
  TEST_EQUALITY(model->get_g_index("Dummy g(5)"),std::make_pair(5,0));
  TEST_EQUALITY(model->get_g_index("Constraint: 2*T_left=T_right"),std::make_pair(6,0));
}
