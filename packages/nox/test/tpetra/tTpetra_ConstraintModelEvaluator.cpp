//@HEADER
// ************************************************************************
//
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER
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

#include "Thyra_TpetraThyraWrappers.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "ME_Tpetra_1DFEM.hpp"

#include "LOCA_Tpetra_ConstraintModelEvaluator.hpp"
#include "NOX_Thyra_Vector.H"
#include "NOX_Thyra_MultiVector.H"
#include "NOX_TpetraTypedefs.hpp"

TEUCHOS_UNIT_TEST(NOX_Tpetra_1DFEM, Responses)
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

  LOCA::ParameterVector p_vec;
  p_vec.addParameter("k",1.0);
  p_vec.addParameter("T_left",1.0);
  std::vector<std::string> g_names;
  g_names.push_back("Constraint: T_right=2");
  g_names.push_back("Constraint: 2*T_left=T_right");
  auto x_thyra = ::Thyra::createMember(model->get_x_space(),"x");
  NOX::Thyra::Vector x(x_thyra);
  LOCA::MultiContinuation::ConstraintModelEvaluator constraints(model,p_vec,g_names,x);

  // Set inputs
  x.init(1.0);
  constraints.setX(x);
  constraints.setParam(0,1.0);
  constraints.setParam(1,1.2);

  // Compute constraints
  out << "Checking constraints:\n";
  constraints.computeConstraints();
  auto g = constraints.getConstraints();
  auto tol = std::numeric_limits<Scalar>::epsilon()*100.0;
  out << "g(0) = " << g(0,0) << std::endl;
  TEST_FLOATING_EQUALITY(g(0,0),-1.0,tol);
  out << "g(1) = " << g(1,0) << std::endl;
  TEST_FLOATING_EQUALITY(g(1,0),1.0,tol);

  // Compute DX
  out << "\nChecking DgDx:\n";
  constraints.computeDX();
  auto input = x.createMultiVector(1,NOX::ShapeCopy);
  input->init(1.0);
  auto result = x.createMultiVector(2,NOX::ShapeCopy);
  NOX::Abstract::MultiVector::DenseMatrix b(2,2);
  b(0,0) = 1.0;
  b(0,1) = 0.0;
  b(1,0) = 0.0;
  b(1,1) = 1.0;
  constraints.addDX(Teuchos::NO_TRANS,1.0,b,0.0,*result);
  auto result_thyra = Teuchos::rcp_dynamic_cast<NOX::Thyra::MultiVector>(result)->getThyraMultiVector();
  auto DgDx = converter::getTpetraMultiVector(result_thyra);
  DgDx->sync_host();
  DgDx->describe(out,Teuchos::VERB_EXTREME);
  auto DgDx_host = DgDx->getLocalViewHost();
  for (size_t i=0; i < DgDx_host.extent(0); ++i) {
    if ( (comm->getRank() == 0) && (i == 0) ) {
      TEST_FLOATING_EQUALITY(DgDx_host(0,0),0.0,tol);
      TEST_FLOATING_EQUALITY(DgDx_host(0,1),2.0,tol);
    }
    else if ( (comm->getRank() == (comm->getSize()-1)) && (i == DgDx_host.extent(0)-1) ) {
      TEST_FLOATING_EQUALITY(DgDx_host(DgDx_host.extent(0)-1,0),1.0,tol);
      TEST_FLOATING_EQUALITY(DgDx_host(DgDx_host.extent(0)-1,1),-1.0,tol);
    }
    else {
      TEST_FLOATING_EQUALITY(DgDx_host(i,0),0.0,tol);
      TEST_FLOATING_EQUALITY(DgDx_host(i,1),0.0,tol);
    }
  }
  Teuchos::Array<NOX::TMultiVector::mag_type> norms(1);
  DgDx->norm2(norms);
  TEST_FLOATING_EQUALITY(norms[0],1.0,tol);

  // Compute g and DgDp
  out << "\nChecking DgDp:\n";
  NOX::Abstract::MultiVector::DenseMatrix dgdp(2,3,true); // first col is g
  std::vector<int> paramIDs(2);
  paramIDs[0] = 0;
  paramIDs[1] = 1;
  constraints.computeDP(paramIDs,dgdp,false);
  TEST_FLOATING_EQUALITY(dgdp(0,0),-1.0,tol);
  TEST_FLOATING_EQUALITY(dgdp(1,0),1.0,tol);
  TEST_FLOATING_EQUALITY(dgdp(0,1),0.0,tol);
  TEST_FLOATING_EQUALITY(dgdp(0,2),0.0,tol);
  TEST_FLOATING_EQUALITY(dgdp(1,1),0.0,tol);
  TEST_FLOATING_EQUALITY(dgdp(1,2),0.0,tol);

  // Test clone/copy
  out << "\nChecking clone and copy methods:\n";
  TEST_ASSERT(constraints.isConstraints());
  TEST_ASSERT(constraints.isDX());

  // Shape copy does not copy values, not equal
  auto constraints_clone = constraints.clone(NOX::ShapeCopy);
  TEST_ASSERT(constraints_clone->isConstraints() != constraints.isConstraints());
  TEST_ASSERT(constraints_clone->isDX() != constraints.isDX());

  // Copy method should now make them equal
  constraints_clone->copy(constraints);
  TEST_ASSERT(constraints_clone->isConstraints() == constraints.isConstraints());
  TEST_ASSERT(constraints_clone->isDX() == constraints.isDX());

  // Let's change the solution so everyhting is invalid in the clone
  constraints_clone->setParam(0,0.0);
  TEST_ASSERT(constraints_clone->isConstraints() != constraints.isConstraints());
  TEST_ASSERT(constraints_clone->isDX() != constraints.isDX());

  // Deep copy should make everything equal
  constraints_clone = constraints.clone(NOX::DeepCopy);
  TEST_ASSERT(constraints_clone->isConstraints() == constraints.isConstraints());
  TEST_ASSERT(constraints_clone->isDX() == constraints.isDX());
}
