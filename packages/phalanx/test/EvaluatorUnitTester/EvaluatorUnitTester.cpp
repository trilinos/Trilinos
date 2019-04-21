// @HEADER
// ************************************************************************
//
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                    Copyright 2008 Sandia Corporation
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER

#include "Phalanx_KokkosDeviceTypes.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_FieldTag_Tag.hpp"
#include "Phalanx_FieldManager.hpp"
#include "Phalanx_TypeStrings.hpp"
#include "Phalanx_DimTag.hpp"
#include "Phalanx_Evaluator_UnmanagedFieldDummy.hpp"
#include "Phalanx_Evaluator_UnitTester.hpp"
#include "Phalanx_MDField_UnmanagedAllocator.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include <limits>

#include "MyTraits.hpp"

PHX_DIM_TAG_DECLARATION(CELL)
PHX_DIM_TAG_IMPLEMENTATION(CELL)

PHX_DIM_TAG_DECLARATION(QP)
PHX_DIM_TAG_IMPLEMENTATION(QP)

PHX_DIM_TAG_DECLARATION(DIM)
PHX_DIM_TAG_IMPLEMENTATION(DIM)

PHX_DIM_TAG_DECLARATION(R)
PHX_DIM_TAG_IMPLEMENTATION(R)

// requires the dim tags defined above
#include "SimpleEvaluator.hpp"
#include "DuplicateFieldEvaluator.hpp"
#include "AllRanksEvaluator.hpp"

TEUCHOS_UNIT_TEST(evaluator_unit_tester, simple)
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;
  using EvalType = MyTraits::Residual;
  using Scalar = EvalType::ScalarT;

  const int num_cells = 10;
  const int num_qp = 8;
  const int num_dim = 3;
  RCP<MDALayout<CELL,QP>> adl = rcp(new MDALayout<CELL,QP>(num_cells,num_qp));
  RCP<MDALayout<CELL,QP>> bdl = adl;
  RCP<MDALayout<CELL,QP,DIM>> cdl = rcp(new MDALayout<CELL,QP,DIM>(num_cells,num_qp,num_dim));

  RCP<SimpleEvaluator<EvalType,MyTraits>> e = rcp(new SimpleEvaluator<EvalType,MyTraits>(adl,bdl,cdl));

  MDField<Scalar,CELL,QP> b = allocateUnmanagedMDField<Scalar,CELL,QP>("b",bdl);
  b.deep_copy(2.0);
  MDField<Scalar,CELL,QP,DIM> c = allocateUnmanagedMDField<Scalar,CELL,QP,DIM>("c",cdl);
  c.deep_copy(3.0);
  Kokkos::fence();

  PHX::EvaluatorUnitTester<EvalType,MyTraits> tester;
  tester.setEvaluatorToTest(e);
  tester.setDependentFieldValues(b);
  tester.setDependentFieldValues(c);
  tester.testEvaluator(num_cells,num_cells,num_cells,num_cells);
  Kokkos::fence();

  MDField<Scalar,CELL,QP> pyrite_a = allocateUnmanagedMDField<Scalar,CELL,QP>("a",adl);
  pyrite_a.deep_copy(36.0);
  const Scalar tol = 1000.0 * std::numeric_limits<Scalar>::epsilon();
  tester.checkFloatValues2(pyrite_a,tol,success,out);
}

TEUCHOS_UNIT_TEST(evaluator_unit_tester, DuplicateField)
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;
  using EvalType = MyTraits::Residual;
  using Scalar = EvalType::ScalarT;

  const int num_cells = 10;
  const int num_qp = 8;
  const int num_dim = 3;
  RCP<MDALayout<CELL,QP>> adl = rcp(new MDALayout<CELL,QP>(num_cells,num_qp));
  RCP<MDALayout<CELL,QP>> bdl = adl;
  RCP<MDALayout<CELL,QP,DIM>> cdl = rcp(new MDALayout<CELL,QP,DIM>(num_cells,num_qp,num_dim));

  RCP<DuplicateFieldEvaluator<EvalType,MyTraits>> e = 
    rcp(new DuplicateFieldEvaluator<EvalType,MyTraits>(adl,bdl,cdl));

  MDField<Scalar,CELL,QP> b = allocateUnmanagedMDField<Scalar,CELL,QP>("b",bdl);
  b.deep_copy(2.0);
  MDField<Scalar,CELL,QP,DIM> c = allocateUnmanagedMDField<Scalar,CELL,QP,DIM>("c",cdl);
  c.deep_copy(3.0);
  Kokkos::fence();

  PHX::EvaluatorUnitTester<EvalType,MyTraits> tester;
  tester.setEvaluatorToTest(e);
  tester.setDependentFieldValues(b);
  tester.setDependentFieldValues(c);
  tester.testEvaluator(num_cells,num_cells,num_cells,num_cells);
  Kokkos::fence();

  MDField<Scalar,CELL,QP> pyrite_a = allocateUnmanagedMDField<Scalar,CELL,QP>("a",adl);
  pyrite_a.deep_copy(36.0);
  const Scalar tol = 1000.0 * std::numeric_limits<Scalar>::epsilon();
  tester.checkFloatValues2(pyrite_a,tol,success,out);
}

TEUCHOS_UNIT_TEST(evaluator_unit_tester, AllRanks)
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;
  using EvalType = MyTraits::Residual;
  using Scalar = EvalType::ScalarT;

  const int r1 = 2;
  const int r2 = 2;
  const int r3 = 2;
  const int r4 = 2;
  const int r5 = 2;
  const int r6 = 2;
  RCP<MDALayout<R>> dl1 = rcp(new MDALayout<R>(r1));
  RCP<MDALayout<R,R>> dl2 = rcp(new MDALayout<R,R>(r1,r2));
  RCP<MDALayout<R,R,R>> dl3 = rcp(new MDALayout<R,R,R>(r1,r2,r3));
  RCP<MDALayout<R,R,R,R>> dl4 = rcp(new MDALayout<R,R,R,R>(r1,r2,r3,r4));
  RCP<MDALayout<R,R,R,R,R>> dl5 = rcp(new MDALayout<R,R,R,R,R>(r1,r2,r3,r4,r5));
  RCP<MDALayout<R,R,R,R,R,R>> dl6 = rcp(new MDALayout<R,R,R,R,R,R>(r1,r2,r3,r4,r5,r6));
  
  RCP<AllRanksEvaluator<EvalType,MyTraits>> e = 
    rcp(new AllRanksEvaluator<EvalType,MyTraits>(dl1,dl2,dl3,dl4,dl5,dl6));

  MDField<Scalar,R> f1 = allocateUnmanagedMDField<Scalar,R>("f1",dl1);
  MDField<Scalar,R,R> f2 = allocateUnmanagedMDField<Scalar,R,R>("f2",dl2);
  MDField<Scalar,R,R,R> f3 = allocateUnmanagedMDField<Scalar,R,R,R>("f3",dl3);
  MDField<Scalar,R,R,R,R> f4 = allocateUnmanagedMDField<Scalar,R,R,R,R>("f4",dl4);
  MDField<Scalar,R,R,R,R,R> f5 = allocateUnmanagedMDField<Scalar,R,R,R,R,R>("f5",dl5);
  MDField<Scalar,R,R,R,R,R,R> f6 = allocateUnmanagedMDField<Scalar,R,R,R,R,R,R>("f6",dl6);
  f1.deep_copy(1.0);
  f2.deep_copy(2.0);
  f3.deep_copy(3.0);
  f4.deep_copy(4.0);
  f5.deep_copy(5.0);
  f6.deep_copy(6.0);
  Kokkos::fence();

  PHX::EvaluatorUnitTester<EvalType,MyTraits> tester;
  tester.setEvaluatorToTest(e);
  tester.setDependentFieldValues(f1);
  tester.setDependentFieldValues(f2);
  tester.setDependentFieldValues(f3);
  tester.setDependentFieldValues(f4);
  tester.setDependentFieldValues(f5);
  tester.setDependentFieldValues(f6);
  tester.testEvaluator(r1,r1,r1,r1);
  Kokkos::fence();

  MDField<Scalar,R> gx1 = allocateUnmanagedMDField<Scalar,R>("x1",dl1);
  MDField<Scalar,R,R> gx2 = allocateUnmanagedMDField<Scalar,R,R>("x2",dl2);
  MDField<Scalar,R,R,R> gx3 = allocateUnmanagedMDField<Scalar,R,R,R>("x3",dl3);
  MDField<Scalar,R,R,R,R> gx4 = allocateUnmanagedMDField<Scalar,R,R,R,R>("x4",dl4);
  MDField<Scalar,R,R,R,R,R> gx5 = allocateUnmanagedMDField<Scalar,R,R,R,R,R>("x5",dl5);
  MDField<Scalar,R,R,R,R,R,R> gx6 = allocateUnmanagedMDField<Scalar,R,R,R,R,R,R>("x6",dl6);
  gx1.deep_copy(1.0);
  gx2.deep_copy(4.0);
  gx3.deep_copy(9.0);
  gx4.deep_copy(16.0);
  gx5.deep_copy(25.0);
  gx6.deep_copy(36.0);
  const Scalar tol = 1000.0 * std::numeric_limits<Scalar>::epsilon();
  tester.checkFloatValues1(gx1,tol,success,out);
  tester.checkFloatValues2(gx2,tol,success,out);
  tester.checkFloatValues3(gx3,tol,success,out);
  tester.checkFloatValues4(gx4,tol,success,out);
  tester.checkFloatValues5(gx5,tol,success,out);
  tester.checkFloatValues6(gx6,tol,success,out);
}
