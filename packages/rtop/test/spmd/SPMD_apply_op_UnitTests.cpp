// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "RTOpPack_Types.hpp"

#include "Teuchos_UnitTestHarness.hpp"

#include "RTOpPack_SPMD_apply_op.hpp"
#include "RTOpPack_ROpSum.hpp"
#include "RTOpPack_TOpAssignScalar.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_as.hpp"


namespace RTOpPack {


//
// Helpers
//


int g_localDim = 2;
bool g_dumpRTOps = false;


TEUCHOS_STATIC_SETUP()
{
  Teuchos::UnitTestRepository::getCLP().setOption(
    "localDim", &g_localDim, "Number of elements in a process" );
  Teuchos::UnitTestRepository::getCLP().setOption(
    "dump-rtops", "no-dump-rtops", &g_dumpRTOps, "Set if RTOps are dumped or not." );
}


using Teuchos::as;
using Teuchos::inoutArg;

#define PRINT_VAR(varName) \
  out << #varName" = " << (varName) << "\n"


Ordinal computeLocalOffset(const Teuchos::Comm<Ordinal> &comm, const Ordinal localDim)
{
  Ordinal localOffset = 0;
  Teuchos::scan<Ordinal>(comm, Teuchos::REDUCE_SUM, localDim, inoutArg(localOffset));
  localOffset -= localDim;
  return localOffset;
}


template<typename Scalar>
RTOpPack::SubVectorView<Scalar> getLocalSubVectorView(
  const Ordinal localOffset, const Ordinal localDim, const Scalar val)
{
  Teuchos::ArrayRCP<Scalar> x_dat(localDim);
  std::fill_n(x_dat.begin(), localDim, val);
  return RTOpPack::SubVectorView<Scalar>(
    localOffset, localDim, x_dat, 1);
}


template<typename Scalar>
RTOpPack::SubMultiVectorView<Scalar> getLocalSubMultiVectorView(
  const Ordinal localOffset, const Ordinal localDim,
  const Ordinal numCols, const Scalar val)
{
  const Ordinal totalLen = localDim*numCols;
  Teuchos::ArrayRCP<Scalar> x_dat(totalLen);
  std::fill_n(x_dat.begin(), totalLen, val);
  return RTOpPack::SubMultiVectorView<Scalar>(
    localOffset, localDim, 0, numCols, x_dat, localDim);
}


template<typename Scalar>
void assertMultiVectorSums(const Teuchos::Comm<Ordinal> &comm,
  const ConstSubMultiVectorView<Scalar> &mv,
  const Ordinal localDim, const Scalar val, const int numProcsFact,
  Teuchos::FancyOStream &out, bool &success)
{
  const Ordinal numCols = mv.numSubCols();
  RTOpPack::ROpSum<Scalar> sumOp;
  Array<RCP<RTOpPack::ReductTarget> > sumTargets_store(numCols);
  Array<RTOpPack::ReductTarget*> sumTargets(numCols);
  for (int j = 0; j < numCols; ++j) {
    sumTargets_store[j] = sumOp.reduct_obj_create();
    sumTargets[j] = sumTargets_store[j].getRawPtr();
  }
  RTOpPack::SPMD_apply_op<Scalar>(&comm, sumOp, numCols, 1, &mv, 0, 0,
    sumTargets.getRawPtr());

  PRINT_VAR(localDim);
  PRINT_VAR(val);
  PRINT_VAR(numProcsFact);
  for (int j = 0; j < numCols; ++j) {
    PRINT_VAR(j);
    Scalar sum_mv_j = sumOp(*sumTargets[j]);
    TEST_EQUALITY(sum_mv_j, as<Scalar>(localDim*numProcsFact)*val);
  }
}


//
// Unit tests
//


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( SPMD_apply_op, vec_args_1_0_sum, Scalar )
{

  const RCP<const Teuchos::Comm<Ordinal> > comm =
    Teuchos::DefaultComm<Ordinal>::getComm();

  //const int procRank = rank(*comm);
  const Ordinal localDim = g_localDim;
  const Ordinal localOffset = computeLocalOffset(*comm, localDim);
  const Scalar val = 1.1;

  if (g_dumpRTOps) {
    RTOpPack::set_SPMD_apply_op_dump_out(rcpFromRef(out));
  }

  RTOpPack::SubVectorView<Scalar> x =
    getLocalSubVectorView<Scalar>(localOffset, localDim, val);

  RTOpPack::ROpSum<Scalar> sumOp;
  RCP<RTOpPack::ReductTarget> sumTarget = sumOp.reduct_obj_create();
  RTOpPack::SPMD_apply_op<Scalar>(&*comm, sumOp, 1, &x, 0, 0, &*sumTarget);
  Scalar sum_x = sumOp(*sumTarget);

  TEST_EQUALITY(sum_x, as<Scalar>(localDim * comm->getSize()) * val )

  RTOpPack::set_SPMD_apply_op_dump_out(Teuchos::null);

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES(
  SPMD_apply_op, vec_args_1_0_sum)


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( SPMD_apply_op, vec_args_1_0_sum_zero_p0, Scalar )
{

  const RCP<const Teuchos::Comm<Ordinal> > comm =
    Teuchos::DefaultComm<Ordinal>::getComm();

  const int procRank = rank(*comm);
  const Ordinal localDim = (procRank == 0 ? 0 : g_localDim);
  const Ordinal localOffset = computeLocalOffset(*comm, localDim);
  const Scalar val = 1.1;

  if (g_dumpRTOps) {
    RTOpPack::set_SPMD_apply_op_dump_out(rcpFromRef(out));
  }

  RTOpPack::SubVectorView<Scalar> x =
    getLocalSubVectorView<Scalar>(localOffset, localDim, val);

  RTOpPack::ROpSum<Scalar> sumOp;
  RCP<RTOpPack::ReductTarget> sumTarget = sumOp.reduct_obj_create();
  RTOpPack::SPMD_apply_op<Scalar>(&*comm, sumOp, 1, &x, 0, 0, &*sumTarget);
  Scalar sum_x = sumOp(*sumTarget);

  TEST_EQUALITY(sum_x, as<Scalar>(g_localDim * (comm->getSize()-1)) * val);

  RTOpPack::set_SPMD_apply_op_dump_out(Teuchos::null);

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES(
  SPMD_apply_op, vec_args_1_0_sum_zero_p0)


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( SPMD_apply_op, vec_args_1_0_sum_zero_p1, Scalar )
{

  const RCP<const Teuchos::Comm<Ordinal> > comm =
    Teuchos::DefaultComm<Ordinal>::getComm();

  const int procRank = rank(*comm);
  const Ordinal localDim = (procRank == 1 ? 0 : g_localDim);
  PRINT_VAR(g_localDim);
  PRINT_VAR(localDim);
  const Ordinal localOffset = computeLocalOffset(*comm, localDim);
  const Scalar val = 1.1;
  PRINT_VAR(val);

  if (g_dumpRTOps) {
    RTOpPack::set_SPMD_apply_op_dump_out(rcpFromRef(out));
  }

  RTOpPack::SubVectorView<Scalar> x =
    getLocalSubVectorView<Scalar>(localOffset, localDim, val);

  RTOpPack::ROpSum<Scalar> sumOp;
  RCP<RTOpPack::ReductTarget> sumTarget = sumOp.reduct_obj_create();
  RTOpPack::SPMD_apply_op<Scalar>(&*comm, sumOp, 1, &x, 0, 0, &*sumTarget);
  Scalar sum_x = sumOp(*sumTarget);

  const Ordinal procFactor = (comm->getSize() == 1 ? 1 : (comm->getSize()-1));
  PRINT_VAR(procFactor);
  TEST_EQUALITY(sum_x, as<Scalar>(g_localDim * procFactor) * val);

  RTOpPack::set_SPMD_apply_op_dump_out(Teuchos::null);

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES(
  SPMD_apply_op, vec_args_1_0_sum_zero_p1)


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( SPMD_apply_op, multivec_args_1_0_sum, Scalar )
{

  const RCP<const Teuchos::Comm<Ordinal> > comm =
    Teuchos::DefaultComm<Ordinal>::getComm();

  //const int procRank = rank(*comm);

  const Ordinal localDim = g_localDim;
  PRINT_VAR(localDim);
  const Ordinal numCols = 3;
  PRINT_VAR(numCols);
  const Ordinal localOffset = computeLocalOffset(*comm, localDim);
  PRINT_VAR(localOffset);
  const Scalar val = 1.1;
  PRINT_VAR(val);

  RTOpPack::SubMultiVectorView<Scalar> mv =
    getLocalSubMultiVectorView<Scalar>(localOffset, localDim, numCols, val);

  assertMultiVectorSums<Scalar>(*comm, mv, g_localDim, val, comm->getSize(),
    out, success);

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES(
  SPMD_apply_op, multivec_args_1_0_sum)


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( SPMD_apply_op, multivec_args_1_0_sum_zero_p0, Scalar )
{

  const RCP<const Teuchos::Comm<Ordinal> > comm =
    Teuchos::DefaultComm<Ordinal>::getComm();

  const int procRank = rank(*comm);

  const Ordinal localDim = (procRank == 0 ? 0 : g_localDim);
  PRINT_VAR(localDim);
  const Ordinal numCols = 3;
  PRINT_VAR(numCols);
  const Ordinal localOffset = computeLocalOffset(*comm, localDim);
  PRINT_VAR(localOffset);
  const Scalar val = 1.1;
  PRINT_VAR(val);

  if (g_dumpRTOps) {
    RTOpPack::set_SPMD_apply_op_dump_out(rcpFromRef(out));
  }

  RTOpPack::SubMultiVectorView<Scalar> mv =
    getLocalSubMultiVectorView<Scalar>(localOffset, localDim, numCols, val);

  assertMultiVectorSums<Scalar>(*comm, mv, g_localDim, val, comm->getSize()-1,
    out, success);

  RTOpPack::set_SPMD_apply_op_dump_out(Teuchos::null);

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES(
  SPMD_apply_op, multivec_args_1_0_sum_zero_p0)


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( SPMD_apply_op, multivec_args_0_1_assign, Scalar )
{

  const RCP<const Teuchos::Comm<Ordinal> > comm =
    Teuchos::DefaultComm<Ordinal>::getComm();

  //const int procRank = rank(*comm);

  const Ordinal localDim = g_localDim;
  PRINT_VAR(localDim);
  const Ordinal numCols = 3;
  PRINT_VAR(numCols);
  const Ordinal localOffset = computeLocalOffset(*comm, localDim);
  PRINT_VAR(localOffset);
  const Scalar val_init = -0.1;
  PRINT_VAR(val_init);

  if (g_dumpRTOps) {
    RTOpPack::set_SPMD_apply_op_dump_out(rcpFromRef(out));
  }

  RTOpPack::SubMultiVectorView<Scalar> mv =
    getLocalSubMultiVectorView<Scalar>(localOffset, localDim, numCols, val_init);

  const Scalar val = 1.2;
  PRINT_VAR(val);
  RTOpPack::TOpAssignScalar<Scalar> assignOp(val);
  RTOpPack::SPMD_apply_op<Scalar>(&*comm, assignOp, numCols, 0, 0, 1, &mv, 0);

  assertMultiVectorSums<Scalar>(*comm, mv, g_localDim, val, comm->getSize(),
    out, success);

  RTOpPack::set_SPMD_apply_op_dump_out(Teuchos::null);

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES(
  SPMD_apply_op, multivec_args_0_1_assign)


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( SPMD_apply_op, multivec_args_0_1_assign_zero_p0, Scalar )
{

  const RCP<const Teuchos::Comm<Ordinal> > comm =
    Teuchos::DefaultComm<Ordinal>::getComm();

  if (g_dumpRTOps) {
    RTOpPack::set_SPMD_apply_op_dump_out(rcpFromRef(out));
  }

  const int procRank = rank(*comm);

  const Ordinal localDim = (procRank == 0 ? 0 : g_localDim);
  PRINT_VAR(localDim);
  const Ordinal numCols = 3;
  PRINT_VAR(numCols);
  const Ordinal localOffset = computeLocalOffset(*comm, localDim);
  PRINT_VAR(localOffset);
  const Scalar val_init = -0.1;
  PRINT_VAR(val_init);

  RTOpPack::SubMultiVectorView<Scalar> mv =
    getLocalSubMultiVectorView<Scalar>(localOffset, localDim, numCols, val_init);

  const Scalar val = 1.2;
  PRINT_VAR(val);
  RTOpPack::TOpAssignScalar<Scalar> assignOp(val);
  RTOpPack::SPMD_apply_op<Scalar>(&*comm, assignOp, numCols, 0, 0, 1, &mv, 0);

  assertMultiVectorSums<Scalar>(*comm, mv, g_localDim, val, comm->getSize()-1,
    out, success);

  RTOpPack::set_SPMD_apply_op_dump_out(Teuchos::null);

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES(
  SPMD_apply_op, multivec_args_0_1_assign_zero_p0)




} // namespace RTOpPack
