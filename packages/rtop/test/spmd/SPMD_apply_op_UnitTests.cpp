/*
// @HEADER
// ***********************************************************************
// 
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//                Copyright (2006) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/


#include "RTOpPack_Types.hpp"

#include "Teuchos_UnitTestHarness.hpp"

#include "RTOpPack_SPMD_apply_op.hpp"
#include "RTOpPack_ROpSum.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_as.hpp"


namespace RTOpPack {


//
// Helpers
//


int g_localDim = 2;


TEUCHOS_STATIC_SETUP()
{
  Teuchos::UnitTestRepository::getCLP().setOption(
    "localDim", &g_localDim, "Number of elements in a process" );
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
  const Ordinal numCols,
  const Scalar val)
{
  const Ordinal totalLen = localDim*numCols;
  Teuchos::ArrayRCP<Scalar> x_dat(totalLen);
  std::fill_n(x_dat.begin(), totalLen, val);
  return RTOpPack::SubMultiVectorView<Scalar>(
    localOffset, localDim, 0, numCols, x_dat, localDim);
}


//
// Unit tests
//


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( SPMD_apply_op, vec_args_1_0_reduce, Scalar )
{

  const RCP<const Teuchos::Comm<Ordinal> > comm =
    Teuchos::DefaultComm<Ordinal>::getComm();

  //const int procRank = rank(*comm);

  const Ordinal localDim = g_localDim;

  const Ordinal localOffset = computeLocalOffset(*comm, localDim);
  
  const Scalar val = 1.1;

  RTOpPack::SubVectorView<Scalar> x =
    getLocalSubVectorView<Scalar>(localOffset, localDim, val);

  RTOpPack::ROpSum<Scalar> sumOp;
  RCP<RTOpPack::ReductTarget> sumTarget = sumOp.reduct_obj_create();
  RTOpPack::SPMD_apply_op<Scalar>(&*comm, sumOp, 1, &x, 0, 0, &*sumTarget);
  Scalar sum_x = sumOp(*sumTarget);

  TEST_EQUALITY(sum_x, as<Scalar>(localDim * val * comm->getSize()))

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES(
  SPMD_apply_op, vec_args_1_0_reduce)


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( SPMD_apply_op, vec_args_1_0_reduce_zero_p0, Scalar )
{

  const RCP<const Teuchos::Comm<Ordinal> > comm =
    Teuchos::DefaultComm<Ordinal>::getComm();

  const int procRank = rank(*comm);

  const Ordinal localDim = (procRank == 0 ? 0 : g_localDim);

  const Ordinal localOffset = computeLocalOffset(*comm, localDim);
  
  const Scalar val = 1.1;

  RTOpPack::SubVectorView<Scalar> x =
    getLocalSubVectorView<Scalar>(localOffset, localDim, val);

  RTOpPack::ROpSum<Scalar> sumOp;
  RCP<RTOpPack::ReductTarget> sumTarget = sumOp.reduct_obj_create();
  RTOpPack::SPMD_apply_op<Scalar>(&*comm, sumOp, 1, &x, 0, 0, &*sumTarget);
  Scalar sum_x = sumOp(*sumTarget);

  TEST_EQUALITY(sum_x, as<Scalar>(g_localDim * val * (comm->getSize()-1)))

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES(
  SPMD_apply_op, vec_args_1_0_reduce_zero_p0)


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( SPMD_apply_op, vec_args_1_0_reduce_zero_p1, Scalar )
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

  RTOpPack::SubVectorView<Scalar> x =
    getLocalSubVectorView<Scalar>(localOffset, localDim, val);

  RTOpPack::ROpSum<Scalar> sumOp;
  RCP<RTOpPack::ReductTarget> sumTarget = sumOp.reduct_obj_create();
  RTOpPack::SPMD_apply_op<Scalar>(&*comm, sumOp, 1, &x, 0, 0, &*sumTarget);
  Scalar sum_x = sumOp(*sumTarget);

  const Ordinal procFactor = (comm->getSize() == 1 ? 1 : (comm->getSize()-1));
  PRINT_VAR(procFactor);
  TEST_EQUALITY(sum_x, as<Scalar>(g_localDim * val * procFactor))

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES(
  SPMD_apply_op, vec_args_1_0_reduce_zero_p1)


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( SPMD_apply_op, multivec_args_1_0_reduce, Scalar )
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

  RTOpPack::ROpSum<Scalar> sumOp;
  Array<RCP<RTOpPack::ReductTarget> > sumTargets_store(numCols);
  Array<RTOpPack::ReductTarget*> sumTargets(numCols);
  for (int j = 0; j < numCols; ++j) {
    sumTargets_store[j] = sumOp.reduct_obj_create();
    sumTargets[j] = sumTargets_store[j].getRawPtr();
  }
  RTOpPack::SPMD_apply_op<Scalar>(&*comm, sumOp, numCols, 1, &mv, 0, 0,
    sumTargets.getRawPtr());

  for (int j = 0; j < numCols; ++j) {
    Scalar sum_mv_j = sumOp(*sumTargets[j]);
    TEST_EQUALITY(sum_mv_j, as<Scalar>(localDim*val* comm->getSize()));
  }

}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES(
  SPMD_apply_op, multivec_args_1_0_reduce)


} // namespace RTOpPack
