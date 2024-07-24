// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DefaultClusteredSpmdProductVectorSpace.hpp"
#include "Thyra_TestingTools.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Thyra_VectorSpaceTester.hpp"


namespace {


//
// Helper code and declarations
//


using Teuchos::as;
using Teuchos::null;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::get_extra_data;
using Teuchos::ScalarTraits;
using Thyra::VectorSpaceBase;
using Thyra::VectorBase;
using Thyra::MultiVectorBase;
using Thyra::createMember;
using Thyra::createMembers;
using Thyra::DefaultSpmdVectorSpace;
using Thyra::defaultSpmdVectorSpace;
using Thyra::locallyReplicatedDefaultSpmdVectorSpace;
typedef Thyra::Ordinal Ordinal;


const int g_localDim = 4; // ToDo: Make variable!


template<class Scalar>
RCP<VectorSpaceBase<Scalar> > 
createSpmdVectorSpace(const Teuchos_Ordinal localDim)
{
  return defaultSpmdVectorSpace<Scalar>(
    Teuchos::DefaultComm<Teuchos_Ordinal>::getComm(),
    localDim, -1 );
}


//
// Unit Tests
//


TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultClusteredSpmdProductVectorSpace,
  defaultConstruct, Scalar )
{

  RCP<const Thyra::DefaultClusteredSpmdProductVectorSpace<Scalar> > cpvs =
    Teuchos::rcp(new Thyra::DefaultClusteredSpmdProductVectorSpace<Scalar>());
  TEST_EQUALITY(cpvs->intraClusterComm(), null);
  TEST_EQUALITY(cpvs->clusterRootRank(), as<Ordinal>(-1));
  TEST_EQUALITY(cpvs->interClusterComm(), null);
  TEST_EQUALITY(cpvs->clusterSubDim(), as<Ordinal>(-1));
  TEST_EQUALITY(cpvs->clusterOffset(), as<Ordinal>(-1));
  TEST_EQUALITY(cpvs->numBlocks(), as<int>(0));
  TEST_EQUALITY(cpvs->dim(), as<Ordinal>(0));
  
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultClusteredSpmdProductVectorSpace,
  defaultConstruct )



TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL( DefaultClusteredSpmdProductVectorSpace,
  parallelConstruct, Scalar )
{
  const RCP<const Teuchos::Comm<Ordinal> > intraComm =
    Teuchos::DefaultComm<Teuchos_Ordinal>::getComm();
  const RCP<const Teuchos::Comm<Ordinal> > interComm = null;
  const RCP<const Thyra::VectorSpaceBase<Scalar> > vs =
    defaultSpmdVectorSpace<Scalar>(intraComm, g_localDim, -1);
  const RCP<const Thyra::DefaultClusteredSpmdProductVectorSpace<Scalar> > cpvs =
    Teuchos::rcp(
      new Thyra::DefaultClusteredSpmdProductVectorSpace<Scalar>(
        intraComm,
        0,
        interComm,
        1,
        &vs
        )
      );
  TEST_EQUALITY(cpvs->intraClusterComm(), intraComm);
  TEST_EQUALITY(cpvs->clusterRootRank(), as<Ordinal>(0));
  TEST_EQUALITY(cpvs->interClusterComm(), interComm);
  //TEST_EQUALITY(cpvs->clusterSubDim(), as<Ordinal>(-1));
  //TEST_EQUALITY(cpvs->clusterOffset(), as<Ordinal>(-1));
  TEST_EQUALITY(cpvs->numBlocks(), as<int>(1));
  //TEST_EQUALITY(cpvs->dim(), as<Ordinal>(-1));

  TEST_ASSERT(cpvs->isCompatible(*cpvs));
  const RCP<const Thyra::VectorSpaceBase<Scalar> > cpvs_clone = cpvs->clone();
  TEST_ASSERT(cpvs_clone->isCompatible(*cpvs));
  TEST_ASSERT(!cpvs->isCompatible(*vs)); // Needs to be fixed!
  TEST_ASSERT(vs->isCompatible(*cpvs));

  //Thyra::VectorSpaceTester<Scalar> vecSpcTester;
  //TEST_ASSERT(vecSpcTester.check(*cpvs, &out));
}

TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT_SCALAR_TYPES( DefaultClusteredSpmdProductVectorSpace,
  parallelConstruct )


// ToDo: Add a *lot* more testing!


} // namespace
