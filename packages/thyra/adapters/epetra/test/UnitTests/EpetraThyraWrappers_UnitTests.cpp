// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_SpmdVectorSpaceBase.hpp"
#include "Thyra_VectorSpaceTester.hpp"
#include "Thyra_TestingTools.hpp"
#include "Epetra_Vector.h"
#include "Teuchos_GlobalMPISession.hpp"

#include "EpetraThyraAdaptersTestHelpers.hpp"

#include "Teuchos_UnitTestHarness.hpp"


namespace {


using Teuchos::ptrFromRef;
using Teuchos::Comm;
using Teuchos::rcp_dynamic_cast;
using Teuchos::rcp;


void createEpetraVsAndMap(const Thyra::Ordinal localDim_in,
  const Ptr<RCP<const Thyra::VectorSpaceBase<double> > > &vs,
  const Ptr<RCP<const Epetra_Map> > &epetra_map,
  const int emptyProcRootRank = -1
  )
{
  const RCP<const Epetra_Comm> epetra_comm = getEpetraComm();
  const int procRank = epetra_comm->MyPID();
  const Thyra::Ordinal localDim = (procRank == emptyProcRootRank ? 0 : localDim_in);
  *epetra_map = rcp(new Epetra_Map(-1, as<int>(localDim), 0, *epetra_comm));
  *vs =  Thyra::create_VectorSpace(*epetra_map);
}


void runVectorSpaceTesterTest(const int emptyProc, 
  Teuchos::FancyOStream &out, bool &success)
{
  using Thyra::VectorSpaceBase;
  using Thyra::SpmdVectorSpaceBase;
  using Thyra::MultiVectorBase;
  
  RCP<const VectorSpaceBase<double> > vs;
  RCP<const Epetra_Map> epetra_map;
  createEpetraVsAndMap(g_localDim, outArg(vs), outArg(epetra_map), emptyProc);
  const int numProcs = epetra_map->Comm().NumProc();

  if (emptyProc >= numProcs) {
    out << "emptyProc = " << emptyProc << " >= numProcs = " << numProcs
        << ": Skipping this test case!\n";
    return;
  }

  const Ordinal dimMultiplier = (emptyProc < 0 ? numProcs : numProcs-1);
  
  TEST_EQUALITY(vs->dim(), g_localDim * dimMultiplier);

  const RCP<const SpmdVectorSpaceBase<double> > spmd_vs =
    rcp_dynamic_cast<const SpmdVectorSpaceBase<double> >(vs, true);

  TEST_EQUALITY(spmd_vs->localSubDim(), as<int>(epetra_map->NumMyElements()));
  
  Thyra::VectorSpaceTester<double> vectorSpaceTester;
  const double tol = 100.0 * Teuchos::ScalarTraits<double>::eps();
  vectorSpaceTester.warning_tol((0.1)*tol);
  vectorSpaceTester.error_tol(tol);
  vectorSpaceTester.show_all_tests(g_show_all_tests);
  vectorSpaceTester.dump_all(g_dumpAll);
  TEST_ASSERT(vectorSpaceTester.check(*vs, &out));
}


void runCreateVectorUnitTest(const int emptyProc, 
  Teuchos::FancyOStream &out, bool &success)
{
  using Thyra::VectorBase;
  using Thyra::VectorSpaceBase;
  using Thyra::MultiVectorBase;
  
  RCP<const VectorSpaceBase<double> > vs;
  RCP<const Epetra_Map> epetra_map;
  createEpetraVsAndMap(g_localDim, outArg(vs), outArg(epetra_map), emptyProc);
  const int numProcs = epetra_map->Comm().NumProc();

  if (emptyProc >= numProcs) {
    out << "emptyProc = " << emptyProc << " >= numProcs = " << numProcs
        << ": Skipping this test case!\n";
    return;
  }

  const RCP<Epetra_Vector> epetra_vec = rcp(new Epetra_Vector(*epetra_map));
  const RCP<VectorBase<double> > thyra_vec = Thyra::create_Vector(epetra_vec, vs);
  const RCP<Epetra_Vector> epetra_vec2 =
    Thyra::get_Epetra_Vector(*epetra_map, thyra_vec);
  TEST_EQUALITY(epetra_vec, epetra_vec2);

  const RCP<const Epetra_Vector> const_epetra_vec2 =
    Thyra::get_Epetra_Vector(*epetra_map, thyra_vec.getConst());
  TEST_EQUALITY(epetra_vec.getConst(), const_epetra_vec2);

}


//
// Unit Tests
//


TEUCHOS_UNIT_TEST( EpetraThyraWrappers, comm )
{
  typedef Teuchos::GlobalMPISession GMPIS;
  const RCP<const Epetra_Comm> epetra_comm = getEpetraComm();
  TEST_EQUALITY(epetra_comm->NumProc(), GMPIS::getNProc());
  TEST_EQUALITY(epetra_comm->MyPID(), GMPIS::getRank());
  const RCP<const Comm<Ordinal> > comm = Thyra::create_Comm(epetra_comm);
  TEST_EQUALITY(comm->getSize(), GMPIS::getNProc());
  TEST_EQUALITY(comm->getRank(), GMPIS::getRank());
  const RCP<const Epetra_Comm> epetra_comm2 = Thyra::get_Epetra_Comm(*comm);
  TEST_EQUALITY(epetra_comm2->NumProc(), GMPIS::getNProc());
  TEST_EQUALITY(epetra_comm2->MyPID(), GMPIS::getRank());
}


TEUCHOS_UNIT_TEST( EpetraThyraWrappers, vectorSpaceTester )
{
  runVectorSpaceTesterTest(-1, out, success);
}


TEUCHOS_UNIT_TEST( EpetraThyraWrappers, vectorSpaceTester_empty_p0 )
{
  runVectorSpaceTesterTest(0, out, success);
}


TEUCHOS_UNIT_TEST( EpetraThyraWrappers, vectorSpaceTester_empty_p1 )
{
  runVectorSpaceTesterTest(1, out, success);
}


TEUCHOS_UNIT_TEST( EpetraThyraWrappers, vectorSpaceTester_empty_pLast )
{
  runVectorSpaceTesterTest(Teuchos::GlobalMPISession::getNProc()-1, out, success);
}


TEUCHOS_UNIT_TEST( EpetraThyraWrappers, createAndViewVector_empty_p0 )
{
  runCreateVectorUnitTest(0, out, success);
}


TEUCHOS_UNIT_TEST( EpetraThyraWrappers, createAndViewVector_empty_p1 )
{
  runCreateVectorUnitTest(1, out, success);
}


TEUCHOS_UNIT_TEST( EpetraThyraWrappers, createAndViewVector_empty_pLast )
{
  runCreateVectorUnitTest(Teuchos::GlobalMPISession::getNProc()-1, out, success);
}


// Test Thyra::create_MultiVector() (const and nonconst) and all views with empty processes


TEUCHOS_UNIT_TEST( EpetraThyraWrappers, get_Epetra_MultiVector_singleBlockProductVector )
{
using Thyra::VectorSpaceBase;
using Thyra::MultiVectorBase;

  RCP<const VectorSpaceBase<double> > vs;
  RCP<const Epetra_Map> epetra_map;
  createEpetraVsAndMap(g_localDim, outArg(vs), outArg(epetra_map));

  const RCP<const VectorSpaceBase<double> > pvs = Thyra::productVectorSpace(vs, 1);

  const RCP<MultiVectorBase<double> > pmv = Thyra::createMembers(pvs, 1);

  const double alpha = 3.5;
  Thyra::assign<double>( pmv.ptr(), alpha );

  const RCP<Epetra_MultiVector> epetra_mv =
    Thyra::get_Epetra_MultiVector(*epetra_map, pmv);

  const RCP<MultiVectorBase<double> > mv2 =
    Thyra::create_MultiVector(epetra_mv, pvs);

  Thyra::testRelNormDiffErr<double>(
    "*pmv->col(0)", *pmv->col(0),
    "*mv2->col(0)", *mv2->col(0),
    "max-error", 0.0,
    "max-warning", 0.0,
    &out
    );
   
}


} // namespace
