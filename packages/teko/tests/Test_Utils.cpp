// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Test_Utils.hpp"

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_RCP.hpp"

// Thyra includes
#include "Thyra_SolveSupportTypes.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_VectorSpaceBase.hpp"
#include "Thyra_ProductVectorSpaceBase.hpp"
#include "Thyra_DefaultProductVector.hpp"
#include "Thyra_DefaultProductMultiVector.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_DefaultLinearOpSource.hpp"
#include "Thyra_PreconditionerFactoryHelpers.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_OperatorVectorTypes.hpp"
#include "Thyra_VectorStdOps.hpp"

#include "Teko_Utilities.hpp"

#include "Tpetra_Vector.hpp"
#include "Tpetra_Core.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"

#include <iostream>

using namespace Teuchos;

namespace Teko {
namespace Test {

const RCP<const Thyra::LinearOpBase<ST> > build2x2(const RCP<const Teuchos::Comm<int> > comm, ST a,
                                                   ST b, ST c, ST d) {
  RCP<Tpetra::Map<LO, GO, NT> > map = rcp(new Tpetra::Map<LO, GO, NT>(2, 0, comm));

  GO indices[2];
  ST row0[2];
  ST row1[2];

  indices[0] = 0;
  indices[1] = 1;

  // build a CrsMatrix
  RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > blk = Tpetra::createCrsMatrix<ST, LO, GO, NT>(map, 2);
  row0[0]                                     = a;
  row0[1]                                     = b;  // do a transpose here!
  row1[0]                                     = c;
  row1[1]                                     = d;
  blk->insertGlobalValues(0, Teuchos::ArrayView<GO>(indices, 2), Teuchos::ArrayView<ST>(row0, 2));
  blk->insertGlobalValues(1, Teuchos::ArrayView<GO>(indices, 2), Teuchos::ArrayView<ST>(row1, 2));
  blk->fillComplete();

  return Thyra::tpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(blk->getDomainMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(blk->getRangeMap()), blk);
}

const RCP<const Thyra::MultiVectorBase<ST> > BlockVector(
    const Tpetra::Vector<ST, LO, GO, NT>& tu, const Tpetra::Vector<ST, LO, GO, NT>& tv,
    const RCP<const Thyra::VectorSpaceBase<ST> >& vs) {
  typedef RCP<const Thyra::MultiVectorBase<ST> > Vector;

  const RCP<const Thyra::ProductVectorSpaceBase<ST> > pvs =
      rcp_dynamic_cast<const Thyra::ProductVectorSpaceBase<ST> >(vs);

  RCP<const Tpetra::MultiVector<ST, LO, GO, NT> > mtu = rcpFromRef(tu);
  RCP<const Tpetra::MultiVector<ST, LO, GO, NT> > mtv = rcpFromRef(tv);
  const Vector u = Thyra::createConstMultiVector<ST, LO, GO, NT>(mtu, pvs->getBlock(0));
  const Vector v = Thyra::createConstMultiVector<ST, LO, GO, NT>(mtv, pvs->getBlock(1));

  // build rhs: this is ugly...in 2 steps
  // (i). allocate space for rhs "Product" vector, this is the range of A
  // (j). Need to do a dynamic cast to set the blocks (VectorBase doesn't know about setBlock)
  // const RCP<Thyra::MultiVectorBase<double> > rhs = Thyra::createMembers(vs,1);
  // rcp_dynamic_cast<Thyra::DefaultProductMultiVector<double> >(rhs)->setBlock(0,u);
  // rcp_dynamic_cast<Thyra::DefaultProductMultiVector<double> >(rhs)->setBlock(1,v);
  std::vector<RCP<Thyra::MultiVectorBase<ST> > > blocks;
  blocks.push_back(rcp_const_cast<Thyra::MultiVectorBase<ST> >(u));
  blocks.push_back(rcp_const_cast<Thyra::MultiVectorBase<ST> >(v));

  return buildBlockedMultiVector(blocks);
}

void Print(std::ostream& os, const std::string& s,
           const RCP<const Thyra::MultiVectorBase<double> >& v) {
  Thyra::ConstDetachedVectorView<double> view(v->col(0));

  os << s << " = ";
  for (int i = 0; i < v->range()->dim(); i++) os << view[i] << " ";
  os << std::endl;
}

RCP<Thyra::VectorBase<double> > BuildIVector(int j,
                                             const RCP<const Thyra::VectorSpaceBase<double> >& vs) {
  RCP<Thyra::VectorBase<double> > v = Thyra::createMember(vs);
  Thyra::assign<double>(v.ptr(), 0.0);

  {
    Thyra::DetachedVectorView<double> view(v);
    view[j] = 1.0;
  }

  return v;
}

void HardExtract(std::ostream& os, const RCP<const Thyra::LinearOpBase<double> >& A) {
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;

  int ncols = A->domain()->dim();
  int nrows = A->range()->dim();

  // print the matrix the hard way
  for (int i = 0; i < nrows; i++) {
    const RCP<Thyra::VectorBase<double> > u = Teko::Test::BuildIVector(i, A->range());
    for (int j = 0; j < ncols; j++) {
      const RCP<const Thyra::VectorBase<double> > ej = Teko::Test::BuildIVector(j, A->domain());
      const RCP<Thyra::VectorBase<double> > v        = Thyra::createMember(A->range());
      ::Thyra::apply(*A, Thyra::NOTRANS, *ej, v.ptr());

      // this for some reason messes up valgrind (2/9/09)
      double value = A->range()->scalarProd(*u, *v);

      os << value << "  ";
    }
    os << std::endl;
  }
}

// add two Thyra vectors
const Teuchos::RCP<Thyra::MultiVectorBase<double> > Add(
    const Teuchos::RCP<const Thyra::MultiVectorBase<double> >& x,
    const Teuchos::RCP<const Thyra::MultiVectorBase<double> >& y) {
  return Add(1.0, x, 1.0, y);
}

const Teuchos::RCP<Thyra::MultiVectorBase<double> > Add(
    double ax, const Teuchos::RCP<const Thyra::MultiVectorBase<double> >& x, double ay,
    const Teuchos::RCP<const Thyra::MultiVectorBase<double> >& y) {
  const Teuchos::RCP<Thyra::MultiVectorBase<double> > z = Thyra::createMembers(x->range(), 1);

  Thyra::linear_combination<double>(Teuchos::tuple<double>(ax, ay),
                                    Teuchos::tuple(Teuchos::ptrInArg(*x), Teuchos::ptrInArg(*y)),
                                    0.0, z.ptr());

  return z;
}

// compute ||x-y||_2
double Difference(const Teuchos::RCP<const Thyra::MultiVectorBase<double> >& x,
                  const Teuchos::RCP<const Thyra::MultiVectorBase<double> >& y) {
  Teuchos::RCP<const Thyra::MultiVectorBase<double> > z = Add(1.0, x, -1.0, y);

  return Thyra::norm_2(*z->col(0));
}

// construct a diagonal matrix

const Teuchos::RCP<const Thyra::LinearOpBase<double> > DiagMatrix_tpetra(GO cnt, ST* vec,
                                                                         std::string label) {
  const RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const RCP<const Tpetra::Map<LO, GO, NT> > map =
      rcp(new const Tpetra::Map<LO, GO, NT>(cnt, 0, comm));
  const RCP<Tpetra::CrsMatrix<ST, LO, GO, NT> > ptrF =
      Tpetra::createCrsMatrix<ST, LO, GO, NT>(map, 1);

  // construct a diagonal matrix
  GO indices[1];
  ST values[1];
  for (int i = 0; i < cnt; i++) {
    indices[0] = i;
    values[0]  = vec[i];
    ptrF->insertGlobalValues(i, Teuchos::ArrayView<GO>(indices, 1),
                             Teuchos::ArrayView<ST>(values, 1));
  }
  ptrF->fillComplete();

  // return thyra object
  return Thyra::tpetraLinearOp<ST, LO, GO, NT>(
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(ptrF->getDomainMap()),
      Thyra::tpetraVectorSpace<ST, LO, GO, NT>(ptrF->getRangeMap()), ptrF);
}

// declare static allocation
std::list<std::pair<Teuchos::RCP<UnitTest>, std::string> > UnitTest::testList;
Teuchos::RCP<const Teuchos::Comm<int> > UnitTest::comm_tpetra_;

void UnitTest::AddTest(const Teuchos::RCP<UnitTest>& ut, const std::string& name) {
  // add a unit test and string to the list
  testList.push_back(std::make_pair(ut, name));
}

bool UnitTest::RunTests_tpetra(int verbosity, std::ostream& stdstrm, std::ostream& failstrm) {
  bool allPassed = true;
  int count = 0, numfailed = 0;
  std::list<std::pair<Teuchos::RCP<UnitTest>, std::string> >::iterator itr;

  bool isParallel = GetComm_tpetra()->getSize() > 1;

  // loop over the tests and run each
  for (itr = testList.begin(); itr != testList.end(); ++itr) {
    int localrun  = 0;
    int localfail = 0;

    // skip a test if its not parallel
    if (isParallel && not(itr->first)->isParallel()) continue;

    stdstrm << "Running test \"" << itr->second << (verbosity >= 1 ? "\"\n" : "\" ... ");

    // run the tests
    (itr->first)->initializeTest();
    bool status = 0 == (localfail = (itr->first)->runTest(verbosity, stdstrm, failstrm, localrun));

    // output some stuff to the standard stream
    if (verbosity >= 1) {
      stdstrm << "Test \"" << itr->second << "\" completed ... ";
      if (status)
        stdstrm << Teko::Test::toString(status) << " (" << localrun << ")" << std::endl;
      else
        stdstrm << Teko::Test::toString(status) << " (" << localfail << ")" << std::endl;
    }

    allPassed &= status;
    if (not status) numfailed += localfail;
    count += localrun;
  }

  // output status
  stdstrm << std::endl;
  stdstrm << "Tests Passed: " << count - numfailed << ", Tests Failed: " << numfailed << std::endl;
  stdstrm << "(Incidently, you want no failures)" << std::endl;

  return allPassed;
}

Teuchos::RCP<const Teuchos::Comm<int> > UnitTest::GetComm_tpetra() { return comm_tpetra_; }

void UnitTest::SetComm_tpetra(const Teuchos::RCP<const Teuchos::Comm<int> >& c) {
  comm_tpetra_ = c;
}

void UnitTest::ClearTests() {
  testList.clear();
  comm_tpetra_ = Teuchos::null;
}

bool UnitTest::CheckParallelBools_tpetra(bool myBool, int& failPID) {
  char myInt = myBool ? 1 : 0;
  std::vector<char> bools(GetComm_tpetra()->getSize());

  GetComm_tpetra()->gatherAll(1, &myInt, GetComm_tpetra()->getSize(), &bools[0]);

  failPID = -1;

  bool result = true;
  for (unsigned int i = 0; i < bools.size(); i++) {
    result &= bools[i] == 1 ? true : false;
    if (bools[i] != 1) failPID = i;
  }

  return result;
}

}  // namespace Test
}  // namespace Teko
