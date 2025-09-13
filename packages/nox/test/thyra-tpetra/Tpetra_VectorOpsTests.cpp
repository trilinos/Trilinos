// @HEADER
// *****************************************************************************
//            NOX: An Object-Oriented Nonlinear Solver Package
//
// Copyright 2002 NTESS and the NOX contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_RCP.hpp"

#include "Tpetra_Core.hpp"
#include "Tpetra_Vector.hpp"
#include "NOX_Thyra_Vector.H"
#include "Thyra_TpetraThyraWrappers.hpp"

// Typedefs

typedef Tpetra::Vector<>::scalar_type Scalar;
typedef Tpetra::Vector<>::local_ordinal_type LO;
typedef Tpetra::Vector<>::global_ordinal_type GO;
typedef Tpetra::Vector<>::node_type Node;

typedef Teuchos::ScalarTraits<Scalar> ST;
typedef Tpetra::Map<LO,GO,Node> Map;
typedef Tpetra::Vector<Scalar,LO,GO,Node> TV;
typedef Thyra::VectorBase<Scalar> TVB;
typedef NOX::Thyra::Vector NTV;
typedef typename TV::mag_type mag_type;

#if defined HAVE_TPETRACORE_CUDA
#define NUM_LOCAL 1000000
#else
#define NUM_LOCAL 1000
#endif
const std::size_t numLocalElements = NUM_LOCAL;

//Routines for checking solution

bool checkVectors(const Teuchos::RCP<TV>& a,
                  const Teuchos::RCP<TV>& b,
                  const mag_type& expectedNorm,
                  Teuchos::FancyOStream& out)
{
  bool success = true;
  ST::magnitudeType tol = 1.0e-14;
  b->update(-1.0*ST::one(), *a, ST::one());
  TEUCHOS_TEST_EQUALITY(b->norm2() < tol, true, out, success);
  TEUCHOS_TEST_FLOATING_EQUALITY(a->norm2(), expectedNorm, tol, out, success);
  return success;
}

template <class T>
bool checkReduction(const T& a,
                    const T& b,
                    const T& expectedVal,
                    Teuchos::FancyOStream& out)
{
  bool success = true;
  typename Teuchos::ScalarTraits<T>::magnitudeType tol = 1.0e-14;
  TEUCHOS_TEST_FLOATING_EQUALITY(expectedVal, a, tol, out, success);
  TEUCHOS_TEST_FLOATING_EQUALITY(expectedVal, b, tol, out, success);
  return success;
}

// Unit tests

TEUCHOS_UNIT_TEST(Tpetra_VectorOps, CopyConstructor)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;

  // Create Tpetra vectors
  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<TV> x = Teuchos::rcp(new TV(map));

  x->putScalar(ST::one());

  // Create new NOX vector with copy of Thyra vector
  Teuchos::RCP<TVB> x_thyra = Thyra::createVector(x);
  Teuchos::RCP<NTV> y = Teuchos::rcp(new NTV(*x_thyra));

  // Check for correct answer
  typedef Thyra::TpetraOperatorVectorExtraction<Scalar, LO, GO, Node> TOVE;
  mag_type ans = static_cast<mag_type>(ST::squareroot(numGlobalElements));
  success = checkVectors(x, TOVE::getTpetraVector(y->getThyraRCPVector()), ans, out);
}

TEUCHOS_UNIT_TEST(Tpetra_VectorOps, OperatorEquals)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;

  // Create Tpetra vectors
  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<TV> x = Teuchos::rcp(new TV(map));
  Teuchos::RCP<TV> y = Teuchos::rcp(new TV(map));

  x->putScalar(ST::one());
  y->putScalar(ST::zero());

  // Copy x into y through NOX interface
  Teuchos::RCP<TVB> x_thyra = Thyra::createVector(x);
  Teuchos::RCP<NTV> x_nox = Teuchos::rcp(new NTV(x_thyra));
  Teuchos::RCP<TVB> y_thyra = Thyra::createVector(y);
  Teuchos::RCP<NTV> y_nox = Teuchos::rcp(new NTV(y_thyra));

  *y_nox = *x_nox;

  // Check for correct answer
  mag_type ans = static_cast<mag_type>(ST::squareroot(numGlobalElements));;
  success = checkVectors(x, y, ans, out);
}

TEUCHOS_UNIT_TEST(Tpetra_VectorOps, Init)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;

  // Create Tpetra vectors
  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<TV> x = Teuchos::rcp(new TV(map));
  Teuchos::RCP<TV> y = Teuchos::rcp(new TV(map));

  // Perform operation with Tpetra directly
  x->putScalar(2.0*ST::one());

  // Do the same thing through NOX/Thyra interface
  Teuchos::RCP<TVB> y_thyra = Thyra::createVector(y);
  Teuchos::RCP<NTV> y_nox = Teuchos::rcp(new NTV(y_thyra));

  y_nox->init(2.0*ST::one());

  // Check for correct answer
  mag_type ans = static_cast<mag_type>(ST::squareroot(4.0*numGlobalElements));;
  success = checkVectors(x, y, ans, out);
}

TEUCHOS_UNIT_TEST(Tpetra_VectorOps, Random)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;

  // Create Tpetra vectors
  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<TV> x = Teuchos::rcp(new TV(map));
  Teuchos::RCP<TV> y = Teuchos::rcp(new TV(map));

  // Perform operation with Tpetra directly
  x->randomize();

  // Do the same thing through NOX/Thyra interface
  Teuchos::RCP<TVB> y_thyra = Thyra::createVector(y);
  Teuchos::RCP<NTV> y_nox = Teuchos::rcp(new NTV(y_thyra));

  y_nox->random();

  // Probaby can't expect these to get the same answer
}

TEUCHOS_UNIT_TEST(Tpetra_VectorOps, Abs)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;

  // Create Tpetra vectors
  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<TV> x = Teuchos::rcp(new TV(map));
  Teuchos::RCP<TV> y = Teuchos::rcp(new TV(map));

  x->putScalar(-1.0*ST::one());
  y->putScalar(-1.0*ST::one());

  // Perform operation with Tpetra directly
  x->abs(*x);

  // Do the same thing through NOX/Thyra interface
  Teuchos::RCP<TVB> y_thyra = Thyra::createVector(y);
  Teuchos::RCP<NTV> y_nox = Teuchos::rcp(new NTV(y_thyra));

  y_nox->abs(*y_nox);

  // Check for correct answer
  mag_type ans = static_cast<mag_type>(ST::squareroot(numGlobalElements));;
  success = checkVectors(x, y, ans, out);
}

TEUCHOS_UNIT_TEST(Tpetra_VectorOps, Reciprocal)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;

  // Create Tpetra vectors
  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<TV> x = Teuchos::rcp(new TV(map));
  Teuchos::RCP<TV> y = Teuchos::rcp(new TV(map));

  x->putScalar(2.0*ST::one());
  y->putScalar(2.0*ST::one());

  // Perform operation with Tpetra directly
  x->reciprocal(*x);

  // Do the same thing through NOX/Thyra interface
  Teuchos::RCP<TVB> y_thyra = Thyra::createVector(y);
  Teuchos::RCP<NTV> y_nox = Teuchos::rcp(new NTV(y_thyra));

  y_nox->reciprocal(*y_nox);

  // Check for correct answer
  mag_type ans = static_cast<mag_type>(ST::squareroot(0.25*numGlobalElements));;
  success = checkVectors(x, y, ans, out);
}

TEUCHOS_UNIT_TEST(Tpetra_VectorOps, Scale_Scalar)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;

  // Create Tpetra vectors
  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<TV> x = Teuchos::rcp(new TV(map));
  Teuchos::RCP<TV> y = Teuchos::rcp(new TV(map));

  x->putScalar(ST::one());
  y->putScalar(ST::one());

  // Perform operation with Tpetra directly
  x->scale(2.0*ST::one());

  // Do the same thing through NOX/Thyra interface
  Teuchos::RCP<TVB> y_thyra = Thyra::createVector(y);
  Teuchos::RCP<NTV> y_nox = Teuchos::rcp(new NTV(y_thyra));

  y_nox->scale(2.0*ST::one());

  // Check for correct answer
  mag_type ans = static_cast<mag_type>(ST::squareroot(4.0*numGlobalElements));;
  success = checkVectors(x, y, ans, out);
}

TEUCHOS_UNIT_TEST(Tpetra_VectorOps, Scale_Vector)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;

  // Create Tpetra vectors
  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<TV> x = Teuchos::rcp(new TV(map));
  Teuchos::RCP<TV> y = Teuchos::rcp(new TV(map));
  Teuchos::RCP<TV> z = Teuchos::rcp(new TV(map));

  x->putScalar(ST::one());
  y->putScalar(ST::one());
  z->putScalar(2.0*ST::one());

  // Perform operation with Tpetra directly
  x->elementWiseMultiply(ST::one(), *z, *x, ST::zero());

  // Do the same thing through NOX/Thyra interface
  Teuchos::RCP<TVB> y_thyra = Thyra::createVector(y);
  Teuchos::RCP<NTV> y_nox = Teuchos::rcp(new NTV(y_thyra));
  Teuchos::RCP<TVB> z_thyra = Thyra::createVector(z);
  Teuchos::RCP<NTV> z_nox = Teuchos::rcp(new NTV(z_thyra));

  y_nox->scale(*z_nox);

  // Check for correct answer
  mag_type ans = static_cast<mag_type>(ST::squareroot(4.0*numGlobalElements));;
  success = checkVectors(x, y, ans, out);
}

TEUCHOS_UNIT_TEST(Tpetra_VectorOps, Update_1)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;

  // Create Tpetra vectors
  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<TV> x = Teuchos::rcp(new TV(map));
  Teuchos::RCP<TV> y = Teuchos::rcp(new TV(map));
  Teuchos::RCP<TV> z = Teuchos::rcp(new TV(map));

  x->putScalar(ST::one());
  y->putScalar(ST::one());
  z->putScalar(ST::one());

  // Perform operation with Tpetra directly
  x->update(2.0*ST::one(), *z, ST::one());

  // Do the same thing through NOX/Thyra interface
  Teuchos::RCP<TVB> y_thyra = Thyra::createVector(y);
  Teuchos::RCP<NTV> y_nox = Teuchos::rcp(new NTV(y_thyra));
  Teuchos::RCP<TVB> z_thyra = Thyra::createVector(z);
  Teuchos::RCP<NTV> z_nox = Teuchos::rcp(new NTV(z_thyra));

  y_nox->update(2.0*ST::one(), *z_nox, ST::one());

  // Check for correct answer
  mag_type ans = static_cast<mag_type>(ST::squareroot(9.0*numGlobalElements));;
  success = checkVectors(x, y, ans, out);
}

TEUCHOS_UNIT_TEST(Tpetra_VectorOps, Update_2)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;

  // Create Tpetra vectors
  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<TV> w = Teuchos::rcp(new TV(map));
  Teuchos::RCP<TV> x = Teuchos::rcp(new TV(map));
  Teuchos::RCP<TV> y = Teuchos::rcp(new TV(map));
  Teuchos::RCP<TV> z = Teuchos::rcp(new TV(map));

  w->putScalar(ST::one());
  x->putScalar(ST::one());
  y->putScalar(ST::one());
  z->putScalar(ST::one());

  // Perform operation with Tpetra directly
  x->update(2.0*ST::one(), *w, 2.0*ST::one(), *z, ST::one());

  // Do the same thing through NOX/Thyra interface
  Teuchos::RCP<TVB> w_thyra = Thyra::createVector(w);
  Teuchos::RCP<NTV> w_nox = Teuchos::rcp(new NTV(w_thyra));
  Teuchos::RCP<TVB> y_thyra = Thyra::createVector(y);
  Teuchos::RCP<NTV> y_nox = Teuchos::rcp(new NTV(y_thyra));
  Teuchos::RCP<TVB> z_thyra = Thyra::createVector(z);
  Teuchos::RCP<NTV> z_nox = Teuchos::rcp(new NTV(z_thyra));

  y_nox->update(2.0*ST::one(), *w_nox, 2.0*ST::one(), *z_nox, ST::one());

  // Check for correct answer
  mag_type ans = static_cast<mag_type>(ST::squareroot(25.0*numGlobalElements));;
  success = checkVectors(x, y, ans, out);
}

TEUCHOS_UNIT_TEST(Tpetra_VectorOps, Norm_1)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;

  // Create Tpetra vectors
  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<TV> x = Teuchos::rcp(new TV(map));

  x->putScalar(ST::one());

  // Perform operation with Tpetra directly
  ST::magnitudeType norm_tpetra = x->norm1();

  // Do the same thing through NOX/Thyra interface
  Teuchos::RCP<TVB> x_thyra = Thyra::createVector(x);
  Teuchos::RCP<NTV> x_nox = Teuchos::rcp(new NTV(x_thyra));

  ST::magnitudeType norm_nox = x_nox->norm(NOX::Abstract::Vector::OneNorm);

  // Check for correct answer
  ST::magnitudeType ans = static_cast<mag_type>(numGlobalElements);;
  success = checkReduction(norm_tpetra, norm_nox, ans, out);
}

TEUCHOS_UNIT_TEST(Tpetra_VectorOps, Norm_2)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;

  // Create Tpetra vectors
  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<TV> x = Teuchos::rcp(new TV(map));

  x->putScalar(ST::one());

  // Perform operation with Tpetra directly
  ST::magnitudeType norm_tpetra = x->norm2();

  // Do the same thing through NOX/Thyra interface
  Teuchos::RCP<TVB> x_thyra = Thyra::createVector(x);
  Teuchos::RCP<NTV> x_nox = Teuchos::rcp(new NTV(x_thyra));

  ST::magnitudeType norm_nox = x_nox->norm(NOX::Abstract::Vector::TwoNorm);

  // Check for correct answer
  ST::magnitudeType ans = static_cast<mag_type>(ST::squareroot(numGlobalElements));;
  success = checkReduction(norm_tpetra, norm_nox, ans, out);
}

TEUCHOS_UNIT_TEST(Tpetra_VectorOps, Norm_Inf)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;

  // Create Tpetra vectors
  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<TV> x = Teuchos::rcp(new TV(map));

  x->putScalar(ST::one());

  // Perform operation with Tpetra directly
  ST::magnitudeType norm_tpetra = x->normInf();

  // Do the same thing through NOX/Thyra interface
  Teuchos::RCP<TVB> x_thyra = Thyra::createVector(x);
  Teuchos::RCP<NTV> x_nox = Teuchos::rcp(new NTV(x_thyra));

  ST::magnitudeType norm_nox = x_nox->norm(NOX::Abstract::Vector::MaxNorm);

  // Check for correct answer
  ST::magnitudeType ans = static_cast<mag_type>(1.0);;
  success = checkReduction(norm_tpetra, norm_nox, ans, out);
}

TEUCHOS_UNIT_TEST(Tpetra_VectorOps, Norm_Weighted)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;

  // Create Tpetra vectors
  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<TV> x = Teuchos::rcp(new TV(map));
  Teuchos::RCP<TV> y = Teuchos::rcp(new TV(map));

  x->putScalar(ST::one());
  y->putScalar(2.0*ST::one());

  // Tpetra weighted vector norm is deprecated, so need to do weighting manually
  ST::magnitudeType norm_tpetra;
  {
   Teuchos::RCP<TV> copy = Teuchos::rcp(new TV(*x, Teuchos::Copy));
   copy->elementWiseMultiply(ST::one(), *y, *copy, ST::zero());
   norm_tpetra = ST::magnitude(ST::squareroot(copy->dot(*x)));
  }

  // Do the same thing through NOX/Thyra interface
  Teuchos::RCP<TVB> x_thyra = Thyra::createVector(x);
  Teuchos::RCP<NTV> x_nox = Teuchos::rcp(new NTV(x_thyra));
  Teuchos::RCP<TVB> y_thyra = Thyra::createVector(y);
  Teuchos::RCP<NTV> y_nox = Teuchos::rcp(new NTV(y_thyra));

  ST::magnitudeType norm_nox = x_nox->norm(*y_nox);

  // Check for correct answer
  ST::magnitudeType ans = static_cast<mag_type>(ST::squareroot(2.0*numGlobalElements));;
  success = checkReduction(norm_tpetra, norm_nox, ans, out);
}

TEUCHOS_UNIT_TEST(Tpetra_VectorOps, InnerProduct)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;

  // Create Tpetra vectors
  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<TV> x = Teuchos::rcp(new TV(map));
  Teuchos::RCP<TV> y = Teuchos::rcp(new TV(map));

  x->putScalar(ST::one());
  y->putScalar(2.0*ST::one());

  // Perform operation with Tpetra directly
  Scalar dot_tpetra = x->dot(*y);

  // Do the same thing through NOX/Thyra interface
  Teuchos::RCP<TVB> x_thyra = Thyra::createVector(x);
  Teuchos::RCP<NTV> x_nox = Teuchos::rcp(new NTV(x_thyra));
  Teuchos::RCP<TVB> y_thyra = Thyra::createVector(y);
  Teuchos::RCP<NTV> y_nox = Teuchos::rcp(new NTV(y_thyra));

  Scalar dot_nox = x_nox->innerProduct(*y_nox);

  // Check for correct answer
  Scalar ans = static_cast<mag_type>(2.0*numGlobalElements);;
  success = checkReduction(dot_tpetra, dot_nox, ans, out);
}

