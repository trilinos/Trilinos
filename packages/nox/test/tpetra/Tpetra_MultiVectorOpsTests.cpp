#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_RCP.hpp"

#include "Tpetra_Core.hpp"
#include "Tpetra_MultiVector.hpp"
#include "NOX_Thyra_MultiVector.H"
#include "Thyra_TpetraThyraWrappers.hpp"

// Typedefs and constants

typedef Tpetra::MultiVector<>::scalar_type Scalar;
typedef Tpetra::MultiVector<>::local_ordinal_type LO;
typedef Tpetra::MultiVector<>::global_ordinal_type GO;
typedef Tpetra::MultiVector<>::node_type Node;

typedef Teuchos::ScalarTraits<Scalar> ST;
typedef Tpetra::Map<LO,GO,Node> Map;
typedef Tpetra::MultiVector<Scalar,LO,GO,Node> TMV;
typedef Thyra::MultiVectorBase<Scalar> TMVB;
typedef NOX::Thyra::MultiVector NTMV;

#if defined HAVE_TPETRACORE_CUDA
#define NUM_LOCAL 1000000
#else
#define NUM_LOCAL 1000
#endif
const std::size_t numLocalElements = NUM_LOCAL;
const std::size_t numCols = 10;

//Routines for checking solution

bool checkEqualMultiVectors(const Teuchos::RCP<TMV>& a,
                            const Teuchos::RCP<TMV>& b,
                            Teuchos::FancyOStream& out)
{
  TEUCHOS_ASSERT((a->getNumVectors() == numCols) &&  (b->getNumVectors() == numCols));
  Scalar one = ST::one();
  b->update(-1.0*one, *a, one);
  Teuchos::Tuple<ST::magnitudeType, numCols> norms;
  b->norm2(norms);

  bool success = true;
  ST::magnitudeType tol = 1.0e-14;
  ST::magnitudeType zero = ST::magnitude(ST::zero());
  for (auto norm = norms.begin(); norm != norms.end(); ++norm) {
    TEUCHOS_TEST_FLOATING_EQUALITY(*norm, zero, tol, out, success);
  }
  return success;
}

template <class T>
bool checkEqualReductions(const Teuchos::ArrayView<T>& a,
                          const Teuchos::ArrayView<T>& b,
                          Teuchos::FancyOStream& out)
{
  TEUCHOS_ASSERT(a.size() == b.size());

  bool success = true;
  typename Teuchos::ScalarTraits<T>::magnitudeType tol = 1.0e-14;
  auto a_itr = a.begin();
  auto b_itr = b.begin();
  for (; a_itr != a.end(); ++a_itr, ++b_itr) {
    TEUCHOS_TEST_FLOATING_EQUALITY(*a_itr, *b_itr, tol, out, success);
  }
  return success;
}

// Unit Tests

TEUCHOS_UNIT_TEST(Tpetra_MultiVectorOps, CopyConstructor)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;

  // Create Tpetra multivectors
  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<TMV> x = Teuchos::rcp(new TMV(map, numCols));

  x->putScalar(ST::one());

  // Create new NOX multivector with copy of Thyra multivector
  Teuchos::RCP<TMVB> x_thyra  = Thyra::createMultiVector(x);
  Teuchos::RCP<NTMV> y = Teuchos::rcp(new NTMV(*x_thyra));

  // Check for correct answer
  typedef Thyra::TpetraOperatorVectorExtraction<Scalar, LO, GO, Node> TOVE;
  success = checkEqualMultiVectors(x, TOVE::getTpetraMultiVector(y->getThyraMultiVector()), out);
}

TEUCHOS_UNIT_TEST(Tpetra_MultiVectorOps, OperatorEquals)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;

  // Create Tpetra multivectors
  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<TMV> x = Teuchos::rcp(new TMV(map, numCols));
  Teuchos::RCP<TMV> y = Teuchos::rcp(new TMV(map, numCols));

  x->putScalar(ST::one());
  y->putScalar(ST::zero());

  // Copy x into y through NOX interface
  Teuchos::RCP<TMVB> x_thyra  = Thyra::createMultiVector(x);
  Teuchos::RCP<NTMV> x_nox = Teuchos::rcp(new NTMV(x_thyra));
  Teuchos::RCP<TMVB> y_thyra  = Thyra::createMultiVector(y);
  Teuchos::RCP<NTMV> y_nox = Teuchos::rcp(new NTMV(y_thyra));

  // This is causing a link error for some reason???
  //*y_nox = *x_nox;

  // Check for correct answer
  //success = checkEqualMultiVectors(x, y, out);
}

TEUCHOS_UNIT_TEST(Tpetra_MultiVectorOps, Init)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;

  // Create Tpetra multivectors
  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<TMV> x = Teuchos::rcp(new TMV(map, numCols));
  Teuchos::RCP<TMV> y = Teuchos::rcp(new TMV(map, numCols));

  // Perform operation with Tpetra directly
  x->putScalar(2.0*ST::one());

  // Do the same thing through NOX/Thyra interface
  Teuchos::RCP<TMVB> y_thyra  = Thyra::createMultiVector(y);
  Teuchos::RCP<NTMV> y_nox = Teuchos::rcp(new NTMV(y_thyra));

  y_nox->init(2.0*ST::one());

  // Check for correct answer
  success = checkEqualMultiVectors(x, y, out);
}

TEUCHOS_UNIT_TEST(Tpetra_MultiVectorOps, Random)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;

  // Create Tpetra multivectors
  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<TMV> x = Teuchos::rcp(new TMV(map, numCols));
  Teuchos::RCP<TMV> y = Teuchos::rcp(new TMV(map, numCols));

  // Perform operation with Tpetra directly
  x->randomize();

  // Do the same thing through NOX/Thyra interface
  Teuchos::RCP<TMVB> y_thyra  = Thyra::createMultiVector(y);
  Teuchos::RCP<NTMV> y_nox = Teuchos::rcp(new NTMV(y_thyra));

  y_nox->random();

  // Check for correct answer
  // Probaby can't expect these to get the same answer
  //success = checkEqualMultiVectors(x, y, out);
}

TEUCHOS_UNIT_TEST(Tpetra_MultiVectorOps, Scale)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;

  // Create Tpetra multivectors
  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<TMV> x = Teuchos::rcp(new TMV(map, numCols));
  Teuchos::RCP<TMV> y = Teuchos::rcp(new TMV(map, numCols));

  x->putScalar(ST::one());
  y->putScalar(ST::one());

  // Perform operation with Tpetra directly
  x->scale(2.0*ST::one());

  // Do the same thing through NOX/Thyra interface
  Teuchos::RCP<TMVB> y_thyra  = Thyra::createMultiVector(y);
  Teuchos::RCP<NTMV> y_nox = Teuchos::rcp(new NTMV(y_thyra));

  y_nox->scale(2.0*ST::one());

  // Check for correct answer
  success = checkEqualMultiVectors(x, y, out);
}

TEUCHOS_UNIT_TEST(Tpetra_MultiVectorOps, Update_1)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;

  // Create Tpetra multivectors
  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<TMV> x = Teuchos::rcp(new TMV(map, numCols));
  Teuchos::RCP<TMV> y = Teuchos::rcp(new TMV(map, numCols));
  Teuchos::RCP<TMV> z = Teuchos::rcp(new TMV(map, numCols));

  x->putScalar(ST::one());
  y->putScalar(ST::one());
  z->putScalar(ST::one());

  // Perform operation with Tpetra directly
  x->update(2.0*ST::one(), *z, ST::one());

  // Do the same thing through NOX/Thyra interface
  Teuchos::RCP<TMVB> y_thyra  = Thyra::createMultiVector(y);
  Teuchos::RCP<NTMV> y_nox = Teuchos::rcp(new NTMV(y_thyra));
  Teuchos::RCP<TMVB> z_thyra  = Thyra::createMultiVector(z);
  Teuchos::RCP<NTMV> z_nox = Teuchos::rcp(new NTMV(z_thyra));


  y_nox->update(2.0*ST::one(), *z_nox, ST::one());

  // Check for correct answer
  success = checkEqualMultiVectors(x, y, out);
}

TEUCHOS_UNIT_TEST(Tpetra_MultiVectorOps, Update_2)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;

  // Create Tpetra multivectors
  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<TMV> w = Teuchos::rcp(new TMV(map, numCols));
  Teuchos::RCP<TMV> x = Teuchos::rcp(new TMV(map, numCols));
  Teuchos::RCP<TMV> y = Teuchos::rcp(new TMV(map, numCols));
  Teuchos::RCP<TMV> z = Teuchos::rcp(new TMV(map, numCols));

  w->putScalar(ST::one());
  x->putScalar(ST::one());
  y->putScalar(ST::one());
  z->putScalar(ST::one());

  // Perform operation with Tpetra directly
  x->update(2.0*ST::one(), *w, 2.0*ST::one(), *z, ST::one());

  // Do the same thing through NOX/Thyra interface
  Teuchos::RCP<TMVB> w_thyra  = Thyra::createMultiVector(w);
  Teuchos::RCP<NTMV> w_nox = Teuchos::rcp(new NTMV(w_thyra));
  Teuchos::RCP<TMVB> y_thyra  = Thyra::createMultiVector(y);
  Teuchos::RCP<NTMV> y_nox = Teuchos::rcp(new NTMV(y_thyra));
  Teuchos::RCP<TMVB> z_thyra  = Thyra::createMultiVector(z);
  Teuchos::RCP<NTMV> z_nox = Teuchos::rcp(new NTMV(z_thyra));


  y_nox->update(2.0*ST::one(), *w_nox, 2.0*ST::one(), *z_nox, ST::one());

  // Check for correct answer
  success = checkEqualMultiVectors(x, y, out);
}

TEUCHOS_UNIT_TEST(Tpetra_MultiVectorOps, Update_3)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;

  // Create Tpetra multivectors
  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<TMV> x = Teuchos::rcp(new TMV(map, numCols));
  Teuchos::RCP<TMV> y = Teuchos::rcp(new TMV(map, numCols));
  Teuchos::RCP<TMV> z = Teuchos::rcp(new TMV(map, 5));

  x->putScalar(ST::one());
  y->putScalar(ST::one());
  z->putScalar(ST::one());

  // Create matrix for the update
  NOX::Abstract::MultiVector::DenseMatrix mat(5, numCols);
  for (int j = 0; j < numCols; ++j) {
    for (int i = 0; i < 5; ++i) {
      mat(i,j) = j;
    }
  }

  // Create locally replicated MultiVector with this data
  // It seems that Teuchos::SerialDenseMatrix is column-major...
  Teuchos::RCP<const Map> map_replicated = Teuchos::rcp(new const Map(5, 0, comm, Tpetra::LocallyReplicated));
  Teuchos::ArrayView<Scalar> mat_view = Teuchos::arrayView(mat.values(), 5*numCols);
  Teuchos::RCP<TMV> mat_mv = Teuchos::rcp(new TMV(map_replicated, mat_view, 5, numCols));

  // Perform operation with Tpetra directly
  x->multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, ST::one(), *z, *mat_mv, ST::one());
  x->describe(out, Teuchos::VERB_EXTREME);

  // Do the same thing through NOX/Thyra interface
  Teuchos::RCP<TMVB> y_thyra  = Thyra::createMultiVector(y);
  Teuchos::RCP<NTMV> y_nox = Teuchos::rcp(new NTMV(y_thyra));
  Teuchos::RCP<TMVB> z_thyra  = Thyra::createMultiVector(z);
  Teuchos::RCP<NTMV> z_nox = Teuchos::rcp(new NTMV(z_thyra));

  y_nox->update(Teuchos::NO_TRANS, ST::one(), *z_nox, mat, ST::one());

  // Check for correct answer
  success = checkEqualMultiVectors(x, y, out);
}

TEUCHOS_UNIT_TEST(Tpetra_MultiVectorOps, Norm_1)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;

  // Create Tpetra multivectors
  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<TMV> x = Teuchos::rcp(new TMV(map, numCols));

  x->putScalar(ST::one());

  // Perform operation with Tpetra directly
  Teuchos::Tuple<ST::magnitudeType, numCols> norms_tpetra;
  x->norm1(norms_tpetra);

  // Do the same thing through NOX/Thyra interface
  Teuchos::RCP<TMVB> x_thyra  = Thyra::createMultiVector(x);
  Teuchos::RCP<NTMV> x_nox = Teuchos::rcp(new NTMV(x_thyra));

  std::vector<ST::magnitudeType> norms_nox(numCols);
  x_nox->norm(norms_nox, NOX::Abstract::Vector::OneNorm);

  // Check for correct answer
  success = checkEqualReductions(norms_tpetra, Teuchos::arrayViewFromVector(norms_nox), out);
}

TEUCHOS_UNIT_TEST(Tpetra_MultiVectorOps, Norm_2)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;

  // Create Tpetra multivectors
  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<TMV> x = Teuchos::rcp(new TMV(map, numCols));

  x->putScalar(ST::one());

  // Perform operation with Tpetra directly
  Teuchos::Tuple<ST::magnitudeType, numCols> norms_tpetra;
  x->norm2(norms_tpetra);

  // Do the same thing through NOX/Thyra interface
  Teuchos::RCP<TMVB> x_thyra  = Thyra::createMultiVector(x);
  Teuchos::RCP<NTMV> x_nox = Teuchos::rcp(new NTMV(x_thyra));

  std::vector<ST::magnitudeType> norms_nox(numCols);
  x_nox->norm(norms_nox, NOX::Abstract::Vector::TwoNorm);

  // Check for correct answer
  success = checkEqualReductions(norms_tpetra, Teuchos::arrayViewFromVector(norms_nox), out);
}

TEUCHOS_UNIT_TEST(Tpetra_MultiVectorOps, Norm_Inf)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;

  // Create Tpetra multivectors
  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<TMV> x = Teuchos::rcp(new TMV(map, numCols));

  x->putScalar(ST::one());

  // Perform operation with Tpetra directly
  Teuchos::Tuple<ST::magnitudeType, numCols> norms_tpetra;
  x->normInf(norms_tpetra);

  // Do the same thing through NOX/Thyra interface
  Teuchos::RCP<TMVB> x_thyra  = Thyra::createMultiVector(x);
  Teuchos::RCP<NTMV> x_nox = Teuchos::rcp(new NTMV(x_thyra));

  std::vector<ST::magnitudeType> norms_nox(numCols);
  x_nox->norm(norms_nox, NOX::Abstract::Vector::MaxNorm);

  // Check for correct answer
  success = checkEqualReductions(norms_tpetra, Teuchos::arrayViewFromVector(norms_nox), out);
}

TEUCHOS_UNIT_TEST(Tpetra_MultiVectorOps, Multiply)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;

  // Create Tpetra multivectors
  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<TMV> x = Teuchos::rcp(new TMV(map, numCols));
  Teuchos::RCP<TMV> y = Teuchos::rcp(new TMV(map, numCols));

  for (int i = 0; i < numCols; ++i) {
    x->getVectorNonConst(i)->putScalar(i*ST::one());
    y->getVectorNonConst(i)->putScalar(2.0*i*ST::one());
  }

  // Perform operation with Tpetra directly
  NOX::Abstract::MultiVector::DenseMatrix dots_tpetra(numCols, numCols);
  for (int j = 0; j < numCols; ++j) {
    for (int i = 0; i < numCols; ++i) {
      dots_tpetra(i, j) = x->getVector(i)->dot(*(y->getVector(j)));
    }
  }

  // Do the same thing through NOX/Thyra interface
  Teuchos::RCP<TMVB> x_thyra  = Thyra::createMultiVector(x);
  Teuchos::RCP<NTMV> x_nox = Teuchos::rcp(new NTMV(x_thyra));
  Teuchos::RCP<TMVB> y_thyra  = Thyra::createMultiVector(y);
  Teuchos::RCP<NTMV> y_nox = Teuchos::rcp(new NTMV(y_thyra));

  NOX::Abstract::MultiVector::DenseMatrix dots_nox(numCols, numCols);
  x_nox->multiply(ST::one(), *y_nox, dots_nox); 

  // Check for correct answer
  success = checkEqualReductions(Teuchos::arrayView(dots_tpetra.values(), numCols*numCols),
    Teuchos::arrayView(dots_nox.values(), numCols*numCols), out);
}

