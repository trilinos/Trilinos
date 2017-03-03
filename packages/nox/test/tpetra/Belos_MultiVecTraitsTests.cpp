#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_RCP.hpp"

#include "Tpetra_Core.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"
#include "BelosThyraAdapter.hpp"

// Typedefs and constants

typedef Tpetra::MultiVector<>::scalar_type Scalar;
typedef Tpetra::MultiVector<>::local_ordinal_type LO;
typedef Tpetra::MultiVector<>::global_ordinal_type GO;
typedef Tpetra::MultiVector<>::node_type Node;

typedef Teuchos::ScalarTraits<Scalar> ST;
typedef Tpetra::Map<LO,GO,Node> Map;
typedef Tpetra::Vector<Scalar,LO,GO,Node> TV;
typedef Tpetra::MultiVector<Scalar,LO,GO,Node> TMV;
typedef Thyra::VectorSpaceBase<Scalar> TVSB;
typedef Thyra::VectorBase<Scalar> TVB;
typedef Thyra::MultiVectorBase<Scalar> TMVB;
typedef Thyra::TpetraOperatorVectorExtraction<Scalar,LO,GO,Node> TOVE;
typedef Belos::MultiVecTraits<Scalar, TMVB> BMVT;
typedef typename TV::mag_type mag_type;

#if defined HAVE_TPETRACORE_CUDA
#define NUM_LOCAL 1000000
#else
#define NUM_LOCAL 1000
#endif
const std::size_t numLocalElements = NUM_LOCAL;

//Routines for checking solution

void checkMultiVectors(const Teuchos::RCP<TMV>& a,
                       const Teuchos::RCP<TMV>& b,
                       const Teuchos::ArrayView<mag_type>& expectedNorms,
                       Teuchos::FancyOStream& out,
                       bool& success)
{
  TEUCHOS_ASSERT(a->getNumVectors() == b->getNumVectors());
  TEUCHOS_ASSERT(a->getNumVectors() == expectedNorms.size());
  Scalar one = ST::one();
  b->update(-1.0*one, *a, one);
  std::vector<ST::magnitudeType> norms(expectedNorms.size());
  b->norm2(Teuchos::arrayViewFromVector(norms));

  ST::magnitudeType tol = 1.0e-14;
  for (auto norm = norms.begin(); norm != norms.end(); ++norm) {
    TEUCHOS_TEST_EQUALITY(*norm < tol, true, out, success);
  }
  a->norm2(Teuchos::arrayViewFromVector(norms));
  auto normIter = norms.begin();
  auto expNormIter = expectedNorms.begin();
  for (; normIter != norms.end(); ++normIter, ++expNormIter) {
    TEUCHOS_TEST_FLOATING_EQUALITY(*normIter, *expNormIter, tol, out, success);
  }
}

template <class T>
void checkReductions(const Teuchos::ArrayView<T>& a,
                     const Teuchos::ArrayView<T>& b,
                     const Teuchos::ArrayView<T>& expectedVals,
                     Teuchos::FancyOStream& out,
                     bool& success)
{
  TEUCHOS_ASSERT(a.size() == b.size());
  TEUCHOS_ASSERT(a.size() == expectedVals.size());

  typename Teuchos::ScalarTraits<T>::magnitudeType tol = 1.0e-14;
  auto aItr = a.begin();
  auto bItr = b.begin();
  auto expItr = expectedVals.begin();
  for (; aItr != a.end(); ++aItr, ++bItr, ++expItr) {
    TEUCHOS_TEST_FLOATING_EQUALITY(*aItr, *bItr, tol, out, success);
    TEUCHOS_TEST_FLOATING_EQUALITY(*aItr, *expItr, tol, out, success);
  }
}

// Unit Tests

TEUCHOS_UNIT_TEST(Belos_MultiVectorTraits, ThyraTpetraVector_Basic)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;
  const Scalar one = ST::one();

  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<const TVSB> space = Thyra::createVectorSpace<Scalar,LO,GO,Node>(map);

  Teuchos::RCP<TVB> x = Thyra::createMember(space);
  Thyra::assign(x.ptr(), one);
  Teuchos::RCP<TV> x_tpetra = TOVE::getTpetraVector(x);

  // Check attributes
  TEUCHOS_TEST_EQUALITY(BMVT::GetGlobalLength(*x), x_tpetra->getGlobalLength(), out, success);
  TEUCHOS_TEST_EQUALITY(BMVT::GetNumberVecs(*x), x_tpetra->getNumVectors(), out, success);
  // Should HasConstantStride exist?
  //TEUCHOS_TEST_EQUALITY(BMVT::HasConstantStride(*x), x_tpetra->isConstantStride(), out, success);
}

TEUCHOS_UNIT_TEST(Belos_MultiVectorTraits, ThyraTpetraVector_CloneCopy)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;
  const Scalar one = ST::one();

  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<const TVSB> space = Thyra::createVectorSpace<Scalar,LO,GO,Node>(map);

  Teuchos::RCP<TVB> x = Thyra::createMember(space);
  Thyra::assign(x.ptr(), one);
  Teuchos::RCP<TV> x_tpetra = TOVE::getTpetraVector(x);

  Teuchos::RCP<TMVB> y = BMVT::CloneCopy(*x);
  Teuchos::RCP<TMV> y_tpetra = TOVE::getTpetraMultiVector(y);

  mag_type val = static_cast<mag_type>(ST::squareroot(numGlobalElements));
  auto ans = Teuchos::tuple(val);
  checkMultiVectors(x_tpetra, y_tpetra, ans, out, success);
}

TEUCHOS_UNIT_TEST(Belos_MultiVectorTraits, ThyraTpetraVector_CloneViewNonConst)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;
  const std::size_t numCols = 1;
  const Scalar one = ST::one();

  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<const TVSB> space = Thyra::createVectorSpace<Scalar,LO,GO,Node>(map);
  Teuchos::RCP<TVB> x = Thyra::createMember(space);
  Thyra::assign(x.ptr(), one);
  Teuchos::RCP<TV> x_tpetra = TOVE::getTpetraVector(x);

  Teuchos::RCP<TMVB> y = BMVT::CloneViewNonConst(*x, Teuchos::Range1D());
  Teuchos::RCP<TMV> y_tpetra = TOVE::getTpetraMultiVector(y);
  // y is view of x, so changing y should change x as well
  y_tpetra->scale(static_cast<Scalar>(2.0));

  mag_type val = static_cast<mag_type>(2.0*ST::squareroot(numGlobalElements));
  auto ans = Teuchos::tuple(val);
  checkMultiVectors(x_tpetra, y_tpetra->subCopy(Teuchos::Range1D(0,numCols-1)), ans, out, success);
}


TEUCHOS_UNIT_TEST(Belos_MultiVectorTraits, ThyraTpetraVector_MvTimesMatAddMv)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;
  const std::size_t numCols = 1;
  const Scalar one = ST::one();

  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<const TVSB> space = Thyra::createVectorSpace<Scalar,LO,GO,Node>(map);

  Teuchos::RCP<TVB> x = Thyra::createMember(space);
  Teuchos::RCP<TVB> y = Thyra::createMember(space);
  Teuchos::RCP<TVB> z = Thyra::createMember(space);
  Thyra::assign(x.ptr(), one);
  Thyra::assign(y.ptr(), one);
  Thyra::assign(z.ptr(), one);
  Teuchos::RCP<TV> x_tpetra = TOVE::getTpetraVector(x);
  Teuchos::RCP<TV> y_tpetra = TOVE::getTpetraVector(y);
  Teuchos::RCP<TV> z_tpetra = TOVE::getTpetraVector(z);

  Teuchos::SerialDenseMatrix<int, Scalar> mat(numCols, numCols, true);
  for (int i = 0; i < mat.numRows(); ++i) {
      mat(i,i) = static_cast<Scalar>(i+2);
  }
  BMVT::MvTimesMatAddMv(one, *z, mat, one, *x);
  for (std::size_t i = 0; i < numCols; ++i) {
    y_tpetra->getVectorNonConst(i)->update(mat(i,i), *z_tpetra->getVector(i), one);
  }

  mag_type val = static_cast<mag_type>(3.0*ST::squareroot(numGlobalElements));
  auto ans = Teuchos::tuple(val);
  checkMultiVectors(x_tpetra, y_tpetra, ans, out, success);
}

TEUCHOS_UNIT_TEST(Belos_MultiVectorTraits, ThyraTpetraVector_MvAddMv)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;
  const Scalar one = ST::one();

  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<const TVSB> space = Thyra::createVectorSpace<Scalar,LO,GO,Node>(map);

  Teuchos::RCP<TVB> x = Thyra::createMember(space);
  Teuchos::RCP<TVB> y = Thyra::createMember(space);
  Teuchos::RCP<TVB> z = Thyra::createMember(space);
  Thyra::assign(x.ptr(), one);
  Thyra::assign(y.ptr(), one);
  Thyra::assign(z.ptr(), one);
  Teuchos::RCP<TV> x_tpetra = TOVE::getTpetraVector(x);
  Teuchos::RCP<TV> y_tpetra = TOVE::getTpetraVector(y);
  Teuchos::RCP<TV> z_tpetra = TOVE::getTpetraVector(z);

  BMVT::MvAddMv(2.0*one, *z, 2.0*one, *x, *x);
  y_tpetra->update(2.0*one, *z_tpetra, 2.0*one);

  mag_type val = static_cast<mag_type>(4.0*ST::squareroot(numGlobalElements));
  auto ans = Teuchos::tuple(val);
  checkMultiVectors(x_tpetra, y_tpetra, ans, out, success);
}

TEUCHOS_UNIT_TEST(Belos_MultiVectorTraits, ThyraTpetraVector_MvScaleScalar)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;
  const Scalar one = ST::one();

  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<const TVSB> space = Thyra::createVectorSpace<Scalar,LO,GO,Node>(map);

  Teuchos::RCP<TVB> x = Thyra::createMember(space);
  Teuchos::RCP<TVB> y = Thyra::createMember(space);
  Thyra::assign(x.ptr(), one);
  Thyra::assign(y.ptr(), one);
  Teuchos::RCP<TV> x_tpetra = TOVE::getTpetraVector(x);
  Teuchos::RCP<TV> y_tpetra = TOVE::getTpetraVector(y);

  BMVT::MvScale(*x, 2.0*one);
  y_tpetra->scale(2.0*one);

  mag_type val = static_cast<mag_type>(2.0*ST::squareroot(numGlobalElements));
  auto ans = Teuchos::tuple(val);
  checkMultiVectors(x_tpetra, y_tpetra, ans, out, success);
}

TEUCHOS_UNIT_TEST(Belos_MultiVectorTraits, ThyraTpetraVector_MvScaleVector)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;
  const std::size_t numCols = 1;
  const Scalar one = ST::one();

  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<const TVSB> space = Thyra::createVectorSpace<Scalar,LO,GO,Node>(map);

  Teuchos::RCP<TVB> x = Thyra::createMember(space);
  Teuchos::RCP<TVB> y = Thyra::createMember(space);
  Thyra::assign(x.ptr(), one);
  Thyra::assign(y.ptr(), one);
  Teuchos::RCP<TV> x_tpetra = TOVE::getTpetraVector(x);
  Teuchos::RCP<TV> y_tpetra = TOVE::getTpetraVector(y);

  std::vector<Scalar> alpha(numCols);
  for (std::vector<Scalar>::size_type i = 0; i < alpha.size(); ++i)
    alpha[i] = static_cast<Scalar>(i+2);
  BMVT::MvScale(*x, alpha);
  y_tpetra->scale(Teuchos::arrayViewFromVector(alpha));

  mag_type val = static_cast<mag_type>(2.0*ST::squareroot(numGlobalElements));
  auto ans = Teuchos::tuple(val);
  checkMultiVectors(x_tpetra, y_tpetra, ans, out, success);
}

TEUCHOS_UNIT_TEST(Belos_MultiVectorTraits, ThyraTpetraVector_MvTransMv)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;
  const std::size_t numCols = 1;
  const Scalar one = ST::one();

  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<const TVSB> space = Thyra::createVectorSpace<Scalar,LO,GO,Node>(map);

  Teuchos::RCP<TVB> x = Thyra::createMember(space);
  Teuchos::RCP<TVB> y = Thyra::createMember(space);
  Thyra::assign(x.ptr(), one);
  Thyra::assign(y.ptr(), one);
  Teuchos::RCP<TV> x_tpetra = TOVE::getTpetraVector(x);
  Teuchos::RCP<TV> y_tpetra = TOVE::getTpetraVector(y);

  Teuchos::SerialDenseMatrix<int, Scalar> mat_belos(numCols, numCols);
  Teuchos::SerialDenseMatrix<int, Scalar> mat_tpetra(numCols, numCols);
  BMVT::MvTransMv(one, *y, *x, mat_belos);
  for (int i = 0; i < mat_tpetra.numCols(); ++ i) {
    for (int j = 0; j < mat_tpetra.numRows(); ++j) {
      mat_tpetra(i,j) = x_tpetra->getVector(i)->dot(*y_tpetra->getVector(j));
    }
  }

  Scalar val = static_cast<Scalar>(numGlobalElements);
  Teuchos::Array<Scalar> ans(numCols*numCols, val);
  checkReductions(Teuchos::arrayView(mat_belos.values(), numCols*numCols),
    Teuchos::arrayView(mat_tpetra.values(), numCols*numCols), ans(), out, success);
}

TEUCHOS_UNIT_TEST(Belos_MultiVectorTraits, ThyraTpetraVector_MvDot)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;
  const std::size_t numCols = 1;
  const Scalar one = ST::one();

  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<const TVSB> space = Thyra::createVectorSpace<Scalar,LO,GO,Node>(map);

  Teuchos::RCP<TVB> x = Thyra::createMember(space);
  Teuchos::RCP<TVB> y = Thyra::createMember(space);
  Thyra::assign(x.ptr(), one);
  Thyra::assign(y.ptr(), 5.0*one);
  Teuchos::RCP<TV> x_tpetra = TOVE::getTpetraVector(x);
  Teuchos::RCP<TV> y_tpetra = TOVE::getTpetraVector(y);

  std::vector<Scalar> dot_belos(numCols);
  Teuchos::Tuple<Scalar, numCols> dot_tpetra;
  BMVT::MvDot(*x, *y, dot_belos);
  x_tpetra->dot(*y_tpetra, dot_tpetra);

  Scalar val = static_cast<Scalar>(5.0*numGlobalElements);
  auto ans = Teuchos::tuple(val);
  checkReductions(Teuchos::arrayViewFromVector(dot_belos), dot_tpetra, ans, out, success);
}

TEUCHOS_UNIT_TEST(Belos_MultiVectorTraits, ThyraTpetraVector_MvNorm)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;
  const std::size_t numCols = 1;
  const Scalar one = ST::one();

  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<const TVSB> space = Thyra::createVectorSpace<Scalar,LO,GO,Node>(map);

  Teuchos::RCP<TVB> x = Thyra::createMember(space);
  Thyra::assign(x.ptr(), one);
  Teuchos::RCP<TV> x_tpetra = TOVE::getTpetraVector(x);

  std::vector<Scalar> norm_belos(numCols);
  Teuchos::Tuple<Scalar, numCols> norm_tpetra;
  BMVT::MvNorm(*x, norm_belos, Belos::TwoNorm);
  x_tpetra->norm2(norm_tpetra);

  mag_type val = static_cast<mag_type>(ST::squareroot(numGlobalElements));
  auto ans = Teuchos::tuple(val);
  checkReductions(Teuchos::arrayViewFromVector(norm_belos), norm_tpetra, ans, out, success);
}

TEUCHOS_UNIT_TEST(Belos_MultiVectorTraits, ThyraTpetraVector_Assign)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;
  const Scalar one = ST::one();

  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<const TVSB> space = Thyra::createVectorSpace<Scalar,LO,GO,Node>(map);

  Teuchos::RCP<TVB> x = Thyra::createMember(space);
  Teuchos::RCP<TVB> y = Thyra::createMember(space);
  Thyra::assign(x.ptr(), one);
  Thyra::assign(y.ptr(), 5.0*one);
  Teuchos::RCP<TV> x_tpetra = TOVE::getTpetraVector(x);
  Teuchos::RCP<TV> y_tpetra = TOVE::getTpetraVector(y);

  BMVT::Assign(*y, *x);

  mag_type val = static_cast<mag_type>(5.0*ST::squareroot(numGlobalElements));
  auto ans = Teuchos::tuple(val);
  checkMultiVectors(x_tpetra, y_tpetra, ans, out, success);
}

TEUCHOS_UNIT_TEST(Belos_MultiVectorTraits, ThyraTpetraVector_MvRandom)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;
  const  Scalar one = ST::one();

  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<const TVSB> space = Thyra::createVectorSpace<Scalar,LO,GO,Node>(map);

  Teuchos::RCP<TVB> x = Thyra::createMember(space);
  Teuchos::RCP<TVB> y = Thyra::createMember(space);
  //Teuchos::RCP<TV> x_tpetra = TOVE::getTpetraVector(x);
  Teuchos::RCP<TV> y_tpetra = TOVE::getTpetraVector(y);

  BMVT::MvRandom(*x);
  y_tpetra->randomize(-one, one);

  //checkMultiVectors(x_tpetra, y_tpetra, ans, out, success);
}

TEUCHOS_UNIT_TEST(Belos_MultiVectorTraits, ThyraTpetraVector_MvInit)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;
  const Scalar one = ST::one();

  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<const TVSB> space = Thyra::createVectorSpace<Scalar,LO,GO,Node>(map);

  Teuchos::RCP<TVB> x = Thyra::createMember(space);
  Teuchos::RCP<TVB> y = Thyra::createMember(space);
  Teuchos::RCP<TV> x_tpetra = TOVE::getTpetraVector(x);
  Teuchos::RCP<TV> y_tpetra = TOVE::getTpetraVector(y);

  BMVT::MvInit(*x, 2.0*one);
  y_tpetra->putScalar(2.0*one);

  mag_type val = static_cast<mag_type>(2.0*ST::squareroot(numGlobalElements));
  auto ans = Teuchos::tuple(val);
  checkMultiVectors(x_tpetra, y_tpetra, ans, out, success);
}

TEUCHOS_UNIT_TEST(Belos_MultiVectorTraits, ThyraTpetraMultiVector_Basic)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;
  const std::size_t numCols = 5;
  const Scalar one = ST::one();

  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<const TVSB> space = Thyra::createVectorSpace<Scalar,LO,GO,Node>(map);

  Teuchos::RCP<TMVB> x = Thyra::createMembers(space, numCols);
  Thyra::assign(x.ptr(), one);
  Teuchos::RCP<TMV> x_tpetra = TOVE::getTpetraMultiVector(x);

  // Check attributes
  TEUCHOS_TEST_EQUALITY(BMVT::GetGlobalLength(*x), x_tpetra->getGlobalLength(), out, success);
  TEUCHOS_TEST_EQUALITY(BMVT::GetNumberVecs(*x), x_tpetra->getNumVectors(), out, success);
  // Should HasConstantStride exist?
  //TEUCHOS_TEST_EQUALITY(BMVT::HasConstantStride(*x), x_tpetra->isConstantStride(), out, success);
}

TEUCHOS_UNIT_TEST(Belos_MultiVectorTraits, ThyraTpetraMultiVector_CloneCopy)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;
  const std::size_t numCols = 5;
  const Scalar one = ST::one();

  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<const TVSB> space = Thyra::createVectorSpace<Scalar,LO,GO,Node>(map);

  Teuchos::RCP<TMVB> x = Thyra::createMembers(space, numCols);
  Thyra::assign(x.ptr(), one);
  Teuchos::RCP<TMV> x_tpetra = TOVE::getTpetraMultiVector(x);

  Teuchos::RCP<TMVB> y = BMVT::CloneCopy(*x);
  Teuchos::RCP<TMV> y_tpetra = TOVE::getTpetraMultiVector(y);

  mag_type val = static_cast<mag_type>(ST::squareroot(numGlobalElements));
  auto ans = Teuchos::tuple(val, val, val, val, val);
  checkMultiVectors(x_tpetra, y_tpetra, ans, out, success);
}

TEUCHOS_UNIT_TEST(Belos_MultiVectorTraits, ThyraTpetraMultiVector_CloneViewNonConst)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;
  const std::size_t numCols = 5;
  const Scalar one = ST::one();

  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<const TVSB> space = Thyra::createVectorSpace<Scalar,LO,GO,Node>(map);

  Teuchos::RCP<TMVB> x = Thyra::createMembers(space, numCols);
  Thyra::assign(x.ptr(), one);
  Teuchos::RCP<TMV> x_tpetra = TOVE::getTpetraMultiVector(x);

  Teuchos::RCP<TMVB> y = BMVT::CloneViewNonConst(*x, Teuchos::Range1D());
  Teuchos::RCP<TMV> y_tpetra = TOVE::getTpetraMultiVector(y);
  // y is view of x, so changing y should change x as well
  y_tpetra->scale(static_cast<Scalar>(2.0));

  mag_type val = static_cast<mag_type>(2.0*ST::squareroot(numGlobalElements));
  auto ans = Teuchos::tuple(val, val, val, val, val);
    checkMultiVectors(x_tpetra, x_tpetra->subCopy(Teuchos::Range1D(0,numCols-1)), ans, out, success);
}

TEUCHOS_UNIT_TEST(Belos_MultiVectorTraits, ThyraTpetraMultiVector_MvTimesMatAddMv)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;
  const std::size_t numCols = 5;
  const Scalar one = ST::one();

  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<const TVSB> space = Thyra::createVectorSpace<Scalar,LO,GO,Node>(map);

  Teuchos::RCP<TMVB> x = Thyra::createMembers(space, numCols);
  Teuchos::RCP<TMVB> y = Thyra::createMembers(space, numCols);
  Teuchos::RCP<TMVB> z = Thyra::createMembers(space, numCols);
  Thyra::assign(x.ptr(), one);
  Thyra::assign(y.ptr(), one);
  Thyra::assign(z.ptr(), one);
  Teuchos::RCP<TMV> x_tpetra = TOVE::getTpetraMultiVector(x);
  Teuchos::RCP<TMV> y_tpetra = TOVE::getTpetraMultiVector(y);
  Teuchos::RCP<TMV> z_tpetra = TOVE::getTpetraMultiVector(z);

  Thyra::assign(x.ptr(), one);
  Teuchos::SerialDenseMatrix<int, Scalar> mat(numCols, numCols, true);
  for (int i = 0; i < mat.numRows(); ++i) {
      mat(i,i) = static_cast<Scalar>(i+2);
  }
  BMVT::MvTimesMatAddMv(one, *z, mat, one, *x);
  for (std::size_t i = 0; i < numCols; ++i) {
    y_tpetra->getVectorNonConst(i)->update(mat(i,i), *z_tpetra->getVector(i), one);
  }

  mag_type val = static_cast<mag_type>(ST::squareroot(numGlobalElements));
  auto ans = Teuchos::tuple(3.0*val, 4.0*val, 5.0*val, 6.0*val, 7.0*val);
  checkMultiVectors(x_tpetra, y_tpetra, ans, out, success);
}

TEUCHOS_UNIT_TEST(Belos_MultiVectorTraits, ThyraTpetraMultiVector_MvAddMv)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;
  const std::size_t numCols = 5;
  const Scalar one = ST::one();

  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<const TVSB> space = Thyra::createVectorSpace<Scalar,LO,GO,Node>(map);

  Teuchos::RCP<TMVB> x = Thyra::createMembers(space, numCols);
  Teuchos::RCP<TMVB> y = Thyra::createMembers(space, numCols);
  Teuchos::RCP<TMVB> z = Thyra::createMembers(space, numCols);
  Thyra::assign(x.ptr(), one);
  Thyra::assign(y.ptr(), one);
  Thyra::assign(z.ptr(), one);
  Teuchos::RCP<TMV> x_tpetra = TOVE::getTpetraMultiVector(x);
  Teuchos::RCP<TMV> y_tpetra = TOVE::getTpetraMultiVector(y);
  Teuchos::RCP<TMV> z_tpetra = TOVE::getTpetraMultiVector(z);

  BMVT::MvAddMv(2.0*one, *z, 2.0*one, *x, *x);
  y_tpetra->update(2.0*one, *z_tpetra, 2.0*one);

  mag_type val = static_cast<mag_type>(4.0*ST::squareroot(numGlobalElements));
  auto ans = Teuchos::tuple(val, val, val, val, val);
  checkMultiVectors(x_tpetra, y_tpetra, ans, out, success);
}

TEUCHOS_UNIT_TEST(Belos_MultiVectorTraits, ThyraTpetraMultiVector_MvScaleScalar)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;
  const std::size_t numCols = 5;
  const Scalar one = ST::one();

  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<const TVSB> space = Thyra::createVectorSpace<Scalar,LO,GO,Node>(map);

  Teuchos::RCP<TMVB> x = Thyra::createMembers(space, numCols);
  Teuchos::RCP<TMVB> y = Thyra::createMembers(space, numCols);
  Thyra::assign(x.ptr(), one);
  Thyra::assign(y.ptr(), one);
  Teuchos::RCP<TMV> x_tpetra = TOVE::getTpetraMultiVector(x);
  Teuchos::RCP<TMV> y_tpetra = TOVE::getTpetraMultiVector(y);

  BMVT::MvScale(*x, 2.0*one);
  y_tpetra->scale(2.0*one);

  mag_type val = static_cast<mag_type>(2.0*ST::squareroot(numGlobalElements));
  auto ans = Teuchos::tuple(val, val, val, val, val);
  checkMultiVectors(x_tpetra, y_tpetra, ans, out, success);
}

TEUCHOS_UNIT_TEST(Belos_MultiVectorTraits, ThyraTpetraMultiVector_MvScaleVector)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;
  const std::size_t numCols = 5;
  const Scalar one = ST::one();

  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<const TVSB> space = Thyra::createVectorSpace<Scalar,LO,GO,Node>(map);

  Teuchos::RCP<TMVB> x = Thyra::createMembers(space, numCols);
  Teuchos::RCP<TMVB> y = Thyra::createMembers(space, numCols);
  Thyra::assign(x.ptr(), one);
  Thyra::assign(y.ptr(), one);
  Teuchos::RCP<TMV> x_tpetra = TOVE::getTpetraMultiVector(x);
  Teuchos::RCP<TMV> y_tpetra = TOVE::getTpetraMultiVector(y);

  std::vector<Scalar> alpha(numCols);
  for (std::vector<Scalar>::size_type i = 0; i < alpha.size(); ++i)
    alpha[i] = static_cast<Scalar>(i+2);
  BMVT::MvScale(*x, alpha);
  y_tpetra->scale(Teuchos::arrayViewFromVector(alpha));

  mag_type val = static_cast<mag_type>(ST::squareroot(numGlobalElements));
  auto ans = Teuchos::tuple(2.0*val, 3.0*val, 4.0*val, 5.0*val, 6.0*val);
  checkMultiVectors(x_tpetra, y_tpetra, ans, out, success);
}

TEUCHOS_UNIT_TEST(Belos_MultiVectorTraits, ThyraTpetraMultiVector_MvTransMv)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;
  const std::size_t numCols = 5;
  const Scalar one = ST::one();

  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<const TVSB> space = Thyra::createVectorSpace<Scalar,LO,GO,Node>(map);

  Teuchos::RCP<TMVB> x = Thyra::createMembers(space, numCols);
  Teuchos::RCP<TMVB> y = Thyra::createMembers(space, numCols);
  Thyra::assign(x.ptr(), one);
  Thyra::assign(y.ptr(), one);
  Teuchos::RCP<TMV> x_tpetra = TOVE::getTpetraMultiVector(x);
  Teuchos::RCP<TMV> y_tpetra = TOVE::getTpetraMultiVector(y);

  Teuchos::SerialDenseMatrix<int, Scalar> mat_belos(numCols, numCols);
  Teuchos::SerialDenseMatrix<int, Scalar> mat_tpetra(numCols, numCols);
  BMVT::MvTransMv(one, *y, *x, mat_belos);
  for (int i = 0; i < mat_tpetra.numCols(); ++ i) {
    for (int j = 0; j < mat_tpetra.numRows(); ++j) {
      mat_tpetra(i,j) = x_tpetra->getVector(i)->dot(*y_tpetra->getVector(j));
    }
  }

  Scalar val = static_cast<Scalar>(numGlobalElements);
  Teuchos::Array<Scalar> ans(numCols*numCols, val);
  checkReductions(Teuchos::arrayView(mat_belos.values(), numCols*numCols),
    Teuchos::arrayView(mat_tpetra.values(), numCols*numCols), ans(), out, success);
}

TEUCHOS_UNIT_TEST(Belos_MultiVectorTraits, ThyraTpetraMultiVector_MvDot)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;
  const std::size_t numCols = 5;
  const Scalar one = ST::one();

  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<const TVSB> space = Thyra::createVectorSpace<Scalar,LO,GO,Node>(map);

  Teuchos::RCP<TMVB> x = Thyra::createMembers(space, numCols);
  Teuchos::RCP<TMVB> y = Thyra::createMembers(space, numCols);
  Thyra::assign(x.ptr(), one);
  Thyra::assign(y.ptr(), 5.0*one);
  Teuchos::RCP<TMV> x_tpetra = TOVE::getTpetraMultiVector(x);
  Teuchos::RCP<TMV> y_tpetra = TOVE::getTpetraMultiVector(y);

  std::vector<Scalar> dot_belos(numCols);
  Teuchos::Tuple<Scalar, numCols> dot_tpetra;
  BMVT::MvDot(*x, *y, dot_belos);
  x_tpetra->dot(*y_tpetra, dot_tpetra);

  Scalar val = static_cast<Scalar>(5.0*numGlobalElements);
  auto ans = Teuchos::tuple(val, val, val, val, val);
  checkReductions(Teuchos::arrayViewFromVector(dot_belos), dot_tpetra, ans, out, success);
}

TEUCHOS_UNIT_TEST(Belos_MultiVectorTraits, ThyraTpetraMultiVector_MvNorm)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;
  const std::size_t numCols = 5;
  const Scalar one = ST::one();

  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<const TVSB> space = Thyra::createVectorSpace<Scalar,LO,GO,Node>(map);

  Teuchos::RCP<TMVB> x = Thyra::createMembers(space, numCols);
  Thyra::assign(x.ptr(), one);
  Teuchos::RCP<TMV> x_tpetra = TOVE::getTpetraMultiVector(x);

  std::vector<Scalar> norm_belos(numCols);
  Teuchos::Tuple<Scalar, numCols> norm_tpetra;
  BMVT::MvNorm(*x, norm_belos, Belos::TwoNorm);
  x_tpetra->norm2(norm_tpetra);

  mag_type val = static_cast<mag_type>(ST::squareroot(numGlobalElements));
  auto ans = Teuchos::tuple(val, val, val, val, val);
  checkReductions(Teuchos::arrayViewFromVector(norm_belos), norm_tpetra, ans, out, success);
}

TEUCHOS_UNIT_TEST(Belos_MultiVectorTraits, ThyraTpetraMultiVector_Assign)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;
  const std::size_t numCols = 5;
  const Scalar one = ST::one();

  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<const TVSB> space = Thyra::createVectorSpace<Scalar,LO,GO,Node>(map);

  Teuchos::RCP<TMVB> x = Thyra::createMembers(space, numCols);
  Teuchos::RCP<TMVB> y = Thyra::createMembers(space, numCols);
  Thyra::assign(x.ptr(), one);
  Thyra::assign(y.ptr(), 5.0*one);
  Teuchos::RCP<TMV> x_tpetra = TOVE::getTpetraMultiVector(x);
  Teuchos::RCP<TMV> y_tpetra = TOVE::getTpetraMultiVector(y);

  BMVT::Assign(*y, *x);

  mag_type val = static_cast<mag_type>(5.0*ST::squareroot(numGlobalElements));
  auto ans = Teuchos::tuple(val, val, val, val, val);
  checkMultiVectors(x_tpetra, y_tpetra, ans, out, success);
}

TEUCHOS_UNIT_TEST(Belos_MultiVectorTraits, ThyraTpetraMultiVector_MvRandom)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;
  const std::size_t numCols = 5;
  const Scalar one = ST::one();

  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<const TVSB> space = Thyra::createVectorSpace<Scalar,LO,GO,Node>(map);

  Teuchos::RCP<TMVB> x = Thyra::createMembers(space, numCols);
  Teuchos::RCP<TMVB> y = Thyra::createMembers(space, numCols);
  Thyra::assign(x.ptr(), one);
  Thyra::assign(y.ptr(), one);
  Teuchos::RCP<TMV> x_tpetra = TOVE::getTpetraMultiVector(x);
  Teuchos::RCP<TMV> y_tpetra = TOVE::getTpetraMultiVector(y);

  BMVT::MvRandom(*x);
  y_tpetra->randomize(-one, one);

  //checkMultiVectors(x_tpetra, y_tpetra, ans, out, success);
}

TEUCHOS_UNIT_TEST(Belos_MultiVectorTraits, ThyraTpetraMultiVector_MvInit)
{
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  const Tpetra::global_size_t numGlobalElements = comm->getSize()*numLocalElements;
  const std::size_t numCols = 5;
  const Scalar one = ST::one();

  Teuchos::RCP<const Map> map = Teuchos::rcp(new const Map(numGlobalElements, numLocalElements, 0, comm));
  Teuchos::RCP<const TVSB> space = Thyra::createVectorSpace<Scalar,LO,GO,Node>(map);

  Teuchos::RCP<TMVB> x = Thyra::createMembers(space, numCols);
  Teuchos::RCP<TMVB> y = Thyra::createMembers(space, numCols);
  Thyra::assign(x.ptr(), one);
  Thyra::assign(y.ptr(), one);
  Teuchos::RCP<TMV> x_tpetra = TOVE::getTpetraMultiVector(x);
  Teuchos::RCP<TMV> y_tpetra = TOVE::getTpetraMultiVector(y);

  BMVT::MvInit(*x, 2.0*one);
  y_tpetra->putScalar(2.0*one);

  mag_type val = static_cast<mag_type>(2.0*ST::squareroot(numGlobalElements));
  auto ans = Teuchos::tuple(val, val, val, val, val);
  checkMultiVectors(x_tpetra, y_tpetra, ans, out, success);
}

