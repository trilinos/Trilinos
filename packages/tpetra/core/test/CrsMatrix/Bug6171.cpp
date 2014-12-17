// Some Macro Magic to ensure that if CUDA and KokkosCompat is enabled
// only the .cu version of this file is actually compiled
#include <Tpetra_config.h>


#include <iostream>
#include <vector>

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Map.hpp>
#include <Teuchos_GlobalMPISession.hpp>

int main (int argc, char *argv[])
{
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;

  typedef int LO; // Local Ordinal
  typedef int GO; // Global Ordinal
  typedef double ST; // data type
  typedef Teuchos::ScalarTraits<ST>::magnitudeType MT;
  typedef Tpetra::global_size_t GST; // Map's constructor needs this

  typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType  NT;
  typedef Tpetra::Map<LO,GO,NT>                                   map_type;
  typedef Tpetra::CrsMatrix<ST,LO,GO,NT>                          crs_matrix_type;
  typedef Tpetra::Vector<ST,LO,GO,NT>                             vec_type;

  Teuchos::GlobalMPISession mpiSession (&argc, &argv, NULL);
  RCP<const Teuchos::Comm<int> > comm =
    Tpetra::DefaultPlatform::getDefaultPlatform ().getComm ();
  TEUCHOS_TEST_FOR_EXCEPTION(
    comm->getSize () != 2, std::logic_error,
    "This test must be run with exactly 2 MPI processes.");

  const int myPID = comm->getRank();

  std::vector<GO> globalIDs;
  std::vector< std::vector<LO> > indices;
  std::vector< std::vector<ST> > values;
  if (myPID == 0) {
    globalIDs.push_back(0);
    globalIDs.push_back(1);
    globalIDs.push_back(2);
    values.resize(3);
    values[0].resize(2, 1);
    values[1].resize(2, 1);
    values[2].resize(2, 1);
    indices.resize(3);
    indices[0].push_back(0); indices[0].push_back(4);
    indices[1].push_back(0); indices[1].push_back(1);
    indices[2].push_back(1); indices[2].push_back(3);
  }
  else {
    globalIDs.push_back(3);
    globalIDs.push_back(4);
    values.resize(2);
    values[0].resize(2, 1);
    values[1].resize(2, 1);
    indices.resize(2);
    indices[0].push_back(2); indices[0].push_back(3);
    indices[1].push_back(2); indices[1].push_back(4);
  }

  const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
  const GO indexBase = 0;
  RCP<const map_type> rowMap =
    rcp (new map_type (INVALID, ArrayView<GO> (globalIDs), indexBase, comm));
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! rowMap->isOneToOne (), std::logic_error,
    "In this test, the row Map is supposed to be one to one.");

  RCP<crs_matrix_type> matrix = rcp (new crs_matrix_type (rowMap, 0));
  for (size_t i = 0; i < static_cast<size_t> (globalIDs.size ()); ++i) {
    matrix->insertGlobalValues (globalIDs[i],
                                ArrayView<const GO> (indices[i]),
                                ArrayView<const ST> (values[i]));
  }
  matrix->fillComplete ();

  vec_type x(rowMap), y(rowMap);
  x.putScalar (1.0);
  y.putScalar (1.0);
  const MT normBefore = y.norm2 ();
  if (myPID == 0) {
    std::cout << "norm of y before matrix->apply = " << normBefore << std::endl;
  }
  ST alpha (0.0), beta (1.0);
  matrix->apply (x, y, Teuchos::TRANS, alpha, beta);
  const MT normAfter = y.norm2 ();
  if (myPID == 0) {
    std::cout << "norm of y after matrix->apply = " << normAfter << std::endl;
  }

  TEUCHOS_TEST_FOR_EXCEPTION(
    normBefore != normAfter, std::logic_error,
    "normBefore = " << normBefore << " != normAfter = " << normAfter << ".");

  std::cout << "End Result: TEST PASSED" << std::endl;
  return EXIT_SUCCESS;
}

