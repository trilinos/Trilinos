// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Testing CrsMatrix accessors
#include "Tpetra_TestingUtilities.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "TpetraUtils_MatrixGenerator.hpp"
#include <type_traits> // std::is_same

//
// UNIT TESTS
//

namespace {

template <typename Scalar, typename LO, typename GO, typename Node>
Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LO, GO, Node> >
buildMatrix(const Teuchos::RCP<const Teuchos::Comm<int> > &comm)
{
  using matrix_t = Tpetra::CrsMatrix<Scalar,LO,GO,Node>;
  using map_t = Tpetra::Map<LO,GO,Node>;

  const size_t nRows = 10;
  const Tpetra::global_size_t dummy = 
        Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
  Teuchos::RCP<const map_t> rowmap = rcp(new map_t(dummy, nRows, 0, comm));

  // create a matrix with val(row,col) = row gid
  const int nPerRow = 5;
  Teuchos::RCP<matrix_t> matrix = rcp(new matrix_t(rowmap, nPerRow));
  Scalar vals[nPerRow];
  GO colinds[nPerRow];

  for (size_t i = 0; i < nRows; i++)
  {
    LO nnz = 0;
    GO gid = rowmap->getGlobalElement(i);

    if (gid-2 >= rowmap->getMinAllGlobalIndex()) {
      colinds[nnz] = gid-2;
      vals[nnz] = gid;
      nnz++;
    }
    if (gid-1 >= rowmap->getMinAllGlobalIndex()) {
      colinds[nnz] = gid-1;
      vals[nnz] = gid;
      nnz++;
    }
    colinds[nnz] = gid;
    vals[nnz] = gid;
    nnz++;
    if (gid+1 <= rowmap->getMaxAllGlobalIndex()) {
      colinds[nnz] = gid+1;
      vals[nnz] = gid;
      nnz++;
    }
    if (gid+2 <= rowmap->getMaxAllGlobalIndex()) {
      colinds[nnz] = gid+2;
      vals[nnz] = gid;
      nnz++;
    }
    matrix->insertGlobalValues(gid, nnz, vals, colinds);
  }
  matrix->fillComplete();
  return matrix;
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CrsMatrix, localValues, LO, GO, Scalar, Node)
{
  using matrix_t = Tpetra::CrsMatrix<Scalar,LO,GO,Node>;
  using impl_scalar_t = typename matrix_t::impl_scalar_type;


  auto comm = Tpetra::getDefaultComm();
  auto me = comm->getRank();
  auto np = comm->getSize();

  auto mat = buildMatrix<Scalar, LO, GO, Node>(comm);

  using host_range_policy = 
        Kokkos::RangePolicy<Kokkos::HostSpace::execution_space, LO>;
  using dev_range_policy = 
        Kokkos::RangePolicy<typename Node::device_type::execution_space, LO>;

  // check the host values view
  {
    int checkHostValuesErrors = 0;
    auto offs = mat->getLocalRowPtrsHost();
    auto inds = mat->getLocalIndicesHost();
    auto view = mat->getLocalValuesHost(Tpetra::Access::ReadOnly);
    auto map = mat->getRowMap();
    Kokkos::parallel_reduce(
      "checkHostValues", 
      host_range_policy(0, mat->getLocalNumRows()),
      [&] (const LO i, int &lerr) {
        GO gid = map->getGlobalElement(i);
        for (size_t j = offs[i]; j < offs[i+1]; j++) {
          if (view[j] != impl_scalar_t(gid)) lerr++;
        } 
      },
      checkHostValuesErrors);
    TEST_EQUALITY(checkHostValuesErrors, 0);
  }

  // check the device values view
  {
    int checkDeviceValuesErrors = 0;
    auto offs = mat->getLocalRowPtrsDevice();
    auto inds = mat->getLocalIndicesDevice();
    auto view = mat->getLocalValuesDevice(Tpetra::Access::ReadOnly);
    auto lclMap = mat->getRowMap()->getLocalMap();
    Kokkos::parallel_reduce(
      "checkDeviceValues", 
      dev_range_policy(0, mat->getLocalNumRows()),
      KOKKOS_LAMBDA (const LO i, int &lerr) {
        GO gid = lclMap.getGlobalElement(i);
        for (size_t j = offs[i]; j < offs[i+1]; j++) {
          if (view[j] != impl_scalar_t(gid)) lerr++;
        } 
      },
      checkDeviceValuesErrors);
    TEST_EQUALITY(checkDeviceValuesErrors, 0);
  }

  // modify the views on device; check them on host
  {
    {
      auto view = mat->getLocalValuesDevice(Tpetra::Access::ReadWrite);
      Kokkos::parallel_for(
        "modifyValuesOnDevice", 
        dev_range_policy(0, mat->getLocalNumEntries()),
        KOKKOS_LAMBDA (const LO i) { view[i] *= 2; }
        );
    }

    {
      int checkDeviceModifiedValuesOnHostErrors = 0;
      auto offs = mat->getLocalRowPtrsHost();
      auto inds = mat->getLocalIndicesHost();
      auto view = mat->getLocalValuesHost(Tpetra::Access::ReadOnly);
      auto map = mat->getRowMap();
      Kokkos::parallel_reduce(
        "checkDeviceModifiedValuesOnHost", 
        host_range_policy(0, mat->getLocalNumRows()),
        [&] (const LO i, int &lerr) {
          GO gid = map->getGlobalElement(i);
          for (size_t j = offs[i]; j < offs[i+1]; j++) {
            if (view[j] != impl_scalar_t(2 * gid)) lerr++;
          } 
        },
        checkDeviceModifiedValuesOnHostErrors);
      TEST_EQUALITY(checkDeviceModifiedValuesOnHostErrors, 0);
    }
  }

  // modify the views on host; check them on device
  {
    {
      auto view = mat->getLocalValuesHost(Tpetra::Access::ReadWrite);
      Kokkos::parallel_for(
        "modifyValuesOnHost", 
        host_range_policy(0, mat->getLocalNumEntries()),
        [&] (const LO i) { view[i] *= 2; }
        );
    }
    {
      int checkHostModifiedValuesOnDeviceErrors = 0;
      auto offs = mat->getLocalRowPtrsDevice();
      auto inds = mat->getLocalIndicesDevice();
      auto view = mat->getLocalValuesDevice(Tpetra::Access::ReadOnly);
      auto lclMap = mat->getRowMap()->getLocalMap();
      Kokkos::parallel_reduce(
        "checkHostModifiedValuesOnDevice", 
        dev_range_policy(0, mat->getLocalNumRows()),
        KOKKOS_LAMBDA (const LO i, int &lerr) {
          GO gid = lclMap.getGlobalElement(i);
          for (size_t j = offs[i]; j < offs[i+1]; j++) {
            if (view[j] != impl_scalar_t(4 * gid)) lerr++;
          } 
        },
        checkHostModifiedValuesOnDeviceErrors);
      TEST_EQUALITY(checkHostModifiedValuesOnDeviceErrors, 0);
    }
  }

  // overwrite the views on device; check them on host
  {
    {
      auto view = mat->getLocalValuesDevice(Tpetra::Access::OverwriteAll);
      Kokkos::parallel_for(
        "overwriteValuesOnDevice", 
        dev_range_policy(0, mat->getLocalNumEntries()),
        KOKKOS_LAMBDA (const LO i) { view[i] = me; }
        );
    }

    {
      int checkDeviceOverwrittenValuesOnHostErrors = 0;
      auto view = mat->getLocalValuesHost(Tpetra::Access::ReadOnly);
      Kokkos::parallel_reduce(
        "checkDeviceOverwrittenValuesOnHost", 
        host_range_policy(0, mat->getLocalNumEntries()),
        [&] (const LO i, int &lerr) {
          if (view[i] != impl_scalar_t(me)) lerr++;
        },
        checkDeviceOverwrittenValuesOnHostErrors);
      TEST_EQUALITY(checkDeviceOverwrittenValuesOnHostErrors, 0);
    }
  }

  // overwrite the views on host; check them on device
  {
    {
      auto view = mat->getLocalValuesHost(Tpetra::Access::OverwriteAll);
      Kokkos::parallel_for(
        "overwriteValuesOnHost", 
        host_range_policy(0, mat->getLocalNumEntries()),
        [&] (const LO i) { view[i] = np; }
        );
    }
    {
      int checkHostOverwrittenValuesOnDeviceErrors = 0;
      auto view = mat->getLocalValuesDevice(Tpetra::Access::ReadOnly);
      Kokkos::parallel_reduce(
        "checkHostOverwrittenValuesOnDevice", 
        dev_range_policy(0, mat->getLocalNumEntries()),
        KOKKOS_LAMBDA (const LO i, int &lerr) {
          if (view[i] != impl_scalar_t(np)) lerr++;
        },
        checkHostOverwrittenValuesOnDeviceErrors);
      TEST_EQUALITY(checkHostOverwrittenValuesOnDeviceErrors, 0);
    }
  }
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CrsMatrix, copyOrView, LO, GO, Scalar, Node)
{
  // Test CrsMatrix constructor that accepts Teuchos' copyOrView argument
  auto comm = Tpetra::getDefaultComm();

  using host_range_policy = 
        Kokkos::RangePolicy<Kokkos::HostSpace::execution_space, LO>;

  using matrix_t = Tpetra::CrsMatrix<Scalar,LO,GO,Node>;
  Teuchos::RCP<matrix_t> origMat = buildMatrix<Scalar, LO, GO, Node>(comm);
  origMat->setAllToScalar(1965);

  {
    // Create a deep copy, change its values, and 
    // confirm that the matrices have different values
    matrix_t copyMat(*origMat, Teuchos::Copy);

    copyMat.setAllToScalar(2020);

    auto origVals = origMat->getLocalValuesHost(Tpetra::Access::ReadOnly);
    auto copyVals = copyMat.getLocalValuesHost(Tpetra::Access::ReadOnly);
    int incorrectMatchingValues = 0;
    Kokkos::parallel_reduce(
      "incorrectMatchingValues",
      host_range_policy(0, origMat->getLocalNumEntries()),
      [&] (const LO i, int &lerr) { if (origVals[i] == copyVals[i]) lerr++; },
      incorrectMatchingValues);
    TEST_EQUALITY(incorrectMatchingValues, 0);
  }
    
  {
    // Create a view copy, change its values, and 
    // confirm that the matrices have different values
    matrix_t viewMat(*origMat, Teuchos::View);

    viewMat.setAllToScalar(2020);
  
    auto origVals = origMat->getLocalValuesHost(Tpetra::Access::ReadOnly);
    auto viewVals = viewMat.getLocalValuesHost(Tpetra::Access::ReadOnly);
    int incorrectDifferingValues = 0;
    Kokkos::parallel_reduce(
      "incorrectMatchingValues",
      host_range_policy(0, origMat->getLocalNumEntries()),
      [&] (const LO i, int &lerr) { if (origVals[i] != viewVals[i]) lerr++; },
      incorrectDifferingValues);
    TEST_EQUALITY(incorrectDifferingValues, 0);
  }
}

//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, localValues, LO, GO, SCALAR, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, copyOrView, LO, GO, SCALAR, NODE ) 

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_SLGN( UNIT_TEST_GROUP )

}
