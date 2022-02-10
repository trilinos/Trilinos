/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER
*/

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

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CrsMatrix, localValues, LO, GO, Scalar, Node)
{
  using matrix_t = Tpetra::CrsMatrix<Scalar,LO,GO,Node>;
  using map_t = Tpetra::Map<LO,GO,Node>;
  using impl_scalar_t = typename matrix_t::impl_scalar_type;

  const Tpetra::global_size_t dummy = 
        Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();

  auto comm = Tpetra::getDefaultComm();
  auto me = comm->getRank();
  auto np = comm->getSize();

  const size_t nRows = 10;
  Teuchos::RCP<const map_t> map = rcp(new map_t(dummy, nRows, 0, comm));

  // create a matrix with val(row,col) = row gid
  const int nPerRow = 5;
  matrix_t mat(map, nPerRow);
  Scalar vals[nPerRow];
  GO colinds[nPerRow];

  for (size_t i = 0; i < nRows; i++)
  {
    LO nnz = 0;
    GO gid = map->getGlobalElement(i);

    if (gid-2 >= map->getMinAllGlobalIndex()) {
      colinds[nnz] = gid-2;
      vals[nnz] = gid;
      nnz++;
    }
    if (gid-1 >= map->getMinAllGlobalIndex()) {
      colinds[nnz] = gid-1;
      vals[nnz] = gid;
      nnz++;
    }
    colinds[nnz] = gid;
    vals[nnz] = gid;
    nnz++;
    if (gid+1 <= map->getMaxAllGlobalIndex()) {
      colinds[nnz] = gid+1;
      vals[nnz] = gid;
      nnz++;
    }
    if (gid+2 <= map->getMaxAllGlobalIndex()) {
      colinds[nnz] = gid+2;
      vals[nnz] = gid;
      nnz++;
    }
    mat.insertGlobalValues(gid, nnz, vals, colinds);
  }
  mat.fillComplete();

  using host_range_policy = 
        Kokkos::RangePolicy<Kokkos::HostSpace::execution_space, LO>;
  using dev_range_policy = 
        Kokkos::RangePolicy<typename Node::device_type::execution_space, LO>;

  // check the host values view
  {
    int checkHostValuesErrors = 0;
    auto offs = mat.getLocalRowPtrsHost();
    auto inds = mat.getLocalIndicesHost();
    auto view = mat.getLocalValuesHost(Tpetra::Access::ReadOnly);
    Kokkos::parallel_reduce(
      "checkHostValues", 
      host_range_policy(0, nRows),
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
    auto offs = mat.getLocalRowPtrsDevice();
    auto inds = mat.getLocalIndicesDevice();
    auto view = mat.getLocalValuesDevice(Tpetra::Access::ReadOnly);
    auto lclMap = mat.getRowMap()->getLocalMap();
    Kokkos::parallel_reduce(
      "checkDeviceValues", 
      dev_range_policy(0, nRows),
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
      auto view = mat.getLocalValuesDevice(Tpetra::Access::ReadWrite);
      Kokkos::parallel_for(
        "modifyValuesOnDevice", 
        dev_range_policy(0, mat.getNodeNumEntries()),
        KOKKOS_LAMBDA (const LO i) { view[i] *= 2; }
        );
    }

    {
      int checkDeviceModifiedValuesOnHostErrors = 0;
      auto offs = mat.getLocalRowPtrsHost();
      auto inds = mat.getLocalIndicesHost();
      auto view = mat.getLocalValuesHost(Tpetra::Access::ReadOnly);
      Kokkos::parallel_reduce(
        "checkDeviceModifiedValuesOnHost", 
        host_range_policy(0, nRows),
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
      auto view = mat.getLocalValuesHost(Tpetra::Access::ReadWrite);
      Kokkos::parallel_for(
        "modifyValuesOnHost", 
        host_range_policy(0, mat.getNodeNumEntries()),
        [&] (const LO i) { view[i] *= 2; }
        );
    }
    {
      int checkHostModifiedValuesOnDeviceErrors = 0;
      auto offs = mat.getLocalRowPtrsDevice();
      auto inds = mat.getLocalIndicesDevice();
      auto view = mat.getLocalValuesDevice(Tpetra::Access::ReadOnly);
      auto lclMap = mat.getRowMap()->getLocalMap();
      Kokkos::parallel_reduce(
        "checkHostModifiedValuesOnDevice", 
        dev_range_policy(0, nRows),
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
      auto view = mat.getLocalValuesDevice(Tpetra::Access::OverwriteAll);
      Kokkos::parallel_for(
        "overwriteValuesOnDevice", 
        dev_range_policy(0, mat.getNodeNumEntries()),
        KOKKOS_LAMBDA (const LO i) { view[i] = me; }
        );
    }

    {
      int checkDeviceOverwrittenValuesOnHostErrors = 0;
      auto view = mat.getLocalValuesHost(Tpetra::Access::ReadOnly);
      Kokkos::parallel_reduce(
        "checkDeviceOverwrittenValuesOnHost", 
        host_range_policy(0, mat.getNodeNumEntries()),
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
      auto view = mat.getLocalValuesHost(Tpetra::Access::OverwriteAll);
      Kokkos::parallel_for(
        "overwriteValuesOnHost", 
        host_range_policy(0, mat.getNodeNumEntries()),
        [&] (const LO i) { view[i] = np; }
        );
    }
    {
      int checkHostOverwrittenValuesOnDeviceErrors = 0;
      auto view = mat.getLocalValuesDevice(Tpetra::Access::ReadOnly);
      Kokkos::parallel_reduce(
        "checkHostOverwrittenValuesOnDevice", 
        dev_range_policy(0, mat.getNodeNumEntries()),
        KOKKOS_LAMBDA (const LO i, int &lerr) {
          if (view[i] != impl_scalar_t(np)) lerr++;
        },
        checkHostOverwrittenValuesOnDeviceErrors);
      TEST_EQUALITY(checkHostOverwrittenValuesOnDeviceErrors, 0);
    }
  }
}


//
// INSTANTIATIONS
//

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
      TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( CrsMatrix, localValues, LO, GO, SCALAR, NODE ) 

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_SLGN( UNIT_TEST_GROUP )

}
