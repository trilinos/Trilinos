// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <iostream>
#include <algorithm>  // shuffle
#include <vector>     // vector
#include <random>

#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_IO.hpp>

// Teuchos
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_StackedTimer.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraUtils.hpp>
#include <Galeri_XpetraMaps.hpp>
// MueLu
#include "MueLu.hpp"
#include "MueLu_TestHelpers.hpp"
#include "MueLu_PerfModels.hpp"
#include "MueLu_PerfModelReporter.hpp"
#include <MatrixLoad.hpp>

#include "Xpetra_TpetraMultiVector.hpp"
#include "Xpetra_TpetraImport.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_MultiVector.hpp"
#include "KokkosSparse_spmv.hpp"
#include "Kokkos_Core.hpp"

#if defined(HAVE_MUELU_CUSPARSE)
#include "cublas_v2.h"
#include "cusparse.h"
#endif

#if defined(HAVE_MUELU_HYPRE)
#include "HYPRE_IJ_mv.h"
#include "HYPRE_parcsr_ls.h"
#include "HYPRE_parcsr_mv.h"
#include "HYPRE.h"
#endif

#if defined(HAVE_MUELU_PETSC)
#include "petscksp.h"
#endif

Teuchos::RCP<Teuchos::StackedTimer> stacked_timer;
Teuchos::RCP<Teuchos::TimeMonitor> globalTimeMonitor;

// =========================================================================
// Support Routines
// =========================================================================

template <class View1, class View2>
inline void copy_view_n(const int n, const View1 x1, View2 x2) {
  Kokkos::parallel_for(
      n, KOKKOS_LAMBDA(const size_t i) {
        x2[i] = x1[i];
      });
}

template <class View1, class View2>
inline void copy_view(const View1 x1, View2 x2) {
  copy_view_n(x1.extent(0), x1, x2);
}

template <class V1, class V2>
void print_crs_graph(std::string name, const V1 rowptr, const V2 colind) {
  printf("%s rowptr[%d] = ", name.c_str(), rowptr.extent(0));
  for (size_t i = 0; i < rowptr.extent(0); i++)
    printf(" %d", (int)rowptr[i]);
  printf("\n%s colind[%d] = ", name.c_str(), colind.extent(0));
  for (size_t i = 0; i < colind.extent(0); i++)
    printf(" %d", (int)colind[i]);
  printf("\n");
}

//==============================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MakeContiguousMaps(const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Matrix,
                        Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> >& ContiguousRowMap,
                        Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> >& ContiguousColumnMap,
                        Teuchos::RCP<const Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> >& ContiguousDomainMap) {
  using Teuchos::rcp;
  using Teuchos::RCP;
  using LO             = LocalOrdinal;
  using GO             = GlobalOrdinal;
  using NO             = Node;
  using map_type       = Tpetra::Map<LO, GO, NO>;
  using import_type    = Tpetra::Import<LO, GO, NO>;
  using go_vector_type = Tpetra::Vector<GO, LO, GO, NO>;

  // Must create GloballyContiguous DomainMap (which is a permutation of Matrix_'s
  // DomainMap) and the corresponding permuted ColumnMap.
  //   Epetra_GID  --------->   LID   ----------> HYPRE_GID
  //           via DomainMap.LID()       via GloballyContiguousDomainMap.GID()
  RCP<const map_type> RowMap      = Matrix.getRowMap();
  RCP<const map_type> DomainMap   = Matrix.getDomainMap();
  RCP<const map_type> ColumnMap   = Matrix.getColMap();
  RCP<const import_type> importer = Matrix.getGraph()->getImporter();

  if (RowMap->isContiguous()) {
    // If the row map is linear, then we can just use the row map as is.
    ContiguousRowMap = RowMap;
  } else {
    // The row map isn't linear, so we need a new row map
    ContiguousRowMap = rcp(new map_type(Matrix.getRowMap()->getGlobalNumElements(),
                                        Matrix.getRowMap()->getLocalNumElements(), 0, Matrix.getRowMap()->getComm()));
  }

  if (DomainMap->isContiguous()) {
    // If the domain map is linear, then we can just use the column map as is.
    ContiguousDomainMap = DomainMap;
    ContiguousColumnMap = ColumnMap;
  } else {
    // The domain map isn't linear, so we need a new domain map
    ContiguousDomainMap = rcp(new map_type(DomainMap->getGlobalNumElements(),
                                           DomainMap->getLocalNumElements(), 0, DomainMap->getComm()));
    if (importer) {
      // If there's an importer then we can use it to get a new column map
      go_vector_type MyGIDsHYPRE(DomainMap, ContiguousDomainMap->getLocalElementList());

      // import the HYPRE GIDs
      go_vector_type ColGIDsHYPRE(ColumnMap);
      ColGIDsHYPRE.doImport(MyGIDsHYPRE, *importer, Tpetra::INSERT);

      // Make a HYPRE numbering-based column map.
      ContiguousColumnMap = rcp(new map_type(ColumnMap->getGlobalNumElements(), ColGIDsHYPRE.getDataNonConst()(), 0, ColumnMap->getComm()));
    } else {
      // The problem has matching domain/column maps, and somehow the domain map isn't linear, so just use the new domain map
      ContiguousColumnMap = rcp(new map_type(ColumnMap->getGlobalNumElements(), ContiguousDomainMap->getLocalElementList(), 0, ColumnMap->getComm()));
    }
  }
}

// =========================================================================
// PETSC Testing
// =========================================================================
#if defined(HAVE_MUELU_PETSC) && defined(HAVE_MPI)

/*#define PETSC_CHK_ERR(x)  {                   \
  PetscCall(x); \
  }*/

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
class Petsc_SpmV_Pack {
  typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> crs_matrix_type;
  typedef Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> vector_type;
  typedef Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> map_type;

 public:
  Petsc_SpmV_Pack(const crs_matrix_type& A,
                  const vector_type& X,
                  vector_type& Y) {
    using LO             = LocalOrdinal;
    const MPI_Comm& comm = *Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >(A.getRowMap()->getComm())->getRawMpiComm();
    LO Nx                = (LO)X.getMap()->getLocalNumElements();
    LO Ny                = (LO)Y.getMap()->getLocalNumElements();

    // NOTE:  Petsc requires contiguous GIDs for the row map
    Teuchos::RCP<const map_type> c_row_map, c_col_map, c_domain_map;
    MakeContiguousMaps(A, c_row_map, c_col_map, c_domain_map);

    // Note: PETSC appears to favor local indices for vector insertion
    std::vector<PetscInt> l_indices(std::max(Nx, Ny));
    for (int i = 0; i < std::max(Nx, Ny); i++)
      l_indices[i] = i;

    // x vector
    VecCreate(comm, &x_p);
    VecSetType(x_p, VECMPI);
    VecSetSizes(x_p, Nx, X.getMap()->getGlobalNumElements());
    VecSetValues(x_p, Nx, l_indices.data(), X.getData(0).getRawPtr(), INSERT_VALUES);
    VecAssemblyBegin(x_p);
    VecAssemblyEnd(x_p);

    // y vector
    VecCreate(comm, &y_p);
    VecSetType(y_p, VECMPI);
    VecSetSizes(y_p, Ny, Y.getMap()->getGlobalNumElements());
    VecSetValues(y_p, Ny, l_indices.data(), Y.getData(0).getRawPtr(), INSERT_VALUES);
    VecAssemblyBegin(y_p);
    VecAssemblyEnd(y_p);

    // A matrix
    // NOTE: This is an overallocation, but that's OK
    PetscInt max_nnz = A.getLocalMaxNumRowEntries();
    MatCreateAIJ(comm, c_row_map->getLocalNumElements(), c_domain_map->getLocalNumElements(), PETSC_DECIDE, PETSC_DECIDE,
                 max_nnz, NULL, max_nnz, NULL, &A_p);

    std::vector<PetscInt> new_indices(max_nnz);
    for (LO i = 0; i < (LO)A.getLocalNumRows(); i++) {
      typename crs_matrix_type::values_host_view_type values;
      typename crs_matrix_type::local_inds_host_view_type indices;
      A.getLocalRowView(i, indices, values);
      for (LO j = 0; j < (LO)indices.extent(0); j++) {
        new_indices[j] = c_col_map->getGlobalElement(indices(j));
      }
      PetscInt GlobalRow[1];
      PetscInt numEntries = (PetscInt)indices.extent(0);
      GlobalRow[0]        = c_row_map->getGlobalElement(i);

      MatSetValues(A_p, 1, GlobalRow, numEntries, new_indices.data(), values.data(), INSERT_VALUES);
    }
    MatAssemblyBegin(A_p, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A_p, MAT_FINAL_ASSEMBLY);
  }

  ~Petsc_SpmV_Pack() {
    VecDestroy(&x_p);
    VecDestroy(&y_p);
    MatDestroy(&A_p);
  }

  bool spmv(const Scalar alpha, const Scalar beta) {
    int rv = MatMult(A_p, x_p, y_p);
    Kokkos::fence();
    return (rv != 0);
  }

 private:
  Mat A_p;
  Vec x_p, y_p;
};

#endif

// =========================================================================
// HYPRE Testing
// =========================================================================
#if defined(HAVE_MUELU_HYPRE) && defined(HAVE_MPI)

#define HYPRE_CHK_ERR(x)                                                              \
  {                                                                                   \
    if (x != 0) throw std::runtime_error("ERROR: HYPRE returned non-zero exit code"); \
  }

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
class HYPRE_SpmV_Pack {
  typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> crs_matrix_type;
  typedef Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> vector_type;
  typedef Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node> map_type;

 public:
  HYPRE_SpmV_Pack(const crs_matrix_type& A,
                  const vector_type& X,
                  vector_type& Y) {
    using LO = LocalOrdinal;
    using GO = GlobalOrdinal;

    // Time to copy the matrix / vector over to HYPRE datastructures
    const MPI_Comm& comm = *Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >(A.getRowMap()->getComm())->getRawMpiComm();

    // NOTE:  Hypre requires contiguous GIDs for the row map
    Teuchos::RCP<const map_type> c_row_map, c_col_map, c_domain_map;
    MakeContiguousMaps(A, c_row_map, c_col_map, c_domain_map);

    // Create matrix
    int row_lo = c_row_map->getMinGlobalIndex();
    int row_hi = c_row_map->getMaxGlobalIndex();
    int dom_lo = c_domain_map->getMinGlobalIndex();
    int dom_hi = c_domain_map->getMaxGlobalIndex();
    HYPRE_CHK_ERR(HYPRE_IJMatrixCreate(comm, row_lo, row_hi, dom_lo, dom_hi, &ij_matrix));
    HYPRE_CHK_ERR(HYPRE_IJMatrixSetObjectType(ij_matrix, HYPRE_PARCSR));
    HYPRE_CHK_ERR(HYPRE_IJMatrixInitialize(ij_matrix));

    // Fill matrix
    std::vector<GO> new_indices(A.getLocalMaxNumRowEntries());
    for (LO i = 0; i < (LO)A.getLocalNumRows(); i++) {
      typename crs_matrix_type::values_host_view_type values;
      typename crs_matrix_type::local_inds_host_view_type indices;
      A.getLocalRowView(i, indices, values);
      for (LO j = 0; j < (LO)indices.extent(0); j++) {
        new_indices[j] = c_col_map->getGlobalElement(indices(j));
      }
      GO GlobalRow[1];
      GO numEntries = (GO)indices.extent(0);
      GlobalRow[0]  = c_row_map->getGlobalElement(i);
      HYPRE_CHK_ERR(HYPRE_IJMatrixSetValues(ij_matrix, 1, &numEntries, GlobalRow, new_indices.data(), values.data()));
    }
    HYPRE_CHK_ERR(HYPRE_IJMatrixAssemble(ij_matrix));
    HYPRE_CHK_ERR(HYPRE_IJMatrixGetObject(ij_matrix, (void**)&parcsr_matrix));

    // Now the x vector
    GO* dom_indices = const_cast<GO*>(c_domain_map->getLocalElementList().getRawPtr());
    HYPRE_CHK_ERR(HYPRE_IJVectorCreate(comm, dom_lo, dom_hi, &x_ij));
    HYPRE_CHK_ERR(HYPRE_IJVectorSetObjectType(x_ij, HYPRE_PARCSR));
    HYPRE_CHK_ERR(HYPRE_IJVectorInitialize(x_ij));
    HYPRE_CHK_ERR(HYPRE_IJVectorSetValues(x_ij, X.getLocalLength(), dom_indices, const_cast<vector_type*>(&X)->getDataNonConst(0).getRawPtr()));
    HYPRE_CHK_ERR(HYPRE_IJVectorAssemble(x_ij));
    HYPRE_CHK_ERR(HYPRE_IJVectorGetObject(x_ij, (void**)&x_par));

    // Now the y vector
    GO* row_indices = const_cast<GO*>(c_row_map->getLocalElementList().getRawPtr());
    HYPRE_CHK_ERR(HYPRE_IJVectorCreate(comm, row_lo, row_hi, &y_ij));
    HYPRE_CHK_ERR(HYPRE_IJVectorSetObjectType(y_ij, HYPRE_PARCSR));
    HYPRE_CHK_ERR(HYPRE_IJVectorInitialize(y_ij));
    HYPRE_CHK_ERR(HYPRE_IJVectorSetValues(y_ij, Y.getLocalLength(), row_indices, Y.getDataNonConst(0).getRawPtr()));
    HYPRE_CHK_ERR(HYPRE_IJVectorAssemble(y_ij));
    HYPRE_CHK_ERR(HYPRE_IJVectorGetObject(y_ij, (void**)&y_par));
  }

  ~HYPRE_SpmV_Pack() {
    HYPRE_IJMatrixDestroy(ij_matrix);
    HYPRE_IJVectorDestroy(x_ij);
    HYPRE_IJVectorDestroy(y_ij);
  }

  bool spmv(const Scalar alpha, const Scalar beta) {
    int rv = HYPRE_ParCSRMatrixMatvec(alpha, parcsr_matrix, x_par, beta, y_par);
    Kokkos::fence();
    return (rv != 0);
  }

 private:
  // NOTE:  The par matrix/vectors are just pointers to the internal data of the ij guys.  They do not need to
  // be deallocated on their own.
  HYPRE_IJMatrix ij_matrix;
  HYPRE_ParCSRMatrix parcsr_matrix;

  HYPRE_IJVector x_ij, y_ij;
  HYPRE_ParVector x_par, y_par;
};

#endif

// =========================================================================
// MAGMASparse Testing
// =========================================================================
#if defined(HAVE_MUELU_MAGMASPARSE)
#include <magma_v2.h>
#include <magmasparse.h>
template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
class MagmaSparse_SpmV_Pack {
  typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> crs_matrix_type;
  typedef Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> vector_type;

 public:
  MagmaSparse_SpmV_Pack(const crs_matrix_type& A,
                        const vector_type& X,
                        vector_type& Y) {}

  ~MagmaSparse_SpmV_Pack(){};

  bool spmv(const Scalar alpha, const Scalar beta) { return (true); }
};

template <typename LocalOrdinal, typename GlobalOrdinal>
class MagmaSparse_SpmV_Pack<double, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosCudaWrapperNode> {
  // typedefs shared among other TPLs
  typedef double Scalar;
  typedef typename Tpetra::KokkosCompat::KokkosCudaWrapperNode Node;
  typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> crs_matrix_type;
  typedef Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> vector_type;
  typedef typename crs_matrix_type::local_matrix_host_type KCRS;
  typedef typename KCRS::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::non_const_type lno_view_t;
  typedef typename graph_t::row_map_type::const_type c_lno_view_t;
  typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
  typedef typename graph_t::entries_type::const_type c_lno_nnz_view_t;
  typedef typename KCRS::values_type::non_const_type scalar_view_t;
  typedef typename Node::device_type device_type;

  typedef typename Kokkos::View<int*,
                                typename lno_nnz_view_t::array_layout,
                                typename lno_nnz_view_t::device_type>
      int_array_type;

 public:
  MagmaSparse_SpmV_Pack(const crs_matrix_type& A,
                        const vector_type& X,
                        vector_type& Y) {
    // data access common to other TPLs
    const KCRS& Amat          = A.getLocalMatrixHost();
    c_lno_view_t Arowptr      = Amat.graph.row_map;
    c_lno_nnz_view_t Acolind  = Amat.graph.entries;
    const scalar_view_t Avals = Amat.values;

    // type conversion
    Arowptr_int = int_array_type("Arowptr", Arowptr.extent(0));
    // copy the ordinals into the local view (type conversion)
    copy_view(Arowptr, Arowptr_int);

    m      = static_cast<int>(A.getLocalNumRows());
    n      = static_cast<int>(A.getLocalNumCols());
    nnz    = static_cast<int>(Acolind.extent(0));
    vals   = reinterpret_cast<Scalar*>(Avals.data());
    cols   = const_cast<int*>(reinterpret_cast<const int*>(Acolind.data()));
    rowptr = reinterpret_cast<int*>(Arowptr_int.data());

    auto X_lcl = X.template getLocalView<device_type>();
    auto Y_lcl = Y.template getLocalView<device_type>();
    x          = reinterpret_cast<Scalar*>(X_lcl.data());
    y          = reinterpret_cast<Scalar*>(Y_lcl.data());

    magma_init();
    int device;
    magma_getdevice(&device);
    magma_queue_create(device, &queue);

    // create the magma data container
    magma_dev_Acrs = {Magma_CSR};
    magma_Acrs     = {Magma_CSR};
    magma_dev_x    = {Magma_DENSE};
    magma_dev_y    = {Magma_DENSE};

    //
    magma_dvset_dev(m, 1,          // assume mx1 size
                    x,             // ptr to data on the device
                    &magma_dev_x,  // magma vector to populate
                    queue);

    magma_dvset_dev(m, 1,          // assume mx1 size
                    y,             // ptr to data on the device
                    &magma_dev_y,  // magma vector to populate
                    queue);
    magma_dcsrset(m, n, rowptr, cols, vals, &magma_Acrs, queue);
    magma_dmtransfer(magma_Acrs, &magma_dev_Acrs, Magma_DEV, Magma_DEV, queue);
  }

  ~MagmaSparse_SpmV_Pack() {
    // we dont' own the data... not sure what magma does here...
    magma_dmfree(&magma_dev_x, queue);
    magma_dmfree(&magma_dev_y, queue);
    magma_dmfree(&magma_dev_Acrs, queue);

    magma_finalize();
  };

  bool spmv(const Scalar alpha, const Scalar beta) {
    magma_d_spmv(1.0, magma_dev_Acrs, magma_dev_x, 0.0, magma_dev_y, queue);
    Kokkos::fence();
  }

 private:
  magma_d_matrix magma_dev_Acrs;
  magma_d_matrix magma_Acrs;
  magma_d_matrix magma_dev_x;
  magma_d_matrix magma_dev_y;
  magma_queue_t queue;

  int m        = -1;
  int n        = -1;
  int nnz      = -1;
  Scalar* vals = nullptr;  // aliased
  int* cols    = nullptr;  // aliased
  int* rowptr  = nullptr;  // copied
  Scalar* x    = nullptr;  // aliased
  Scalar* y    = nullptr;  // aliased

  // handles to the copied data
  int_array_type Arowptr_int;
};
#endif

// =========================================================================
// CuSparse Testing
// =========================================================================
#if defined(HAVE_MUELU_CUSPARSE)

#define CHECK_CUDA(func)                                         \
  {                                                              \
    cudaError_t status = (func);                                 \
    if (status != cudaSuccess) {                                 \
      printf("CUDA API failed at line %d with error: %s (%d)\n", \
             __LINE__, cudaGetErrorString(status), status);      \
    }                                                            \
  }

#define CHECK_CUSPARSE(func)                                         \
  {                                                                  \
    cusparseStatus_t status = (func);                                \
    if (status != CUSPARSE_STATUS_SUCCESS) {                         \
      printf("CUSPARSE API failed at line %d with error: %s (%d)\n", \
             __LINE__, cusparseGetErrorString(status), status);      \
    }                                                                \
  }

template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
class CuSparse_SpmV_Pack {
  typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> crs_matrix_type;
  typedef Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> vector_type;

 public:
  CuSparse_SpmV_Pack(const crs_matrix_type& A,
                     const vector_type& X,
                     vector_type& Y) {}

  ~CuSparse_SpmV_Pack(){};

  cusparseStatus_t spmv(const Scalar alpha, const Scalar beta) { return CUSPARSE_STATUS_SUCCESS; }
};

template <typename LocalOrdinal, typename GlobalOrdinal>
class CuSparse_SpmV_Pack<double, LocalOrdinal, GlobalOrdinal, Tpetra::KokkosCompat::KokkosCudaWrapperNode> {
  // typedefs shared among other TPLs
  typedef double Scalar;
  typedef typename Tpetra::KokkosCompat::KokkosCudaWrapperNode Node;
  typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> crs_matrix_type;
  typedef Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> vector_type;
  typedef typename crs_matrix_type::local_matrix_device_type KCRS;
  typedef typename KCRS::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::non_const_type lno_view_t;
  typedef typename graph_t::row_map_type::const_type c_lno_view_t;
  typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
  typedef typename graph_t::entries_type::const_type c_lno_nnz_view_t;
  typedef typename KCRS::values_type::non_const_type scalar_view_t;
  typedef typename Node::device_type device_type;

  typedef typename Kokkos::View<int*,
                                typename lno_nnz_view_t::array_layout,
                                typename lno_nnz_view_t::device_type>
      cusparse_int_type;

 public:
  CuSparse_SpmV_Pack(const crs_matrix_type& A,
                     const vector_type& X,
                     vector_type& Y) {
    // data access common to other TPLs
    const KCRS& Amat          = A.getLocalMatrixDevice();
    c_lno_view_t Arowptr      = Amat.graph.row_map;
    c_lno_nnz_view_t Acolind  = Amat.graph.entries;
    const scalar_view_t Avals = Amat.values;

    Arowptr_cusparse = cusparse_int_type("Arowptr", Arowptr.extent(0));
    Acolind_cusparse = cusparse_int_type("Acolind", Acolind.extent(0));
    // copy the ordinals into the local view (type conversion)
    copy_view(Arowptr, Arowptr_cusparse);
    copy_view(Acolind, Acolind_cusparse);

    m      = static_cast<int>(A.getLocalNumRows());
    n      = static_cast<int>(A.getLocalNumCols());
    nnz    = static_cast<int>(Acolind_cusparse.extent(0));
    vals   = reinterpret_cast<Scalar*>(Avals.data());
    cols   = reinterpret_cast<int*>(Acolind_cusparse.data());
    rowptr = reinterpret_cast<int*>(Arowptr_cusparse.data());

    auto X_lcl = X.getLocalViewDevice(Tpetra::Access::ReadOnly);
    auto Y_lcl = Y.getLocalViewDevice(Tpetra::Access::ReadWrite);
    x          = const_cast<Scalar*>(reinterpret_cast<const Scalar*>(X_lcl.data()));
    y          = reinterpret_cast<Scalar*>(Y_lcl.data());

    /* Get handle to the CUBLAS context */
    cublasCreate(&cublasHandle);
    /* Get handle to the CUSPARSE context */
    cusparseCreate(&cusparseHandle);

    CHECK_CUSPARSE(cusparseCreateCsr(&descrA, m, n, nnz,
                                     rowptr, cols, vals,
                                     CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
                                     CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F));

    CHECK_CUSPARSE(cusparseCreateDnVec(&vecY, m, y, CUDA_R_64F));
    CHECK_CUSPARSE(cusparseCreateDnVec(&vecX, n, x, CUDA_R_64F));
  }

  ~CuSparse_SpmV_Pack() {
    CHECK_CUSPARSE(cusparseDestroySpMat(descrA));
    CHECK_CUSPARSE(cusparseDestroyDnVec(vecX));
    CHECK_CUSPARSE(cusparseDestroyDnVec(vecY));
    cusparseDestroy(cusparseHandle);
    cublasDestroy(cublasHandle);
  }

  cusparseStatus_t spmv(const Scalar alpha, const Scalar beta) {
    // compute: y = alpha*Ax + beta*y

#if CUSPARSE_VERSION >= 11201
    cusparseSpMVAlg_t alg = CUSPARSE_SPMV_ALG_DEFAULT;
#else
    cusparseSpMVAlg_t alg = CUSPARSE_MV_ALG_DEFAULT;
#endif

    size_t bufferSize;
    CHECK_CUSPARSE(cusparseSpMV_bufferSize(cusparseHandle,
                                           transA,
                                           &alpha,
                                           descrA,
                                           vecX,
                                           &beta,
                                           vecY,
                                           CUDA_R_64F,
                                           alg,
                                           &bufferSize));

    void* dBuffer = NULL;
    CHECK_CUDA(cudaMalloc(&dBuffer, bufferSize));

    cusparseStatus_t rc = cusparseSpMV(cusparseHandle,
                                       transA,
                                       &alpha,
                                       descrA,
                                       vecX,
                                       &beta,
                                       vecY,
                                       CUDA_R_64F,
                                       alg,
                                       dBuffer);

    CHECK_CUDA(cudaFree(dBuffer));
    Kokkos::fence();
    return (rc);
  }

 private:
  cublasHandle_t cublasHandle     = 0;
  cusparseHandle_t cusparseHandle = 0;
  cusparseSpMatDescr_t descrA     = 0;
  cusparseDnVecDescr_t vecX = 0, vecY = 0;

  // CUSPARSE_OPERATION_NON_TRANSPOSE
  // CUSPARSE_OPERATION_TRANSPOSE
  // CUSPARSE_OPERATION_CONJUGATE_TRANSPOSE
  cusparseOperation_t transA = CUSPARSE_OPERATION_NON_TRANSPOSE;
  int m                      = -1;
  int n                      = -1;
  int nnz                    = -1;
  Scalar* vals               = nullptr;  // aliased
  int* cols                  = nullptr;  // copied
  int* rowptr                = nullptr;  // copied
  Scalar* x                  = nullptr;  // aliased
  Scalar* y                  = nullptr;  // aliased

  // handles to the copied data
  cusparse_int_type Arowptr_cusparse;
  cusparse_int_type Acolind_cusparse;
};
#endif

// =========================================================================
// MKL Testing
// =========================================================================
#if defined(HAVE_MUELU_MKL)
#include "mkl.h"

std::string mkl_error(sparse_status_t code) {
  switch (code) {
    case SPARSE_STATUS_SUCCESS:
      return std::string("Success");
    case SPARSE_STATUS_NOT_INITIALIZED:
      return std::string("Empty handle or matrix array");
    case SPARSE_STATUS_ALLOC_FAILED:
      return std::string("Memory allocation failed");
    case SPARSE_STATUS_INVALID_VALUE:
      return std::string("Input contains an invalid value");
    case SPARSE_STATUS_EXECUTION_FAILED:
      return std::string("Execution failed");
    case SPARSE_STATUS_INTERNAL_ERROR:
      return std::string("Internal error");
    case SPARSE_STATUS_NOT_SUPPORTED:
      return std::string("Operation not supported");
  };
}

struct matrix_descr mkl_descr;

void MV_MKL(sparse_matrix_t& AMKL, double* x, double* y) {
  // sparse_status_t mkl_sparse_d_mv (sparse_operation_t operation, double alpha, const sparse_matrix_t A, struct matrix_descr descr, const double *x, double beta, double *y);

  mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0, AMKL, mkl_descr, x, 0.0, y);
  Kokkos::fence();
}

#endif

// =========================================================================
// Tpetra Kernel Testing
// =========================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MV_Tpetra(const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A, const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x, Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& y) {
  A.apply(x, y);
  Kokkos::fence();
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void MV_KK(const Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& A, const Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& x, Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& y) {
  typedef typename Node::device_type device_type;
  const auto& AK = A.getLocalMatrixDevice();
  auto X_lcl     = x.getLocalViewDevice(Tpetra::Access::ReadOnly);
  auto Y_lcl     = y.getLocalViewDevice(Tpetra::Access::OverwriteAll);
  KokkosSparse::spmv(KokkosSparse::NoTranspose, Teuchos::ScalarTraits<Scalar>::one(), AK, X_lcl, Teuchos::ScalarTraits<Scalar>::zero(), Y_lcl);
  Kokkos::fence();
}

// =========================================================================
// =========================================================================
// =========================================================================
template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor& clp, Xpetra::UnderlyingLib& lib, int argc, char* argv[]) {
#include <MueLu_UseShortNames.hpp>

  using std::endl;
  using Teuchos::RCP;  // reference count pointers
  using Teuchos::rcp;
  using Teuchos::TimeMonitor;

  bool success = false;
  bool verbose = true;
    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
    int numProc                         = comm->getSize();

    // =========================================================================
    // Convenient definitions
    // =========================================================================
    //    SC zero = Teuchos::ScalarTraits<SC>::zero(), one = Teuchos::ScalarTraits<SC>::one();

    // Instead of checking each time for rank, create a rank 0 stream
    RCP<Teuchos::FancyOStream> pOut = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    Teuchos::FancyOStream& out      = *pOut;
    out.setOutputToRootOnly(0);

    // =========================================================================
    // Parameters initialization
    // =========================================================================
    GO nx = 100, ny = 100, nz = 50;
    Galeri::Xpetra::Parameters<GO> galeriParameters(clp, nx, ny, nz, "Laplace2D");  // manage parameters of the test case
    Xpetra::Parameters xpetraParameters(clp);                                       // manage parameters of Xpetra

    bool binaryFormat = false;
    clp.setOption("binary", "ascii", &binaryFormat, "print timings to screen");
    std::string matrixFile;
    clp.setOption("matrixfile", &matrixFile, "matrix market file containing matrix");
    std::string rowMapFile;
    clp.setOption("rowmap", &rowMapFile, "map data file");
    std::string colMapFile;
    clp.setOption("colmap", &colMapFile, "colmap data file");
    std::string domainMapFile;
    clp.setOption("domainmap", &domainMapFile, "domainmap data file");
    std::string rangeMapFile;
    clp.setOption("rangemap", &rangeMapFile, "rangemap data file");

    bool printTimings = true;
    clp.setOption("timings", "notimings", &printTimings, "print timings to screen");
    int nrepeat = 1;
    clp.setOption("nrepeat", &nrepeat, "repeat the experiment N times");

    bool describeMatrix = true;
    clp.setOption("showmatrix", "noshowmatrix", &describeMatrix, "describe matrix");
    bool useStackedTimer = false;
    clp.setOption("stackedtimer", "nostackedtimer", &useStackedTimer, "use stacked timer");
    bool verboseModel = false;
    clp.setOption("verbosemodel", "noverbosemodel", &verboseModel, "use stacked verbose performance model");



    std::ostringstream galeriStream;
    std::string rhsFile, coordFile, coordMapFile, nullFile, materialFile;  // unused
    typedef typename Teuchos::ScalarTraits<SC>::magnitudeType real_type;
    typedef Xpetra::MultiVector<real_type, LO, GO, NO> RealValuedMultiVector;
    RCP<RealValuedMultiVector> coordinates;
    RCP<MultiVector> nullspace, material, x, b;
    RCP<Matrix> A;
    RCP<const Map> map;

    switch (clp.parse(argc, argv)) {
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED: return EXIT_SUCCESS;
      case Teuchos::CommandLineProcessor::PARSE_ERROR:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL: break;
    }

    // Load the matrix off disk (or generate it via Galeri), assuming only one right hand side is loaded.
    MatrixLoad<SC, LO, GO, NO>(comm, lib, binaryFormat, matrixFile, rhsFile, rowMapFile, colMapFile,
                               domainMapFile, rangeMapFile, coordFile, coordMapFile, nullFile, materialFile,
                               map, A, coordinates, nullspace, material, x, b, 1,
                               galeriParameters, xpetraParameters, galeriStream);


    out << "========================================================" << endl
        << xpetraParameters
        // << matrixParameters
        << "========================================================" << endl
        << "Template Types:" << endl
        << "  Scalar:        " << Teuchos::demangleName(typeid(SC).name()) << endl
        << "  LocalOrdinal:  " << Teuchos::demangleName(typeid(LO).name()) << endl
        << "  GlobalOrdinal: " << Teuchos::demangleName(typeid(GO).name()) << endl
        << "  Node:          " << Teuchos::demangleName(typeid(NO).name()) << endl
        << "Sizes:" << endl
        << "  Scalar:        " << sizeof(SC) << endl
        << "  LocalOrdinal:  " << sizeof(LO) << endl
        << "  GlobalOrdinal: " << sizeof(GO) << endl
        << "========================================================" << endl
        << "Matrix:        " << Teuchos::demangleName(typeid(Matrix).name()) << endl
        << "Vector:        " << Teuchos::demangleName(typeid(MultiVector).name()) << endl
        << "Hierarchy:     " << Teuchos::demangleName(typeid(Hierarchy).name()) << endl
        << "========================================================" << endl
        << " MPI Ranks:    " << numProc << endl
        << "========================================================" << endl;

#if defined(HAVE_TPETRA_INST_OPENMP)
    out << "Tpetra::KokkosCompat::KokkosOpenMPWrapperNode::execution_space().concurrency() = " << Tpetra::KokkosCompat::KokkosOpenMPWrapperNode::execution_space().concurrency() << endl
        << "========================================================" << endl;
#endif

    // =========================================================================
    // Problem construction
    // =========================================================================

    RCP<MultiVector> y          = Xpetra::VectorFactory<SC, LO, GO, Node>::Build(A->getRowMap());
    RCP<MultiVector> y_baseline = Xpetra::VectorFactory<SC, LO, GO, Node>::Build(A->getRowMap());
    x->putScalar(Teuchos::ScalarTraits<Scalar>::one());
    y->putScalar(Teuchos::ScalarTraits<Scalar>::nan());
    y_baseline->putScalar(Teuchos::ScalarTraits<Scalar>::nan());

    typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> crs_matrix_type;
    typedef Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> vector_type;
    RCP<const crs_matrix_type> At;
    vector_type xt;
    vector_type yt;

    At                         = Utilities::Op2TpetraCrs(A);

}  // main

//- -- --------------------------------------------------------
#define MUELU_AUTOMATIC_TEST_ETI_NAME main_
#include "MueLu_Test_ETI.hpp"

int main(int argc, char* argv[]) {
  return Automatic_Test_ETI(argc, argv);
}
