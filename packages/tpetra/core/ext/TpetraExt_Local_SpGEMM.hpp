// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRAEXT_LOCAL_SPGEMM_HPP
#define TPETRAEXT_LOCAL_SPGEMM_HPP

#include "KokkosSparse_spgemm.hpp"

namespace Tpetra::MMdetails {

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node,
          class LocalOrdinalViewType>
void mult_A_B_newmatrix_local_serial(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
                                     CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Bview,
                                     const LocalOrdinalViewType& targetMapToOrigRow_dev,
                                     const LocalOrdinalViewType& targetMapToImportRow_dev,
                                     const LocalOrdinalViewType& Bcol2Ccol_dev,
                                     const LocalOrdinalViewType& Icol2Ccol_dev,
                                     CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C,
                                     Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Node>> Cimport,
                                     const std::string& label,
                                     const Teuchos::RCP<Teuchos::ParameterList>& params) {
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;

  Tpetra::Details::ProfilingRegion MM("TpetraExt: MMM: Newmatrix Core");

  // Lots and lots of typedefs
  typedef typename Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_host_type KCRS;
  typedef typename KCRS::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::const_type c_lno_view_t;
  typedef typename graph_t::row_map_type::non_const_type lno_view_t;
  typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
  typedef typename KCRS::values_type::non_const_type scalar_view_t;

  typedef Scalar SC;
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;
  typedef Node NO;
  typedef Map<LO, GO, NO> map_type;
  const size_t ST_INVALID = Teuchos::OrdinalTraits<LO>::invalid();
  const LO LO_INVALID     = Teuchos::OrdinalTraits<LO>::invalid();
  const SC SC_ZERO        = Teuchos::ScalarTraits<Scalar>::zero();

  auto targetMapToOrigRow =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),
                                          targetMapToOrigRow_dev);
  auto targetMapToImportRow =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),
                                          targetMapToImportRow_dev);
  auto Bcol2Ccol =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),
                                          Bcol2Ccol_dev);
  auto Icol2Ccol =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),
                                          Icol2Ccol_dev);

  // Sizes
  RCP<const map_type> Ccolmap = C.getColMap();
  size_t m                    = Aview.origMatrix->getLocalNumRows();
  size_t n                    = Ccolmap->getLocalNumElements();
  size_t b_max_nnz_per_row    = Bview.origMatrix->getLocalMaxNumRowEntries();

  // Grab the  Kokkos::SparseCrsMatrices & inner stuff
  const KCRS Amat = Aview.origMatrix->getLocalMatrixHost();
  const KCRS Bmat = Bview.origMatrix->getLocalMatrixHost();

  c_lno_view_t Arowptr = Amat.graph.row_map, Browptr = Bmat.graph.row_map;
  const lno_nnz_view_t Acolind = Amat.graph.entries, Bcolind = Bmat.graph.entries;
  const scalar_view_t Avals = Amat.values, Bvals = Bmat.values;

  c_lno_view_t Irowptr;
  lno_nnz_view_t Icolind;
  scalar_view_t Ivals;
  if (!Bview.importMatrix.is_null()) {
    auto lclB         = Bview.importMatrix->getLocalMatrixHost();
    Irowptr           = lclB.graph.row_map;
    Icolind           = lclB.graph.entries;
    Ivals             = lclB.values;
    b_max_nnz_per_row = std::max(b_max_nnz_per_row, Bview.importMatrix->getLocalMaxNumRowEntries());
  }

  // Classic csr assembly (low memory edition)
  //
  // mfh 27 Sep 2016: C_estimate_nnz does not promise an upper bound.
  // The method loops over rows of A, and may resize after processing
  // each row.  Chris Siefert says that this reflects experience in
  // ML; for the non-threaded case, ML found it faster to spend less
  // effort on estimation and risk an occasional reallocation.
  size_t CSR_alloc = std::max(C_estimate_nnz(*Aview.origMatrix, *Bview.origMatrix), n);
  lno_view_t Crowptr(Kokkos::ViewAllocateWithoutInitializing("Crowptr"), m + 1);
  lno_nnz_view_t Ccolind(Kokkos::ViewAllocateWithoutInitializing("Ccolind"), CSR_alloc);
  scalar_view_t Cvals(Kokkos::ViewAllocateWithoutInitializing("Cvals"), CSR_alloc);

  // mfh 27 Sep 2016: The c_status array is an implementation detail
  // of the local sparse matrix-matrix multiply routine.

  // The status array will contain the index into colind where this entry was last deposited.
  //   c_status[i] <  CSR_ip - not in the row yet
  //   c_status[i] >= CSR_ip - this is the entry where you can find the data
  // We start with this filled with INVALID's indicating that there are no entries yet.
  // Sadly, this complicates the code due to the fact that size_t's are unsigned.
  size_t INVALID = Teuchos::OrdinalTraits<size_t>::invalid();
  std::vector<size_t> c_status(n, ST_INVALID);

  // mfh 27 Sep 2016: Here is the local sparse matrix-matrix multiply
  // routine.  The routine computes C := A * (B_local + B_remote).
  //
  // For column index Aik in row i of A, targetMapToOrigRow[Aik] tells
  // you whether the corresponding row of B belongs to B_local
  // ("orig") or B_remote ("Import").

  // For each row of A/C
  size_t CSR_ip = 0, OLD_ip = 0;
  for (size_t i = 0; i < m; i++) {
    // mfh 27 Sep 2016: m is the number of rows in the input matrix A
    // on the calling process.
    Crowptr[i] = CSR_ip;

    // mfh 27 Sep 2016: For each entry of A in the current row of A
    for (size_t k = Arowptr[i]; k < Arowptr[i + 1]; k++) {
      LO Aik        = Acolind[k];  // local column index of current entry of A
      const SC Aval = Avals[k];    // value of current entry of A
      if (Aval == SC_ZERO)
        continue;  // skip explicitly stored zero values in A

      if (targetMapToOrigRow[Aik] != LO_INVALID) {
        // mfh 27 Sep 2016: If the entry of targetMapToOrigRow
        // corresponding to the current entry of A is populated, then
        // the corresponding row of B is in B_local (i.e., it lives on
        // the calling process).

        // Local matrix
        size_t Bk = static_cast<size_t>(targetMapToOrigRow[Aik]);

        // mfh 27 Sep 2016: Go through all entries in that row of B_local.
        for (size_t j = Browptr[Bk]; j < Browptr[Bk + 1]; ++j) {
          LO Bkj = Bcolind[j];
          LO Cij = Bcol2Ccol[Bkj];

          if (c_status[Cij] == INVALID || c_status[Cij] < OLD_ip) {
            // New entry
            c_status[Cij]   = CSR_ip;
            Ccolind[CSR_ip] = Cij;
            Cvals[CSR_ip]   = Aval * Bvals[j];
            CSR_ip++;

          } else {
            Cvals[c_status[Cij]] += Aval * Bvals[j];
          }
        }

      } else {
        // mfh 27 Sep 2016: If the entry of targetMapToOrigRow
        // corresponding to the current entry of A NOT populated (has
        // a flag "invalid" value), then the corresponding row of B is
        // in B_local (i.e., it lives on the calling process).

        // Remote matrix
        size_t Ik = static_cast<size_t>(targetMapToImportRow[Aik]);
        for (size_t j = Irowptr[Ik]; j < Irowptr[Ik + 1]; ++j) {
          LO Ikj = Icolind[j];
          LO Cij = Icol2Ccol[Ikj];

          if (c_status[Cij] == INVALID || c_status[Cij] < OLD_ip) {
            // New entry
            c_status[Cij]   = CSR_ip;
            Ccolind[CSR_ip] = Cij;
            Cvals[CSR_ip]   = Aval * Ivals[j];
            CSR_ip++;
          } else {
            Cvals[c_status[Cij]] += Aval * Ivals[j];
          }
        }
      }
    }

    // Resize for next pass if needed
    if (i + 1 < m && CSR_ip + std::min(n, (Arowptr[i + 2] - Arowptr[i + 1]) * b_max_nnz_per_row) > CSR_alloc) {
      CSR_alloc *= 2;
      Kokkos::resize(Ccolind, CSR_alloc);
      Kokkos::resize(Cvals, CSR_alloc);
    }
    OLD_ip = CSR_ip;
  }

  Crowptr[m] = CSR_ip;

  // Downward resize
  Kokkos::resize(Ccolind, CSR_ip);
  Kokkos::resize(Cvals, CSR_ip);

  {
    Tpetra::Details::ProfilingRegion MM3("TpetraExt: MMM: Newmatrix Final Sort");

    // Final sort & set of CRS arrays
    if (params.is_null() || params->get("sort entries", Details::MatrixTraits<Scalar, LocalOrdinal, GlobalOrdinal, Node>::spgemmNeedsSortedInputs()))
      Import_Util::sortCrsEntries(Crowptr, Ccolind, Cvals);
    C.setAllValues(Crowptr, Ccolind, Cvals);
  }

  Tpetra::Details::ProfilingRegion MM4("TpetraExt: MMM: Newmatrix ESCC");
  {
    // Final FillComplete
    //
    // mfh 27 Sep 2016: So-called "expert static fill complete" bypasses
    // Import (from domain Map to column Map) construction (which costs
    // lots of communication) by taking the previously constructed
    // Import object.  We should be able to do this without interfering
    // with the implementation of the local part of sparse matrix-matrix
    // multply above.
    RCP<Teuchos::ParameterList> labelList = rcp(new Teuchos::ParameterList);
    labelList->set("Timer Label", label);
    if (!params.is_null()) labelList->set("compute global constants", params->get("compute global constants", true));
    RCP<const Export<LO, GO, NO>> dummyExport;
    C.expertStaticFillComplete(Bview.origMatrix->getDomainMap(), Aview.origMatrix->getRangeMap(), Cimport, dummyExport, labelList);
  }
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node,
          class LocalOrdinalViewType>
void mult_A_B_newmatrix_local_KokkosKernels(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
                                            CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Bview,
                                            const LocalOrdinalViewType& Acol2Brow,
                                            const LocalOrdinalViewType& Acol2Irow,
                                            const LocalOrdinalViewType& Bcol2Ccol,
                                            const LocalOrdinalViewType& Icol2Ccol,
                                            CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C,
                                            Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Node>> Cimport,
                                            const std::string& label,
                                            const Teuchos::RCP<Teuchos::ParameterList>& params) {
  using Teuchos::RCP;
  using Teuchos::rcp;
  RCP<Tpetra::Details::ProfilingRegion> MM = rcp(new Tpetra::Details::ProfilingRegion("TpetraExt: MMM: Newmatrix"));

  // Node-specific code
  std::string myalg("SPGEMM_KK_MEMORY");
  int team_work_size;
  std::string nodeLowerCase;
  std::string nodeUpperCase;

  // This code works for serial and openmp backends as well.
  // Taking out the static asserts is all that is needed.
#ifdef KOKKOS_ENABLE_SERIAL
  static_assert(!std::is_same_v<Node, Tpetra::KokkosCompat::KokkosSerialWrapperNode>,
                "This is the default implementation of mult_A_B_newmatrix_local_KokkosKernels. We should never get here but instead use the overload for serial node.");
  if constexpr (std::is_same_v<Node, Tpetra::KokkosCompat::KokkosSerialWrapperNode>) {
    nodeLowerCase  = "serial";
    nodeUpperCase  = "Serial";
    team_work_size = 1;
  }
#endif
#ifdef KOKKOS_ENABLE_OPENMP
  static_assert(!std::is_same_v<Node, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>,
                "This is the default implementation of mult_A_B_newmatrix_local_KokkosKernels. We should never get here but instead use the overload for openmp node.");
  if constexpr (std::is_same_v<Node, Tpetra::KokkosCompat::KokkosOpenMPWrapperNode>) {
    nodeLowerCase  = "openmp";
    nodeUpperCase  = "OpenMP";
    team_work_size = 16;
  }
#endif
#ifdef KOKKOS_ENABLE_CUDA
  if constexpr (std::is_same_v<Node, Tpetra::KokkosCompat::KokkosCudaWrapperNode>) {
    nodeLowerCase  = "cuda";
    nodeUpperCase  = "Cuda";
    team_work_size = 16;
  }
#endif
#ifdef KOKKOS_ENABLE_HIP
  if constexpr (std::is_same_v<Node, Tpetra::KokkosCompat::KokkosHipWrapperNode>) {
    nodeLowerCase  = "hip";
    nodeUpperCase  = "HIP";
    team_work_size = 16;
  }
#endif
#ifdef KOKKOS_ENABLE_SYCL
  if constexpr (std::is_same_v<Node, Tpetra::KokkosCompat::KokkosSYCLWrapperNode>) {
    nodeLowerCase  = "sycl";
    nodeUpperCase  = "SYCL";
    team_work_size = 16;
  }
#endif

  // Lots and lots of typedefs
  using KCRS           = typename Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_device_type;
  using device_t       = typename KCRS::device_type;
  using graph_t        = typename KCRS::StaticCrsGraphType;
  using lno_view_t     = typename graph_t::row_map_type::non_const_type;
  using int_view_t     = Kokkos::View<int*, typename lno_view_t::array_layout, typename lno_view_t::memory_space, typename lno_view_t::memory_traits>;
  using c_lno_view_t   = typename graph_t::row_map_type::const_type;
  using lno_nnz_view_t = typename graph_t::entries_type::non_const_type;
  using scalar_view_t  = typename KCRS::values_type::non_const_type;
  // typedef typename graph_t::row_map_type::const_type lno_view_t_const;

  // Options
  if (!params.is_null()) {
    if (params->isParameter(nodeLowerCase + ": algorithm"))
      myalg = params->get(nodeLowerCase + ": algorithm", myalg);
    if (params->isParameter(nodeLowerCase + ": team work size"))
      team_work_size = params->get(nodeLowerCase + ": team work size", team_work_size);
  }

  // KokkosKernelsHandle
  using KernelHandle = KokkosKernels::Experimental::KokkosKernelsHandle<
      typename lno_view_t::const_value_type, typename lno_nnz_view_t::const_value_type, typename scalar_view_t::const_value_type,
      typename device_t::execution_space, typename device_t::memory_space, typename device_t::memory_space>;
  using IntKernelHandle = KokkosKernels::Experimental::KokkosKernelsHandle<
      typename int_view_t::const_value_type, typename lno_nnz_view_t::const_value_type, typename scalar_view_t::const_value_type,
      typename device_t::execution_space, typename device_t::memory_space, typename device_t::memory_space>;

  // Grab the  Kokkos::SparseCrsMatrices
  const KCRS Amat = Aview.origMatrix->getLocalMatrixDevice();
  const KCRS Bmat = Bview.origMatrix->getLocalMatrixDevice();

  c_lno_view_t Arowptr         = Amat.graph.row_map;
  c_lno_view_t Browptr         = Bmat.graph.row_map;
  const lno_nnz_view_t Acolind = Amat.graph.entries;
  const lno_nnz_view_t Bcolind = Bmat.graph.entries;
  const scalar_view_t Avals    = Amat.values;
  const scalar_view_t Bvals    = Bmat.values;

  // Get the algorithm mode
  std::string alg = nodeUpperCase + std::string(" algorithm");
  //  printf("DEBUG: Using kernel: %s\n",myalg.c_str());
  if (!params.is_null() && params->isParameter(alg)) myalg = params->get(alg, myalg);
  KokkosSparse::SPGEMMAlgorithm alg_enum = KokkosSparse::StringToSPGEMMAlgorithm(myalg);

  // Merge the B and Bimport matrices
  KCRS Bmerged = Tpetra::MMdetails::merge_matrices(Aview, Bview, Acol2Brow, Acol2Irow, Bcol2Ccol, Icol2Ccol, C.getColMap()->getLocalNumElements());

  if constexpr (Details::MatrixTraits<Scalar, LocalOrdinal, GlobalOrdinal, Node>::spgemmNeedsSortedInputs()) {
    if (!KokkosSparse::isCrsGraphSorted(Bmerged.graph.row_map, Bmerged.graph.entries)) {
      KokkosSparse::sort_crs_matrix(Bmerged);
    }
  }

  MM = Teuchos::null;
  MM = rcp(new Tpetra::Details::ProfilingRegion("TpetraExt: MMM: Newmatrix Core"));

  // Do the multiply on whatever we've got
  typename KernelHandle::nnz_lno_t AnumRows = Amat.numRows();
  typename KernelHandle::nnz_lno_t BnumRows = Bmerged.numRows();
  typename KernelHandle::nnz_lno_t BnumCols = Bmerged.numCols();

  // regardless of whether integer row ptrs are used, need to ultimately produce the row pointer, entries, and values of the expected types
  lno_view_t row_mapC(Kokkos::ViewAllocateWithoutInitializing("non_const_lno_row"), AnumRows + 1);
  lno_nnz_view_t entriesC;
  scalar_view_t valuesC;

  Tpetra::Details::IntRowPtrHelper<decltype(Bmerged)> irph(Bmerged.nnz(), Bmerged.graph.row_map);
  const bool useIntRowptrs =
      irph.shouldUseIntRowptrs() &&
      Aview.origMatrix->getApplyHelper()->shouldUseIntRowptrs();

  if (useIntRowptrs) {
    IntKernelHandle kh;
    kh.create_spgemm_handle(alg_enum);
    kh.set_team_work_size(team_work_size);

    int_view_t int_row_mapC(Kokkos::ViewAllocateWithoutInitializing("non_const_int_row"), AnumRows + 1);

    auto Aint = Aview.origMatrix->getApplyHelper()->getIntRowptrMatrix(Amat);
    auto Bint = irph.getIntRowptrMatrix(Bmerged);

    {
      Tpetra::Details::ProfilingRegion MM2("TpetraExt: MMM: Newmatrix KokkosKernels symbolic int");
      KokkosSparse::spgemm_symbolic(&kh, AnumRows, BnumRows, BnumCols, Aint.graph.row_map, Aint.graph.entries, false, Bint.graph.row_map, Bint.graph.entries, false, int_row_mapC);
    }

    Tpetra::Details::ProfilingRegion MM2("TpetraExt: MMM: Newmatrix KokkosKernels numeric int");

    size_t c_nnz_size = kh.get_spgemm_handle()->get_c_nnz();
    if (c_nnz_size) {
      entriesC = lno_nnz_view_t(Kokkos::ViewAllocateWithoutInitializing("entriesC"), c_nnz_size);
      valuesC  = scalar_view_t(Kokkos::ViewAllocateWithoutInitializing("valuesC"), c_nnz_size);
    }
    KokkosSparse::spgemm_numeric(&kh, AnumRows, BnumRows, BnumCols, Aint.graph.row_map, Aint.graph.entries, Aint.values, false, Bint.graph.row_map, Bint.graph.entries, Bint.values, false, int_row_mapC, entriesC, valuesC);
    // transfer the integer rowptrs back to the correct rowptr type
    Kokkos::parallel_for(
        int_row_mapC.size(), KOKKOS_LAMBDA(int i) { row_mapC(i) = int_row_mapC(i); });
    kh.destroy_spgemm_handle();

  } else {
    KernelHandle kh;
    kh.create_spgemm_handle(alg_enum);
    kh.set_team_work_size(team_work_size);

    {
      Tpetra::Details::ProfilingRegion MM2("TpetraExt: Jacobi: Newmatrix KokkosKernels symbolic non-int");
      KokkosSparse::spgemm_symbolic(&kh, AnumRows, BnumRows, BnumCols, Amat.graph.row_map, Amat.graph.entries, false, Bmerged.graph.row_map, Bmerged.graph.entries, false, row_mapC);
    }

    Tpetra::Details::ProfilingRegion MM2("TpetraExt: Jacobi: Newmatrix KokkosKernels numeric non-int");
    size_t c_nnz_size = kh.get_spgemm_handle()->get_c_nnz();
    if (c_nnz_size) {
      entriesC = lno_nnz_view_t(Kokkos::ViewAllocateWithoutInitializing("entriesC"), c_nnz_size);
      valuesC  = scalar_view_t(Kokkos::ViewAllocateWithoutInitializing("valuesC"), c_nnz_size);
    }

    KokkosSparse::spgemm_numeric(&kh, AnumRows, BnumRows, BnumCols, Amat.graph.row_map, Amat.graph.entries, Amat.values, false, Bmerged.graph.row_map, Bmerged.graph.entries, Bmerged.values, false, row_mapC, entriesC, valuesC);

    kh.destroy_spgemm_handle();
  }

  MM = Teuchos::null;
  MM = rcp(new Tpetra::Details::ProfilingRegion("TpetraExt: MMM: Newmatrix Sort"));

  // Sort & set values
  if (params.is_null() || params->get("sort entries", true))
    Import_Util::sortCrsEntries(row_mapC, entriesC, valuesC);
  C.setAllValues(row_mapC, entriesC, valuesC);

  MM = Teuchos::null;
  MM = rcp(new Tpetra::Details::ProfilingRegion("TpetraExt: MMM: Newmatrix ESFC"));

  // Final Fillcomplete
  RCP<Teuchos::ParameterList> labelList = rcp(new Teuchos::ParameterList);
  labelList->set("Timer Label", label);
  if (!params.is_null()) labelList->set("compute global constants", params->get("compute global constants", true));
  RCP<const Export<LocalOrdinal, GlobalOrdinal, Node>> dummyExport;
  C.expertStaticFillComplete(Bview.origMatrix->getDomainMap(), Aview.origMatrix->getRangeMap(), Cimport, dummyExport, labelList);
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node,
          class LocalOrdinalViewType>
void mult_A_B_reuse_local_serial(CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Aview,
                                 CrsMatrixStruct<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Bview,
                                 const LocalOrdinalViewType& targetMapToOrigRow_dev,
                                 const LocalOrdinalViewType& targetMapToImportRow_dev,
                                 const LocalOrdinalViewType& Bcol2Ccol_dev,
                                 const LocalOrdinalViewType& Icol2Ccol_dev,
                                 CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>& C,
                                 Teuchos::RCP<const Import<LocalOrdinal, GlobalOrdinal, Node>> Cimport,
                                 const std::string& label,
                                 const Teuchos::RCP<Teuchos::ParameterList>& params) {
  using Teuchos::RCP;
  using Teuchos::rcp;

  // Lots and lots of typedefs
  typedef typename Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_host_type KCRS;
  typedef typename KCRS::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::const_type c_lno_view_t;
  typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
  typedef typename KCRS::values_type::non_const_type scalar_view_t;

  typedef Scalar SC;
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;
  typedef Node NO;
  typedef Map<LO, GO, NO> map_type;
  const size_t ST_INVALID = Teuchos::OrdinalTraits<LO>::invalid();
  const LO LO_INVALID     = Teuchos::OrdinalTraits<LO>::invalid();
  const SC SC_ZERO        = Teuchos::ScalarTraits<Scalar>::zero();

  Tpetra::Details::ProfilingRegion MM("TpetraExt: MMM: Reuse SerialCore");
  // Since this is being run on Cuda, we need to fence because the below code will use UVM
  // typename graph_t::execution_space().fence();

  // KDDKDD UVM Without UVM, need to copy targetMap arrays to host.
  // KDDKDD UVM Ideally, this function would run on device and use
  // KDDKDD UVM KokkosKernels instead of this host implementation.
  auto targetMapToOrigRow =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),
                                          targetMapToOrigRow_dev);
  auto targetMapToImportRow =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),
                                          targetMapToImportRow_dev);
  auto Bcol2Ccol =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),
                                          Bcol2Ccol_dev);
  auto Icol2Ccol =
      Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(),
                                          Icol2Ccol_dev);

  // Sizes
  RCP<const map_type> Ccolmap = C.getColMap();
  size_t m                    = Aview.origMatrix->getLocalNumRows();
  size_t n                    = Ccolmap->getLocalNumElements();

  // Grab the  Kokkos::SparseCrsMatrices & inner stuff
  const KCRS Amat = Aview.origMatrix->getLocalMatrixHost();
  const KCRS Bmat = Bview.origMatrix->getLocalMatrixHost();
  const KCRS Cmat = C.getLocalMatrixHost();

  c_lno_view_t Arowptr         = Amat.graph.row_map,
               Browptr         = Bmat.graph.row_map,
               Crowptr         = Cmat.graph.row_map;
  const lno_nnz_view_t Acolind = Amat.graph.entries,
                       Bcolind = Bmat.graph.entries,
                       Ccolind = Cmat.graph.entries;
  const scalar_view_t Avals = Amat.values, Bvals = Bmat.values;
  scalar_view_t Cvals = Cmat.values;

  c_lno_view_t Irowptr;
  lno_nnz_view_t Icolind;
  scalar_view_t Ivals;
  if (!Bview.importMatrix.is_null()) {
    auto lclB = Bview.importMatrix->getLocalMatrixHost();
    Irowptr   = lclB.graph.row_map;
    Icolind   = lclB.graph.entries;
    Ivals     = lclB.values;
  }

  // Classic csr assembly (low memory edition)
  // mfh 27 Sep 2016: The c_status array is an implementation detail
  // of the local sparse matrix-matrix multiply routine.

  // The status array will contain the index into colind where this entry was last deposited.
  //   c_status[i] <  CSR_ip - not in the row yet
  //   c_status[i] >= CSR_ip - this is the entry where you can find the data
  // We start with this filled with INVALID's indicating that there are no entries yet.
  // Sadly, this complicates the code due to the fact that size_t's are unsigned.
  std::vector<size_t> c_status(n, ST_INVALID);

  // For each row of A/C
  size_t CSR_ip = 0, OLD_ip = 0;
  for (size_t i = 0; i < m; i++) {
    // First fill the c_status array w/ locations where we're allowed to
    // generate nonzeros for this row
    OLD_ip = Crowptr[i];
    CSR_ip = Crowptr[i + 1];
    for (size_t k = OLD_ip; k < CSR_ip; k++) {
      c_status[Ccolind[k]] = k;

      // Reset values in the row of C
      Cvals[k] = SC_ZERO;
    }

    for (size_t k = Arowptr[i]; k < Arowptr[i + 1]; k++) {
      LO Aik        = Acolind[k];
      const SC Aval = Avals[k];
      if (Aval == SC_ZERO)
        continue;

      if (targetMapToOrigRow[Aik] != LO_INVALID) {
        // Local matrix
        size_t Bk = static_cast<size_t>(targetMapToOrigRow[Aik]);

        for (size_t j = Browptr[Bk]; j < Browptr[Bk + 1]; ++j) {
          LO Bkj = Bcolind[j];
          LO Cij = Bcol2Ccol[Bkj];

          TEUCHOS_TEST_FOR_EXCEPTION(c_status[Cij] < OLD_ip || c_status[Cij] >= CSR_ip,
                                     std::runtime_error, "Trying to insert a new entry (" << i << "," << Cij << ") into a static graph "
                                                                                          << "(c_status = " << c_status[Cij] << " of [" << OLD_ip << "," << CSR_ip << "))");

          Cvals[c_status[Cij]] += Aval * Bvals[j];
        }

      } else {
        // Remote matrix
        size_t Ik = static_cast<size_t>(targetMapToImportRow[Aik]);
        for (size_t j = Irowptr[Ik]; j < Irowptr[Ik + 1]; ++j) {
          LO Ikj = Icolind[j];
          LO Cij = Icol2Ccol[Ikj];

          TEUCHOS_TEST_FOR_EXCEPTION(c_status[Cij] < OLD_ip || c_status[Cij] >= CSR_ip,
                                     std::runtime_error, "Trying to insert a new entry (" << i << "," << Cij << ") into a static graph "
                                                                                          << "(c_status = " << c_status[Cij] << " of [" << OLD_ip << "," << CSR_ip << "))");

          Cvals[c_status[Cij]] += Aval * Ivals[j];
        }
      }
    }
  }

  {
    Tpetra::Details::ProfilingRegion MM3("TpetraExt: MMM: Reuse ESFC");
    C.fillComplete(C.getDomainMap(), C.getRangeMap());
  }
}

}  // namespace Tpetra::MMdetails

#endif
