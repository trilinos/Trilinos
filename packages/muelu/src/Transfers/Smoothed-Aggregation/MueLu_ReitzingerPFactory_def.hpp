// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_REITZINGERPFACTORY_DEF_HPP
#define MUELU_REITZINGERPFACTORY_DEF_HPP

#include <Xpetra_MapFactory.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_CrsMatrix.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MatrixMatrix.hpp>
#include <Xpetra_MultiVector.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Import.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
// #include <Xpetra_IO.hpp>

#include "MueLu_ReitzingerPFactory_decl.hpp"

#include <Teuchos_ScalarTraits.hpp>

#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_ImportUtils.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> ReitzingerPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());

#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
  SET_VALID_ENTRY("repartition: enable");
  SET_VALID_ENTRY("repartition: use subcommunicators");
  SET_VALID_ENTRY("tentative: calculate qr");
  SET_VALID_ENTRY("tentative: constant column sums");
#undef SET_VALID_ENTRY

  validParamList->set<RCP<const FactoryBase> >("D0", Teuchos::null, "Generating factory of the matrix D0");
  validParamList->set<RCP<const FactoryBase> >("NodeAggMatrix", Teuchos::null, "Generating factory of the matrix NodeAggMatrix");
  validParamList->set<RCP<const FactoryBase> >("Pnodal", Teuchos::null, "Generating factory of the matrix P");

  // Make sure we don't recursively validate options for the matrixmatrix kernels
  ParameterList norecurse;
  norecurse.disableRecursiveValidation();
  validParamList->set<ParameterList>("matrixmatrix: kernel params", norecurse, "MatrixMatrix kernel parameters");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ReitzingerPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& fineLevel, Level& coarseLevel) const {
  Input(fineLevel, "D0");
  Input(coarseLevel, "NodeAggMatrix");
  Input(coarseLevel, "Pnodal");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ReitzingerPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& fineLevel, Level& coarseLevel) const {
  return BuildP(fineLevel, coarseLevel);
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ReitzingerPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildP(Level& fineLevel, Level& coarseLevel) const {
  FactoryMonitor m(*this, "Build", coarseLevel);

  using XMM               = Xpetra::MatrixMatrix<SC, LO, GO, NO>;
  using local_matrix_type = typename Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::local_matrix_type;
  using rowptr_type       = typename local_matrix_type::row_map_type::non_const_type;
  using colidx_type       = typename local_matrix_type::index_type::non_const_type;
  using values_type       = typename local_matrix_type::values_type::non_const_type;

  using impl_scalar_type = typename Matrix::impl_scalar_type;
  using ATS              = KokkosKernels::ArithTraits<impl_scalar_type>;
  using mag_type         = typename KokkosKernels::ArithTraits<impl_scalar_type>::magnitudeType;
  using magATS           = KokkosKernels::ArithTraits<mag_type>;

  using execution_space = typename Node::execution_space;
  using memory_space    = typename Node::memory_space;

  const auto one_Scalar      = Teuchos::ScalarTraits<Scalar>::one();
  const auto one_impl_scalar = ATS::one();
  const auto zero_LO         = KokkosKernels::ArithTraits<LocalOrdinal>::zero();
  const auto one_LO          = KokkosKernels::ArithTraits<LocalOrdinal>::one();
  const auto one_mag         = magATS::one();
  const auto eps_mag         = magATS::epsilon();
  const auto INVALID_GO      = Teuchos::OrdinalTraits<GlobalOrdinal>::invalid();

  // Using a nodal prolongator Pn and the discrete gradient matrix D0, this factory constructs
  // a coarse discrete gradient matrix D0H and an edge prolongator Pe such that the commuting
  // relationship
  //
  //  D0 * Pn = Pe * D0H
  //
  // holds.

  // The construction of the coarse discrete gradient works as follows.
  // We create edges between aggregates that contain at least one pair of connected nodes.
  // This boils down to computing the matrix
  //
  //  Z := (D0*Pn)^T * (D0 * Pn).
  //
  // If Z_ij != 0 we create an edge e between the nodal aggregates i and j.
  // Z is clearly symmetric. We only create a single edge between i and j.
  // In the distributed case, we also need to decide which rank owns the coarse edge e.
  // If both endpoints i and j live on process proc0 then proc0 should obiously own the edge e.
  // If i lives on proc0 and j lives on proc1, we tie-break based on the rule
  //
  //  min{proc0, proc1} if proc0+proc1 is odd,
  //  max{proc0, proc1} if proc0+proc1 is even.
  //
  // The orientation of the edge (encoded in the values 1 and -1) is determined by the GIDs of the endpoints.
  // All edges point from smaller GID to larger GID, i.e. i < j.

  // We also perform detection of boundary condtions and add additional edges to the coarse discrete gradient.
  // We detect all edges in D0 that only connect to a single node.

  Teuchos::FancyOStream& out0 = GetBlackHole();
  const ParameterList& pL     = GetParameterList();

  bool update_communicators = pL.get<bool>("repartition: enable") && pL.get<bool>("repartition: use subcommunicators");

  RCP<Matrix> D0 = Get<RCP<Matrix> >(fineLevel, "D0");
  RCP<Matrix> Pn = Get<RCP<Matrix> >(coarseLevel, "Pnodal");

  // This needs to be an Operator because if NodeMatrix gets repartitioned away, we get an Operator on the level
  RCP<Operator> CoarseNodeMatrix = Get<RCP<Operator> >(coarseLevel, "NodeAggMatrix");

  // Matrix matrix params
  RCP<ParameterList> mm_params = rcp(new ParameterList);
  if (pL.isSublist("matrixmatrix: kernel params"))
    mm_params->sublist("matrixmatrix: kernel params") = pL.sublist("matrixmatrix: kernel params");

  {  // Check that Pn is piecewise constant
    // TODO: Should this be a debug-only check?

    auto vec_ones = VectorFactory::Build(Pn->getDomainMap(), false);
    vec_ones->putScalar(one_Scalar);
    auto vec_rowsums = VectorFactory::Build(Pn->getRangeMap(), false);
    Pn->apply(*vec_ones, *vec_rowsums, Teuchos::NO_TRANS);

    auto lclPn      = Pn->getLocalMatrixDevice();
    auto lclRowSums = vec_rowsums->getLocalViewDevice(Tpetra::Access::ReadOnly);

    bool all_entries_ok = true;
    Kokkos::parallel_reduce(
        Kokkos::RangePolicy<execution_space>(0, lclPn.numRows()), KOKKOS_LAMBDA(const LocalOrdinal rlid, bool& entries_ok) {
      // rowsums are 1
      entries_ok = entries_ok && (ATS::magnitude(lclRowSums(rlid, 0) - one_impl_scalar) < eps_mag);

      // all nonzero entries are 1
      auto row = lclPn.rowConst(rlid);
      for (LocalOrdinal k = 0; k < row.length; ++k) {
        entries_ok = entries_ok && (ATS::magnitude(row.value(k)-one_impl_scalar) < eps_mag);

      } }, Kokkos::LAnd<bool>(all_entries_ok));

    TEUCHOS_TEST_FOR_EXCEPTION(!all_entries_ok, std::runtime_error, "The prolongator needs to be piecewise constant and all entries need to be 1.");
  }

  RCP<Matrix> D0_Pn;
  RCP<Matrix> D0H;
  LocalOrdinal numCoarseEdges          = 0;
  LocalOrdinal numCoarseRegularEdges   = 0;
  LocalOrdinal numCoarseDirichletEdges = 0;
  auto isDirichletFineEdge             = Xpetra::VectorFactory<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>::Build(D0->getRowMap(), false);
  auto numFineEdges                    = isDirichletFineEdge->getMap()->getLocalNumElements();
  {
    // Construct D0*Pn and Z := (D0*Pn)^T * (D0*Pn)
    RCP<Matrix> dummy;
    D0_Pn         = XMM::Multiply(*D0, false, *Pn, false, dummy, GetOStream(Runtime0), true, true);
    RCP<Matrix> Z = XMM::Multiply(*D0_Pn, true, *D0_Pn, false, dummy, GetOStream(Runtime0), true, true);

    auto rowMap       = Z->getRowMap();
    auto colMap       = Z->getColMap();
    auto lclRowMap    = rowMap->getLocalMap();
    auto lclColMap    = colMap->getLocalMap();
    auto lclZ         = Z->getLocalMatrixDevice();
    auto numLocalRows = lclZ.numRows();

#ifdef HAVE_MUELU_DEBUG
    TEUCHOS_ASSERT(Utilities::MapsAreNested(*rowMap, *colMap));
#endif

    auto importer = Z->getCrsGraph()->getImporter();
    // TODO: replace with Kokkos once PID data lives on device
    Teuchos::Array<int> Z_col_pids;
    Kokkos::View<int*, memory_space> Z_col_pids_d;
    if (!importer.is_null()) {
      MueLu::ImportUtils<LO, GO, NO> utils;
      utils.getPids(*importer, Z_col_pids, false);
      Kokkos::View<int*, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged> > Z_col_pids_h(Z_col_pids.data(), Z_col_pids.size());
      Z_col_pids_d = Kokkos::View<int*, memory_space>("Z_col_pids_d", Z_col_pids.size());
      Kokkos::deep_copy(Z_col_pids_d, Z_col_pids_h);
    }

    int myProcId = rowMap->getComm()->getRank();

    // Tie-break criterion for owner of coarse edges
    auto tie_break = KOKKOS_LAMBDA(int proc0, int proc1) {
      if ((proc0 + proc1) % 2 == 1) {
        return Kokkos::min(proc0, proc1);
      } else {
        return Kokkos::max(proc0, proc1);
      }
    };

    // Utility function to determine whether we need to add a coarse edge for entry
    // (rlid, clid) of Z.
    auto add_edge = KOKKOS_LAMBDA(LocalOrdinal rlid, LocalOrdinal clid) {
      if (clid < numLocalRows) {
        // Both row and column index are local, this process owns the new edge.
        // Only create one edge between rlid and clid and ignore the transposed entry.
        return (rlid < clid);
      } else {
        // Column index is nonlocal. Need to decide if this process owns the new edge.
        int otherProcId = Z_col_pids_d(clid);
        int owner       = tie_break(myProcId, otherProcId);
        return (owner == myProcId);
      }
    };

    // Count up how many coarse regular edges we are creating.
    Kokkos::parallel_reduce(
        Kokkos::RangePolicy<execution_space>(0, numLocalRows), KOKKOS_LAMBDA(const LocalOrdinal rlid, LocalOrdinal& ne) {
          auto row = lclZ.rowConst(rlid);
          // Loop over entries in row of Z
          for (LocalOrdinal k = 0; k < row.length; ++k) {
            auto clid = row.colidx(k);
            if (add_edge(rlid, clid))
              ++ne;
          }
        },
        numCoarseRegularEdges);

    // Mark as singleParents any D0 edges with only one node (so these are
    // edges that connect an interior node with a Dirichlet node).
    // isDirichletFineEdge is 1 for edges with a single endpoint and 0 otherwise.
    using LOMatrix = Tpetra::CrsMatrix<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>;
    RCP<LOMatrix> D0_LocalOrdinal;
    {
      // We want to do a apply using D0 on vectors with Scalar=LocalOrdinal.
      // Something like this could work and would not require any memory allocations, but it requires
      // "convert" to be ETI'd for all possible scalar types.

      // toTpetra(D0)->template convert<LocalOrdinal>()->apply(*toTpetra(oneVec), *toTpetra(isDirichletFineEdge), Teuchos::NO_TRANS);

      using lo_local_matrix_type = typename LOMatrix::local_matrix_device_type;

      auto lclGraph = D0->getCrsGraph()->getLocalGraphDevice();
      Kokkos::View<LocalOrdinal*, memory_space> values(Kokkos::ViewAllocateWithoutInitializing("values_LocalOrdinal"), D0->getLocalNumEntries());
      {
        auto values_scalar = D0->getLocalMatrixDevice().values;
        Kokkos::parallel_for(
            "MueLu::ReitzingerPFactory::convert", Kokkos::RangePolicy<execution_space>(0, values.extent(0)), KOKKOS_LAMBDA(const size_t i) {
              if (values_scalar(i) == one_impl_scalar)
                values(i) = one_LO;
              else if (values_scalar(i) == -one_impl_scalar)
                values(i) = -one_LO;
              else
                Kokkos::abort("D0 contains bad values");
            });
      }
      lo_local_matrix_type lclMatrix("D0_LocalOrdinal", D0->getLocalMatrixDevice().numCols(), values, lclGraph);

      D0_LocalOrdinal = rcp(new LOMatrix(lclMatrix, toTpetra(D0->getRowMap()), toTpetra(D0->getColMap()), toTpetra(D0->getDomainMap()), toTpetra(D0->getRangeMap())));
    }
    {
      auto oneVec = Xpetra::VectorFactory<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>::Build(D0->getDomainMap(), false);
      oneVec->putScalar(KokkosKernels::ArithTraits<LocalOrdinal>::one());
      D0_LocalOrdinal->apply(*toTpetra(oneVec), *toTpetra(isDirichletFineEdge), Teuchos::NO_TRANS);

      auto lcl_isDirichletFineEdge = isDirichletFineEdge->getLocalViewDevice(Tpetra::Access::ReadWrite);
      auto lcl_D0_Pn               = D0_Pn->getLocalMatrixDevice();
      Kokkos::parallel_for(
          Kokkos::RangePolicy<execution_space>(0, lcl_isDirichletFineEdge.extent(0)), KOKKOS_LAMBDA(const LocalOrdinal i) {
            if (ATS::magnitude(ATS::magnitude(lcl_isDirichletFineEdge(i, 0)) - one_mag) > eps_mag) {
              // This is a regular edge with two end points.
              lcl_isDirichletFineEdge(i, 0) = zero_LO;
            } else {
              // This is an edge with one endpoint.
              lcl_isDirichletFineEdge(i, 0) = one_LO;
            }
          });
    }

    // Count the number of fine Dirichlet edges that are connected to every coarse nodal aggregate via
    // the graph of D0*Pn.
    auto numberConnectedFineDirichletEdgesToCoarseNode = Xpetra::VectorFactory<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node>::Build(D0_Pn->getDomainMap(), false);
    {
      auto abs_D0_Pn = LOMatrix(toTpetra(D0_Pn->getCrsGraph()));
      abs_D0_Pn.fillComplete(toTpetra(D0_Pn->getDomainMap()), toTpetra(D0_Pn->getRangeMap()));
      abs_D0_Pn.setAllToScalar(KokkosKernels::ArithTraits<LocalOrdinal>::one());
      abs_D0_Pn.apply(*toTpetra(isDirichletFineEdge), *toTpetra(numberConnectedFineDirichletEdgesToCoarseNode), Teuchos::TRANS);
    }

    // Count local Dirichlet coarse edges
    {
      auto lcl_numberConnectedFineDirichletEdgesToCoarseNode = numberConnectedFineDirichletEdgesToCoarseNode->getLocalViewDevice(Tpetra::Access::ReadOnly);

      Kokkos::parallel_reduce(
          Kokkos::RangePolicy<execution_space>(0, lcl_numberConnectedFineDirichletEdgesToCoarseNode.extent(0)),
          KOKKOS_LAMBDA(const LocalOrdinal i, LocalOrdinal& ne) {
            if (ATS::magnitude(lcl_numberConnectedFineDirichletEdgesToCoarseNode(i, 0)) > eps_mag) {
              ++ne;
            }
          },
          numCoarseDirichletEdges);
    }

    if (IsPrint(Statistics0)) {
      LocalOrdinal numGlobalRegularEdges;
      LocalOrdinal numGlobalDirichletEdges;
      MueLu_sumAll(rowMap->getComm(), numCoarseRegularEdges, numGlobalRegularEdges);
      MueLu_sumAll(rowMap->getComm(), numCoarseDirichletEdges, numGlobalDirichletEdges);
      GetOStream(Statistics0) << "regular edges: " << numGlobalRegularEdges << ", Dirichlet edges: " << numGlobalDirichletEdges << std::endl;
    }

    numCoarseEdges = numCoarseRegularEdges + numCoarseDirichletEdges;
    rowptr_type rowptr(Kokkos::ViewAllocateWithoutInitializing("rowptr D0H"), numCoarseEdges + 1);
    // 2 entries per regular edge, 1 entry per Dirichlet edge
    LocalOrdinal nnz = 2 * numCoarseRegularEdges + numCoarseDirichletEdges;
    colidx_type colidx(Kokkos::ViewAllocateWithoutInitializing("colidx D0H"), nnz);
    values_type values(Kokkos::ViewAllocateWithoutInitializing("values D0H"), nnz);

    // Fill regular edges
    Kokkos::parallel_scan(
        Kokkos::RangePolicy<execution_space>(0, numLocalRows),
        KOKKOS_LAMBDA(const LocalOrdinal rlid, LocalOrdinal& ne, const bool update) {
          auto row = lclZ.rowConst(rlid);
          if (!update) {
            // First pass: figure out offsets for entries.
            for (LocalOrdinal k = 0; k < row.length; ++k) {
              auto clid = row.colidx(k);
              if (add_edge(rlid, clid))
                ++ne;
            }
          } else {
            // Second pass: enter entries

            // initialize
            if (rlid == 0)
              rowptr(rlid) = 0;

            auto rgid  = lclRowMap.getGlobalElement(rlid);
            auto rclid = lclColMap.getLocalElement(rgid);

            // loop over entries in row of Z
            for (LocalOrdinal k = 0; k < row.length; ++k) {
              auto clid = row.colidx(k);
              if (add_edge(rlid, clid)) {
                auto cgid = lclColMap.getGlobalElement(clid);
                // enter the two end-points of the edge, orient edge based on GIDs of nodal endpoints
                colidx(2 * ne)     = rclid;
                colidx(2 * ne + 1) = clid;
                if (rgid < cgid) {
                  values(2 * ne)     = -one_impl_scalar;
                  values(2 * ne + 1) = one_impl_scalar;
                } else {
                  values(2 * ne)     = one_impl_scalar;
                  values(2 * ne + 1) = -one_impl_scalar;
                }
                ++ne;
                rowptr(ne) = 2 * ne;
              }
            }
          }
        });

    // Fill Dirichlet edges
    // Create one coarse Dirichlet edge for every nodal aggregate that is connected to at least one fine Dirichlet edge.
    {
      auto lcl_numberConnectedFineDirichletEdgesToCoarseNode = numberConnectedFineDirichletEdgesToCoarseNode->getLocalViewDevice(Tpetra::Access::ReadOnly);
      Kokkos::parallel_scan(
          Kokkos::RangePolicy<execution_space>(0, lcl_numberConnectedFineDirichletEdgesToCoarseNode.extent(0)),
          KOKKOS_LAMBDA(const LocalOrdinal agg_lid, LocalOrdinal& ne, const bool update) {
            if (ATS::magnitude(lcl_numberConnectedFineDirichletEdgesToCoarseNode(agg_lid, 0)) > eps_mag) {
              if (!update) {
                // First pass: figure out offsets
                ++ne;
              } else {
                // Second pass: fill
                colidx(2 * numCoarseRegularEdges + ne) = agg_lid;
                values(2 * numCoarseRegularEdges + ne) = one_impl_scalar;
                ++ne;
                rowptr(numCoarseRegularEdges + ne) = 2 * numCoarseRegularEdges + ne;
              }
            }
          });
    }

    auto D0H_rowmap = MapFactory::Build(rowMap->lib(), INVALID_GO, numCoarseEdges, 0, rowMap->getComm());
    auto lclD0H     = local_matrix_type("D0H", numCoarseEdges, colMap->getLocalNumElements(), nnz, values, rowptr, colidx);

    // Construct distributed matrix
    D0H = MatrixFactory::Build(lclD0H, D0H_rowmap, colMap, Z->getDomainMap(), D0H_rowmap);
  }

  // Create the Pe matrix, but with the extra entries.  From ML's notes:
  /* The general idea is that the matrix                              */
  /*                        T_h P_n T_H^*                             */
  /* is almost Pe. If we make sure that P_n contains 1's and -1's, the*/
  /* matrix triple product will yield a matrix with +/- 1 and +/- 2's.*/
  /* If we remove all the 1's and divide the 2's by 2. we arrive at Pe*/

  RCP<Matrix> D0_Pn_D0HT;
  {
    SubFactoryMonitor m2(*this, "Generate Pe (pre-fix)", coarseLevel);
#if 0
    {
      // If you're concerned about processor / rank mismatches, this debugging code might help
      int rank =  D0->getRowMap()->getComm()->getRank();
      int fine_level = fineLevel.GetLevelID();
      printf("[%d] Level %d Checkpoint #2 Pn = %d/%d/%d/%d D0c = %d/%d/%d/%d D0 = %d/%d/%d/%d\n",rank,fine_level,
             Pn->getRangeMap()->getComm()->getSize(),
             Pn->getRowMap()->getComm()->getSize(),
             Pn->getColMap()->getComm()->getSize(),
             Pn->getDomainMap()->getComm()->getSize(),
             D0H->getRangeMap()->getComm()->getSize(),
             D0H->getRowMap()->getComm()->getSize(),
             D0H->getColMap()->getComm()->getSize(),
             D0H->getDomainMap()->getComm()->getSize(),
             D0->getRangeMap()->getComm()->getSize(),
             D0->getRowMap()->getComm()->getSize(),
             D0->getColMap()->getComm()->getSize(),
             D0->getDomainMap()->getComm()->getSize());
      fflush(stdout);
      D0->getRowMap()->getComm()->barrier();
    }
#endif
    RCP<Matrix> dummy;
    RCP<Matrix> Pn_D0cT = XMM::Multiply(*Pn, false, *D0H, true, dummy, out0, true, true, "Pn*D0c'", mm_params);

    // We don't want this guy getting accidently used later
    if (!mm_params.is_null()) mm_params->remove("importer", false);

    D0_Pn_D0HT = XMM::Multiply(*D0, false, *Pn_D0cT, false, dummy, out0, true, true, "D0*(Pn*D0c')", mm_params);

    // TODO: Something like this *might* work.  But this specifically, doesn't
    // Pe = XMM::Multiply(*D0_Pn_nonghosted,false,*D0H,true,dummy,out0,true,true,"(D0*Pn)*D0c'",mm_params);
  }

  RCP<Matrix> Pe;
  {
    auto lcl_D0_Pn               = D0_Pn->getLocalMatrixDevice();
    auto lcl_D0_Pn_D0HT          = D0_Pn_D0HT->getLocalMatrixDevice();
    auto lcl_isDirichletFineEdge = isDirichletFineEdge->getLocalViewDevice(Tpetra::Access::ReadOnly);

    auto lcl_colmap_D0_Pn_D0HT = D0_Pn_D0HT->getColMap()->getLocalMap();

    const auto half = one_impl_scalar / (one_impl_scalar + one_impl_scalar);

    // overallocate by 1 to allow for easier counting
    rowptr_type Pe_rowptr("Pe_rowptr", numFineEdges + 2);

    // count entries per row
    Kokkos::parallel_for(
        "Pe_count_entries", Kokkos::RangePolicy<execution_space>(0, numFineEdges), KOKKOS_LAMBDA(const LocalOrdinal fineEdge) {
          if (lcl_isDirichletFineEdge(fineEdge, 0) != one_LO) {
            // regular fine edge
            auto row = lcl_D0_Pn_D0HT.rowConst(fineEdge);
            for (int k = 0; k < row.length; ++k) {
              auto val = row.value(k);
              // filter out entries +-1 and 0
              if (!((ATS::magnitude(val - one_impl_scalar) < eps_mag) || (ATS::magnitude(val + one_impl_scalar) < eps_mag) || (ATS::magnitude(val) < eps_mag))) {
                // add entry (fineEdge, clid) -> val/2.
                ++Pe_rowptr(fineEdge + 2);
              }
            }
          } else {
            // Dirichlet interior fine edge
            ++Pe_rowptr(fineEdge + 2);
          }
        });

    // prefix sum
    LocalOrdinal Pe_nnz;
    Kokkos::parallel_scan(
        "Pe_prefix_sum", Kokkos::RangePolicy<execution_space>(0, numFineEdges + 2), KOKKOS_LAMBDA(const LocalOrdinal rlid, LocalOrdinal& nnz, const bool update) {
          nnz += Pe_rowptr(rlid);
          if (update) {
            Pe_rowptr(rlid) = nnz;
          }
        },
        Pe_nnz);

    // allocate view for indices and values
    colidx_type Pe_colidx("Pe_colidx", Pe_nnz);
    values_type Pe_values("Pe_values", Pe_nnz);

    // We build the mapping from coarse nodes lids wrt column map of D0_Pn to coarse edge gids.
    RCP<GOVector> map_coarseNodes_colMap_D0_Pn_to_coarseEdges;
    {
      auto map_coarseEdges_rowMap_D0H_to_coarseEdges = Xpetra::VectorFactory<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node>::Build(D0H->getRowMap());
      {
        auto lcl_map_coarseEdges_rowMap_D0H_to_coarseEdges = map_coarseEdges_rowMap_D0H_to_coarseEdges->getLocalViewDevice(Tpetra::Access::OverwriteAll);
        auto lclMap                                        = D0H->getRowMap()->getLocalMap();
        Kokkos::parallel_for(
            Kokkos::RangePolicy<execution_space>(numCoarseRegularEdges, numCoarseEdges), KOKKOS_LAMBDA(const LocalOrdinal coarseEdge) {
              lcl_map_coarseEdges_rowMap_D0H_to_coarseEdges(coarseEdge, 0) = lclMap.getGlobalElement(coarseEdge);
            });
      }
      auto map_coarseNodes_domainMap_D0_Pn_to_coarseEdges = Xpetra::VectorFactory<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node>::Build(D0H->getDomainMap());
      {
        // We want to do a transpose apply using D0H on vectors with Scalar=GlobalOrdinal.
        // Something like this could work and would not require any memory allocations, but it requires
        // "convert" to be ETI'd for all possible scalar types.

        // toTpetra(D0H)->template convert<GlobalOrdinal>()->apply(*toTpetra(map_coarseEdges_rowMap_D0H_to_coarseEdges), *toTpetra(map_coarseNodes_domainMap_D0_Pn_to_coarseEdges), Teuchos::TRANS);

        using GOMatrix             = Tpetra::CrsMatrix<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node>;
        using go_local_matrix_type = typename GOMatrix::local_matrix_device_type;

        auto lclGraph = D0H->getCrsGraph()->getLocalGraphDevice();
        typename go_local_matrix_type::values_type::non_const_type ones("ones_GlobalOrdinal", D0H->getLocalNumEntries());
        const auto one_GO = KokkosKernels::ArithTraits<typename go_local_matrix_type::values_type::value_type>::one();
        Kokkos::deep_copy(ones, one_GO);

        go_local_matrix_type lclMatrix("D0H_GlobalOrdinal", D0H->getLocalMatrixDevice().numCols(), ones, lclGraph);

        auto D0H_GlobalOrdinal = GOMatrix(lclMatrix, toTpetra(D0H->getRowMap()), toTpetra(D0H->getColMap()), toTpetra(D0H->getDomainMap()), toTpetra(D0H->getRangeMap()));
        D0H_GlobalOrdinal.apply(*toTpetra(map_coarseEdges_rowMap_D0H_to_coarseEdges), *toTpetra(map_coarseNodes_domainMap_D0_Pn_to_coarseEdges), Teuchos::TRANS);
      }

      auto importer = D0_Pn->getCrsGraph()->getImporter();
      if (!importer.is_null()) {
        map_coarseNodes_colMap_D0_Pn_to_coarseEdges = Xpetra::VectorFactory<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node>::Build(D0_Pn->getColMap());
        map_coarseNodes_colMap_D0_Pn_to_coarseEdges->doImport(*map_coarseNodes_domainMap_D0_Pn_to_coarseEdges, *importer, Xpetra::INSERT);
      } else {
        map_coarseNodes_colMap_D0_Pn_to_coarseEdges = map_coarseNodes_domainMap_D0_Pn_to_coarseEdges;
      }
    }
    {
      auto lcl_map_coarseNodes_colMap_D0_Pn_to_coarseEdges = map_coarseNodes_colMap_D0_Pn_to_coarseEdges->getLocalViewDevice(Tpetra::Access::ReadOnly);

      // fill
      Kokkos::parallel_for(
          "Pe_fill", Kokkos::RangePolicy<execution_space>(0, numFineEdges), KOKKOS_LAMBDA(const LocalOrdinal fineEdge_lid) {
            if (lcl_isDirichletFineEdge(fineEdge_lid, 0) != one_LO) {
              // regular fine edge
              auto row = lcl_D0_Pn_D0HT.rowConst(fineEdge_lid);
              for (int k = 0; k < row.length; ++k) {
                auto val = row.value(k);
                if (!((ATS::magnitude(val - one_impl_scalar) < eps_mag) || (ATS::magnitude(val + one_impl_scalar) < eps_mag) || (ATS::magnitude(val) < eps_mag))) {
                  auto clid = row.colidx(k);
                  // add entry (fineEdge_lid, clid) -> val/2.
                  auto offset       = Pe_rowptr(fineEdge_lid + 1);
                  Pe_colidx(offset) = clid;
                  Pe_values(offset) = val * half;
                  ++Pe_rowptr(fineEdge_lid + 1);
                }
              }
            } else {
              // Dirichlet interior fine edge
              // Only one nonzero entry in row of D0_Pn: (fineEdge_lid, coarseNode_lid_D0_Pn) -> val_D0_Pn
              for (auto offset_D0_Pn = lcl_D0_Pn.graph.row_map(fineEdge_lid); offset_D0_Pn < lcl_D0_Pn.graph.row_map(fineEdge_lid + 1); ++offset_D0_Pn) {
                LocalOrdinal coarseNode_lid_D0_Pn = lcl_D0_Pn.graph.entries(offset_D0_Pn);
                impl_scalar_type val_D0_Pn        = lcl_D0_Pn.values(offset_D0_Pn);
                if (ATS::magnitude(val_D0_Pn) > eps_mag) {
                  GlobalOrdinal coarseEdge_gid = lcl_map_coarseNodes_colMap_D0_Pn_to_coarseEdges(coarseNode_lid_D0_Pn, 0);

                  auto coarseEdge_lid_D0_Pn_D0HT = lcl_colmap_D0_Pn_D0HT.getLocalElement(coarseEdge_gid);

                  // We rely on the fact that all coarse interior edges have been created with value 1.
                  const auto val_D0H = one_impl_scalar;

                  // add entry (fineEdge_lid, coarseEdge) -> val_D0_Pn/val_D0H to edge prolongator
                  auto offset_Pe       = Pe_rowptr(fineEdge_lid + 1);
                  Pe_colidx(offset_Pe) = coarseEdge_lid_D0_Pn_D0HT;
                  Pe_values(offset_Pe) = val_D0_Pn / val_D0H;
                  ++Pe_rowptr(fineEdge_lid + 1);
                  break;
                }
              }
            }
          });
    }
    auto lclPe = local_matrix_type("Pe", numFineEdges, D0_Pn_D0HT->getColMap()->getLocalNumElements(), Pe_nnz, Pe_values, Kokkos::subview(Pe_rowptr, Kokkos::make_pair((decltype(numFineEdges))0, numFineEdges + 1)), Pe_colidx);

    // Construct distributed matrix
    Pe = MatrixFactory::Build(lclPe, D0->getRowMap(), D0_Pn_D0HT->getColMap(), D0H->getRangeMap(), D0->getRangeMap());
  }

  /* Check commuting property */
  CheckCommutingProperty(*Pe, *D0H, *D0, *Pn);

  /*  If we're repartitioning here, we need to cut down the communicators */
  // NOTE: We need to do this *after* checking the commuting property, since
  // that's going to need to fineLevel's communicators, not the repartitioned ones
  if (update_communicators) {
    // NOTE: We can only do D0 here.  We have to do Ke_coarse=(Re Ke_fine Pe) in RebalanceAcFactory
    RCP<const Teuchos::Comm<int> > newComm;
    if (!CoarseNodeMatrix.is_null()) newComm = CoarseNodeMatrix->getDomainMap()->getComm();
    RCP<const Map> newMap = MapFactory::copyMapWithNewComm(D0H->getRowMap(), newComm);
    D0H->removeEmptyProcessesInPlace(newMap);

    // The "in place" still leaves a dummy matrix here.  That needs to go
    if (newMap.is_null()) D0H = Teuchos::null;

    Set(coarseLevel, "InPlaceMap", newMap);
  }
  /* Set output on the level */
  Set(coarseLevel, "P", Pe);
  Set(coarseLevel, "Ptent", Pe);

  Set(coarseLevel, "D0", D0H);

  /* This needs to be kept for the smoothers */
  coarseLevel.Set("D0", D0H, NoFactory::get());
  coarseLevel.AddKeepFlag("D0", NoFactory::get(), MueLu::Final);
  coarseLevel.RemoveKeepFlag("D0", NoFactory::get(), MueLu::UserData);

#if 0
  {
    int numProcs = Pe->getRowMap()->getComm()->getSize();
    char fname[80];

    sprintf(fname, "Pe_%d_%d.mat", numProcs, fineLevel.GetLevelID());
    Xpetra::IO<SC, LO, GO, NO>::Write(fname, *Pe);
    sprintf(fname, "Pn_%d_%d.mat", numProcs, fineLevel.GetLevelID());
    Xpetra::IO<SC, LO, GO, NO>::Write(fname, *Pn);
    if (!D0H.is_null()) {
      sprintf(fname, "D0c_%d_%d.mat", numProcs, fineLevel.GetLevelID());
      Xpetra::IO<SC, LO, GO, NO>::Write(fname, *D0H);
    }
    sprintf(fname, "D0f_%d_%d.mat", numProcs, fineLevel.GetLevelID());
    Xpetra::IO<SC, LO, GO, NO>::Write(fname, *D0);
  }
#endif

}  // end Build

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ReitzingerPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
    CheckCommutingProperty(const Matrix& Pe, const Matrix& D0_c, const Matrix& D0_f, const Matrix& Pn) const {
  if (IsPrint(Statistics0)) {
    using XMM = MatrixMatrix;
    auto one  = Teuchos::ScalarTraits<SC>::one();

    RCP<Matrix> dummy;
    RCP<Matrix> left  = XMM::Multiply(Pe, false, D0_c, false, dummy, GetOStream(Runtime0));
    RCP<Matrix> right = XMM::Multiply(D0_f, false, Pn, false, dummy, GetOStream(Runtime0));

    RCP<Matrix> summation;
    XMM::TwoMatrixAdd(*left, false, one, *right, false, -one, summation, GetOStream(Runtime0));
    summation->fillComplete(left->getDomainMap(), left->getRangeMap());

    auto norm = summation->getFrobeniusNorm();
    GetOStream(Statistics0) << "CheckCommutingProperty: || Pe D0_c - D0_f Pn || = " << norm << std::endl;
  }

}  // end CheckCommutingProperty

}  // namespace MueLu

#define MUELU_REITZINGERPFACTORY_SHORT
#endif  // MUELU_REITZINGERPFACTORY_DEF_HPP
