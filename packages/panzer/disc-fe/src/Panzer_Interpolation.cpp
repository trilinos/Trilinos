// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Panzer_Interpolation.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_IntrepidOrientation.hpp"
#include "Intrepid2_OrientationTools.hpp"
#include "Intrepid2_LagrangianInterpolation.hpp"
#ifdef PANZER_HAVE_EPETRA_STACK
#include "Epetra_MpiComm.h"
#endif

// #define PANZER_INTERPOLATION_DEBUG_OUTPUT = 1

namespace panzer {

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
removeSmallEntries(Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >& A,
                   typename Teuchos::ScalarTraits<Scalar>::magnitudeType tol) {
  using crs_matrix   = Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using row_ptr_type = typename crs_matrix::local_graph_device_type::row_map_type::non_const_type;
  using col_idx_type = typename crs_matrix::local_graph_device_type::entries_type::non_const_type;
  using vals_type    = typename crs_matrix::local_matrix_device_type::values_type;

  using ATS = Kokkos::ArithTraits<Scalar>;
  using impl_SC  = typename ATS::val_type;
  using impl_ATS = Kokkos::ArithTraits<impl_SC>;

  auto lclA = A->getLocalMatrixDevice();

  auto rowptr = row_ptr_type("rowptr", lclA.numRows() + 1);

  Kokkos::parallel_for(
      "removeSmallEntries::rowptr1",
      Kokkos::RangePolicy<LocalOrdinal>(0, lclA.numRows()),
      KOKKOS_LAMBDA(const LocalOrdinal rlid) {
        auto row = lclA.row(rlid);
        for (LocalOrdinal k = 0; k < row.length; ++k) {
          if (impl_ATS::magnitude(row.value(k)) > tol) {
            rowptr(rlid + 1) += 1;
          }
        }
      });
  LocalOrdinal nnz;
  Kokkos::parallel_scan(
      "removeSmallEntries::rowptr2",
      Kokkos::RangePolicy<LocalOrdinal>(0, lclA.numRows()),
      KOKKOS_LAMBDA(const LocalOrdinal rlid, LocalOrdinal& partial_nnz, bool is_final) {
        partial_nnz += rowptr(rlid + 1);
        if (is_final)
          rowptr(rlid + 1) = partial_nnz;
      },
      nnz);

  auto idx  = col_idx_type("idx", nnz);
  auto vals = vals_type("vals", nnz);

  Kokkos::parallel_for(
      "removeSmallEntries::indicesValues",
      Kokkos::RangePolicy<LocalOrdinal>(0, lclA.numRows()),
      KOKKOS_LAMBDA(const LocalOrdinal rlid) {
        auto row = lclA.row(rlid);
        auto I   = rowptr(rlid);
        for (LocalOrdinal k = 0; k < row.length; ++k) {
          if (impl_ATS::magnitude(row.value(k)) > tol) {
            idx(I)  = row.colidx(k);
            vals(I) = row.value(k);
            I += 1;
          }
        }
      });
  Kokkos::fence();

  auto newA = Teuchos::rcp(new crs_matrix(A->getRowMap(), A->getColMap(), rowptr, idx, vals));
  newA->fillComplete(A->getDomainMap(),
                     A->getRangeMap());
  return newA;
}


Teuchos::RCP<Thyra::LinearOpBase<double>> buildInterpolation(const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits>> &linObjFactory,
                                                             const std::string &domain_basis_name, const std::string &range_basis_name,
                                                             Intrepid2::EOperator op, size_t worksetSize,
                                                             const bool matrixFree)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;

  using Scalar = double;

  using tpetraBlockedLinObjFactory = typename panzer::BlockedTpetraLinearObjFactory<panzer::Traits, Scalar, LocalOrdinal, GlobalOrdinal>;
#ifdef PANZER_HAVE_EPETRA_STACK
  using epetraBlockedLinObjFactory = typename panzer::BlockedEpetraLinearObjFactory<panzer::Traits, LocalOrdinal>;
#endif
  using UGI = panzer::GlobalIndexer;

  // must be able to cast to a block linear object factory
  RCP<const tpetraBlockedLinObjFactory > tblof = rcp_dynamic_cast<const tpetraBlockedLinObjFactory >(linObjFactory);
#ifdef PANZER_HAVE_EPETRA_STACK
  RCP<const epetraBlockedLinObjFactory > eblof = rcp_dynamic_cast<const epetraBlockedLinObjFactory >(linObjFactory);
#endif

  RCP<const panzer::BlockedDOFManager> blockedDOFMngr;
  if (tblof != Teuchos::null) {
    blockedDOFMngr = tblof->getGlobalIndexer();
#ifdef PANZER_HAVE_EPETRA_STACK
  } else if (eblof != Teuchos::null) {
    blockedDOFMngr = eblof->getGlobalIndexer();
#endif
  } else {
    TEUCHOS_ASSERT(false);
  }

  // get global indexers for domain and range dofs
  std::vector<RCP<UGI> > fieldDOFMngrs = blockedDOFMngr->getFieldDOFManagers();
  int domainFieldNum = blockedDOFMngr->getFieldNum(domain_basis_name);
  int rangeFieldNum = blockedDOFMngr->getFieldNum(range_basis_name);
  int domainBlockIndex = blockedDOFMngr->getFieldBlock(domainFieldNum);
  int rangeBlockIndex = blockedDOFMngr->getFieldBlock(rangeFieldNum);
  RCP<panzer::DOFManager> domain_ugi = rcp_dynamic_cast<panzer::DOFManager>(blockedDOFMngr->getFieldDOFManagers()[domainBlockIndex], true);
  RCP<panzer::DOFManager> range_ugi = rcp_dynamic_cast<panzer::DOFManager>(blockedDOFMngr->getFieldDOFManagers()[rangeBlockIndex], true);

  RCP<const panzer::ConnManager> conn = blockedDOFMngr->getConnManager();

  return buildInterpolation(conn, domain_ugi, range_ugi, domain_basis_name, range_basis_name, op, worksetSize,
                            /*forceVectorial=*/false,
                            /*useTpetra=*/tblof != Teuchos::null,
                            matrixFree);
}


Teuchos::RCP<Thyra::LinearOpBase<double> > buildInterpolation(const Teuchos::RCP<const panzer::ConnManager> &conn,
                                                              const Teuchos::RCP<panzer::DOFManager> &domain_ugi,
                                                              const Teuchos::RCP<panzer::DOFManager> &range_ugi,
                                                              const std::string& domain_basis_name,
                                                              const std::string& range_basis_name,
                                                              Intrepid2::EOperator op,
                                                              size_t worksetSize,
                                                              const bool force_vectorial,
                                                              const bool useTpetra,
                                                              const bool matrixFree)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;

  using Scalar = double;

  using STS = Teuchos::ScalarTraits<Scalar>;
  using KAT = Kokkos::ArithTraits<Scalar>;
  using OT  = Teuchos::OrdinalTraits<GlobalOrdinal>;

  using DeviceSpace = PHX::Device;
  using HostSpace = Kokkos::HostSpace;
  using ots = Intrepid2::OrientationTools<DeviceSpace>;
  using li = Intrepid2::LagrangianInterpolation<DeviceSpace>;
  using DynRankDeviceView = Kokkos::DynRankView<double, DeviceSpace>;

  using tp_graph = Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal>;
  using tp_matrix = Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal>;
  using tp_map = Tpetra::Map<LocalOrdinal, GlobalOrdinal>;
#ifdef PANZER_HAVE_EPETRA_STACK
  using ep_linObjContainer = panzer::BlockedEpetraLinearObjContainer;
  using ep_matrix = Epetra_CrsMatrix;
  using ep_map = Epetra_Map;
#endif

  if (matrixFree) {
    TEUCHOS_ASSERT(useTpetra);
    TEUCHOS_ASSERT(!force_vectorial);
    auto mfOp = rcp(new MatrixFreeInterpolationOp<Scalar,LocalOrdinal,GlobalOrdinal>(conn, domain_ugi, range_ugi, domain_basis_name, range_basis_name, op, worksetSize));
    return Thyra::tpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,typename tp_matrix::node_type>(Thyra::createVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal>(mfOp->getRangeMap()),
                                                                                                  Thyra::createVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal>(mfOp->getDomainMap()),
                                                                                                  mfOp);
  }

  // get the domain and range bases
  auto domain_fieldPattern = domain_ugi->getFieldPattern(domain_basis_name);
  auto domain_basis = rcp_dynamic_cast<const panzer::Intrepid2FieldPattern>(domain_fieldPattern,true)->getIntrepidBasis();
  auto range_fieldPattern = range_ugi->getFieldPattern(range_basis_name);
  auto range_basis = rcp_dynamic_cast<const panzer::Intrepid2FieldPattern>(range_fieldPattern,true)->getIntrepidBasis();

  // cardinalities
  const size_t domainCardinality = domain_basis->getCardinality();
  const size_t rangeCardinality = range_basis->getCardinality();

  const int dim = range_basis->getBaseCellTopology().getDimension();

  if (op == Intrepid2::OPERATOR_VALUE) {
    TEUCHOS_ASSERT(rangeCardinality >= domainCardinality);
    TEUCHOS_ASSERT_EQUALITY(domain_basis->getFunctionSpace(), range_basis->getFunctionSpace());
  }

  // Create the global interp matrix.
  RCP<const tp_map> tp_rangemap;
  RCP<const tp_map> tp_domainmap;
  RCP<const tp_map> tp_rowmap;
  RCP<const tp_map> tp_colmap;
  RCP<tp_matrix> tp_interp_matrix;
  typename tp_matrix::local_matrix_device_type lcl_tp_interp_matrix;
#ifdef PANZER_HAVE_EPETRA_STACK
  RCP<const ep_map> ep_rangemap;
  RCP<const ep_map> ep_domainmap;
  RCP<const ep_map> ep_rowmap;
  RCP<const ep_map> ep_colmap;
  RCP<ep_matrix> ep_interp_matrix;
#endif

  auto rangeElementLIDs_d = range_ugi->getLIDs();
  auto domainElementLIDs_d = domain_ugi->getLIDs();

  RCP<Thyra::LinearOpBase<Scalar> > thyra_interp;
  size_t maxNumElementsPerBlock = 0;
  LocalOrdinal minLocalIndex = 0;
  LocalOrdinal maxLocalIndex = 0;
  if (useTpetra) {
    // build maps
    std::vector<GlobalOrdinal> gids;
    range_ugi->getOwnedIndices(gids);
    tp_rowmap = rcp(new tp_map(OT::invalid(), gids.data(), static_cast<LocalOrdinal>(gids.size()), OT::zero(), range_ugi->getComm()));
    tp_rangemap = tp_rowmap;
    domain_ugi->getOwnedIndices(gids);
    tp_domainmap = rcp(new tp_map(OT::invalid(), gids.data(), static_cast<LocalOrdinal>(gids.size()), OT::zero(), domain_ugi->getComm()));
    domain_ugi->getOwnedAndGhostedIndices(gids);
    tp_colmap = rcp(new tp_map(OT::invalid(), gids.data(), static_cast<LocalOrdinal>(gids.size()), OT::zero(), domain_ugi->getComm()));

    minLocalIndex = tp_rowmap->getMinLocalIndex();
    maxLocalIndex = tp_rowmap->getMaxLocalIndex();

    // estimate number of entries per row
    // This is an upper bound, as we are counting dofs that are on shared nodes, edges, faces more than once.
    using dv = Kokkos::DualView<size_t*, typename tp_graph::device_type>;
    dv numEntriesPerRow("numEntriesPerRow", tp_rowmap->getLocalNumElements());
    {
      auto numEntriesPerRow_d = numEntriesPerRow.view_device();

      // loop over element blocks
      std::vector<std::string> elementBlockIds;
      range_ugi->getElementBlockIds(elementBlockIds);
      for(std::size_t blockIter = 0; blockIter < elementBlockIds.size(); ++blockIter) {

        // loop over elements
        std::vector<int> elementIds = range_ugi->getElementBlock(elementBlockIds[blockIter]);
        Kokkos::View<int *, HostSpace> elementIds_h(elementIds.data(), elementIds.size());
        Kokkos::View<int *, DeviceSpace> elementIds_d("elementIds_d", elementIds_h.extent(0));
        Kokkos::deep_copy(elementIds_d, elementIds_h);
        maxNumElementsPerBlock = std::max(maxNumElementsPerBlock, elementIds.size());

        Kokkos::parallel_for("MiniEM_Interpolation::numEntriesPerRow",
                             Kokkos::RangePolicy<size_t, typename tp_matrix::node_type::execution_space>(0, elementIds.size()),
                             KOKKOS_LAMBDA(const size_t elemIter) {
                               auto elemId = elementIds_d(elemIter);

                               // get IDs for range dofs
                               auto rangeLIDs_d = Kokkos::subview(rangeElementLIDs_d, elemId, Kokkos::ALL());

                               // loop over range LIDs
                               for(size_t rangeIter = 0; rangeIter < rangeCardinality; ++rangeIter) {
                                 const LocalOrdinal range_row = rangeLIDs_d(rangeIter);
                                 const bool isOwned = ((minLocalIndex <= range_row) && (range_row <= maxLocalIndex));
                                 if (isOwned)
				   Kokkos::atomic_add(&numEntriesPerRow_d(range_row), domainCardinality);
                               } //end range LID loop
                             });
      } // blocks loop
      numEntriesPerRow.template modify<typename dv::t_dev>();
      numEntriesPerRow.template sync<typename dv::t_host>();
    }

    // Set up graph
    auto tp_interp_graph = rcp(new tp_graph(tp_rowmap, tp_colmap, numEntriesPerRow));

    { // This runs on host
      Kokkos::View<LocalOrdinal**, HostSpace> rangeElementLIDs_h("rangeElementLIDs_h", rangeElementLIDs_d.extent(0), rangeCardinality);
      Kokkos::View<LocalOrdinal**, HostSpace> domainElementLIDs_h("domainElementLIDs_h", domainElementLIDs_d.extent(0), domainCardinality);
      Kokkos::deep_copy(rangeElementLIDs_h, rangeElementLIDs_d);
      Kokkos::deep_copy(domainElementLIDs_h, domainElementLIDs_d);

      // loop over element blocks
      std::vector<std::string> elementBlockIds;
      range_ugi->getElementBlockIds(elementBlockIds);
      for(std::size_t blockIter = 0; blockIter < elementBlockIds.size(); ++blockIter) {

        // loop over elements
        std::vector<int> elementIds = range_ugi->getElementBlock(elementBlockIds[blockIter]);
        maxNumElementsPerBlock = std::max(maxNumElementsPerBlock, elementIds.size());
        for(std::size_t elemIter = 0; elemIter < elementIds.size(); ++elemIter) {
          auto elemId = elementIds[elemIter];

          // get IDs for range dofs
          auto rangeLIDs_h = Kokkos::subview(rangeElementLIDs_h, elemId, Kokkos::ALL());
          auto domainLIDs_h = Kokkos::subview(domainElementLIDs_h, elemId, Kokkos::ALL());

          // loop over range LIDs
          for(size_t rangeIter = 0; rangeIter < rangeCardinality; ++rangeIter) {
            const LocalOrdinal range_row = rangeLIDs_h(rangeIter);
            const bool isOwned = ((minLocalIndex <= range_row) && (range_row <= maxLocalIndex));
            if (isOwned) {
              Teuchos::ArrayView<LocalOrdinal> domainLIDs_av = Teuchos::ArrayView<LocalOrdinal>(domainLIDs_h.data(), domainLIDs_h.extent_int(0));
              tp_interp_graph->insertLocalIndices(range_row, domainLIDs_av);
            }
          } //end range LID loop
        } // elements loop
      } // blocks loop
    }

    RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList());
    pl->set("Optimize Storage", true);
    tp_interp_graph->fillComplete(tp_domainmap, tp_rangemap, pl);

    tp_interp_matrix = rcp(new tp_matrix(tp_interp_graph));
    lcl_tp_interp_matrix = tp_interp_matrix->getLocalMatrixDevice();
  }
#ifdef PANZER_HAVE_EPETRA_STACK
  else {

    const RCP<const Teuchos::MpiComm<int> > mpiComm = Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >(range_ugi->getComm());
    auto ep_comm = Teuchos::rcp(new Epetra_MpiComm(*mpiComm->getRawMpiComm()));
    std::vector<GlobalOrdinal> gids;
    range_ugi->getOwnedIndices(gids);
    ep_rowmap = rcp(new ep_map(OT::invalid(), static_cast<LocalOrdinal>(gids.size()), gids.data(), OT::zero(), *ep_comm));
    ep_rangemap = ep_rowmap;
    domain_ugi->getOwnedIndices(gids);
    ep_domainmap = rcp(new ep_map(OT::invalid(), static_cast<LocalOrdinal>(gids.size()), gids.data(), OT::zero(), *ep_comm));
    domain_ugi->getOwnedAndGhostedIndices(gids);
    ep_colmap = rcp(new ep_map(OT::invalid(), static_cast<LocalOrdinal>(gids.size()), gids.data(), OT::zero(), *ep_comm));

    {
      // loop over element blocks
      std::vector<std::string> elementBlockIds;
      range_ugi->getElementBlockIds(elementBlockIds);
      for(std::size_t blockIter = 0; blockIter < elementBlockIds.size(); ++blockIter) {

        // loop over elements
        std::vector<int> elementIds = range_ugi->getElementBlock(elementBlockIds[blockIter]);
        maxNumElementsPerBlock = std::max(maxNumElementsPerBlock, elementIds.size());
      }
    }

    // TODO: Fix this.
    size_t nnzPerRowEstimate = 25*domainCardinality;

    ep_interp_matrix = rcp(new ep_matrix(Copy, *ep_rowmap, *ep_colmap, static_cast<int>(nnzPerRowEstimate), /*StaticProfile=*/true));

    RCP<const Thyra::LinearOpBase<double> > th_ep_interp = Thyra::epetraLinearOp(ep_interp_matrix,
                                                                                 Thyra::NOTRANS,
                                                                                 Thyra::EPETRA_OP_APPLY_APPLY,
                                                                                 Thyra::EPETRA_OP_ADJOINT_SUPPORTED,
                                                                                 Thyra::create_VectorSpace(ep_rangemap),
                                                                                 Thyra::create_VectorSpace(ep_domainmap));
    thyra_interp = Teuchos::rcp_const_cast<Thyra::LinearOpBase<double> >(th_ep_interp);
  }
#endif

  // allocate some views
  int numCells;
  if (maxNumElementsPerBlock > 0)
    numCells = static_cast<int>(std::min(maxNumElementsPerBlock, worksetSize));
  else
    numCells = static_cast<int>(worksetSize);

  DynRankDeviceView range_dofCoords_d("range_dofCoords_d", rangeCardinality, dim);
  DynRankDeviceView basisCoeffsLIOriented_d("basisCoeffsLIOriented_d", numCells, rangeCardinality, domainCardinality);

  // the ranks of these depend on dimension
  DynRankDeviceView valuesAtDofCoordsNonOriented_d;
  DynRankDeviceView valuesAtDofCoordsOriented_d;

  if (!force_vectorial) {
    // Let Intrepid2 give us the correctly dimensioned view, then build one with +1 ranks and extent(0) == numCells
    auto temp = domain_basis->allocateOutputView(static_cast<int>(rangeCardinality), op);

    // These view have dimensions
    //  numCells, numFields=domainCardinality, numPoints=rangeCardinality, (spatialDim)
    //
    if (temp.rank() == 3) {
      valuesAtDofCoordsNonOriented_d = DynRankDeviceView("valuesAtDofCoordsNonOriented_d", temp.extent(0), temp.extent(1), temp.extent(2));
      valuesAtDofCoordsOriented_d    = DynRankDeviceView("valuesAtDofCoordsOriented_d", numCells, temp.extent(0), temp.extent(1), temp.extent(2));
    } else {
      valuesAtDofCoordsNonOriented_d = DynRankDeviceView("valuesAtDofCoordsNonOriented_d", temp.extent(0), temp.extent(1));
      valuesAtDofCoordsOriented_d    = DynRankDeviceView("valuesAtDofCoordsOriented_d", numCells, temp.extent(0), temp.extent(1));
    }
  } else {
    valuesAtDofCoordsNonOriented_d = DynRankDeviceView("valuesAtDofCoordsNonOriented_d", domainCardinality, rangeCardinality, dim);
    valuesAtDofCoordsOriented_d    = DynRankDeviceView("valuesAtDofCoordsOriented_d", numCells, domainCardinality, rangeCardinality, dim);
  }

  int fieldRank = Intrepid2::getFieldRank(range_basis->getFunctionSpace());
  TEUCHOS_ASSERT((fieldRank == 0) || (fieldRank == 1));

  auto entryFilterTol = 1e5*Teuchos::ScalarTraits<typename STS::magnitudeType>::eps();

  // range dof coordinates
  range_basis->getDofCoords(range_dofCoords_d);

  // compute values of op * (domain basis) at range dof coords on reference element
  domain_basis->getValues(valuesAtDofCoordsNonOriented_d, range_dofCoords_d, op);

  // get block ids
  std::vector<std::string> elementBlockIds;
  range_ugi->getElementBlockIds(elementBlockIds);

  // get orientations for all blocks
  std::map<std::string, std::vector<Intrepid2::Orientation> > orientations;
  buildIntrepidOrientations(elementBlockIds, *conn, orientations);

  // loop over element blocks
  for(std::size_t blockIter = 0; blockIter < elementBlockIds.size(); ++blockIter) {

    auto rangeOffsets_d = range_ugi->getGIDFieldOffsetsKokkos(elementBlockIds[blockIter], 0);
    auto domainOffsets_d = domain_ugi->getGIDFieldOffsetsKokkos(elementBlockIds[blockIter], 0);
#ifdef PANZER_INTERPOLATION_DEBUG_OUTPUT
    std::cout << "rangeOffsets_d" << std::endl;
    for (int i = 0; i < rangeCardinality; i++)
      std::cout << rangeOffsets_d(i) << " ";
    std::cout << std::endl;
    std::cout << "domainOffsets_d" << std::endl;
    for (int i = 0; i < domainCardinality; i++)
      std::cout << domainOffsets_d(i) << " ";
    std::cout << std::endl;
#endif

    // get element ids
    Kokkos::View<int *, DeviceSpace> elementIds_d;
    {
      std::vector<int> elementIds = range_ugi->getElementBlock(elementBlockIds[blockIter]);
      Kokkos::View<int *, HostSpace> elementIds_h(elementIds.data(), elementIds.size());
      elementIds_d = Kokkos::View<int *, DeviceSpace>("elementIds_d", elementIds_h.extent(0));
      Kokkos::deep_copy(elementIds_d, elementIds_h);
    }

    // get element orientations
    typename Kokkos::DynRankView<Intrepid2::Orientation,DeviceSpace> elemOrts_d ("elemOrts_d", elementIds_d.extent(0));
    {
      // copy orientations to device
      auto blockOrientations = orientations[elementBlockIds[blockIter]];
      Kokkos::View<Intrepid2::Orientation*, HostSpace> elemOrts_h(blockOrientations.data(), blockOrientations.size());
      Kokkos::deep_copy(elemOrts_d, elemOrts_h);
    }

    // loop over element worksets
    for(std::size_t elemIter = 0; elemIter < elementIds_d.extent(0); elemIter += numCells) {

      int endCellRange =
        std::min(numCells, Teuchos::as<int>(elementIds_d.extent(0)) -
                           Teuchos::as<int>(elemIter));

      // get subviews on workset
      auto ortsRange = Kokkos::make_pair(elemIter, std::min(elemIter + numCells, elemOrts_d.extent(0)));
      auto elemOrtsWorkset_d = Kokkos::subview(elemOrts_d, ortsRange);
      // Last workset might be shorter.
      auto worksetRange = Kokkos::make_pair(0, endCellRange);
      auto valuesAtDofCoordsOrientedWorkset_d = Kokkos::subview(valuesAtDofCoordsOriented_d, worksetRange, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
      auto basisCoeffsLIOrientedWorkset_d     = Kokkos::subview(basisCoeffsLIOriented_d,     worksetRange, Kokkos::ALL(), Kokkos::ALL());

      // apply orientations for domain basis
      // shuffles things in the second dimension, i.e. wrt domain basis
      ots::modifyBasisByOrientation(valuesAtDofCoordsOrientedWorkset_d,
                                    valuesAtDofCoordsNonOriented_d,
                                    elemOrtsWorkset_d,
                                    domain_basis.get());

      // get basis coefficients of domain basis functions wrt range basis
      for(size_t domainIter=0; domainIter<domainCardinality; domainIter++)
        // Get basis coeffs wrt range basis on reference element.
        // basisCoeffsLI has dimensions (numCells, numFields=rangeCardinality, domainCardinality)
        li::getBasisCoeffs(Kokkos::subview(basisCoeffsLIOrientedWorkset_d, Kokkos::ALL(), Kokkos::ALL(), domainIter),
                           Kokkos::subview(valuesAtDofCoordsOrientedWorkset_d, Kokkos::ALL(), domainIter, Kokkos::ALL(), Kokkos::ALL()),
                           range_basis.get(),
                           elemOrtsWorkset_d);

#ifdef PANZER_HAVE_EPETRA_STACK
      if (!useTpetra) { // Epetra fill

	Kokkos::View<LocalOrdinal*,DeviceSpace> indices_d("indices", domainCardinality);
	Kokkos::View<Scalar*,      DeviceSpace> values_d ("values",  domainCardinality);


        for (int cellNo = 0; cellNo < endCellRange; cellNo++) {
          auto elemId = elementIds_d(elemIter+cellNo);

          // get IDs for range and domain dofs
          auto rangeLIDs_d = Kokkos::subview(rangeElementLIDs_d, elemId, Kokkos::ALL());
          auto domainLIDs_d = Kokkos::subview(domainElementLIDs_d, elemId, Kokkos::ALL());

          // loop over range LIDs
          for(size_t rangeIter = 0; rangeIter < rangeCardinality; ++rangeIter) {
            const LocalOrdinal range_row = rangeLIDs_d(rangeOffsets_d(rangeIter));
            const bool isOwned = ep_rowmap->MyLID(range_row);
            if (isOwned) {
              // filter entries for zeros
              LocalOrdinal rowNNZ = 0;
              for(size_t domainIter = 0; domainIter < domainCardinality; domainIter++) {
                Scalar val = basisCoeffsLIOriented_d(cellNo, rangeIter, domainIter);
                if (KAT::magnitude(val) > entryFilterTol) {
                  indices_d(rowNNZ) = domainLIDs_d(domainOffsets_d(domainIter));
                  values_d(rowNNZ) = val;
                  ++rowNNZ;
                }
              }

              int ret = ep_interp_matrix->ReplaceMyValues(range_row, rowNNZ, values_d.data(), indices_d.data());
              if (ret != 0) {
                ret = ep_interp_matrix->InsertMyValues(range_row, rowNNZ, values_d.data(), indices_d.data());
                TEUCHOS_ASSERT(ret == 0);
              }
            } //end if owned
          } // end range LID loop
        } //end workset loop
      } // Epetra fill
      else
#endif
      { // Tpetra fill
        Kokkos::parallel_for(
                             "MiniEM_Interpolation::worksetLoop",
                             Kokkos::RangePolicy<int, typename tp_matrix::node_type::execution_space>(0, endCellRange),
                             KOKKOS_LAMBDA(const int cellNo) {
                               auto elemId = elementIds_d(elemIter+cellNo);

                               // get IDs for range and domain dofs
                               auto rangeLIDs_d = Kokkos::subview(rangeElementLIDs_d, elemId, Kokkos::ALL());
                               auto domainLIDs_d = Kokkos::subview(domainElementLIDs_d, elemId, Kokkos::ALL());

#ifdef PANZER_INTERPOLATION_DEBUG_OUTPUT
                               std::cout << "\n" << elemOrts_d(elemIter+cellNo).to_string() << std::endl;
                               std::cout << "rangeLIDs" << std::endl;
                               for (int i = 0; i < rangeCardinality; i++)
                                 std::cout << rangeLIDs_d(i) << " ";
                               std::cout << std::endl << "domainLIDs" << std::endl;
                               for (int i = 0; i < domainCardinality; i++)
                                 std::cout << domainLIDs_d(i) << " ";
                               std::cout << std::endl;
#endif
                               // loop over range LIDs
                               for(size_t rangeIter = 0; rangeIter < rangeCardinality; ++rangeIter) {
                                 const LocalOrdinal range_row = rangeLIDs_d(rangeOffsets_d(rangeIter));
                                 const bool isOwned = ((minLocalIndex <= range_row) && (range_row <= maxLocalIndex));
                                 if (isOwned) {
                                   // filter entries for zeros
                                   for(size_t domainIter = 0; domainIter < domainCardinality; domainIter++) {
                                     Scalar val = basisCoeffsLIOriented_d(cellNo, rangeIter, domainIter);
                                     if (KAT::magnitude(val) > entryFilterTol) {

#if defined(PANZER_INTERPOLATION_DEBUG_OUTPUT) || defined(PANZER_DEBUG)
                                       {
                                         // Check that there is no entry yet or that we are overwriting it with the same value
                                         auto row = lcl_tp_interp_matrix.rowConst(range_row);
                                         for(LocalOrdinal kk = 0; kk<row.length; ++kk)
                                           if (row.colidx(kk) == domainLIDs_d(domainOffsets_d(domainIter)))
                                             if (!(KAT::magnitude(row.value(kk)-val) < entryFilterTol || KAT::magnitude(row.value(kk)) < entryFilterTol)) {
#ifdef PANZER_INTERPOLATION_DEBUG_OUTPUT
                                               std::cout << "Replacing (" << range_row << "," << row.colidx(kk) << ") = " << row.value(kk) << " with " << val << std::endl;
#endif
#ifdef PANZER_DEBUG
                                               TEUCHOS_ASSERT(false);
#endif
                                             }
                                       }
#endif
#ifdef PANZER_INTERPOLATION_DEBUG_OUTPUT
                                       std::cout << "Setting (" << range_row << "," << domainLIDs_d(domainOffsets_d(domainIter)) << ") = " << val << std::endl;
#endif
				       lcl_tp_interp_matrix.replaceValues(range_row, &(domainLIDs_d(domainOffsets_d(domainIter))), 1, &val, /*is_sorted=*/false, /*force_atomic=*/true);
                                     }
                                   } //end if owned
                                 } // isOwned
                               } // end range LID loop
                             }); //end workset loop
      } // Tpetra fill
    } //end element loop
  } //end element block loop

  if (useTpetra) {
    tp_interp_matrix->fillComplete(tp_domainmap, tp_rangemap);

    if (op != Intrepid2::OPERATOR_VALUE) {
      // Discrete gradient, curl, etc actually live on edges, faces, etc, so we have a lot of zeros in the matrix.
      tp_interp_matrix = removeSmallEntries(tp_interp_matrix, entryFilterTol);
    }

    thyra_interp = Thyra::tpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,typename tp_matrix::node_type>(Thyra::createVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal>(tp_rangemap),
                                                                                                          Thyra::createVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal>(tp_domainmap),
                                                                                                          tp_interp_matrix);
  }
#ifdef PANZER_HAVE_EPETRA_STACK
  else
    ep_interp_matrix->FillComplete(*ep_domainmap, *ep_rangemap);
#endif

  return thyra_interp;
}

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MatrixFreeInterpolationOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  MatrixFreeInterpolationOp(const Teuchos::RCP<const panzer::ConnManager> &conn,
                            const Teuchos::RCP<panzer::DOFManager> &_domain_ugi,
                            const Teuchos::RCP<panzer::DOFManager> &_range_ugi,
                            const std::string& _domain_basis_name,
                            const std::string& _range_basis_name,
                            Intrepid2::EOperator _op,
                            size_t _worksetSize) :
    name(""),
    domain_basis_name(_domain_basis_name),
    range_basis_name(_range_basis_name),
    op(_op),
    worksetSize(_worksetSize),
    domain_ugi(_domain_ugi),
    range_ugi(_range_ugi)
  {

    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_dynamic_cast;

    using OT  = Teuchos::OrdinalTraits<GlobalOrdinal>;

    // typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal> tp_matrix;
    typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal> tp_map;

    RCP<const tp_map> tp_rangemap;
    RCP<const tp_map> tp_domainmap;
    RCP<const tp_map> tp_rowmap;
    RCP<const tp_map> tp_colmap;

    {
      // build maps
      std::vector<GlobalOrdinal> gids;
      range_ugi->getOwnedIndices(gids);
      tp_rowmap = rcp(new tp_map(OT::invalid(), gids.data(), gids.size(), OT::zero(), range_ugi->getComm()));
      tp_rangemap = tp_rowmap;
      domain_ugi->getOwnedIndices(gids);
      tp_domainmap = rcp(new tp_map(OT::invalid(), gids.data(), gids.size(), OT::zero(), domain_ugi->getComm()));
      domain_ugi->getOwnedAndGhostedIndices(gids);
      tp_colmap = rcp(new tp_map(OT::invalid(), gids.data(), gids.size(), OT::zero(), domain_ugi->getComm()));
    }

    domainMap_ = tp_domainmap;
    rangeMap_ = tp_rangemap;
    columnMap_ = tp_colmap;
    import_ = rcp(new Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node>(domainMap_, columnMap_));

    allocateColumnMapVector(1);

    precomputeOwnersAndOrientations(conn);

  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MatrixFreeInterpolationOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  allocateColumnMapVector(size_t numVectors) {
    colmapMV_ = rcp(new Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>(columnMap_, numVectors));
  }

  // Pre-compute elements that own a DoF
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MatrixFreeInterpolationOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  precomputeOwnersAndOrientations(const Teuchos::RCP<const panzer::ConnManager> &conn) {
    size_t maxNumElementsPerBlock = 0;
    int numCells;
    if (maxNumElementsPerBlock > 0)
      numCells = std::min(maxNumElementsPerBlock, worksetSize);
    else
      numCells = worksetSize;

    LocalOrdinal lclTargetSize = getRangeMap()->getLocalNumElements();
    owner_d_ = Kokkos::View<LocalOrdinal*,DeviceSpace>("owner", lclTargetSize);

    auto rangeLIDs_d = range_ugi->getLIDs();

    auto owner_d = owner_d_;

    // loop over element blocks
    std::vector<std::string> elementBlockIds;
    range_ugi->getElementBlockIds(elementBlockIds);
    for(std::size_t blockIter = 0; blockIter < elementBlockIds.size(); ++blockIter) {

      // loop over element worksets
      std::vector<int> elementIds = range_ugi->getElementBlock(elementBlockIds[blockIter]);
      Kokkos::View<int*,DeviceSpace>::HostMirror elementIds_h(elementIds.data(), elementIds.size());
      Kokkos::View<int*,DeviceSpace> elementIds_d("elementIds_d", elementIds_h.extent(0));
      Kokkos::deep_copy(elementIds_d, elementIds_h);
      for(std::size_t elemIter = 0; elemIter < elementIds_d.extent(0); elemIter += numCells) {
        using range_type = Kokkos::RangePolicy<LocalOrdinal, DeviceSpace>;
        Kokkos::parallel_for("miniEM::MatrixFreeInterpolationOp::cellLoop",
                             range_type(elemIter, std::min(elemIter+numCells,
                                                           elementIds_d.extent(0))),
                             KOKKOS_LAMBDA(const LocalOrdinal cellNo2) {
                               auto elemId = elementIds_d(cellNo2);

                               // loop over range LIDs
                               for(size_t rangeIter = 0; rangeIter < rangeLIDs_d.extent(1); ++rangeIter) {
                                 LocalOrdinal range_row = rangeLIDs_d(elemId, rangeIter);
                                 if (range_row < lclTargetSize)
                                   owner_d(range_row) = elemId;
                               }
                             });
      }
    }

    // get orientations for all blocks
    std::map<std::string, std::vector<Intrepid2::Orientation> > orientations;
    buildIntrepidOrientations(elementBlockIds, *conn, orientations);

    using HostSpace = Kokkos::HostSpace;

    // copy orientations to device
    for (auto const& orientation : orientations) {
      auto blockOrientations = orientation.second;
      orientations_[orientation.first] = typename Kokkos::DynRankView<Intrepid2::Orientation,DeviceSpace>("elemOrts_d", blockOrientations.size());
      Kokkos::View<Intrepid2::Orientation*, HostSpace> elemOrts_h(blockOrientations.data(), blockOrientations.size());
      Kokkos::deep_copy(orientations_[orientation.first], elemOrts_h);
    }
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MatrixFreeInterpolationOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  apply (const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X,
         Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y,
         Teuchos::ETransp mode,
         Scalar alpha,
         Scalar beta) const {

    if (mode == Teuchos::NO_TRANS) {
      applyNonTransposed(X, Y, alpha, beta);
      return;
    } else if (mode == Teuchos::TRANS) {
      applyTransposed(X, Y, alpha, beta);
      return;
    } else
      TEUCHOS_ASSERT(false);
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MatrixFreeInterpolationOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  applyNonTransposed (const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X,
                      Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y,
                      Scalar alpha,
                      Scalar beta) const {

    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_dynamic_cast;
    using range_type = Kokkos::RangePolicy<LocalOrdinal, DeviceSpace>;

    using ots = Intrepid2::OrientationTools<DeviceSpace>;
    using li = Intrepid2::LagrangianInterpolation<DeviceSpace>;
    using DynRankDeviceView = Kokkos::DynRankView<Scalar,DeviceSpace>;
    using view_t = typename Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::dual_view_type::t_dev;
    using const_view_t = typename Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::dual_view_type::t_dev::const_type;

    Teuchos::TimeMonitor tm(*Teuchos::TimeMonitor::getNewTimer(std::string("Mini-EM: matrix-free apply no_trans ") + name));

    using TST = Teuchos::ScalarTraits<Scalar>;
    const Scalar ZERO = TST::zero();
    if (beta == ZERO)
      Y.putScalar(ZERO);
    colmapMV_->doImport(X, *import_,Tpetra::INSERT);

    const_view_t lclX = colmapMV_->getLocalViewDevice(Tpetra::Access::ReadOnly);
    view_t lclY = Y.getLocalViewDevice(Tpetra::Access::ReadWrite);
    size_t numVectors = lclY.extent(1);
    LocalOrdinal lclTargetSize = getRangeMap()->getLocalNumElements();

    // get the domain and range bases
    auto domain_fieldPattern = domain_ugi->getFieldPattern(domain_basis_name);
    auto domain_basis = rcp_dynamic_cast<const panzer::Intrepid2FieldPattern>(domain_fieldPattern,true)->getIntrepidBasis();
    auto range_fieldPattern = range_ugi->getFieldPattern(range_basis_name);
    auto range_basis = rcp_dynamic_cast<const panzer::Intrepid2FieldPattern>(range_fieldPattern,true)->getIntrepidBasis();

    // cardinalities
    const size_t domainCardinality = domain_basis->getCardinality();
    const size_t rangeCardinality = range_basis->getCardinality();

    size_t maxNumElementsPerBlock = 0;

    const int dim = range_basis->getBaseCellTopology().getDimension();

    if (op == Intrepid2::OPERATOR_VALUE) {
      TEUCHOS_ASSERT(rangeCardinality >= domainCardinality);
      TEUCHOS_ASSERT_EQUALITY(domain_basis->getFunctionSpace(), range_basis->getFunctionSpace());
    }

    // allocate some views
    int numCells;
    if (maxNumElementsPerBlock > 0)
      numCells = std::min(maxNumElementsPerBlock, worksetSize);
    else
      numCells = worksetSize;
    DynRankDeviceView range_dofCoords_d("range_dofCoords_d", rangeCardinality, dim);
    DynRankDeviceView basisCoeffsLIOriented_d("basisCoeffsLIOriented_d", numCells, rangeCardinality, domainCardinality);

    // the ranks of these depend on dimension
    DynRankDeviceView valuesAtDofCoordsNonOriented_d;
    DynRankDeviceView valuesAtDofCoordsOriented_d;
    DynRankDeviceView reducedValuesAtDofCoordsOriented_d;

    {
      // Let Intrepid2 give us the correctly dimensioned view, then build one with +1 ranks and extent(0) == numCells
      auto temp = domain_basis->allocateOutputView(rangeCardinality, op);

      // These view have dimensions
      //  numCells, numFields=domainCardinality, numPoints=rangeCardinality, (spatialDim)
      //
      if (temp.rank() == 3) {
        valuesAtDofCoordsNonOriented_d     = DynRankDeviceView("valuesAtDofCoordsNonOriented_d", temp.extent(0), temp.extent(1), temp.extent(2));
        valuesAtDofCoordsOriented_d        = DynRankDeviceView("valuesAtDofCoordsOriented_d", numCells, temp.extent(0), temp.extent(1), temp.extent(2));
        reducedValuesAtDofCoordsOriented_d = DynRankDeviceView("reducedValuesAtDofCoordsOriented_d", numCells, temp.extent(1), temp.extent(2), numVectors);
      } else {
        valuesAtDofCoordsNonOriented_d     = DynRankDeviceView("valuesAtDofCoordsNonOriented_d", temp.extent(0), temp.extent(1));
        valuesAtDofCoordsOriented_d        = DynRankDeviceView("valuesAtDofCoordsOriented_d", numCells, temp.extent(0), temp.extent(1));
        reducedValuesAtDofCoordsOriented_d = DynRankDeviceView("reducedValuesAtDofCoordsOriented_d", numCells, temp.extent(1), numVectors);
      }
    }

    int fieldRank = Intrepid2::getFieldRank(range_basis->getFunctionSpace());
    TEUCHOS_ASSERT((fieldRank == 0) || (fieldRank == 1));

    auto rangeLIDs_d = range_ugi->getLIDs();
    auto domainLIDs_d = domain_ugi->getLIDs();

    // range dof coordinates
    range_basis->getDofCoords(range_dofCoords_d);

    // compute values of op * (domain basis) at range dof coords on reference element
    domain_basis->getValues(valuesAtDofCoordsNonOriented_d, range_dofCoords_d, op);

    // loop over element blocks
    std::vector<std::string> elementBlockIds;
    range_ugi->getElementBlockIds(elementBlockIds);
    for(std::size_t blockIter = 0; blockIter < elementBlockIds.size(); ++blockIter) {

      auto rangeOffsets_d = range_ugi->getGIDFieldOffsetsKokkos(elementBlockIds[blockIter], 0);
      auto domainOffsets_d = domain_ugi->getGIDFieldOffsetsKokkos(elementBlockIds[blockIter], 0);

      // get element ids
      Kokkos::View<int *, DeviceSpace> elementIds_d;
      {
        std::vector<int> elementIds = range_ugi->getElementBlock(elementBlockIds[blockIter]);
        Kokkos::View<int *, Kokkos::HostSpace> elementIds_h(elementIds.data(), elementIds.size());
        elementIds_d = Kokkos::View<int *, DeviceSpace>("elementIds_d", elementIds_h.extent(0));
        Kokkos::deep_copy(elementIds_d, elementIds_h);
      }

      // get element orientations
      auto elemOrts_d = orientations_.at(elementBlockIds[blockIter]);

      for(std::size_t elemIter = 0; elemIter < elementIds_d.extent(0); elemIter += numCells) {

        int endCellRange =
          std::min(numCells, Teuchos::as<int>(elementIds_d.extent(0)) -
                             Teuchos::as<int>(elemIter));

        // get subviews on workset
        auto ortsRange = Kokkos::make_pair(elemIter, std::min(elemIter + numCells, elemOrts_d.extent(0)));
        auto elemOrtsWorkset_d = Kokkos::subview(elemOrts_d, ortsRange);
        // Last workset might be shorter.
        auto worksetRange = Kokkos::make_pair(0, endCellRange);
        auto valuesAtDofCoordsOrientedWorkset_d = Kokkos::subview(valuesAtDofCoordsOriented_d, worksetRange, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
        auto basisCoeffsLIOrientedWorkset_d     = Kokkos::subview(basisCoeffsLIOriented_d,     worksetRange, Kokkos::ALL(), Kokkos::ALL());

        // apply orientations for LO basis
        // shuffles things in the second dimension, i.e. wrt LO basis
        ots::modifyBasisByOrientation(valuesAtDofCoordsOrientedWorkset_d,
                                      valuesAtDofCoordsNonOriented_d,
                                      elemOrtsWorkset_d,
                                      domain_basis.get());

        Kokkos::deep_copy(reducedValuesAtDofCoordsOriented_d, 0.0);

        if (reducedValuesAtDofCoordsOriented_d.rank() == 4) {
          Kokkos::parallel_for("miniEM:MatrixFreeInterpolationOp:cellLoop1",
                               range_type(0, std::min(numCells, elementIds_d.extent_int(0)-Teuchos::as<int>(elemIter))),
                               KOKKOS_LAMBDA(const LocalOrdinal cellNo) {
                                 LocalOrdinal cellNo2 = elemIter+cellNo;
                                 LocalOrdinal elemId = elementIds_d(cellNo2);
                                 for(size_t domainIter=0; domainIter<domainCardinality; domainIter++) {
                                   LocalOrdinal J = domainLIDs_d(elemId, domainOffsets_d(domainIter));
                                   for(size_t rangeIter=0; rangeIter<rangeCardinality; rangeIter++) {
                                     for(size_t d=0; d<valuesAtDofCoordsOriented_d.extent(3); d++) {
                                       Scalar val = valuesAtDofCoordsOriented_d(cellNo, domainIter, rangeIter, d);
                                       for (size_t j = 0; j<numVectors; ++j)
                                         reducedValuesAtDofCoordsOriented_d(cellNo, rangeIter, d, j) += val * lclX(J, j);
                                     }
                                   }
                                 }
          });

          for (size_t j = 0; j<numVectors; ++j)
            li::getBasisCoeffs(Kokkos::subview(basisCoeffsLIOrientedWorkset_d, Kokkos::ALL(), Kokkos::ALL(), j),
                               Kokkos::subview(reducedValuesAtDofCoordsOriented_d, worksetRange, Kokkos::ALL(), Kokkos::ALL(), j),
                               range_basis.get(),
                               elemOrtsWorkset_d
                               );
        } else {
          Kokkos::parallel_for("miniEM:MatrixFreeInterpolationOp:cellLoop1",
                               range_type(0, std::min(numCells, elementIds_d.extent_int(0)-Teuchos::as<int>(elemIter))),
                               KOKKOS_LAMBDA(const LocalOrdinal cellNo) {
                                 LocalOrdinal cellNo2 = elemIter+cellNo;
                                 LocalOrdinal elemId = elementIds_d(cellNo2);
                                 for(size_t domainIter=0; domainIter<domainCardinality; domainIter++) {
                                   LocalOrdinal J = domainLIDs_d(elemId, domainOffsets_d(domainIter));
                                   for(size_t rangeIter=0; rangeIter<rangeCardinality; rangeIter++) {
                                     Scalar val = valuesAtDofCoordsOriented_d(cellNo, domainIter, rangeIter);
                                     for (size_t j = 0; j<numVectors; ++j)
                                       reducedValuesAtDofCoordsOriented_d(cellNo, rangeIter, j) += val * lclX(J, j);
                                   }
                                 }
          });

          for (size_t j = 0; j<numVectors; ++j)
            li::getBasisCoeffs(Kokkos::subview(basisCoeffsLIOrientedWorkset_d, Kokkos::ALL(), Kokkos::ALL(), j),
                               Kokkos::subview(reducedValuesAtDofCoordsOriented_d, worksetRange, Kokkos::ALL(), j),
                               range_basis.get(),
                               elemOrtsWorkset_d
                               );
        }

        auto owner_d = owner_d_;

        Kokkos::parallel_for("miniEM::MatrixFreeInterpolationOp::cellLoop2",
                             range_type(elemIter, std::min(elemIter+numCells,
                                                           elementIds_d.extent(0))),
                             KOKKOS_LAMBDA(const LocalOrdinal cellNo2) {
                               LocalOrdinal cellNo = cellNo2-elemIter;
                               LocalOrdinal elemId = elementIds_d(cellNo2);

                               // loop over range LIDs
                               for(size_t rangeIter = 0; rangeIter < rangeLIDs_d.extent(1); ++rangeIter) {
                                 LocalOrdinal range_row = rangeLIDs_d(elemId, rangeOffsets_d(rangeIter));

                                 // if owned
                                 if ((range_row < lclTargetSize) && (owner_d(range_row) == elemId)) {

                                   for (size_t j = 0; j<numVectors; ++j) {
                                     Scalar val = basisCoeffsLIOriented_d(cellNo, rangeIter, j);
                                     lclY(range_row,j) = beta*lclY(range_row,j) + alpha*val;
                                   }
                                 } // end if owned
                               } // end range LID loop
                             }); // end element loop
      } //end workset loop
    } //end element block loop

  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MatrixFreeInterpolationOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  applyTransposed (const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X,
                   Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y,
                   Scalar alpha,
                   Scalar beta) const {

    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_dynamic_cast;
    using range_type = Kokkos::RangePolicy<LocalOrdinal, DeviceSpace>;

    typedef Intrepid2::OrientationTools<DeviceSpace> ots;
    typedef Intrepid2::LagrangianInterpolation<DeviceSpace> li;
    typedef Kokkos::DynRankView<Scalar,DeviceSpace> DynRankDeviceView;
    using TST = Teuchos::ScalarTraits<Scalar>;
    const Scalar ZERO = TST::zero();
    using view_t = typename Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::dual_view_type::t_dev;
    using const_view_t = typename Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::dual_view_type::t_dev::const_type;

    Teuchos::TimeMonitor tm(*Teuchos::TimeMonitor::getNewTimer(std::string("Mini-EM: matrix-free apply trans ") + name));

    const_view_t lclX = X.getLocalViewDevice(Tpetra::Access::ReadOnly);
    colmapMV_->putScalar(ZERO);
    view_t lclYtemp = colmapMV_->getLocalViewDevice(Tpetra::Access::ReadWrite);

    // get the domain and range bases
    auto domain_fieldPattern = domain_ugi->getFieldPattern(domain_basis_name);
    auto domain_basis = rcp_dynamic_cast<const panzer::Intrepid2FieldPattern>(domain_fieldPattern,true)->getIntrepidBasis();
    auto range_fieldPattern = range_ugi->getFieldPattern(range_basis_name);
    auto range_basis = rcp_dynamic_cast<const panzer::Intrepid2FieldPattern>(range_fieldPattern,true)->getIntrepidBasis();

    // cardinalities
    const size_t domainCardinality = domain_basis->getCardinality();
    const size_t rangeCardinality = range_basis->getCardinality();

    size_t maxNumElementsPerBlock = 0;

    const int dim = range_basis->getBaseCellTopology().getDimension();

    if (op == Intrepid2::OPERATOR_VALUE) {
      TEUCHOS_ASSERT(rangeCardinality >= domainCardinality);
      TEUCHOS_ASSERT_EQUALITY(domain_basis->getFunctionSpace(), range_basis->getFunctionSpace());
    }

    // allocate some views
    int numCells;
    if (maxNumElementsPerBlock > 0)
      numCells = std::min(maxNumElementsPerBlock, worksetSize);
    else
      numCells = worksetSize;
    DynRankDeviceView range_dofCoords_d("range_dofCoords_d", rangeCardinality, dim);
    DynRankDeviceView basisCoeffsLIOriented_d("basisCoeffsLIOriented_d", numCells, rangeCardinality, domainCardinality);

    // the ranks of these depend on dimension
    DynRankDeviceView valuesAtDofCoordsNonOriented_d;
    DynRankDeviceView valuesAtDofCoordsOriented_d;

    {
      // Let Intrepid2 give us the correctly dimensioned view, then build one with +1 ranks and extent(0) == numCells
      auto temp = domain_basis->allocateOutputView(rangeCardinality, op);

      // These view have dimensions
      //  numCells, numFields=domainCardinality, numPoints=rangeCardinality, (spatialDim)
      //
      if (temp.rank() == 3) {
        valuesAtDofCoordsNonOriented_d = DynRankDeviceView("valuesAtDofCoordsNonOriented_d", temp.extent(0), temp.extent(1), temp.extent(2));
        valuesAtDofCoordsOriented_d    = DynRankDeviceView("valuesAtDofCoordsOriented_d", numCells, temp.extent(0), temp.extent(1), temp.extent(2));
      } else {
        valuesAtDofCoordsNonOriented_d = DynRankDeviceView("valuesAtDofCoordsNonOriented_d", temp.extent(0), temp.extent(1));
        valuesAtDofCoordsOriented_d    = DynRankDeviceView("valuesAtDofCoordsOriented_d", numCells, temp.extent(0), temp.extent(1));
      }
    }

    int fieldRank = Intrepid2::getFieldRank(range_basis->getFunctionSpace());
    TEUCHOS_ASSERT((fieldRank == 0) || (fieldRank == 1));

    auto rangeLIDs_d = range_ugi->getLIDs();
    auto domainLIDs_d = domain_ugi->getLIDs();
    Kokkos::fence();

    // range dof coordinates
    range_basis->getDofCoords(range_dofCoords_d);

    // compute values of op * (domain basis) at range dof coords on reference element
    domain_basis->getValues(valuesAtDofCoordsNonOriented_d, range_dofCoords_d, op);

    // loop over element blocks
    std::vector<std::string> elementBlockIds;
    range_ugi->getElementBlockIds(elementBlockIds);
    for(std::size_t blockIter = 0; blockIter < elementBlockIds.size(); ++blockIter) {

      auto rangeOffsets_d = range_ugi->getGIDFieldOffsetsKokkos(elementBlockIds[blockIter], 0);
      auto domainOffsets_d = domain_ugi->getGIDFieldOffsetsKokkos(elementBlockIds[blockIter], 0);

      // get element ids
      Kokkos::View<int *, DeviceSpace> elementIds_d;
      {
        std::vector<int> elementIds = range_ugi->getElementBlock(elementBlockIds[blockIter]);
        Kokkos::View<int *, Kokkos::HostSpace> elementIds_h(elementIds.data(), elementIds.size());
        elementIds_d = Kokkos::View<int *, DeviceSpace>("elementIds_d", elementIds_h.extent(0));
        Kokkos::deep_copy(elementIds_d, elementIds_h);
      }

      // get element orientations
      auto elemOrts_d = orientations_.at(elementBlockIds[blockIter]);

      for(std::size_t elemIter = 0; elemIter < elementIds_d.extent(0); elemIter += numCells) {

        int endCellRange =
          std::min(numCells, Teuchos::as<int>(elementIds_d.extent(0)) -
                             Teuchos::as<int>(elemIter));

        // get subviews on workset
        auto ortsRange = Kokkos::make_pair(elemIter, std::min(elemIter + numCells, elemOrts_d.extent(0)));
        auto elemOrtsWorkset_d = Kokkos::subview(elemOrts_d, ortsRange);
        // Last workset might be shorter.
        auto worksetRange = Kokkos::make_pair(0, endCellRange);
        auto valuesAtDofCoordsOrientedWorkset_d = Kokkos::subview(valuesAtDofCoordsOriented_d, worksetRange, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
        auto basisCoeffsLIOrientedWorkset_d     = Kokkos::subview(basisCoeffsLIOriented_d,     worksetRange, Kokkos::ALL(), Kokkos::ALL());

        // apply orientations for domain basis
        // shuffles things in the second dimension, i.e. wrt domain basis
        ots::modifyBasisByOrientation(valuesAtDofCoordsOrientedWorkset_d,
                                      valuesAtDofCoordsNonOriented_d,
                                      elemOrtsWorkset_d,
                                      domain_basis.get());
        Kokkos::fence();

        // get basis coefficients of domain basis functions wrt range basis
        for(size_t domainIter=0; domainIter<domainCardinality; domainIter++)
          // Get basis coeffs wrt range basis on reference element.
          // basisCoeffsLI has dimensions (numCells, numFields=rangeCardinality, domainCardinality)
          li::getBasisCoeffs(Kokkos::subview(basisCoeffsLIOrientedWorkset_d, Kokkos::ALL(), Kokkos::ALL(), domainIter),
                             Kokkos::subview(valuesAtDofCoordsOrientedWorkset_d, Kokkos::ALL(), domainIter, Kokkos::ALL(), Kokkos::ALL()),
                             range_basis.get(), elemOrtsWorkset_d);
        Kokkos::fence();

        auto owner_d = owner_d_;


        Kokkos::parallel_for("miniEM::MatrixFreeInterpolationOp::cellLoop",
                             range_type(elemIter, std::min(elemIter+numCells,
                                                           elementIds_d.extent(0))),
                             KOKKOS_LAMBDA(const LocalOrdinal cellNo2) {
                               LocalOrdinal cellNo = cellNo2-elemIter;
                               LocalOrdinal elemId = elementIds_d(cellNo2);

                               // loop over range LIDs
                               for(size_t rangeIter = 0; rangeIter < rangeLIDs_d.extent(1); ++rangeIter) {
                                 LocalOrdinal range_row = rangeLIDs_d(elemId, rangeOffsets_d(rangeIter));

                                 // if owned
                                 if ((range_row < (LocalOrdinal) lclX.extent(0)) && (owner_d(range_row) == elemId)) {

                                   for(size_t domainIter = 0; domainIter < domainLIDs_d.extent(1); domainIter++) {
                                     LocalOrdinal J = domainLIDs_d(elemId, domainOffsets_d(domainIter));
                                     Scalar val = basisCoeffsLIOriented_d(cellNo, rangeIter, domainIter);
                                     for (size_t j = 0; j<lclYtemp.extent(1); ++j)
                                       Kokkos::atomic_add(&lclYtemp(J,j), alpha*val*lclX(range_row,j));
                                   }
                                 } //end if owned
                               } //end range LID loop
                             }); // end element loop
      } //end workset loop
    } //end element block loop
    Kokkos::fence();

    if (beta == ZERO)
      Y.putScalar(ZERO);
    else
      Y.scale(beta);
    Y.doExport(*colmapMV_, *import_, Tpetra::ADD_ASSIGN);
    Kokkos::fence();
  }


}
