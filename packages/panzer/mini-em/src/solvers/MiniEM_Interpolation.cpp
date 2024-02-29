#include "MiniEM_Interpolation.hpp"


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


Teko::LinearOp buildInterpolation(const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > &linObjFactory,
                                  const std::string& lo_basis_name,
                                  const std::string& ho_basis_name,
                                  Intrepid2::EOperator op,
                                  size_t worksetSize)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;

  using Scalar = double;
  using LocalOrdinal = int;
  using GlobalOrdinal = panzer::GlobalOrdinal;

  using STS = Teuchos::ScalarTraits<Scalar>;
  using KAT = Kokkos::ArithTraits<Scalar>;
  using OT  = Teuchos::OrdinalTraits<GlobalOrdinal>;

  using tpetraBlockedLinObjFactory = typename panzer::BlockedTpetraLinearObjFactory<panzer::Traits, Scalar, LocalOrdinal, GlobalOrdinal>;
#ifdef PANZER_HAVE_EPETRA_STACK
  using epetraBlockedLinObjFactory = typename panzer::BlockedEpetraLinearObjFactory<panzer::Traits, LocalOrdinal>;
#endif
  using UGI = panzer::GlobalIndexer;
  using DeviceSpace = PHX::Device;
  using HostSpace = Kokkos::HostSpace;
  using ots = Intrepid2::OrientationTools<DeviceSpace>;
  using li = Intrepid2::LagrangianInterpolation<DeviceSpace>;
  using DynRankDeviceView = Kokkos::DynRankView<double, DeviceSpace>;

  // must be able to cast to a block linear object factory
  RCP<const tpetraBlockedLinObjFactory > tblof = rcp_dynamic_cast<const tpetraBlockedLinObjFactory >(linObjFactory);
#ifdef PANZER_HAVE_EPETRA_STACK
  RCP<const epetraBlockedLinObjFactory > eblof = rcp_dynamic_cast<const epetraBlockedLinObjFactory >(linObjFactory);
#endif

  using tp_graph = Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal>;
  using tp_matrix = Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal>;
  using tp_map = Tpetra::Map<LocalOrdinal, GlobalOrdinal>;
#ifdef PANZER_HAVE_EPETRA_STACK
  using ep_linObjContainer = panzer::BlockedEpetraLinearObjContainer;
  using ep_matrix = Epetra_CrsMatrix;
  using ep_map = Epetra_Map;
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

  // get global indexers for LO and HO dofs
  std::vector<RCP<UGI> > fieldDOFMngrs = blockedDOFMngr->getFieldDOFManagers();
  int loFieldNum = blockedDOFMngr->getFieldNum(lo_basis_name);
  int hoFieldNum = blockedDOFMngr->getFieldNum(ho_basis_name);
  int loBlockIndex = blockedDOFMngr->getFieldBlock(loFieldNum);
  int hoBlockIndex = blockedDOFMngr->getFieldBlock(hoFieldNum);
  RCP<panzer::DOFManager> lo_ugi = rcp_dynamic_cast<panzer::DOFManager>(blockedDOFMngr->getFieldDOFManagers()[loBlockIndex],true);
  RCP<panzer::DOFManager> ho_ugi = rcp_dynamic_cast<panzer::DOFManager>(blockedDOFMngr->getFieldDOFManagers()[hoBlockIndex],true);

  // get the LO and HO bases
  auto lo_fieldPattern = lo_ugi->getFieldPattern(lo_basis_name);
  auto lo_basis = rcp_dynamic_cast<const panzer::Intrepid2FieldPattern>(lo_fieldPattern,true)->getIntrepidBasis();
  auto ho_fieldPattern = ho_ugi->getFieldPattern(ho_basis_name);
  auto ho_basis = rcp_dynamic_cast<const panzer::Intrepid2FieldPattern>(ho_fieldPattern,true)->getIntrepidBasis();

  // cardinalities
  const size_t loCardinality = lo_basis->getCardinality();
  const size_t hoCardinality = ho_basis->getCardinality();

  auto hoLIDs_h = Kokkos::View<LocalOrdinal*,HostSpace>("hoLIDs_h",hoCardinality);
  auto loLIDs_h = Kokkos::View<LocalOrdinal*,HostSpace>("loLIDs_h",loCardinality);
  size_t maxNumElementsPerBlock = 0;

  // Create the global interp matrix.
  // The operator maps from LO (domain) to HO (range)
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

  auto hoElementLIDs_d = ho_ugi->getLIDs();
  auto loElementLIDs_d = lo_ugi->getLIDs();

  RCP<Thyra::LinearOpBase<Scalar> > thyra_interp;
  LocalOrdinal minLocalIndex = 0;
  LocalOrdinal maxLocalIndex = 0;
  if (tblof != Teuchos::null) {
    // build maps
    std::vector<GlobalOrdinal> gids;
    ho_ugi->getOwnedIndices(gids);
    tp_rowmap = rcp(new tp_map(OT::invalid(), gids.data(), static_cast<LocalOrdinal>(gids.size()), OT::zero(), ho_ugi->getComm()));
    tp_rangemap = tp_rowmap;
    lo_ugi->getOwnedIndices(gids);
    tp_domainmap = rcp(new tp_map(OT::invalid(), gids.data(), static_cast<LocalOrdinal>(gids.size()), OT::zero(), lo_ugi->getComm()));
    lo_ugi->getOwnedAndGhostedIndices(gids);
    tp_colmap = rcp(new tp_map(OT::invalid(), gids.data(), static_cast<LocalOrdinal>(gids.size()), OT::zero(), lo_ugi->getComm()));

    minLocalIndex = tp_rowmap->getMinLocalIndex();
    maxLocalIndex = tp_rowmap->getMaxLocalIndex();

    // estimate number of entries per row
    // This is an upper bound, as we are counting dofs that are on shared nodes, edges, faces more than once.
    // Kokkos::View<size_t*,HostSpace> numEntriesPerRow("numEntriesPerRow", tp_rowmap->getLocalNumElements());
    Kokkos::DualView<size_t*, DeviceSpace> numEntriesPerRow("numEntriesPerRow", tp_rowmap->getLocalNumElements());

    {
      auto numEntriesPerRow_d = numEntriesPerRow.view_device();

      // loop over element blocks
      std::vector<std::string> elementBlockIds;
      blockedDOFMngr->getElementBlockIds(elementBlockIds);
      for(std::size_t blockIter = 0; blockIter < elementBlockIds.size(); ++blockIter) {

        // loop over elements
        std::vector<int> elementIds = ho_ugi->getElementBlock(elementBlockIds[blockIter]);
        Kokkos::View<int *, HostSpace> elementIds_h(elementIds.data(), elementIds.size());
        Kokkos::View<int *, DeviceSpace> elementIds_d("elementIds_d", elementIds_h.extent(0));
        Kokkos::deep_copy(elementIds_d, elementIds_h);
        maxNumElementsPerBlock = std::max(maxNumElementsPerBlock, elementIds.size());

        Kokkos::parallel_for("MiniEM_Interpolation::numEntriesPerRow",
                             Kokkos::RangePolicy<size_t, typename tp_matrix::node_type::execution_space>(0, elementIds.size()),
                             KOKKOS_LAMBDA(const size_t elemIter) {
                               auto elemId = elementIds_d(elemIter);

                               // get IDs for HO dofs
                               auto hoLIDs_d = Kokkos::subview(hoElementLIDs_d, elemId, Kokkos::ALL());

                               // loop over HO LIDs
                               for(size_t hoIter = 0; hoIter < hoCardinality; ++hoIter) {
                                 const LocalOrdinal ho_row = hoLIDs_d(hoIter);
                                 const bool isOwned = ((minLocalIndex <= ho_row) && (ho_row <= maxLocalIndex));
                                 if (isOwned)
				   Kokkos::atomic_add(&numEntriesPerRow_d(ho_row), loCardinality);
                               } //end HO LID loop
                             });
      } // blocks loop
      numEntriesPerRow.template modify<DeviceSpace>();
      numEntriesPerRow.template sync<HostSpace>();
    }

    // Set up graph
    auto tp_interp_graph = rcp(new tp_graph(tp_rowmap, tp_colmap, numEntriesPerRow));

    { // This runs on host
      Kokkos::View<LocalOrdinal**, HostSpace> hoElementLIDs_h("hoElementLIDs_h", hoElementLIDs_d.extent(0), hoCardinality);
      Kokkos::View<LocalOrdinal**, HostSpace> loElementLIDs_h("loElementLIDs_h", loElementLIDs_d.extent(0), loCardinality);
      Kokkos::deep_copy(hoElementLIDs_h, hoElementLIDs_d);
      Kokkos::deep_copy(loElementLIDs_h, loElementLIDs_d);

      // loop over element blocks
      std::vector<std::string> elementBlockIds;
      blockedDOFMngr->getElementBlockIds(elementBlockIds);
      for(std::size_t blockIter = 0; blockIter < elementBlockIds.size(); ++blockIter) {

        // loop over elements
        std::vector<int> elementIds = ho_ugi->getElementBlock(elementBlockIds[blockIter]);
        maxNumElementsPerBlock = std::max(maxNumElementsPerBlock, elementIds.size());
        for(std::size_t elemIter = 0; elemIter < elementIds.size(); ++elemIter) {
          auto elemId = elementIds[elemIter];

          // get IDs for HO dofs
          auto hoLIDs_h = Kokkos::subview(hoElementLIDs_h, elemId, Kokkos::ALL());
          auto loLIDs_h = Kokkos::subview(loElementLIDs_h, elemId, Kokkos::ALL());

          // loop over HO LIDs
          for(size_t hoIter = 0; hoIter < hoCardinality; ++hoIter) {
            const LocalOrdinal ho_row = hoLIDs_h(hoIter);
            const bool isOwned = ((minLocalIndex <= ho_row) && (ho_row <= maxLocalIndex));
            if (isOwned) {
              Teuchos::ArrayView<LocalOrdinal> loLIDs_av = Teuchos::ArrayView<LocalOrdinal>(loLIDs_h.data(), loLIDs_h.extent_int(0));
              tp_interp_graph->insertLocalIndices(ho_row, loLIDs_av);
            }
          } //end HO LID loop
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
  else if (eblof != Teuchos::null) {
    RCP<panzer::GlobalEvaluationData> dataObject
      = rcp(new panzer::LOCPair_GlobalEvaluationData(eblof,panzer::LinearObjContainer::Mat));
    RCP<panzer::LinearObjContainer> global_loc
      = rcp_dynamic_cast<panzer::LOCPair_GlobalEvaluationData>(dataObject,true)->getGlobalLOC();
    RCP<panzer::LinearObjContainer> ghosted_loc
      = rcp_dynamic_cast<panzer::LOCPair_GlobalEvaluationData>(dataObject,true)->getGhostedLOC();

    RCP<ep_linObjContainer> global_eloc = rcp_dynamic_cast<ep_linObjContainer>(global_loc,true);
    RCP<ep_linObjContainer> ghosted_eloc = rcp_dynamic_cast<ep_linObjContainer>(ghosted_loc,true);

    ep_rangemap  = global_eloc->getMapForBlock(hoBlockIndex);
    ep_domainmap = global_eloc->getMapForBlock(loBlockIndex);
    ep_rowmap    = global_eloc->getMapForBlock(hoBlockIndex);
    ep_colmap    = ghosted_eloc->getMapForBlock(loBlockIndex);

    {
      // loop over element blocks
      std::vector<std::string> elementBlockIds;
      blockedDOFMngr->getElementBlockIds(elementBlockIds);
      for(std::size_t blockIter = 0; blockIter < elementBlockIds.size(); ++blockIter) {

        // loop over elements
        std::vector<int> elementIds = ho_ugi->getElementBlock(elementBlockIds[blockIter]);
        maxNumElementsPerBlock = std::max(maxNumElementsPerBlock, elementIds.size());
      }
    }

    // TODO: Fix this.
    size_t nnzPerRowEstimate = 25*loCardinality;

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

  RCP<const panzer::ConnManager> conn = blockedDOFMngr->getConnManager();

  // assume only one cell toplogy per rank
  std::vector<shards::CellTopology> topologies;
  conn->getElementBlockTopologies(topologies);
  shards::CellTopology topology = topologies[0];
  int dim = static_cast<int>(topology.getDimension());
  // num vertices in an element
  const int numElemVertices = static_cast<int>(topology.getVertexCount());

  // set up a node only conn manager
  auto node_basis = panzer::createIntrepid2Basis<DeviceSpace,Scalar,Scalar>("HGrad",1,topology);
  auto node_fieldPattern = rcp(new panzer::Intrepid2FieldPattern(node_basis));
  RCP<panzer::ConnManager> node_conn = Teuchos::rcp_dynamic_cast<panzer::ConnManager>(conn->noConnectivityClone(),true);
  node_conn->buildConnectivity(*node_fieldPattern);

  if (op == Intrepid2::OPERATOR_VALUE) {
    TEUCHOS_ASSERT(hoCardinality >= loCardinality);
    TEUCHOS_ASSERT_EQUALITY(lo_basis->getFunctionSpace(), ho_basis->getFunctionSpace());
  }

  // allocate some views
  int numCells;
  if (maxNumElementsPerBlock > 0)
    numCells = static_cast<int>(std::min(maxNumElementsPerBlock, worksetSize));
  else
    numCells = static_cast<int>(worksetSize);
  DynRankDeviceView             ho_dofCoords_d("ho_dofCoords_d", hoCardinality, dim);
  DynRankDeviceView             basisCoeffsLIOriented_d("basisCoeffsLIOriented_d", numCells, hoCardinality, loCardinality);
  typename Kokkos::DynRankView<Intrepid2::Orientation,DeviceSpace>     elemOrts_d ("elemOrts_d",  numCells);
  typename Kokkos::DynRankView<GlobalOrdinal, DeviceSpace>::HostMirror elemNodes_h("elemNodes_h", numCells, numElemVertices);
  typename Kokkos::DynRankView<GlobalOrdinal, DeviceSpace>             elemNodes_d("elemNodes_d", numCells, numElemVertices);

  // the ranks of these depend on dimension
  DynRankDeviceView ho_dofCoeffs_d;
  DynRankDeviceView valuesAtDofCoordsNonOriented_d;
  DynRankDeviceView valuesAtDofCoordsOriented_d;

  {
    // Let Intrepid2 give us the correctly dimensioned view, then build one with +1 ranks and extent(0) == numCells
    auto temp = lo_basis->allocateOutputView(static_cast<int>(hoCardinality), op);

    // These view have dimensions
    //  numCells, numFields=loCardinality, numPoints=hoCardinality, (spatialDim)
    //
    if (temp.rank() == 3) {
      valuesAtDofCoordsNonOriented_d = DynRankDeviceView("valuesAtDofCoordsNonOriented_d", temp.extent(0), temp.extent(1), temp.extent(2));
      valuesAtDofCoordsOriented_d    = DynRankDeviceView("valuesAtDofCoordsOriented_d", numCells, temp.extent(0), temp.extent(1), temp.extent(2));
    } else {
      valuesAtDofCoordsNonOriented_d = DynRankDeviceView("valuesAtDofCoordsNonOriented_d", temp.extent(0), temp.extent(1));
      valuesAtDofCoordsOriented_d    = DynRankDeviceView("valuesAtDofCoordsOriented_d", numCells, temp.extent(0), temp.extent(1));
    }
  }

  int fieldRank = Intrepid2::getFieldRank(ho_basis->getFunctionSpace());
  TEUCHOS_ASSERT((fieldRank == 0) || (fieldRank == 1));

  auto entryFilterTol = 100*Teuchos::ScalarTraits<typename STS::magnitudeType>::eps();

  // loop over element blocks
  std::vector<std::string> elementBlockIds;
  blockedDOFMngr->getElementBlockIds(elementBlockIds);
  for(std::size_t blockIter = 0; blockIter < elementBlockIds.size(); ++blockIter) {

    // loop over element worksets
    std::vector<int> elementIds = ho_ugi->getElementBlock(elementBlockIds[blockIter]);
    Kokkos::View<int *, HostSpace> elementIds_h(elementIds.data(), elementIds.size());
    Kokkos::View<int *, DeviceSpace> elementIds_d("elementIds_d", elementIds_h.extent(0));
    Kokkos::deep_copy(elementIds_d, elementIds_h);

    for(std::size_t elemIter = 0; elemIter < elementIds.size(); elemIter += numCells) {

      // get element orientations
      for (int cellNo = 0; cellNo < numCells; cellNo++) {
        if (elemIter+cellNo >= elementIds.size())
          continue;
        const GlobalOrdinal* node_ids = node_conn->getConnectivity(elementIds[elemIter+cellNo]);
        for(int i = 0; i < numElemVertices; i++)
          elemNodes_h(cellNo, i) = node_ids[i];
      }
      Kokkos::deep_copy(elemNodes_d, elemNodes_h);

      ots::getOrientation(elemOrts_d, elemNodes_d, topology);

      // HO dof coordinates and coefficients
      ho_basis->getDofCoords(ho_dofCoords_d);

      // compute values of op * (LO basis) at HO dof coords on reference element
      lo_basis->getValues(valuesAtDofCoordsNonOriented_d, ho_dofCoords_d, op);

      // apply orientations for LO basis
      // shuffles things in the second dimension, i.e. wrt LO basis
      ots::modifyBasisByOrientation(valuesAtDofCoordsOriented_d,
                                    valuesAtDofCoordsNonOriented_d,
                                    elemOrts_d,
                                    lo_basis.get());

      //get basis coefficients of LO basis functions wrt HO basis
      for(size_t loIter=0; loIter<loCardinality; loIter++)
        // Get basis coeffs wrt HO basis on reference element.
        // basisCoeffsLI has dimensions (numCells, numFields=hoCardinality, loCardinality)
        li::getBasisCoeffs(Kokkos::subview(basisCoeffsLIOriented_d, Kokkos::ALL(), Kokkos::ALL(), loIter),
                           Kokkos::subview(valuesAtDofCoordsOriented_d, Kokkos::ALL(), loIter, Kokkos::ALL(), Kokkos::ALL()),
                           ho_basis.get(), elemOrts_d);

      int endCellRange =
          std::min(numCells, Teuchos::as<int>(elementIds.size()) -
                                 Teuchos::as<int>(elemIter));

#ifdef PANZER_HAVE_EPETRA_STACK
      if (tblof.is_null()) { // Epetra fill

	Kokkos::View<LocalOrdinal*,DeviceSpace> indices_d("indices", loCardinality);
	Kokkos::View<Scalar*,      DeviceSpace> values_d ("values",  loCardinality);

	
        for (int cellNo = 0; cellNo < endCellRange; cellNo++) {
          auto elemId = elementIds_d(elemIter+cellNo);

          // get IDs for HO and LO dofs
          auto hoLIDs_d = Kokkos::subview(hoElementLIDs_d, elemId, Kokkos::ALL());
          auto loLIDs_d = Kokkos::subview(loElementLIDs_d, elemId, Kokkos::ALL());

          // loop over HO LIDs
          for(size_t hoIter = 0; hoIter < hoCardinality; ++hoIter) {
            const LocalOrdinal ho_row = hoLIDs_d(hoIter);
            const bool isOwned = ep_rowmap->MyLID(ho_row);
            if (isOwned) {
              // filter entries for zeros
              LocalOrdinal rowNNZ = 0;
              for(size_t loIter = 0; loIter < loCardinality; loIter++) {
                Scalar val = basisCoeffsLIOriented_d(cellNo, hoIter, loIter);
                if (KAT::magnitude(val) > entryFilterTol) {
                  indices_d(rowNNZ) = loLIDs_d(loIter);
                  values_d(rowNNZ) = val;
                  ++rowNNZ;
                }
              }

              int ret = ep_interp_matrix->ReplaceMyValues(ho_row, rowNNZ, values_d.data(), indices_d.data());
              if (ret != 0) {
                ret = ep_interp_matrix->InsertMyValues(ho_row, rowNNZ, values_d.data(), indices_d.data());
                TEUCHOS_ASSERT(ret == 0);
              }
            } //end if owned
          } // end HO LID loop
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

                               // get IDs for HO and LO dofs
                               auto hoLIDs_d = Kokkos::subview(hoElementLIDs_d, elemId, Kokkos::ALL());
                               auto loLIDs_d = Kokkos::subview(loElementLIDs_d, elemId, Kokkos::ALL());

                               // loop over HO LIDs
                               for(size_t hoIter = 0; hoIter < hoCardinality; ++hoIter) {
                                 const LocalOrdinal ho_row = hoLIDs_d(hoIter);
                                 const bool isOwned = ((minLocalIndex <= ho_row) && (ho_row <= maxLocalIndex));
                                 if (isOwned) {
                                   // filter entries for zeros
                                   for(size_t loIter = 0; loIter < loCardinality; loIter++) {
                                     Scalar val = basisCoeffsLIOriented_d(cellNo, hoIter, loIter);
                                     if (KAT::magnitude(val) > entryFilterTol)
				       lcl_tp_interp_matrix.replaceValues(ho_row, &(loLIDs_d(loIter)), 1, &val, /*is_sorted=*/false, /*force_atomic=*/true);
                                   } //end if owned
                                 } // isOwned
                               } // end HO LID loop
                             }); //end workset loop
      } // Tpetra fill
    } //end element loop
  } //end element block loop

  if (tblof != Teuchos::null) {
    tp_interp_matrix->fillComplete(tp_domainmap, tp_rangemap);

    if (op != Intrepid2::OPERATOR_VALUE) {
      // Discrete gradient, curl, etc actually live on edges, faces, etc, so we have a lot of zeros in the matrix.
      tp_interp_matrix = removeSmallEntries(tp_interp_matrix, entryFilterTol);
    }

    thyra_interp = Thyra::tpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,typename tp_matrix::node_type>(Thyra::createVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal>(tp_rangemap),
                                                                                                          Thyra::createVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal>(tp_domainmap),
                                                                                                          tp_interp_matrix);
#if 0
    // compare the sparse matrix version and the matrix-free apply
    auto mfOp = rcp(new mini_em::MatrixFreeInterpolationOp<Scalar,LocalOrdinal,GlobalOrdinal>("test", linObjFactory, lo_basis_name, ho_basis_name, op, worksetSize));
    auto thyra_mfOp = Thyra::tpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,typename tp_matrix::node_type>(Thyra::createVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal>(mfOp->getRangeMap()),
                                                                                                             Thyra::createVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal>(mfOp->getDomainMap()),
                                                                                                             mfOp);
    {
      auto testX  = rcp(new Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal>(tp_domainmap, 1));
      auto testY1 = rcp(new Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal>(tp_rangemap, 1));
      auto testY2 = rcp(new Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal>(tp_rangemap, 1));
      testX->randomize();
      testY1->putScalar(1.);
      testY2->putScalar(1.);

      tp_interp_matrix->apply(*testX, *testY1, Teuchos::NO_TRANS, 3.0, 2.0);
      mfOp->apply(*testX, *testY2, Teuchos::NO_TRANS, 3.0, 2.0);
      testY1->update(-1.0,*testY2,1.0);
      std::cout << "norm difference for 3 * M * X + 2 Y: " << testY1->getVector(0)->norm2() << std::endl;

      tp_interp_matrix->apply(*testX, *testY1);
      mfOp->apply(*testX, *testY2);
      testY1->update(-1.0,*testY2,1.0);
      std::cout << "norm difference for M * X: " << testY1->getVector(0)->norm2() << std::endl;

      testX  = rcp(new Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal>(tp_rangemap, 1));
      testY1 = rcp(new Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal>(tp_domainmap, 1));
      testY2 = rcp(new Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal>(tp_domainmap, 1));
      testX->randomize();
      testY1->putScalar(1.);
      testY2->putScalar(1.);

      tp_interp_matrix->apply(*testX, *testY1, Teuchos::TRANS, 3.0, 2.0);
      mfOp->apply(*testX, *testY2, Teuchos::TRANS, 3.0, 2.0);
      testY1->update(-1.0,*testY2,1.0);
      std::cout << "norm difference for 3 * M^T * X + 2 Y: " << testY1->getVector(0)->norm2() << std::endl;

      tp_interp_matrix->apply(*testX, *testY1, Teuchos::TRANS);
      mfOp->apply(*testX, *testY2, Teuchos::TRANS);

      static int counter = 0;
      Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal> >::writeDenseFile("X_" + std::to_string(counter)+".mm", *testX);
      Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal> >::writeDenseFile("Y1_" + std::to_string(counter)+".mm", *testY1);
      Tpetra::MatrixMarket::Writer<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal> >::writeDenseFile("Y2_" + std::to_string(counter)+".mm", *testY2);
      ++counter;

      testY1->update(-1.0,*testY2,1.0);
      std::cout << "norm difference for M^T * X: " << testY1->getVector(0)->norm2() << std::endl;
    }

#endif

  }
#ifdef PANZER_HAVE_EPETRA_STACK
  else
    ep_interp_matrix->FillComplete(*ep_domainmap, *ep_rangemap);
#endif

  return thyra_interp;
}

void addInterpolationToRequestHandler(const std::string& name,
                                      const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > &linObjFactory,
                                      const Teuchos::RCP<Teko::RequestHandler> & reqHandler,
                                      const std::string& lo_basis_name,
                                      const std::string& ho_basis_name,
                                      Intrepid2::EOperator op,
                                      const bool waitForRequest,
                                      const bool dump,
                                      const size_t worksetSize,
                                      const bool matrixFree) {

  // add interpolation callback to request handler
  reqHandler->addRequestCallback(Teuchos::rcp(new InterpolationRequestCallback(name, linObjFactory, lo_basis_name, ho_basis_name, op, waitForRequest, dump, worksetSize, matrixFree)));
}


InterpolationRequestCallback::
InterpolationRequestCallback(const std::string& name,
                             const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > &linObjFactory,
                             const std::string& lo_basis_name,
                             const std::string& ho_basis_name,
                             Intrepid2::EOperator op,
                             const bool waitForRequest,
                             const bool dump,
                             const size_t worksetSize,
                             const bool matrixFree)
  : name_(name), linObjFactory_(linObjFactory), lo_basis_name_(lo_basis_name), ho_basis_name_(ho_basis_name), op_(op), dump_(dump), worksetSize_(worksetSize), matrixFree_(matrixFree)
{
  if (!waitForRequest)
    build();
}

void
InterpolationRequestCallback::
build() {
  if (!matrixFree_) {
    Teuchos::TimeMonitor tm(*Teuchos::TimeMonitor::getNewTimer(std::string("Mini-EM: assemble ") + name_));
    interp_ = buildInterpolation(linObjFactory_, lo_basis_name_, ho_basis_name_, op_, worksetSize_);
  } else {
    using Scalar = double;
    using LocalOrdinal = int;
    using GlobalOrdinal = panzer::GlobalOrdinal;
    using tp_matrix = Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal>;

    Teuchos::TimeMonitor tm(*Teuchos::TimeMonitor::getNewTimer(std::string("Mini-EM: matrix-free setup ") + name_));
    auto mfOp = rcp(new mini_em::MatrixFreeInterpolationOp<Scalar,LocalOrdinal,GlobalOrdinal>(name_, linObjFactory_, lo_basis_name_, ho_basis_name_, op_, worksetSize_));
    interp_ = Thyra::tpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,typename tp_matrix::node_type>(Thyra::createVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal>(mfOp->getRangeMap()),
                                                                                                     Thyra::createVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal>(mfOp->getDomainMap()),
                                                                                                     mfOp);
  }
  if (dump_ && !matrixFree_) {
    std::string filename = name_ + ".mm";

    mini_em::writeOut(filename, *interp_);
  }
}


bool
InterpolationRequestCallback::
handlesRequest(const Teko::RequestMesg & rm) {
  std::string name = rm.getName();

  return (name==name_);
}


Teko::LinearOp
InterpolationRequestCallback::
request(const Teko::RequestMesg & rm) {
  TEUCHOS_ASSERT(handlesRequest(rm));
  std::string name = rm.getName();

  if(name==name_) {
    if (interp_.is_null()) {
      build();
    }
    return interp_;
  } else
    TEUCHOS_ASSERT(false);
}

void
InterpolationRequestCallback::
preRequest(const Teko::RequestMesg & rm) {
    // checking for its existance is as good as pre requesting
    TEUCHOS_ASSERT(handlesRequest(rm));
  }
