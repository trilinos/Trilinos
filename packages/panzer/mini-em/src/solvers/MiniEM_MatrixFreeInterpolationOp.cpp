#include "MiniEM_MatrixFreeInterpolationOp.hpp"
#include <Tpetra_Operator.hpp>
#include "Panzer_LOCPair_GlobalEvaluationData.hpp"
#include "Panzer_IntrepidOrientation.hpp"
#include "Panzer_IntrepidBasisFactory.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_BlockedTpetraLinearObjFactory.hpp"
#include "Intrepid2_OrientationTools.hpp"
#include "Intrepid2_LagrangianInterpolation.hpp"
#include "MiniEM_Utils.hpp"

namespace mini_em {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  MatrixFreeInterpolationOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  MatrixFreeInterpolationOp(const std::string& _name,
                            Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > _linObjFactory,
                            const std::string& _lo_basis_name,
                            const std::string& _ho_basis_name,
                            Intrepid2::EOperator _op,
                            size_t _worksetSize) :
    name(_name),
    linObjFactory(_linObjFactory),
    lo_basis_name(_lo_basis_name),
    ho_basis_name(_ho_basis_name),
    op(_op),
    worksetSize(_worksetSize)
  {

    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_dynamic_cast;

    using OT  = Teuchos::OrdinalTraits<GlobalOrdinal>;

    typedef typename panzer::BlockedTpetraLinearObjFactory<panzer::Traits,Scalar,LocalOrdinal,GlobalOrdinal> tpetraBlockedLinObjFactory;
    typedef panzer::GlobalIndexer UGI;

    // must be able to cast to a block linear object factory
    RCP<const tpetraBlockedLinObjFactory > tblof = rcp_dynamic_cast<const tpetraBlockedLinObjFactory >(linObjFactory);

    // typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal> tp_matrix;
    typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal> tp_map;

    blockedDOFMngr = tblof->getGlobalIndexer();

    // get global indexers for LO and HO dofs
    std::vector<RCP<UGI> > fieldDOFMngrs = blockedDOFMngr->getFieldDOFManagers();
    int loFieldNum = blockedDOFMngr->getFieldNum(lo_basis_name);
    int hoFieldNum = blockedDOFMngr->getFieldNum(ho_basis_name);
    int loBlockIndex = blockedDOFMngr->getFieldBlock(loFieldNum);
    int hoBlockIndex = blockedDOFMngr->getFieldBlock(hoFieldNum);
    lo_ugi = rcp_dynamic_cast<panzer::DOFManager>(blockedDOFMngr->getFieldDOFManagers()[loBlockIndex],true);
    ho_ugi = rcp_dynamic_cast<panzer::DOFManager>(blockedDOFMngr->getFieldDOFManagers()[hoBlockIndex],true);

    // The operator maps from LO (domain) to HO (range)
    RCP<const tp_map> tp_rangemap, tp_domainmap, tp_rowmap, tp_colmap;

    {
      // build maps
      std::vector<GlobalOrdinal> gids;
      ho_ugi->getOwnedIndices(gids);
      tp_rowmap = rcp(new tp_map(OT::invalid(), gids.data(), gids.size(), OT::zero(), ho_ugi->getComm()));
      tp_rangemap = tp_rowmap;
      lo_ugi->getOwnedIndices(gids);
      tp_domainmap = rcp(new tp_map(OT::invalid(), gids.data(), gids.size(), OT::zero(), lo_ugi->getComm()));
      lo_ugi->getOwnedAndGhostedIndices(gids);
      tp_colmap = rcp(new tp_map(OT::invalid(), gids.data(), gids.size(), OT::zero(), lo_ugi->getComm()));
    }

    domainMap_ = tp_domainmap;
    rangeMap_ = tp_rangemap;
    columnMap_ = tp_colmap;
    import_ = rcp(new Tpetra::Import<LocalOrdinal,GlobalOrdinal,Node>(domainMap_, columnMap_));

    allocateColumnMapVector(1);

    precomputeOwners();

    setupNodeOnlyConnManager();

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
  precomputeOwners() {
    size_t maxNumElementsPerBlock = 0;
    int numCells;
    if (maxNumElementsPerBlock > 0)
      numCells = std::min(maxNumElementsPerBlock, worksetSize);
    else
      numCells = worksetSize;

    LocalOrdinal lclTargetSize = getRangeMap()->getLocalNumElements();
    owner_d_ = Kokkos::View<LocalOrdinal*,DeviceSpace>("owner", lclTargetSize);

    auto hoLIDs_d = ho_ugi->getLIDs();

    auto owner_d = owner_d_;

    // loop over element blocks
    std::vector<std::string> elementBlockIds;
    blockedDOFMngr->getElementBlockIds(elementBlockIds);
    for(std::size_t blockIter = 0; blockIter < elementBlockIds.size(); ++blockIter) {

      // loop over element worksets
      std::vector<int> elementIds = ho_ugi->getElementBlock(elementBlockIds[blockIter]);
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

                               // loop over HO LIDs
                               for(size_t hoIter = 0; hoIter < hoLIDs_d.extent(1); ++hoIter) {
                                 LocalOrdinal ho_row = hoLIDs_d(elemId, hoIter);
                                 if (ho_row < lclTargetSize)
                                   owner_d(ho_row) = elemId;
                               }
                             });
      }
    }
  }


  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void
  MatrixFreeInterpolationOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  setupNodeOnlyConnManager() {
    Teuchos::RCP<const panzer::ConnManager> conn = blockedDOFMngr->getConnManager();

    // assume only one cell toplogy per rank
    std::vector<shards::CellTopology> topologies;
    conn->getElementBlockTopologies(topologies);
    shards::CellTopology topology = topologies[0];
    // set up a node only conn manager
    auto node_basis = panzer::createIntrepid2Basis<DeviceSpace,Scalar,Scalar>("HGrad",1,topology);
    auto node_fieldPattern = Teuchos::rcp(new panzer::Intrepid2FieldPattern(node_basis));
    node_conn_ = Teuchos::rcp_dynamic_cast<panzer_stk::STKConnManager>(conn->noConnectivityClone(),true);
    node_conn_->buildConnectivity(*node_fieldPattern);
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

    typedef Intrepid2::OrientationTools<DeviceSpace> ots;
    typedef Intrepid2::Experimental::LagrangianInterpolation<DeviceSpace> li;
    typedef Kokkos::DynRankView<Scalar,DeviceSpace> DynRankDeviceView;

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

    // get the LO and HO bases
    auto lo_fieldPattern = lo_ugi->getFieldPattern(lo_basis_name);
    auto lo_basis = rcp_dynamic_cast<const panzer::Intrepid2FieldPattern>(lo_fieldPattern,true)->getIntrepidBasis();
    auto ho_fieldPattern = ho_ugi->getFieldPattern(ho_basis_name);
    auto ho_basis = rcp_dynamic_cast<const panzer::Intrepid2FieldPattern>(ho_fieldPattern,true)->getIntrepidBasis();

    // cardinalities
    const size_t loCardinality = lo_basis->getCardinality();
    const size_t hoCardinality = ho_basis->getCardinality();

    size_t maxNumElementsPerBlock = 0;

    RCP<const panzer::ConnManager> conn = blockedDOFMngr->getConnManager();

    // assume only one cell toplogy per rank
    std::vector<shards::CellTopology> topologies;
    conn->getElementBlockTopologies(topologies);
    shards::CellTopology topology = topologies[0];
    int dim = topology.getDimension();
    // num vertices in an element
    const int numElemVertices = topology.getVertexCount();

    // // set up a node only conn manager
    // auto node_basis = panzer::createIntrepid2Basis<DeviceSpace,Scalar,Scalar>("HGrad",1,topology);
    // auto node_fieldPattern = rcp(new panzer::Intrepid2FieldPattern(node_basis));
    // RCP<panzer_stk::STKConnManager> node_conn = Teuchos::rcp_dynamic_cast<panzer_stk::STKConnManager>(conn->noConnectivityClone(),true);
    // node_conn->buildConnectivity(*node_fieldPattern);

    if (op == Intrepid2::OPERATOR_VALUE) {
      TEUCHOS_ASSERT(hoCardinality >= loCardinality);
      TEUCHOS_ASSERT_EQUALITY(lo_basis->getFunctionSpace(), ho_basis->getFunctionSpace());
    }

    // allocate some views
    int numCells;
    if (maxNumElementsPerBlock > 0)
      numCells = std::min(maxNumElementsPerBlock, worksetSize);
    else
      numCells = worksetSize;
    DynRankDeviceView ho_dofCoords_d("ho_dofCoords_d", numCells, hoCardinality, dim);
    DynRankDeviceView basisCoeffsLIOriented_d("basisCoeffsLIOriented_d", numCells, hoCardinality, loCardinality);

    typename Kokkos::DynRankView<Intrepid2::Orientation,DeviceSpace> elemOrts_d ("elemOrts_d",  numCells);
    typename Kokkos::DynRankView<GlobalOrdinal, DeviceSpace>         elemNodes_d("elemNodes_d", numCells, numElemVertices);

    // the ranks of these depend on dimension
    DynRankDeviceView ho_dofCoeffs_d;
    DynRankDeviceView valuesAtDofCoordsNonOriented_d;
    DynRankDeviceView valuesAtDofCoordsOriented_d;
    DynRankDeviceView reducedValuesAtDofCoordsOriented_d;

    {
      // Let Intrepid2 give us the correctly dimensioned view, then build one with +1 ranks and extent(0) == numCells
      auto temp = lo_basis->allocateOutputView(hoCardinality, op);

      // These view have dimensions
      //  numCells, numFields=loCardinality, numPoints=hoCardinality, (spatialDim)
      //
      if (temp.rank() == 3) {
        valuesAtDofCoordsNonOriented_d     = DynRankDeviceView("valuesAtDofCoordsNonOriented_d", numCells, temp.extent(0), temp.extent(1), temp.extent(2));
        valuesAtDofCoordsOriented_d        = DynRankDeviceView("valuesAtDofCoordsOriented_d", numCells, temp.extent(0), temp.extent(1), temp.extent(2));
        reducedValuesAtDofCoordsOriented_d = DynRankDeviceView("reducedValuesAtDofCoordsOriented_d", numCells, temp.extent(1), temp.extent(2), numVectors);
      } else {
        valuesAtDofCoordsNonOriented_d     = DynRankDeviceView("valuesAtDofCoordsNonOriented_d", numCells, temp.extent(0), temp.extent(1));
        valuesAtDofCoordsOriented_d        = DynRankDeviceView("valuesAtDofCoordsOriented_d", numCells, temp.extent(0), temp.extent(1));
        reducedValuesAtDofCoordsOriented_d = DynRankDeviceView("reducedValuesAtDofCoordsOriented_d", numCells, temp.extent(1), 1, numVectors);
      }
    }

    int fieldRank = Intrepid2::getFieldRank(ho_basis->getFunctionSpace());
    TEUCHOS_ASSERT((fieldRank == 0) || (fieldRank == 1));

    const bool isVectorBasis = (fieldRank == 1);
    if (isVectorBasis) {
      ho_dofCoeffs_d = DynRankDeviceView("ho_dofCoeffs_d", numCells, hoCardinality, dim);
    } else {
      ho_dofCoeffs_d = DynRankDeviceView("ho_dofCoeffs_d", numCells, hoCardinality);
    }

    auto hoLIDs_d = ho_ugi->getLIDs();
    auto loLIDs_d = lo_ugi->getLIDs();

    auto node_connectivity_h = node_conn_->getConnectivityView();
    Kokkos::View<GlobalOrdinal*,DeviceSpace> node_connectivity_d("node_connectivity_d", node_connectivity_h.extent(0));
    Kokkos::deep_copy(node_connectivity_d, node_connectivity_h);
    auto node_connectivitySize_h = node_conn_->getConnectivitySizeView();
    Kokkos::View<LocalOrdinal*,DeviceSpace> node_connectivitySize_d("node_connectivitySize_d", node_connectivitySize_h.extent(0));
    Kokkos::deep_copy(node_connectivitySize_d, node_connectivitySize_h);
    auto node_elementLidToConn_h = node_conn_->getElementLidToConnView();
    Kokkos::View<LocalOrdinal*,DeviceSpace> node_elementLidToConn_d("node_elementLidToConn_d", node_elementLidToConn_h.extent(0));
    Kokkos::deep_copy(node_elementLidToConn_d, node_elementLidToConn_h);

    // loop over element blocks
    std::vector<std::string> elementBlockIds;
    blockedDOFMngr->getElementBlockIds(elementBlockIds);
    for(std::size_t blockIter = 0; blockIter < elementBlockIds.size(); ++blockIter) {

      // loop over element worksets
      std::vector<int> elementIds = ho_ugi->getElementBlock(elementBlockIds[blockIter]);
      Kokkos::View<int*,DeviceSpace>::HostMirror elementIds_h(elementIds.data(), elementIds.size());
      Kokkos::View<int*,DeviceSpace> elementIds_d("elementIds_d", elementIds_h.extent(0));
      Kokkos::deep_copy(elementIds_d, elementIds_h);
      for(std::size_t elemIter = 0; elemIter < elementIds_d.extent(0); elemIter += numCells) {

        // get element orientations
        Kokkos::parallel_for("miniEM:MatrixFreeInterpolationOp::connectivity",
                             range_type(0, std::min(numCells, elementIds_d.extent_int(0)-Teuchos::as<int>(elemIter))),
                             KOKKOS_LAMBDA(const LocalOrdinal cellNo) {
                               LocalOrdinal cellNo2 = elemIter+cellNo;
                               LocalOrdinal elementID = elementIds_d(cellNo2);
                               LocalOrdinal k = node_elementLidToConn_d(elementID);
                               for(int i = 0; i < node_connectivitySize_d(elementID); i++)
                                 elemNodes_d(cellNo, i) = node_connectivity_d(k+i);
                             });
        ots::getOrientation(elemOrts_d, elemNodes_d, topology);

        // HO dof coordinates and coefficients
        li::getDofCoordsAndCoeffs(ho_dofCoords_d, ho_dofCoeffs_d, ho_basis.get(), elemOrts_d);

        // compute values of op * (LO basis) at HO dof coords on reference element
        // TODO: Once this is supported by Intrepid2, make this a parallel_for.
        for (int cellNo = 0; cellNo < numCells; cellNo++)
          lo_basis->getValues(Kokkos::subview(valuesAtDofCoordsNonOriented_d, cellNo, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL()),
                              Kokkos::subview(ho_dofCoords_d, cellNo, Kokkos::ALL(), Kokkos::ALL()),
                              op);

        // apply orientations for LO basis
        // shuffles things in the second dimension, i.e. wrt LO basis
        ots::modifyBasisByOrientation(valuesAtDofCoordsOriented_d,
                                      valuesAtDofCoordsNonOriented_d,
                                      elemOrts_d,
                                      lo_basis.get());

        Kokkos::deep_copy(reducedValuesAtDofCoordsOriented_d, 0.0);
        Kokkos::parallel_for("miniEM:MatrixFreeInterpolationOp:cellLoop1",
                             range_type(0, std::min(numCells, elementIds_d.extent_int(0)-Teuchos::as<int>(elemIter))),
                             KOKKOS_LAMBDA(const LocalOrdinal cellNo) {
                               LocalOrdinal cellNo2 = elemIter+cellNo;
                               LocalOrdinal elemId = elementIds_d(cellNo2);
                               for(size_t loIter=0; loIter<loCardinality; loIter++) {
                                 LocalOrdinal J = loLIDs_d(elemId, loIter);
                                 for(size_t hoIter=0; hoIter<hoCardinality; hoIter++) {
                                   for(size_t d=0; d<valuesAtDofCoordsOriented_d.extent(3); d++) {
                                     Scalar val = valuesAtDofCoordsOriented_d(cellNo, loIter, hoIter, d);
                                     for (size_t j = 0; j<numVectors; ++j)
                                       reducedValuesAtDofCoordsOriented_d(cellNo, hoIter, d, j) += val * lclX(J, j);
                                   }
                                 }
                               }
                             });


        for (size_t j = 0; j<numVectors; ++j)
          li::getBasisCoeffs(Kokkos::subview(basisCoeffsLIOriented_d, Kokkos::ALL(), Kokkos::ALL(), j),
                             Kokkos::subview(reducedValuesAtDofCoordsOriented_d, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL(), j),
                             ho_dofCoeffs_d);

        auto owner_d = owner_d_;

        Kokkos::parallel_for("miniEM::MatrixFreeInterpolationOp::cellLoop2",
                             range_type(elemIter, std::min(elemIter+numCells,
                                                           elementIds_d.extent(0))),
                             KOKKOS_LAMBDA(const LocalOrdinal cellNo2) {
                               LocalOrdinal cellNo = cellNo2-elemIter;
                               LocalOrdinal elemId = elementIds_d(cellNo2);

                               // loop over HO LIDs
                               for(size_t hoIter = 0; hoIter < hoLIDs_d.extent(1); ++hoIter) {
                                 LocalOrdinal ho_row = hoLIDs_d(elemId, hoIter);

                                 // if owned
                                 if ((ho_row < lclTargetSize) && (owner_d(ho_row) == elemId)) {

                                   for (size_t j = 0; j<numVectors; ++j) {
                                     Scalar val = basisCoeffsLIOriented_d(cellNo, hoIter, j);
                                     lclY(ho_row,j) = beta*lclY(ho_row,j) + alpha*val;
                                   }
                                 } // end if owned
                               } // end HO LID loop
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
    typedef Intrepid2::Experimental::LagrangianInterpolation<DeviceSpace> li;
    typedef Kokkos::DynRankView<Scalar,DeviceSpace> DynRankDeviceView;

    using TST = Teuchos::ScalarTraits<Scalar>;
    const Scalar ZERO = TST::zero();
    using view_t = typename Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::dual_view_type::t_dev;
    using const_view_t = typename Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::dual_view_type::t_dev::const_type;

    Teuchos::TimeMonitor tm(*Teuchos::TimeMonitor::getNewTimer(std::string("Mini-EM: matrix-free apply trans ") + name));

    const_view_t lclX = X.getLocalViewDevice(Tpetra::Access::ReadOnly);
    colmapMV_->putScalar(ZERO);
    view_t lclYtemp = colmapMV_->getLocalViewDevice(Tpetra::Access::ReadWrite);

    // get the LO and HO bases
    auto lo_fieldPattern = lo_ugi->getFieldPattern(lo_basis_name);
    auto lo_basis = rcp_dynamic_cast<const panzer::Intrepid2FieldPattern>(lo_fieldPattern,true)->getIntrepidBasis();
    auto ho_fieldPattern = ho_ugi->getFieldPattern(ho_basis_name);
    auto ho_basis = rcp_dynamic_cast<const panzer::Intrepid2FieldPattern>(ho_fieldPattern,true)->getIntrepidBasis();

    // cardinalities
    const size_t loCardinality = lo_basis->getCardinality();
    const size_t hoCardinality = ho_basis->getCardinality();

    size_t maxNumElementsPerBlock = 0;

    RCP<const panzer::ConnManager> conn = blockedDOFMngr->getConnManager();

    // assume only one cell toplogy per rank
    std::vector<shards::CellTopology> topologies;
    conn->getElementBlockTopologies(topologies);
    shards::CellTopology topology = topologies[0];
    int dim = topology.getDimension();
    // num vertices in an element
    const int numElemVertices = topology.getVertexCount();

    // // set up a node only conn manager
    // auto node_basis = panzer::createIntrepid2Basis<DeviceSpace,Scalar,Scalar>("HGrad",1,topology);
    // auto node_fieldPattern = rcp(new panzer::Intrepid2FieldPattern(node_basis));
    // RCP<panzer_stk::STKConnManager> node_conn = Teuchos::rcp_dynamic_cast<panzer_stk::STKConnManager>(conn->noConnectivityClone(),true);
    // node_conn->buildConnectivity(*node_fieldPattern);

    if (op == Intrepid2::OPERATOR_VALUE) {
      TEUCHOS_ASSERT(hoCardinality >= loCardinality);
      TEUCHOS_ASSERT_EQUALITY(lo_basis->getFunctionSpace(), ho_basis->getFunctionSpace());
    }

    // allocate some views
    int numCells;
    if (maxNumElementsPerBlock > 0)
      numCells = std::min(maxNumElementsPerBlock, worksetSize);
    else
      numCells = worksetSize;
    DynRankDeviceView ho_dofCoords_d("ho_dofCoords_d", numCells, hoCardinality, dim);
    DynRankDeviceView basisCoeffsLIOriented_d("basisCoeffsLIOriented_d", numCells, hoCardinality, loCardinality);

    typename Kokkos::DynRankView<Intrepid2::Orientation,DeviceSpace> elemOrts_d ("elemOrts_d",  numCells);
    typename Kokkos::DynRankView<GlobalOrdinal, DeviceSpace>         elemNodes_d("elemNodes_d", numCells, numElemVertices);

    // the ranks of these depend on dimension
    DynRankDeviceView ho_dofCoeffs_d;
    DynRankDeviceView valuesAtDofCoordsNonOriented_d;
    DynRankDeviceView valuesAtDofCoordsOriented_d;

    {
      // Let Intrepid2 give us the correctly dimensioned view, then build one with +1 ranks and extent(0) == numCells
      auto temp = lo_basis->allocateOutputView(hoCardinality, op);

      // These view have dimensions
      //  numCells, numFields=loCardinality, numPoints=hoCardinality, (spatialDim)
      //
      if (temp.rank() == 3) {
        valuesAtDofCoordsNonOriented_d = DynRankDeviceView("valuesAtDofCoordsNonOriented_d", numCells, temp.extent(0), temp.extent(1), temp.extent(2));
        valuesAtDofCoordsOriented_d    = DynRankDeviceView("valuesAtDofCoordsNonOriented_d", numCells, temp.extent(0), temp.extent(1), temp.extent(2));
      } else {
        valuesAtDofCoordsNonOriented_d = DynRankDeviceView("valuesAtDofCoordsNonOriented_d", numCells, temp.extent(0), temp.extent(1));
        valuesAtDofCoordsOriented_d    = DynRankDeviceView("valuesAtDofCoordsNonOriented_d", numCells, temp.extent(0), temp.extent(1));
      }
    }

    int fieldRank = Intrepid2::getFieldRank(ho_basis->getFunctionSpace());
    TEUCHOS_ASSERT((fieldRank == 0) || (fieldRank == 1));

    const bool isVectorBasis = (fieldRank == 1);
    if (isVectorBasis) {
      ho_dofCoeffs_d = DynRankDeviceView("ho_dofCoeffs_d", numCells, hoCardinality, dim);
    } else {
      ho_dofCoeffs_d = DynRankDeviceView("ho_dofCoeffs_d", numCells, hoCardinality);
    }

    auto hoLIDs_d = ho_ugi->getLIDs();
    auto loLIDs_d = lo_ugi->getLIDs();
    Kokkos::fence();

    // loop over element blocks
    std::vector<std::string> elementBlockIds;
    blockedDOFMngr->getElementBlockIds(elementBlockIds);
    for(std::size_t blockIter = 0; blockIter < elementBlockIds.size(); ++blockIter) {

      // loop over element worksets
      std::vector<int> elementIds = ho_ugi->getElementBlock(elementBlockIds[blockIter]);
      Kokkos::View<int*,DeviceSpace>::HostMirror elementIds_h(elementIds.data(), elementIds.size());
      Kokkos::View<int*,DeviceSpace> elementIds_d("elementIds_d", elementIds_h.extent(0));
      Kokkos::deep_copy(elementIds_d, elementIds_h);
      Kokkos::fence();
      for(std::size_t elemIter = 0; elemIter < elementIds_d.extent(0); elemIter += numCells) {

        // get element orientations
        auto node_connectivity_h = node_conn_->getConnectivityView();
        Kokkos::View<GlobalOrdinal*,DeviceSpace> node_connectivity_d("node_connectivity_d", node_connectivity_h.extent(0));
        Kokkos::deep_copy(node_connectivity_d, node_connectivity_h);
        auto node_connectivitySize_h = node_conn_->getConnectivitySizeView();
        Kokkos::View<LocalOrdinal*,DeviceSpace> node_connectivitySize_d("node_connectivitySize_d", node_connectivitySize_h.extent(0));
        Kokkos::deep_copy(node_connectivitySize_d, node_connectivitySize_h);
        auto node_elementLidToConn_h = node_conn_->getElementLidToConnView();
        Kokkos::View<LocalOrdinal*,DeviceSpace> node_elementLidToConn_d("node_elementLidToConn_d", node_elementLidToConn_h.extent(0));
        Kokkos::deep_copy(node_elementLidToConn_d, node_elementLidToConn_h);
        Kokkos::fence();

        Kokkos::parallel_for("miniEM:MatrixFreeInterpolationOp::connectivity",
                             range_type(0, std::min(numCells, elementIds_d.extent_int(0)-Teuchos::as<int>(elemIter))),
                             KOKKOS_LAMBDA(const LocalOrdinal cellNo) {
                               LocalOrdinal cellNo2 = elemIter+cellNo;
                               LocalOrdinal elementID = elementIds_d(cellNo2);
                               LocalOrdinal k = node_elementLidToConn_d(elementID);
                               for(int i = 0; i < node_connectivitySize_d(elementID); i++)
                                 elemNodes_d(cellNo, i) = node_connectivity_d(k+i);
                             });
        Kokkos::fence();
        ots::getOrientation(elemOrts_d, elemNodes_d, topology);
        Kokkos::fence();

        // HO dof coordinates and coefficients
        li::getDofCoordsAndCoeffs(ho_dofCoords_d, ho_dofCoeffs_d, ho_basis.get(), elemOrts_d);
        Kokkos::fence();

        // compute values of op * (LO basis) at HO dof coords on reference element
        // TODO: Once this is supported by Intrepid2, make this a parallel_for.
        for (int cellNo = 0; cellNo < numCells; cellNo++)
          lo_basis->getValues(Kokkos::subview(valuesAtDofCoordsNonOriented_d, cellNo, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL()),
                              Kokkos::subview(ho_dofCoords_d, cellNo, Kokkos::ALL(), Kokkos::ALL()),
                              op);
        Kokkos::fence();

        // apply orientations for LO basis
        // shuffles things in the second dimension, i.e. wrt LO basis
        ots::modifyBasisByOrientation(valuesAtDofCoordsOriented_d,
                                      valuesAtDofCoordsNonOriented_d,
                                      elemOrts_d,
                                      lo_basis.get());
        Kokkos::fence();

        // get basis coefficients of LO basis functions wrt HO basis
        for(size_t loIter=0; loIter<loCardinality; loIter++)
          // Get basis coeffs wrt HO basis on reference element.
          // basisCoeffsLI has dimensions (numCells, numFields=hoCardinality, loCardinality)
          li::getBasisCoeffs(Kokkos::subview(basisCoeffsLIOriented_d, Kokkos::ALL(), Kokkos::ALL(), loIter),
                             Kokkos::subview(valuesAtDofCoordsOriented_d, Kokkos::ALL(), loIter, Kokkos::ALL(), Kokkos::ALL()),
                             ho_dofCoeffs_d);
        Kokkos::fence();

        auto owner_d = owner_d_;


        Kokkos::parallel_for("miniEM::MatrixFreeInterpolationOp::cellLoop",
                             range_type(elemIter, std::min(elemIter+numCells,
                                                           elementIds_d.extent(0))),
                             KOKKOS_LAMBDA(const LocalOrdinal cellNo2) {
                               LocalOrdinal cellNo = cellNo2-elemIter;
                               LocalOrdinal elemId = elementIds_d(cellNo2);

                               // loop over HO LIDs
                               for(size_t hoIter = 0; hoIter < hoLIDs_d.extent(1); ++hoIter) {
                                 LocalOrdinal ho_row = hoLIDs_d(elemId, hoIter);

                                 // if owned
                                 if ((ho_row < (LocalOrdinal) lclX.extent(0)) && (owner_d(ho_row) == elemId)) {

                                   for(size_t loIter = 0; loIter < loLIDs_d.extent(1); loIter++) {
                                     LocalOrdinal J = loLIDs_d(elemId, loIter);
                                     Scalar val = basisCoeffsLIOriented_d(cellNo, hoIter, loIter);
                                     for (size_t j = 0; j<lclYtemp.extent(1); ++j)
                                       Kokkos::atomic_add(&lclYtemp(J,j), alpha*val*lclX(ho_row,j));
                                   }
                                 } //end if owned
                               } //end HO LID loop
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
