#ifndef _MiniEM_Interpolation_hpp_
#define _MiniEM_Interpolation_hpp_


#include "Panzer_LOCPair_GlobalEvaluationData.hpp"
#include "Panzer_IntrepidOrientation.hpp"
#include "Panzer_IntrepidBasisFactory.hpp"
#include "Intrepid2_OrientationTools.hpp"
#include "Intrepid2_LagrangianInterpolation.hpp"
#ifdef PANZER_HAVE_EPETRA_STACK
#include "Thyra_EpetraThyraWrappers.hpp"
#endif
#include "MiniEM_Utils.hpp"
#include "MiniEM_MatrixFreeInterpolationOp.hpp"
#include "MiniEM_MatrixFreeInterpolationOp.cpp"


Teko::LinearOp buildInterpolation(const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > linObjFactory,
                                  const std::string& lo_basis_name,
                                  const std::string& ho_basis_name,
                                  Intrepid2::EOperator op=Intrepid2::OPERATOR_VALUE,
                                  size_t worksetSize=1000)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;

  typedef double Scalar;
  typedef int LocalOrdinal;
  typedef panzer::GlobalOrdinal GlobalOrdinal;

  using STS = Teuchos::ScalarTraits<Scalar>;
  using OT  = Teuchos::OrdinalTraits<GlobalOrdinal>;

  typedef typename panzer::BlockedTpetraLinearObjFactory<panzer::Traits,Scalar,LocalOrdinal,GlobalOrdinal> tpetraBlockedLinObjFactory;
#ifdef PANZER_HAVE_EPETRA_STACK
  typedef typename panzer::BlockedEpetraLinearObjFactory<panzer::Traits,LocalOrdinal> epetraBlockedLinObjFactory;
#endif
  typedef panzer::GlobalIndexer UGI;
  typedef PHX::Device DeviceSpace;
  typedef Kokkos::HostSpace HostSpace;
  typedef Intrepid2::OrientationTools<DeviceSpace> ots;
  typedef Intrepid2::Experimental::LagrangianInterpolation<DeviceSpace> li;

  typedef Kokkos::DynRankView<double,DeviceSpace> DynRankDeviceView;

  // must be able to cast to a block linear object factory
  RCP<const tpetraBlockedLinObjFactory > tblof = rcp_dynamic_cast<const tpetraBlockedLinObjFactory >(linObjFactory);
#ifdef PANZER_HAVE_EPETRA_STACK
  RCP<const epetraBlockedLinObjFactory > eblof = rcp_dynamic_cast<const epetraBlockedLinObjFactory >(linObjFactory);
#endif

  typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal> tp_matrix;
  typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal> tp_map;
#ifdef PANZER_HAVE_EPETRA_STACK
  typedef typename panzer::BlockedEpetraLinearObjContainer ep_linObjContainer;
  typedef Epetra_CrsMatrix ep_matrix;
  typedef Epetra_Map ep_map;
#endif

  RCP<const panzer::BlockedDOFManager> blockedDOFMngr;
  if (tblof != Teuchos::null) {
    blockedDOFMngr = tblof->getGlobalIndexer();
#ifdef PANZER_HAVE_EPETRA_STACK
  } else if (eblof != Teuchos::null) {
    TEUCHOS_ASSERT(false);
    // The Epetra code path works, expect for the fact that Epetra
    // does not implement the needed matrix entry insertion. We'd need
    // to handle the overwriting (instead of summing into) of already
    // existing entries by hand.
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
  RCP<const tp_map> tp_rangemap, tp_domainmap, tp_rowmap, tp_colmap;
  RCP<tp_matrix> tp_interp_matrix;
#ifdef PANZER_HAVE_EPETRA_STACK
  RCP<const ep_map> ep_rangemap, ep_domainmap, ep_rowmap, ep_colmap;
  RCP<ep_matrix> ep_interp_matrix;
#endif

  RCP<Thyra::LinearOpBase<Scalar> > thyra_interp;
  if (tblof != Teuchos::null) {
    // build maps
    std::vector<GlobalOrdinal> gids;
    ho_ugi->getOwnedIndices(gids);
    tp_rowmap = rcp(new tp_map(OT::invalid(), gids.data(), gids.size(), OT::zero(), ho_ugi->getComm()));
    tp_rangemap = tp_rowmap;
    lo_ugi->getOwnedIndices(gids);
    tp_domainmap = rcp(new tp_map(OT::invalid(), gids.data(), gids.size(), OT::zero(), lo_ugi->getComm()));
    lo_ugi->getOwnedAndGhostedIndices(gids);
    tp_colmap = rcp(new tp_map(OT::invalid(), gids.data(), gids.size(), OT::zero(), lo_ugi->getComm()));

    // estimate number of entries per row
    // This is an upper bound, as we are counting dofs that are on shared nodes, edges, faces more than once.
    Kokkos::View<size_t*,HostSpace> numEntriesPerRow("numEntriesPerRow", tp_rowmap->getLocalNumElements());
    {
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
          auto hoLIDs_d = ho_ugi->getElementLIDs(elemId);
          Kokkos::deep_copy(hoLIDs_h, hoLIDs_d);

          // TODO: do counting on device

          // loop over HO LIDs
          for(size_t hoIter = 0; hoIter < hoLIDs_h.size(); ++hoIter) {
            const LocalOrdinal ho_row = hoLIDs_h(hoIter);
            const bool isOwned = tp_rowmap->isNodeLocalElement(ho_row);
            if (isOwned)
              numEntriesPerRow(ho_row) += loCardinality;
          } //end HO LID loop
        } // elements loop
      } // blocks loop
    }

    Teuchos::ArrayView<size_t> nEPR = Teuchos::ArrayView<size_t>(numEntriesPerRow.data(), numEntriesPerRow.extent(0));
    tp_interp_matrix = rcp(new tp_matrix(tp_rowmap,tp_colmap,nEPR));
    RCP<const Thyra::VectorSpaceBase<Scalar> > rangeVectorSpace = Thyra::createVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal>(tp_rangemap);
    RCP<const Thyra::VectorSpaceBase<Scalar> > domainVectorSpace = Thyra::createVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal>(tp_domainmap);
    thyra_interp = Thyra::tpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,typename tp_matrix::node_type>(rangeVectorSpace,
                                                                                                          domainVectorSpace,
                                                                                                          tp_interp_matrix);
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

    // TODO: Fix this.
    size_t nnzPerRowEstimate = 25*loCardinality;

    ep_interp_matrix = rcp(new ep_matrix(Copy, *ep_rowmap, *ep_colmap, nnzPerRowEstimate, /*StaticProfile=*/true));

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
  int dim = topology.getDimension();
  // num vertices in an element
  const int numElemVertices = topology.getVertexCount();

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
    numCells = std::min(maxNumElementsPerBlock, worksetSize);
  else
    numCells = worksetSize;
  DynRankDeviceView             ho_dofCoords_d("ho_dofCoords_d", numCells, hoCardinality, dim);
  DynRankDeviceView             basisCoeffsLIOriented_d("basisCoeffsLIOriented_d", numCells, hoCardinality, loCardinality);
  DynRankDeviceView::HostMirror basisCoeffsLIOriented_h = Kokkos::create_mirror_view(basisCoeffsLIOriented_d);
  typename Kokkos::DynRankView<Intrepid2::Orientation,DeviceSpace>     elemOrts_d ("elemOrts_d",  numCells);
  typename Kokkos::DynRankView<GlobalOrdinal, DeviceSpace>::HostMirror elemNodes_h("elemNodes_h", numCells, numElemVertices);
  typename Kokkos::DynRankView<GlobalOrdinal, DeviceSpace>             elemNodes_d("elemNodes_d", numCells, numElemVertices);

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

  Kokkos::View<LocalOrdinal*,HostSpace> indices_h("indices", loCardinality);
  Kokkos::View<Scalar*,      HostSpace> values_h ("values",  loCardinality);

  auto entryFilterTol = 100*Teuchos::ScalarTraits<typename STS::magnitudeType>::eps();

  // loop over element blocks
  std::vector<std::string> elementBlockIds;
  blockedDOFMngr->getElementBlockIds(elementBlockIds);
  for(std::size_t blockIter = 0; blockIter < elementBlockIds.size(); ++blockIter) {

    // loop over element worksets
    std::vector<int> elementIds = ho_ugi->getElementBlock(elementBlockIds[blockIter]);
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

      //get basis coefficients of LO basis functions wrt HO basis
      for(size_t loIter=0; loIter<loCardinality; loIter++)
        // Get basis coeffs wrt HO basis on reference element.
        // basisCoeffsLI has dimensions (numCells, numFields=hoCardinality, loCardinality)
        li::getBasisCoeffs(Kokkos::subview(basisCoeffsLIOriented_d, Kokkos::ALL(), Kokkos::ALL(), loIter),
                           Kokkos::subview(valuesAtDofCoordsOriented_d, Kokkos::ALL(), loIter, Kokkos::ALL(), Kokkos::ALL()),
                           ho_dofCoeffs_d);

      Kokkos::deep_copy(basisCoeffsLIOriented_h, basisCoeffsLIOriented_d);

      for (int cellNo = 0; cellNo < numCells; cellNo++) {
        if (elemIter+cellNo >= elementIds.size())
          continue;
        auto elemId = elementIds[elemIter+cellNo];

        // get IDs for HO and LO dofs
        auto hoLIDs_d = ho_ugi->getElementLIDs(elemId);
        auto loLIDs_d = lo_ugi->getElementLIDs(elemId);
        Kokkos::deep_copy(hoLIDs_h, hoLIDs_d);
        Kokkos::deep_copy(loLIDs_h, loLIDs_d);

        // TODO: do filtering and insert on device

        // loop over HO LIDs
        for(size_t hoIter = 0; hoIter < hoLIDs_h.size(); ++hoIter) {
          LocalOrdinal ho_row = hoLIDs_h(hoIter);
          bool isOwned;
#ifdef PANZER_HAVE_EPETRA_STACK
          if (tblof != Teuchos::null)
            isOwned = tp_rowmap->isNodeLocalElement(ho_row);
          else
            isOwned = ep_rowmap->MyLID(ho_row);
#else
          isOwned = tp_rowmap->isNodeLocalElement(ho_row);
#endif

          if (isOwned) {
            // filter entries for zeros
            size_t rowNNZ = 0;
            for(size_t loIter = 0; loIter < loCardinality; loIter++) {
              Scalar val = basisCoeffsLIOriented_h(cellNo, hoIter, loIter);
              if (STS::magnitude(val) > entryFilterTol) {
                indices_h(rowNNZ) = loLIDs_h(loIter);
                values_h(rowNNZ) = val;
                rowNNZ += 1;
              }
            }

#ifdef PANZER_HAVE_EPETRA_STACK
            if (tblof != Teuchos::null)
              tp_interp_matrix->insertLocalValues(ho_row, rowNNZ, values_h.data(), indices_h.data(), Tpetra::INSERT);
            else
              ep_interp_matrix->InsertMyValues(ho_row, rowNNZ, values_h.data(), indices_h.data());
#else
            tp_interp_matrix->insertLocalValues(ho_row, rowNNZ, values_h.data(), indices_h.data(), Tpetra::INSERT);
#endif
          } //end if owned
        } //end HO LID loop
      } //end workset loop
    } //end element loop
  } //end element block loop


  if (tblof != Teuchos::null) {
    tp_interp_matrix->fillComplete(tp_domainmap, tp_rangemap);

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


class InterpolationRequestCallback : public Teko::RequestCallback<Teko::LinearOp> {
private:

  std::string name_;
  const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > linObjFactory_;
  const std::string lo_basis_name_;
  const std::string ho_basis_name_;
  Intrepid2::EOperator op_;
  const bool dump_;
  Teko::LinearOp interp_;
  const size_t worksetSize_;
  const bool matrixFree_;

public:

  InterpolationRequestCallback(const std::string& name,
                               const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > linObjFactory,
                               const std::string& lo_basis_name,
                               const std::string& ho_basis_name,
                               Intrepid2::EOperator op=Intrepid2::OPERATOR_VALUE,
                               const bool waitForRequest=true,
                               const bool dump=false,
                               const size_t worksetSize=1000,
                               const bool matrixFree=false)
  : name_(name), linObjFactory_(linObjFactory), lo_basis_name_(lo_basis_name), ho_basis_name_(ho_basis_name), op_(op), dump_(dump), worksetSize_(worksetSize), matrixFree_(matrixFree)
  {
    if (!waitForRequest)
      build();
  };

  void build()
  {
    if (!matrixFree_) {
      Teuchos::TimeMonitor tm(*Teuchos::TimeMonitor::getNewTimer(std::string("Mini-EM: assemble ") + name_));
      interp_ = buildInterpolation(linObjFactory_, lo_basis_name_, ho_basis_name_, op_, worksetSize_);
    } else {
      typedef double Scalar;
      typedef int LocalOrdinal;
      typedef panzer::GlobalOrdinal GlobalOrdinal;
      typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal> tp_matrix;

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

  bool handlesRequest(const Teko::RequestMesg & rm)
  {
    std::string name = rm.getName();

    return (name==name_);
  };

  Teko::LinearOp request(const Teko::RequestMesg & rm)
  {
    TEUCHOS_ASSERT(handlesRequest(rm));
    std::string name = rm.getName();

    if(name==name_) {
      if (interp_.is_null()) {
        build();
      }
      return interp_;
    } else
      TEUCHOS_ASSERT(false);
  };

  void preRequest(const Teko::RequestMesg & rm)
  {
    // checking for its existance is as good as pre requesting
    TEUCHOS_ASSERT(handlesRequest(rm));
  };
};




void addInterpolationToRequestHandler(
                                      const std::string& name,
                                      const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > linObjFactory,
                                      const Teuchos::RCP<Teko::RequestHandler> & reqHandler,
                                      const std::string& lo_basis_name,
                                      const std::string& ho_basis_name,
                                      Intrepid2::EOperator op=Intrepid2::OPERATOR_VALUE,
                                      const bool waitForRequest=true,
                                      const bool dump=false,
                                      const size_t worksetSize=1000,
                                      const bool matrixFree=false) {

  // add interpolation callback to request handler
  reqHandler->addRequestCallback(Teuchos::rcp(new InterpolationRequestCallback(name, linObjFactory, lo_basis_name, ho_basis_name, op, waitForRequest, dump, worksetSize, matrixFree)));
}

#endif
