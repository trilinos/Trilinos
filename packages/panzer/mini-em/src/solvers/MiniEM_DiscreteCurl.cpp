#include "MiniEM_DiscreteCurl.hpp"

void addDiscreteCurlToRequestHandler(
                                     const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > linObjFactory,
                                     const Teuchos::RCP<Teko::RequestHandler> & reqHandler)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;

  typedef double Scalar;
  typedef int LocalOrdinal;
  typedef panzer::GlobalOrdinal GlobalOrdinal;

  typedef typename panzer::BlockedTpetraLinearObjFactory<panzer::Traits,Scalar,LocalOrdinal,GlobalOrdinal> tpetraBlockedLinObjFactory;
#ifdef PANZER_HAVE_EPETRA_STACK
  typedef typename panzer::BlockedEpetraLinearObjFactory<panzer::Traits,LocalOrdinal> epetraBlockedLinObjFactory;
#endif
  typedef panzer::GlobalIndexer UGI;
  typedef PHX::Device DeviceSpace;
  typedef Kokkos::HostSpace HostSpace;
  typedef Intrepid2::OrientationTools<HostSpace> ots;
  typedef Intrepid2::Experimental::LagrangianInterpolation<HostSpace> li;

  // use "B_face" and "E_edge" as the DOFs
  std::string face_basis_name = "B_face";
  std::string edge_basis_name = "E_edge";

  // must be able to cast to a block linear object factory
  RCP<const tpetraBlockedLinObjFactory > tblof  = rcp_dynamic_cast<const tpetraBlockedLinObjFactory >(linObjFactory);
#ifdef PANZER_HAVE_EPETRA_STACK
  RCP<const epetraBlockedLinObjFactory > eblof  = rcp_dynamic_cast<const epetraBlockedLinObjFactory >(linObjFactory);
#endif
  if (tblof != Teuchos::null) {
    typedef typename panzer::BlockedTpetraLinearObjContainer<Scalar,LocalOrdinal,GlobalOrdinal> linObjContainer;
    typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,panzer::TpetraNodeType> matrix;
    typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,panzer::TpetraNodeType> map;

    RCP<const panzer::BlockedDOFManager> blockedDOFMngr = tblof->getGlobalIndexer();

    // get global indexers for faces and edges
    std::vector<RCP<UGI> > fieldDOFMngrs = blockedDOFMngr->getFieldDOFManagers();
    int fFieldNum = blockedDOFMngr->getFieldNum(face_basis_name);
    int eFieldNum = blockedDOFMngr->getFieldNum(edge_basis_name);
    int fBlockIndex = blockedDOFMngr->getFieldBlock(fFieldNum);
    int eBlockIndex = blockedDOFMngr->getFieldBlock(eFieldNum);
    RCP<panzer::DOFManager> face_ugi = rcp_dynamic_cast<panzer::DOFManager>(blockedDOFMngr->getFieldDOFManagers()[fBlockIndex],true);
    RCP<panzer::DOFManager> edge_ugi = rcp_dynamic_cast<panzer::DOFManager>(blockedDOFMngr->getFieldDOFManagers()[eBlockIndex],true);

    // extract ghosted and global linear object containers
    RCP<panzer::GlobalEvaluationData> dataObject
      = rcp(new panzer::LOCPair_GlobalEvaluationData(tblof,panzer::LinearObjContainer::Mat));
    RCP<panzer::LinearObjContainer> global_loc
      = rcp_dynamic_cast<panzer::LOCPair_GlobalEvaluationData>(dataObject,true)->getGlobalLOC();
    RCP<panzer::LinearObjContainer> ghosted_loc
      = rcp_dynamic_cast<panzer::LOCPair_GlobalEvaluationData>(dataObject,true)->getGhostedLOC();

    RCP<linObjContainer> global_tloc = rcp_dynamic_cast<linObjContainer>(global_loc,true);
    RCP<linObjContainer> ghosted_tloc = rcp_dynamic_cast<linObjContainer>(ghosted_loc,true);

    RCP<const map> rangemap  = global_tloc->getMapForBlock(fBlockIndex);
    RCP<const map> domainmap = global_tloc->getMapForBlock(eBlockIndex);
    RCP<const map> rowmap    = global_tloc->getMapForBlock(fBlockIndex);
    RCP<const map> colmap    = ghosted_tloc->getMapForBlock(eBlockIndex);

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

    // hexes require scaling
    double area_scaling = 1.0;
    if(numElemVertices == 8)
      area_scaling = 2.0;

    // get the HCurl and HDiv bases
    auto face_fieldPattern = face_ugi->getFieldPattern(face_basis_name);
    auto face_basis = rcp_dynamic_cast<const panzer::Intrepid2FieldPattern>(face_fieldPattern,true)->getIntrepidBasis();
    auto edge_fieldPattern = edge_ugi->getFieldPattern(edge_basis_name);
    auto edge_basis = rcp_dynamic_cast<const panzer::Intrepid2FieldPattern>(edge_fieldPattern,true)->getIntrepidBasis();

    // cardinalities
    int hdivCardinality = face_basis->getCardinality();
    int hcurlCardinality = edge_basis->getCardinality();

    // allocate some view
    Kokkos::DynRankView<double,DeviceSpace> dofCoords("dofCoords", 1, hdivCardinality, dim);
    Kokkos::DynRankView<double,HostSpace> basisCoeffsLI("basisCoeffsLI", 1, hcurlCardinality, hdivCardinality);
    typename Kokkos::DynRankView<Intrepid2::Orientation,HostSpace>::HostMirror elemOrts("elemOrts", 1);
    typename Kokkos::DynRankView<GlobalOrdinal, HostSpace>::HostMirror elemNodes("elemNodes", 1, numElemVertices);
    typename Kokkos::DynRankView<int, HostSpace>::HostMirror fOrt("fOrt", hdivCardinality);
    typename Kokkos::DynRankView<double, HostSpace>::HostMirror ortJacobian("ortJacobian", 2, 2);

    // the ranks of these depend on dimension
    Kokkos::DynRankView<double,HostSpace>   dofCoeffs;
    Kokkos::DynRankView<double,DeviceSpace> refDofCoeffs;
    Kokkos::DynRankView<double,DeviceSpace> curlAtDofCoordsNonOriented_d;
    Kokkos::DynRankView<double,HostSpace>   curlAtDofCoords;
    if(dim==3){
      dofCoeffs                    = Kokkos::DynRankView<double,HostSpace>("dofCoeffs", 1, hdivCardinality,dim);
      refDofCoeffs                 = Kokkos::DynRankView<double,DeviceSpace>("refDofCoeffs", hdivCardinality,dim);
      curlAtDofCoordsNonOriented_d = Kokkos::DynRankView<double,DeviceSpace>("curlAtDofCoordsNonOriented", 1, hcurlCardinality, hdivCardinality, dim);
      curlAtDofCoords              = Kokkos::DynRankView<double,HostSpace>("curlAtDofCoords", 1, hcurlCardinality, hdivCardinality, dim);
    } else {
      dofCoeffs                    = Kokkos::DynRankView<double,HostSpace>("dofCoeffs", 1, hdivCardinality);
      refDofCoeffs                 = Kokkos::DynRankView<double,DeviceSpace>("refDofCoeffs", hdivCardinality);
      curlAtDofCoordsNonOriented_d = Kokkos::DynRankView<double,DeviceSpace>("curlAtDofCoordsNonOriented", 1, hcurlCardinality, hdivCardinality);
      curlAtDofCoords              = Kokkos::DynRankView<double,HostSpace>("curlAtDofCoords", 1, hcurlCardinality, hdivCardinality);
    }
    auto curlAtDofCoordsNonOriented = Kokkos::create_mirror_view(curlAtDofCoordsNonOriented_d);
    face_basis->getDofCoeffs(refDofCoeffs);
    auto dofCoeffs_h = Kokkos::create_mirror_view(dofCoeffs);
    auto refDofCoeffs_h = Kokkos::create_mirror_view(refDofCoeffs);
    Kokkos::deep_copy(refDofCoeffs_h, refDofCoeffs);
    // set up the topology of each face in an element for computing DOF coefficients
    // in 2D coefficients are same as reference coefficients
    shards::CellTopology sub_topologies[hdivCardinality];
    if(dim < 3)
      for(int i = 0; i < hdivCardinality; i++)
        dofCoeffs_h(0,i) = refDofCoeffs_h(i);
    else {
      for(int iface = 0; iface < hdivCardinality; iface++){
        shards::CellTopology sub_topology(topology.getCellTopologyData(dim-1,iface));
        sub_topologies[iface] = sub_topology;
      }
    }

    // compute curls at dof coords
    edge_basis->getValues(Kokkos::subview(curlAtDofCoordsNonOriented_d, 0, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL()),
                          Kokkos::subview(dofCoords, 0, Kokkos::ALL(), Kokkos::ALL()), Intrepid2::OPERATOR_CURL);
    Kokkos::deep_copy(curlAtDofCoordsNonOriented, curlAtDofCoordsNonOriented_d);

    // create the global curl matrix
    RCP<matrix> curl_matrix = rcp(new matrix(rowmap,colmap,hcurlCardinality));

    // get IDs for edges and faces
    auto fLIDs_k = face_ugi->getLIDs();
    auto fLIDs = Kokkos::create_mirror_view(fLIDs_k);
    Kokkos::deep_copy(DeviceSpace::execution_space(), fLIDs, fLIDs_k);
    auto eLIDs_k = edge_ugi->getLIDs();
    auto eLIDs = Kokkos::create_mirror_view(eLIDs_k);
    Kokkos::deep_copy(DeviceSpace::execution_space(), eLIDs, eLIDs_k);

    std::vector<int> edges_on_face(hcurlCardinality,-1);

    // loop over element blocks
    std::vector<std::string> elementBlockIds;
    blockedDOFMngr->getElementBlockIds(elementBlockIds);
    for(std::size_t blockIter = 0; blockIter < elementBlockIds.size(); ++blockIter) {

      // loop over elements
      std::vector<int> elementIds = edge_ugi->getElementBlock(elementBlockIds[blockIter]);
      for(std::size_t elemIter = 0; elemIter < elementIds.size(); ++elemIter){

        int element = elementIds[elemIter];

        // get element orientations
        auto node_ids = node_conn->getConnectivity(element);
        for(int i = 0; i < numElemVertices; i++)
          elemNodes(0,i) = node_ids[i];

        elemOrts(0) = Intrepid2::Orientation::getOrientation(topology,Kokkos::subview(elemNodes,0,Kokkos::ALL()));

        //compute interpolation weights (dofCoeffs)
        if(dim==3){
          elemOrts(0).getFaceOrientation(fOrt.data(),hdivCardinality);
          for(int iface = 0; iface < hdivCardinality; iface++){
            Intrepid2::Impl::OrientationTools::getJacobianOfOrientationMap(ortJacobian, sub_topologies[iface], fOrt(iface));
            auto ortJacobianDet = ortJacobian(0,0)*ortJacobian(1,1)-ortJacobian(1,0)*ortJacobian(0,1);
            for(int idim = 0; idim < dim; idim++)
              dofCoeffs_h(0,iface,idim) = refDofCoeffs_h(iface,idim)*ortJacobianDet;
          }
        }

        //orient basis
        ots::modifyBasisByOrientation(curlAtDofCoords,
                                      curlAtDofCoordsNonOriented,
                                      elemOrts,
                                      edge_basis.get());

        //get basis coefficients (dofs)
        for(int curlIter=0; curlIter<hcurlCardinality; curlIter++)
          li::getBasisCoeffs(Kokkos::subview(basisCoeffsLI,Kokkos::ALL(),curlIter,Kokkos::ALL()),
                             Kokkos::subview(curlAtDofCoords,Kokkos::ALL(),curlIter,Kokkos::ALL(),Kokkos::ALL()),
                             dofCoeffs_h);

        // loop over faces
        for(std::size_t fIter = 0; fIter < fLIDs.extent(1); ++fIter){

          LocalOrdinal fLID = fLIDs(element,fIter);

          // get a list of edges on the face
          if(dim==3)
            edge_fieldPattern->getSubcellClosureIndices(2,fIter,edges_on_face);
          else
            for(int i = 0; i < hcurlCardinality; i++)
              edges_on_face[i] = i;

          const bool isOwned = rowmap->isNodeLocalElement(fLID);

          // if owned
          if(isOwned && curl_matrix->getNumEntriesInLocalRow(fLID) < 1){

            // get values from reference incidence matrix
            std::vector<Scalar> values(hcurlCardinality,0.0);
            for(int curlIter=0; curlIter<hcurlCardinality; curlIter++){
              // only add entries for edges on the face
              bool edge_is_on_face = false;
              for (int i = 0; i < hcurlCardinality; i++) {
                if (edges_on_face[i] == curlIter) {
                  edge_is_on_face = true;
                  break;
                }
              }
              if(edge_is_on_face) {
                // normalize the values
                if(std::abs(basisCoeffsLI(0,curlIter,fIter)) > 1.0e-10)
                  values[curlIter] = basisCoeffsLI(0,curlIter,fIter)*area_scaling;
              }
            }

            // insert values in matrix
            auto indices = Kokkos::subview(eLIDs, element, Kokkos::ALL());
            curl_matrix->insertLocalValues(fLID, values.size(), &values[0], &indices[0]);
            TEUCHOS_ASSERT_EQUALITY(curl_matrix->getNumEntriesInLocalRow(fLID),values.size());

          }//end if owned
        }//end face loop
      }//end element loop
    }//end element block loop
    curl_matrix->fillComplete(domainmap,rangemap);

    RCP<Thyra::LinearOpBase<Scalar> > thyra_curl = Thyra::tpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,typename matrix::node_type>(Thyra::createVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal>(rangemap),
                                                                                                                                       Thyra::createVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal>(domainmap),
                                                                                                                                       curl_matrix);

    // add curl callback to request handler
    reqHandler->addRequestCallback(Teuchos::rcp(new CurlRequestCallback(thyra_curl)));

#ifdef PANZER_HAVE_EPETRA_STACK
  }  else if (eblof != Teuchos::null) {

    typedef typename panzer::BlockedEpetraLinearObjContainer linObjContainer;
    typedef Epetra_CrsMatrix matrix;
    typedef Epetra_Map map;

    RCP<const panzer::BlockedDOFManager> blockedDOFMngr = eblof->getGlobalIndexer();

    // get global indexers for faces and edges
    std::vector<RCP<UGI> > fieldDOFMngrs = blockedDOFMngr->getFieldDOFManagers();
    int fFieldNum = blockedDOFMngr->getFieldNum(face_basis_name);
    int eFieldNum = blockedDOFMngr->getFieldNum(edge_basis_name);
    int fBlockIndex = blockedDOFMngr->getFieldBlock(fFieldNum);
    int eBlockIndex = blockedDOFMngr->getFieldBlock(eFieldNum);
    RCP<panzer::DOFManager> face_ugi = rcp_dynamic_cast<panzer::DOFManager>(blockedDOFMngr->getFieldDOFManagers()[fBlockIndex],true);
    RCP<panzer::DOFManager> edge_ugi = rcp_dynamic_cast<panzer::DOFManager>(blockedDOFMngr->getFieldDOFManagers()[eBlockIndex],true);

    // extract ghosted and global linear object containers
    RCP<panzer::GlobalEvaluationData> dataObject
      = rcp(new panzer::LOCPair_GlobalEvaluationData(eblof,panzer::LinearObjContainer::Mat));
    RCP<panzer::LinearObjContainer> global_loc
      = rcp_dynamic_cast<panzer::LOCPair_GlobalEvaluationData>(dataObject,true)->getGlobalLOC();
    RCP<panzer::LinearObjContainer> ghosted_loc
      = rcp_dynamic_cast<panzer::LOCPair_GlobalEvaluationData>(dataObject,true)->getGhostedLOC();

    RCP<linObjContainer> global_tloc = rcp_dynamic_cast<linObjContainer>(global_loc,true);
    RCP<linObjContainer> ghosted_tloc = rcp_dynamic_cast<linObjContainer>(ghosted_loc,true);

    RCP<const map> rangemap  = global_tloc->getMapForBlock(fBlockIndex);
    RCP<const map> domainmap = global_tloc->getMapForBlock(eBlockIndex);
    RCP<const map> rowmap    = global_tloc->getMapForBlock(fBlockIndex);
    RCP<const map> colmap    = ghosted_tloc->getMapForBlock(eBlockIndex);

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

    // hexes require scaling
    double area_scaling = 1.0;
    if(numElemVertices == 8)
      area_scaling = 2.0;

    // get the HCurl and HDiv bases
    auto face_fieldPattern = face_ugi->getFieldPattern(face_basis_name);
    auto face_basis = rcp_dynamic_cast<const panzer::Intrepid2FieldPattern>(face_fieldPattern,true)->getIntrepidBasis()->getHostBasis();
    auto edge_fieldPattern = edge_ugi->getFieldPattern(edge_basis_name);
    auto edge_basis = rcp_dynamic_cast<const panzer::Intrepid2FieldPattern>(edge_fieldPattern,true)->getIntrepidBasis()->getHostBasis();

    // cardinalities
    int hdivCardinality = face_basis->getCardinality();
    int hcurlCardinality = edge_basis->getCardinality();

    // allocate some view
    Kokkos::DynRankView<double,HostSpace> dofCoords("dofCoords", 1, hdivCardinality, dim);
    Kokkos::DynRankView<double,HostSpace> basisCoeffsLI("basisCoeffsLI", 1, hcurlCardinality, hdivCardinality);
    typename Kokkos::DynRankView<Intrepid2::Orientation,DeviceSpace>::HostMirror elemOrts("elemOrts", 1);
    typename Kokkos::DynRankView<GlobalOrdinal, DeviceSpace>::HostMirror elemNodes("elemNodes", 1, numElemVertices);
    typename Kokkos::DynRankView<int, DeviceSpace>::HostMirror fOrt("fOrt", hdivCardinality);
    typename Kokkos::DynRankView<double, DeviceSpace>::HostMirror ortJacobian("ortJacobian", 2, 2);

    // the ranks of these depend on dimension
    Kokkos::DynRankView<double,HostSpace> dofCoeffs;
    Kokkos::DynRankView<double,HostSpace> refDofCoeffs;
    Kokkos::DynRankView<double,HostSpace> curlAtDofCoordsNonOriented;
    Kokkos::DynRankView<double,HostSpace> curlAtDofCoords;
    if(dim==3){
      dofCoeffs                  = Kokkos::DynRankView<double,HostSpace>("dofCoeffs", 1, hdivCardinality,dim);
      refDofCoeffs               = Kokkos::DynRankView<double,HostSpace>("refDofCoeffs", hdivCardinality,dim);
      curlAtDofCoordsNonOriented = Kokkos::DynRankView<double,HostSpace>("curlAtDofCoordsNonOriented", 1, hcurlCardinality, hdivCardinality, dim);
      curlAtDofCoords            = Kokkos::DynRankView<double,HostSpace>("curlAtDofCoords", 1, hcurlCardinality, hdivCardinality, dim);
    } else {
      dofCoeffs                  = Kokkos::DynRankView<double,HostSpace>("dofCoeffs", 1, hdivCardinality);
      refDofCoeffs               = Kokkos::DynRankView<double,HostSpace>("refDofCoeffs", hdivCardinality);
      curlAtDofCoordsNonOriented = Kokkos::DynRankView<double,HostSpace>("curlAtDofCoordsNonOriented", 1, hcurlCardinality, hdivCardinality);
      curlAtDofCoords            = Kokkos::DynRankView<double,HostSpace>("curlAtDofCoords", 1, hcurlCardinality, hdivCardinality);
    }
    face_basis->getDofCoeffs(refDofCoeffs);

    auto dofCoeffs_h = Kokkos::create_mirror_view(dofCoeffs);
    auto refDofCoeffs_h = Kokkos::create_mirror_view(refDofCoeffs);
    Kokkos::deep_copy(refDofCoeffs_h, refDofCoeffs);
    // set up the topology of each face in an element for computing DOF coefficients
    // in 2D coefficients are same as reference coefficients
    shards::CellTopology sub_topologies[hdivCardinality];
    if(dim < 3)
      for(int i = 0; i < hdivCardinality; i++)
        dofCoeffs_h(0,i) = refDofCoeffs_h(i);
    else {
      for(int iface = 0; iface < hdivCardinality; iface++){
        shards::CellTopology sub_topology(topology.getCellTopologyData(dim-1,iface));
        sub_topologies[iface] = sub_topology;
      }
    }

    // compute curls at dof coords
    edge_basis->getValues(Kokkos::subview(curlAtDofCoordsNonOriented, 0, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL()), Kokkos::subview(dofCoords, 0, Kokkos::ALL(), Kokkos::ALL()), Intrepid2::OPERATOR_CURL);

    // create the global curl matrix
    RCP<matrix> curl_matrix = rcp(new matrix(Copy, *rowmap, *colmap, basisCoeffsLI.extent(1), /*StaticProfile=*/true));

    // loop over element blocks
    std::vector<std::string> elementBlockIds;
    blockedDOFMngr->getElementBlockIds(elementBlockIds);
    for(std::size_t blockIter = 0; blockIter < elementBlockIds.size(); ++blockIter) {

      // loop over elements
      std::vector<int> elementIds = edge_ugi->getElementBlock(elementBlockIds[blockIter]);
      for(std::size_t elemIter = 0; elemIter < elementIds.size(); ++elemIter){

        // get element orientations

        auto node_ids = node_conn->getConnectivity(elementIds[elemIter]);

        for(int i = 0; i < numElemVertices; i++)
          elemNodes(0,i) = node_ids[i];

        elemOrts(0) = Intrepid2::Orientation::getOrientation(topology,Kokkos::subview(elemNodes,0,Kokkos::ALL()));

        //compute interpolation weights (dofCoeffs)
        if(dim==3){
          elemOrts(0).getFaceOrientation(fOrt.data(),hdivCardinality);
          for(int iface = 0; iface < hdivCardinality; iface++){
            Intrepid2::Impl::OrientationTools::getJacobianOfOrientationMap(ortJacobian, sub_topologies[iface], fOrt(iface));
            auto ortJacobianDet = ortJacobian(0,0)*ortJacobian(1,1)-ortJacobian(1,0)*ortJacobian(0,1);
            for(int idim = 0; idim < dim; idim++)
              dofCoeffs_h(0,iface,idim) = refDofCoeffs_h(iface,idim)*ortJacobianDet;
          }
        }
        Kokkos::deep_copy(dofCoeffs, dofCoeffs_h);

        //orient basis
        // Kokkos::DynRankView<Intrepid2::Orientation,DeviceSpace> elemOrts_d("elemOrts_d", 1);
        // Kokkos::deep_copy(elemOrts_d, elemOrts);
        ots::modifyBasisByOrientation(curlAtDofCoords,
                                      curlAtDofCoordsNonOriented,
                                      elemOrts,
                                      edge_basis.get());

        //get basis coefficients (dofs)
        for(int curlIter=0; curlIter<hcurlCardinality; curlIter++)
          li::getBasisCoeffs(Kokkos::subview(basisCoeffsLI,Kokkos::ALL(),curlIter,Kokkos::ALL()), Kokkos::subview(curlAtDofCoords,Kokkos::ALL(),curlIter,Kokkos::ALL(),Kokkos::ALL()), dofCoeffs);
        auto basisCoeffsLI_h = Kokkos::create_mirror_view(basisCoeffsLI);
        Kokkos::deep_copy(basisCoeffsLI_h, basisCoeffsLI);

        // get IDs for edges and faces
        std::vector<GlobalOrdinal> fGIDs;
        face_ugi->getElementGIDs(elementIds[elemIter],fGIDs);
        std::vector<GlobalOrdinal> eGIDs;
        edge_ugi->getElementGIDs(elementIds[elemIter],eGIDs);
        auto eLIDs_k = edge_ugi->getElementLIDs(elementIds[elemIter]);
        auto fLIDs_k = face_ugi->getElementLIDs(elementIds[elemIter]);
        auto eLIDs = Kokkos::create_mirror_view(eLIDs_k);
        auto fLIDs = Kokkos::create_mirror_view(fLIDs_k);
        Kokkos::deep_copy(eLIDs, eLIDs_k);
        Kokkos::deep_copy(fLIDs, fLIDs_k);

        // need to know which faces are owned by this proc
        std::vector<bool> isOwned(fGIDs.size());
        for (size_t face=0; face < fGIDs.size(); ++face)
          isOwned[face] = (rowmap->MyGID(fGIDs[face]));

        // loop over faces
        for(std::size_t fIter = 0; fIter < fLIDs.size(); ++fIter){

          // get a list of edges on the face
          std::vector<int> edges_on_face(hcurlCardinality,-1);
          if(dim==3)
            edge_fieldPattern->getSubcellClosureIndices(2,fIter,edges_on_face);
          else
            for(int i = 0; i < hcurlCardinality; i++)
              edges_on_face[i] = i;

          // if owned
          if(isOwned[fIter] && curl_matrix->NumMyEntries(fLIDs[fIter]) < 1){

            // get values from reference incidence matrix
            std::vector<Scalar> values(hcurlCardinality,0.0);
            for(int curlIter=0; curlIter<hcurlCardinality; curlIter++){
              // only add entries for edges on the face
              auto it = std::find(edges_on_face.begin(), edges_on_face.end(), curlIter);
              if(it!=edges_on_face.end()){
                // normalize the values
                if(std::abs(basisCoeffsLI_h(0,curlIter,fIter)) > 1.0e-10)
                  values[curlIter] = basisCoeffsLI_h(0,curlIter,fIter)*area_scaling;
              }
            }

            // insert values in matrix
            int err = curl_matrix->InsertMyValues(fLIDs[fIter], values.size(), &values[0], &eLIDs[0]);
            TEUCHOS_ASSERT_EQUALITY(err,0);

          }//end if owned
        }//end face loop
      }//end element loop
    }//end element block loop
    curl_matrix->FillComplete(*domainmap, *rangemap);

    RCP<const Thyra::LinearOpBase<Scalar> > thyra_curl = Thyra::epetraLinearOp(curl_matrix,
                                                                               Thyra::NOTRANS,
                                                                               Thyra::EPETRA_OP_APPLY_APPLY,
                                                                               Thyra::EPETRA_OP_ADJOINT_SUPPORTED,
                                                                               Thyra::create_VectorSpace(rangemap),
                                                                               Thyra::create_VectorSpace(domainmap));

    // add curl callback to request handler
    reqHandler->addRequestCallback(Teuchos::rcp(new CurlRequestCallback(thyra_curl)));
#endif
  } else
    TEUCHOS_ASSERT(false);
}
