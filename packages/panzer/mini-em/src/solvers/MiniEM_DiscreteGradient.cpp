#include "MiniEM_DiscreteGradient.hpp"

void addDiscreteGradientToRequestHandler(
                                         const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > linObjFactory,
                                         const Teuchos::RCP<Teko::RequestHandler> & reqHandler)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;

  typedef double Scalar;
#ifdef PANZER_HAVE_EPETRA_STACK
  typedef int LocalOrdinalEpetra;
#endif
  typedef int LocalOrdinalTpetra;
  typedef panzer::GlobalOrdinal GlobalOrdinalTpetra;

  typedef typename panzer::BlockedTpetraLinearObjFactory<panzer::Traits,Scalar,LocalOrdinalTpetra,GlobalOrdinalTpetra> tpetraBlockedLinObjFactory;
#ifdef PANZER_HAVE_EPETRA_STACK
  typedef typename panzer::BlockedEpetraLinearObjFactory<panzer::Traits,LocalOrdinalEpetra> epetraBlockedLinObjFactory;
#endif

  typedef PHX::Device DeviceSpace;

  // use "AUXILIARY_EDGE" and "AUXILIARY_NODE" as the DOFs
  std::string edge_basis_name  = "AUXILIARY_EDGE";
  std::string nodal_basis_name = "AUXILIARY_NODE";

  // must be able to cast to a block linear object factory
  RCP<const tpetraBlockedLinObjFactory > tblof  = rcp_dynamic_cast<const tpetraBlockedLinObjFactory >(linObjFactory);
#ifdef PANZER_HAVE_EPETRA_STACK
  RCP<const epetraBlockedLinObjFactory > eblof  = rcp_dynamic_cast<const epetraBlockedLinObjFactory >(linObjFactory);
#endif
  if (tblof != Teuchos::null) {
    typedef LocalOrdinalTpetra LocalOrdinal;
    typedef GlobalOrdinalTpetra GlobalOrdinal;
    typedef panzer::GlobalIndexer UGI;
    typedef typename panzer::BlockedTpetraLinearObjContainer<Scalar,LocalOrdinal,GlobalOrdinal> linObjContainer;
    typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,panzer::TpetraNodeType> matrix;
    typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal,panzer::TpetraNodeType> map;

    RCP<const panzer::BlockedDOFManager> blockedDOFMngr = tblof->getGlobalIndexer();

    // get global indexers for edges and nodes
    std::vector<RCP<UGI> > fieldDOFMngrs = blockedDOFMngr->getFieldDOFManagers();
    int eFieldNum = blockedDOFMngr->getFieldNum(edge_basis_name);
    int nFieldNum = blockedDOFMngr->getFieldNum(nodal_basis_name);
    int eBlockIndex = blockedDOFMngr->getFieldBlock(eFieldNum);
    int nBlockIndex = blockedDOFMngr->getFieldBlock(nFieldNum);
    RCP<UGI> eUgi = blockedDOFMngr->getFieldDOFManagers()[eBlockIndex];
    RCP<UGI> nUgi = blockedDOFMngr->getFieldDOFManagers()[nBlockIndex];

    // extract ghosted and global linear object containers
    RCP<panzer::GlobalEvaluationData> dataObject
      = rcp(new panzer::LOCPair_GlobalEvaluationData(tblof,panzer::LinearObjContainer::Mat));
    RCP<panzer::LinearObjContainer> global_loc
      = rcp_dynamic_cast<panzer::LOCPair_GlobalEvaluationData>(dataObject,true)->getGlobalLOC();
    RCP<panzer::LinearObjContainer> ghosted_loc
      = rcp_dynamic_cast<panzer::LOCPair_GlobalEvaluationData>(dataObject,true)->getGhostedLOC();

    RCP<linObjContainer> global_tloc = rcp_dynamic_cast<linObjContainer>(global_loc,true);
    RCP<linObjContainer> ghosted_tloc = rcp_dynamic_cast<linObjContainer>(ghosted_loc,true);

    RCP<const map> rangemap  = global_tloc->getMapForBlock(eBlockIndex);
    RCP<const map> domainmap = global_tloc->getMapForBlock(nBlockIndex);
    RCP<const map> rowmap    = global_tloc->getMapForBlock(eBlockIndex);
    RCP<const map> colmap    = ghosted_tloc->getMapForBlock(nBlockIndex);
    RCP<matrix> grad_matrix = rcp(new matrix(rowmap, colmap, 2));

    RCP<const panzer::FieldPattern> field_pattern = blockedDOFMngr->getGeometricFieldPattern();
    shards::CellTopology cell_topology = field_pattern->getCellTopology();
    std::vector<Intrepid2::Orientation> orientations = *panzer::buildIntrepidOrientation(blockedDOFMngr);
    // loop over element blocks
    std::vector<std::string> elementBlockIds;
    blockedDOFMngr->getElementBlockIds(elementBlockIds);
    std::vector<bool> insertedEdges(rowmap->getLocalNumElements(),false);
    auto eLIDs_k = eUgi->getLIDs();
    auto nLIDs_k = nUgi->getLIDs();
    auto eLIDs = Kokkos::create_mirror_view(eLIDs_k);
    auto nLIDs = Kokkos::create_mirror_view(nLIDs_k);
    Kokkos::deep_copy(DeviceSpace::execution_space(), eLIDs, eLIDs_k);
    Kokkos::deep_copy(DeviceSpace::execution_space(), nLIDs, nLIDs_k);

    for(std::size_t blockIter = 0; blockIter < elementBlockIds.size(); ++blockIter) {

      // loop over elements
      std::vector<int> elementIds = blockedDOFMngr->getElementBlock(elementBlockIds[blockIter]);
      for(std::size_t elemIter = 0; elemIter < elementIds.size(); ++elemIter){

        int element = elementIds[elemIter];

        // get element orientations
        std::vector<int> eFieldOffsets = blockedDOFMngr->getGIDFieldOffsets(elementBlockIds[blockIter],eFieldNum);
        std::vector<int> ort(eFieldOffsets.size(),0);
        orientations[element].getEdgeOrientation(&ort[0], eFieldOffsets.size());

        // loop over edges
        for(std::size_t eIter = 0; eIter < eFieldOffsets.size(); ++eIter){

          const bool isOwned = rowmap->isNodeLocalElement(eLIDs(element,eIter));

          if(isOwned && !insertedEdges[eLIDs(element, eIter)]){

            int headIndex = cell_topology.getNodeMap(1, eIter, 0);
            int tailIndex = cell_topology.getNodeMap(1, eIter, 1);

            // assign values based on orientation of edge (-1 for tail, 1 for head)
            double values[2] = {1.0,-1.0};
            if(ort[eIter] == 0)
              {values[0] *= -1.0; values[1] *= -1.0;}

            // get LIDs associated with nodes
            int indices[2] = {nLIDs(element,tailIndex),nLIDs(element,headIndex)};
            grad_matrix->insertLocalValues(eLIDs(element,eIter), 2, values, indices);
            insertedEdges[eLIDs(element,eIter)] = true;
          }//end if
        }//end edge loop
      }//end element loop
    }//end element block loop
    grad_matrix->fillComplete(domainmap,rangemap);

    RCP<Thyra::LinearOpBase<double> > thyra_gradient = Thyra::tpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,typename matrix::node_type>(Thyra::createVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal>(rangemap),
                                                                                                                                           Thyra::createVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal>(domainmap),
                                                                                                                                           grad_matrix);

    // add gradient callback to request handler
    reqHandler->addRequestCallback(Teuchos::rcp(new GradientRequestCallback(thyra_gradient)));
#ifdef PANZER_HAVE_EPETRA_STACK
  } else if (eblof != Teuchos::null) {
    typedef panzer::GlobalIndexer UGI;
    typedef typename panzer::BlockedEpetraLinearObjContainer linObjContainer;
    typedef Epetra_CrsMatrix matrix;
    typedef Epetra_Map map;

    RCP<const panzer::BlockedDOFManager> blockedDOFMngr = eblof->getGlobalIndexer();

    // get global indexers for edges and nodes
    std::vector<RCP<UGI> > fieldDOFMngrs = blockedDOFMngr->getFieldDOFManagers();
    int eFieldNum = blockedDOFMngr->getFieldNum(edge_basis_name);
    int nFieldNum = blockedDOFMngr->getFieldNum(nodal_basis_name);
    int eBlockIndex = blockedDOFMngr->getFieldBlock(eFieldNum);
    int nBlockIndex = blockedDOFMngr->getFieldBlock(nFieldNum);
    RCP<UGI> eUgi = blockedDOFMngr->getFieldDOFManagers()[eBlockIndex];
    RCP<UGI> nUgi = blockedDOFMngr->getFieldDOFManagers()[nBlockIndex];

    // extract ghosted and global linear object containers
    RCP<panzer::GlobalEvaluationData> dataObject
      = rcp(new panzer::LOCPair_GlobalEvaluationData(eblof,panzer::LinearObjContainer::Mat));
    RCP<panzer::LinearObjContainer> global_loc
      = rcp_dynamic_cast<panzer::LOCPair_GlobalEvaluationData>(dataObject,true)->getGlobalLOC();
    RCP<panzer::LinearObjContainer> ghosted_loc
      = rcp_dynamic_cast<panzer::LOCPair_GlobalEvaluationData>(dataObject,true)->getGhostedLOC();

    RCP<linObjContainer> global_tloc = rcp_dynamic_cast<linObjContainer>(global_loc,true);
    RCP<linObjContainer> ghosted_tloc = rcp_dynamic_cast<linObjContainer>(ghosted_loc,true);

    RCP<const map> rangemap  = global_tloc->getMapForBlock(eBlockIndex);
    RCP<const map> domainmap = global_tloc->getMapForBlock(nBlockIndex);
    RCP<const map> rowmap    = global_tloc->getMapForBlock(eBlockIndex);
    RCP<const map> colmap    = ghosted_tloc->getMapForBlock(nBlockIndex);
    RCP<matrix> grad_matrix = rcp(new matrix(Copy, *rowmap, *colmap, 2, /*StaticProfile=*/true));

    RCP<const panzer::FieldPattern> field_pattern = blockedDOFMngr->getGeometricFieldPattern();
    shards::CellTopology cell_topology = field_pattern->getCellTopology();
    std::vector<Intrepid2::Orientation> orientations = *panzer::buildIntrepidOrientation(blockedDOFMngr);
    // loop over element blocks
    std::vector<std::string> elementBlockIds;
    blockedDOFMngr->getElementBlockIds(elementBlockIds);
    std::vector<bool> insertedEdges(rowmap->NumMyElements(),false);
    for(std::size_t blockIter = 0; blockIter < elementBlockIds.size(); ++blockIter) {

      // loop over elements
      std::vector<int> elementIds = blockedDOFMngr->getElementBlock(elementBlockIds[blockIter]);
      for(std::size_t elemIter = 0; elemIter < elementIds.size(); ++elemIter){

        // get IDs for edges and nodes
        std::vector<panzer::GlobalOrdinal> eGIDs;
        eUgi->getElementGIDs(elementIds[elemIter],eGIDs);
        std::vector<panzer::GlobalOrdinal> nGIDs;
        nUgi->getElementGIDs(elementIds[elemIter],nGIDs);
        auto eLIDs_k = eUgi->getElementLIDs(elementIds[elemIter]);
        auto nLIDs_k = nUgi->getElementLIDs(elementIds[elemIter]);
        auto eLIDs = Kokkos::create_mirror_view(eLIDs_k);
        auto nLIDs = Kokkos::create_mirror_view(nLIDs_k);
        Kokkos::deep_copy(eLIDs, eLIDs_k);
        Kokkos::deep_copy(nLIDs, nLIDs_k);

        std::vector<bool> isOwned;
        eUgi->ownedIndices(eGIDs,isOwned);

        // get element orientations
        std::vector<int> eFieldOffsets = blockedDOFMngr->getGIDFieldOffsets(elementBlockIds[blockIter],eFieldNum);
        std::vector<int> ort(eFieldOffsets.size(),0);
        orientations[elementIds[elemIter]].getEdgeOrientation(&ort[0], eFieldOffsets.size());

        // loop over edges
        for(std::size_t eIter = 0; eIter < eFieldOffsets.size(); ++eIter){

          if(isOwned[eIter] && !insertedEdges[eLIDs[eIter]]){

            int headIndex = cell_topology.getNodeMap(1, eIter, 0);
            int tailIndex = cell_topology.getNodeMap(1, eIter, 1);

            // assign values based on orientation of edge (-1 for tail, 1 for head)
            double values[2] = {1.0,-1.0};
            if(ort[eIter] == 0)
              {values[0] *= -1.0; values[1] *= -1.0;}

            // get LIDs associated with nodes
            int indices[2] = {nLIDs[tailIndex],nLIDs[headIndex]};

            // insert values in matrix
            int err = grad_matrix->InsertMyValues(eLIDs[eIter], 2, values, indices);
            TEUCHOS_ASSERT_EQUALITY(err,0);
            insertedEdges[eLIDs[eIter]] = true;
          }//end if
        }//end edge loop
      }//end element loop
    }//end element block loop
    grad_matrix->FillComplete(*domainmap, *rangemap);

    RCP<const Thyra::LinearOpBase<double> > thyra_gradient = Thyra::epetraLinearOp(grad_matrix,
                                                                                   Thyra::NOTRANS,
                                                                                   Thyra::EPETRA_OP_APPLY_APPLY,
                                                                                   Thyra::EPETRA_OP_ADJOINT_SUPPORTED,
                                                                                   Thyra::create_VectorSpace(rangemap),
                                                                                   Thyra::create_VectorSpace(domainmap));

    // add gradient callback to request handler
    reqHandler->addRequestCallback(Teuchos::rcp(new GradientRequestCallback(thyra_gradient)));
#endif
  } else
    TEUCHOS_ASSERT(false);


}
