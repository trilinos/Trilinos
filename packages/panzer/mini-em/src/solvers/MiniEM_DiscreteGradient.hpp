#include "Panzer_LOCPair_GlobalEvaluationData.hpp"
#include "Panzer_IntrepidOrientation.hpp"
#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Tpetra_Import.hpp"

class GradientRequestCallback : public Teko::RequestCallback<Teko::LinearOp> {
private:

  Teko::LinearOp gradient_;

public:

  GradientRequestCallback(const Teko::LinearOp & gradient)
     : gradient_(gradient) {};

  bool handlesRequest(const Teko::RequestMesg & rm)
  {
    std::string name = rm.getName();

    return (name=="Discrete Gradient");
  };

  Teko::LinearOp request(const Teko::RequestMesg & rm)
  {
    TEUCHOS_ASSERT(handlesRequest(rm));
    std::string name = rm.getName();

    if(name=="Discrete Gradient")
      return gradient_;
    else
      TEUCHOS_ASSERT(false);
  };

  void preRequest(const Teko::RequestMesg & rm)
  {
    // checking for its existance is as good as pre requesting
    TEUCHOS_ASSERT(handlesRequest(rm));
  };
};

void addDiscreteGradientToRequestHandler(
       const Teuchos::RCP<const panzer::LinearObjFactory<panzer::Traits> > linObjFactory,
       const Teuchos::RCP<Teko::RequestHandler> & reqHandler)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;

  typedef double Scalar;
  typedef int LocalOrdinalEpetra;
  typedef int GlobalOrdinalEpetra;
  typedef int LocalOrdinalTpetra;
  typedef panzer::Ordinal64 GlobalOrdinalTpetra;

  typedef typename panzer::BlockedTpetraLinearObjFactory<panzer::Traits,Scalar,LocalOrdinalTpetra,GlobalOrdinalTpetra> tpetraBlockedLinObjFactory;
  typedef typename panzer::BlockedEpetraLinearObjFactory<panzer::Traits,LocalOrdinalEpetra> epetraBlockedLinObjFactory;

  // use "AUXILIARY_EDGE" and "AUXILIARY_NODE" as the DOFs
  std::string edge_basis_name  = "AUXILIARY_EDGE";
  std::string nodal_basis_name = "AUXILIARY_NODE";
      
  // must be able to cast to a block linear object factory
  RCP<const tpetraBlockedLinObjFactory > tblof  = rcp_dynamic_cast<const tpetraBlockedLinObjFactory >(linObjFactory);
  RCP<const epetraBlockedLinObjFactory > eblof  = rcp_dynamic_cast<const epetraBlockedLinObjFactory >(linObjFactory);
  if (tblof != Teuchos::null) {
    typedef LocalOrdinalTpetra LocalOrdinal;
    typedef GlobalOrdinalTpetra GlobalOrdinal;
    typedef panzer::UniqueGlobalIndexer<LocalOrdinal,GlobalOrdinal> UGI;
    typedef typename panzer::BlockedTpetraLinearObjContainer<Scalar,LocalOrdinal,GlobalOrdinal> linObjContainer;
    typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal> matrix;
    typedef Tpetra::Map<LocalOrdinal,GlobalOrdinal> map;

    RCP<const panzer::BlockedDOFManager<LocalOrdinal,GlobalOrdinal> > blockedDOFMngr = tblof->getGlobalIndexer();

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
    RCP<matrix> grad_matrix = rcp(new matrix(rowmap, colmap, 2, Tpetra::StaticProfile));

    RCP<const panzer::FieldPattern> field_pattern = blockedDOFMngr->getGeometricFieldPattern();
    shards::CellTopology cell_topology = field_pattern->getCellTopology();
    std::vector<Intrepid2::Orientation> orientations = *panzer::buildIntrepidOrientation(blockedDOFMngr);
    // loop over element blocks
    std::vector<std::string> elementBlockIds;
    blockedDOFMngr->getElementBlockIds(elementBlockIds);
    for(std::size_t blockIter = 0; blockIter < elementBlockIds.size(); ++blockIter) {

      // loop over elements
      std::vector<int> elementIds = blockedDOFMngr->getElementBlock(elementBlockIds[blockIter]);
      for(std::size_t elemIter = 0; elemIter < elementIds.size(); ++elemIter){
 
        // get IDs for edges and nodes
        std::vector<GlobalOrdinal> eGIDs;
        eUgi->getElementGIDs(elementIds[elemIter],eGIDs);
        std::vector<GlobalOrdinal> nGIDs;
        nUgi->getElementGIDs(elementIds[elemIter],nGIDs);
        auto eLIDs = eUgi->getElementLIDs(elementIds[elemIter]);
        auto nLIDs = nUgi->getElementLIDs(elementIds[elemIter]);
 
        std::vector<bool> isOwned;
        eUgi->ownedIndices(eGIDs,isOwned);

        // get element orientations
        std::vector<int> eFieldOffsets = blockedDOFMngr->getGIDFieldOffsets(elementBlockIds[blockIter],eFieldNum);
        std::vector<int> ort(eFieldOffsets.size(),0);
        orientations[elementIds[elemIter]].getEdgeOrientation(&ort[0], eFieldOffsets.size());

        // loop over edges
        for(std::size_t eIter = 0; eIter < eFieldOffsets.size(); ++eIter){
 
          if(isOwned[eIter]){

            int headIndex = cell_topology.getNodeMap(1, eIter, 0);
            int tailIndex = cell_topology.getNodeMap(1, eIter, 1);

            // assign values based on orientation of edge (-1 for tail, 1 for head)
            double values[2] = {1.0,-1.0};
            if(ort[eIter] == 0)
              {values[0] *= -1.0; values[1] *= -1.0;}
 
            // get LIDs associated with nodes
            int indices[2] = {nLIDs[tailIndex],nLIDs[headIndex]};
 
            // insert values in matrix
            LocalOrdinal err = grad_matrix->replaceLocalValues(eLIDs[eIter], 2, values, indices);
            if (err<2)
              grad_matrix->insertLocalValues(eLIDs[eIter], 2, values, indices);
          }//end if
        }//end edge loop
      }//end element loop
    }//end element block loop
    grad_matrix->fillComplete(domainmap,rangemap);

    RCP<Thyra::LinearOpBase<double> > thyra_gradient = Thyra::tpetraLinearOp<Scalar,LocalOrdinal,GlobalOrdinal,typename matrix::node_type>(Thyra::createVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal>(domainmap),
                                                                                                                                           Thyra::createVectorSpace<Scalar,LocalOrdinal,GlobalOrdinal>(rangemap),grad_matrix);

    // add gradient callback to request handler
    reqHandler->addRequestCallback(Teuchos::rcp(new GradientRequestCallback(thyra_gradient)));
  } else if (eblof != Teuchos::null) {
    typedef LocalOrdinalEpetra LocalOrdinal;
    typedef GlobalOrdinalEpetra GlobalOrdinal;
    typedef panzer::UniqueGlobalIndexer<LocalOrdinal,GlobalOrdinal> UGI;
    typedef typename panzer::BlockedEpetraLinearObjContainer linObjContainer;
    typedef Epetra_CrsMatrix matrix;
    typedef Epetra_Map map;

    RCP<const panzer::BlockedDOFManager<LocalOrdinal,GlobalOrdinal> > blockedDOFMngr = eblof->getGlobalIndexer();

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
    for(std::size_t blockIter = 0; blockIter < elementBlockIds.size(); ++blockIter) {

      // loop over elements
      std::vector<int> elementIds = blockedDOFMngr->getElementBlock(elementBlockIds[blockIter]);
      for(std::size_t elemIter = 0; elemIter < elementIds.size(); ++elemIter){

        // get IDs for edges and nodes
        std::vector<GlobalOrdinal> eGIDs;
        eUgi->getElementGIDs(elementIds[elemIter],eGIDs);
        std::vector<GlobalOrdinal> nGIDs;
        nUgi->getElementGIDs(elementIds[elemIter],nGIDs);
        auto eLIDs = eUgi->getElementLIDs(elementIds[elemIter]);
        auto nLIDs = nUgi->getElementLIDs(elementIds[elemIter]);

        std::vector<bool> isOwned;
        eUgi->ownedIndices(eGIDs,isOwned);

        // get element orientations
        std::vector<int> eFieldOffsets = blockedDOFMngr->getGIDFieldOffsets(elementBlockIds[blockIter],eFieldNum);
        std::vector<int> ort(eFieldOffsets.size(),0);
        orientations[elementIds[elemIter]].getEdgeOrientation(&ort[0], eFieldOffsets.size());

        // loop over edges
        for(std::size_t eIter = 0; eIter < eFieldOffsets.size(); ++eIter){

          if(isOwned[eIter]){

            int headIndex = cell_topology.getNodeMap(1, eIter, 0);
            int tailIndex = cell_topology.getNodeMap(1, eIter, 1);

            // assign values based on orientation of edge (-1 for tail, 1 for head)
            double values[2] = {1.0,-1.0};
            if(ort[eIter] == 0)
              {values[0] *= -1.0; values[1] *= -1.0;}

            // get LIDs associated with nodes
            int indices[2] = {nLIDs[tailIndex],nLIDs[headIndex]};

            // insert values in matrix
            int err = grad_matrix->ReplaceMyValues(eLIDs[eIter], 2, values, indices);
            if (err != 0)
              err = grad_matrix->InsertMyValues(eLIDs[eIter], 2, values, indices);
            TEUCHOS_ASSERT_EQUALITY(err,0);
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
  } else
    TEUCHOS_ASSERT(false);
  
}
