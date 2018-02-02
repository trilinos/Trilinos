#include "Panzer_LOCPair_GlobalEvaluationData.hpp"
#include "Panzer_IntrepidOrientation.hpp"

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

  typedef panzer::BlockedTpetraLinearObjFactory<panzer::Traits,double,int,panzer::Ordinal64> BTLOF;
  typedef panzer::BlockedTpetraLinearObjContainer<double,int,panzer::Ordinal64> BTLOC;
  typedef panzer::UniqueGlobalIndexer<int,panzer::Ordinal64> UGI;
  typedef Tpetra::CrsMatrix<double,int,panzer::Ordinal64> Tmatrix;
  typedef Thyra::TpetraLinearOp<double,int,panzer::Ordinal64> ThTLOp;

  // use "AUXILIARY_EDGE" and "AUXILIARY_NODE" as the DOFs
  std::string edge_basis_name  = "AUXILIARY_EDGE";
  std::string nodal_basis_name = "AUXILIARY_NODE";
      
  // must be able to cast to a block linear object factory
  RCP<const BTLOF > blof  = rcp_dynamic_cast<const BTLOF >(linObjFactory,true);
  RCP<const panzer::BlockedDOFManager<int,panzer::Ordinal64> > blockedDOFMngr = blof->getGlobalIndexer();

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
        = rcp(new panzer::LOCPair_GlobalEvaluationData(blof,panzer::LinearObjContainer::Mat));
  RCP<panzer::LinearObjContainer> loc = rcp_dynamic_cast<panzer::LOCPair_GlobalEvaluationData>(dataObject,true)->getGhostedLOC();
  RCP<BTLOC> ghosted_tloc = rcp_dynamic_cast<BTLOC>(loc,true);
  RCP<panzer::LinearObjContainer> global_loc 
     = rcp_dynamic_cast<panzer::LOCPair_GlobalEvaluationData>(dataObject,true)->getGlobalLOC();
  RCP<BTLOC> global_tloc = rcp_dynamic_cast<BTLOC>(global_loc,true);

  // get the node-to-edge block from the ghosted LOC to fill as the gradient
  RCP<Thyra::BlockedLinearOpBase<double> > blocked_op = rcp_dynamic_cast<Thyra::BlockedLinearOpBase<double> >(ghosted_tloc->get_A_th(),true);
  RCP<ThTLOp> thyra_TpetraOp = rcp_dynamic_cast<ThTLOp>(blocked_op->getNonconstBlock(eBlockIndex,nBlockIndex),true);
  RCP<Tmatrix> grad_matrix = rcp_dynamic_cast<Tmatrix>(thyra_TpetraOp->getTpetraOperator(),true);
  grad_matrix->resumeFill();
 
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
      std::vector<panzer::Ordinal64> eGIDs;
      eUgi->getElementGIDs(elementIds[elemIter],eGIDs);
      std::vector<panzer::Ordinal64> nGIDs;
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
          int err = grad_matrix->replaceLocalValues(eLIDs[eIter], 2, values, indices);
          TEUCHOS_ASSERT_EQUALITY(err,2);
        }//end if
      }//end edge loop
    }//end element loop
  }//end element block loop
  grad_matrix->fillComplete();

  // get global linear object from ghosted linear object
  blof->ghostToGlobalContainer(*ghosted_tloc, *global_tloc,panzer::LinearObjContainer::Mat);
  RCP<Thyra::BlockedLinearOpBase<double> > global_blocked_op = rcp_dynamic_cast<Thyra::BlockedLinearOpBase<double> >(global_tloc->get_A_th(),true);
  RCP<Thyra::LinearOpBase<double> > thyra_gradient = global_blocked_op->getNonconstBlock(eBlockIndex,nBlockIndex);
 
  // add gradient callback to request handler
  reqHandler->addRequestCallback(Teuchos::rcp(new GradientRequestCallback(thyra_gradient)));
}
