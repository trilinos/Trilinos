/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>

#include <fei_Lookup_Impl.hpp>

#include <fei_VectorSpace.hpp>

#include <snl_fei_Utils.hpp>
#include <fei_TemplateUtils.hpp>
#include <fei_CommUtils.hpp>

#include <snl_fei_Constraint.hpp>

#include <snl_fei_SubdMsgHandler.hpp>

#undef fei_file
#define fei_file "fei_Lookup_Impl.cpp"
#include <fei_ErrMacros.hpp>

//----------------------------------------------------------------------------
fei::Lookup_Impl::Lookup_Impl(fei::SharedPtr<fei::MatrixGraph> matGraph,
				  int nodeIDType)
  : matGraph_(matGraph),
    ptBlkMap_(NULL),
    vspace_(),
    nodeIDType_(nodeIDType),
    nodenumPairs_(),
    eqnnumPairs_(),
    nodenumSubdomainDB_(),
    databasesBuilt_(false),
    fieldIDs_(),
    fieldSizes_(),
    elemBlockIDs_(0, 4),
    fieldIDs_2D_(),
    workspace_()
{
  vspace_ = matGraph_->getRowSpace();
  ptBlkMap_ = vspace_->getPointBlockMap();

  int err = buildDatabases();
  if (err != 0) {
    voidERReturn;
  }
}

//----------------------------------------------------------------------------
fei::Lookup_Impl::~Lookup_Impl()
{
  fei::destroyValues(nodenumSubdomainDB_);
  nodenumSubdomainDB_.clear();
}

//----------------------------------------------------------------------------
int fei::Lookup_Impl::getEqnNumber(int nodeNumber, int fieldID)
{
  std::map<int,fei::Record*>::iterator
    nnp_iter = nodenumPairs_.find(nodeNumber);

  if (nnp_iter == nodenumPairs_.end()) return(-1);

  fei::Record* node = (*nnp_iter).second;

  std::vector<int>& eqnNums = vspace_->getEqnNumbers();
  int* eqnNumbers = eqnNums.size() > 0 ? &eqnNums[0] : NULL;
    + node->getOffsetIntoEqnNumbers();
  int numInstances = 1;
  int offset = -1;
  node->getFieldMask()->getFieldEqnOffset(fieldID, offset, numInstances);
  return(eqnNumbers[offset]);
}

//----------------------------------------------------------------------------
int fei::Lookup_Impl::getAssociatedNodeNumber(int eqnNumber)
{
  std::map<int,fei::Record*>::iterator
    enp_iter = eqnnumPairs_.find(eqnNumber);

  if (enp_iter == eqnnumPairs_.end()) return(-1);

  fei::Record* node = (*enp_iter).second;

  return( node->getNumber() );
}

//----------------------------------------------------------------------------
int fei::Lookup_Impl::getAssociatedNodeID(int eqnNumber)
{
  std::map<int,fei::Record*>::iterator
    enp_iter = eqnnumPairs_.find(eqnNumber);

  if (enp_iter == eqnnumPairs_.end()) return(-1);

  fei::Record* node = (*enp_iter).second;

  return( node->getID() );
}

bool fei::Lookup_Impl::isInLocalElement(int nodeNumber)
{
  std::map<int,std::vector<int>* >::iterator
    nns_iter = nodenumSubdomainDB_.find(nodeNumber);
  if (nns_iter != nodenumSubdomainDB_.end()) {
    return( true );
  }

  std::map<int,fei::Record*>::iterator
    nnp_iter = nodenumPairs_.find(nodeNumber);

  return(nnp_iter != nodenumPairs_.end() ? true : false);
}

//----------------------------------------------------------------------------
int fei::Lookup_Impl::getOffsetIntoBlkEqn(int blkEqn, int ptEqn)
{
  //assume blkEqn is a node-number, for now.
  std::map<int,fei::Record*>::iterator
    nnp_iter = nodenumPairs_.find(blkEqn);

  if (nnp_iter == nodenumPairs_.end()) return(-1);

  fei::Record* node = (*nnp_iter).second;

  int eqn = vspace_->getEqnNumbers()[node->getOffsetIntoEqnNumbers()];
  return(ptEqn - eqn);
}

//----------------------------------------------------------------------------
int fei::Lookup_Impl::buildDatabases()
{
  if (databasesBuilt_) return(0);

  snl_fei::RecordCollection* collection = NULL;
  int err = vspace_->getRecordCollection(nodeIDType_, collection);
  if (err != 0) {
    //probably means that vspace_ doesn't have 'nodeIDType_', so we'll skip the
    //rest of this function. It's not a problem if vspace_ doesn't have nodeIDType_.
    return(0);
  }

  std::vector<int>& vspcEqnNumbers = vspace_->getEqnNumbers();

  snl_fei::RecordCollection::map_type& rmap = collection->getRecords();

  snl_fei::RecordCollection::map_type::iterator
    r_iter = rmap.begin(),
    r_end = rmap.end();

  for(; r_iter != r_end; ++r_iter) {
    fei::Record* node = (*r_iter).second;

    std::pair<int,fei::Record* > int_node_pair(node->getNumber(), node);

    nodenumPairs_.insert(int_node_pair);

    int numEqns = node->getFieldMask()->getNumIndices();
    int* eqnNumbers = &vspcEqnNumbers[0];
                    + node->getOffsetIntoEqnNumbers();

    for(int eq=0; eq<numEqns; ++eq) {
      std::pair<int,fei::Record* > eqn_node_pair(eqnNumbers[eq], node);
      eqnnumPairs_.insert(eqn_node_pair);
    }
  }

  MPI_Comm comm = matGraph_->getRowSpace()->getCommunicator();

  int numLocalLagrangeConstraints = matGraph_->getLagrangeConstraints().size();

  int numGlobalLagrangeConstraints = 0;
  fei::GlobalSum(comm, numLocalLagrangeConstraints, numGlobalLagrangeConstraints);

  bool noconstraints = numGlobalLagrangeConstraints<1 ? true : false;

  fei::SharedIDs* subdomainIDs = NULL;
  fei::SharedIDs* sharedIDs = NULL;
  int i;
  err = vspace_->getSharedIDs_private(nodeIDType_, sharedIDs);
  if (err < 0) return(-1);

  if (noconstraints == false) {
    subdomainIDs = new fei::SharedIDs;

    snl_fei::SubdMsgHandler subdmsghndlr(collection, sharedIDs, subdomainIDs);

    int idx = snl_fei::binarySearch(nodeIDType_, vspace_->sharedIDTypes_);
    if (idx < 0) ERReturn(-1);

    if ((int)vspace_->ownerPatterns_.size() > idx &&
	(int)vspace_->sharerPatterns_.size() > idx) {
      subdmsghndlr.setSendPattern(vspace_->ownerPatterns_[idx]);
      subdmsghndlr.setRecvPattern(vspace_->sharerPatterns_[idx]);
      CHK_ERR( fei::exchange(comm, &subdmsghndlr) );
    }

    //Now the subdomainIDs object contains a mapping from each shared ID to a
    //list of processors that have that ID in their local subdomain.
    //So what we'll do next is run through the list of IDs in subdomainIDs and
    //for each ID, store the corresponding node-number in a database together
    //with a pointer to a list (vector) of the subdomain-processors.
  }
  else {
    subdomainIDs = sharedIDs;
  }

  int local_proc = fei::localProc(comm);

  fei::SharedIDs::table_type& sdIDTable = subdomainIDs->getSharedIDs();
  fei::SharedIDs::table_type::iterator
    sd_iter = sdIDTable.begin(),
    sd_end  = sdIDTable.end();

  for(i=0; sd_iter != sd_end; ++i, ++sd_iter) {
    int id = (*sd_iter).first;
    fei::SharedIDs::table_type::row_type* procList = (*sd_iter).second;

    fei::Record* node = collection->getRecordWithID(id);
    if (node == NULL) {
      ERReturn(-1);
    }

    std::vector<int>* newarray = new std::vector<int>;
    fei::SharedIDs::table_type::row_type::const_iterator
      p_iter = procList->begin(),
      p_end = procList->end();
    for(; p_iter != p_end; ++p_iter) {
      int proc = *p_iter;
      newarray->push_back(proc);
    }

    if (node->isInLocalSubdomain_) {
      snl_fei::sortedListInsert(local_proc, *newarray);
    }

    nodenumSubdomainDB_.insert(std::pair<int,std::vector<int>*>(node->getNumber(), newarray));
  }

  if (!noconstraints) {
    delete subdomainIDs;
  }

  databasesBuilt_ = true;
  return(0);
}

