/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/


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
  std::map<int,fei::Record<int>*>::iterator
    nnp_iter = nodenumPairs_.find(nodeNumber);

  if (nnp_iter == nodenumPairs_.end()) return(-1);

  fei::Record<int>* node = (*nnp_iter).second;

  std::vector<int>& eqnNums = vspace_->getEqnNumbers();
  int* eqnNumbers = eqnNums.size() > 0 ? &eqnNums[0] : NULL;
  if (eqnNumbers == NULL) {
    throw std::runtime_error("Fatal error in fei::Lookup_Impl::getEqnNumber");
  }
  eqnNumbers += node->getOffsetIntoEqnNumbers();
  int offset = -1;
  node->getFieldMask()->getFieldEqnOffset(fieldID, offset);
  return(eqnNumbers[offset]);
}

//----------------------------------------------------------------------------
int fei::Lookup_Impl::getAssociatedNodeNumber(int eqnNumber)
{
  std::map<int,fei::Record<int>*>::iterator
    enp_iter = eqnnumPairs_.find(eqnNumber);

  if (enp_iter == eqnnumPairs_.end()) return(-1);

  fei::Record<int>* node = (*enp_iter).second;

  return( node->getNumber() );
}

//----------------------------------------------------------------------------
int fei::Lookup_Impl::getAssociatedNodeID(int eqnNumber)
{
  std::map<int,fei::Record<int>*>::iterator
    enp_iter = eqnnumPairs_.find(eqnNumber);

  if (enp_iter == eqnnumPairs_.end()) return(-1);

  fei::Record<int>* node = (*enp_iter).second;

  return( node->getID() );
}

//----------------------------------------------------------------------------
int fei::Lookup_Impl::getAssociatedFieldID(int eqnNumber)
{
  std::map<int,fei::Record<int>*>::iterator
    enp_iter = eqnnumPairs_.find(eqnNumber);

  if (enp_iter == eqnnumPairs_.end()) return(-1);

  fei::Record<int>* node = (*enp_iter).second;

  fei::FieldMask* fm = node->getFieldMask();
  const std::vector<int>& fieldIDs = fm->getFieldIDs();
  const std::vector<int>& fieldSizes = fm->getFieldSizes();

  const std::vector<int>& eqnNumbers = vspace_->getEqnNumbers();

  int baseEqnOffset = node->getOffsetIntoEqnNumbers();
  int numNodalEqns = fm->getNumIndices();

  if (baseEqnOffset + numNodalEqns > (int)eqnNumbers.size()) {
    throw std::runtime_error("fei::Lookup_Impl::getAssociatedFieldID ERROR, nodal eqn offset out of range.");
  }

  int offset = 0;
  int eqn = eqnNumbers[baseEqnOffset];
  while(eqn < eqnNumber && offset < numNodalEqns) {
    eqn = eqnNumbers[baseEqnOffset + ++offset];
  }

  if (eqn != eqnNumber) {
    throw std::runtime_error("fei::Lookup_Impl::getAssociatedFieldID ERROR, eqnNumber not found");
  }

  int fieldSize_total = 0;
  for(size_t i=0; i<fieldSizes.size(); ++i) {
    fieldSize_total += fieldSizes[i];
    if (fieldSize_total > offset) {
      return fieldIDs[i];
    }
  }

  return -1;
}

//----------------------------------------------------------------------------
bool fei::Lookup_Impl::isInLocalElement(int nodeNumber)
{
  std::map<int,std::vector<int>* >::iterator
    nns_iter = nodenumSubdomainDB_.find(nodeNumber);
  if (nns_iter != nodenumSubdomainDB_.end()) {
    return( true );
  }

  std::map<int,fei::Record<int>*>::iterator
    nnp_iter = nodenumPairs_.find(nodeNumber);

  return(nnp_iter != nodenumPairs_.end() ? true : false);
}

//----------------------------------------------------------------------------
int fei::Lookup_Impl::getOffsetIntoBlkEqn(int blkEqn, int ptEqn)
{
  //assume blkEqn is a node-number, for now.
  std::map<int,fei::Record<int>*>::iterator
    nnp_iter = nodenumPairs_.find(blkEqn);

  if (nnp_iter == nodenumPairs_.end()) return(-1);

  fei::Record<int>* node = (*nnp_iter).second;

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

  std::vector<fei::Record<int> >& rvec = collection->getRecords();

  for(size_t i=0; i<rvec.size(); ++i) {
    fei::Record<int>* node = &rvec[i];

    std::pair<int,fei::Record<int>* > int_node_pair(node->getNumber(), node);

    nodenumPairs_.insert(int_node_pair);

    int numEqns = node->getFieldMask()->getNumIndices();
    int* eqnNumbers = &vspcEqnNumbers[0]
                    + node->getOffsetIntoEqnNumbers();

    for(int eq=0; eq<numEqns; ++eq) {
      std::pair<int,fei::Record<int>* > eqn_node_pair(eqnNumbers[eq], node);
      eqnnumPairs_.insert(eqn_node_pair);
    }
  }

  MPI_Comm comm = matGraph_->getRowSpace()->getCommunicator();

  int numLocalLagrangeConstraints = matGraph_->getLagrangeConstraints().size();

  int numGlobalLagrangeConstraints = 0;
  fei::GlobalSum(comm, numLocalLagrangeConstraints, numGlobalLagrangeConstraints);

  bool noconstraints = numGlobalLagrangeConstraints<1 ? true : false;

  fei::SharedIDs<int> subdomainIDs;
  fei::SharedIDs<int>& sharedIDs = vspace_->getSharedIDs(nodeIDType_);

  if (noconstraints == false) {
    snl_fei::SubdMsgHandler subdmsghndlr(collection, &sharedIDs, &subdomainIDs);

    if (vspace_->ownerPatterns_.size() > 0 && vspace_->sharerPatterns_.size() > 0) {
      subdmsghndlr.setSendPattern(vspace_->ownerPatterns_.find(nodeIDType_)->second);
      subdmsghndlr.setRecvPattern(vspace_->sharerPatterns_.find(nodeIDType_)->second);
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

  fei::SharedIDs<int>::map_type& sdIDTable = subdomainIDs.getSharedIDs();
  fei::SharedIDs<int>::map_type::iterator
    sd_iter = sdIDTable.begin(),
    sd_end  = sdIDTable.end();

  for(int i=0; sd_iter != sd_end; ++i, ++sd_iter) {
    int id = sd_iter->first;
    std::set<int>& procList = sd_iter->second;

    fei::Record<int>* node = collection->getRecordWithID(id);
    if (node == NULL) {
      ERReturn(-1);
    }

    std::vector<int>* newarray = new std::vector<int>;
    std::set<int>::const_iterator
      p_iter = procList.begin(), p_end = procList.end();

    for(; p_iter != p_end; ++p_iter) {
      int proc = *p_iter;
      newarray->push_back(proc);
    }

    if (node->isInLocalSubdomain_) {
      fei::sortedListInsert(local_proc, *newarray);
    }

    nodenumSubdomainDB_.insert(std::pair<int,std::vector<int>*>(node->getNumber(), newarray));
  }

  databasesBuilt_ = true;
  return(0);
}

