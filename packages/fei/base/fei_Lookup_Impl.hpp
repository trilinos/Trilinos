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


#ifndef _fei_Lookup_Impl_hpp_
#define _fei_Lookup_Impl_hpp_

#include <fei_macros.hpp>
#include <fei_MatrixGraph.hpp>
#include <fei_Pattern.hpp>
#include <fei_ConnectivityBlock.hpp>
#include <fei_SharedIDs.hpp>
#include <snl_fei_RecordCollection.hpp>
#include <snl_fei_PointBlockMap.hpp>
#include <fei_TemplateUtils.hpp>
#include <fei_Lookup.hpp>

#include <vector>

namespace fei {
  /** An implementation of the Lookup interface.
   */
  class Lookup_Impl : public virtual Lookup {
  public:
    /** Constructor */
    Lookup_Impl(fei::SharedPtr<fei::MatrixGraph> matGraph,
		int nodeIDType);

    /** Destructor */
    virtual ~Lookup_Impl();

    /** Implementation of Lookup:: method */
    int getNumFields()
      {
	return( vspace_->getNumFields() );
      }

    /** Implementation of Lookup:: method */
    int getFieldSize(int fieldID)
      {
	return(vspace_->getFieldSize(fieldID));
      }

    /** Implementation of Lookup:: method */
    const int* getFieldIDsPtr()
      {
	int len = getNumFields();
	if (len < 1) return(NULL);

	vspace_->getFields(fieldIDs_);
	return(&fieldIDs_[0]);
      }

    /** Implementation of Lookup:: method */
    const int* getFieldSizesPtr()
      {
	const int* fieldIDs = getFieldIDsPtr();
	if (fieldIDs == NULL) return(NULL);

	unsigned numFields = fieldIDs_.size();
	fieldSizes_.resize(numFields);
	int* fsPtr = &fieldSizes_[0];
	for(unsigned i=0; i<numFields; ++i) {
	  fsPtr[i] = vspace_->getFieldSize(fieldIDs[i]);
	}
	return(fsPtr);
      }

    /** Implementation of Lookup:: method */
    int getNumElemBlocks()
      { return(matGraph_->getConnectivityBlocks().size()); }

    /** Implementation of Lookup:: method */
    const GlobalID* getElemBlockIDs()
      {
	int err = matGraph_->getConnectivityBlockIDs(elemBlockIDs_);
	return(err==0 ? &elemBlockIDs_[0] : NULL);
      }

    /** Implementation of Lookup:: method */
    void getElemBlockInfo(GlobalID blockID,
                         int& interleaveStrategy, int& lumpingStrategy,
                         int& numElemDOF, int& numElements,
                         int& numNodesPerElem, int& numEqnsPerElem)
    {
      interleaveStrategy = 0; lumpingStrategy = 0;
      numElemDOF = 0;
      const fei::ConnectivityBlock* cblock =
        matGraph_->getConnectivityBlock(blockID);
      numElements = cblock->getConnectivityIDs().size();
      numNodesPerElem = cblock->getRowPattern()->getNumIDs();
      numEqnsPerElem = cblock->getRowPattern()->getNumIndices();
    }

    /** Implementation of Lookup:: method */
    const int* getNumFieldsPerNode(GlobalID blockID)
      {
	const fei::ConnectivityBlock* cblock = matGraph_->getConnectivityBlock(blockID);
	if (cblock==NULL) return(NULL);
	return(cblock->getRowPattern()->getNumFieldsPerID());
      }

    /** Implementation of Lookup:: method */
    const int* const* getFieldIDsTable(GlobalID blockID)
      {
	const fei::ConnectivityBlock* cblock = matGraph_->getConnectivityBlock(blockID);
	if (cblock==NULL) return(NULL);
	int numNodes = cblock->getRowPattern()->getNumIDs();
	const int* numFieldsPerNode = cblock->getRowPattern()->getNumFieldsPerID();
	const int* fieldIDs = cblock->getRowPattern()->getFieldIDs();
	fieldIDs_2D_.resize(numNodes);
	const int** f2dPtr = &fieldIDs_2D_[0];
	int offset = 0;
	for(int i=0; i<numNodes; ++i) {
	  f2dPtr[i] = fieldIDs + offset;
	  offset += numFieldsPerNode[i];
	}
	return(f2dPtr);
      }

    /** Implementation of Lookup:: method */
    int getEqnNumber(int nodeNumber, int fieldID);

    /** Implementation of Lookup:: method */
    int getAssociatedNodeNumber(int eqnNumber);

    int getAssociatedNodeID(int eqnNumber);

    /** Implementation of Lookup:: method */
    int getAssociatedFieldID(int eqnNumber);

    /** Implementation of Lookup:: method */
    bool isInLocalElement(int nodeNumber);

    /** Implementation of Lookup:: method */
    int getNumSubdomains(int nodeNumber)
    {
      std::vector<int>* subdomains = NULL;
      std::map<int,std::vector<int>* >::iterator
        nns_iter = nodenumSubdomainDB_.find(nodeNumber);
      if (nns_iter != nodenumSubdomainDB_.end()) subdomains = (*nns_iter).second;
      return( subdomains==0 ? 0 : subdomains->size() );
    }

    /** Implementation of Lookup:: method */
    int* getSubdomainList(int nodeNumber)
    {
      std::vector<int>* subdomains = NULL;
      std::map<int,std::vector<int>* >::iterator
        nns_iter = nodenumSubdomainDB_.find(nodeNumber);
      if (nns_iter != nodenumSubdomainDB_.end()) subdomains = (*nns_iter).second;

      return( subdomains==0 ? NULL : &(*subdomains)[0] );
    }

    /** Implementation of Lookup:: method */
    int getNumSharedNodes()
    {
      int numShared;
      int err = vspace_->getNumSharedIDs(nodeIDType_, numShared);
      return(err==0 ? numShared : -1);
    }

    /** Implementation of Lookup:: method */
    const int* getSharedNodeNumbers()
    {
      fei::SharedIDs<int>& sharedIDs = vspace_->getSharedIDs(nodeIDType_);

      int numShared = sharedIDs.getSharedIDs().size();
      workspace_.resize(numShared*2);
      int* wkPtr = &workspace_[0];
      fei::copyKeysToArray(sharedIDs.getSharedIDs(), numShared, wkPtr);

      snl_fei::RecordCollection* collection = NULL;
      vspace_->getRecordCollection(nodeIDType_, collection);

      for(int i=0; i<numShared; ++i) {
        fei::Record<int>* node = collection->getRecordWithID(wkPtr[i]);
        if (node == NULL) return NULL;

        wkPtr[numShared+i] = node->getNumber();
      }
      return(wkPtr+numShared);
    }

    /** Implementation of Lookup:: method */
    const int* getSharedNodeProcs(int nodeNumber)
    {
      std::map<int,fei::Record<int>*>::iterator
        nnp_iter = nodenumPairs_.find(nodeNumber);

      if (nnp_iter == nodenumPairs_.end()) return(0);

      fei::Record<int>* node = (*nnp_iter).second;

      const fei::SharedIDs<int>& sharedIDs = vspace_->getSharedIDs(nodeIDType_);

      int shID = node->getID();

      fei::SharedIDs<int>::map_type::const_iterator
        iter = sharedIDs.getSharedIDs().find(shID);
      if (iter == sharedIDs.getSharedIDs().end()) return(NULL);

      const std::set<int>& shprocs = iter->second;

      fei::copySetToVector(shprocs, workspace_);
      return(&workspace_[0]);
    }

    /** Implementation of Lookup:: method */
    int getNumSharingProcs(int nodeNumber)
    {
      std::map<int,fei::Record<int>*>::iterator
        nnp_iter = nodenumPairs_.find(nodeNumber);

      if (nnp_iter == nodenumPairs_.end()) return(0);

      fei::Record<int>* node = (*nnp_iter).second;

      const fei::SharedIDs<int>& sharedIDs = vspace_->getSharedIDs(nodeIDType_);

      int shID = node->getID();

      fei::SharedIDs<int>::map_type::const_iterator
        iter = sharedIDs.getSharedIDs().find(shID);
      if (iter == sharedIDs.getSharedIDs().end()) return(0);

      const std::set<int>& shprocs = iter->second;
      return(shprocs.size());
    }

    /** Implementation of Lookup:: method */
    bool isExactlyBlkEqn(int ptEqn)
      { return( ptBlkMap_->isExactlyBlkEqn(ptEqn) ); }

    /** Implementation of Lookup:: method */
    int ptEqnToBlkEqn(int ptEqn)
      { return( ptBlkMap_->eqnToBlkEqn(ptEqn) ); }

    /** Implementation of Lookup:: method */
    int getOffsetIntoBlkEqn(int blkEqn, int ptEqn);

    /** Implementation of Lookup:: method */
    int getBlkEqnSize(int blkEqn)
    {
      return( ptBlkMap_->getBlkEqnSize(blkEqn) );
    }

  private:
    int buildDatabases();

    fei::SharedPtr<fei::MatrixGraph> matGraph_;
    snl_fei::PointBlockMap* ptBlkMap_;
    fei::SharedPtr<fei::VectorSpace> vspace_;
    int nodeIDType_;

    std::map<int, fei::Record<int>*> nodenumPairs_;
    std::map<int,fei::Record<int>*> eqnnumPairs_;

    std::map<int,std::vector<int>*> nodenumSubdomainDB_;

    bool databasesBuilt_;

    std::vector<int> fieldIDs_;
    std::vector<int> fieldSizes_;
    std::vector<GlobalID> elemBlockIDs_;
    std::vector<const int*> fieldIDs_2D_;
    std::vector<int> workspace_;
  };//class Lookup_Impl
}//namespace fei

#endif // _fei_Lookup_Impl_hpp_
