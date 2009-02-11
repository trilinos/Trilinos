/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

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
    int getAssociatedFieldID(int eqnNumber)
      { return(-1); }

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
	fei::SharedIDs* sharedIDs = NULL;
	int err = vspace_->getSharedIDs_private(nodeIDType_, sharedIDs);
	if (err < 0) return(NULL);

	int numShared = sharedIDs->getSharedIDs().getMap().size();
	workspace_.resize(numShared*2);
	int* wkPtr = &workspace_[0];
	fei::copyKeysToArray(sharedIDs->getSharedIDs().getMap(), numShared, wkPtr);

	snl_fei::RecordCollection* collection = NULL;
	err = vspace_->getRecordCollection(nodeIDType_, collection);
	if (err < 0) return(NULL);

	for(int i=0; i<numShared; ++i) {
	  fei::Record* node = NULL;
	  try {
	    node = collection->getRecordWithID(wkPtr[i]);
	  }
	  catch (std::runtime_error& exc) {
	    return(NULL);
	  }
	  wkPtr[numShared+i] = node->getNumber();
	}
	return(wkPtr+numShared);
      }

    /** Implementation of Lookup:: method */
    const int* getSharedNodeProcs(int nodeNumber)
      {
	std::map<int,fei::Record*>::iterator
	  nnp_iter = nodenumPairs_.find(nodeNumber);

	if (nnp_iter == nodenumPairs_.end()) return(0);

	fei::Record* node = (*nnp_iter).second;

	fei::SharedIDs* sharedIDs = NULL;
	int err = vspace_->getSharedIDs_private(nodeIDType_, sharedIDs);
	if (err < 0) return(NULL);

	int shID = node->getID();

	fei::SharedIDs::table_type::row_type* list = sharedIDs->getSharedIDs().getRow(shID);
	if (list == NULL) return(NULL);

        workspace_.resize(list->size());
	list->copy_to_array(workspace_.size(), &workspace_[0]);
	return(&workspace_[0]);
      }

    /** Implementation of Lookup:: method */
    int getNumSharingProcs(int nodeNumber)
      {
	std::map<int,fei::Record*>::iterator
	  nnp_iter = nodenumPairs_.find(nodeNumber);

	if (nnp_iter == nodenumPairs_.end()) return(0);

	fei::Record* node = (*nnp_iter).second;

	fei::SharedIDs* sharedIDs = NULL;
	int err = vspace_->getSharedIDs_private(nodeIDType_, sharedIDs);
	if (err < 0) return(err);

	int shID = node->getID();

	fei::SharedIDs::table_type::row_type* list = sharedIDs->getSharedIDs().getRow(shID);
	return(list!=NULL ? list->size() : -1);
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

    std::map<int, fei::Record*> nodenumPairs_;
    std::map<int,fei::Record*> eqnnumPairs_;

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
