/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _snl_fei_Broker_FEData_hpp_
#define _snl_fei_Broker_FEData_hpp_

#include <fei_macros.hpp>
#include <fei_mpi.h>
#include <snl_fei_Broker.hpp>
#include <fei_FiniteElementData.hpp>
#include <fei_VectorSpace.hpp>
#include <fei_MatrixGraph.hpp>
#include <fei_Matrix_Impl.hpp>
#include <fei_Pattern.hpp>
#include <fei_Vector_Impl.hpp>
#include <fei_ConnectivityBlock.hpp>
#include <snl_fei_LinearSystem_FEData.hpp>
#include <fei_Lookup_Impl.hpp>

#undef fei_file
#define fei_file "snl_fei_Broker_FEData.hpp"
#include <fei_ErrMacros.hpp>

namespace snl_fei {

  /** Implementation of snl_fei::Broker specialized to broker objects from a
      FiniteElementData instance.
  */
  class Broker_FEData : public snl_fei::Broker {
  public:
    /** Constructor */
    Broker_FEData(fei::SharedPtr<FiniteElementData> feData,
		  fei::SharedPtr<fei::MatrixGraph> matrixGraph,
		  int nodeIDType);

    /** destructor */
    virtual ~Broker_FEData();

    /** Produce an instance of an fei::Vector. This overloading of the
	create() method is for use by Broker implementations that are
	dispensing 'views' of vectors that reside in LinearSystemCore or
	FiniteElementData container implementations. In those cases, there is
	a distinction that must be made between solution-vectors and
	rhs-vectors.

	@param isSolutionVector
     */
    virtual fei::SharedPtr<fei::Vector> createVector(bool isSolutionVector=false)
      {
	fei::SharedPtr<fei::Vector> vptr;
	if (matrixGraph_.get() == NULL) return(vptr);

	if (setStructure() != 0) return(vptr);

        int localsize = matrixGraph_->getRowSpace()->getNumIndices_Owned();
	fei::SharedPtr<fei::Vector> vecptr;
        vecptr.reset(new fei::Vector_Impl<FiniteElementData>(matrixGraph_->getRowSpace(),
                                                   feData_.get(), localsize,
                                                    isSolutionVector));
	return(vecptr);
      }

    /** Produce an instance of an fei::Matrix
     */
    virtual fei::SharedPtr<fei::Matrix> createMatrix()
    {
      fei::SharedPtr<fei::Matrix> mptr;
      if (matrixGraph_.get() == NULL) return(mptr);

      if (setStructure() != 0) return(mptr);
      int localsize = matrixGraph_->getRowSpace()->getNumIndices_Owned();

      bool zeroSharedRows = false;
      fei::SharedPtr<fei::Matrix> matptr;
      matptr.reset(new fei::Matrix_Impl<FiniteElementData>(feData_,
            matrixGraph_, localsize, zeroSharedRows));
      return(matptr);
    }

    /** Produce an instance of an fei::LinearSystem
     */
    virtual fei::SharedPtr<fei::LinearSystem> createLinearSystem()
      {
	fei::SharedPtr<fei::LinearSystem> lsptr;
	if (matrixGraph_.get() == NULL) return(lsptr);

	if (setStructure() != 0) return(lsptr);

	snl_fei::LinearSystem_FEData*
	  linsysfed = new LinearSystem_FEData(feData_,
					      matrixGraph_);
	linsysfed->setLookup(lookup_);
	fei::SharedPtr<fei::LinearSystem> linsysptr(linsysfed);
	return(linsysptr);
      }

    /** Set the MatrixGraph object used by this broker. */
    virtual void setMatrixGraph(fei::SharedPtr<fei::MatrixGraph> matrixGraph)
    {
      matrixGraph_ = matrixGraph;
    }

  private:
    int setStructure()
    {
      if (matrixGraph_.get() == NULL) ERReturn(-1);
      if (setStructure_ == true) return(0);

      lookup_ = new fei::Lookup_Impl(matrixGraph_, nodeIDType_);

      CHK_ERR( feData_->setLookup(*lookup_) );

      fei::SharedPtr<fei::VectorSpace> vspace = matrixGraph_->getRowSpace();

      int numLocalNodes = vspace->getNumOwnedAndSharedIDs(nodeIDType_);

      int numElemBlocks = matrixGraph_->getNumConnectivityBlocks();

      int nodeType = 0;
      snl_fei::RecordCollection* nodeRecords = NULL;
      vspace->getRecordCollection(nodeType, nodeRecords);

      int* intData = new int[numElemBlocks*3];
      int* numElemsPerBlock =       intData;
      int* numNodesPerElem =        intData+numElemBlocks;
      int* elemMatrixSizePerBlock = intData+2*numElemBlocks;
      int i;

      std::vector<int> elemBlockIDs;
      CHK_ERR( matrixGraph_->getConnectivityBlockIDs( elemBlockIDs) );

      for(i=0; i<numElemBlocks; ++i) {
        const fei::ConnectivityBlock* cblock =
          matrixGraph_->getConnectivityBlock(elemBlockIDs[i]);
        if (cblock==NULL) return(-1);
        numElemsPerBlock[i] = cblock->getConnectivityIDs().size();
        numNodesPerElem[i] = cblock->getRowPattern()->getNumIDs();
        elemMatrixSizePerBlock[i] = cblock->getRowPattern()->getNumIndices();
      }

      int numSharedNodes = 0;
      CHK_ERR( vspace->getNumSharedIDs(nodeIDType_, numSharedNodes) );

      int numLagrangeConstraints = matrixGraph_->getLocalNumLagrangeConstraints();

      CHK_ERR( feData_->describeStructure(numElemBlocks,
            numElemsPerBlock,
            numNodesPerElem,
            elemMatrixSizePerBlock,
            numLocalNodes,
            numSharedNodes,
            numLagrangeConstraints) );

      std::map<int,fei::ConnectivityBlock*>::const_iterator
        cdb_iter = matrixGraph_->getConnectivityBlocks().begin();

      std::vector<int> nodeNumbers, numDofPerNode, dof_ids;
      int total_num_dof = 0;
      for(i=0; i<numElemBlocks; ++i, ++cdb_iter) {
        fei::ConnectivityBlock* cblock = (*cdb_iter).second;
        fei::Pattern* pattern = cblock->getRowPattern();

        int numConnectedNodes = pattern->getNumIDs();
        nodeNumbers.resize(numConnectedNodes);
        numDofPerNode.resize(numConnectedNodes);
        int* nodeNumPtr = &nodeNumbers[0];
        int* numDofPtr = &numDofPerNode[0];

        //For the calls to FiniteElementData::setConnectivity, we're going to
        //need a list of num-dof-per-node. So construct that now.
        const int* numFieldsPerID = pattern->getNumFieldsPerID();
        const int* fieldIDs = pattern->getFieldIDs();

        int foffset = 0;
        for(int ii=0; ii<numConnectedNodes; ++ii) {
          int dof = 0;
          for(int f=0; f<numFieldsPerID[ii]; ++f) {
            dof += vspace->getFieldSize(fieldIDs[foffset++]);
          }
          numDofPtr[ii] = dof;
          total_num_dof += dof;
        }

        dof_ids.resize(total_num_dof, 0);
        int* dof_ids_ptr = &dof_ids[0];

        ////
        //Next we'll loop over the connectivity-lists in this block,
        //making a call to FiniteElementData::setConnectivity for each one.

        std::map<int,int>& elemIDs = cblock->getConnectivityIDs();
        int numElems = elemIDs.size();
        int* nodes = &(cblock->getRowConnectivities()[0]);

        int offset = 0;
        for(int elem=0; elem<numElems; ++elem) {
          for(int n=0; n<numConnectedNodes; ++n) {
            fei::Record<int>* node = nodeRecords->getRecordWithLocalID(nodes[offset++]);
            nodeNumPtr[n] = node->getNumber();
          }

          CHK_ERR( feData_->setConnectivity(elemBlockIDs[i], elem,
                numConnectedNodes,
                nodeNumPtr, numDofPtr, dof_ids_ptr));
        }//end for(...numElems...)
      }//end for(...numElemBlocks...)

      delete [] intData;

      setStructure_ = true;
      return(0);
    }

    fei::SharedPtr<FiniteElementData> feData_;
    fei::SharedPtr<fei::MatrixGraph> matrixGraph_;

    int nodeIDType_;
    bool setStructure_;
    bool setMatrixMatrixGraph_;

    fei::Lookup_Impl* lookup_;
  };//class Broker_FEData
}//namespace snl_fei


#endif // _snl_fei_Broker_FEData_hpp_
