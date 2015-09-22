/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _HexBeam_h_
#define _HexBeam_h_

#include <fei_base.hpp>

/**
HexBeam is a data generator used for testing the FEI. It generates data that
represents a cube of hex elements. More precisely, it's a cube that can be
elongated in one dimension, forming an arbitrarily long 3D bar. The dimensions
of the bar are WxWxD. In parallel runs, the bar is divided such that each
processor has roughly D/numprocs "slices" each with WxW elements.

The slices that are separated by processor boundaries share common nodes, which
are "shared nodes" in the global mesh.
 */
class HexBeam {
 public:
  enum { OneD = 0, TwoD = 1, ThreeD = 2 };

  HexBeam(int W, int D, int DofPerNode,
	  int decomp, int numProcs, int localProc);
  virtual ~HexBeam();

  virtual bool inErrorState() { return( inErrorState_ ); }

  virtual int numNodesPerElem() { return( nodesPerElem_ ); }

  virtual int numDofPerNode() { return( dofPerNode_ ); }

  virtual int numLocalElems() { return( localNumElems_ ); }

  virtual int numLocalNodes() { return( localNumNodes_ ); }

  virtual int firstLocalElem() { return( firstLocalElem_ ); }

  virtual int getElemConnectivity(int elemID, int* nodeIDs);

  virtual int getElemStiffnessMatrix(int elemID, double* elemMat);

  virtual int getElemLoadVector(int elemID, double* elemVec);

  virtual int getNumBCNodes();

  virtual int getBCNodes(int numNodes, int* nodeIDs);

  virtual int getBCValues(int numBCNodes, int* offsetsIntoField, double* vals);

  virtual int getNumSharedNodes();

  virtual int getSharedNodes(int numSharedNodes,
		     int*& sharedNodes,
		     int*& numSharingProcsPerNode,
		     int**& sharingProcs);

  virtual int getNumCRs() { return( 0 ); }

  virtual int getNumNodesPerCR() { return( 0 ); }

  virtual int getCRNodes(int** nodeIDs) { (void)nodeIDs; return(0); }

  int W_;
  int D_;
  int decomp_;
  int numProcs_;
  int localProc_;

  int totalNumElems_;
  int totalNumNodes_;
  int localNumElems_;
  int localNumNodes_;
  int firstLocalElem_;
  int firstLocalNode_;

  int numElemsPerSlice_;
  int numNodesPerSlice_;
  int numLocalSlices_;

  bool inErrorState_;

  int nodesPerElem_;
  int dofPerNode_;

  int numLocalDOF_;
  int numGlobalDOF_;
};

namespace HexBeam_Functions {

int init_elem_connectivities(FEI* fei, HexBeam& hexcube);

int init_shared_nodes(FEI* fei, HexBeam& hexcube);

int init_constraints(FEI* fei, HexBeam& hexcube, int& firstLocalCRID);

int load_elem_data(FEI* fei, HexBeam& hexcube);

int load_constraints(FEI* fei, HexBeam& hexcube, int firstLocalCRID);

int load_BC_data(FEI* fei, HexBeam& hexcube);

int print_cube_data(HexBeam& hexcube, int numProcs, int localProc);


int init_elem_connectivities(fei::MatrixGraph* matrixGraph, HexBeam& hexcube);

int init_shared_nodes(fei::MatrixGraph* matrixGraph, HexBeam& hexcube);

int init_constraints(fei::MatrixGraph* matrixGraph, HexBeam& hexcube,
		     int localProc, int& firstLocalCRID);

int init_slave_constraints(fei::MatrixGraph* matrixGraph, HexBeam& hexcube);

int load_elem_data(fei::MatrixGraph* matrixGraph,
		   fei::Matrix* mat,
		   fei::Vector* rhs,
		   HexBeam& hexcube);

int load_constraints(fei::LinearSystem* linSys, HexBeam& hexcube,
		     int firstLocalCRID);

int load_BC_data(fei::LinearSystem* linSys, HexBeam& hexcube);

}//namespace HexBeam_Functions

#endif // _HexBeam_h_
