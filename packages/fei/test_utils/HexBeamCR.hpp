/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _HexBeamCR_hpp_
#define _HexBeamCR_hpp_

#include "fei_macros.hpp"

#include "test_utils/HexBeam.hpp"

/**
HexBeamCR is a data generator used for testing the FEI. It generates data that
represents a beam of hex elements. More precisely, it's a cube that can be
elongated in one dimension, forming an arbitrarily long 3D beam. The dimensions
of the beam are WxWxD. In parallel runs, the beam is divided such that each
processor has roughly D/numprocs "slices" each with WxW elements.

HexBeamCR is the same as HexBeam, except that when slices are made to form
a parallel decomposition, the slices are bound together using
constraint-relations. i.e., the slices that are separated by processor
boundaries do not share common nodes as they do in the case of HexBeam. Instead
the opposing faces have distinct (non-shared) nodes which are tied together
by the constraints.
In addition, HexBeamCR has a constraint-bound slice in the middle of each
processor's partition so that there are constraints even in the case of
a serial run on a single processor.
 */
class HexBeamCR : public HexBeam {
 public:
  enum { OneD = 0, TwoD = 1, ThreeD = 2 };

  HexBeamCR(int W, int D, int DofPerNode,
	  int decomp, int numProcs, int localProc);
  virtual ~HexBeamCR();

  int getElemConnectivity(int elemID, int* nodeIDs);

  int getElemStiffnessMatrix(int elemID, double* elemMat);

  int getElemLoadVector(int elemID, double* elemVec);

  int getNumBCNodes();

  int getBCNodes(int numNodes, int* nodeIDs);

  int getBCGammaValues(int numBCDofs, double* gamma);

  int getNumSharedNodes();

  int getSharedNodes(int numSharedNodes,
		     int*& sharedNodes,
		     int*& numSharingProcsPerNode,
		     int**& sharingProcs);

  int getNumCRs() { return( numLocalCRs_ ); }

  int getNumNodesPerCR(){ return( numNodesPerCR_ ); }

  int getCRNodes(int** nodeIDs);

  int numLocalCRs_;
  int numNodesPerCR_;

  int firstLocalSlice_;
  int localCRslice_;
};

#endif // _HexBeamCR_hpp_
