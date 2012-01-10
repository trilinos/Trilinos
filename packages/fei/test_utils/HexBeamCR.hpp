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
