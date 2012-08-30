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


#ifndef _fei_Vector_core_hpp_
#define _fei_Vector_core_hpp_

#include <fei_iosfwd.hpp>
#include <fei_CSVec.hpp>
#include <fei_VectorSpace.hpp>
#include <fei_Reducer.hpp>
#include <fei_Logger.hpp>
#include <fei_SharedPtr.hpp>
#include <fei_EqnComm.hpp>

namespace fei {

/** Class to provide infrastructure for fei::Vector implementations. */
class Vector_core : protected fei::Logger {
 public:
  /** constructor */
  Vector_core(fei::SharedPtr<fei::VectorSpace> vecSpace, int numLocalEqns);

  /** destructor */
  virtual ~Vector_core();

  /** Retrieve a copy of values from the vector for the specified indices.
      Note that if the specified indices are not local in the underlying
      non-overlapping data decomposition, these values are not guaranteed to
      be correct until after the scatterToOverlap() method has been called.
  */
  int copyOut(int numValues,
	      const int* indices,
	      double* values,
	      int vectorIndex=0) const;

  /** Sum in data for FiniteElementData-specific structure. Power users only. */
  virtual int sumIntoFEVector(int blockID,
			      int connOffset,
			      int numNodes,
			      const int* nodeNumbers,
			      const int* numIndicesPerNode,
            const int* dof_ids,
			      const double* values) = 0;

  /** Another FiniteElementData-specific method. Power users only. */
  virtual int copyOut_FE(int nodeNumber, int dofOffset, double& value) = 0;

  /** Give specified data to underlying vector object. */
  int giveToVector(int numValues,
		   const int* indices,
		   const double* values,
		   bool sumInto=true,
		   int vectorIndex=0);

  /** Scatter locally-owned vector data to remote sharing procs. */
  virtual int scatterToOverlap();

  /** define the extent of the overlapping indices. If the optional arguments
    are not specified, remote eqns will be obtained from the internal
    fei::VectorSpace attribute. */
  void setOverlap(int numRemoteEqns=0, const int* remoteEqns=NULL);

 protected:
  /** Assemble data specified by fieldID, idType, etc. */
  int assembleFieldData(int fieldID,
			int idType,
			int numIDs,
			const int* IDs,
			const double* data,
			bool sumInto=true,
			int vectorIndex=0);

  int assembleFieldDataLocalIDs(int fieldID,
			int idType,
			int numIDs,
			const int* localIDs,
			const double* data,
			bool sumInto=true,
			int vectorIndex=0);

  void setCommSizes();
  /** Gather shared data from remote procs for eqns that are locally-owned. */
  virtual int gatherFromOverlap(bool accumulate = true);

  /** Copy out data specified by fieldID, idType, etc. */
  virtual int copyOutFieldData(int fieldID,
			       int idType,
			       int numIDs,
			       const int* IDs,
			       double* data,
			       int vectorIndex=0);

  /** Review this function. Is it redundant with other functions? */
  virtual int giveToUnderlyingVector(int numValues,
				     const int* indices,
				     const double* values,
				     bool sumInto=true,
				     int vectorIndex=0) = 0;

  /** Review this function. Is it redundant with other functions? */
  virtual int copyOutOfUnderlyingVector(int numValues,
					const int* indices,
					double* values,
					int vectorIndex=0) const = 0;

  /** Establish basic information like sizes etc. */
  /** Write data to specified filename. */
  virtual int writeToFile(const char* filename,
			  bool matrixMarketFormat=true);

  /** Write data to specified ostream. */
  virtual int writeToStream(FEI_OSTREAM& ostrm,
			    bool matrixMarketFormat=true);

  /** Return vector-space that describes that size/layout of this vector. */
  fei::SharedPtr<fei::VectorSpace> get_vector_space() const
    {
      return( vecSpace_ );
    }

  /** Set vector-space that describes the size/layout of this vector. */
  void set_vector_space(fei::SharedPtr<fei::VectorSpace> vspace)
    {
      vecSpace_ = vspace;
    }

  /** Query for first locally-owned vector position. */
  int firstLocalOffset() const { return( firstLocalOffset_ ); }

  /** Query for last locally-owned vector position. */
  int lastLocalOffset() const { return( lastLocalOffset_ ); }

  /** work_indices */
  std::vector<int>& work_indices() { return( work_indices_ ); }
  /** work_indices2 */
  std::vector<int>& work_indices2(){ return( work_indices2_); }

  /** haveFEVector */
  bool haveFEVector() { return( haveFEVector_ ); }
  /** setFEVector */
  void setFEVector(bool flag) {haveFEVector_ = flag; }

  /** remotelyOwned */
  std::vector<CSVec*>& remotelyOwned() { return( remotelyOwned_ ); }
  const std::vector<CSVec*>& remotelyOwned() const { return( remotelyOwned_ ); }
  std::vector<int>& remotelyOwnedProcs() { return( remotelyOwnedProcs_ ); }
  const std::vector<int>& remotelyOwnedProcs() const { return( remotelyOwnedProcs_ ); }

  fei::CSVec* getRemotelyOwned(int proc) {
    std::vector<int>::iterator iter = std::lower_bound(remotelyOwnedProcs_.begin(), remotelyOwnedProcs_.end(), proc);
    fei::CSVec* return_vec = NULL;
    size_t offset = iter - remotelyOwnedProcs_.begin();
    if (iter == remotelyOwnedProcs_.end() || *iter != proc) {
      remotelyOwnedProcs_.insert(iter, proc);
      return_vec = new fei::CSVec;
      remotelyOwned_.insert(remotelyOwned_.begin()+offset, return_vec);
    }
    else {
      return_vec = remotelyOwned_[offset];
    }
 
    return return_vec;
  }

  const fei::CSVec* getRemotelyOwned(int proc) const {
    std::vector<int>::const_iterator iter = std::lower_bound(remotelyOwnedProcs_.begin(), remotelyOwnedProcs_.end(), proc);
    if (iter == remotelyOwnedProcs_.end() || *iter != proc) {
      throw std::runtime_error("failed to find remote-vec for specified processor.");
    }

    size_t offset = iter - remotelyOwnedProcs_.begin();
    return remotelyOwned_[offset];
  }

 protected:
  fei::SharedPtr<fei::EqnComm> eqnComm_;

 private:
  void pack_send_buffers(const std::vector<int>& sendProcs,
                       const std::vector<fei::CSVec*>& remotelyOwned,
                       std::vector<std::vector<char> >& send_chars,
                       bool resize_buffer,
                       bool zeroRemotelyOwnedAfterPacking);

  fei::SharedPtr<fei::VectorSpace> vecSpace_;

  MPI_Comm comm_;

  int firstLocalOffset_, lastLocalOffset_, numLocal_;

  std::vector<int> work_indices_;
  std::vector<int> work_indices2_;

  bool haveFEVector_;

  std::vector<int> remotelyOwnedProcs_;
  std::vector<CSVec*> remotelyOwned_;
  std::vector<int> sendProcs_;
  std::vector<int> recvProcs_;
  std::vector<int> recv_sizes_;
  std::vector<std::vector<char> > recv_chars_;
  std::vector<std::vector<char> > send_chars_;
  bool sendRecvProcsNeedUpdated_;

  bool overlapAlreadySet_;
  std::string dbgprefix_;
};//class Vector_core

}//namespace fei

#endif

