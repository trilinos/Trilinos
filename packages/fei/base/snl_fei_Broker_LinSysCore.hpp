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


#ifndef _snl_fei_Broker_LinSysCore_hpp_
#define _snl_fei_Broker_LinSysCore_hpp_

#include <fei_macros.hpp>
#include <fei_mpi.h>
#include <fei_CommUtils.hpp>
#include <snl_fei_Broker.hpp>
#include <fei_LinearSystemCore.hpp>
#include <fei_VectorSpace.hpp>
#include <fei_Lookup_Impl.hpp>
#include <fei_MatrixGraph.hpp>
#include <fei_SparseRowGraph.hpp>
#include <fei_Vector_Impl.hpp>
#include <fei_Matrix_Impl.hpp>
#include <fei_MatrixReducer.hpp>
#include <fei_VectorReducer.hpp>
#include <fei_Reducer.hpp>
#include <snl_fei_LinearSystem_General.hpp>

#undef fei_file
#define fei_file "snl_fei_Broker_LinSysCore.hpp"
#include <fei_ErrMacros.hpp>

namespace snl_fei {

  /** Implementation of snl_fei::Broker specialized to broker objects from a
      LinearSystemCore instance.
  */
  class Broker_LinSysCore : public snl_fei::Broker {
  public:
    /** Constructor */
    Broker_LinSysCore(fei::SharedPtr<LinearSystemCore> lsc,
                      fei::SharedPtr<fei::MatrixGraph> matrixGraph,
                      fei::SharedPtr<fei::Reducer> reducer,
                      bool blockMatrix);

    /** destructor */
    virtual ~Broker_LinSysCore();

    /** Produce an instance of an fei::Vector. This overloading of the
        createVector() method is for use by Broker implementations that are
        dispensing 'views' of vectors that reside in LinearSystemCore or
        FiniteElementData container implementations. In those cases, there is
        a distinction that must be made between solution-vectors and
        rhs-vectors.

        @param isSolutionVector
     */
    fei::SharedPtr<fei::Vector> createVector(bool isSolutionVector=false)
      {
        fei::SharedPtr<fei::Vector> vptr;
        if (matrixGraph_.get() == NULL) return(vptr);
        if (setGlobalOffsets() != 0) return(vptr);

        fei::SharedPtr<fei::Vector> tmpvec;
        tmpvec.reset(new fei::Vector_Impl<LinearSystemCore >(matrixGraph_->getRowSpace(),
                                                       linsyscore_.get(),
                                                       numLocalEqns_,
                                                       isSolutionVector));

        fei::SharedPtr<fei::Vector> vec;
        if (reducer_.get() != NULL) {
          vec.reset(new fei::VectorReducer(reducer_, tmpvec, isSolutionVector));
        }
        else {
          vec = tmpvec;
        }

        return(vec);
      }

    /** Produce an instance of an fei::Matrix
     */
    fei::SharedPtr<fei::Matrix> createMatrix()
      {
        fei::SharedPtr<fei::Matrix> mptr;
        if (matrixGraph_.get() == NULL) return(mptr);

        if (setMatrixStructure() != 0) return(mptr);

        bool zeroSharedRows = true;
        if (reducer_.get() != NULL) {
          zeroSharedRows = false;
        }

        fei::SharedPtr<fei::Matrix> tmpmat;
        tmpmat.reset(new fei::Matrix_Impl<LinearSystemCore>(linsyscore_, matrixGraph_,
                                                           numLocalEqns_, zeroSharedRows));
        fei::SharedPtr<fei::Matrix> matptr;
        if (reducer_.get() != NULL) {
          matptr.reset(new fei::MatrixReducer(reducer_, tmpmat));
        }
        else {
          matptr = tmpmat;
        }

        return(matptr);
      }

    /** Produce an instance of an fei::LinearSystem
     */
    fei::SharedPtr<fei::LinearSystem> createLinearSystem()
      {
        fei::SharedPtr<fei::LinearSystem> lsptr;
        if (matrixGraph_.get() == NULL) return(lsptr);

        if (setMatrixStructure() != 0) return(lsptr);

        fei::SharedPtr<fei::LinearSystem>
          linSys(new snl_fei::LinearSystem_General(matrixGraph_));
        return(linSys);
      }

    /** Set the MatrixGraph object used by this broker. */
    void setMatrixGraph(fei::SharedPtr<fei::MatrixGraph> matrixGraph)
      { matrixGraph_ = matrixGraph; }

  private:
    int setGlobalOffsets()
      {
        //only set the global offsets once.
        if (setGlobalOffsets_) return(0);

        if (matrixGraph_.get() == NULL) return(-1);

        MPI_Comm comm = matrixGraph_->getRowSpace()->getCommunicator();
        int num_procs = fei::numProcs(comm);
        int local_proc = fei::localProc(comm);

        std::vector<int> globalOffsets;
        std::vector<int> globalBlkOffsets;

        if (reducer_.get() != NULL) {
          int localsize = reducer_->getLocalReducedEqns().size();
          numLocalEqns_ = localsize;
          std::vector<int> lsizes(num_procs, 0);
          std::vector<int> gsizes(num_procs, 0);
          lsizes[local_proc] = localsize;
          fei::GlobalMax(comm, lsizes, gsizes);
          globalOffsets.resize(num_procs+1);
          int offset = 0;
          for(int p=0; p<num_procs; ++p) {
            globalOffsets[p] = offset;
            offset += gsizes[p];
          }
          globalOffsets[num_procs] = offset;
          globalBlkOffsets = globalOffsets;
        }
        else {
          fei::SharedPtr<fei::VectorSpace> vecSpace = 
            matrixGraph_->getRowSpace();

          vecSpace->getGlobalIndexOffsets(globalOffsets);
          vecSpace->getGlobalBlkIndexOffsets(globalBlkOffsets);

          numLocalEqns_ = globalOffsets[local_proc+1]-globalOffsets[local_proc];
        }

        CHK_ERR(linsyscore_->setGlobalOffsets(num_procs+1,
                                              &globalBlkOffsets[0],
                                              &globalOffsets[0],
                                              &globalBlkOffsets[0]));

        setGlobalOffsets_ = true;
        return(0);
      }

    int setMatrixStructure()
      {
        //only set the matrix matrixGraph once.
        if (setMatrixStructure_) return(0);

        if (matrixGraph_.get() == NULL) return(-1);

        lookup_ = new fei::Lookup_Impl(matrixGraph_, 0);

        CHK_ERR( linsyscore_->setLookup(*lookup_) );

        CHK_ERR( setGlobalOffsets() );

        MPI_Comm comm = matrixGraph_->getRowSpace()->getCommunicator();

        fei::SharedPtr<fei::SparseRowGraph> localSRGraph =
          matrixGraph_->createGraph(blockMatrix_);

        std::vector<int>& rowNumbers = localSRGraph->rowNumbers;
        int numLocalRows = rowNumbers.size();
        int* rowOffsets = &(localSRGraph->rowOffsets[0]);
        int numLocalNonzeros = localSRGraph->packedColumnIndices.size();
        int* nonzeros = &(localSRGraph->packedColumnIndices[0]);

        int numGlobalNonzeros = 0;
        fei::GlobalSum(comm, numLocalNonzeros, numGlobalNonzeros);

        std::vector<int*> colPtrs(numLocalRows);
        std::vector<int> ptRowsPerBlkRow(numLocalRows, 1);
        std::vector<int> rowLengths(numLocalRows);
        int* rowLengthsPtr = &rowLengths[0];

        for(int i=0; i<numLocalRows; ++i) {
          colPtrs[i] = &(nonzeros[rowOffsets[i]]);
          rowLengthsPtr[i] = rowOffsets[i+1]-rowOffsets[i];
          if (blockMatrix_ == true) {
            ptRowsPerBlkRow[i] = lookup_->getBlkEqnSize(rowNumbers[i]);
          }
        }

        CHK_ERR( linsyscore_->setMatrixStructure(&colPtrs[0],
                                                 rowLengthsPtr,
                                                 &colPtrs[0],
                                                 rowLengthsPtr,
                                                 &ptRowsPerBlkRow[0]));

        setMatrixStructure_ = true;

        return(0);
      }


    fei::SharedPtr<LinearSystemCore> linsyscore_;
    fei::SharedPtr<fei::MatrixGraph> matrixGraph_;
    fei::SharedPtr<fei::Reducer> reducer_;
    Lookup* lookup_;

    bool setGlobalOffsets_;
    int numLocalEqns_;
    bool setMatrixStructure_;
    bool blockMatrix_;
  };//class Broker_LinSysCore
}//namespace snl_fei

#endif // _snl_fei_Broker_LinSysCore_hpp_
