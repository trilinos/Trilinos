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

#ifndef _fei_AztecDVBR_Matrix_hpp_
#define _fei_AztecDVBR_Matrix_hpp_


#include <fei_SharedPtr.hpp>
//
// This class provides a wrapper for the Aztec DVBR matrix data structure.
//
// Some functions return an int. This will be the ESI_error_code. It
// will be 0 if there were no errors, 1 if an error occurred.
//
namespace fei_trilinos {

class Aztec_BlockMap;

class AztecDVBR_Matrix {
 
  public:
    // Constructor.
    AztecDVBR_Matrix(fei::SharedPtr<Aztec_BlockMap> map);

    //Copy constructor.
    AztecDVBR_Matrix(const AztecDVBR_Matrix& src);

    virtual ~AztecDVBR_Matrix ();

    //query functions.

    int getNumBlocksPerRow(int blkRow, int& nnzBlksPerRow) const;
    int getNumNonzerosPerRow(int blkRow, int& nnzPerRow) const;
    int getNumBlocksPerRow(int* nnzBlksPerRow) const;
    int getNumNonzerosPerRow(int* nnzPerRow) const;
    int getBlockSize(int blkRow, int blkCol, int& ptRows, int& ptCols);

    // Mathematical functions.
    void matvec(const Aztec_LSVector& x, Aztec_LSVector& y) const;

    void put(double s);

    // ... to read matrix.

    int getBlockRow(int blk_row,
                    double* vals,
                    int* blk_col_inds,
                    int num_nz_blocks) const;
    
    // ... to write matrix.

    int putBlockRow(int blk_row,
                    double* vals,
                    int* blk_col_inds,
                    int num_nz_blocks) const;

    int sumIntoBlockRow(int blk_row,
                        double* vals,
                        int* blk_col_inds,
                        int num_nz_blocks) const;
 
    // configuration function.
    void allocate(int* num_nz_blocks, int* blk_col_inds);

    void loadComplete();

    bool isLoaded() const {return(isLoaded_);};
    void setLoaded(bool flag) {isLoaded_ = flag;};
    bool isAllocated() const {return(isAllocated_);};
    void setAllocated(bool flag) {isAllocated_ = flag;};

    AZ_MATRIX* getAZ_MATRIX_Ptr() const {return(Amat_);};

    bool readFromFile(const char *filename);
    bool writeToFile(const char *fileName) const;

    //numRemoteBlocks is the number of block-column-indices on this processor
    //that correspond to block-rows that reside on another processor.
    int getNumRemoteBlocks() {return(numRemoteBlocks_);};
    int* getRemoteBlockIndices() {return(remoteInds_);};
    int* getRemoteBlockSizes() {return(remoteBlockSizes_);};

    int* getUpdate_index() {return(update_index_);};
    int* getData_org() {return(data_org_);};

  private:
    int inUpdate(int globalIndex, int& localIndex) const;

    void readAllocateInfo(FILE* infile, int*& num_nz_blocks, int*& blk_col_inds);
    void readMatrixData(FILE* infile);

    void calcRpntr();
    void calcBpntr(int* nzBlksPerRow);
    void setBindx(int nnzBlks, int* blkColInds);
    void calcIndx(int nnzBlks);

    int getBindxOffset(int blkInd, int bpntrStart, int bpntrEnd) const;

    void calcRemoteInds(int*& remoteInds, int& len);
    void getRemoteBlkSizes(int* remoteBlkSizes, int* remoteInds, int len);
    void insertList(int item, int*& list, int& len);
    void getValuesFromString(char *line, int len, double *values,
                             int lenValues);
    void messageAbort(const char* mesg) const;

    fei::SharedPtr<Aztec_BlockMap> amap_;

    AZ_MATRIX *Amat_;

    int N_update_;
    int* external_;
    int* extern_index_;
    int* update_index_;
    int* data_org_;
    int* orderingUpdate_;

    bool isLoaded_;
    bool isAllocated_;

    int localNNZ_;
    int* nnzPerRow_;

    int numRemoteBlocks_;
    int* remoteInds_;
    int* remoteBlockSizes_;
};

}//namespace fei_trilinos

#endif
