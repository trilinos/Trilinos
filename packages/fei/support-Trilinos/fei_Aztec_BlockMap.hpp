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

#ifndef _fei_Aztec_BlockMap_h_
#define _fei_Aztec_BlockMap_h_


//
// This Aztec_BlockMap class is a wrapper that encapsulates the general
// information needed to describe the layout of an Aztec DVBR matrix.
// It is a companion/support class that goes with the data class wrappers
// Aztec_LSVector and AztecDVBR_Matrix. Aztec_BlockMap inherits from
// Aztec_Map.
//
// Aztec_Map allows the storage and retrieval of information such as
// local and global sizes, the MPI communicator, and the proc_config array.
// Aztec_BlockMap describes the partitioning and layout of a block matrix.
//

namespace fei_trilinos {

class Aztec_BlockMap : public Aztec_Map {
    
  public:
    Aztec_BlockMap(int globalSize, int N_update, const int* update, int localOffset,
                   MPI_Comm comm,
                   int numGlobalBlocks, int numLocalBlocks, const int* blkUpdate,
                   int localBlockOffset, int* blockSizes);

    Aztec_BlockMap(const Aztec_BlockMap& map);       // copy constructor
    virtual ~Aztec_BlockMap(void);

    const int& getNumGlobalBlocks() const {return(numGlobalBlocks_);}
    const int& getNumLocalBlocks() const {return(numLocalBlocks_);}
    const int& getLocalBlockOffset() const {return(localBlockOffset_);}

    const int* getBlockSizes() const {return(blockSizes_);}
    int* getBlockUpdate() {return(blockUpdate_);}

  private:

    void checkInput();

    int numGlobalBlocks_;
    int numLocalBlocks_;
    int localBlockOffset_;
    int* blockSizes_;
    int* blockUpdate_;
};

}// namespace fei_trilinos

#endif

