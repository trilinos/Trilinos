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


#include <fei_trilinos_macros.hpp>
#include <fei_iostream.hpp>

#ifdef HAVE_FEI_AZTECOO

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <fei_mpi.h>

#ifndef FEI_SER

#define AZTEC_MPI AZTEC_MPI
#define AZ_MPI AZ_MPI
#ifndef MPI
#define MPI MPI
#endif

#endif

#include <az_aztec.h>
#include <fei_Aztec_Map.hpp>
#include <fei_Aztec_BlockMap.hpp>

namespace fei_trilinos {

//==============================================================================
Aztec_BlockMap::Aztec_BlockMap(int globalSz, int N_update, const int* update, int localOffs,
                          MPI_Comm comm,
                          int numGlobalBlocks, int numLocalBlocks,
                          const int* blkUpdate,
                          int localBlockOffset, int* blockSizes)

  : Aztec_Map(globalSz, N_update, update, localOffs, comm),
    numGlobalBlocks_(numGlobalBlocks),
    numLocalBlocks_(numLocalBlocks),
    localBlockOffset_(localBlockOffset)
{
    checkInput();

    blockSizes_ = new int[numLocalBlocks_];
    blockUpdate_ = new int[numLocalBlocks_];

    for(int i=0; i<numLocalBlocks_; i++) {
        blockSizes_[i] = blockSizes[i];
        blockUpdate_[i] = blkUpdate[i];
    }
}
 
//==============================================================================
Aztec_BlockMap::Aztec_BlockMap(const Aztec_BlockMap& map)
  : Aztec_Map(map.globalSize(), map.localSize(), &map.update[0], map.localOffset(),
              map.getCommunicator()),
    numGlobalBlocks_(map.numGlobalBlocks_),
    numLocalBlocks_(map.numLocalBlocks_),
    localBlockOffset_(map.localBlockOffset_)
{
    blockSizes_ = new int[numLocalBlocks_];
    blockUpdate_ = new int[numLocalBlocks_];

    for(int i=0; i<numLocalBlocks_; i++) {
        blockSizes_[i] = map.blockSizes_[i];
        blockUpdate_[i] = map.blockUpdate_[i];
    }
}

//==============================================================================
Aztec_BlockMap::~Aztec_BlockMap(void)  {

    delete [] blockSizes_;
    blockSizes_ = NULL;

    delete [] blockUpdate_;
    blockUpdate_ = NULL;

    numGlobalBlocks_ = 0;
    numLocalBlocks_ = 0;
    localBlockOffset_ = 0;
}

//==============================================================================
void Aztec_BlockMap::checkInput() {

    if (numGlobalBlocks_ <= 0) {
       fei::console_out() << "Aztec_BlockMap: ERROR, numGlobalBlocks <= 0. Aborting."
            << FEI_ENDL;
       abort();
    }

    if (localBlockOffset_ < 0) {
       fei::console_out() << "Aztec_BlockMap: ERROR, negative localBlockOffset. Aborting."
            << FEI_ENDL;
       abort();
    }
}

}//namespace fei_trilinos

#endif
//HAVE_FEI_AZTECOO
