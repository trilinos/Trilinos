/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

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


//==============================================================================
Aztec_BlockMap::Aztec_BlockMap(int globalSz, int localSz, int localOffs,
                          MPI_Comm comm,
                          int numGlobalBlocks, int numLocalBlocks,
                          int localBlockOffset, int* blockSizes)

  : Aztec_Map(globalSz, localSz, localOffs, comm),
    numGlobalBlocks_(numGlobalBlocks),
    numLocalBlocks_(numLocalBlocks),
    localBlockOffset_(localBlockOffset)
{
    checkInput();

    blockSizes_ = new int[numLocalBlocks_];

    for(int i=0; i<numLocalBlocks_; i++) {
        blockSizes_[i] = blockSizes[i];
    }
}
 
//==============================================================================
Aztec_BlockMap::Aztec_BlockMap(const Aztec_BlockMap& map)
  : Aztec_Map(map.globalSize(), map.localSize(), map.localOffset(),
              map.getCommunicator()),
    numGlobalBlocks_(map.numGlobalBlocks_),
    numLocalBlocks_(map.numLocalBlocks_),
    localBlockOffset_(map.localBlockOffset_)
{
    blockSizes_ = new int[numLocalBlocks_];

    for(int i=0; i<numLocalBlocks_; i++) {
        blockSizes_[i] = map.blockSizes_[i];
    }
}

//==============================================================================
Aztec_BlockMap::~Aztec_BlockMap(void)  {

    delete [] blockSizes_;
    blockSizes_ = NULL;

    numGlobalBlocks_ = 0;
    numLocalBlocks_ = 0;
    localBlockOffset_ = 0;
}

//==============================================================================
void Aztec_BlockMap::checkInput() {

    if (numGlobalBlocks_ <= 0) {
       FEI_CERR << "Aztec_BlockMap: ERROR, numGlobalBlocks <= 0. Aborting."
            << FEI_ENDL;
       abort();
    }

    if (localBlockOffset_ < 0) {
       FEI_CERR << "Aztec_BlockMap: ERROR, negative localBlockOffset. Aborting."
            << FEI_ENDL;
       abort();
    }
}

#endif
//HAVE_FEI_AZTECOO
