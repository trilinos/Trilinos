/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <assert.h>
#include <fei_iostream.hpp>
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


//==============================================================================
Aztec_Map::Aztec_Map(int globalSz, int localSz,
                     int localOffs, MPI_Comm comm)
  : globalSize_(globalSz),
    localSize_(localSz),
    localOffset_(localOffs),
    comm_(comm)
{
    checkInput();

    proc_config_ = new int[AZ_PROC_SIZE];
    AZ_set_proc_config(proc_config_, comm_);
}
 
//==============================================================================
Aztec_Map::Aztec_Map(const Aztec_Map& map) :
    globalSize_(map.globalSize_),
    localSize_(map.localSize_),
    localOffset_(map.localOffset_),
    comm_(map.comm_)
{
    proc_config_ = new int[AZ_PROC_SIZE];
    AZ_processor_info(proc_config_);
}

//==============================================================================
Aztec_Map::~Aztec_Map(void)  {

    globalSize_ = 0;
    localSize_ = 0;
    localOffset_ = 0;

    delete [] proc_config_;
    proc_config_ = NULL;
}

//==============================================================================
void Aztec_Map::checkInput() {

    if (globalSize_ <= 0) {
        FEI_CERR << "Aztec_Map: ERROR, globalSize <= 0. Aborting." << FEI_ENDL;
        abort();
    }

    if (localSize_ < 0) {
        FEI_CERR << "Aztec_Map: ERROR, localSize negative. Aborting." << FEI_ENDL;
        abort();
    }

    if (localOffset_ < 0) {
        FEI_CERR << "Aztec_Map: ERROR, localOffset negative. Aborting." << FEI_ENDL;
        abort();
    }
}

