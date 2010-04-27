/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

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

namespace fei_trilinos {

//==============================================================================
Aztec_Map::Aztec_Map(int globalSz, int N_update, const int* update,
                     int localOffs, MPI_Comm comm)
 :
    proc_config(AZ_PROC_SIZE),
    update(update, update+N_update),
    external(NULL),
    update_index(NULL),
    extern_index(NULL),
    data_org(NULL),
    orderingUpdate(),
    az_transformed(false),
    globalSize_(globalSz),
    localSize_(N_update),
    localOffset_(localOffs),
    N_update_(N_update),
    comm_(comm)
{
    checkInput();
    AZ_set_proc_config(&proc_config[0], comm_);
}
 
//==============================================================================
Aztec_Map::Aztec_Map(const Aztec_Map& map) :
    proc_config(AZ_PROC_SIZE),
    update(map.update),
    external(NULL),
    update_index(NULL),
    extern_index(NULL),
    data_org(NULL),
    orderingUpdate(map.orderingUpdate),
    az_transformed(map.az_transformed),
    globalSize_(map.globalSize_),
    localSize_(map.localSize_),
    localOffset_(map.localOffset_),
    N_update_(map.localSize_),
    comm_(map.comm_)
{
    AZ_processor_info(&proc_config[0]);
    update_index = (int*)AZ_allocate(N_update_*sizeof(int));
    for(int i=0; i<N_update_; ++i) {
      update_index[i] = map.update_index[i];
    }
    external = (int*)AZ_allocate(data_org[AZ_N_external]*sizeof(int));
    extern_index = (int*)AZ_allocate(data_org[AZ_N_external]*sizeof(int));
    for(int i=0; i<data_org[AZ_N_external]; ++i) {
      external[i] = map.external[i];
      extern_index[i] = map.extern_index[i];
    }
}

//==============================================================================
Aztec_Map::~Aztec_Map()
{
  globalSize_ = 0;
  localSize_ = 0;
  localOffset_ = 0;

  std::free(update_index);
  std::free(external);
  std::free(extern_index);
  std::free(data_org);
}

//==============================================================================
void Aztec_Map::checkInput() {

    if (globalSize_ <= 0) {
        throw std::runtime_error("Aztec_Map: ERROR, globalSize <= 0.");
    }

    if (localSize_ < 0) {
        throw std::runtime_error("Aztec_Map: ERROR, localSize negative.");
    }

    if (localOffset_ < 0) {
        throw std::runtime_error("Aztec_Map: ERROR, localOffset negative.");
    }
}

}//namespace fei_trilinos

#endif
//HAVE_FEI_AZTECOO

