#ifndef _fei_set_shared_ids_hpp_
#define _fei_set_shared_ids_hpp_

/*--------------------------------------------------------------------*/
/*    Copyright 2009 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>
#include <fei_mpi.h>
#include <fei_fwd.hpp>

namespace fei {

/** Given a record-collection, perform inter-processor communication
  to find out which IDs are shared among multiple processors, and
  fill the SharedIDs object with those mappings.
*/
void set_shared_ids(MPI_Comm comm,
                    const snl_fei::RecordCollection& records,
                    fei::SharedIDs<int>& sharedIDs,
                    int lowest_global_id, int highest_global_id);

}

#endif

