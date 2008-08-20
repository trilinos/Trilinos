/*--------------------------------------------------------------------*/
/*    Copyright 2007 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_CommCore_hpp_
#define _fei_CommCore_hpp_

#include "fei_macros.hpp"
#include "fei_mpi.h"

#include <vector>

namespace fei {

  /** core implementation class for basic communication stuff.
    Implementation detail, not for public consumption. */
  class CommCore {
  public:
    /** destructor */
    virtual ~CommCore();

    /** factory-like creation method */
    static CommCore* create();

    /** query for a tag */
    int get_tag();

    /** increment tag attribute */
    void increment_tag();

#ifndef FEI_SER
    /** array of MPI_Request objects */
    std::vector<MPI_Request> mpiReqs_;

    /** array of MPI_Status objects */
    std::vector<MPI_Status> mpiStatuses_;
#endif

    /** integer work-array */
    std::vector<int> tmpIntData_;

  private:
    CommCore();
    int tag_;
  };//class CommCore

} //namespace fei

#endif // _fei_CommCore_hpp_

