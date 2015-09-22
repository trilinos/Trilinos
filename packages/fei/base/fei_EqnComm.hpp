/*--------------------------------------------------------------------*/
/*    Copyright 2007 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_EqnComm_hpp_
#define _fei_EqnComm_hpp_

#include <fei_macros.hpp>
#include <fei_fwd.hpp>
#include <fei_mpi.h>

namespace fei {
class EqnComm {
 public:
  /** constructor */
  EqnComm(MPI_Comm comm, int numLocalEqns);
  EqnComm(MPI_Comm comm, int numLocalEqns, const std::vector<int>& globalOffsets);

  /** destructor */
  virtual ~EqnComm();

  const std::vector<int>& getGlobalOffsets() const;

  int getOwnerProc(int eqn) const;

 private:
  MPI_Comm comm_;
  std::vector<int> globalOffsets_;
};//class EqnComm
}//namespace fei
#endif

