
// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PINTVECTORCOMMUNICATION_H
#define ROL_PINTVECTORCOMMUNICATION_H

#include <ostream>

#include "ROL_Vector.hpp"

namespace ROL {

/** 
 * A virtual base class for abstracting communication of full vectors. Useful for parallel-in-time.
 */ 
template <class Real> 
class PinTVectorCommunication {
public:
  virtual ~PinTVectorCommunication() {}
  /**
   * \brief Send a vector to a neighboring processor.
   */
  virtual void send(MPI_Comm comm,int rank,Vector<Real> & source,int tag=0) const = 0;

  /**
   * \brief Recieve a vector from a neighboring processor. Uses blocking communication. 
   */
  virtual void recv(MPI_Comm comm,int rank,Vector<Real> & dest,int tag=0) const = 0;

  /**
   * \brief Recieve a vector from a neighboring processor. Uses blocking communication. 
   */
  virtual void recvSumInto(MPI_Comm comm,int rank,Vector<Real> & dest,int tag=0) const = 0;

  /**
   * \brief Recieve a vector from a neighboring processor. Uses blocking communication. Optionally, sum into
   */
  void recv(MPI_Comm comm,int rank,Vector<Real> & dest,bool sumInto,int tag=0) const 
  {
    if(sumInto)
      recvSumInto(comm,rank,dest,tag);
    else
      recv(comm,rank,dest,tag);
  }
};

} // namespace ROL

#endif
