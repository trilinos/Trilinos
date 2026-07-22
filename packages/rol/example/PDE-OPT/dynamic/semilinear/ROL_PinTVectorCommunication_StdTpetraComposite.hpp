// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PINTVECTORCOMMUNICATION_STDTPETRACOMPOSITEVECTOR_H
#define ROL_PINTVECTORCOMMUNICATION_STDTPETRACOMPOSITEVECTOR_H

#include "ROL_PinTVectorCommunication.hpp"
#include "ROL_PinTVectorCommunication_StdVector.hpp"
#include "ROL_PinTVectorCommunication_Tpetra.hpp"

/** @ingroup la_group
    \class ROL::PinTVector
    \brief Defines a parallel in time vector.
*/

namespace ROL {

template <class Real> 
class PinTVectorCommunication_StdTpetraComposite : public PinTVectorCommunication<Real> {

  PinTVectorCommunication_StdVector<Real> stdComm_;
  PinTVectorCommunication_Tpetra<Real> tpetraComm_;

public:

  virtual ~PinTVectorCommunication_StdTpetraComposite() {}
 
  /**
   * \brief Send a vector to a neighboring processor.
   */
  void send(MPI_Comm comm,int rank,Vector<Real> & source,int tag=0) const override
  {
    if(isTpetra(source)) {
      tpetraComm_.send(comm,rank,source,tag);
      return;
    }

    if(isStd(source)) {
      stdComm_.send(comm,rank,source,tag);
      return;
    }

    assert(false);
  }

  /**
   * \brief Recieve a vector from a neighboring processor. Uses blocking communication. 
   */
  void recv(MPI_Comm comm,int rank,Vector<Real> & dest,int tag=0) const override
  {
    if(isTpetra(dest)) {
      tpetraComm_.recv(comm,rank,dest,tag);
      return;
    }

    if(isStd(dest)) {
      stdComm_.recv(comm,rank,dest,tag);
      return;
    }

    assert(false);
  }

  /**
   * \brief Recieve a vector from a neighboring processor. Uses blocking communication. 
   */
  void recvSumInto(MPI_Comm comm,int rank,Vector<Real> & dest,int tag=0) const override
  {
    if(isTpetra(dest)) {
      tpetraComm_.recvSumInto(comm,rank,dest,tag);
      return;
    }

    if(isStd(dest)) {
      stdComm_.recvSumInto(comm,rank,dest,tag);
      return;
    }

    assert(false);
  }

  bool isTpetra(const Vector<Real> & vec) const
  {
    auto tpetraVec = dynamicPtrCast<const TpetraMultiVector<Real>>(makePtrFromRef(vec));

    return tpetraVec!=ROL::nullPtr;
  }

  bool isStd(const Vector<Real> & vec) const
  {
    auto stdVec = dynamicPtrCast<const StdVector<Real>>(makePtrFromRef(vec));

    return stdVec!=ROL::nullPtr;
  }
};

} // namespace ROL

#endif
