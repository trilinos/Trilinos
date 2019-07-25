// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
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
