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

#ifndef ROL_PINTVECTORCOMMUNICATION_TPETRA_H
#define ROL_PINTVECTORCOMMUNICATION_TPETRA_H

#include "mpi.h"

#include "ROL_PinTVectorCommunication.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_TpetraMultiVector.hpp"

/** @ingroup la_group
    \class ROL::PinTVector
    \brief Defines a parallel in time vector.
*/

namespace ROL {

template <class Real> 
class PinTVectorCommunication_Tpetra : public PinTVectorCommunication<Real> {
public:
  /**
   * \brief Send a vector to a neighboring processor.
   */
  void send(MPI_Comm comm,int rank,Vector<Real> & source,int tag=0) const override
  {
    auto & tp_source = *dynamic_cast<TpetraMultiVector<Real>&>(source).getVector();
    tp_source.sync_host();

    auto view = tp_source.getLocalViewHost();

    // int myRank = -1;
    // MPI_Comm_rank(comm, &myRank);

    MPI_Send(const_cast<Real*>(&view(0,0)),int(view.extent(0)*view.extent(1)),MPI_DOUBLE,rank,tag,comm);
  }

  /**
   * \brief Recieve a vector from a neighboring processor. Uses blocking communication. 
   */
  void recv(MPI_Comm comm,int rank,Vector<Real> & dest,int tag=0) const override
  {
    auto & tp_dest = *dynamic_cast<TpetraMultiVector<Real>&>(dest).getVector();
    auto view = tp_dest.getLocalViewHost();

    // int myRank = -1;
    // MPI_Comm_rank(comm, &myRank);

    MPI_Recv(&view(0,0),int(view.extent(0)*view.extent(1)),MPI_DOUBLE,rank,tag,comm,MPI_STATUS_IGNORE);

    // tp_source.template sync<Kokkos::DeviceSpace>();
  }

  /**
   * \brief Recieve a vector from a neighboring processor. Uses blocking communication. 
   */
  void recvSumInto(MPI_Comm comm,int rank,Vector<Real> & dest,int tag=0) const override
  {
    auto & tp_dest = *dynamic_cast<TpetraMultiVector<Real>&>(dest).getVector();
    auto view = tp_dest.getLocalViewHost();

    int myRank = -1;
    MPI_Comm_rank(comm, &myRank);

    assert(view.extent(1)==1);

    std::vector<Real> buffer(view.extent(0)*view.extent(1),0.0);
    MPI_Recv(&buffer[0],int(buffer.size()),MPI_DOUBLE,rank,tag,comm,MPI_STATUS_IGNORE);

    for(size_t i=0;i<buffer.size();i++)
      view(i,0) += buffer[i];
  }
};

} // namespace ROL

#endif
