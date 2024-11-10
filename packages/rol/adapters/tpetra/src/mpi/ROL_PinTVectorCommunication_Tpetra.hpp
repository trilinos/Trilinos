// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PINTVECTORCOMMUNICATION_TPETRA_H
#define ROL_PINTVECTORCOMMUNICATION_TPETRA_H

#include "mpi.h"

#include "ROL_PinTVectorCommunication.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_TpetraMultiVector.hpp"
#include "Tpetra_Access.hpp"

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
    auto view = tp_source.getLocalViewHost(Tpetra::Access::ReadOnly);

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
    auto view = tp_dest.getLocalViewHost(Tpetra::Access::ReadWrite);

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
    auto view = tp_dest.getLocalViewHost(Tpetra::Access::ReadWrite);

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
