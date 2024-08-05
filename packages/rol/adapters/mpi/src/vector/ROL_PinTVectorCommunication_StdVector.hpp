// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PINTVECTORCOMMUNICATION_STDVECTOR_H
#define ROL_PINTVECTORCOMMUNICATION_STDVECTOR_H

#include "mpi.h"

#include "ROL_PinTVectorCommunication.hpp"
#include "ROL_StdVector.hpp"

/** @ingroup la_group
    \class ROL::PinTVector
    \brief Defines a parallel in time vector.
*/

namespace ROL {

template <class Real> 
class PinTVectorCommunication_StdVector : public PinTVectorCommunication<Real> {
public:

  virtual ~PinTVectorCommunication_StdVector() {}
 
  /**
   * \brief Send a vector to a neighboring processor.
   */
  void send(MPI_Comm comm,int rank,Vector<Real> & source,int tag=0) const override
  {
    const std::vector<Real> & std_source = *dynamic_cast<const StdVector<Real>&>(source).getVector();

    // int myRank = -1;
    // MPI_Comm_rank(comm, &myRank);

    MPI_Send(const_cast<Real*>(&std_source[0]),int(std_source.size()),MPI_DOUBLE,rank,tag,comm);
  }

  /**
   * \brief Recieve a vector from a neighboring processor. Uses blocking communication. 
   */
  void recv(MPI_Comm comm,int rank,Vector<Real> & dest,int tag=0) const override
  {
    std::vector<Real> & std_dest = *dynamic_cast<StdVector<Real>&>(dest).getVector();

    // int myRank = -1;
    // MPI_Comm_rank(comm, &myRank);

    MPI_Recv(&std_dest[0],int(std_dest.size()),MPI_DOUBLE,rank,tag,comm,MPI_STATUS_IGNORE);
  }

  /**
   * \brief Recieve a vector from a neighboring processor. Uses blocking communication. 
   */
  void recvSumInto(MPI_Comm comm,int rank,Vector<Real> & dest,int tag=0) const override
  {
    std::vector<Real> & std_dest = *dynamic_cast<StdVector<Real>&>(dest).getVector();
    std::vector<Real> buffer(std_dest.size(),0.0);

    // int myRank = -1;
    // MPI_Comm_rank(comm, &myRank);

    MPI_Recv(&buffer[0],int(buffer.size()),MPI_DOUBLE,rank,tag,comm,MPI_STATUS_IGNORE);

    for(size_t i=0;i<std_dest.size();i++)
      std_dest[i] += buffer[i];
  }
};

} // namespace ROL

#endif
