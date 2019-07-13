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
