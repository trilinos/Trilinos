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

#ifndef OED_SPLITCOMM_HPP
#define OED_SPLITCOMM_HPP

#ifdef HAVE_MPI

#include "ROL_Ptr.hpp"
#include <mpi.h>

template<class Ordinal, class Comm>
class OED_SplitComm {
private:
  ROL::Ptr<Comm> comm_design_;
  ROL::Ptr<Comm> comm_sample_;

public:
  OED_SplitComm(int m) {
    Ordinal rank, Ngroups, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    Ordinal M = (m > size ? 1 : m);
    Ngroups = size/M;
    // Instantiate Design Communicator.
    MPI_Comm design_comm;
    MPI_Comm_split(MPI_COMM_WORLD, rank/M, rank, &design_comm);
    Ordinal comRank; // Process rank in linear algebra communicator.
    Ordinal comSize; // Number of processes in linear algebra communicator.
    MPI_Comm_rank(design_comm,&comRank); // Get process rank.
    MPI_Comm_size(design_comm,&comSize); // Get communicator size.
    comm_design_ = ROL::makePtr<Comm>(design_comm); // Wrap as Comm.

    // Determine group ranks for sampling.
    std::vector<Ordinal> granks(Ngroups);
    for (Ordinal i=0;i<Ngroups;i++) granks[i] = comRank + i*M;

    // Build MPI groups for sampling.
    MPI_Group world_comm; // Grab MPI_COMM_WORLD and place in world_comm.
    MPI_Comm_group(MPI_COMM_WORLD,&world_comm);
    MPI_Group group;
    MPI_Group_incl(world_comm,Ngroups,&granks[0],&group);

    // Instantiate Sample Communicator based on group.
    MPI_Comm sample_comm;
    MPI_Comm_create(MPI_COMM_WORLD, group, &sample_comm);
    Ordinal comRank1; // Process rank in sample communicator.
    Ordinal comSize1; // Number of processes in sample communicator.
    MPI_Comm_rank(sample_comm,&comRank1); // Get process rank.
    MPI_Comm_size(sample_comm,&comSize1); // Get communicator size.
    comm_sample_ = ROL::makePtr<Comm>(sample_comm); // Wrap as Comm.
  }

  const ROL::Ptr<Comm> getDesignComm() const {
    return comm_design_;
  }

  const ROL::Ptr<Comm> getSampleComm() const {
    return comm_sample_;
  }
};

#endif // HAVE_MPI
#endif // OED_SPLITCOMM_HPP
