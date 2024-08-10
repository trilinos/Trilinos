// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef OED_SPLITCOMM_HPP
#define OED_SPLITCOMM_HPP

#ifdef HAVE_MPI

#include "ROL_Ptr.hpp"
#include <mpi.h>

namespace OED {

template<class Ordinal, class Comm>
class SplitComm {
private:
  ROL::Ptr<Comm> comm_design_;
  ROL::Ptr<Comm> comm_sample_;

public:
  SplitComm(const int m, Mpi_Comm global_comm = MPI_COMM_WORLD) { // CHECK: ALLOW MPI_COMM_WORLD
    Ordinal rank, Ngroups, size;
    MPI_Comm_rank(global_comm, &rank);
    MPI_Comm_size(global_comm, &size);
    Ordinal M = (m > size ? 1 : m);
    Ngroups = size/M;
    // Instantiate Design Communicator.
    MPI_Comm design_comm;
    MPI_Comm_split(global_comm, rank/M, rank, &design_comm);
    Ordinal comRank; // Process rank in linear algebra communicator.
    Ordinal comSize; // Number of processes in linear algebra communicator.
    MPI_Comm_rank(design_comm,&comRank); // Get process rank.
    MPI_Comm_size(design_comm,&comSize); // Get communicator size.
    comm_design_ = ROL::makePtr<Comm>(design_comm); // Wrap as Comm.

    // Determine group ranks for sampling.
    std::vector<Ordinal> granks(Ngroups);
    for (Ordinal i=0;i<Ngroups;i++) granks[i] = comRank + i*M;

    // Build MPI groups for sampling.
    MPI_Group world_comm; // Grab global_comm and place in world_comm.
    MPI_Comm_group(global_comm,&world_comm);
    MPI_Group group;
    MPI_Group_incl(world_comm,Ngroups,&granks[0],&group);

    // Instantiate Sample Communicator based on group.
    MPI_Comm sample_comm;
    MPI_Comm_create(global_comm, group, &sample_comm);
    Ordinal comRank1; // Process rank in sample communicator.
    Ordinal comSize1; // Number of processes in sample communicator.
    MPI_Comm_rank(sample_comm,&comRank1); // Get process rank.
    MPI_Comm_size(sample_comm,&comSize1); // Get communicator size.
    comm_sample_ = ROL::makePtr<Comm>(sample_comm); // Wrap as Comm.
  }

  const ROL::Ptr<Comm> getDesignComm(void) const {
    return comm_design_;
  }

  const ROL::Ptr<Comm> getSampleComm(void) const {
    return comm_sample_;
  }
};

} // End OED Namespace

#endif // HAVE_MPI
#endif // OED_SPLITCOMM_HPP
