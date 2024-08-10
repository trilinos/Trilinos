// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

void split_comm_world(ROL::Ptr<Teuchos::Comm<int> > & Comm_linalg,
                      ROL::Ptr<Teuchos::Comm<int> > & Comm_sample,
                      int M) {
  int rank, Ngroups, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (M>size) {
    M = 1;
  }
  Ngroups = size/M;
  // Instantiate Linear Algebra Communicator.
  MPI_Comm linalg_comm;
  MPI_Comm_split(MPI_COMM_WORLD, rank/M, rank, &linalg_comm);
  int comRank; // Process rank in linear algebra communicator.
  int comSize; // Number of processes in linear algebra communicator.
  MPI_Comm_rank(linalg_comm,&comRank);     // Get process rank.
  MPI_Comm_size(linalg_comm,&comSize);     // Get communicator size.
  Comm_linalg = ROL::makePtr<Teuchos::MpiComm<int>>(linalg_comm); // Wrap as Teuchos::MpiComm.

  // Determine group ranks for sampling.
  Teuchos::Array<int> granks(Ngroups);
  for (int i=0;i<Ngroups;i++) granks[i] = comRank + i*M;

  // Build MPI groups for sampling.
  MPI_Group world_comm; // Grab MPI_COMM_WORLD and place in world_comm.
  MPI_Comm_group(MPI_COMM_WORLD,&world_comm);
  MPI_Group group;
  MPI_Group_incl(world_comm,Ngroups,&granks[0],&group);

  // Instantiate Sample Communicator based on group.
  MPI_Comm sample_comm;
  MPI_Comm_create(MPI_COMM_WORLD, group, &sample_comm);
  int comRank1; // Process rank in sample communicator.
  int comSize1; // Number of processes in sample communicator.
  MPI_Comm_rank(sample_comm,&comRank1);    // Get process rank.
  MPI_Comm_size(sample_comm,&comSize1);    // Get communicator size.
  Comm_sample = ROL::makePtr<Teuchos::MpiComm<int>>(sample_comm); // Wrap as Teuchos::MpiComm.
}

