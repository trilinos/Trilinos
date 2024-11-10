// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

void createBrickMap1D(int numGlobalElements, std::vector<int>& ind, Teuchos::RCP<const Teuchos::Comm<int> > comm) {
  // INPUT: numGlobalElements = number of finite difference nodes along x
  // INPUT: ind = vector containing indices of the rows owned by the current MPI processor
  // INPUT: comm = MPI communicator (MPI_COMM_WORLD)

  int mypid = comm->getRank();

  ind.reserve(static_cast<int>(numGlobalElements / (comm->getSize() - 1) + 1));
  if (mypid != 0 && mypid != comm->getSize() - 1)
    for (int i = 0; i <= (static_cast<int>(numGlobalElements / (comm->getSize() - 1))) - 1; ++i)
      ind.emplace_back((mypid - 1) * static_cast<int>(numGlobalElements / (comm->getSize() - 1)) + i);

  if (mypid == comm->getSize() - 1)
    for (int i = (mypid - 1) * static_cast<int>(numGlobalElements / (comm->getSize() - 1)); i != numGlobalElements; ++i)
      ind.emplace_back(i);
}

void createBrickMap2D(int nx, int brick_sizex, int brick_sizey, std::vector<int>& ind, Teuchos::RCP<const Teuchos::Comm<int> > comm) {
  // INPUT: nx = number of finite difference nodes along x
  // INPUT: brick_sizex = size of a brick aggregate along x-direction
  // INPUT: brick_sizey = size of a brick aggregate along y-direction
  // INPUT: ind = vector containing indices of the rows owned by the current MPI processor
  // INPUT: comm = MPI communicator (MPI_COMM_WORLD)

  int mypid = comm->getRank();

  ind.reserve(brick_sizex * brick_sizey);

  int ndx = nx / brick_sizex;

  if (mypid != 0) {
    int grid_row = std::ceil(static_cast<double>(mypid) / ndx);
    int ypos     = grid_row;
    int xpos;

    if (0 != mypid % ndx)
      xpos = mypid % ndx;
    else
      xpos = ndx;

    int preliminary = nx * brick_sizey * (ypos - 1);

    for (int row = 0; row < brick_sizey; ++row)
      for (int col = brick_sizex * (xpos - 1) + 1; col <= brick_sizex * xpos; ++col)
        ind.emplace_back(preliminary + row * nx + col - 1);
  }
}

void createBrickMap3D(int nx, int ny, int brick_sizex, int brick_sizey, int brick_sizez, std::vector<int>& ind, Teuchos::RCP<const Teuchos::Comm<int> > comm) {
  // INPUT: nx = number of finite difference nodes along x
  // INPUT: ny = number of finite difference nodes along y
  // INPUT: brick_sizex = size of a brick aggregate along x-direction
  // INPUT: brick_sizey = size of a brick aggregate along y-direction
  // INPUT: ind = vector containing indices of the rows owned by the current MPI processor
  // INPUT: comm = MPI communicator (MPI_COMM_WORLD)

  int mypid = comm->getRank();

  ind.reserve(brick_sizex * brick_sizey * brick_sizez);

  // determine the number of subdomains along x and y directions
  int ndx = nx / brick_sizex;
  int ndy = ny / brick_sizey;

  if (mypid != 0) {
    int grid_plane = std::ceil(static_cast<double>(mypid) / (ndx * ndy));
    int plane_id   = mypid % (ndx * ndy);

    if (0 == plane_id)
      plane_id = ndx * ndy;

    int plane_row = std::ceil(static_cast<double>(plane_id) / ndx);
    int ypos      = plane_row;
    int xpos;

    if (0 != plane_id % ndx)
      xpos = plane_id % ndx;
    else
      xpos = ndx;

    int preliminary = nx * ny * brick_sizez * (grid_plane - 1) + nx * brick_sizey * (ypos - 1);

    for (int l = 0; l < brick_sizez; ++l)
      for (int row = 0; row < brick_sizey; ++row)
        for (int col = brick_sizex * (xpos - 1) + 1; col <= brick_sizex * xpos; ++col)
          ind.emplace_back(preliminary + l * nx * ny + row * nx + col - 1);
  }
}
