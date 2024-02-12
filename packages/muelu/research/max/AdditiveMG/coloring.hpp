#include <math.h>

int coloring2D(int mypid, int ndx) {
  // INPUT: mypid is the MPI id of the processor
  // INPUT: ndx is the nuber of domains along x-direction
  // INPUT: ndy is the nuber of domains along y-direction
  // OUTPUT: color label associated with the current MPI processor (=subdomain)

  // detect the row of the domain grid where the current subdomain is located
  int grid_row = std::ceil(static_cast<double>(mypid) / ndx);

  // detect the y-coordinate in the domain grid associated with the current subdomain
  int ypos = grid_row % 3;
  int xpos;
  int color = -1;

  // detect the x-coordinate in the domain grid associated with current subdomain
  if (0 != mypid % ndx)
    xpos = static_cast<int>(mypid - std::floor(static_cast<double>(mypid) / ndx) * ndx) % 3;
  else
    xpos = (mypid - ((ypos - 1) * ndx)) % 3;

  // use x and y coordinates to determine the color of the current subdomain
  if (xpos > 0 && ypos > 0)
    color = (ypos - 1) * 3 + xpos;
  else if (xpos > 0 && ypos == 0)
    color = 6 + xpos;
  else if (xpos == 0 && ypos > 0)
    color = ypos * 3;
  else
    color = 9;

  // TEUCHOS_TEST_FOR_EXCEPT( color < 1 );
  return color - 1;
}

int coloring3D(int mypid, int ndx, int ndy) {
  // INPUT: mypid is the MPI id of the processor
  // INPUT: ndx is the nuber of domains along x-direction
  // INPUT: ndy is the nuber of domains along y-direction
  // OUTPUT: color label associated with the current MPI processor (=subdomain)

  // detect the plane of the domain grid where the current subdomain resides
  int grid_plane = std::ceil(static_cast<double>(mypid) / (ndx * ndy));

  // On the given plane, find the local id of the current subdomain
  int plane_id = mypid % (ndx * ndy);

  // detect the row on the current two-dimensional grid where the current subdomain is located
  int plane_row = std::ceil(static_cast<double>(plane_id) / ndx);

  // detect the y-coordinate on the given plane associated with the current subdomain
  int ypos = plane_row % 3;
  int xpos;
  int label_plane;
  int color = -1;

  // detect the x-coordinate on the given plane associated with the current subdomain
  if (0 != plane_id % ndx)
    xpos = static_cast<int>(mypid - std::floor(static_cast<double>(plane_id) / ndx) * ndx) % 3;
  else
    xpos = (plane_id - ((ypos - 1) * ndx)) % 3;

  // find a two-dimensional coloring to give to the current subdomain according to its position
  // in the plane where it resides
  if (xpos > 0 && ypos > 0)
    label_plane = (ypos - 1) * 3 + xpos;
  else if (xpos > 0 && ypos == 0)
    label_plane = 6 + xpos;
  else if (xpos == 0 && ypos > 0)
    label_plane = ypos * 3;
  else
    label_plane = 9;

  // use the two-dimensional coloring to determine the three-dimensional coloring
  if (grid_plane % 3 == 1)
    color = label_plane;

  else if (grid_plane % 3 == 2)
    color = 9 + label_plane;

  else if (grid_plane % 3 == 0)
    color = 18 + label_plane;

  // TEUCHOS_TEST_FOR_EXCEPT( color < 1 );
  return color - 1;
}
