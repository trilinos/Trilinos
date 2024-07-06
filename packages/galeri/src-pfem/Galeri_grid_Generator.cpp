// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_RefCountPtr.hpp"

#include "Galeri_grid_Generator.h"
#include "Galeri_grid_Triangle.h"
#include "Galeri_grid_Hex.h"

// ===========================================================================
void Galeri::grid::Generator::
getSquareWithTriangles(Epetra_Comm& comm, 
                       const int numGlobalElementsX, const int numGlobalElementsY,
                       const int numDomainsX, const int numDomainsY,
                       Galeri::grid::Loadable& domain, Galeri::grid::Loadable& boundary)
{
  getSquare(comm, numGlobalElementsX, numGlobalElementsY,
            numDomainsX, numDomainsY, domain, boundary, "Triangle");
}

// ===========================================================================
void Galeri::grid::Generator::
getSquareWithQuads(Epetra_Comm& comm, 
                   const int numGlobalElementsX, const int numGlobalElementsY,
                   const int numDomainsX, const int numDomainsY,
                   Galeri::grid::Loadable& domain, Galeri::grid::Loadable& boundary)
{
  getSquare(comm, numGlobalElementsX, numGlobalElementsY,
            numDomainsX, numDomainsY, domain, boundary, "Quad");
}

// ===========================================================================
void Galeri::grid::Generator::
getSquare(Epetra_Comm& comm, 
          const int numGlobalElementsX, const int numGlobalElementsY,
          const int numDomainsX, const int numDomainsY,
          Galeri::grid::Loadable& domain, Galeri::grid::Loadable& boundary,
          const std::string what)
{
  TEUCHOS_TEST_FOR_EXCEPTION(numDomainsX * numDomainsY != comm.NumProc(), std::logic_error,
                     "the number of processor should equal numDomainsX * numDomainsY"
                     << ", now numProcs = " << comm.NumProc()
                     << " and numDomainsX * numDomainsY = " << numDomainsX * numDomainsY);

  TEUCHOS_TEST_FOR_EXCEPTION(numGlobalElementsX % numDomainsX != 0, std::logic_error,
                     "numGlobalElementsX must be a multiple of numDomainsX");

  TEUCHOS_TEST_FOR_EXCEPTION(numGlobalElementsY % numDomainsY != 0, std::logic_error,
                     "numGlobalElementsY must be a multiple of numDomainsY");

  double lx = 1.0;
  double ly = 1.0;

  // these are the global number of elements and vertices
  int numGlobalElements = numGlobalElementsX * numGlobalElementsY;
  if (what == "Triangle") numGlobalElements *= 2;
  int numGlobalVertices = (numGlobalElementsX + 1) * (numGlobalElementsY + 1);

  int numGlobalVerticesX = numGlobalElementsX + 1;
  int numGlobalVerticesY = numGlobalElementsY + 1;

  // these are the mesh sizes, hx and hy
  double deltax = lx / numGlobalElementsX;
  double deltay = ly / numGlobalElementsY;

  // (px, py) are the coordinates of this processor.
  int px = comm.MyPID() % numDomainsX;
  int py = comm.MyPID() / numDomainsX;

  // (numMyElementsX, numMyElementsY) are the number of elements
  // in the square assigned to this processor, and
  // (numMyVerticesX, numMyVerticesY) the number of vertices.
  int numMyElementsX = numGlobalElementsX / numDomainsX;
  int numMyElementsY = numGlobalElementsY / numDomainsY;

  int numMyVerticesX = numMyElementsX + 1;
  int numMyVerticesY = numMyElementsY + 1;

  // (sx, sy) are the coordinates of the first element of this processor.
  int sx = px * numMyElementsX;
  int sy = py * numMyElementsY;

  // and these are the number of vertices and elements assigned
  // to this processor.
  int numMyElements = numMyElementsX * numMyElementsY;
  if (what == "Triangle") numMyElements *= 2;
  int numMyVertices = (numMyElementsX + 1) * (numMyElementsY + 1);

  Triangle triangle;

  domain.initialize(comm, numGlobalElements, numMyElements, triangle);

  int elementOffset = numMyElements * comm.MyPID();
  int vertexOffset  = px * numMyElementsX + py * numMyElementsY * numGlobalVerticesX;

  int count = 0;
  if (what == "Triangle")
  {
    for (int iy = 0; iy < numMyElementsY; ++iy)
    {
      for (int ix = 0; ix < numMyElementsX; ++ix)
      {
        int GEID = elementOffset + count++;
        int GVID = vertexOffset + ix + iy * numGlobalVerticesX;

        domain.setGlobalConnectivity(GEID, 0, GVID);
        domain.setGlobalConnectivity(GEID, 1, GVID + 1);
        domain.setGlobalConnectivity(GEID, 2, GVID + 2 + numGlobalElementsX);

        GEID = elementOffset + count++;
        domain.setGlobalConnectivity(GEID, 0, GVID + 2 + numGlobalElementsX);
        domain.setGlobalConnectivity(GEID, 1, GVID + 1 + numGlobalElementsX);
        domain.setGlobalConnectivity(GEID, 2, GVID);
      }
    }
  }
  else
  {
    for (int iy = 0; iy < numMyElementsY; ++iy)
    {
      for (int ix = 0; ix < numMyElementsX; ++ix)
      {
        int GEID = elementOffset + count++;
        int GVID = vertexOffset + ix + iy * numGlobalVerticesX;

        domain.setGlobalConnectivity(GEID, 0, GVID);
        domain.setGlobalConnectivity(GEID, 1, GVID + 1);
        domain.setGlobalConnectivity(GEID, 2, GVID + 2 + numGlobalElementsX);
        domain.setGlobalConnectivity(GEID, 3, GVID + 1 + numGlobalElementsX);
      }
    }
  }

  domain.freezeConnectivity();

  for (int iy = 0; iy < numMyVerticesY; ++iy)
  {
    for (int ix = 0; ix < numMyVerticesX; ++ix)
    {
      int GVID = vertexOffset + ix + iy * numGlobalVerticesX;

      domain.setGlobalCoordinates(GVID, 0, (sx + ix) * deltax);
      domain.setGlobalCoordinates(GVID, 1, (sy + iy) * deltay);
    }
  }

  domain.freezeCoordinates();

  // now build boundary faces

  int numMyBoundaries = 0;

  if (py == 0)               numMyBoundaries += numMyElementsX;
  if (py == numDomainsY - 1) numMyBoundaries += numMyElementsX;

  if (px == 0)               numMyBoundaries += numMyElementsY;
  if (px == numDomainsX - 1) numMyBoundaries += numMyElementsY;

  int pos = 0;
  vector<int> list(numMyBoundaries);

  if (py == 0)
  {
    int offset = px * numMyElementsX;

    for (int i = 0; i < numMyElementsX; ++i)
      list[pos++] = offset + i;
  }

  if (px == numDomainsX - 1)
  {
    int offset = numGlobalElementsX + py * numMyElementsY;

    for (int i = 0; i < numMyElementsY; ++i)
      list[pos++] = offset + i;
  }

  if (py == numDomainsY - 1)
  {
    int offset = numGlobalElementsX + numGlobalElementsY + px * numMyElementsX;

    for (int i = 0; i < numMyElementsX; ++i)
      list[pos++] = offset + i;
  }

  if (px == 0)
  {
    int offset = 2 * numGlobalElementsX + numGlobalElementsY + py * numMyElementsY;

    for (int i = 0; i < numMyElementsY; ++i)
      list[pos++] = offset + i;
  }

  TEUCHOS_TEST_FOR_EXCEPTION(pos != numMyBoundaries, std::logic_error,
                     "internal error in boundary list definition, " 
                     << pos << " vs. " << numMyBoundaries);

  Segment segment;
  boundary.initialize(comm, -1, numMyBoundaries, segment, &list[0]);

  // now insert the actual vertices in the grid

  if (py == 0)
  {
    int offset = px * numMyElementsX;

    for (int i = 0; i < numMyElementsX; ++i)
    {
      boundary.setGlobalConnectivity(offset + i, 0, offset + i);
      boundary.setGlobalConnectivity(offset + i, 1, offset + i + 1);
    }
  }

  if (px == numDomainsX - 1)
  {
    int offset = numGlobalVerticesX * py * numMyElementsY + numGlobalElementsX;
    int offset2 = numGlobalElementsX + py * numMyElementsY;

    for (int i = 0; i < numMyElementsY; ++i)
    {
      boundary.setGlobalConnectivity(offset2 + i, 0, offset + i * numGlobalVerticesX);
      boundary.setGlobalConnectivity(offset2 + i, 1, offset + (i + 1) * numGlobalVerticesX);
    }
  }

  if (py == numDomainsY - 1)
  {
    int offset = numGlobalVerticesX * numGlobalElementsY + px * numMyElementsX;
    int offset2 = numGlobalElementsX + numGlobalElementsY + px * numMyElementsX;

    for (int i = 0; i < numMyElementsX; ++i)
    {
      boundary.setGlobalConnectivity(offset2 + i, 0, offset + i);
      boundary.setGlobalConnectivity(offset2 + i, 1, offset + i + 1);
    }
  }

  if (px == 0)
  {
    int offset = numGlobalVerticesX * py * numMyElementsY;
    int offset2 = 2 * numGlobalElementsX + numGlobalElementsY + py * numMyElementsY;

    for (int i = 0; i < numMyElementsY; ++i)
    {
      boundary.setGlobalConnectivity(offset2 + i, 0, offset + i * numGlobalVerticesX);
      boundary.setGlobalConnectivity(offset2 + i, 1, offset + (i + 1) * numGlobalVerticesX);
    }
  }

  boundary.freezeConnectivity();

  if (py == 0)
  {
    int offset = px * numMyElementsX + 1;

    for (int i = 0; i < numMyElementsX + 1; ++i)
    {
      boundary.setGlobalCoordinates(offset + i, 0, deltax * (offset + i));
      boundary.setGlobalCoordinates(offset + i, 1, 0.0);
    }
  }

  if (px == numDomainsX - 1)
  {
    int offset = numGlobalVerticesX + py * numMyElementsY - 1;
    int offset2 = px * numMyElementsX;

    for (int i = 0; i < numMyElementsY + 1; ++i)
    {
      boundary.setGlobalCoordinates(offset + i * numGlobalVerticesX, 0, lx);
      boundary.setGlobalCoordinates(offset + i * numGlobalVerticesX, 1, deltay * (offset2 + i));
    }
  }

  if (py == numDomainsY - 1)
  {
    int offset = px * numMyElementsX;
    int offset2 = numGlobalVerticesX * numGlobalElementsY + px * numMyElementsX;

    for (int i = 0; i < numMyElementsX + 1; ++i)
    {
      boundary.setGlobalCoordinates(offset2 + i, 0, deltax * (offset + i));
      boundary.setGlobalCoordinates(offset2 + i, 1, ly);
    }
  }

  if (px == 0)
  {
    int offset = numGlobalVerticesX * py * numMyElementsY;
    int offset2 = py * numMyElementsX;

    for (int i = 0; i < numMyElementsY + 1; ++i)
    {
      boundary.setGlobalCoordinates(offset + i * numGlobalVerticesX, 0, 0.0);
      boundary.setGlobalCoordinates(offset + i * numGlobalVerticesX, 1, deltay * (offset2 + i));
    }
  }

  boundary.freezeCoordinates();
}

// ===========================================================================
void Galeri::grid::Generator::
getCubeWithHexs(Epetra_Comm& comm, 
                const int numGlobalElementsX, const int numGlobalElementsY, const int numGlobalElementsZ,
                const int numDomainsX, const int numDomainsY, const int numDomainsZ,
                Galeri::grid::Loadable& domain, Galeri::grid::Loadable& boundary)
{
  map<char, int> numGlobalElements;
  numGlobalElements['x'] = numGlobalElementsX;
  numGlobalElements['y'] = numGlobalElementsY;
  numGlobalElements['z'] = numGlobalElementsZ;
  numGlobalElements['p'] = numGlobalElements['x'] * numGlobalElements['y'];
  numGlobalElements['q'] = numGlobalElements['y'] * numGlobalElements['z'];
  numGlobalElements['r'] = numGlobalElements['x'] * numGlobalElements['z'];
  numGlobalElements['a'] = numGlobalElements['x'] * numGlobalElements['y'] * numGlobalElements['z'];

  map<char, int> numGlobalVertices;
  numGlobalVertices['x'] = numGlobalElements['x'] + 1;
  numGlobalVertices['y'] = numGlobalElements['y'] + 1;
  numGlobalVertices['z'] = numGlobalElements['z'] + 1;
  numGlobalVertices['p'] = numGlobalVertices['x'] * numGlobalVertices['y'];
  numGlobalVertices['q'] = numGlobalVertices['y'] * numGlobalVertices['z'];
  numGlobalVertices['r'] = numGlobalVertices['x'] * numGlobalVertices['z'];
  numGlobalVertices['a'] = numGlobalVertices['x'] * numGlobalVertices['y'] * numGlobalVertices['z'];

  map<char, double> length;
  length['x'] = 1.0;
  length['y'] = 1.0;
  length['z'] = 1.0;

  map<char, int> numDomains;
  numDomains['x'] = numDomainsX;
  numDomains['y'] = numDomainsY;
  numDomains['z'] = numDomainsZ;
  numDomains['p'] = numDomains['x'] * numDomains['y'];

  map<char, int> numMyElements;
  numMyElements['x'] = numGlobalElements['x'] / numDomains['x'];
  numMyElements['y'] = numGlobalElements['y'] / numDomains['y'];
  numMyElements['z'] = numGlobalElements['z'] / numDomains['z'];
  numMyElements['a'] = numMyElements['x'] * numMyElements['y'] * numMyElements['z'];

  map<char, int> numMyVertices;
  numMyVertices['x'] = numMyElements['x'] + 1;
  numMyVertices['y'] = numMyElements['y'] + 1;
  numMyVertices['z'] = numMyElements['z'] + 1;
  numMyVertices['a'] = numMyVertices['x'] * numMyVertices['y'] * numMyVertices['z'];

  map<char, int> pos;
  pos['z'] = comm.MyPID() / numDomains['p'];
  pos['y'] = (comm.MyPID() - pos['z'] * numDomains['p']) / numDomains['x'];
  pos['x'] = (comm.MyPID() - pos['z'] * numDomains['p']) % numDomains['x'];

  vector<int> list;

  for (int iz = 0; iz < numMyElements['z']; ++iz)
  {
    for (int iy = 0; iy < numMyElements['y']; ++iy)
    {
      for (int ix = 0; ix < numMyElements['x']; ++ix)
      {
        int IX = ix + numMyElements['x'] * pos['x'];
        int IY = iy + numMyElements['y'] * pos['y'];
        int IZ = iz + numMyElements['z'] * pos['z'];

        int GEID = IZ * numGlobalElements['p'] + IY * numGlobalElements['x'] + IX;
        list.push_back(GEID);
      }
    }
  }

  Hex hex;
  domain.initialize(comm, numGlobalElements['a'], list.size(), hex, &list[0]);

  for (int iz = 0; iz < numMyElements['z']; ++iz)
  {
    for (int iy = 0; iy < numMyElements['y']; ++iy)
    {
      for (int ix = 0; ix < numMyElements['x']; ++ix)
      {
        int IX = ix + numMyElements['x'] * pos['x'];
        int IY = iy + numMyElements['y'] * pos['y'];
        int IZ = iz + numMyElements['z'] * pos['z'];

        int GEID = IZ * numGlobalElements['p'] + IY * numGlobalElements['x'] + IX;
        int offset = IZ * numGlobalVertices['p'] + IY * numGlobalVertices['x'] + IX;

        domain.setGlobalConnectivity(GEID, 0, offset);
        domain.setGlobalConnectivity(GEID, 1, offset + 1);
        domain.setGlobalConnectivity(GEID, 2, offset + numGlobalVertices['x'] + 1);
        domain.setGlobalConnectivity(GEID, 3, offset + numGlobalVertices['x']);
        domain.setGlobalConnectivity(GEID, 4, offset + numGlobalVertices['p']);
        domain.setGlobalConnectivity(GEID, 5, offset + numGlobalVertices['p'] + 1);
        domain.setGlobalConnectivity(GEID, 6, offset + numGlobalVertices['p'] + numGlobalVertices['x'] + 1);
        domain.setGlobalConnectivity(GEID, 7, offset + numGlobalVertices['p'] + numGlobalVertices['x']);
      }
    }
  }

  domain.freezeConnectivity();

  double hx = length['x'] / numGlobalElements['x'];
  double hy = length['y'] / numGlobalElements['y'];
  double hz = length['z'] / numGlobalElements['z'];

  for (int iz = 0; iz < numMyVertices['z']; ++iz)
  {
    for (int iy = 0; iy < numMyVertices['y']; ++iy)
    {
      for (int ix = 0; ix < numMyVertices['x']; ++ix)
      {
        int IX = ix + numMyElements['x'] * pos['x'];
        int IY = iy + numMyElements['y'] * pos['y'];
        int IZ = iz + numMyElements['z'] * pos['z'];

        int GVID = IZ * numGlobalVertices['p'] + IY * numGlobalVertices['x'] + IX;
        domain.setGlobalCoordinates(GVID, 0, hx * IX);
        domain.setGlobalCoordinates(GVID, 1, hy * IY);
        domain.setGlobalCoordinates(GVID, 2, hz * IZ);
      }
    }
  }

  domain.freezeCoordinates();

  // ===================== //
  // now fix the boudaries //
  // ===================== //

  int numMyBoundaries = 0;

  if (pos['z'] == 0)                   numMyBoundaries += numMyVertices['p'];
  if (pos['z'] == numDomains['z'] - 1) numMyBoundaries += numMyVertices['p'];
  if (pos['x'] == 0)                   numMyBoundaries += numMyVertices['q'];
  if (pos['x'] == numDomains['y'] - 1) numMyBoundaries += numMyVertices['q'];
  if (pos['y'] == 0)                   numMyBoundaries += numMyVertices['r'];
  if (pos['y'] == numDomains['z'] - 1) numMyBoundaries += numMyVertices['r'];

  list.resize(0);

  // bottom 
  if (pos['z'] == 0)
  {
    for (int iy = 0; iy < numMyVertices['y']; ++iy)
    {
      for (int ix = 0; ix < numMyVertices['x']; ++ix)
      {
        int IX = ix + numMyElements['x'] * pos['x'];
        int IY = iy + numMyElements['y'] * pos['y'];

        int GBID = IY * numGlobalVertices['x'] + IX;
        list.push_back(GBID);
      }
    }
  }

  // top 
  if (pos['z'] == numDomains['x'] - 1)
  {
    for (int iy = 0; iy < numMyVertices['y']; ++iy)
    {
      for (int ix = 0; ix < numMyVertices['x']; ++ix)
      {
        int IX = ix + numMyElements['x'] * pos['x'];
        int IY = iy + numMyElements['y'] * pos['y'];

        int GBID = IY * numGlobalVertices['x'] + IX + numGlobalVertices['p'] * numGlobalElements['z'];
        list.push_back(GBID);
      }
    }
  }

  // front 
  for (int iz = 0; iz < numMyVertices['z']; ++iz)
  {
    for (int ix = 0; ix < numMyVertices['x']; ++ix)
    {
      int IX = ix + numMyElements['x'] * pos['x'];
      int IZ = iz + numMyElements['z'] * pos['z'];

      int GBID = IZ * numGlobalVertices['p'] + IX;
      list.push_back(GBID);
    }
  }

  // rear 
  for (int iz = 0; iz < numMyVertices['z']; ++iz)
  {
    for (int ix = 0; ix < numMyVertices['x']; ++ix)
    {
      int IX = ix + numMyElements['x'] * pos['x'];
      int IZ = iz + numMyElements['z'] * pos['z'];

      int GBID = (IZ + 1) * numGlobalVertices['p'] - numGlobalVertices['x'] + IX;
      list.push_back(GBID);
    }
  }

  // left 
  for (int iz = 0; iz < numMyVertices['z']; ++iz)
  {
    for (int iy = 0; iy < numMyVertices['y']; ++iy)
    {
      int IY = iy + numMyElements['y'] * pos['y'];
      int IZ = iz + numMyElements['z'] * pos['z'];

      int GBID = IZ * numGlobalVertices['p'] + IY * numGlobalVertices['x'];
      list.push_back(GBID);
    }
  }

  // right 
  for (int iz = 0; iz < numMyVertices['z']; ++iz)
  {
    for (int iy = 0; iy < numMyVertices['y']; ++iy)
    {
      int IY = iy + numMyElements['y'] * pos['y'];
      int IZ = iz + numMyElements['z'] * pos['z'];

      int GBID = IZ * numGlobalVertices['p'] + (IY + 1) * numGlobalVertices['x'] - 1;
      list.push_back(GBID);
    }
  }

  Point point;
  boundary.initialize(comm, -1, list.size(), point, &list[0]);

  // ===================== //
  // now insert vertex IDs //
  // ===================== //

  // bottom 
  if (pos['z'] == 0)
  {
    for (int iy = 0; iy < numMyVertices['y']; ++iy)
    {
      for (int ix = 0; ix < numMyVertices['x']; ++ix)
      {
        int IX = ix + numMyElements['x'] * pos['x'];
        int IY = iy + numMyElements['y'] * pos['y'];

        int GBID = IY * numGlobalVertices['x'] + IX;
        boundary.setGlobalConnectivity(GBID, 0, GBID);
      }
    }
  }

  // top 
  if (pos['z'] == numDomains['x'] - 1)
  {
    for (int iy = 0; iy < numMyVertices['y']; ++iy)
    {
      for (int ix = 0; ix < numMyVertices['x']; ++ix)
      {
        int IX = ix + numMyElements['x'] * pos['x'];
        int IY = iy + numMyElements['y'] * pos['y'];

        int GBID = IY * numGlobalVertices['x'] + IX + numGlobalVertices['p'] * numGlobalElements['z'];
        boundary.setGlobalConnectivity(GBID, 0, GBID);
      }
    }
  }

  // front 
  for (int iz = 0; iz < numMyVertices['z']; ++iz)
  {
    for (int ix = 0; ix < numMyVertices['x']; ++ix)
    {
      int IX = ix + numMyElements['x'] * pos['x'];
      int IZ = iz + numMyElements['z'] * pos['z'];

      int GBID = IZ * numGlobalVertices['p'] + IX;
      boundary.setGlobalConnectivity(GBID, 0, GBID);
    }
  }

  // rear 
  for (int iz = 0; iz < numMyVertices['z']; ++iz)
  {
    for (int ix = 0; ix < numMyVertices['x']; ++ix)
    {
      int IX = ix + numMyElements['x'] * pos['x'];
      int IZ = iz + numMyElements['z'] * pos['z'];

      int GBID = (IZ + 1) * numGlobalVertices['p'] - numGlobalVertices['x'] + IX;
      boundary.setGlobalConnectivity(GBID, 0, GBID);
    }
  }

  // left 
  for (int iz = 0; iz < numMyVertices['z']; ++iz)
  {
    for (int iy = 0; iy < numMyVertices['y']; ++iy)
    {
      int IY = iy + numMyElements['y'] * pos['y'];
      int IZ = iz + numMyElements['z'] * pos['z'];

      int GBID = IZ * numGlobalVertices['p'] + IY * numGlobalVertices['x'];
      boundary.setGlobalConnectivity(GBID, 0, GBID);
    }
  }

  // right 
  for (int iz = 0; iz < numMyVertices['z']; ++iz)
  {
    for (int iy = 0; iy < numMyVertices['y']; ++iy)
    {
      int IY = iy + numMyElements['y'] * pos['y'];
      int IZ = iz + numMyElements['z'] * pos['z'];

      int GBID = IZ * numGlobalVertices['p'] + (IY + 1) * numGlobalVertices['x'] - 1;
      boundary.setGlobalConnectivity(GBID, 0, GBID);
    }
  }

  boundary.freezeConnectivity();

  // bottom
  if (pos['z'] == 0)
  {
    for (int iy = 0; iy < numMyVertices['y']; ++iy)
    {
      for (int ix = 0; ix < numMyVertices['x']; ++ix)
      {
        int IX = ix + numMyElements['x'] * pos['x'];
        int IY = iy + numMyElements['y'] * pos['y'];

        int GBID = IY * numGlobalVertices['x'] + IX;

        boundary.setGlobalCoordinates(GBID, 0, hx * IX);
        boundary.setGlobalCoordinates(GBID, 1, hy * IY);
        boundary.setGlobalCoordinates(GBID, 2, 0.0);
      }
    }
  }

  // top
  if (pos['z'] == numDomains['x'] - 1)
  {
    for (int iy = 0; iy < numMyVertices['y']; ++iy)
    {
      for (int ix = 0; ix < numMyVertices['x']; ++ix)
      {
        int IX = ix + numMyElements['x'] * pos['x'];
        int IY = iy + numMyElements['y'] * pos['y'];

        int GBID = IY * numGlobalVertices['x'] + IX + numGlobalVertices['p'] * numGlobalElements['z'];

        boundary.setGlobalCoordinates(GBID, 0, hx * IX);
        boundary.setGlobalCoordinates(GBID, 1, hy * IY);
        boundary.setGlobalCoordinates(GBID, 2, length['z']);
      }
    }
  }

  // front
  for (int iz = 0; iz < numMyVertices['z']; ++iz)
  {
    for (int ix = 0; ix < numMyVertices['x']; ++ix)
    {
      int IX = ix + numMyElements['x'] * pos['x'];
      int IZ = iz + numMyElements['z'] * pos['z'];

      int GBID = IZ * numGlobalVertices['p'] + IX;

      boundary.setGlobalCoordinates(GBID, 0, hx * IX);
      boundary.setGlobalCoordinates(GBID, 1, 0.0);
      boundary.setGlobalCoordinates(GBID, 2, hz * IZ);
    }
  }

  // rear
  for (int iz = 0; iz < numMyVertices['z']; ++iz)
  {
    for (int ix = 0; ix < numMyVertices['x']; ++ix)
    {
      int IX = ix + numMyElements['x'] * pos['x'];
      int IZ = iz + numMyElements['z'] * pos['z'];

      int GBID = (IZ + 1) * numGlobalVertices['p'] - numGlobalVertices['x'] + IX;

      boundary.setGlobalCoordinates(GBID, 0, hx * IX);
      boundary.setGlobalCoordinates(GBID, 1, length['y']);
      boundary.setGlobalCoordinates(GBID, 2, hz * IZ);
    }
  }

  // left 
  for (int iz = 0; iz < numMyVertices['z']; ++iz)
  {
    for (int iy = 0; iy < numMyVertices['y']; ++iy)
    {
      int IY = iy + numMyElements['y'] * pos['y'];
      int IZ = iz + numMyElements['z'] * pos['z'];

      int GBID = IZ * numGlobalVertices['p'] + IY * numGlobalVertices['x'];

      boundary.setGlobalCoordinates(GBID, 0, 0.0);
      boundary.setGlobalCoordinates(GBID, 1, hy * IY);
      boundary.setGlobalCoordinates(GBID, 2, hz * IZ);
    }
  }

  // right 
  for (int iz = 0; iz < numMyVertices['z']; ++iz)
  {
    for (int iy = 0; iy < numMyVertices['y']; ++iy)
    {
      int IY = iy + numMyElements['y'] * pos['y'];
      int IZ = iz + numMyElements['z'] * pos['z'];

      int GBID = IZ * numGlobalVertices['p'] + (IY + 1) * numGlobalVertices['x'] - 1;

      boundary.setGlobalCoordinates(GBID, 0, length['x']);
      boundary.setGlobalCoordinates(GBID, 1, hy * IY);
      boundary.setGlobalCoordinates(GBID, 2, hz * IZ);
    }
  }

  boundary.freezeCoordinates();
}
