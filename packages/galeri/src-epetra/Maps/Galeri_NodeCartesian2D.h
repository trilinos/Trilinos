// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef GALERI_NODECARTESIAN2D_H
#define GALERI_NODECARTESIAN2D_H
#include "Galeri_Exception.h"
#include "Galeri_Utils.h"
#ifndef HAVE_MPI
#include "Epetra_Comm.h"
#else
#include "Epetra_MpiComm.h"
#endif
#include "Epetra_Map.h"

namespace Galeri {

namespace Maps {

template<typename int_type>
inline
Epetra_Map* 
TNodeCartesian2D(const Epetra_Comm& Comm, const Epetra_Comm & NodeComm, const int MyNodeID,
                const int nx, const int ny,
                const int ndx, const int ndy, const int px, const int py )
{
  if (nx <= 0 || ny <= 0 || (ndx > nx) || (ndy > ny) || px <= 0 || py <= 0)
    throw(Exception(__FILE__, __LINE__,
                    "Incorrect input parameter to Maps::NodeCartesian2D()",
                    "nx = " + toString(nx) +
                    ", ny = " + toString(ny) +
                    ", ndx = " + toString(ndx) +
                    ", ndy = " + toString(ndy) +
                    ", px = " + toString(px) +
                    ", py = " + toString(py)));
  int NodePID=NodeComm.MyPID();
  // Compute nodal position
  int nxpid = MyNodeID % ndx;
  int nypid = MyNodeID / ndx;

  // Compute position w/i node
  int pxpid = NodePID % px;
  int pypid = NodePID / px;

  // Compute effective PID, mx, my
  int mx = ndx*px;
  int my = ndy*py;
  int ppn=px*py;
  int EffectivePID = nypid*ndx*ppn + nxpid*px + pypid*mx + pxpid;

  // Build the map based on the Effective PID - code copied from Galeri_NodeCartesian2D.h
  int xpid = EffectivePID % mx;
  int ypid = EffectivePID / mx;

  int PerProcSmallXDir = (int) (((double) nx)/((double) mx));
  int PerProcSmallYDir = (int) (((double) ny)/((double) my));
  int NBigXDir         = nx - PerProcSmallXDir*mx;
  int NBigYDir         = ny - PerProcSmallYDir*my;
  int startx,starty,endx,endy;
  
  if (xpid < NBigXDir) startx =  xpid*(PerProcSmallXDir+1);
  else                 startx = (xpid-NBigXDir)*PerProcSmallXDir + 
                                      NBigXDir*(PerProcSmallXDir+1);
  endx = startx + PerProcSmallXDir;
  if (xpid < NBigXDir) endx++;

  if (ypid < NBigYDir) starty =  ypid*(PerProcSmallYDir+1);
  else                 starty = (ypid-NBigYDir)*PerProcSmallYDir +
                                      NBigYDir*(PerProcSmallYDir+1);
  endy = starty + PerProcSmallYDir;
  if ( ypid < NBigYDir) endy++;

  int NumMyElements = (endx - startx) * (endy - starty);
  std::vector<int_type> MyGlobalElements(NumMyElements);
  int count = 0;

  for (int i = startx ; i < endx ; ++i) 
    for (int j = starty ; j < endy ; ++j) 
      MyGlobalElements[count++] = i + j * nx;

#if !defined(EPETRA_NO_32BIT_GLOBAL_INDICES) || !defined(EPETRA_NO_64BIT_GLOBAL_INDICES)
  return(new Epetra_Map(nx * ny,  NumMyElements, &MyGlobalElements[0],
                        0, Comm));
#else
  return(new Epetra_Map);
#endif

  
} // TNodeCartesian2D()

#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
Epetra_Map* 
NodeCartesian2D(const Epetra_Comm& Comm, const Epetra_Comm & NodeComm, const int MyNodeID,
                const int nx, const int ny,
                const int ndx, const int ndy, const int px, const int py ) {
  return TNodeCartesian2D<int>(Comm, NodeComm, MyNodeID, nx, ny, ndx, ndy, px, py);
}
#endif

#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
Epetra_Map* 
NodeCartesian2D64(const Epetra_Comm& Comm, const Epetra_Comm & NodeComm, const int MyNodeID,
                const int nx, const int ny,
                const int ndx, const int ndy, const int px, const int py ) {
  return TNodeCartesian2D<long long>(Comm, NodeComm, MyNodeID, nx, ny, ndx, ndy, px, py);
}
#endif

} // namespace BlockMaps
} // namespace Galeri
#endif
