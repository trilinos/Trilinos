// @HEADER
// ************************************************************************
//
//           Galeri: Finite Element and Matrix Generation Package
//                 Copyright (2006) ETHZ/Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
//
// Questions about Galeri? Contact Marzio Sala (marzio.sala _AT_ gmail.com)
//
// ************************************************************************
// @HEADER

#ifndef GALERI_NODECARTESIAN2D_H
#define GALERI_NODECARTESIAN2D_H
#include "Galeri_Exception.h"
#include "Galeri_Utils.h"
#include "Epetra_Comm.h"
#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"

namespace Galeri {

namespace Maps {

inline
Epetra_Map* 
NodeCartesian2D(const Epetra_Comm& Comm, const Epetra_Comm & NodeComm, const int MyNodeID,
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
  vector<int> MyGlobalElements(NumMyElements);
  int count = 0;

  for (int i = startx ; i < endx ; ++i) 
    for (int j = starty ; j < endy ; ++j) 
      MyGlobalElements[count++] = i + j * nx;

  return(new Epetra_Map(nx * ny,  NumMyElements, &MyGlobalElements[0],
                        0, Comm));

  
} // NodeCartesian2D()

} // namespace BlockMaps
} // namespace Galeri
#endif
