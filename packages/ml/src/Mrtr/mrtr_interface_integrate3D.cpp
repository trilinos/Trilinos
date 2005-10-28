/*
#@HEADER
# ************************************************************************
#
#                 Copyright (2002) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
# Questions? Contact Jonathan Hu (jhu@sandia.gov) or Ray Tuminaro 
# (rstumin@sandia.gov).
#
# ************************************************************************
#@HEADER
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
#ifdef TRILINOS_PACKAGE

#include <ctime>
#include <vector>

#include "mrtr_interface.H"
#include "mrtr_utils.H"
#include "mrtr_pnode.H"
#include "mrtr_segment.H"
#include "mrtr_integrator.H"
#include "mrtr_projector.H"
#include "mrtr_overlap.H"

#include "Epetra_SerialDenseMatrix.h"

/*----------------------------------------------------------------------*
 |  make mortar integration of master/slave side in 3D (2D interface)   |
 *----------------------------------------------------------------------*/
bool MRTR::Interface::Integrate_3D(Epetra_CrsMatrix& M,
                                   Epetra_CrsMatrix& D)
{ 
  if (!IsComplete())
  {
    if (gcomm_.MyPID()==0)
      cout << "***ERR*** MRTR::Interface::Integrate_3D:\n"
           << "***ERR*** Complete() not called on interface " << Id_ << "\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  if (!lComm()) return true;

  // get the sides
  int mside = MortarSide();
  int sside = OtherSide(mside);

  
  // loop over all segments of slave side
  map<int,MRTR::Segment*>::iterator scurr;
  for (scurr=rseg_[sside].begin(); scurr!=rseg_[sside].end(); ++scurr)
  {
    // the segment to be integrated
    MRTR::Segment* actsseg = scurr->second;

#if 1
    cout << "\nActive sseg id " << actsseg->Id() << "\n\n";
#endif

    // check whether I own at least one of the nodes on this slave segment
    int nnode = actsseg->Nnode();
    MRTR::Node** nodes = actsseg->Nodes();
    bool foundone = false;
    for (int i=0; i<nnode; ++i)
      if (NodePID(nodes[i]->Id()) == lComm()->MyPID())
      {
        foundone = true;
        break;
      }
    // if none of the nodes belongs to me, do nothing on this segment
    if (!foundone) continue;
    
    // loop over all segments on the master side
    map<int,MRTR::Segment*>::iterator mcurr;
    for (mcurr=rseg_[mside].begin(); mcurr!=rseg_[mside].end(); ++mcurr)    
    {
      MRTR::Segment* actmseg = mcurr->second;
      
#if 1
    cout << "Active mseg id " << actmseg->Id() << endl;
#endif
      // if there is an overlap, integrate the pair
      // (whether there is an overlap or not will be checked inside)
      Integrate_3D_Section(*actsseg,*actmseg,M,D);
      
    } // for (mcurr=rseg_[mside].begin(); mcurr!=rseg_[mside].end(); ++mcurr)  
  } // for (scurr=rseg_[sside].begin(); scurr!=rseg_[sside].end(); ++scurr)

  return true;
}

/*----------------------------------------------------------------------*
 | integrate the master/slave side's contribution from the overlap      |
 | of 2 segments (3D version) IF there is an overlap                    |
 *----------------------------------------------------------------------*/
bool MRTR::Interface::Integrate_3D_Section(MRTR::Segment& sseg, 
                                           MRTR::Segment& mseg,
                                           Epetra_CrsMatrix& M,
                                           Epetra_CrsMatrix& D)
{ 
  bool ok;
  
  // if one of the segments is quadratic, we have to do something here
  if (sseg.Type()!=MRTR::Segment::seg_BiLinearTri || mseg.Type()!=MRTR::Segment::seg_BiLinearTri)
  {
    cout << "***ERR*** MRTR::Interface::Integrate_3D_Section:\n"
         << "***ERR*** Integration of other then bilinear triangles not yet implemented\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }

  // first determine whether there is an overlap between sseg and mseg
  // for this purpose, the 'overlapper' class is used
  MRTR::Overlap overlap(sseg,mseg,*this);

  // determine the overlap triangulation
  ok = overlap.ComputeOverlap();
  if (!ok) // there is no overlap
    return true;
  

  return false;
}





















#endif // TRILINOS_PACKAGE
