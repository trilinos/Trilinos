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

#include "mrtr_overlap.H"
#include "mrtr_projector.H"
#include "mrtr_node.H"
#include "mrtr_segment.H"
#include "mrtr_interface.H"
#include "mrtr_utils.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 10/05|
 *----------------------------------------------------------------------*/
MRTR::Overlap::Overlap(MRTR::Segment& sseg, MRTR::Segment& mseg, MRTR::Interface& inter) :
inter_(inter),
sseg_(sseg),
mseg_(mseg),
overlap_(false)
{
  if (sseg.Type()!=MRTR::Segment::seg_BiLinearTri || mseg.Type()!=MRTR::Segment::seg_BiLinearTri)
  {
    cout << "***ERR*** MRTR::Overlap::Overlap:\n"
         << "***ERR*** Overlap of other then bilinear triangles not yet implemented\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  // first we project the master segment's nodes onto the slave segment
  for (int i=0; i<4; ++i) in_[i] = false;
               nnode_ = mseg.Nnode();
  MRTR::Node** mnode  = mseg.Nodes();
  MRTR::Projector projector(inter_.IsOneDimensional());
  for (int i=0; i<nnode_; ++i)
  {
    // project node i onto sseg
    projector.ProjectNodetoSegment_SegmentNormal(*mnode[i],sseg,xi_[i]);
    // check whether i is inside sseg
    if (xi_[i][0]<=1. && xi_[i][1]<=abs(1.-xi_[i][0]) && xi_[i][0]>=0. && xi_[i][1]>=0.)
    {
      overlap_ = true;
      in_[i]   = true;
    }
  }
  cout << "-----in: " << in_[0] << "/" << in_[1] << "/" << in_[2] << endl;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 10/05|
 *----------------------------------------------------------------------*/
MRTR::Overlap::~Overlap()
{
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MRTR::Overlap::HaveOverlap()
{
  // if one of the master nodes projects into the slave element, there is
  // an overlap for sure
  for (int i=0; i<nnode_; ++i)
    if (in_[i])
      return true;

  // none of the master segments project into the slave element
  // The slave nodes might project into the master element
  double xi[4][2];
  MRTR::Node** snode  = sseg_.Nodes();
  MRTR::Projector projector(inter_.IsOneDimensional());
  for (int i=0; i<nnode_; ++i)
  {
    // project slave node i onto mseg
    projector.ProjectNodetoSegment_NodalNormal(*snode[i],mseg_,xi[i]);
    // check whether i is inside sseg
    if (xi[i][0]<=1. && xi[i][1]<=abs(1.-xi[i][0]) && xi[i][0]>=0. && xi[i][1]>=0.)
    {
      cout << "-----Slave node projects into master\n";
      overlap_ = true;
      return true;
    }
  }
  
  

  return false;
}



















































#endif // TRILINOS_PACKAGE
