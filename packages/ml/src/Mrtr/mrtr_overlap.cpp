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
#include "mrtr_point.H"
#include "mrtr_edge.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 10/05|
 *----------------------------------------------------------------------*/
MRTR::Overlap::Overlap(MRTR::Segment& sseg, MRTR::Segment& mseg, MRTR::Interface& inter) :
inter_(inter),
sseg_(sseg),
mseg_(mseg),
overlap_(false),
havemxi_(false),
havesxi_(false),
haveline_(false)
{
  if (sseg.Type()!=MRTR::Segment::seg_BiLinearTri || mseg.Type()!=MRTR::Segment::seg_BiLinearTri)
  {
    cout << "***ERR*** MRTR::Overlap::Overlap:\n"
         << "***ERR*** Overlap of other then bilinear triangles not yet implemented\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 10/05|
 *----------------------------------------------------------------------*/
MRTR::Overlap::~Overlap()
{
}

/*----------------------------------------------------------------------*
 |  determine whether 2 triangles overlap (private)          mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MRTR::Overlap::intersect_line2D(int s, int m, double* alpha, double* beta, double* xi)
{
  if (!haveline_)
  {
    cout << "***ERR*** MRTR::Overlap::intersect_line2D:\n"
         << "***ERR*** line information has to be build before intersecting them\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  
  double A[2][2];
  double x[2];
  double b[2];
  A[0][0] = sline_[s][2] - sline_[s][0];
  A[1][0] = sline_[s][3] - sline_[s][1];
  
  A[0][1] = mline_[m][0] - mline_[m][2];
  A[1][1] = mline_[m][1] - mline_[m][3];

  b[0] = mline_[m][0] - sline_[s][0];
  b[1] = mline_[m][1] - sline_[s][1];

  MRTR::solve22(A,x,b);
  
  if (alpha)
    *alpha = x[0];
    
  if (beta)
    *beta = x[1];
    
  if (xi)
  {
    xi[0] = sline_[s][0] + x[0]*(sline_[s][2] - sline_[s][0]);
    xi[1] = sline_[s][1] + x[0]*(sline_[s][3] - sline_[s][1]);
  }

  if ( 0. <= x[0] && x[0] <= 1. && 0. <= x[1] && x[1] <= 1. )
  {
    if (xi)
    {
      double xitest[2];
      xitest[0] = mline_[m][0] + x[1]*(mline_[m][2] - mline_[m][0]);
      xitest[1] = mline_[m][1] + x[1]*(mline_[m][3] - mline_[m][1]);
      if ( abs(xi[0]-xitest[0])>1.0e-10 || abs(xi[1]-xitest[1])>1.0e-10)
      {
        cout << " xi     " << xi[0] << "/" << xi[1] << endl; 
        cout << " xitest " << xitest[0] << "/" << xitest[1] << endl; 
        cout << "***WRN*** MRTR::Overlap::intersect_line2D:\n"
             << "***WRN*** weired result in intersecting lines\n"
             << "***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      }
    }
    return true; 
  }
  return false;
}

/*----------------------------------------------------------------------*
 |  build line info from triangles (private)                 mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MRTR::Overlap::build_lines()
{
  if (!havemxi_)
  {
    cout << "***ERR*** MRTR::Overlap::build_lines:\n"
         << "***ERR*** projection of master element nodes has to be done before\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }

  // slave segment, line 0
  sseg_.LocalCoordinatesOfNode(0,&sline_[0][0]);
  sseg_.LocalCoordinatesOfNode(1,&sline_[0][2]);
  // slave segment, line 1
  sseg_.LocalCoordinatesOfNode(1,&sline_[1][0]);
  sseg_.LocalCoordinatesOfNode(2,&sline_[1][2]);
  // slave segment, line 2
  sseg_.LocalCoordinatesOfNode(2,&sline_[2][0]);
  sseg_.LocalCoordinatesOfNode(0,&sline_[2][2]);
  
  // master segment, line 0
  mline_[0][0] = mxi_[0][0];
  mline_[0][1] = mxi_[0][1];
  mline_[0][2] = mxi_[1][0];
  mline_[0][3] = mxi_[1][1];
  // master segment, line 1
  mline_[1][0] = mxi_[1][0];
  mline_[1][1] = mxi_[1][1];
  mline_[1][2] = mxi_[2][0];
  mline_[1][3] = mxi_[2][1];
  // master segment, line 2
  mline_[2][0] = mxi_[2][0];
  mline_[2][1] = mxi_[2][1];
  mline_[2][2] = mxi_[0][0];
  mline_[2][3] = mxi_[0][1];

  return true;
}

/*----------------------------------------------------------------------*
 |  project master nodes onto slave element (private)        mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MRTR::Overlap::build_mxi()
{
  // project the master segment's nodes onto the slave segment
  for (int i=0; i<4; ++i) min_[i] = false;
               nnode_ = mseg_.Nnode();
  MRTR::Node** mnode  = mseg_.Nodes();
  MRTR::Projector projector(inter_.IsOneDimensional());
  for (int i=0; i<nnode_; ++i)
  {
    // project node i onto sseg
    projector.ProjectNodetoSegment_SegmentNormal(*mnode[i],sseg_,mxi_[i]);
    // check whether i is inside sseg
    if (mxi_[i][0]<=1. && mxi_[i][1]<=abs(1.-mxi_[i][0]) && mxi_[i][0]>=0. && mxi_[i][1]>=0.)
    {
      min_[i]   = true;
      cout << "OVERLAP: master node in: " << i << endl;
    }
  }
  havemxi_ = true;
  return true;
}

/*----------------------------------------------------------------------*
 |  project slave nodes onto master element (private)        mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MRTR::Overlap::build_sxi()
{
  for (int i=0; i<4; ++i) sin_[i] = false;
  MRTR::Node** snode  = sseg_.Nodes();
  MRTR::Projector projector(inter_.IsOneDimensional());
  for (int i=0; i<nnode_; ++i)
  {
    // project slave node i onto mseg
    projector.ProjectNodetoSegment_NodalNormal(*snode[i],mseg_,sxi_[i]);
    // check whether i is inside sseg
    if (sxi_[i][0]<=1. && sxi_[i][1]<=abs(1.-sxi_[i][0]) && sxi_[i][0]>=0. && sxi_[i][1]>=0.)
    {
      sin_[i] = true;
      cout << "OVERLAP: slave node in: " << i << endl;
    }
  }
  havesxi_ = true;
  return true;
}

/*----------------------------------------------------------------------*
 |  get a view from a point (private)                        mwgee 10/05|
 *----------------------------------------------------------------------*/
MRTR::Point* MRTR::Overlap::PointView(int id)
{
  map<int,MRTR::Point*>::iterator curr = p_.find(id);
  if (curr != p_.end())
    return curr->second;
  else
    return NULL;
}

/*----------------------------------------------------------------------*
 |  get a view from a edge (private)                         mwgee 10/05|
 *----------------------------------------------------------------------*/
MRTR::Edge* MRTR::Overlap::EdgeView(int id)
{
  map<int,MRTR::Edge*>::iterator curr = e_.find(id);
  if (curr != e_.end())
    return curr->second;
  else
    return NULL;
}

/*----------------------------------------------------------------------*
 |  add point to polygon (private)                           mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MRTR::Overlap::AddPoint(MRTR::Point* p)
{
  map<int,MRTR::Point*>::iterator curr = p_.find(p->Id());
  if (curr != p_.end())
  {
    delete curr->second;
    curr->second = p;
  }
  else
    p_.insert(pair<int,MRTR::Point*>(p->Id(),p));
  return true;
}

/*----------------------------------------------------------------------*
 |  add edge to polygon (private)                            mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MRTR::Overlap::AddEdge(MRTR::Edge* e)
{
  map<int,MRTR::Edge*>::iterator curr = e_.find(e->Id());
  if (curr != e_.end())
  {
    delete curr->second;
    curr->second = e;
  }
  else
    e_.insert(pair<int,MRTR::Edge*>(e->Id(),e));
  return true;
}


/*----------------------------------------------------------------------*
 |  determine whether 2 triangles overlap (public)           mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MRTR::Overlap::HaveOverlap()
{
  // maybe we computed this before
  if (overlap_) 
    return true;
  
  // first we project the master segment's nodes onto the slave segment
  for (int i=0; i<4; ++i) min_[i] = false;
               nnode_ = mseg_.Nnode();
  MRTR::Node** mnode  = mseg_.Nodes();
  MRTR::Projector projector(inter_.IsOneDimensional());
  for (int i=0; i<nnode_; ++i)
  {
    // project node i onto sseg
    projector.ProjectNodetoSegment_SegmentNormal(*mnode[i],sseg_,mxi_[i]);
    // check whether i is inside sseg
    if (mxi_[i][0]<=1. && mxi_[i][1]<=abs(1.-mxi_[i][0]) && mxi_[i][0]>=0. && mxi_[i][1]>=0.)
    {
      overlap_ = true;
      min_[i]   = true;
      //cout << "OVERLAP: master node in: " << i << endl;
      return true;
    }
  }
  havemxi_ = true;

  // none of the master segments project into the slave element
  // The slave nodes might project into the master element
  for (int i=0; i<4; ++i) sin_[i] = false;
  MRTR::Node** snode  = sseg_.Nodes();
  for (int i=0; i<nnode_; ++i)
  {
    // project slave node i onto mseg
    projector.ProjectNodetoSegment_NodalNormal(*snode[i],mseg_,sxi_[i]);
    // check whether i is inside sseg
    if (sxi_[i][0]<=1. && sxi_[i][1]<=abs(1.-sxi_[i][0]) && sxi_[i][0]>=0. && sxi_[i][1]>=0.)
    {
      sin_[i] = true;
      //cout << "OVERLAP: slave node in: " << i << endl;
      overlap_ = true;
      return true;
    }
  }
  havesxi_ = true;
  
  // determine whether the triangles' lines intersect
  
  // build the information about sseg_ and mseg_'s edges
  haveline_ = build_lines();
  
  // loop lines of slave triangles and intersect them with master lines
  // if any intersection is found, then segments overlap
  for (int i=0; i<3; ++i)
  {
    for (int j=0; j<3; ++j)
    {
      double eta[2];
      bool ok = intersect_line2D(i,j,NULL,NULL,eta);
      if (ok) 
      {
        //cout << "OVERLAP: sline " << i << " intersects mline " << j << endl;
        overlap_ = true;
        return true;
      }
    }
  }
  
  // there doesn't seem to be an intersection between these 2 triangles
  return false;
}

/*----------------------------------------------------------------------*
 |  compute the overlap (public)                             mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MRTR::Overlap::ComputeOverlap()
{
  // project master nodes onto slave segment if not already done
  if (!havemxi_)
    havemxi_ = build_mxi();

  // project slave nodes onto master segment if not already done
  if (!havesxi_)
    havesxi_ = build_sxi();
  
  // compute information about the edges of the slave/master triangle
  if (!haveline_)
    haveline_ = build_lines();

#if 1
  // debugging output
  cout << "OVERLAP: Slave Segment\n" << sseg_;
  cout << "OVERLAP: Master Segment\n" << mseg_;
  cout << "node 0: " << mxi_[0][0] << "/" << mxi_[0][1] << endl
       << "     1: " << mxi_[1][0] << "/" << mxi_[1][1] << endl
       << "     2: " << mxi_[2][0] << "/" << mxi_[2][1] << endl;
#endif

  // create master nodes and edges and push them on the polygon
  // do nodes
  MRTR::Point* point;
  point = new MRTR::Point(0,&mxi_[0][0]);
  AddPoint(point);
  point = new MRTR::Point(1,&mxi_[1][0]);
  AddPoint(point);
  point = new MRTR::Point(2,&mxi_[2][0]);
  AddPoint(point);
  
  // do edges
  MRTR::Edge* edge;
  edge = new MRTR::Edge(0,0,1);
  AddEdge(edge);    
  edge = new MRTR::Edge(1,1,2);    
  AddEdge(edge);    
  edge = new MRTR::Edge(2,2,0);    
  AddEdge(edge);    
  
  

  return false;
}
















































#endif // TRILINOS_PACKAGE
