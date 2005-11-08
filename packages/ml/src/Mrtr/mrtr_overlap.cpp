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
#include "mrtr_pnode.H"
#include "mrtr_segment.H"
#include "mrtr_segment_bilineartri.H"
#include "mrtr_interface.H"
#include "mrtr_utils.H"

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
  p_.clear();
  s_.clear();
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 10/05|
 *----------------------------------------------------------------------*/
MRTR::Overlap::~Overlap()
{
  // destroy the point map
  p_.clear();
  // destroy the segment map
  s_.clear();
}

/*----------------------------------------------------------------------*
 |  copy a polygon of points (private)                       mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MRTR::Overlap::CopyPointPolygon(map<int,RefCountPtr<MRTR::Point> >& from, map<int,RefCountPtr<MRTR::Point> >& to)
{
  map<int,RefCountPtr<MRTR::Point> >::iterator pcurr;
  for (pcurr=from.begin(); pcurr != from.end(); ++pcurr)
    if (pcurr->second != null)
    {
      RefCountPtr<MRTR::Point> tmp = rcp(new MRTR::Point(pcurr->second->Id(),pcurr->second->Xi()));
      to.insert(pair<int,RefCountPtr<MRTR::Point> >(tmp->Id(),tmp));
    }
  return true;
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
               nnode_ = mseg_.Nnode();
  MRTR::Node** mnode  = mseg_.Nodes();
  MRTR::Projector projector(inter_.IsOneDimensional());
  for (int i=0; i<nnode_; ++i)
  {
    // project node i onto sseg
    projector.ProjectNodetoSegment_SegmentNormal(*mnode[i],sseg_,mxi_[i]);
#if 0
    // check whether i is inside sseg
    if (mxi_[i][0]<=1. && mxi_[i][1]<=abs(1.-mxi_[i][0]) && mxi_[i][0]>=0. && mxi_[i][1]>=0.)
    {
      cout << "OVERLAP: master node in: " << i << endl;
    }
#endif    
  }
  havemxi_ = true;
  return true;
}

/*----------------------------------------------------------------------*
 |  project slave nodes onto master element (private)        mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MRTR::Overlap::build_sxi()
{
  // this makes coords of snodes on the sseg local coord system
  sxi_[0][0] = 0.;
  sxi_[0][1] = 0.;
  sxi_[1][0] = 1.;
  sxi_[1][1] = 0.;
  sxi_[2][0] = 0.;
  sxi_[2][1] = 1.;

  havesxi_ = true;
  return true;
}

/*----------------------------------------------------------------------*
 |  outward normal to sseg's edges in local coords (private) mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MRTR::Overlap::build_normal()
{
  sn_[0][0] = 0.;
  sn_[0][1] = -1.;
  sn_[1][0] = 1.;
  sn_[1][1] = 1.;
  sn_[2][0] = -1.;
  sn_[2][1] = 0.;
  return true;
}

/*----------------------------------------------------------------------*
 |  compute the overlap (public)                             mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MRTR::Overlap::ComputeOverlap()
{

  // project master nodes onto slave segment if not already done
  if (!havemxi_)
    havemxi_ = build_mxi();

  // perform a quick test
  bool ok = QuickOverlapTest();
  if (!ok)
    return false;

  // build array of local sxi coords of nodes
  havesxi_ = build_sxi();

  // build outward normal of edges of sseg (in local coords)
  build_normal();
  
  // compute information about the edges of the slave/master triangle
  if (!haveline_)
    haveline_ = build_lines();

  // perform clipping algorithm
  ok = Clipelements();
  if (!ok)
    return false;
  
  // make a triangulation of the overlap polygon
  ok = Triangulization();
  if (!ok)
    return false;

  return true;
}

/*----------------------------------------------------------------------*
 |  find intersection (private)                              mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MRTR::Overlap::Clip_Intersect(const double* N,const double* PE,const double* P0,const double* P1,double* xi)
{
  double P1P0[2];
  P1P0[0] = P1[0] - P0[0];
  P1P0[1] = P1[1] - P0[1];

  // determine the denominator first
  double denom = -(N[0]*P1P0[0] + N[1]*P1P0[1]);
  
  // if the denom is zero, then lines are parallel, no intersection
  if (fabs(denom)<1.0e-10)
    return false;
    
  double alpha = (N[0]*(P0[0]-PE[0]) + N[1]*(P0[1]-PE[1]))/denom;
  
  // alpha is the line parameter of the line P0 - p1
  // if it's outside 0 <= alpha <= 1 there is no intersection
  if (alpha<0.0 || 1.0<alpha)
    return false;
    
  // Compute the coord xi of the intersecting point
  xi[0] = P0[0] + alpha*P1P0[0];
  xi[1] = P0[1] + alpha*P1P0[1];

  //cout << "OVERLAP Clip_Intersect: found intersection xi " << xi[0] << "/" << xi[1] << endl;
  return true;
}

/*----------------------------------------------------------------------*
 |  test point (private)                                     mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MRTR::Overlap::Clip_TestPoint(const double* N, const double* PE, const double* P)
{
  double PPE[2];
  PPE[0] = P[0] - PE[0];
  PPE[1] = P[1] - PE[1];

  double dotproduct = PPE[0]*N[0]+PPE[1]*N[1];
  //cout << "OVERLAP Clip_TestPoint: dotproduct " << dotproduct << endl;

  if (dotproduct>1.0e-10)
    return false;
  else
    return true;
}

/*----------------------------------------------------------------------*
 |  find parameterization alpha for point on line (private)  mwgee 10/05|
 *----------------------------------------------------------------------*/
double MRTR::Overlap::Clip_ParameterPointOnLine(const double* P0,const double* P1,const double* P)
{
  double dist1 = sqrt( (P[0]-P0[0])*(P[0]-P0[0])+(P[1]-P0[1])*(P[1]-P0[1]) );
  double dist2 = sqrt( (P1[0]-P0[0])*(P1[0]-P0[0])+(P1[1]-P0[1])*(P1[1]-P0[1]) );
  if (dist2<1.0e-10) 
  {
    cout << "***ERR*** MRTR::Overlap::Clip_ParameterPointOnLine:\n"
         << "***ERR*** edge length too small, division by near zero\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  return (dist1/dist2);
}

/*----------------------------------------------------------------------*
 |  add point (private)                                      mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MRTR::Overlap::AddSegment(int id, MRTR::Segment* seg)
{
  RefCountPtr<MRTR::Segment> tmp = rcp(seg);
  s_.insert(pair<int,RefCountPtr<MRTR::Segment> >(id,tmp));
  return true;
}

/*----------------------------------------------------------------------*
 |  add point (private)                                      mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MRTR::Overlap::AddPointtoPolygon(const int id,const double* P)
{
  // check whether this point is already in there
  map<int,RefCountPtr<MRTR::Point> >::iterator curr = p_.find(id);
  // it's there
  if (curr != p_.end())
    curr->second->SetXi(P);
  else
  {
    //cout << "OVERLAP Clip_AddPointtoPolygon: added point " << id << endl;
    RefCountPtr<MRTR::Point> p = rcp(new MRTR::Point(id,P));
    p_.insert(pair<int,RefCountPtr<MRTR::Point> >(id,p));
  }
  return true;
}

/*----------------------------------------------------------------------*
 |  add point (private)                                      mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MRTR::Overlap::AddPointtoPolygon(map<int,RefCountPtr<MRTR::Point> >& p,const int id,const double* P)
{
  RefCountPtr<MRTR::Point> point = rcp(new MRTR::Point(id,P));
  p.insert(pair<int,RefCountPtr<MRTR::Point> >(id,point));
  return true;
}

/*----------------------------------------------------------------------*
 |  remove point (private)                                      mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MRTR::Overlap::RemovePointfromPolygon(const int id,const double* P)
{
  // check whether this point is in there
  map<int,RefCountPtr<MRTR::Point> >::iterator curr = p_.find(id);
  if (curr != p_.end())
  {
    curr->second = null;
    p_.erase(id);
    //cout << "OVERLAP Clip_RemovePointfromPolygon: removed point " << id << endl;
    return true;
  }
  else
  {
    //cout << "OVERLAP Clip_RemovePointfromPolygon: do nothing\n";
    return false;
  }
}

/*----------------------------------------------------------------------*
 |  get point view (private)                                 mwgee 10/05|
 *----------------------------------------------------------------------*/
void MRTR::Overlap::PointView(vector<RefCountPtr<MRTR::Point> >& points)
{
  // allocate vector of ptrs
  points.resize(SizePointPolygon());
  
  // get the point views
  int count=0;
  map<int,RefCountPtr<MRTR::Point> >::iterator pcurr;
  for (pcurr=p_.begin(); pcurr != p_.end(); ++pcurr)
  {
    points[count] = pcurr->second;
    ++count;
  }
  if (count != SizePointPolygon())
  {
    cout << "***ERR*** MRTR::Overlap::PointView:\n"
         << "***ERR*** number of point wrong\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  return;
}

/*----------------------------------------------------------------------*
 |  get segment view (protected)                             mwgee 11/05|
 *----------------------------------------------------------------------*/
void MRTR::Overlap::SegmentView(vector<RefCountPtr<MRTR::Segment> >& segs)
{
  // allocate vector of ptrs
  segs.resize(Nseg());
  
  // get the segment views
  int count=0;
  map<int,RefCountPtr<MRTR::Segment> >::iterator curr;
  for (curr=s_.begin(); curr != s_.end(); ++curr)
  {
    segs[count] = curr->second;
    ++count;
  }
  if (count != Nseg())
  {
    cout << "***ERR*** MRTR::Overlap::SegmentView:\n"
         << "***ERR*** number of segments wrong\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  return;
}

/*----------------------------------------------------------------------*
 |  get point view (private)                                 mwgee 10/05|
 *----------------------------------------------------------------------*/
void MRTR::Overlap::PointView(map<int,RefCountPtr<MRTR::Point> >& p,
                                       vector<RefCountPtr<MRTR::Point> >& points)
{
  // allocate vector of ptrs
  int np = p.size();
  points.resize(np);
  
  // get the point views
  int count=0;
  map<int,RefCountPtr<MRTR::Point> >::iterator pcurr;
  for (pcurr=p.begin(); pcurr != p.end(); ++pcurr)
  {
    points[count] = pcurr->second;
    ++count;
  }
  if (count != np)
  {
    cout << "***ERR*** MRTR::Overlap::PointView:\n"
         << "***ERR*** number of point wrong\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  return;
}

/*----------------------------------------------------------------------*
 |  get point view (private)                                 mwgee 11/05|
 *----------------------------------------------------------------------*/
void MRTR::Overlap::PointView(vector<MRTR::Point*>& p,const int* nodeids,const int np)
{
  p.resize(np);
  
  for (int i=0; i<np; ++i)
  {
    map<int,RefCountPtr<MRTR::Point> >::iterator pcurr = p_.find(nodeids[i]);
    if (pcurr==p_.end())
    {
      cout << "***ERR*** MRTR::Overlap::PointView:\n"
           << "***ERR*** cannot find point " << nodeids[i] << "\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);
    }
    p[i] = pcurr->second.get();
  }
  return;
}

/*----------------------------------------------------------------------*
 |  perform clipping algorithm (private)                     mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MRTR::Overlap::Clipelements()
{
  if (!havemxi_ || !havesxi_ || !haveline_)
  {
    cout << "***ERR*** MRTR::Overlap::Clipelements:\n"
         << "***ERR*** initialization of Overlap class missing\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  
  // put all mseg nodes in polygon
  AddPointtoPolygon(100+0,&mline_[0][0]);
  AddPointtoPolygon(100+1,&mline_[1][0]);
  AddPointtoPolygon(100+2,&mline_[2][0]);
  
  //===========================================================================

  // loop edges of slave segment and clip the master edges
  // add all intersection points that can be found
  for (int clipedge=0; clipedge<3; ++clipedge)
  {
    // point on that clip edge (dim 2)
    double* PE = &sline_[clipedge][0];
    // the outward normal to the clip edge (dim 2)
    double* N = &sn_[clipedge][0];

    // loop edges of master segment and clip them against the clip edge
    // add all intersections
    for (int medge=0; medge<3; ++medge)
    {
      bool ok;
      
      // start point of medge with id 100+medge
      double* P0 = &mline_[medge][0];
      //int id0 = 100+medge;
      
      // end point of medge has id 100+(medge+1 < 3)
      double* P1 = &mline_[medge][2];
      //int id1 = medge+1;
      //if (id1==3) id1 = 100+0;
      //else        id1 += 100;
      
      // find intersection between medge and clipedge
      // if there's an intersection, put it in the polygon
      // the id is 10*clipedge+medge
      double xi[2];
      ok = Clip_Intersect(N,PE,P0,P1,xi);
      if (ok)
      {
        //cout << "OVERLAP Clipelements: adding intersection point " << 10*clipedge+medge << endl;
        AddPointtoPolygon(10*clipedge+medge,xi);
      }
    } // for (int medge=0; medge<3; ++medge)
  } // for (int clipedge=0; clipedge<3; ++clipedge)

  //===========================================================================

  // loop all clipedges and clip all points incl intersections 
  // that are in the polygon
  // throw away al point that are not inside
  {
    int           np    = SizePointPolygon();
    
    vector<RefCountPtr<MRTR::Point> > point;
    PointView(point);
    int p;
    for (p=0; p<np; ++p)
    {
      bool ok = true;
      //cout << "OVERLAP Clipelements: now clipping point " << point[p]->Id() << endl;
      const double* P = point[p]->Xi();
      for (int clipedge=0; clipedge<3; ++clipedge)
      {
        // point on that clip edge (dim 2)
        double* PE = &sline_[clipedge][0];
        // the outward normal to the clip edge (dim 2)
        double* N = &sn_[clipedge][0];
        
        // clip point P against this edge
        ok = Clip_TestPoint(N,PE,P);
        // leave point in
        if (ok) continue;
        // remove this point
        else
        {
          ok = false;
          break;
        }
      } // for (int clipedge=0; clipedge<3; ++clipedge)
      
      // if not found inside, remove point
      if (!ok)
      {
        //cout << "OVERLAP Clipelements: removing point " << point[p]->Id() << endl;
        RemovePointfromPolygon(point[p]->Id(),P);
      }
            
    } // for (int p=0; p<np; ++p)
    point.clear();
  }

  //===========================================================================

  // see whether there are points left in the polygon
  // if there are no points in the polygon, the sseg could still be completely inside
  // the mseg, now test for that case
  {
    int np = SizePointPolygon();
#if 0
    if (np) cout << "OVERLAP Clipelements: # point in polygon after clipping " << np << endl;
#endif
    if  (!np)
    {
      // we have to project 1 slave node onto the master element to make decision
      // that sseg is not completely inside mseg
      MRTR::Node** snode = sseg_.Nodes();
      MRTR::Projector projector(inter_.IsOneDimensional());
      double xi[2];
      projector.ProjectNodetoSegment_NodalNormal(*snode[0],mseg_,xi);
      if (xi[0]<=1. && xi[1]<=abs(1.-xi[0]) && xi[0]>=0. && xi[1]>=0.)
      {
        cout << "OVERLAP !!!!!! sseg is in mseg though clipping is zero!\n";
        exit(EXIT_FAILURE);
      }
      else 
        return false;
    }
  }

  //===========================================================================
  //===========================================================================
  //===========================================================================

#if 0
  // make printout of the polygon so far
  {
    int np    = SizePointPolygon();
    vector<RefCountPtr<MRTR::Point> > point; PointView(point);
    for (int p=0; p<np; ++p)
    {
      cout << "OVERLAP Clipelements: point " << setw(3) << point[p]->Id() 
           << " xi " << point[p]->Xi()[0] 
           << "/" << point[p]->Xi()[1] << endl;
    }
    point.clear();
  }
#endif

  //===========================================================================

  // count how many corner nodes of mseg are in and how many
  // intersections there are
  int np = SizePointPolygon();
  vector<RefCountPtr<MRTR::Point> > point; PointView(point);
  int nmcorner=0;
  int nminter=0;
  for (int p=0; p<np; ++p)
  {
    if (point[p]->Id() < 100) ++nminter;
    else                      ++nmcorner;
  }

  // maybe its a good idea to collapse intersections with nodes (that are in)
  // that are really close to each other

  // there are 12 cases watched out for
  // cases can be partially determined by looking at nminter and nmcorner
  MRTR::Overlap::Polygon_Type type = MRTR::Overlap::poly_none;

#if 0  
  cout << "OVERLAP Clipelements: corners in " << nmcorner << " intersections " << nminter << endl;
#endif

  //---------------------------------------------------------------------------
  // cases 1_2
  if (nmcorner==1 && nminter==2)
  {
    int pair[3]; pair[0] = pair[1] = pair[2] = 0;
    for (int p=0; p<np; ++p)
    {
      if (point[p]->Id()>=100) continue;
      int ten = MRTR::digit_ten(point[p]->Id());
      if (ten != 0 && ten != 1 && ten != 2)
      {
        cout << "***ERR*** MRTR::Overlap::Clipelements:\n"
             << "***ERR*** weird overflow\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        exit(EXIT_FAILURE);
      }
      pair[ten] += 1;
    }    
    //// case 1_2_same (1 corner in, 2 intersections same edge)
    for (int i=0; i<3; ++i)
      if (pair[i]==2)
      {
        type = MRTR::Overlap::poly_1_2_same;
        break;
      }
    //// case 1_2_dif (1 corner in, 2 intersections different edges)
    if (type == MRTR::Overlap::poly_none)
      type = MRTR::Overlap::poly_1_2_dif;
  }



  //---------------------------------------------------------------------------
  // cases 2_2
  else if (nmcorner==2 && nminter==2)
  {
    int pair[3]; pair[0] = pair[1] = pair[2] = 0;
    for (int p=0; p<np; ++p)
    {
      if (point[p]->Id()>=100) continue;
      int ten = MRTR::digit_ten(point[p]->Id());
      if (ten != 0 && ten != 1 && ten != 2)
      {
        cout << "***ERR*** MRTR::Overlap::Clipelements:\n"
             << "***ERR*** weird overflow\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        exit(EXIT_FAILURE);
      }
      pair[ten] += 1;
    }    
    //// case 2_2_same (2 corner in, 2 intersections same edge)
    for (int i=0; i<3; ++i)
      if (pair[i]==2)
      {
        type = MRTR::Overlap::poly_2_2_same;
        break;
      }
    
    //// case 2_2_dif (2 corner in, 2 intersections different edges)
    if (type == MRTR::Overlap::poly_none)
      type = MRTR::Overlap::poly_2_2_dif;
  }



  //---------------------------------------------------------------------------
  // cases 0_4
  else if (nmcorner==0 && nminter==4)
  {
    int pair[3]; pair[0] = pair[1] = pair[2] = 0;
    for (int p=0; p<np; ++p)
    {
      if (point[p]->Id()>=100) continue;
      int ten = MRTR::digit_ten(point[p]->Id());
      if (ten != 0 && ten != 1 && ten != 2)
      {
        cout << "***ERR*** MRTR::Overlap::Clipelements:\n"
             << "***ERR*** weird overflow\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        exit(EXIT_FAILURE);
      }
      pair[ten] += 1;
    }    
    //// case 0_4_same (pairs of 2 intersections on same edge)
    int count=0;
    for (int i=0; i<3; ++i)
      if (pair[i]==2) ++count;
    if (count==2)
      type = MRTR::Overlap::poly_0_4_same;
    
    //// case 0_4_dif (pair of 2 intersections on same edge, other 2 on different edges)
    else if (count==1)
      type = MRTR::Overlap::poly_0_4_dif;
      
    //// there are no pairs
    else
    {
      cout << "***ERR*** MRTR::Overlap::Clipelements:\n"
           << "***ERR*** unknown case of polygon\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);
    }
  }



  //---------------------------------------------------------------------------
  // case 1_4 (1 corner in, 4 intersections) 
  else if (nmcorner==1 && nminter==4)
  {
    int pair[3]; pair[0] = pair[1] = pair[2] = 0;
    for (int p=0; p<np; ++p)
    {
      if (point[p]->Id()>=100) continue;
      int ten = MRTR::digit_ten(point[p]->Id());
      if (ten != 0 && ten != 1 && ten != 2)
      {
        cout << "***ERR*** MRTR::Overlap::Clipelements:\n"
             << "***ERR*** weird overflow\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        exit(EXIT_FAILURE);
      }
      pair[ten] += 1;
    }    
    int count=0;
    for (int i=0; i<3; ++i)
      if (pair[i]==2) ++count;
    //// case 1_4_same (1 corner in, pairs of 2 intersections on same edge)
    if (count==2)
      type = MRTR::Overlap::poly_1_4_same;
    //// case 1_4_dif (1 corner in, 1 pair of intersections on same edge)
    else if (count==1) 
      type = MRTR::Overlap::poly_1_4_dif;
    //// there are no pairs
    else
    {
      cout << "***ERR*** MRTR::Overlap::Clipelements:\n"
           << "***ERR*** unknown case of polygon\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);
    }
  }



  //---------------------------------------------------------------------------
  // case 0_2 (2 intersections on different edges)
  else if (nmcorner==0 && nminter==2)
  {
    // intersections have to be on different edges
    int pair[3]; pair[0] = pair[1] = pair[2] = 0;
    for (int p=0; p<np; ++p)
    {
      if (point[p]->Id()>=100) continue;
      int ten = MRTR::digit_ten(point[p]->Id());
      if (ten != 0 && ten != 1 && ten != 2)
      {
        cout << "***ERR*** MRTR::Overlap::Clipelements:\n"
             << "***ERR*** weird overflow\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        exit(EXIT_FAILURE);
      }
      pair[ten] += 1;
    }    
    int count=0;
    for (int i=0; i<3; ++i)
      if (pair[i]==2) ++count;
    if (count!=0)
    {
      cout << "***ERR*** MRTR::Overlap::Clipelements:\n"
           << "***ERR*** unknown case of polygon\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);
    }
    else
      type = MRTR::Overlap::poly_0_2;
  }



  //---------------------------------------------------------------------------
  // case 0_6 (0 corners in, 6 intersections) 
  else if (nmcorner==0 && nminter==6)
    type = MRTR::Overlap::poly_0_6;



  //---------------------------------------------------------------------------
  // case 3_0 (3 corners in, mseg completely inside sseg) 
  else if (nmcorner==3 && nminter==0)
    type = MRTR::Overlap::poly_3_0;



  //---------------------------------------------------------------------------
  // case 0_0 (no corners in, sseg completely inside mseg)
  else if (nmcorner==0 && nminter==0)
    type = MRTR::Overlap::poly_0_0;



  //---------------------------------------------------------------------------
  // unknown case
  else
  {
    cout << "***ERR*** MRTR::Overlap::Clipelements:\n"
         << "***ERR*** unknown case of polygon\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }

  //===========================================================================

  // finish the polygon
  switch(type)
  {
  //---------------------------------------------------------------------------
  case poly_1_2_same:
    Clip_FixPolygon_1_2_same(np,&point[0]);
  break;  
  //---------------------------------------------------------------------------
  case poly_1_2_dif:
    Clip_FixPolygon_1_2_dif(np,&point[0]);
  break;  
  //---------------------------------------------------------------------------
  case poly_1_4_same:
    Clip_FixPolygon_1_4_same(np,&point[0]);
  break;  
  //---------------------------------------------------------------------------
  case poly_1_4_dif:
    Clip_FixPolygon_1_4_dif(np,&point[0]);
  break;  
  //---------------------------------------------------------------------------
  case poly_2_2_same:
    Clip_FixPolygon_2_2_same(np,&point[0]);
  break;  
  //---------------------------------------------------------------------------
  case poly_2_2_dif:
    Clip_FixPolygon_2_2_dif(np,&point[0]);
  break;  
  //---------------------------------------------------------------------------
  case poly_3_0:
    Clip_FixPolygon_3_0(np,&point[0]);
  break;  
  //---------------------------------------------------------------------------
  case poly_0_0:
    Clip_FixPolygon_0_0(np,&point[0]);
  break;  
  //---------------------------------------------------------------------------
  case poly_0_4_same:
    Clip_FixPolygon_0_4_same(np,&point[0]);
  break;  
  //---------------------------------------------------------------------------
  case poly_0_4_dif:
    Clip_FixPolygon_0_4_dif(np,&point[0]);
  break;  
  //---------------------------------------------------------------------------
  case poly_0_2:
    Clip_FixPolygon_0_2(np,&point[0]);
  break;  
  //---------------------------------------------------------------------------
  case poly_0_6:
    Clip_FixPolygon_0_6(np,&point[0]);
  break;  
  //---------------------------------------------------------------------------
  case MRTR::Overlap::poly_none:
    cout << "***ERR*** MRTR::Overlap::Clipelements:\n"
         << "***ERR*** unknown case of polygon\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  break;
  default:
    cout << "***ERR*** MRTR::Overlap::Clipelements:\n"
         << "***ERR*** unknown case of polygon\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  break;
  }

  point.clear();
  return true;
}


/*----------------------------------------------------------------------*
 |  fix polygon in the case 1_2_dif (private)                mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MRTR::Overlap::Clip_FixPolygon_1_2_dif(int np, RefCountPtr<MRTR::Point>* point)
{
  //cout << "OVERLAP This is case 1_2_dif\n";
  
  // identify the slave edges with the 2 intersections
  int sedge[2];
  int count=0;
  for (int p=0; p<np; ++p)
  {
    if (point[p]->Id() >= 100) continue;
    sedge[count] = MRTR::digit_ten(point[p]->Id());
    ++count;
  }
  
  // identify the slave node between them
  int prevedge = -1;
  int spoint = -1;
  int nextedge = -1;
  // order the edges
  if (sedge[0]>sedge[1])
    MRTR::swap(sedge[0],sedge[1]);
  if (sedge[0]==0 && sedge[1]==1) 
  {
    spoint = 1;
    prevedge = 0;
    nextedge = 1;
  }
  else if (sedge[0]==1 && sedge[1]==2)
  {
    spoint = 2;
    prevedge = 1;
    nextedge = 2;
  }
  else if (sedge[0]==0 && sedge[1]==2)
  {
    spoint=0;
    prevedge = 2;
    nextedge = 0;
  }
  else
  {
    cout << "***ERR*** MRTR::Overlap::Clip_FixPolygon_1_2_dif:\n"
         << "***ERR*** can't find slave point\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }

  // add the spoint to the polygon
  AddPointtoPolygon(1000+spoint,&sxi_[spoint][0]);
  
  // find the convex hull for all points in polygon
  ConvexHull(p_);

  
#if 0
  // printout the polygon
  map<int,MRTR::Point*>::iterator pcurr;
  for (pcurr=p_.begin(); pcurr != p_.end(); ++pcurr)
    if (pcurr->second)
      cout << "OVERLAP finished polygon point " << pcurr->second->Id()
           << " xi " << pcurr->second->Xi()[0] << " / " << pcurr->second->Xi()[1] << endl;
#endif  
  return true;
}


/*----------------------------------------------------------------------*
 |  fix polygon in the case 1_2_same (private)                mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MRTR::Overlap::Clip_FixPolygon_1_2_same(int np, RefCountPtr<MRTR::Point>* point)
{
  //cout << "OVERLAP This is case 1_2_same\n";

  // make sure there's two intersections and one internal
  RefCountPtr<MRTR::Point> sedgepoint[2];
  int count1=0;
  int count2=0;
  for (int p=0; p<np; ++p)
  {
    if (point[p]->Id() >= 100) 
    {
      ++count2;
      continue;
    }
    sedgepoint[count1] = point[p];
    ++count1;
  }
  if (count1 != 2 || count2 != 1)
  {
    cout << "***ERR*** MRTR::Overlap::Clip_FixPolygon_1_2_same:\n"
         << "***ERR*** number of points wrong\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }

  // find the convex hull for all points in polygon
  ConvexHull(p_);

  
#if 0
  // printout the polygon
  map<int,MRTR::Point*>::iterator pcurr;
  for (pcurr=p_.begin(); pcurr != p_.end(); ++pcurr)
    if (pcurr->second)
      cout << "OVERLAP finished polygon point " << pcurr->second->Id()
           << " xi " << pcurr->second->Xi()[0] << " / " << pcurr->second->Xi()[1] << endl;
#endif  
  return true;
}

/*----------------------------------------------------------------------*
 |  fix polygon in the case 2_2_same (private)               mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MRTR::Overlap::Clip_FixPolygon_2_2_same(int np, RefCountPtr<MRTR::Point>* point)
{
  //cout << "OVERLAP This is case 2_2_same\n";

  // identify the slave edge
  int sedge;
  RefCountPtr<MRTR::Point> sedgepoint[2];
  RefCountPtr<MRTR::Point> inpoint[2];
  int count1=0;
  int count2=0;
  for (int p=0; p<np; ++p)
  {
    if (point[p]->Id() >= 100) 
    {
      inpoint[count2] = point[p];
      ++count2;
    }
    else
    {
      sedge = MRTR::digit_ten(point[p]->Id());
      sedgepoint[count1] = point[p];
      ++count1;
    }
  }
  if (count1 != 2 || count2 !=2)
  {
    cout << "***ERR*** MRTR::Overlap::Clip_FixPolygon_2_2_same:\n"
         << "***ERR*** number of intersections or number of inner nodes wrong\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  
  // find the convex hull for all points in polygon
  ConvexHull(p_);
  
  
#if 0
  // printout the polygon
  map<int,MRTR::Point*>::iterator pcurr;
  for (pcurr=p_.begin(); pcurr != p_.end(); ++pcurr)
    if (pcurr->second)
      cout << "OVERLAP finished polygon point " << pcurr->second->Id()
           << " xi " << pcurr->second->Xi()[0] << " / " << pcurr->second->Xi()[1] << endl;
#endif  
  return true;
}


/*----------------------------------------------------------------------*
 |  fix polygon in the case 2_2_dif (private)                mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MRTR::Overlap::Clip_FixPolygon_2_2_dif(int np, RefCountPtr<MRTR::Point>* point)
{
  //cout << "OVERLAP This is case 2_2_dif\n";

  // identify the slave edge
  int sedge[2];
  RefCountPtr<MRTR::Point> sedgepoint[2];
  RefCountPtr<MRTR::Point> inpoint[2];
  int count1=0;
  int count2=0;
  for (int p=0; p<np; ++p)
  {
    if (point[p]->Id() >= 100) 
    {
      inpoint[count2] = point[p];
      ++count2;
    }
    else
    {
      sedge[count1] = MRTR::digit_ten(point[p]->Id());
      sedgepoint[count1] = point[p];
      ++count1;
    }
  }
  if (count1 != 2 || count2 !=2)
  {
    cout << "***ERR*** MRTR::Overlap::Clip_FixPolygon_2_2_dif:\n"
         << "***ERR*** number of intersections or number of inner nodes wrong\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }

  // identify the slave node between them
  int spoint = -1;
  // order the edges
  if (sedge[0]>sedge[1])
    MRTR::swap(sedge[0],sedge[1]);
  if (sedge[0]==0 && sedge[1]==1) 
    spoint = 1;
  else if (sedge[0]==1 && sedge[1]==2)
    spoint = 2;
  else if (sedge[0]==0 && sedge[1]==2)
    spoint=0;
  else
  {
    cout << "***ERR*** MRTR::Overlap::Clip_FixPolygon_2_2_dif:\n"
         << "***ERR*** can't find slave point\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }

  // add the slave node to the polygon
  AddPointtoPolygon(1000+spoint,&sxi_[spoint][0]);  
  
  // find the convex hull for all points in polygon
  ConvexHull(p_);




#if 0
  // printout the polygon
  map<int,MRTR::Point*>::iterator pcurr;
  for (pcurr=p_.begin(); pcurr != p_.end(); ++pcurr)
    if (pcurr->second)
      cout << "OVERLAP finished polygon point " << pcurr->second->Id()
           << " xi " << pcurr->second->Xi()[0] << " / " << pcurr->second->Xi()[1] << endl;
#endif  
  return true;
}


/*----------------------------------------------------------------------*
 |  fix polygon in the case 1_4_same (private)               mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MRTR::Overlap::Clip_FixPolygon_1_4_same(int np, RefCountPtr<MRTR::Point>* point)
{
  //cout << "OVERLAP This is case 1_4_same\n";

  // identify the slave edges
  int sedge[4];
  RefCountPtr<MRTR::Point> sedgepoint[4];
  RefCountPtr<MRTR::Point> inpoint;
  int count1=0;
  int count2=0;
  for (int p=0; p<np; ++p)
  {
    if (point[p]->Id() >= 100) 
    {
      inpoint = point[p];
      ++count2;
    }
    else
    {
      sedge[count1] = MRTR::digit_ten(point[p]->Id());
      sedgepoint[count1] = point[p];
      ++count1;
    }
  }
  if (count1 != 4 || count2 !=1)
  {
    cout << "***ERR*** MRTR::Overlap::Clip_FixPolygon_1_4_same:\n"
         << "***ERR*** number of intersections or number of inner nodes wrong\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  
  // find the convex hull for all points in polygon
  ConvexHull(p_);

  
#if 0
  // printout the polygon
  map<int,MRTR::Point*>::iterator pcurr;
  for (pcurr=p_.begin(); pcurr != p_.end(); ++pcurr)
    if (pcurr->second)
      cout << "OVERLAP finished polygon point " << pcurr->second->Id()
           << " xi " << pcurr->second->Xi()[0] << " / " << pcurr->second->Xi()[1] << endl;
#endif  
  return true;
}


/*----------------------------------------------------------------------*
 |  fix polygon in the case 1_4_dif (private)                mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MRTR::Overlap::Clip_FixPolygon_1_4_dif(int np, RefCountPtr<MRTR::Point>* point)
{
  //cout << "OVERLAP This is case 1_4_dif\n";

  // identify the slave edges
  int sedge[4];
  RefCountPtr<MRTR::Point> sedgepoint[4];
  RefCountPtr<MRTR::Point> inpoint;
  int count1=0;
  int count2=0;
  for (int p=0; p<np; ++p)
  {
    if (point[p]->Id() >= 100) 
    {
      inpoint = point[p];
      ++count2;
    }
    else
    {
      sedge[count1] = MRTR::digit_ten(point[p]->Id());
      sedgepoint[count1] = point[p];
      ++count1;
    }
  }
  if (count1 != 4 || count2 !=1)
  {
    cout << "***ERR*** MRTR::Overlap::Clip_FixPolygon_1_4_dif:\n"
         << "***ERR*** number of intersections or number of inner nodes wrong\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  
  // find the 2 intersection points that are 'lonely' on an edge
  int lonelies[2];
  count1=0;
  for (int i=0; i<4; ++i)
  {
    bool foundagain=false;
    for (int j=0; j<4; ++j)
    {
      if (j==i) continue;
      if (sedge[j]==sedge[i])
      {
        foundagain = true;
        break;
      }
    }
    if (!foundagain)
    {
      lonelies[count1] = i;
      ++count1;
    }
  }
  if (count1 != 2)
  {
    cout << "***ERR*** MRTR::Overlap::Clip_FixPolygon_1_4_dif:\n"
         << "***ERR*** overflow appeared\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  
  // identify the slave node between the 2 lonely points
  int spoint = -1;
  // get the lonely edges
  int e[2];
  e[0] = sedge[lonelies[0]];
  e[1] = sedge[lonelies[1]];
  // order the edges
  if (e[0]>e[1])
    MRTR::swap(e[0],e[1]);
  if (e[0]==0 && e[1]==1) 
    spoint = 1;
  else if (e[0]==1 && e[1]==2)
    spoint = 2;
  else if (e[0]==0 && e[1]==2)
    spoint=0;
  else
  {
    cout << "***ERR*** MRTR::Overlap::Clip_FixPolygon_1_4_dif:\n"
         << "***ERR*** can't find slave point\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  
  // put the slave point in
  AddPointtoPolygon(1000+spoint,&sxi_[spoint][0]);  

  // find the convex hull for all points in polygon
  ConvexHull(p_);
  
#if 0
  // printout the polygon
  map<int,MRTR::Point*>::iterator pcurr;
  for (pcurr=p_.begin(); pcurr != p_.end(); ++pcurr)
    if (pcurr->second)
      cout << "OVERLAP finished polygon point " << pcurr->second->Id()
           << " xi " << pcurr->second->Xi()[0] << " / " << pcurr->second->Xi()[1] << endl;
#endif  
  return true;
}

/*----------------------------------------------------------------------*
 |  fix polygon in the case 3_0 (private)                    mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MRTR::Overlap::Clip_FixPolygon_3_0(int np, RefCountPtr<MRTR::Point>* point)
{
  //cout << "OVERLAP This is case 3_0\n";

  // make sure there's no intersections and three internal
  int count=0;
  for (int p=0; p<np; ++p)
  {
    if (point[p]->Id() >= 100) 
    {
      ++count;
      continue;
    }
  }
  if (count != 3 )
  {
    cout << "***ERR*** MRTR::Overlap::Clip_FixPolygon_3_0:\n"
         << "***ERR*** number of points wrong\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }

  // find the convex hull for all points in polygon
  ConvexHull(p_);

  
#if 0
  // printout the polygon
  map<int,MRTR::Point*>::iterator pcurr;
  for (pcurr=p_.begin(); pcurr != p_.end(); ++pcurr)
    if (pcurr->second)
      cout << "OVERLAP finished polygon point " << pcurr->second->Id()
           << " xi " << pcurr->second->Xi()[0] << " / " << pcurr->second->Xi()[1] << endl;
#endif  
  return true;
}

/*----------------------------------------------------------------------*
 |  fix polygon in the case 0_0 (private)                    mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MRTR::Overlap::Clip_FixPolygon_0_0(int np, RefCountPtr<MRTR::Point>* point)
{
  //cout << "OVERLAP This is case 0_0\n";

  // make sure there's no points in polygon
  if (np)
  {
    cout << "***ERR*** MRTR::Overlap::Clip_FixPolygon_0_0:\n"
         << "***ERR*** number of points wrong\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  
  
  // all slave points are in the master segment, add them
  AddPointtoPolygon(1000+0,&sxi_[0][0]);  
  AddPointtoPolygon(1000+1,&sxi_[1][0]);  
  AddPointtoPolygon(1000+2,&sxi_[2][0]);  
  
  
  // find the convex hull for all points in polygon
  ConvexHull(p_);

  
#if 0
  // printout the polygon
  map<int,MRTR::Point*>::iterator pcurr;
  for (pcurr=p_.begin(); pcurr != p_.end(); ++pcurr)
    if (pcurr->second)
      cout << "OVERLAP finished polygon point " << pcurr->second->Id()
           << " xi " << pcurr->second->Xi()[0] << " / " << pcurr->second->Xi()[1] << endl;
#endif  
  return true;
}

/*----------------------------------------------------------------------*
 |  fix polygon in the case 0_4_same (private)                mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MRTR::Overlap::Clip_FixPolygon_0_4_same(int np, RefCountPtr<MRTR::Point>* point)
{
  //cout << "OVERLAP This is case 0_4_same\n";

  // make sure there's two intersections and one internal
  int count1=0;
  int count2=0;
  for (int p=0; p<np; ++p)
  {
    if (point[p]->Id() >= 100) 
    {
      ++count2;
      continue;
    }
    else ++count1;
  }
  if (count1 != 4 || count2 != 0)
  {
    cout << "***ERR*** MRTR::Overlap::Clip_FixPolygon_0_4_same:\n"
         << "***ERR*** number of points wrong\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }

  // find the convex hull for all points in polygon
  ConvexHull(p_);

  
#if 0
  // printout the polygon
  map<int,MRTR::Point*>::iterator pcurr;
  for (pcurr=p_.begin(); pcurr != p_.end(); ++pcurr)
    if (pcurr->second)
      cout << "OVERLAP finished polygon point " << pcurr->second->Id()
           << " xi " << pcurr->second->Xi()[0] << " / " << pcurr->second->Xi()[1] << endl;
#endif  
  return true;
}

/*----------------------------------------------------------------------*
 |  fix polygon in the case 0_4_dif (private)                mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MRTR::Overlap::Clip_FixPolygon_0_4_dif(int np, RefCountPtr<MRTR::Point>* point)
{
  //cout << "OVERLAP This is case 0_4_dif\n";

  // identify the slave edges
  int sedge[4];
  RefCountPtr<MRTR::Point> sedgepoint[4];
  int count1=0;
  int count2=0;
  for (int p=0; p<np; ++p)
  {
    if (point[p]->Id() >= 100) 
      ++count2;
    else
    {
      sedge[count1] = MRTR::digit_ten(point[p]->Id());
      sedgepoint[count1] = point[p];
      ++count1;
    }
  }
  if (count1 != 4)
  {
    cout << "***ERR*** MRTR::Overlap::Clip_FixPolygon_0_4_dif:\n"
         << "***ERR*** number of intersections wrong\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  
  // find the 2 intersection points that are 'lonely' on an edge
  int lonelies[2];
  count1=0;
  for (int i=0; i<4; ++i)
  {
    bool foundagain=false;
    for (int j=0; j<4; ++j)
    {
      if (j==i) continue;
      if (sedge[j]==sedge[i])
      {
        foundagain = true;
        break;
      }
    }
    if (!foundagain)
    {
      lonelies[count1] = i;
      ++count1;
    }
  }
  if (count1 != 2)
  {
    cout << "***ERR*** MRTR::Overlap::Clip_FixPolygon_0_4_dif:\n"
         << "***ERR*** overflow appeared\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  
  // identify the slave node between the 2 lonely points
  int spoint = -1;
  // get the lonely edges
  int e[2];
  e[0] = sedge[lonelies[0]];
  e[1] = sedge[lonelies[1]];
  // order the edges
  if (e[0]>e[1])
    MRTR::swap(e[0],e[1]);
  if (e[0]==0 && e[1]==1) 
    spoint = 1;
  else if (e[0]==1 && e[1]==2)
    spoint = 2;
  else if (e[0]==0 && e[1]==2)
    spoint=0;
  else
  {
    cout << "***ERR*** MRTR::Overlap::Clip_FixPolygon_0_4_dif:\n"
         << "***ERR*** can't find slave point\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  
  // put the slave point in
  AddPointtoPolygon(1000+spoint,&sxi_[spoint][0]);  

  // find the convex hull for all points in polygon
  ConvexHull(p_);

  
#if 0
  // printout the polygon
  map<int,MRTR::Point*>::iterator pcurr;
  for (pcurr=p_.begin(); pcurr != p_.end(); ++pcurr)
    if (pcurr->second)
      cout << "OVERLAP finished polygon point " << pcurr->second->Id()
           << " xi " << pcurr->second->Xi()[0] << " / " << pcurr->second->Xi()[1] << endl;
#endif  
  return true;
}

/*----------------------------------------------------------------------*
 |  fix polygon in the case 0_2 (private)                    mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MRTR::Overlap::Clip_FixPolygon_0_2(int np, RefCountPtr<MRTR::Point>* point)
{
  //cout << "OVERLAP This is case 0_2\n";

  // identify the slave edges
  int sedge[2];
  RefCountPtr<MRTR::Point> sedgepoint[2];
  int count1=0;
  int count2=0;
  for (int p=0; p<np; ++p)
  {
    if (point[p]->Id() >= 100) 
      ++count2;
    else
    {
      sedge[count1] = MRTR::digit_ten(point[p]->Id());
      sedgepoint[count1] = point[p];
      ++count1;
    }
  }
  if (count1 != 2)
  {
    cout << "***ERR*** MRTR::Overlap::Clip_FixPolygon_0_2:\n"
         << "***ERR*** number of intersections wrong\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  
  // the 2 intersections have to be on different edges
  if (sedge[0]==sedge[1])
  {
    cout << "***ERR*** MRTR::Overlap::Clip_FixPolygon_0_2:\n"
         << "***ERR*** weird: 2 intersecions on same edge and nothing else!\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  
    // identify the slave node between the 2 edge
  int spoint = -1;
  // order the edges
  if (sedge[0]>sedge[1])
    MRTR::swap(sedge[0],sedge[1]);
  if (sedge[0]==0 && sedge[1]==1) 
    spoint = 1;
  else if (sedge[0]==1 && sedge[1]==2)
    spoint = 2;
  else if (sedge[0]==0 && sedge[1]==2)
    spoint=0;
  else
  {
    cout << "***ERR*** MRTR::Overlap::Clip_FixPolygon_0_2:\n"
         << "***ERR*** can't find slave point\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  
  // the polygon can now either contain spoint or the other 2 edges of sseg
  // the only way to find out is to project spoint onto mseg
  // if its in : put it in
  //        out: put the other 2 in
  MRTR::Node** snode = sseg_.Nodes();
  MRTR::Projector projector(inter_.IsOneDimensional());
  double xi[2];
  bool isin = false;
  projector.ProjectNodetoSegment_NodalNormal(*snode[spoint],mseg_,xi);
  if (xi[0]<=1. && xi[1]<=abs(1.-xi[0]) && xi[0]>=0. && xi[1]>=0.)
    isin = true;

  if (isin)
  {
    AddPointtoPolygon(1000+spoint,&sxi_[spoint][0]);  
  }
  else
  {
    if (spoint==0)
    {
      AddPointtoPolygon(1000+1,&sxi_[1][0]);  
      AddPointtoPolygon(1000+2,&sxi_[2][0]);  
    }
    else if (spoint==1)
    {
      AddPointtoPolygon(1000+0,&sxi_[0][0]);  
      AddPointtoPolygon(1000+2,&sxi_[2][0]);  
    }
    else if (spoint==2)
    {
      AddPointtoPolygon(1000+0,&sxi_[0][0]);  
      AddPointtoPolygon(1000+1,&sxi_[1][0]);  
    }
    else
    {
      cout << "***ERR*** MRTR::Overlap::Clip_FixPolygon_0_2:\n"
           << "***ERR*** point " << spoint << " out of scope (0,1,2)\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);
    }
  }  
  


  // find the convex hull for all points in polygon
  ConvexHull(p_);

  
#if 0
  // printout the polygon
  map<int,MRTR::Point*>::iterator pcurr;
  for (pcurr=p_.begin(); pcurr != p_.end(); ++pcurr)
    if (pcurr->second)
      cout << "OVERLAP finished polygon point " << pcurr->second->Id()
           << " xi " << pcurr->second->Xi()[0] << " / " << pcurr->second->Xi()[1] << endl;
#endif  
  return true;
}

/*----------------------------------------------------------------------*
 |  fix polygon in the case 0_6 (private)                    mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MRTR::Overlap::Clip_FixPolygon_0_6(int np, RefCountPtr<MRTR::Point>* point)
{
  //cout << "OVERLAP This is case 0_6\n";

  // make sure there's 6 intersections and no internal
  int count1=0;
  int count2=0;
  for (int p=0; p<np; ++p)
  {
    if (point[p]->Id() >= 100) 
    {
      ++count2;
      continue;
    }
    else ++count1;
  }
  if (count1 != 6 || count2 != 0)
  {
    cout << "***ERR*** MRTR::Overlap::Clip_FixPolygon_0_6:\n"
         << "***ERR*** number of points wrong\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }

  // find the convex hull for all points in polygon
  ConvexHull(p_);

  
#if 0
  // printout the polygon
  map<int,MRTR::Point*>::iterator pcurr;
  for (pcurr=p_.begin(); pcurr != p_.end(); ++pcurr)
    if (pcurr->second)
      cout << "OVERLAP finished polygon point " << pcurr->second->Id()
           << " xi " << pcurr->second->Xi()[0] << " / " << pcurr->second->Xi()[1] << endl;
#endif  
  return true;
}

/*----------------------------------------------------------------------*
 |  make a triangulization of a polygon (private)             mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MRTR::Overlap::Triangulization()
{
  // we have a polygon that is in clockwise order at this moment

  // if more than 3 points, find a center point
  int np = SizePointPolygon();
  if (np<3)
  {
    cout << "***ERR*** MRTR::Overlap::Triangulization:\n"
         << "***ERR*** # point in polygon < 3 ... very strange!!!\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  if (np>3) // we have to add a center point
  {
    double xi[2]; xi[0] = xi[1] = 0.0;
    vector<RefCountPtr<MRTR::Point> > points; PointView(points);
    for (int i=0; i<np; ++i)
    {
      const double* pxi = points[i]->Xi();
      xi[0] += pxi[0];
      xi[1] += pxi[1];
    }
    xi[0] /= np;
    xi[1] /= np;
    // create a point -1 as center point and add it to the polygon
    AddPointtoPolygon(-1,xi);
    points.clear();
  } // if (np>3)
  
  np = SizePointPolygon();
  vector<RefCountPtr<MRTR::Point> >  points; PointView(points);

  // create a MRTR::Node for every point
  int dof[3]; dof[0] = dof[1] = dof[2] = -1;
  
  // find real world coords for all points
  // find real world normal for all points
  // note that the polygon is in slave segment parameter space and is
  // completely contained in the slave segment. we can therefore use
  // slave segment values to interpolate polgon point values
  for (int i=0; i<np; ++i)
  {
    double x[3]; x[0] = x[1] = x[2] = 0.0;
    double n[3]; n[0] = n[1] = n[2] = 0.0;
    double val[sseg_.Nnode()];
    sseg_.EvaluateFunction(0,points[i]->Xi(),val,sseg_.Nnode(),NULL);
    MRTR::Node** snodes = sseg_.Nodes();
    for (int j=0; j<sseg_.Nnode(); ++j)
      for (int k=0; k<3; ++k)
      {
        x[k] += val[j]*snodes[j]->X()[k];
        n[k] += val[j]*snodes[j]->N()[k];
      }
    double length = sqrt(MRTR::dot(n,n,3));
    for (int j=0; j<3; ++j) n[j] /= length;
    // create a node with this coords and normal;
    MRTR::Node* node = new MRTR::Node(points[i]->Id(),x,3,dof);
    node->SetN(n);
    // set node in point
    points[i]->SetNode(node);

#if 0
    cout << *points[i];
#endif

  }  

  // find projection values for all points in polygon on mseg
  {
    double mxi[2];
    MRTR::Projector projector(inter_.IsOneDimensional());
    for (int i=0; i<np; ++i)
    {
      RefCountPtr<MRTR::Node> node = points[i]->Node();
      projector.ProjectNodetoSegment_NodalNormal(*node,mseg_,mxi);
      // create a projected node and set it in node
      MRTR::ProjectedNode* pnode = new MRTR::ProjectedNode(*node,mxi,&mseg_);
      node->SetProjectedNode(pnode);
#if 0
      if (mxi[0]<=1. && mxi[1]<=abs(1.-mxi[0]) && mxi[0]>=0. && mxi[1]>=0.)
      cout << "OVERLAP: point " << points[i]->Id() << " is in mseg, mxi " << mxi[0] << " / " << mxi[1] << endl;
      else
      cout << "OVERLAP: point " << points[i]->Id() << " is NOT in mseg, mxi " << mxi[0] << " / " << mxi[1] << endl;
#endif

    }
  }

  // every point has 3 values of the 3 shape functions of a triangle:
  // 3 values from function 0 from sseg
  // 3 values from function 1 from sseg
  // 3 values from function 0 from mseg
  for (int i=0; i<np; ++i)
  {
    double val[3];
    // evaluate function 0 from sseg
    sseg_.EvaluateFunction(0,points[i]->Xi(),val,3,NULL);
    points[i]->StoreFunctionValues(0,val,3);
    // evaluate function 1 from sseg
    sseg_.EvaluateFunction(1,points[i]->Xi(),val,3,NULL);
    points[i]->StoreFunctionValues(1,val,3);
    // evaluate function 0 from mseg
    mseg_.EvaluateFunction(0,points[i]->Node()->GetProjectedNode()->Xi(),val,3,NULL);
    points[i]->StoreFunctionValues(2,val,3);
  }
  
  // create the triangle elements from polygon with centerpoint
  // In case of np==3, there is no centerpoint
  if (np>3)
  {
    // the polygon is in clockwise order, center point is points[0] and has
    // id = -1
    if (points[0]->Id() != -1)
    {
      cout << "***ERR*** MRTR::Overlap::Triangulization:\n"
           << "***ERR*** points[0]->Id() is not -1\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);
    }
    int nodeid[3];
    MRTR::Segment_BiLinearTri* tmp;
    RefCountPtr<MRTR::Function_LinearTri> func = rcp(new Function_LinearTri());
    for (int i=2; i<np; ++i)
    {
      // there are np-1 triangles
      // triangle ids go from 0 to np-2
      // a triangle is defined by nodes 0, i, i-1
      nodeid[0]   = points[0]->Id();
      nodeid[1]   = points[i]->Id();
      nodeid[2] = points[i-1]->Id();
      tmp = new MRTR::Segment_BiLinearTri(i-2,3,nodeid);
      // set a linear shape function to this triangle
      tmp->SetFunction(0,func);
      // add triangle the *this class
      AddSegment(tmp->Id(),tmp);
    }
    // add the last triangle defined by nodes 0, 1, np-1 separately
    nodeid[0]   = points[0]->Id();
    nodeid[1]   = points[1]->Id();
    nodeid[2] = points[np-1]->Id();
    tmp = new MRTR::Segment_BiLinearTri(np-2,3,nodeid);
    // set a linear shape function to this triangle
    tmp->SetFunction(0,func);
    // add triangle to the *this class
    AddSegment(tmp->Id(),tmp);
  }
  else if (np==3) // single triangle without centerpoint
  {
    int nodeid[3];
    MRTR::Segment_BiLinearTri* tmp;
    RefCountPtr<MRTR::Function_LinearTri> func = rcp(new Function_LinearTri());
    nodeid[0] = points[0]->Id();
    nodeid[1] = points[1]->Id();
    nodeid[2] = points[2]->Id();
    tmp = new MRTR::Segment_BiLinearTri(0,3,nodeid);
    // set a linear shape function to this triangle
    tmp->SetFunction(0,func);
    // add triangle the *this class
    AddSegment(tmp->Id(),tmp);
  }
  else
  {
    cout << "***ERR*** MRTR::Overlap::Triangulization:\n"
         << "***ERR*** # point in polygon < 3 ... very strange!!!\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  
  // create ptr topology between triangle and nodes in points
  vector<MRTR::Node*> nodes(np);
  for (int i=0; i<np; ++i) nodes[i] = points[i]->Node().get();
  
  // loop segments and set ptr to nodes in them
  map<int,RefCountPtr<MRTR::Segment> >::iterator curr;
  for (curr=s_.begin(); curr != s_.end(); ++curr)
    curr->second->GetPtrstoNodes(nodes);
  
  nodes.clear();
  
  points.clear();
  
  return true;
}

/*----------------------------------------------------------------------*
 |  perform a quick search (private)                         mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MRTR::Overlap::QuickOverlapTest()
{
  // we need the projection of the master element on the slave element 
  // to do this test, otherwise we assmue we passed it
  if (!havemxi_) 
    return true;
  
  // find the max and min coords of slave and master points
  double mmax[2]; 
  mmax[0] = mxi_[0][0];
  mmax[1] = mxi_[0][1];
  double mmin[2]; 
  mmin[0] = mxi_[0][0];
  mmin[1] = mxi_[0][1];
  for (int point=1; point<3; ++point)
    for (int dim=0; dim<2; ++dim)
    {
      if (mmax[dim] < mxi_[point][dim]) 
        mmax[dim] = mxi_[point][dim];
      if (mmin[dim] > mxi_[point][dim]) 
        mmin[dim] = mxi_[point][dim];
    }


  bool isin = false;
  // check xi0 direction
  if (  (0.0 < mmin[0] && mmin[0] < 1.0) ||
        (0.0 < mmax[0] && mmax[0] < 1.0) ||
        (mmin[0] < 0.0 && 1.0 < mmax[0])  )
  isin = true;
  
  if (!isin) 
    return false;

  // check xi1 direction
  if (  (0.0 < mmin[1] && mmin[1] < 1.0) ||
        (0.0 < mmax[1] && mmax[1] < 1.0) ||
        (mmin[1] < 0.0 && 1.0 < mmax[1])  )
  isin = true;
  
  if (!isin)
    return false;

  return true;
}










































#endif // TRILINOS_PACKAGE
