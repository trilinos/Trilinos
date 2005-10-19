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
  e_.clear();
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 10/05|
 *----------------------------------------------------------------------*/
MRTR::Overlap::~Overlap()
{
  map<int,MRTR::Point*>::iterator pcurr;
  for (pcurr=p_.begin(); pcurr != p_.end(); ++pcurr)
    if (pcurr->second)
      delete pcurr->second;
  p_.clear();

  map<int,MRTR::Edge*>::iterator ecurr;
  for (ecurr=e_.begin(); ecurr != e_.end(); ++ecurr)
    if (ecurr->second)
      delete ecurr->second;
  e_.clear();
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
    // check whether i is inside sseg
    if (mxi_[i][0]<=1. && mxi_[i][1]<=abs(1.-mxi_[i][0]) && mxi_[i][0]>=0. && mxi_[i][1]>=0.)
    {
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

  // project slave nodes onto master segment if not already done
  if (!havesxi_)
  {
    havesxi_ = build_sxi();
               build_normal();
  }
  
  // compute information about the edges of the slave/master triangle
  if (!haveline_)
    haveline_ = build_lines();

#if 1
  // debugging output
  cout << "OVERLAP: Slave  " << sseg_;
  cout << "OVERLAP: Master " << mseg_;
  cout << "node 0: " << mxi_[0][0] << "/" << mxi_[0][1] << endl
       << "     1: " << mxi_[1][0] << "/" << mxi_[1][1] << endl
       << "     2: " << mxi_[2][0] << "/" << mxi_[2][1] << endl;
#endif

  
  // perform clipping algorithm
  bool ok = Clipelements();
  if (!ok)
    return false;
  

  return false;
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
  Clip_AddPointtoPolygon(100+0,&mline_[0][0]);
  Clip_AddPointtoPolygon(100+1,&mline_[1][0]);
  Clip_AddPointtoPolygon(100+2,&mline_[2][0]);
  
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
        Clip_AddPointtoPolygon(10*clipedge+medge,xi);
      }
    } // for (int medge=0; medge<3; ++medge)
  } // for (int clipedge=0; clipedge<3; ++clipedge)

  //===========================================================================

  // loop all clipedges and clip all points incl intersections 
  // that are in the polygon
  // throw away al point that are not inside
  {
    int           np    = Clip_SizePolygon();
    MRTR::Point** point = Clip_PointView();
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
        Clip_RemovePointfromPolygon(point[p]->Id(),P);
      }
            
    } // for (int p=0; p<np; ++p)
    delete [] point; point = NULL;
  }

  //===========================================================================

  // see whether there are points left in the polygon
  // if there are no points in the polygon, the sseg could still be completely inside
  // the mseg, now test for that case
  {
    int np = Clip_SizePolygon();
    if (np) cout << "OVERLAP Clipelements: # point in polygon after clipping " << np << endl;
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

#if 1
  // make printout of the polygon so far
  {
    int np    = Clip_SizePolygon();
    MRTR::Point** point = Clip_PointView();
    for (int p=0; p<np; ++p)
    {
      cout << "OVERLAP Clipelements: point " << setw(3) << point[p]->Id() 
           << " xi " << point[p]->Xi()[0] 
           << "/" << point[p]->Xi()[1] << endl;
    }
    delete [] point; point = NULL;
  }
#endif

  //===========================================================================

  // count how many corner nodes of mseg are in and how many
  // intersections there are
  int np    = Clip_SizePolygon();
  MRTR::Point** point = Clip_PointView();
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
  
  cout << "OVERLAP Clipelements: corners in " << nmcorner << " intersections " << nminter << endl;

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
  break;  
  //---------------------------------------------------------------------------
  case poly_1_2_dif:
  break;  
  //---------------------------------------------------------------------------
  case poly_1_4_same:
  break;  
  //---------------------------------------------------------------------------
  case poly_1_4_dif:
  break;  
  //---------------------------------------------------------------------------
  case poly_2_2_same:
  break;  
  //---------------------------------------------------------------------------
  case poly_2_2_dif:
  break;  
  //---------------------------------------------------------------------------
  case poly_3_0:
  break;  
  //---------------------------------------------------------------------------
  case poly_0_0:
  break;  
  //---------------------------------------------------------------------------
  case poly_0_4_same:
  break;  
  //---------------------------------------------------------------------------
  case poly_0_4_dif:
  break;  
  //---------------------------------------------------------------------------
  case poly_0_2:
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




  if (point) delete point; point = NULL;
  return true;
}


/*----------------------------------------------------------------------*
 |  find intersection (private)                              mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MRTR::Overlap::Clip_Intersect(double* N,double* PE,double* P0,double* P1,double* xi)
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
 |  add point (private)                                      mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MRTR::Overlap::Clip_AddPointtoPolygon(int id, double* P)
{
  // check whether this point is already in there
  map<int,MRTR::Point*>::iterator curr = p_.find(id);
  // it's there
  if (curr != p_.end())
    curr->second->SetXi(P);
  else
  {
    //cout << "OVERLAP Clip_AddPointtoPolygon: added point " << id << endl;
    MRTR::Point* p = new MRTR::Point(id,P);
    p_.insert(pair<int,MRTR::Point*>(id,p));
  }
  return true;
}

/*----------------------------------------------------------------------*
 |  remove point (private)                                      mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MRTR::Overlap::Clip_RemovePointfromPolygon(const int id,const double* P)
{
  // check whether this point is in there
  MRTR::Point* p;
  map<int,MRTR::Point*>::iterator curr = p_.find(id);
  if (curr != p_.end())
  {
    if (curr->second) delete curr->second;
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
 |  add edge (private)                                       mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MRTR::Overlap::Clip_AddEdgetoPolygon(int id, int node0, int node1)
{
  // check whether this edge is already in there
  map<int,MRTR::Edge*>::iterator curr = e_.find(id);
  // it's there
  if (curr != e_.end())
    curr->second->SetNodes(node0,node1);
  else
  {
    cout << "OVERLAP Clip_AddEdgetoPolygon: added edge " << id << " between nodes " << node0 << " - " << node1 << endl;
    MRTR::Edge* e = new MRTR::Edge(id,node0,node1);
    e_.insert(pair<int,MRTR::Edge*>(id,e));
  }
  return true;
}

/*----------------------------------------------------------------------*
 |  get edge view (private)                                  mwgee 10/05|
 *----------------------------------------------------------------------*/
MRTR::Edge* MRTR::Overlap::GetEdgeViewfromPolygon(int id)
{
  // check whether this edge is already in there
  map<int,MRTR::Edge*>::iterator curr = e_.find(id);
  // it's there
  if (curr != e_.end())
    return curr->second;
  else
    return NULL;
}

/*----------------------------------------------------------------------*
 |  get point view (private)                                 mwgee 10/05|
 *----------------------------------------------------------------------*/
MRTR::Point** MRTR::Overlap::Clip_PointView()
{
  // allocate vector of ptrs
  MRTR::Point** points;
  if (Clip_SizePolygon())
    points = new MRTR::Point*[Clip_SizePolygon()];
  else 
    return NULL;
  
  // get the point views
  int count=0;
  map<int,MRTR::Point*>::iterator pcurr;
  for (pcurr=p_.begin(); pcurr != p_.end(); ++pcurr)
  {
    points[count] = pcurr->second;
    ++count;
  }
  if (count != Clip_SizePolygon())
  {
    cout << "***ERR*** MRTR::Overlap::Clip_PointView:\n"
         << "***ERR*** number of point wrong\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  return points;
}














































#endif // TRILINOS_PACKAGE
