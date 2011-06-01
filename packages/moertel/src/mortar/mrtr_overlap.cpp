/*
#@HEADER
# ************************************************************************
#
#                          Moertel FE Package
#                 Copyright (2006) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? Contact Glen Hansen (Glen.Hansen@inl.gov)
#
# ************************************************************************
#@HEADER
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
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
MOERTEL::Overlap::Overlap(MOERTEL::Segment& sseg, MOERTEL::Segment& mseg, 
                          MOERTEL::Interface& inter, bool exactvalues, int outlevel) :
inter_(inter),
sseg_(sseg),
mseg_(mseg),
outputlevel_(outlevel),
overlap_(false),
havemxi_(false),
havesxi_(false),
havelines_(false),
havesxim_(false),
havelinem_(false),
exactvalues_(exactvalues)
{
  if ( (sseg.Type()!=MOERTEL::Segment::seg_BiLinearTri &&
        sseg.Type()!=MOERTEL::Segment::seg_BiLinearQuad    )  || 
       (mseg.Type()!=MOERTEL::Segment::seg_BiLinearTri && 
        mseg.Type()!=MOERTEL::Segment::seg_BiLinearQuad    )  )
  {
			std::stringstream oss;
		oss << "***ERR*** MOERTEL::Overlap::Overlap:\n"
         << "***ERR*** Overlap of other then bilinear triangles/quads not yet implemented\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
		throw ReportError(oss);
  }
  p_.clear();
  s_.clear();

}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 10/05|
 *----------------------------------------------------------------------*/
MOERTEL::Overlap::~Overlap()
{
  // destroy the point map
  p_.clear();
  // destroy the segment map
  s_.clear();

}

/*----------------------------------------------------------------------*
 |  build line info from triangles/quads (private)           mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MOERTEL::Overlap::build_lines_s()
{
  if (!havemxi_)
  {
	std::cout << "***ERR*** MOERTEL::Overlap::build_lines:\n"
         << "***ERR*** projection of master element nodes has to be done before\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }

  if (sseg_.Nnode()==3) {

    // slave segment, line 0
    sseg_.LocalCoordinatesOfNode(0,&sline_[0][0]);
    sseg_.LocalCoordinatesOfNode(1,&sline_[0][2]);

    // slave segment, line 1
    sseg_.LocalCoordinatesOfNode(1,&sline_[1][0]);
    sseg_.LocalCoordinatesOfNode(2,&sline_[1][2]);

    // slave segment, line 2
    sseg_.LocalCoordinatesOfNode(2,&sline_[2][0]);
    sseg_.LocalCoordinatesOfNode(0,&sline_[2][2]);

  }
  else if (sseg_.Nnode()==4) {

    // slave segment, line 0
    sseg_.LocalCoordinatesOfNode(0,&sline_[0][0]);
    sseg_.LocalCoordinatesOfNode(1,&sline_[0][2]);
    // slave segment, line 1
    sseg_.LocalCoordinatesOfNode(1,&sline_[1][0]);
    sseg_.LocalCoordinatesOfNode(2,&sline_[1][2]);
    // slave segment, line 2
    sseg_.LocalCoordinatesOfNode(2,&sline_[2][0]);
    sseg_.LocalCoordinatesOfNode(3,&sline_[2][2]);
    // slave segment, line 3
    sseg_.LocalCoordinatesOfNode(3,&sline_[3][0]);
    sseg_.LocalCoordinatesOfNode(0,&sline_[3][2]);
  }
  else {

			std::stringstream oss;
				oss << "***ERR*** MOERTEL::Overlap::build_lines_s:\n"
					<< "***ERR*** # slave node " << sseg_.Nnode()
					<< " not supported\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
			throw ReportError(oss);
  }

/*
	mxi_[i][j] - the i index is the ith corner node of the master segment. The j index selects either the
	\xi (j = 0) or the \eta (j = 1) coordinate of the projection of the node into sseg.

	mline_[i][j] - the i index selects the outside edge of the master segment. There are
	3 edges for triangular segments, and 4 for quads. j = 0 is the \xi coord in sseg of the beginning of the edge,
	j = 1 is the \eta coord, j = 2 is the \xi coord in sseg of the end point of the edge,
	and j = 3 is the \eta coord of the end point of the edge.
*/

  if (mseg_.Nnode()==3) {

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
  }
  else if (mseg_.Nnode()==4)
  {
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
    mline_[2][2] = mxi_[3][0];
    mline_[2][3] = mxi_[3][1];

    // master segment, line 3
    mline_[3][0] = mxi_[3][0];
    mline_[3][1] = mxi_[3][1];
    mline_[3][2] = mxi_[0][0];
    mline_[3][3] = mxi_[0][1];

  }
  else {

			std::stringstream oss;
				oss << "***ERR*** MOERTEL::Overlap::build_lines_s:\n"
					<< "***ERR*** # master node " << mseg_.Nnode()
					<< " not supported\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
			throw ReportError(oss);

  }
  
  return true;

}

/*----------------------------------------------------------------------*
 |  build line info from triangles (private)                 mwgee 11/05|
 *----------------------------------------------------------------------*/
bool MOERTEL::Overlap::build_lines_m()
{
  if (!havesxim_)
  {
	std::cout << "***ERR*** MOERTEL::Overlap::build_lines:\n"
         << "***ERR*** projection of slave element nodes has to be done before\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }

  if (mseg_.Nnode()==3)
  {
    // master segment, line 0
    mseg_.LocalCoordinatesOfNode(0,&mlinem_[0][0]);
    mseg_.LocalCoordinatesOfNode(1,&mlinem_[0][2]);
    // master segment, line 1
    mseg_.LocalCoordinatesOfNode(1,&mlinem_[1][0]);
    mseg_.LocalCoordinatesOfNode(2,&mlinem_[1][2]);
    // master segment, line 2
    mseg_.LocalCoordinatesOfNode(2,&mlinem_[2][0]);
    mseg_.LocalCoordinatesOfNode(0,&mlinem_[2][2]);
  }
  else if (mseg_.Nnode()==4)
  {
    // master segment, line 0
    mseg_.LocalCoordinatesOfNode(0,&mlinem_[0][0]);
    mseg_.LocalCoordinatesOfNode(1,&mlinem_[0][2]);
    // master segment, line 1
    mseg_.LocalCoordinatesOfNode(1,&mlinem_[1][0]);
    mseg_.LocalCoordinatesOfNode(2,&mlinem_[1][2]);
    // master segment, line 2
    mseg_.LocalCoordinatesOfNode(2,&mlinem_[2][0]);
    mseg_.LocalCoordinatesOfNode(3,&mlinem_[2][2]);
    // master segment, line 3
    mseg_.LocalCoordinatesOfNode(3,&mlinem_[3][0]);
    mseg_.LocalCoordinatesOfNode(0,&mlinem_[3][2]);
  }
  else {

			std::stringstream oss;
				oss << "***ERR*** MOERTEL::Overlap::build_lines_s:\n"
					<< "***ERR*** # master node " << sseg_.Nnode()
					<< " not supported\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
			throw ReportError(oss);
  }

  if (sseg_.Nnode()==3) {

    // slave segment, line 0
    slinem_[0][0] = sxim_[0][0];
    slinem_[0][1] = sxim_[0][1];
    slinem_[0][2] = sxim_[1][0];
    slinem_[0][3] = sxim_[1][1];
    // slave segment, line 1
    slinem_[1][0] = sxim_[1][0];
    slinem_[1][1] = sxim_[1][1];
    slinem_[1][2] = sxim_[2][0];
    slinem_[1][3] = sxim_[2][1];
    // slave segment, line 2
    slinem_[2][0] = sxim_[2][0];
    slinem_[2][1] = sxim_[2][1];
    slinem_[2][2] = sxim_[0][0];
    slinem_[2][3] = sxim_[0][1];
  }
   else if  (sseg_.Nnode()==4)
   {
    // slave segment, line 0
    slinem_[0][0] = sxim_[0][0];
    slinem_[0][1] = sxim_[0][1];
    slinem_[0][2] = sxim_[1][0];
    slinem_[0][3] = sxim_[1][1];
    // slave segment, line 1
    slinem_[1][0] = sxim_[1][0];
    slinem_[1][1] = sxim_[1][1];
    slinem_[1][2] = sxim_[2][0];
    slinem_[1][3] = sxim_[2][1];
    // slave segment, line 2
    slinem_[2][0] = sxim_[2][0];
    slinem_[2][1] = sxim_[2][1];
    slinem_[2][2] = sxim_[3][0];
    slinem_[2][3] = sxim_[3][1];
    // slave segment, line 3
    slinem_[3][0] = sxim_[3][0];
    slinem_[3][1] = sxim_[3][1];
    slinem_[3][2] = sxim_[0][0];
    slinem_[3][3] = sxim_[0][1];
  }
  else {

			std::stringstream oss;
				oss << "***ERR*** MOERTEL::Overlap::build_lines_s:\n"
					<< "***ERR*** # slave node " << mseg_.Nnode()
					<< " not supported\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
			throw ReportError(oss);
  }
  
  
  return true;
}

/*----------------------------------------------------------------------*
 |  project master nodes onto slave element (private)        mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MOERTEL::Overlap::build_mxi()
{
  // project the master segment's nodes onto the slave segment
  MOERTEL::Node** mnode  = mseg_.Nodes();
  MOERTEL::Projector projector(inter_.IsOneDimensional(),OutLevel());
  double gap;
  for (int i=0; i<mseg_.Nnode(); ++i)
  {
    // project node i onto sseg

    projector.ProjectNodetoSegment_SegmentNormal(*mnode[i],sseg_,mxi_[i],gap);

/*
	mxi_[i] output. mxi_[i][0] is the \xi coordinate of the projection in sseg_, and
	mxi_[i][1] is the \eta coordinate of the projection in sseg_.
*/

#if 0
    // check whether i is inside sseg
    if (mxi_[i][0]<=1. && mxi_[i][1]<=abs(1.-mxi_[i][0]) && mxi_[i][0]>=0. && mxi_[i][1]>=0.)
    {
	  std::cout << "OVERLAP: master node in: " << mnode[i]->Id() << endl;
    }
#endif    

  }
  havemxi_ = true;
  return true;
}

/*----------------------------------------------------------------------*
 |  project slave nodes onto master element (private)        mwgee 11/05|
 *----------------------------------------------------------------------*/
bool MOERTEL::Overlap::build_sxim()
{
  // project the slave segment's nodes onto the master segment
  int nsnode = sseg_.Nnode();
  MOERTEL::Node** snode  = sseg_.Nodes();
  MOERTEL::Projector projector(inter_.IsOneDimensional(),OutLevel());
  double gap;
  for (int i=0; i<nsnode; ++i)
  {
    // project node i onto sseg
    projector.ProjectNodetoSegment_NodalNormal(*snode[i],mseg_,sxim_[i],gap);
#if 0
    // check whether i is inside sseg
    if (sxim_[i][0]<=1. && sxim_[i][1]<=abs(1.-sxim_[i][0]) && sxim_[i][0]>=0. && sxim_[i][1]>=0.)
    {
	  std::cout << "OVERLAP: slave node in: " << snode[i]->Id() << endl;
    }
#endif    
  }
  havesxim_ = true;
  return true;
}

/*----------------------------------------------------------------------*
 |  get coords of slave nodes (private)                      mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MOERTEL::Overlap::build_sxi()
{
  // this makes coords of snodes on the sseg local coord system
  if (sseg_.Nnode()==3)
  {
    sxi_[0][0] = 0.;
    sxi_[0][1] = 0.;
    sxi_[1][0] = 1.;
    sxi_[1][1] = 0.;
    sxi_[2][0] = 0.;
    sxi_[2][1] = 1.;
  }
  else if (sseg_.Nnode()==4)
  {
    sxi_[0][0] = -1.;
    sxi_[0][1] = -1.;
    sxi_[1][0] =  1.;
    sxi_[1][1] = -1.;
    sxi_[2][0] =  1.;
    sxi_[2][1] =  1.;
    sxi_[3][0] = -1.;
    sxi_[3][1] =  1.;
  }
  else
  {
			std::stringstream oss;
		oss << "***ERR*** MOERTEL::Overlap::build_sxi:\n"
         << "***ERR*** # slave node " << sseg_.Nnode() << " not supported\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
		throw ReportError(oss);
  }

  havesxi_ = true;
  return true;
}

/*----------------------------------------------------------------------*
 |  outward normal to sseg's edges in local coords (private) mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MOERTEL::Overlap::build_normal()
{
  if (sseg_.Nnode()==3)
  {
    sn_[0][0] = 0.;
    sn_[0][1] = -1.;
    sn_[1][0] = 1.;
    sn_[1][1] = 1.;
    sn_[2][0] = -1.;
    sn_[2][1] = 0.;
  }
  else if (sseg_.Nnode()==4)
  {
    sn_[0][0] = 0.;
    sn_[0][1] = -1.;
    sn_[1][0] = 1.;
    sn_[1][1] = 0.;
    sn_[2][0] = 0.;
    sn_[2][1] = 1.;
    sn_[3][0] = -1.;
    sn_[3][1] = 0.;
  }
  else
  {
			std::stringstream oss;
		oss << "***ERR*** MOERTEL::Overlap::build_normal:\n"
         << "***ERR*** # slave node " << sseg_.Nnode() << " not supported\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
		throw ReportError(oss);
  }
  if (mseg_.Nnode()==3)
  {
    mn_[0][0] = 0.;
    mn_[0][1] = -1.;
    mn_[1][0] = 1.;
    mn_[1][1] = 1.;
    mn_[2][0] = -1.;
    mn_[2][1] = 0.;
  }
  else if (mseg_.Nnode()==4)
  {
    mn_[0][0] = 0.;
    mn_[0][1] = -1.;
    mn_[1][0] = 1.;
    mn_[1][1] = 0.;
    mn_[2][0] = 0.;
    mn_[2][1] = 1.;
    mn_[3][0] = -1.;
    mn_[3][1] = 0.;
  }
  else
  {
			std::stringstream oss;
		oss << "***ERR*** MOERTEL::Overlap::build_normal:\n"
         << "***ERR*** # master node " << sseg_.Nnode() << " not supported\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
		throw ReportError(oss);
  }
  
  return true;
}

/*----------------------------------------------------------------------*
 |  compute the overlap (public)                             mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MOERTEL::Overlap::ComputeOverlap() {

  // Returning false means there is no overlap between mseg and sseg

  bool ok = false;

  // If we pass the quick overlap test, that means that the master and slave segments are close
  // enough to each other that they could interact. 

  ok = QuickOverlapTest();

  if (!ok)
  {
    return false; // cannot interact
  }

  // project master nodes onto slave segment if not already done
  // The master nodes that define the master segment are projected into the slave
  // segment using:
  //  ProjectNodetoSegment_SegmentNormal(*mnode[i],sseg_,mxi_[i],gap);

  if (!havemxi_)
    havemxi_ = build_mxi();

  // build array of local sxi coords of nodes
  // These do not need to be projected as they are the coords of the slave nodes in the
  // slave segment

  havesxi_ = build_sxi();

  // project slave nodes onto master segment
  // These nodal coordinates are projected from the slave segment into the master, using
  //  ProjectNodetoSegment_NodalNormal(*snode[i],mseg_,sxim_[i],gap);

  havesxim_ = build_sxim();
  
  // build outward normal of edges of sseg (in local coords)
  build_normal();
  
  // compute information about the edges of the slave/master polygons
  // in slave polygon coord system

  if (!havelines_)
    havelines_ = build_lines_s();

  // compute information about the edges of the slave/master polygons
  // in master polygon coord system
  if (!havelinem_)
    havelinem_ = build_lines_m();

  // perform clipping algorithm
  // If this returns false, it means that there is no effective overlap between the master and slave
  // segments
  
//  ok = Clipelements();
  ok = ClipelementsSH(); // Use Sutherland-Hodgman algorithm
  
  if (!ok)
    return false;
  
  // make a triangulation of the overlap polygon
  ok = Triangulation();
  if (!ok)
    return false;

  return true;
}

/*----------------------------------------------------------------------*
 |  perform clipping algorithm (private)                     mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MOERTEL::Overlap::Clipelements() {

  if (!havemxi_ || !havesxi_ || !havelines_ || !havesxim_ || !havelinem_){

			std::stringstream oss;
				oss << "***ERR*** MOERTEL::Overlap::Clipelements:\n"
         << "***ERR*** initialization of Overlap class missing\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
			throw ReportError(oss);

  }
  
  const int nmnode = mseg_.Nnode();
  const int nsnode = sseg_.Nnode();
  
  // put all mseg nodes in polygon
  // These coordinates are the corners of the master segment projected into the slave segment's
  // coordinate system



  for (int i=0; i<nmnode; ++i)
    AddPointtoPolygon(100+i,&mline_[i][0]);
  
  //===========================================================================

  // loop edges of slave segment and clip the master edges
  // add all intersection points that can be found

  for (int clipedge=0; clipedge<nsnode; ++clipedge) {

    // point on that clip edge (dim 2). These are the edges of sseg in its \xi \eta space. Note
    // that if sseg is a quad, we are clipping against its parametric space which is a square
    // -1 <= \xi <= 1 and -1 <= \eta <= 1

    double* PE = &sline_[clipedge][0];

    // the outward normal to the clip edge (dim 2)

    double* N = &sn_[clipedge][0];

    // loop edges of master segment and clip them against the clip edge
    // add all intersections

    for (int medge=0; medge<nmnode; ++medge) {

      bool ok;

      // start point of medge with id 100+medge
      double* P0 = &mline_[medge][0];
      //int id0 = 100+medge;
      
      // end point of medge has id 100+(medge+1 < 3)
      double* P1 = &mline_[medge][2];
      
      // find intersection between medge and clipedge
      // if there's an intersection, put it in the polygon
      // the id is 10*clipedge+medge

      double xi[2];
      ok = Clip_Intersect(N,PE,P0,P1,xi);

      if (ok) {
	  
        //std::cout << "OVERLAP Clipelements: adding intersection point " << 10*clipedge+medge << endl;
// GAH - ASSUME less than 10 medges and 10*clipedge+medge < 100
        AddPointtoPolygon(10*clipedge+medge,xi);
      }
    } // for (int medge=0; medge<3; ++medge)
  } // for (int clipedge=0; clipedge<3; ++clipedge)
  
  //===========================================================================

  // loop all clipedges and clip all points incl intersections 
  // that are in the polygon
  // throw away all points that are not inside

  // Again, these are the edges of sseg in its \xi \eta space. Note
  // that if sseg is a quad, we are clipping against its parametric space which is a square
  // -1 <= \xi <= 1 and -1 <= \eta <= 1

  {
    int np = SizePointPolygon();
    
	std::vector<Teuchos::RCP<MOERTEL::Point> > point;
    PointView(point);
    int p;

    for (p=0; p<np; ++p) {

      bool ok = true;
      //std::cout << "OVERLAP Clipelements: now clipping point " << point[p]->Id() << endl;
      const double* P = point[p]->Xi();

      for (int clipedge=0; clipedge<nsnode; ++clipedge) {

        // point on that clip edge (dim 2)
        double* PE = &sline_[clipedge][0];
        // the outward normal to the clip edge (dim 2)
        double* N = &sn_[clipedge][0];
        
        // clip point P against this edge
// GAH - EPSILON P must be epsilon (in a cos(theta) sense) in front of the clip edge (on the normal side)
//					to be removed from the polygon.

        ok = Clip_TestPoint(N,PE,P,1.0e-10);

        if (ok) {

			continue; // leave point in

        }
        else {		// remove this point
      
        //std::cout << "OVERLAP Clipelements: removing point " << point[p]->Id() << endl;
        RemovePointfromPolygon(point[p]->Id(),P);
          break;

      }
            
      } // for (int clipedge=0; clipedge<3; ++clipedge)
      
    } // for (int p=0; p<np; ++p)

    point.clear();

  }

  //===========================================================================
  // loop all slave nodes and clip against master element.
  // put those slave nodes that are inside master element into polygon
  // Note that this works in master segment coords and the node that is put in
  // is put in with slave segment coords as the polygon is completely in
  // slave segment coords

  // master point have ids 100,101,102
  // slave points have ids 1000,1001,1002
  // edge intersections have ids sedge*10+medge
  
  // try only to put slave points in if the polygon is not empty

  int np = SizePointPolygon();
  if (np)
  {
    for (int i=0; i<nsnode; ++i)
    {
      bool ok = true;
      // get the slave point in master coords
      double* P = sxim_[i];
      // loop master clip edges
      for (int clipedge=0; clipedge<nmnode; ++clipedge)
      {
        // point on master clipedge
        double* PE = &mlinem_[clipedge][0];
        // the outward normal to the clip edge (dim 2)
        double* N = &mn_[clipedge][0];
        
        // clip point P against this edge
// GAH - EPSILON clip test point
// GAH - EPSILON P must be epsilon (in a cos(theta) sense) in front of the clip edge (on the normal side)
        ok = Clip_TestPoint(N,PE,P,1.0e-5);
        // put point in
        if (ok) continue;
        else
        {
          ok = false;
          break;
        }
      } // for (int clipedge=0; clipedge<3; ++clipedge)
      // don't put point in
      if (!ok)
        continue;
      else
      {
        //std::cout << "OVERLAP Clipelements: inserting slave point " << 1000+i << " xi="
        //     << sxi_[i][0] << "/" << sxi_[i][1] << endl;
        AddPointtoPolygon(1000+i,sxi_[i]);
      }
    } // for (int i=0; i<3; ++i)
  }

  //===========================================================================
  //===========================================================================
  //===========================================================================
  
#if 0
  // make printout of the polygon so far
  {
    int np    = SizePointPolygon();
	std::vector<Teuchos::RCP<MOERTEL::Point> > point; PointView(point);
    for (int p=0; p<np; ++p)
    {
      std::cout << "OVERLAP Clipelements: point " << std::setw(3) << point[p]->Id()
           << " xi " << point[p]->Xi()[0] 
           << "/" << point[p]->Xi()[1] << std::endl;
    }
    point.clear();
  }
#endif


  //===========================================================================
  // count how many corner nodes of mseg are in and how many
  // intersections there are
  np = SizePointPolygon();

  //===========================================================================
  // if there are no points, there is still the chance that all slave points
  // are inside the master element
  if (!np)
  {
    for (int i=0; i<3; ++i)
    {
      bool ok = true;
      // get the slave point in master coords
      double* P = sxim_[i];
      
      for (int clipedge=0; clipedge<3; ++clipedge)
      {
        // point on master clipedge
        double* PE = &mlinem_[clipedge][0];
        // the outward normal to the clip edge (dim 2)
        double* N = &mn_[clipedge][0];
        
        // clip point P against this edge
// GAH - EPSILON P must be epsilon (in a cos(theta) sense) in front of the clip edge (on the normal side)
        ok = Clip_TestPoint(N,PE,P,1.0e-5);
        // put point in
        if (ok) continue;
        else
        {
          ok = false;
          break;
        }
      } // for (int clipedge=0; clipedge<3; ++clipedge)
      // don't put point in
      if (!ok)
        continue;
      else
      {
        //std::cout << "OVERLAP Clipelements: inserting slave point " << 1000+i << " xi="
        //     << sxi_[i][0] << "/" << sxi_[i][1] << endl;
        AddPointtoPolygon(1000+i,sxi_[i]);
      }
      
    } // for (int i=0; i<3; ++i)
    
    //=========================================================================
    // check again how many points there are inside
    np = SizePointPolygon();
    if (np>2); // all slave points are in master segment, just continue
    else if (np && np <= 2)
    {
      if (inter_.OutLevel()>8)
        std::cout << "MOERTEL: ***WRN*** MOERTEL::Overlap::Clipelements:\n"
             << "MOERTEL: ***WRN*** " << np << " slave nodes seem to be in master segment but no overlap detected\n"
             << "MOERTEL: ***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
    }
    else // with none or less then 3 points in we assume no overlap
      return false;
  }

  //===========================================================================
  // maybe its a good idea to collapse intersections with nodes (that are in)
  // that are really close to each other
  if (p_.size()>3)
    CollapsePoints(p_,1.0e-5);

#if 0
  // make printout of the polygon so far
  {
	std::cout << "--------------------------------------------\n";
    int np    = SizePointPolygon();
	std::vector<Teuchos::RCP<MOERTEL::Point> > point; PointView(point);
    for (int p=0; p<np; ++p)
    {
      std::cout << "OVERLAP Clipelements: point " << std::setw(3) << point[p]->Id()
           << " xi " << point[p]->Xi()[0] 
           << "/" << point[p]->Xi()[1] << std::endl;
    }
    point.clear();
  }
#endif

  //===========================================================================
  // check again
  np = SizePointPolygon();
  if (np && np<3)
  {
    if (OutLevel()>8)
      std::cout << "MOERTEL: ***WRN*** MOERTEL::Overlap::Clipelements:\n"
           << "MOERTEL: ***WRN*** " << np << " nodes in polygon but could not detect overlap\n"
           << "MOERTEL: ***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  else if (!np)
    return false;

  //===========================================================================
  // finish the polygon
  ConvexHull(p_);  
  
#if 0
  // make printout of the polygon so far
  {
    int np    = SizePointPolygon();
    std::vector<Teuchos::RCP<MOERTEL::Point> > point; PointView(point);
    for (int p=0; p<np; ++p)
    {
      std::cout << "OVERLAP Clipelements: point " << std::setw(3) << point[p]->Id()
           << " xi " << point[p]->Xi()[0]
           << "/" << point[p]->Xi()[1] << std::endl;
    }
    point.clear();
  }
#endif
  
  return true;
}

/*---------------------------------------------------------------------------------*
 |  Clipping elements using Sutherland-Hodgman algorithm (1974) (private)  gah 10/10|
 *---------------------------------------------------------------------------------*/

bool MOERTEL::Overlap::buildPoly(std::vector<double>& source_xi, std::vector<double>& source_eta,
	std::vector<double>& target_xi, std::vector<double>& target_eta, double *PE, double *N) {

	bool s_in, p_in;
	double point_s[2], point_p[2];
	double eps = 1.0e-10; // GAH EPSILON
	int index = -1;
	bool ok;

	// Make sure they are empty
	
	target_xi.clear();
	target_eta.clear();

	// Find the first s inside the polygon

	for(unsigned int i = 0; i < source_xi.size(); i++){
        
		// clip point s against this edge (PE = point on edge, N is outward normal of edge)

		point_s[0] = source_xi[i];
		point_s[1] = source_eta[i];

		s_in = Clip_TestPoint(N, PE, point_s, eps); // true if point is inside edge

		if(s_in){ // if point s is inside the clipedge

			index = i; // placeholder
			break;

		}
	}

	if(index < 0){ // There were no points found inside the clipedge

		return false;

	}

	// Note that we are at index i and s_in is true

	// Go to the end of the input poly data

	for(unsigned int i = static_cast<unsigned int>(index); i < source_xi.size(); i++){

		if(i + 1 >= source_xi.size() - 1){ // p is at the beginning of the list

			point_p[0] = source_xi[0];
			point_p[1] = source_eta[0];

		}
		else {

			point_p[0] = source_xi[i + 1];
			point_p[1] = source_eta[i + 1];

		}

		p_in = Clip_TestPoint(N, PE, point_p, eps); // true if point is inside edge

		if(s_in && p_in){  // both s and p are in polygon

			// Store p in new polygon

			target_xi.push_back(point_p[0]);
			target_eta.push_back(point_p[1]);

		}
		else if(s_in && !p_in){  // s is in, and p is out

			double xi[2];
			ok = Guarded_Clip_Intersect(N, PE, point_s, point_p, xi);

			if(!ok){ // Cannot find the intersection between edge and clip edge

				std::stringstream oss;
					oss << "***ERR*** MOERTEL::Overlap::buildPoly:\n"
						<< "***ERR*** Cannot find the intersection between edge and clip edge\n"
						<< "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
				throw ReportError(oss);
			}
			
			// Store i in new polygon

			target_xi.push_back(xi[0]);
			target_eta.push_back(xi[1]);

		}
		else if(!s_in && p_in){ // s is out, p is in

			double xi[2];
			ok = Guarded_Clip_Intersect(N, PE, point_s, point_p, xi);

			if(!ok){ // Cannot find the intersection between edge and clip edge
				std::stringstream oss;
					oss << "***ERR*** MOERTEL::Overlap::buildPoly:\n"
						<< "***ERR*** Cannot find the intersection between edge and clip edge\n"
						<< "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
				throw ReportError(oss);
			}
			
			// Store i in new polygon

			target_xi.push_back(xi[0]);
			target_eta.push_back(xi[1]);
			
			// Store p in new polygon

			target_xi.push_back(point_p[0]);
			target_eta.push_back(point_p[1]);

		}

		// Do nothing if neither s nor p are inside clip edge
		
		// Advance s

		s_in = p_in;
		point_s[0] = point_p[0];
		point_s[1] = point_p[1];

	}

	// Note that here s_in, and point_s are set to the 0th point in the source data
	// We now need to march back up to index
	
	for(int i = 0; i < index; i++){

		point_p[0] = source_xi[i + 1];
		point_p[1] = source_eta[i + 1];

		p_in = Clip_TestPoint(N, PE, point_p, eps); // true if point is inside edge

		if(s_in && p_in){  // both s and p are in polygon

			// Store p in new polygon

			target_xi.push_back(point_p[0]);
			target_eta.push_back(point_p[1]);

		}
		else if(s_in && !p_in){  // s is in, and p is out

			double xi[2];
			ok = Guarded_Clip_Intersect(N, PE, point_s, point_p, xi);

			if(!ok){ // Cannot find the intersection between edge and clip edge
				std::stringstream oss;
					oss << "***ERR*** MOERTEL::Overlap::buildPoly:\n"
						<< "***ERR*** Cannot find the intersection between edge and clip edge\n"
						<< "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
				throw ReportError(oss);
			}
			
			// Store i in new polygon

			target_xi.push_back(xi[0]);
			target_eta.push_back(xi[1]);

		}
		else if(!s_in && p_in){ // s is out, p is in

			double xi[2];
			ok = Guarded_Clip_Intersect(N, PE, point_s, point_p, xi);

			if(!ok){ // Cannot find the intersection between edge and clip edge
				std::stringstream oss;
					oss << "***ERR*** MOERTEL::Overlap::buildPoly:\n"
						<< "***ERR*** Cannot find the intersection between edge and clip edge\n"
						<< "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
				throw ReportError(oss);
			}
			
			// Store i in new polygon

			target_xi.push_back(xi[0]);
			target_eta.push_back(xi[1]);
			
			// Store p in new polygon

			target_xi.push_back(point_p[0]);
			target_eta.push_back(point_p[1]);

		}

		// Do nothing if neither s nor p are inside clip edge
		
		// Advance s

		s_in = p_in;
		point_s[0] = point_p[0];
		point_s[1] = point_p[1];

	}

	return true;  // If we have not thrown an error or returned false, return the polygon

}

bool MOERTEL::Overlap::ClipelementsSH() {

  if (!havemxi_ || !havesxi_ || !havelines_ || !havesxim_ || !havelinem_){

			std::stringstream oss;
				oss << "***ERR*** MOERTEL::Overlap::Clipelements:\n"
					<< "***ERR*** initialization of Overlap class missing\n"
					<< "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
			throw ReportError(oss);

  }
  
  const int nmnode = mseg_.Nnode();
  const int nsnode = sseg_.Nnode();
  bool ok = true;
  double eps = 1.0e-10; // GAH EPSILON


  // I am reading http://cs.fit.edu/~wds/classes/graphics/Clip/clip/clip.html and 
  // http://en.wikipedia.org/wiki/Sutherland–Hodgman_algorithm as I write this code.

  // I assume that the master segment is convex. It should be if this is a mesh. Secondly, the
  // slave segment is a square in its parametric space -1 <= \xi <= 1 and -1 <= \eta <= 1. It need
  // only be convex.

  // For the Sutherland-Hodgman algorithm, the slave segment is used for the clip polygon.

  // Start with an input list of polygon vertices, in the slave coordinate system. Note that these points
  // are ordered around the polygon, as i is the line number of the boundaries of the master seg.

  std::vector<double> s_poly_xi, s_poly_eta, t_poly_xi, t_poly_eta;

  // Put all the corners of the master segment into the polygon

  for (int i=0; i<nmnode; ++i) { // loop over the corners of the master seg

	s_poly_xi.push_back(mline_[i][0]); 
	s_poly_eta.push_back(mline_[i][1]);

  }

  // Clip that poly against the slave poly one edge at a time

  for (int clipedge = 0; clipedge < nsnode; ++clipedge) {

	// point on that clip edge (dim 2)
	double* PE = &sline_[clipedge][0];

	// the outward normal to the clip edge (dim 2)
	double* N = &sn_[clipedge][0];

	ok = buildPoly(s_poly_xi, s_poly_eta, t_poly_xi, t_poly_eta, PE, N);

	if(!ok){

		// there is no intersection between polys. However, it is possible that the
		// slave poly is completely contained within the master
		
		break;

	}

	// The target vectors now become the source vectors
	
	s_poly_xi = t_poly_xi;
	s_poly_eta = t_poly_eta;

  } // for (int clipedge=0; clipedge<3; ++clipedge)

  if(ok){ // We have a target polygon

		double xi[2];

		// Compress the polygon by removing adjacent points less than epsilon apart

		for(unsigned int p = t_poly_xi.size() - 1; p >= 1; p--){

			xi[0] = t_poly_xi[p] - t_poly_xi[p - 1];
			xi[1] = t_poly_eta[p] - t_poly_eta[p - 1];

			double dist = MOERTEL::length(xi, 2);

			if(dist <= eps){ // remove the point p

				t_poly_xi.erase(t_poly_xi.begin() + p);
				t_poly_eta.erase(t_poly_eta.begin() + p);

			}
		}

		// Check the last point too

		if(t_poly_xi.size() >= 2){
		
			xi[0] = t_poly_xi[t_poly_xi.size() - 1] - t_poly_xi[0];
			xi[1] = t_poly_eta[t_poly_xi.size() - 1] - t_poly_eta[0];

			double dist = MOERTEL::length(xi, 2);

			if(dist <= eps){ // remove the point p

				t_poly_xi.erase(t_poly_xi.end() - 1);
				t_poly_eta.erase(t_poly_eta.end() - 1);

			}
		}

		// Do we still have a polygon? If not, just move on

		if(t_poly_xi.size() < 3){

			return false;

		}

		// Store it in the polygon

		for (unsigned int p = 0; p < t_poly_xi.size(); ++p) {

			xi[0] = t_poly_xi[p];
			xi[1] = t_poly_eta[p];

			AddPointtoPolygon(p, xi);

		}

#if 0
		// make printout of the polygon so far
		{
		int np    = SizePointPolygon();
		std::vector<Teuchos::RCP<MOERTEL::Point> > point; PointView(point);

		std::cout << "Master is in slave" << std::endl;

		for (int p=0; p<np; ++p) {

			std::cout << "OVERLAP Clipelements: point " << std::setw(3) << point[p]->Id()
						<< " xi " << point[p]->Xi()[0]
						<< "/" << point[p]->Xi()[1] << std::endl;
		}
		point.clear();
		}
#endif
  
		return true;

  }

  //===========================================================================

  // We come down here if there is no polygon from the above. This could mean that the slave segment is
  // completely inside the master segment.
  
	std::vector<int> s_node_id;

	for (int i=0; i<nsnode; ++i) {

		bool ok = true;

		// get the slave point in master coords
		double* P = sxim_[i];

		// loop master clip edges
		for (int clipedge=0; clipedge<nmnode; ++clipedge) {

			// point on master clipedge
			double* PE = &mlinem_[clipedge][0];

			// the outward normal to the clip edge (dim 2)
			double* N = &mn_[clipedge][0];
        
			// clip point P against this edge
			// GAH - EPSILON clip test point
			ok = Clip_TestPoint(N,PE,P,1.0e-5);

			// put point in
			if (ok) 
				continue;
			else {
				ok = false;
				break;
			}
		} // for (int clipedge=0; clipedge<3; ++clipedge)

		// We will be here, with ok == true only if the point is inside ALL clip edges

		// don't put point in
		if (!ok)
			continue;
		else { // Point is inside ALL clip edges

			s_node_id.push_back(i);
		}
	} // for (int i=0; i<3; ++i)

	if(s_node_id.size() < static_cast<unsigned int>(nsnode)){ // Slave poly does not lie within master either. 
															// There is no overlap. Move on.

/* The reasoning here is that zero of the master polygon was found in the slave (no nodes). If the slave is completely
 * within the master, then all 4 nodes need to cleanly show up inside the master. If they do not, chances are only
 * a corner or edge of the master is touched.
 */

		return false;

	}

	// Slave is completely in master. Put the slave in the polygon
	
	for(unsigned int i = 0; i < s_node_id.size(); i++){

		//std::cout << "OVERLAP Clipelements: inserting slave point " << 1000+i << " xi="
		//     << sxi_[i][0] << "/" << sxi_[i][1] << endl;
		AddPointtoPolygon(i,sxi_[i]);

	}

#if 0
		// make printout of the polygon so far
		{
		int np    = SizePointPolygon();
		std::vector<Teuchos::RCP<MOERTEL::Point> > point; PointView(point);

		std::cout << "Slave is in master" << std::endl;

		for (int p=0; p<np; ++p) {

			std::cout << "OVERLAP Clipelements: point " << std::setw(3) << point[p]->Id()
						<< " xi " << point[p]->Xi()[0]
						<< "/" << point[p]->Xi()[1] << std::endl;
		}
		point.clear();
		}
#endif

	return true;
}



/*----------------------------------------------------------------------*
 |  make a triangulization of a polygon (private)             mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MOERTEL::Overlap::Triangulation()
{
  bool ok = false;
  // we have a polygon that is in clockwise order at this moment

  // if more than 3 points, find a center point
  int np = SizePointPolygon();
  if (np<3)
  {
		std::stringstream oss;
			oss << "***ERR*** MOERTEL::Overlap::Triangulization:\n"
         << "***ERR*** # point in polygon < 3 ... very strange!!!\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
		throw ReportError(oss);

  }
  if (np>3) // we have to add a center point
  {
    double xi[2];
	std::vector<Teuchos::RCP<MOERTEL::Point> > points; PointView(points);

#if 1 // compute center point xi as centroid
    ok = Centroid(xi,points,np);    

	if(!ok) return false; // Area of polygon is zero

    //std::cout << "Centroid xi : " << xi[0] << " / " << xi[1] << endl; fflush(stdout);
#endif
    
#if 0 // compute center point as nodal avergage
    xi[0] = xi[1] = 0.0;
    for (int i=0; i<np; ++i)
    {
      const double* pxi = points[i]->Xi();
      xi[0] += pxi[0];
      xi[1] += pxi[1];
    }
    xi[0] /= np;
    xi[1] /= np;
    //std::cout << "Nodal    xi : " << xi[0] << " / " << xi[1] << endl << endl; fflush(stdout);
#endif
    
    // create a point -1 as center point and add it to the polygon
    AddPointtoPolygon(-1,xi);
    points.clear();
  } // if (np>3)
  
  np = SizePointPolygon();
  std::vector<Teuchos::RCP<MOERTEL::Point> >  points; PointView(points);

  // create a MOERTEL::Node for every point
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
    double val[20];
    sseg_.EvaluateFunction(0,points[i]->Xi(),val,sseg_.Nnode(),NULL);
    MOERTEL::Node** snodes = sseg_.Nodes();
    for (int j=0; j<sseg_.Nnode(); ++j)
      for (int k=0; k<3; ++k)
      {
        x[k] += val[j]*snodes[j]->X()[k];
        n[k] += val[j]*snodes[j]->N()[k];
      }
    const double length = MOERTEL::length(n,3);
    for (int j=0; j<3; ++j) n[j] /= length;
    // create a node with this coords and normal;
    MOERTEL::Node* node = new MOERTEL::Node(points[i]->Id(),x,3,dof,false,OutLevel());
    node->SetN(n);
    // set node in point
    points[i]->SetNode(node);
#if 0
	std::cout << *points[i];
#endif
  }  

  // find projection values for all points in polygon on mseg
  {
    double mxi[2];
	double gap;
    MOERTEL::Projector projector(inter_.IsOneDimensional(),OutLevel());
    for (int i=0; i<np; ++i)
    {
	  Teuchos::RCP<MOERTEL::Node> node = points[i]->Node();
      projector.ProjectNodetoSegment_NodalNormal(*node,mseg_,mxi,gap);
      // create a projected node and set it in node
      MOERTEL::ProjectedNode* pnode = new MOERTEL::ProjectedNode(*node,mxi,&mseg_);
      node->SetProjectedNode(pnode);
	  node->SetGap(gap);
#if 0
	  std::cout << "-------------------------------------------------------\n";
      if (mseg_.Nnode()==3)
      {
        if (mxi[0]<=1. && mxi[1]<=abs(1.-mxi[0]) && mxi[0]>=0. && mxi[1]>=0.)
          std::cout << "OVERLAP: point " << points[i]->Id() << " is in mseg, mxi " << mxi[0] << " / " << mxi[1] << endl;
        else
          std::cout << "OVERLAP: point " << points[i]->Id() << " is NOT in mseg, mxi " << mxi[0] << " / " << mxi[1] << endl;
      }
      else if (mseg_.Nnode()==4)
      {
        if (mxi[0]<=1.001 && mxi[0]>=-1.001 && mxi[1]<=1.001 && mxi[1]>=-1.001)
          std::cout << "OVERLAP: point " << points[i]->Id() << " is in mseg, mxi " << mxi[0] << " / " << mxi[1] << endl;
        else
          std::cout << "OVERLAP: point " << points[i]->Id() << " is NOT in mseg, mxi " << mxi[0] << " / " << mxi[1] << endl;
      }
      else
      {
			std::stringstream oss;
			oss << "***ERR*** MOERTEL::Overlap::Triangulization:\n"
             << "***ERR*** # nodes " << mseg_.Nnode() << " of master segment is unknown\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
			throw ReportError(oss);
      }
	  std::cout << "-------------------------------------------------------\n";
#endif

    }
  }

  // if we plan to interpolate function values at the gaussian points we need the
  // function values at the points
  if (exactvalues_==false)
  {
    // every point has 3 values of the 3 shape functions of the elements:
    // max 4 values from function 0 from sseg
    // max 4 values from function 1 from sseg
    // max 4 values from function 0 from mseg
    for (int i=0; i<np; ++i)
    {
      double val[20];
      int nmval = mseg_.Nnode();
      int nsval = sseg_.Nnode();
      // evaluate function 0 from sseg
      sseg_.EvaluateFunction(0,points[i]->Xi(),val,nsval,NULL);
      points[i]->StoreFunctionValues(0,val,nsval);
      // evaluate function 1 from sseg
      sseg_.EvaluateFunction(1,points[i]->Xi(),val,nsval,NULL);
      points[i]->StoreFunctionValues(1,val,nsval);
      // evaluate function 0 from mseg
      mseg_.EvaluateFunction(0,points[i]->Node()->GetProjectedNode()->Xi(),val,nmval,NULL);
      points[i]->StoreFunctionValues(2,val,nmval);
    }
  }
  
  // create the triangle elements from polygon with centerpoint
  // In case of np==3, there is no centerpoint
  if (np>3)
  {
    // the polygon is in clockwise order, center point is points[0] and has
    // id = -1
    if (points[0]->Id() != -1)
    {
			std::stringstream oss;
		oss << "***ERR*** MOERTEL::Overlap::Triangulization:\n"
           << "***ERR*** points[0]->Id() is not -1\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
			throw ReportError(oss);
    }
    int nodeid[3];
    MOERTEL::Segment_BiLinearTri* tmp;
    MOERTEL::Function_LinearTri* func = new Function_LinearTri(OutLevel());
    for (int i=2; i<np; ++i)
    {
      // there are np-1 triangles
      // triangle ids go from 0 to np-2
      // a triangle is defined by nodes 0, i, i-1
      // *this class takes ownership of triangle
      nodeid[0]   = points[0]->Id();
      nodeid[1]   = points[i]->Id();
      nodeid[2]   = points[i-1]->Id();
      tmp = new MOERTEL::Segment_BiLinearTri(i-2,3,nodeid,OutLevel());
      // set a linear shape function to this triangle
      tmp->SetFunction(0,func);
      // add triangle to the *this class
      AddSegment(tmp->Id(),tmp);
    }
    // add the last triangle defined by nodes 0, 1, np-1 separately
    // *this class takes ownership of triangle
    nodeid[0]   = points[0]->Id();
    nodeid[1]   = points[1]->Id();
    nodeid[2]   = points[np-1]->Id();
    tmp = new MOERTEL::Segment_BiLinearTri(np-2,3,nodeid,OutLevel());
    // set a linear shape function to this triangle
    tmp->SetFunction(0,func);
    // add triangle to the *this class
    AddSegment(tmp->Id(),tmp);
    if (func) delete func; func = NULL;
  }
  else if (np==3) // single triangle without centerpoint
  {
    int nodeid[3];
    MOERTEL::Segment_BiLinearTri* tmp;
    MOERTEL::Function_LinearTri* func = new Function_LinearTri(OutLevel());
    nodeid[0] = points[0]->Id();
    nodeid[1] = points[2]->Id();
    nodeid[2] = points[1]->Id();
    // *this class takes ownership of triangle
    tmp = new MOERTEL::Segment_BiLinearTri(0,3,nodeid,OutLevel());
    // set a linear shape function to this triangle
    tmp->SetFunction(0,func);
    // add triangle the *this class
    AddSegment(tmp->Id(),tmp);
    if (func) delete func; func = NULL;
  }
  else
  {
			std::stringstream oss;
		oss << "***ERR*** MOERTEL::Overlap::Triangulization:\n"
         << "***ERR*** # point in polygon < 3 ... very strange!!!\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
			throw ReportError(oss);
  }
  
  // create ptr topology between triangle and nodes in points
  std::vector<MOERTEL::Node*> nodes(np);
  for (int i=0; i<np; ++i) nodes[i] = points[i]->Node().get();
  
  // loop segments and set ptr to nodes in them
  std::map<int,Teuchos::RCP<MOERTEL::Segment> >::iterator curr;
  for (curr=s_.begin(); curr != s_.end(); ++curr)
    curr->second->GetPtrstoNodes(nodes);
  
  nodes.clear();
  
  points.clear();
  
  return true;
}
