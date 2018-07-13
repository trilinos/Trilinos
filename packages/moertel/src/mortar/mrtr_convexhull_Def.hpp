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
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Questions? Contact Glen Hansen (gahanse@sandia.gov)
#
# ************************************************************************
#@HEADER
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
#include "mrtr_overlap.hpp"
#include "mrtr_projector.H"
#include "mrtr_node.H"
#include "mrtr_segment.H"
#include "mrtr_interface.H"
#include "mrtr_utils.H"

/*----------------------------------------------------------------------*
 |  create a convexhull of a set of points (private)         mwgee 10/05|
 *----------------------------------------------------------------------*/
template <class IFace>
bool MOERTEL::Overlap<IFace>::ConvexHull(std::map<int,Teuchos::RCP<MOERTEL::Point> >& p)
{
  // # points
  int np = p.size();

  // sort points by xi[0] coordinate
  std::vector<Teuchos::RCP<MOERTEL::Point> > points; PointView(p,points);
  std::vector<double> dlistx(np);
  std::vector<double> dlisty(np);
  std::vector<int> list2(np);
  for (int i=0; i<np; ++i)
  {
    dlistx[i] = points[i]->Xi()[0];
    list2[i] = i;
  }
  MOERTEL::sort(&dlistx[0],np,&list2[0]);

  // get assoc. y-list
  for (int i=0; i<np; ++i)
    dlisty[i] = points[list2[i]]->Xi()[1];

  // sort list again for those points where x is identical
  for (int i=0; i<np-1; ++i)
    if (dlistx[i]==dlistx[i+1])
      if (dlisty[i]>dlisty[i+1])
      {
        MOERTEL::swap(dlistx[i],dlistx[i+1]);
        MOERTEL::swap(dlisty[i],dlisty[i+1]);
        MOERTEL::swap(list2[i],list2[i+1]);
      }  
  
  // create a new polygon and put points in in sorted order
  std::map<int,Teuchos::RCP<MOERTEL::Point> > newp;
  for (int i=0; i<np; ++i)
  {
	Teuchos::RCP<MOERTEL::Point> tmp = Teuchos::rcp(new MOERTEL::Point(i,points[list2[i]]->Xi(),OutLevel()));
    newp.insert(std::pair<int,Teuchos::RCP<MOERTEL::Point> >(i,tmp));
  }
  
  // delete the old one
  p.clear();
  // copy new one over
  CopyPointPolygon(newp,p);
  // destroy the new one
  newp.clear();
  
  points.clear();

  std::map<int,Teuchos::RCP<MOERTEL::Point> >::iterator pcurr;

#if 0
  // printout the polygon
  std::cout << "Input polygon:\n";
  for (pcurr=p.begin(); pcurr != p.end(); ++pcurr)
    if (pcurr->second != Teuchos::null)
      std::cout << *(pcurr->second);
#endif  

  //===========================================================================
  // build the upper hull
  std::map<int,Teuchos::RCP<MOERTEL::Point> > upper;
  PointView(p,points);
  // put in the first 2 points
  AddPointtoPolygon(upper,0,points[0]->Xi()); //std::cout << *points[0];
  AddPointtoPolygon(upper,1,points[1]->Xi()); //std::cout << *points[1];
  //---------------------------------------------------------------------------
  for (int i=2; i<np; ++i)
  {
    // add point[i] to upper
    AddPointtoPolygon(upper,i,points[i]->Xi()); //std::cout << *points[i];

    // find whether we still have a convex hull
    while (upper.size()>2 && !MakeRightTurnUpper(i,upper))
      RemovePointBefore(i,upper);
#if 0
  // printout the current upper hull
  std::map<int,MOERTEL::Point*>::iterator pcurr;
  for (pcurr=upper.begin(); pcurr != upper.end(); ++pcurr)
    if (pcurr->second)
      std::cout << *(pcurr->second);
#endif  
  } // for (int i=2; i<np; ++i)

  //===========================================================================
  // build the lower hull
  std::map<int,Teuchos::RCP<MOERTEL::Point> > lower;
  // put in the first 2 points
  AddPointtoPolygon(lower,np-1,points[np-1]->Xi()); //std::cout << *points[np-1];
  AddPointtoPolygon(lower,np-2,points[np-2]->Xi()); //std::cout << *points[np-2];
  //---------------------------------------------------------------------------
  for (int i=np-3; i>=0; --i)
  {
    // add point[i] to lower
    AddPointtoPolygon(lower,i,points[i]->Xi()); //std::cout << *points[i];

    // find whether we still have a convex hull
    while (lower.size()>2 && !MakeRightTurnLower(i,lower))
      RemovePointAfter(i,lower);
#if 0
  // printout the current lower hull
  std::map<int,Teuchos::RCP<MOERTEL::Point> >::iterator pcurr;
  for (pcurr=lower.begin(); pcurr != lower.end(); ++pcurr)
    if (pcurr->second)
      std::cout << *(pcurr->second);
#endif  
  } // for (int i=np-3; i>=0; --i)  

  //===========================================================================
  points.clear();


  //===========================================================================
  // join upper and lower hull
  // note not to put in duplicate start and end point
  std::map<int,Teuchos::RCP<MOERTEL::Point> > finalp;
  
  // put upper hull in
  int i=0;
  for (pcurr=upper.begin(); pcurr != upper.end(); ++pcurr)
  {
    //std::cout << *(pcurr->second);
    AddPointtoPolygon(finalp,i,pcurr->second->Xi());
    ++i;
  }
  
  // put lower hull in, skip first and last point
  pcurr = lower.end();
  --pcurr; --pcurr;
  for (; pcurr != lower.begin(); --pcurr)
  {
    //std::cout << *(pcurr->second);
    AddPointtoPolygon(finalp,i,pcurr->second->Xi());
    ++i;
  }

#if 0
  // printout the polygon
  std::cout << "--------------------------------------------\n";
  std::cout << "Final polygon:\n";
  for (pcurr=finalp.begin(); pcurr != finalp.end(); ++pcurr)
    if (pcurr->second != Teuchos::null)
      std::cout << *(pcurr->second);
#endif  

  // compare size of convex hull polygon newp to size of p, all
  // nodes must be part of the convex hull
  if (finalp.size() != p.size())
  {
    if (OutLevel()>8)
    std::cout << "MOERTEL: ***WRN*** MOERTEL::Overlap::ConvexHull:\n"
         << "MOERTEL: ***WRN*** size of convex hull " << finalp.size() << " not # nodes " << p.size() << std::endl
         << "MOERTEL: ***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
  }
  
  // copy the polygon over to the input map p
  p.clear();
  upper.clear();
  lower.clear();
  CopyPointPolygon(finalp,p);
  finalp.clear();

  return true;
}


/*----------------------------------------------------------------------*
 |  test whether three points make a right turn (private)    mwgee 10/05|
 *----------------------------------------------------------------------*/
template <class IFace>
bool MOERTEL::Overlap<IFace>::MakeRightTurnUpper(int i,std::map<int,Teuchos::RCP<MOERTEL::Point> >& hull)
{
  // note:
  // point i for sure exists as it was added as last point
  // the points i-1 and i-2 do not necessary have ids i-1 and i-2, they 
  // are just the 2 point BEFORE i (could have any id < i)
  std::map<int,Teuchos::RCP<MOERTEL::Point> >::iterator curr = hull.find(i);

  if (curr==hull.end()) {

	  std::stringstream oss;
			oss << "MOERTEL: ***ERR*** MOERTEL::Overlap::MakeRightTurn:\n"
         << "MOERTEL: ***ERR*** cannot find point " << i << " in convex hull\n"
         << "MOERTEL: ***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	  throw ReportError(oss);
  }

  Teuchos::RCP<MOERTEL::Point> point = curr->second; //std::cout << *point;
  curr--;
  Teuchos::RCP<MOERTEL::Point> pointm1 = curr->second; //std::cout << *pointm1;
  curr--;
  Teuchos::RCP<MOERTEL::Point> pointm2 = curr->second; //std::cout << *pointm2;
  double N[2];
  N[0] =  pointm1->Xi()[1] - pointm2->Xi()[1];
  N[1] = -(pointm1->Xi()[0] - pointm2->Xi()[0]);
  double P[2];
  P[0] = point->Xi()[0] - pointm1->Xi()[0];
  P[1] = point->Xi()[1] - pointm1->Xi()[1];
  
  double dotp = MOERTEL::dot(N,P,2);
  if (dotp>=0.0000)
  {
    //std::cout << "Makes a right\n";
    return true;
  }
  else 
  {
    //std::cout << "Makes NO right\n";
    return false;
  }
}

/*----------------------------------------------------------------------*
 |  test whether three points make a right turn (private)    mwgee 10/05|
 *----------------------------------------------------------------------*/
template <class IFace>
bool MOERTEL::Overlap<IFace>::MakeRightTurnLower(int i,std::map<int,Teuchos::RCP<MOERTEL::Point> >& hull)
{
  // note:
  // point i for sure exists as it was added as last point
  // the points i-1 and i-2 do not necessary have ids i-1 and i-2, they 
  // are just the 2 point BEFORE i (could have any id < i)
  std::map<int,Teuchos::RCP<MOERTEL::Point> >::iterator curr = hull.find(i);

  if (curr==hull.end()) {

	  std::stringstream oss;
			oss << "MOERTEL: ***ERR*** MOERTEL::Overlap::MakeRightTurn:\n"
				<< "MOERTEL: ***ERR*** cannot find point " << i << " in convex hull\n"
				<< "MOERTEL: ***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	  throw ReportError(oss);

  }

  Teuchos::RCP<MOERTEL::Point> point = curr->second; //std::cout << *point;
  curr++;
  Teuchos::RCP<MOERTEL::Point> pointm1 = curr->second; //std::cout << *pointm1;
  curr++;
  Teuchos::RCP<MOERTEL::Point> pointm2 = curr->second; //std::cout << *pointm2;
  double N[2];
  N[0] =  pointm1->Xi()[1] - pointm2->Xi()[1];
  N[1] = -(pointm1->Xi()[0] - pointm2->Xi()[0]);
  double P[2];
  P[0] = point->Xi()[0] - pointm1->Xi()[0];
  P[1] = point->Xi()[1] - pointm1->Xi()[1];
  
  double dotp = MOERTEL::dot(N,P,2);
  if (dotp>=0.0000)
  {
    //std::cout << "Makes a right\n";
    return true;
  }
  else 
  {
    //std::cout << "Makes NO right\n";
    return false;
  }
}

/*----------------------------------------------------------------------*
 |  test whether three points make a right turn (private)    mwgee 10/05|
 *----------------------------------------------------------------------*/
template <class IFace>
void MOERTEL::Overlap<IFace>::RemovePointBefore(int i,std::map<int,Teuchos::RCP<MOERTEL::Point> >& hull)
{
  // note:
  // point i for sure exists as it was added as last point
  // the points i-1 does not necessary have id i-1  
  std::map<int,Teuchos::RCP<MOERTEL::Point> >::iterator curr = hull.find(i);

  if (curr==hull.end()) {

	  std::stringstream oss;
			oss << "MOERTEL: ***ERR*** MOERTEL::Overlap::RemovePointBefore:\n"
				<< "MOERTEL: ***ERR*** cannot find point " << i << " in convex hull\n"
				<< "MOERTEL: ***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	  throw ReportError(oss);
  }

  curr--;
  //std::cout << "Erasing point " << curr->first << " from hull\n";
  hull.erase(curr->first);
  return;
}

/*----------------------------------------------------------------------*
 |  test whether three points make a right turn (private)    mwgee 10/05|
 *----------------------------------------------------------------------*/
template <class IFace>
void MOERTEL::Overlap<IFace>::RemovePointAfter(int i,std::map<int,Teuchos::RCP<MOERTEL::Point> >& hull)
{
  // note:
  // point i for sure exists as it was added as last point
  // the points i-1 does not necessary have id i-1  
  std::map<int,Teuchos::RCP<MOERTEL::Point> >::iterator curr = hull.find(i);

  if (curr==hull.end()) {

	  std::stringstream oss;
			oss << "MOERTEL: ***ERR*** MOERTEL::Overlap::RemovePointBefore:\n"
				<< "MOERTEL: ***ERR*** cannot find point " << i << " in convex hull\n"
				<< "MOERTEL: ***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	  throw ReportError(oss);
  }

  curr++;
  //std::cout << "Erasing point " << curr->first << " from hull\n";
  hull.erase(curr->first);
  return;
}


/*----------------------------------------------------------------------*
 |  collapse points that are really close to one point       mwgee 11/05|
 *----------------------------------------------------------------------*/
template <class IFace>
bool MOERTEL::Overlap<IFace>::CollapsePoints(std::map<int,Teuchos::RCP<MOERTEL::Point> >& p,
                                   const double eps)
{
  // we don't want to collapse on a polygon that has just three or less points
  if (p.size() <= 3)
    return true;
  
  // # points
  int np = p.size();

  // get a view of all points
  std::vector<Teuchos::RCP<MOERTEL::Point> > points; 
  PointView(p,points);

  // create a new std::map for collapsed points
  std::map<int,Teuchos::RCP<MOERTEL::Point> > pnew;

  // create  vector holding points to collapse
  std::vector<Teuchos::RCP<MOERTEL::Point> > collapse(points.size());

  // loop points and compare coords
  for (int i=0; i<np; ++i)
  {
    if (points[i] == Teuchos::null)
      continue;
    
    // put point i into vector with collapse points
    // std::cout << "Adding " << i << " to collapse\n";
    collapse[0] = points[i];
    int count = 1;
    
    for (int j=i+1; j<np; ++j)
    {
      if (points[j] == Teuchos::null)
        continue;
      double xi1xi2[2];
      xi1xi2[0] = points[j]->Xi()[0] - points[i]->Xi()[0];
      xi1xi2[1] = points[j]->Xi()[1] - points[i]->Xi()[1];
      double dist = MOERTEL::length(xi1xi2,2);
      // std::cout << "distance between " << i << " and " << j << " : " << dist << std::endl;
      if (dist<eps)
      {
        // std::cout << "Adding " << j << " to collapse\n";
        // add point2 to collapse vector
        collapse[count] = points[j];
        ++count;
        points[j] = Teuchos::null;
      }
    }
    
    // loop all nodes in the collapse vector and put in the one with an
    // id above 100 (if there is any)
    bool foundit = false;
    //if (count>1)
    //  std::cout << "Collapsing " << count << " nodes\n";
    for (int j=0; j<count; ++j)
    {
      if (collapse[j]->Id()>=100)
      {
        AddPointtoPolygon(pnew,collapse[j]->Id(),collapse[j]->Xi());
        foundit = true;
        break;
      }
    }
    if (!foundit) // there is no point with id >= 100
    {
      AddPointtoPolygon(pnew,collapse[0]->Id(),collapse[0]->Xi());
    }
  } // for (int i=0; i<np; ++i)

  points.clear();
  
  // the new polygon is supposed to have at least three points, otherwise its
  // no better than the old uncollapsed one and we'd rather keep the old one
  // ^^^ this comment is inaccurate ^^^
  // soon after this method is called, there is a check for a degenerate overlap
  // by asking how many things are in p.  If p is less than three, the overlap
  // is skipped, otherwise the computation continues.  After that, an area weighted
  // sum is computed on the polygon.  If the area is 0, this generates NaNs.
  // By commenting out this if block, the check for degeneracy will work as
  // claimed in the comments.
  /* */
  if (pnew.size() < 3)
  {
    pnew.clear();
    collapse.clear();
    return true;
  }
  /* */

  p.clear();
  CopyPointPolygon(pnew,p);
  pnew.clear();  

  return true;
}
