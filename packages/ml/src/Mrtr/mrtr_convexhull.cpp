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
 |  create a convexhull of a set of points (private)         mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MRTR::Overlap::ConvexHull(map<int,RefCountPtr<MRTR::Point> >& p)
{
  // # points
  int np = p.size();

  // sort points by xi[0] coordinate
  vector<RefCountPtr<MRTR::Point> > points; PointView(p,points);
  vector<double> dlistx(np);
  vector<double> dlisty(np);
  vector<int> list2(np);
  for (int i=0; i<np; ++i)
  {
    dlistx[i] = points[i]->Xi()[0];
    list2[i] = i;
  }
  MRTR::sort(&dlistx[0],np,&list2[0]);

  // get assoc. y-list
  for (int i=0; i<np; ++i)
    dlisty[i] = points[list2[i]]->Xi()[1];

  // sort list again for those points where x is identical
  for (int i=0; i<np-1; ++i)
    if (dlistx[i]==dlistx[i+1])
      if (dlisty[i]>dlisty[i+1])
      {
        MRTR::swap(dlistx[i],dlistx[i+1]);
        MRTR::swap(dlisty[i],dlisty[i+1]);
        MRTR::swap(list2[i],list2[i+1]);
      }  
  
  // create a new polygon and put points in in sorted order
  map<int,RefCountPtr<MRTR::Point> > newp;
  for (int i=0; i<np; ++i)
  {
    RefCountPtr<MRTR::Point> tmp = rcp(new MRTR::Point(i,points[list2[i]]->Xi()));
    newp.insert(pair<int,RefCountPtr<MRTR::Point> >(i,tmp));
  }
  
  // delete the old one
  p.clear();
  // copy new one over
  CopyPointPolygon(newp,p);
  // destroy the new one
  newp.clear();
  
  points.clear();

  map<int,RefCountPtr<MRTR::Point> >::iterator pcurr;

#if 0
  // printout the polygon
  for (pcurr=p.begin(); pcurr != p.end(); ++pcurr)
    if (pcurr->second)
      cout << *(pcurr->second);
#endif  

  //===========================================================================
  // build the upper hull
  map<int,RefCountPtr<MRTR::Point> > upper;
  PointView(p,points);
  // put in the first 2 points
  AddPointtoPolygon(upper,0,points[0]->Xi()); //cout << *points[0];
  AddPointtoPolygon(upper,1,points[1]->Xi()); //cout << *points[1];
  //---------------------------------------------------------------------------
  for (int i=2; i<np; ++i)
  {
    // add point[i] to upper
    AddPointtoPolygon(upper,i,points[i]->Xi()); //cout << *points[i];

    // find whether we still have a convex hull
    while (upper.size()>2 && !MakeRightTurnUpper(i,upper))
      RemovePointBefore(i,upper);
#if 0
  // printout the current upper hull
  map<int,MRTR::Point*>::iterator pcurr;
  for (pcurr=upper.begin(); pcurr != upper.end(); ++pcurr)
    if (pcurr->second)
      cout << *(pcurr->second);
#endif  
  } // for (int i=2; i<np; ++i)

  //===========================================================================
  // build the lower hull
  map<int,RefCountPtr<MRTR::Point> > lower;
  // put in the first 2 points
  AddPointtoPolygon(lower,np-1,points[np-1]->Xi()); //cout << *points[np-1];
  AddPointtoPolygon(lower,np-2,points[np-2]->Xi()); //cout << *points[np-2];
  //---------------------------------------------------------------------------
  for (int i=np-3; i>=0; --i)
  {
    // add point[i] to lower
    AddPointtoPolygon(lower,i,points[i]->Xi()); //cout << *points[i];

    // find whether we still have a convex hull
    while (lower.size()>2 && !MakeRightTurnLower(i,lower))
      RemovePointAfter(i,lower);
#if 0
  // printout the current lower hull
  map<int,RefCountPtr<MRTR::Point> >::iterator pcurr;
  for (pcurr=lower.begin(); pcurr != lower.end(); ++pcurr)
    if (pcurr->second)
      cout << *(pcurr->second);
#endif  
  } // for (int i=np-3; i>=0; --i)  

  //===========================================================================
  points.clear();


  //===========================================================================
  // join upper and lower hull
  // note not to put in duplicate start and end point
  map<int,RefCountPtr<MRTR::Point> > finalp;
  
  // put upper hull in
  int i=0;
  for (pcurr=upper.begin(); pcurr != upper.end(); ++pcurr)
  {
    //cout << *(pcurr->second);
    AddPointtoPolygon(finalp,i,pcurr->second->Xi());
    ++i;
  }
  
  // put lower hull in, skip first and last point
  pcurr = lower.end();
  --pcurr; --pcurr;
  for (; pcurr != lower.begin(); --pcurr)
  {
    //cout << *(pcurr->second);
    AddPointtoPolygon(finalp,i,pcurr->second->Xi());
    ++i;
  }

  // compare size of convex hull polygon newp to size of p, all
  // nodes must be part of the convex hull
  if (finalp.size() != p.size())
  {
    cout << "***ERR*** MRTR::Overlap::ConvexHull:\n"
         << "***ERR*** size of convex hull " << finalp.size() << " not # nodes " << p.size() << endl
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
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
bool MRTR::Overlap::MakeRightTurnUpper(int i,map<int,RefCountPtr<MRTR::Point> >& hull)
{
  // note:
  // point i for sure exists as it was added as last point
  // the points i-1 and i-2 do not necessary have ids i-1 and i-2, they 
  // are just the 2 point BEFORE i (could have any id < i)
  map<int,RefCountPtr<MRTR::Point> >::iterator curr = hull.find(i);
  if (curr==hull.end())
  {
    cout << "***ERR*** MRTR::Overlap::MakeRightTurn:\n"
         << "***ERR*** cannot find point " << i << " in convex hull\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  RefCountPtr<MRTR::Point> point = curr->second; //cout << *point;
  curr--;
  RefCountPtr<MRTR::Point> pointm1 = curr->second; //cout << *pointm1;
  curr--;
  RefCountPtr<MRTR::Point> pointm2 = curr->second; //cout << *pointm2;
  double N[2];
  N[0] =  pointm1->Xi()[1] - pointm2->Xi()[1];
  N[1] = -(pointm1->Xi()[0] - pointm2->Xi()[0]);
  double P[2];
  P[0] = point->Xi()[0] - pointm1->Xi()[0];
  P[1] = point->Xi()[1] - pointm1->Xi()[1];
  
  double dotp = MRTR::dot(N,P,2);
  if (dotp>=0.0000)
  {
    //cout << "Makes a right\n";
    return true;
  }
  else 
  {
    //cout << "Makes NO right\n";
    return false;
  }
}

/*----------------------------------------------------------------------*
 |  test whether three points make a right turn (private)    mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MRTR::Overlap::MakeRightTurnLower(int i,map<int,RefCountPtr<MRTR::Point> >& hull)
{
  // note:
  // point i for sure exists as it was added as last point
  // the points i-1 and i-2 do not necessary have ids i-1 and i-2, they 
  // are just the 2 point BEFORE i (could have any id < i)
  map<int,RefCountPtr<MRTR::Point> >::iterator curr = hull.find(i);
  if (curr==hull.end())
  {
    cout << "***ERR*** MRTR::Overlap::MakeRightTurn:\n"
         << "***ERR*** cannot find point " << i << " in convex hull\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  RefCountPtr<MRTR::Point> point = curr->second; //cout << *point;
  curr++;
  RefCountPtr<MRTR::Point> pointm1 = curr->second; //cout << *pointm1;
  curr++;
  RefCountPtr<MRTR::Point> pointm2 = curr->second; //cout << *pointm2;
  double N[2];
  N[0] =  pointm1->Xi()[1] - pointm2->Xi()[1];
  N[1] = -(pointm1->Xi()[0] - pointm2->Xi()[0]);
  double P[2];
  P[0] = point->Xi()[0] - pointm1->Xi()[0];
  P[1] = point->Xi()[1] - pointm1->Xi()[1];
  
  double dotp = MRTR::dot(N,P,2);
  if (dotp>=0.0000)
  {
    //cout << "Makes a right\n";
    return true;
  }
  else 
  {
    //cout << "Makes NO right\n";
    return false;
  }
}

/*----------------------------------------------------------------------*
 |  test whether three points make a right turn (private)    mwgee 10/05|
 *----------------------------------------------------------------------*/
void MRTR::Overlap::RemovePointBefore(int i,map<int,RefCountPtr<MRTR::Point> >& hull)
{
  // note:
  // point i for sure exists as it was added as last point
  // the points i-1 does not necessary have id i-1  
  map<int,RefCountPtr<MRTR::Point> >::iterator curr = hull.find(i);
  if (curr==hull.end())
  {
    cout << "***ERR*** MRTR::Overlap::RemovePointBefore:\n"
         << "***ERR*** cannot find point " << i << " in convex hull\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  curr--;
  //cout << "Erasing point " << curr->first << " from hull\n";
  hull.erase(curr->first);
  return;
}

/*----------------------------------------------------------------------*
 |  test whether three points make a right turn (private)    mwgee 10/05|
 *----------------------------------------------------------------------*/
void MRTR::Overlap::RemovePointAfter(int i,map<int,RefCountPtr<MRTR::Point> >& hull)
{
  // note:
  // point i for sure exists as it was added as last point
  // the points i-1 does not necessary have id i-1  
  map<int,RefCountPtr<MRTR::Point> >::iterator curr = hull.find(i);
  if (curr==hull.end())
  {
    cout << "***ERR*** MRTR::Overlap::RemovePointBefore:\n"
         << "***ERR*** cannot find point " << i << " in convex hull\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  curr++;
  //cout << "Erasing point " << curr->first << " from hull\n";
  hull.erase(curr->first);
  return;
}































#endif // TRILINOS_PACKAGE
