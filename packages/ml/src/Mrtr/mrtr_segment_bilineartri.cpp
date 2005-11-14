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

#include "mrtr_segment_bilineartri.H"
#include "mrtr_interface.H"
#include "mrtr_utils.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 10/05|
 |  id               (in)  a unique segment id                          |
 |  nnode            (in)  number of nodes on this segment              |
 |  nodeId           (in)  unique node ids of nodes on this segment     |
 |                         nodeIds have to be sorted on input such that |
 |                         1D case: end nodes of segments first going   |
 |                                  through segment in mathematical     |
 |                                  positive sense. This is             |
 |                                  important to compute the direction  |
 |                                  of the outward normal of the segment|
 |                         2D case: corner nodes of segment first in    |
 |                                  counterclockwise order              |
 |                                  This is                             |
 |                                  important to compute the direction  |
 |                                  of the outward normal of the segment|
 *----------------------------------------------------------------------*/
MRTR::Segment_BiLinearTri::Segment_BiLinearTri(int id, int nnode, int* nodeId) :
MRTR::Segment(id,nnode,nodeId)
{
  stype_ = MRTR::Segment::seg_BiLinearTri;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 10/05|
 |  This constructor should not be used by the user, it is used         |
 |  internally together with Pack/Unpack for communication              |
 *----------------------------------------------------------------------*/
MRTR::Segment_BiLinearTri::Segment_BiLinearTri() :
MRTR::Segment()
{
  stype_ = MRTR::Segment::seg_none;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 10/05|
 *----------------------------------------------------------------------*/
MRTR::Segment_BiLinearTri::Segment_BiLinearTri(MRTR::Segment_BiLinearTri& old) :
MRTR::Segment(old)
{
  // all date lives in the base class and is copied in MRTR::Segment(old)
}

/*----------------------------------------------------------------------*
 | pack all data in this segment into a vector               mwgee 10/05|
 *----------------------------------------------------------------------*/
int* MRTR::Segment_BiLinearTri::Pack(int* size)
{ 
  // note: first there has to be the size and second there has to be the type
  // *size = *size  + stype_ + Id_ + nodeId_.size() + nnode*sizeof(int) + Nfunctions() + 2*Nfunctions()*sizeof(int)
     *size =  1     +	1    +  1  +   1            +	nodeId_.size()  +      1       +      2*Nfunctions();
  int* pack = new int[*size];
  
  int count=0;
  
  pack[count++] = *size;
  pack[count++] = (int)stype_;
  pack[count++] = Id_;
  pack[count++] = nodeId_.size();
  for (int i=0; i<nodeId_.size(); ++i) 
    pack[count++] = nodeId_[i];
  pack[count++] = Nfunctions();
  
  map<int,RefCountPtr<MRTR::Function> >::iterator curr;
  for (curr = functions_.begin(); curr != functions_.end(); ++curr)
  {
    pack[count++] = curr->first;
    pack[count++] = curr->second->Type();
  }
  
  if (count != *size)
  {
    cout << "***ERR*** MRTR::Segment_BiLinearTri::Pack:\n"
         << "***ERR*** mismatch in packing size\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);     
  }
  
  return pack;
}

/*----------------------------------------------------------------------*
 | unpack all data from a vector into this class             mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MRTR::Segment_BiLinearTri::UnPack(int* pack)
{ 
  // note: first there has to be the size and second there has to be the type
  int count = 0;
  int size  = pack[count++];
  stype_    = (MRTR::Segment::SegmentType)pack[count++];
  Id_       = pack[count++];
  nodeId_.resize(pack[count++]);
  for (int i=0; i<nodeId_.size(); ++i)
    nodeId_[i] = pack[count++];
  
  int nfunc = pack[count++];
  
  for (int i=0; i<nfunc; ++i)
  {
    int id   = pack[count++];
    int type = pack[count++];
    MRTR::Function* func = MRTR::AllocateFunction((MRTR::Function::FunctionType)type);
    RefCountPtr<MRTR::Function> rcptrfunc = rcp(func);
    functions_.insert(pair<int,RefCountPtr<MRTR::Function> >(id,rcptrfunc));
  }
  
  if (count != size)
  {
    cout << "***ERR*** MRTR::Segment_BiLinearTri::UnPack:\n"
         << "***ERR*** mismatch in packing size\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);     
  }
  
  return true;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 10/05|
 *----------------------------------------------------------------------*/
MRTR::Segment_BiLinearTri::~Segment_BiLinearTri()
{ 
  // data held in base class is destroyed by base class destructor
}

/*----------------------------------------------------------------------*
 |  clone this segment (public)                              mwgee 10/05|
 *----------------------------------------------------------------------*/
MRTR::Segment* MRTR::Segment_BiLinearTri::Clone()
{ 
  MRTR::Segment_BiLinearTri* newseg = new MRTR::Segment_BiLinearTri(*this);
  return (newseg);
}

/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 10/05|
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const MRTR::Segment_BiLinearTri& seg)
{
  seg.Print(); 
  return os;
}

/*----------------------------------------------------------------------*
 | build an outward normal at a node adjacent to this        mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MRTR::Segment_BiLinearTri::LocalCoordinatesOfNode(int lid, double* xi)
{ 
  if (lid==0)
  {
    xi[0] = 0.0;
    xi[1] = 0.0;
  }
  else if (lid==1)
  {
    xi[0] = 1.0;
    xi[1] = 0.0;
  }
  else if (lid==2)
  {
    xi[0] = 0.0;
    xi[1] = 1.0;
  }
  else
  {
    cout << "***ERR*** MRTR::Segment_BiLinearTri::LocalCoordinatesOfNode:\n"
         << "***ERR*** local node number " << lid << " out of range (0..2)\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);     
  }
  return true;
}

/*----------------------------------------------------------------------*
 | build basis vectors and metric at given point xi          mwgee 10/05|
 *----------------------------------------------------------------------*/
double MRTR::Segment_BiLinearTri::Metric(double* xi, double g[], double G[][3])
{ 
  cout << "***ERR*** MRTR::Segment_BiLinearTri::Metric:\n"
       << "***ERR*** not impl.\n"
       << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
  exit(EXIT_FAILURE);     
  return 0.0;
}

/*----------------------------------------------------------------------*
 | build an outward normal at a node adjacent to this        mwgee 10/05|
 | returns allocated vector of length 3 with outward normal             |
 *----------------------------------------------------------------------*/
double* MRTR::Segment_BiLinearTri::BuildNormal(double* xi)
{ 
  // linear triangles are planar, so we don't care were exactly to build the normal

  // build base vectors
  double g1[3];
  double g2[3];
  for (int i=0; i<3; ++i)
  {
    g1[i] = Nodes()[1]->X()[i] - Nodes()[0]->X()[i];  
    g2[i] = Nodes()[2]->X()[i] - Nodes()[0]->X()[i]; 
  }
  
  // build normal as their cross product
  double* n = new double[3];
  
  MRTR::cross(n,g1,g2);

  return n;
}

/*----------------------------------------------------------------------*
 | compute the length (Area) of this segment                 mwgee 10/05|
 *----------------------------------------------------------------------*/
double MRTR::Segment_BiLinearTri::Area()
{ 
  double xi[2];
  xi[0] = xi[1] = 0.0;
  
  double* n = BuildNormal(xi);

  double a = 0.0;
  for (int i=0; i<3; ++i)
    a+= n[i]*n[i];
  
  delete [] n;

  return (sqrt(a)/2.0);
}



#endif // TRILINOS_PACKAGE
