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
#include "mrtr_segment_bilinearquad.H"
#include "mrtr_interface.H"

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
MOERTEL::Segment_BiLinearQuad::Segment_BiLinearQuad(int id, int nnode, int* nodeId, int out) :
MOERTEL::Segment(id,nnode,nodeId,out)
{
  stype_ = MOERTEL::Segment::seg_BiLinearQuad;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 10/05|
 |  This constructor should not be used by the user, it is used         |
 |  internally together with Pack/Unpack for communication              |
 *----------------------------------------------------------------------*/
MOERTEL::Segment_BiLinearQuad::Segment_BiLinearQuad(int out) :
MOERTEL::Segment(out)
{
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 10/05|
 *----------------------------------------------------------------------*/
MOERTEL::Segment_BiLinearQuad::Segment_BiLinearQuad(MOERTEL::Segment_BiLinearQuad& old) :
MOERTEL::Segment(old)
{
  // all date lives in the base class and is copied in MOERTEL::Segment(old)
}

/*----------------------------------------------------------------------*
 | pack all data in this segment into a vector               mwgee 10/05|
 *----------------------------------------------------------------------*/
int* MOERTEL::Segment_BiLinearQuad::Pack(int* size)
{ 
  // note: first there has to be the size and second there has to be the type
  // *size = *size  + stype_ + Id_ + nodeId_.size() + nnode*sizeof(int) + Nfunctions() + 2*Nfunctions()*sizeof(int)
     *size =  1     +	1    +  1  +   1            +	nodeId_.size()  +      1       +      2*Nfunctions();
  
  int* pack = new int[*size];
  int count = 0;
  
  pack[count++] = *size;
  pack[count++] = (int)stype_;
  pack[count++] = Id_;
  pack[count++] = nodeId_.size();
  for (int i=0; i<(int)nodeId_.size(); ++i) 
    pack[count++] = nodeId_[i];
  pack[count++] = Nfunctions();
  std::map<int,Teuchos::RCP<MOERTEL::Function> >::iterator curr;
  for (curr = functions_.begin(); curr != functions_.end(); ++curr)
  {
    pack[count++] = curr->first;
    pack[count++] = curr->second->Type();
  }
  
  if (count != *size)
  {
	  std::stringstream oss;
    oss << "***ERR*** MOERTEL::Segment_BiLinearQuad::Pack:\n"
         << "***ERR*** mismatch in packing size\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	throw ReportError(oss);
  }
  
  return pack;
}

/*----------------------------------------------------------------------*
 | unpack all data from a vector into this class             mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MOERTEL::Segment_BiLinearQuad::UnPack(int* pack)
{ 
  // note: first there has to be the size and second there has to be the type
  int count = 0;
  int size  = pack[count++];
  stype_    = (MOERTEL::Segment::SegmentType)pack[count++];
  Id_       = pack[count++];
  nodeId_.resize(pack[count++]);
  for (int i=0; i<(int)nodeId_.size(); ++i)
    nodeId_[i] = pack[count++];
  int nfunc = pack[count++];
  for (int i=0; i<nfunc; ++i)
  {
    int id   = pack[count++];
    int type = pack[count++];
    MOERTEL::Function* func = MOERTEL::AllocateFunction((MOERTEL::Function::FunctionType)type,OutLevel());
	Teuchos::RCP<MOERTEL::Function> rcptrfunc = Teuchos::rcp(func);
    functions_.insert(std::pair<int,Teuchos::RCP<MOERTEL::Function> >(id,rcptrfunc));
  }
  
  if (count != size)
  {
	  std::stringstream oss;
    oss << "***ERR*** MOERTEL::Segment_BiLinearQuad::UnPack:\n"
         << "***ERR*** mismatch in packing size\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	throw ReportError(oss);
  }
  
  return true;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 10/05|
 *----------------------------------------------------------------------*/
MOERTEL::Segment_BiLinearQuad::~Segment_BiLinearQuad()
{ 
  // data held in base class is destroyed by base class destructor
}

/*----------------------------------------------------------------------*
 |  clone this segment (public)                              mwgee 10/05|
 *----------------------------------------------------------------------*/
MOERTEL::Segment* MOERTEL::Segment_BiLinearQuad::Clone()
{ 
  MOERTEL::Segment_BiLinearQuad* newseg = new MOERTEL::Segment_BiLinearQuad(*this);
  return (newseg);
}

/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 10/05|
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const MOERTEL::Segment_BiLinearQuad& seg)
{
  seg.Print(); 
  return os;
}

/*----------------------------------------------------------------------*
 | build an outward normal at a node adjacent to this        mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MOERTEL::Segment_BiLinearQuad::LocalCoordinatesOfNode(int lid, double* xi)
{ 
  if (lid==0)
  {
    xi[0] = -1.;
    xi[1] = -1.;
  }
  else if (lid==1)
  {
    xi[0] =  1.;
    xi[1] = -1.;
  }
  else if (lid==2)
  {
    xi[0] =  1.;
    xi[1] =  1.;
  }
  else if (lid==3)
  {
    xi[0] = -1.;
    xi[1] =  1.;
  }
  else
  {
	  std::stringstream oss;
    oss << "***ERR*** MOERTEL::Segment_BiLinearQuad::LocalCoordinatesOfNode:\n"
         << "***ERR*** local node number " << lid << " out of range (0..3)\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	throw ReportError(oss);
  }
  return true;
}

/*----------------------------------------------------------------------*
 | build an outward normal at a node adjacent to this        mwgee 10/05|
 | returns allocated vector of length 3 with outward normal             |
 *----------------------------------------------------------------------*/
double* MOERTEL::Segment_BiLinearQuad::BuildNormal(double* xi)
{ 
  // A bilinear quad in 3D can be warped, so it does matter where
  // to build the normal
  double G[3][3];
  Metric(xi,NULL,G);
  
  // the normal is G[2]
  double* n = new double[3];
  
  for (int i=0; i<3; ++i)
    n[i] = G[2][i];
  
  return n;
}

/*----------------------------------------------------------------------*
 | build basis vectors and metric at given point xi          mwgee 10/05|
 *----------------------------------------------------------------------*/
double MOERTEL::Segment_BiLinearQuad::Metric(double* xi, double g[], double G[][3])
{ 
  // get nodal coords;
  const double* x[4];
  for (int i=0; i<4; ++i) x[i] = nodeptr_[i]->X();
  
  // get shape functions and derivatives at xi
  double val[4];
  double deriv[8];
  EvaluateFunction(0,xi,val,4,deriv);
  
  // Build kovariant metric G1 and G2 = partial x / partial theta sup i
  for (int i=0; i<2; ++i)
    for (int dim=0; dim<3; ++dim)
    {
      G[i][dim] = 0.0;
      for (int node=0; node<4; ++node)
        G[i][dim] += deriv[node*2+i] * x[node][dim];
    }
  
  // build G3 as cross product of G1 x G2
  MOERTEL::cross(G[2],G[0],G[1]);
  
  // dA at this point is length of G[3] or |G1 x G2|
  double dA = MOERTEL::length(G[2],3);  
  return dA;
}

/*----------------------------------------------------------------------*
 | compute the length (Area) of this segment                 mwgee 10/05|
 *----------------------------------------------------------------------*/
double MOERTEL::Segment_BiLinearQuad::Area()
{ 
  double coord[4][2];
  double sqrtthreeinv = 1./(sqrt(3.));
  coord[0][0] = -sqrtthreeinv;
  coord[0][1] = -sqrtthreeinv;
  coord[1][0] =  sqrtthreeinv;
  coord[1][1] = -sqrtthreeinv;
  coord[2][0] =  sqrtthreeinv;
  coord[2][1] =  sqrtthreeinv;
  coord[3][0] = -sqrtthreeinv;
  coord[3][1] =  sqrtthreeinv;
  double A = 0.0;

  // create an integrator to get the gaussian points
  for (int gp=0; gp<4; ++gp)
  {
    double G[3][3];
    A += Metric(coord[gp],NULL,G);
  }
  return A;
}
