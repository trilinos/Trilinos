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
#include "mrtr_segment_linear1D.H"
#include "mrtr_interface.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 06/05|
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
MOERTEL::Segment_Linear1D::Segment_Linear1D(int id, int nnode, int* nodeId, int out) :
MOERTEL::Segment(id,nnode,nodeId,out)
{
  stype_ = MOERTEL::Segment::seg_Linear1D;
}

MOERTEL::Segment_Linear1D::Segment_Linear1D(int id, const std::vector<int>& nodev, int out) :
MOERTEL::Segment(id, nodev, out)
{
  stype_ = MOERTEL::Segment::seg_Linear1D;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 07/05|
 |  This constructor should not be used by the user, it is used         |
 |  internally together with Pack/Unpack for communication              |
 *----------------------------------------------------------------------*/
MOERTEL::Segment_Linear1D::Segment_Linear1D(int out) :
MOERTEL::Segment(out)
{
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 06/05|
 *----------------------------------------------------------------------*/
MOERTEL::Segment_Linear1D::Segment_Linear1D(MOERTEL::Segment_Linear1D& old) :
MOERTEL::Segment(old)
{
  // all date lives in the base class and is copied in MOERTEL::Segment(old)
}

/*----------------------------------------------------------------------*
 | pack all data in this segment into a vector               mwgee 07/05|
 *----------------------------------------------------------------------*/
int* MOERTEL::Segment_Linear1D::Pack(int* size)
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
    oss << "***ERR*** MOERTEL::Segment_Linear1D::Pack:\n"
         << "***ERR*** mismatch in packing size\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	throw ReportError(oss);
  }
  
  return pack;
}

/*----------------------------------------------------------------------*
 | unpack all data from a vector into this class             mwgee 07/05|
 *----------------------------------------------------------------------*/
bool MOERTEL::Segment_Linear1D::UnPack(int* pack)
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
    oss << "***ERR*** MOERTEL::Segment_Linear1D::UnPack:\n"
         << "***ERR*** mismatch in packing size\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	throw ReportError(oss);
  }
  
  return true;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 06/05|
 *----------------------------------------------------------------------*/
MOERTEL::Segment_Linear1D::~Segment_Linear1D()
{ 
  // data held in base class is destroyed by base class destructor
}

/*----------------------------------------------------------------------*
 |  clone this segment (public)                              mwgee 06/05|
 *----------------------------------------------------------------------*/
MOERTEL::Segment* MOERTEL::Segment_Linear1D::Clone()
{ 
  MOERTEL::Segment_Linear1D* newseg = new MOERTEL::Segment_Linear1D(*this);
  return (newseg);
}

/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 06/05|
 *----------------------------------------------------------------------*/
std::ostream& operator << (std::ostream& os, const MOERTEL::Segment_Linear1D& seg)
{
  seg.Print(); 
  return os;
}

/*----------------------------------------------------------------------*
 | build an outward normal at a node adjacent to this        mwgee 07/05|
 *----------------------------------------------------------------------*/
bool MOERTEL::Segment_Linear1D::LocalCoordinatesOfNode(int lid, double* xi)
{ 
  if (lid==0) xi[0] = -1.0;
  else if (lid==1) xi[0] = 1.0;
  else 
  {
	  std::stringstream oss;
    oss << "***ERR*** MOERTEL::Segment_Linear1D::LocalCoordinatesOfNode:\n"
  	 << "***ERR*** Segment " << Id() << ": node number " << lid << " out of range\n"
  	 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	throw ReportError(oss);
  }
  return true;
}

/*----------------------------------------------------------------------*
 | build basis vectors and metric at given point xi          mwgee 07/05|
 *----------------------------------------------------------------------*/
double MOERTEL::Segment_Linear1D::Metric(double* xi, double g[], double G[][3])
{ 
  // get nodal coordinates
  const double* x[2];
  x[0] = nodeptr_[0]->XCoords();
  x[1] = nodeptr_[1]->XCoords();
  
  // get shape functions
  double val[2];
  double deriv[2];
  functions_[0]->EvaluateFunction(*this,xi,val,2,deriv);

  double glocal[2];
  double* gl;
  if (g) gl = g;
  else   gl = glocal;
  
  // build covariant basis vector g = partial x / partial theta sup i
  for (int dim=0; dim<2; ++dim)
  {
    gl[dim] = 0.0;
    for (int node=0; node<2; ++node)
      gl[dim] += deriv[node] * x[node][dim];
  }
  
  // build metric tensor G sub ij = g sub i dot g sub j
  // in this 1D case, it's a scalar
  if (G)
  {
    G[0][0] = 0;
    for (int i=0; i<2; ++i) G[0][0] += gl[i]*gl[i];
  }

  // FIXME: look at file shell8/s8_tvmr.c & shell8/s8_static_keug.c 
  // build dA as g1 cross g2
  // in this linear 1D case, it's just the length of g
  double dl = sqrt(gl[0]*gl[0]+gl[1]*gl[1]);
  
  return dl;
}

/*----------------------------------------------------------------------*
 | build an outward normal at a node adjacent to this        mwgee 07/05|
 | returns allocated vector of length 3 with outward normal             |
 *----------------------------------------------------------------------*/
double* MOERTEL::Segment_Linear1D::BuildNormal(double* xi)
{ 
  // build the metric vectors at this local coordinates xi
  double g[3]; for (int i=0; i<3; ++i) g[i] = 0.0;
  Metric(xi,g,NULL);
  
  // in 3D, the outward normal is g1 cross g2, in 2D, the normal is
  // n1 = g2 and n2 = -g1
  double* n = new double[3];
  n[0] = g[1];
  n[1] = -g[0];
  n[2] = 0.0;
  double length = sqrt(n[0]*n[0]+n[1]*n[1]);
  n[0] /= length;
  n[1] /= length;
  return n;
}

/*----------------------------------------------------------------------*
 | compute the length (Area) of this segment                 mwgee 07/05|
 *----------------------------------------------------------------------*/
double MOERTEL::Segment_Linear1D::Area()
{ 
  // get nodal coordinates
  const double* x[2];
  x[0] = nodeptr_[0]->XCoords();
  x[1] = nodeptr_[1]->XCoords();

  // build vector from x[0] to x[1]
  double tangent[2];
  tangent[0] = x[1][0] - x[0][0];
  tangent[1] = x[1][1] - x[0][1];

  double length = sqrt(tangent[0]*tangent[0]+tangent[1]*tangent[1]);
  return length;
}
