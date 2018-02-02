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
MOERTEL::Segment_BiLinearTri::Segment_BiLinearTri(int id, int nnode, int* nodeId, int out) :
MOERTEL::Segment(id,nnode,nodeId,out)
{
  stype_ = MOERTEL::Segment::seg_BiLinearTri;
}

MOERTEL::Segment_BiLinearTri::Segment_BiLinearTri(int id, const std::vector<int>& nodev, int out) :
MOERTEL::Segment(id, nodev, out)
{
  stype_ = MOERTEL::Segment::seg_BiLinearTri;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 10/05|
 |  This constructor should not be used by the user, it is used         |
 |  internally together with Pack/Unpack for communication              |
 *----------------------------------------------------------------------*/
MOERTEL::Segment_BiLinearTri::Segment_BiLinearTri(int out) :
MOERTEL::Segment(out)
{
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 10/05|
 *----------------------------------------------------------------------*/
MOERTEL::Segment_BiLinearTri::Segment_BiLinearTri(MOERTEL::Segment_BiLinearTri& old) :
MOERTEL::Segment(old)
{
  // all date lives in the base class and is copied in MOERTEL::Segment(old)
}

/*----------------------------------------------------------------------*
 | pack all data in this segment into a vector               mwgee 10/05|
 *----------------------------------------------------------------------*/
int* MOERTEL::Segment_BiLinearTri::Pack(int* size)
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
    oss << "***ERR*** MOERTEL::Segment_BiLinearTri::Pack:\n"
         << "***ERR*** mismatch in packing size\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	throw ReportError(oss);
  }
  
  return pack;
}

/*----------------------------------------------------------------------*
 | unpack all data from a vector into this class             mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MOERTEL::Segment_BiLinearTri::UnPack(int* pack)
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
    oss << "***ERR*** MOERTEL::Segment_BiLinearTri::UnPack:\n"
         << "***ERR*** mismatch in packing size\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	throw ReportError(oss);
  }
  
  return true;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 10/05|
 *----------------------------------------------------------------------*/
MOERTEL::Segment_BiLinearTri::~Segment_BiLinearTri()
{ 
  // data held in base class is destroyed by base class destructor
}

/*----------------------------------------------------------------------*
 |  clone this segment (public)                              mwgee 10/05|
 *----------------------------------------------------------------------*/
MOERTEL::Segment* MOERTEL::Segment_BiLinearTri::Clone()
{ 
  MOERTEL::Segment_BiLinearTri* newseg = new MOERTEL::Segment_BiLinearTri(*this);
  return (newseg);
}

/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 10/05|
 *----------------------------------------------------------------------*/
std::ostream& operator << (std::ostream& os, const MOERTEL::Segment_BiLinearTri& seg)
{
  seg.Print(); 
  return os;
}

/*----------------------------------------------------------------------*
 | build an outward normal at a node adjacent to this        mwgee 10/05|
 *----------------------------------------------------------------------*/
bool MOERTEL::Segment_BiLinearTri::LocalCoordinatesOfNode(int lid, double* xi)
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
	  std::stringstream oss;
    oss << "***ERR*** MOERTEL::Segment_BiLinearTri::LocalCoordinatesOfNode:\n"
         << "***ERR*** local node number " << lid << " out of range (0..2)\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	throw ReportError(oss);
  }
  return true;
}

/*----------------------------------------------------------------------*
 | build basis vectors and metric at given point xi          mwgee 10/05|
 *----------------------------------------------------------------------*/
double MOERTEL::Segment_BiLinearTri::Metric(double* xi, double g[], double G[][3])
{ 
	  std::stringstream oss;
  oss << "***ERR*** MOERTEL::Segment_BiLinearTri::Metric:\n"
       << "***ERR*** not impl.\n"
       << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	throw ReportError(oss);
  return 0.0;
}

/*----------------------------------------------------------------------*
 | build an outward normal at a node adjacent to this        mwgee 10/05|
 | returns allocated vector of length 3 with outward normal             |
 *----------------------------------------------------------------------*/
double* MOERTEL::Segment_BiLinearTri::BuildNormal(double* xi)
{ 
  // linear triangles are planar, so we don't care were exactly to build the normal

  // build base vectors
  double g1[3];
  double g2[3];
  for (int i=0; i<3; ++i)
  {
    g1[i] = Nodes()[1]->XCoords()[i] - Nodes()[0]->XCoords()[i];  
    g2[i] = Nodes()[2]->XCoords()[i] - Nodes()[0]->XCoords()[i]; 
  }
  
  // build normal as their cross product
  double* n = new double[3];
  
  MOERTEL::cross(n,g1,g2);

  return n;
}

/*----------------------------------------------------------------------*
 | compute the length (Area) of this segment                 mwgee 10/05|
 *----------------------------------------------------------------------*/
double MOERTEL::Segment_BiLinearTri::Area()
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
