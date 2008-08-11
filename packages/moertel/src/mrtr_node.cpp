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
# Questions? Contact Michael Gee (mwgee@sandia.gov)
#
# ************************************************************************
#@HEADER
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
#include "mrtr_node.H"
#include "mrtr_interface.H"
#include "mrtr_pnode.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 06/05|
 *----------------------------------------------------------------------*/
MOERTEL::Node::Node(int Id, const double* x, int ndof, const int* dof, bool isonboundary, int out) :
Id_(Id),
outputlevel_(out),
iscorner_(false),
isonboundary_(isonboundary),
Drow_(null),
Mrow_(null),
Mmodrow_(null)
{
  seg_.resize(0);
  segptr_.resize(0);
  
  for (int i=0; i<3; ++i) 
  {
    x_[i] = x[i];
    n_[i] = 0.0;
  }
  
  dof_.resize(ndof);
  for (int i=0; i<ndof; ++i) dof_[i] = dof[i];
  
  LMdof_.resize(0);
  pnode_.resize(0);
  supportedby_.clear();
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 07/05|
 |  This constructor should not be used by the user, it is used         |
 |  used internally                                                     |
 *----------------------------------------------------------------------*/
MOERTEL::Node::Node(int out) :
Id_(-1),
outputlevel_(out),
iscorner_(false),
isonboundary_(false),
Drow_(null),
Mrow_(null),
Mmodrow_(null)
{
  seg_.resize(0);
  segptr_.resize(0);
  
  dof_.clear();
  LMdof_.resize(0);
  
  for (int i=0; i<3; ++i)
  { 
    n_[i] = 0.0;
    x_[i] = 0.0;
  }

  pnode_.resize(0);
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 06/05|
 *----------------------------------------------------------------------*/
MOERTEL::Node::Node(const MOERTEL::Node& old) :
supportedby_(old.supportedby_)
{
  Id_ = old.Id();
  outputlevel_ = old.outputlevel_;
  iscorner_ = old.iscorner_;
  isonboundary_ = old.isonboundary_;
  
  for (int i=0; i<3; ++i) 
  {
    x_[i] = old.x_[i];
    n_[i] = old.n_[i];
  }
  
  if (old.dof_.size())
    dof_ = old.dof_;
  else
    dof_.resize(0);
    
  if (old.LMdof_.size())
    LMdof_ = old.LMdof_;
  else
    LMdof_.resize(0);  
  
  if (old.seg_.size())
  {
    seg_.resize(old.seg_.size());
    seg_ = old.seg_;
  }
  else
    seg_.resize(0);
  
  if (old.segptr_.size())
  {
    segptr_ = old.segptr_;
  }
  else
    segptr_.resize(0);
  
  pnode_.resize(old.pnode_.size());  
  for (int i=0; i<(int)pnode_.size(); ++i)
    if (old.pnode_[i].get() != NULL)
    {
      pnode_[i] = rcp(new MOERTEL::ProjectedNode(*(old.pnode_[i])));
    }
  
  if (old.Drow_ != null)
  {
    map<int,double>* tmp = new map<int,double>(*(old.Drow_));
    Drow_ = rcp(tmp);
  }
  else
    Drow_ = null;    

  if (old.Mrow_ != null)
  {
    map<int,double>* tmp = new map<int,double>(*(old.Mrow_));
    Mrow_ = rcp(tmp);
  }
  else
    Mrow_ = null;    
}

/*----------------------------------------------------------------------*
 | pack data in this node into a vector                      mwgee 07/05|
 *----------------------------------------------------------------------*/
double* MOERTEL::Node::Pack(int* size)
{
  // *size = *size + Id_ + x_[3] + n_[3] + dof_.size() + ndof_*sizeof(double) + seg_.size() + nseg_*sizeof(double) + iscorner_ + isonboundary_
     *size = 1     + 1  +  3     + 3     +  1          + dof_.size()                + 1     + seg_.size()          + 1         + 1;
  double* pack = new double[*size];
  int count = 0;
  
  pack[count++] = (double)(*size);
  pack[count++] = (double)Id_;
  for (int i=0; i<3; ++i)
    pack[count++] = x_[i];
  for (int i=0; i<3; ++i)
    pack[count++] = n_[i];
  pack[count++] = (double)dof_.size();
  for (int i=0; i<(int)dof_.size(); ++i)
    pack[count++] = (double)dof_[i];
  pack[count++] = (double)seg_.size();
  for (int i=0; i<(int)seg_.size(); ++i)
    pack[count++] = (double)seg_[i];
  pack[count++] = (double)(iscorner_);
  pack[count++] = (double)(isonboundary_);

  if (count != *size)
  {
    cout << "***ERR*** MOERTEL::Node::Pack:\n"
         << "***ERR*** mismatch in packing size\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);     
  }
  
  return pack;
}

/*----------------------------------------------------------------------*
 | unpack data from a vector in this class                   mwgee 07/05|
 *----------------------------------------------------------------------*/
bool MOERTEL::Node::UnPack(double* pack)
{
  int count = 0;
  int size  = (int)pack[count++];
  Id_       = (int)pack[count++];
  for (int i=0; i<3; ++i)
    x_[i] = pack[count++];
  for (int i=0; i<3; ++i)
    n_[i] = pack[count++];
  dof_.resize((int)pack[count++]);  
  for (int i=0; i<(int)dof_.size(); ++i)
    dof_[i] = (int)pack[count++];
  seg_.resize((int)pack[count++]);
  for (int i=0; i<(int)seg_.size(); ++i)
    seg_[i] = (int)pack[count++];
  iscorner_ = (bool)pack[count++];
  isonboundary_ = (bool)pack[count++];
    
  if (count != size)
  {
    cout << "***ERR*** MOERTEL::Node::UnPack:\n"
         << "***ERR*** mismatch in packing size\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);     
  }
  
  return true;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 06/05|
 *----------------------------------------------------------------------*/
MOERTEL::Node::~Node()
{
  dof_.clear();
  LMdof_.clear();
  seg_.clear();
  segptr_.clear();
  pnode_.clear();
  Drow_ = null;
  Mrow_ = null;
  Mmodrow_ = null;
  supportedby_.clear();
}

/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 06/05|
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const MOERTEL::Node& node)
{ 
  node.Print();
  return (os);
}

/*----------------------------------------------------------------------*
 |  print node                                               mwgee 06/05|
 *----------------------------------------------------------------------*/
bool MOERTEL::Node::Print() const
{ 
  cout << "Node " << setw(6) << Id_ << "\tCoords ";
  for (int i=0; i<3; ++i)
    cout << setw(12) << x_[i] << " ";

  cout << "Normal " ;
  for (int i=0; i<3; ++i)
    cout << setw(12) << n_[i] << " ";

  cout << "#Dofs " << dof_.size() << " Dofs ";
  for (int i=0; i<(int)dof_.size(); ++i)
    cout << dof_[i] << " ";

  cout << "#LMDofs " << LMdof_.size();
  if (LMdof_.size())
  {
    cout << " LMDofs ";
    for (int i=0; i<(int)LMdof_.size(); ++i)
      cout << LMdof_[i] << " ";
  }
  
  if (IsCorner())
    cout << " is shared among 1D interfaces";

  if (IsOnBoundary())
  {
    cout << " is boundary of 2D-interface";
    cout << ", member of " << NSupportSet() << " support sets";
  }

  cout << endl;
  return true;
}

/*----------------------------------------------------------------------*
 |  set lagrange multiplier dof id                           mwgee 07/05|
 *----------------------------------------------------------------------*/
bool MOERTEL::Node::SetLagrangeMultiplierId(int LMId)
{ 
  // first check whether this dof has been set before
  // if so, do nothing
  for (int i=0; i<(int)LMdof_.size(); ++i)
    if (LMdof_[i] == LMId)
      return true;
      
  // resize the vector to tak the new dof
  LMdof_.resize(LMdof_.size()+1);
  
  // put in the new dof
  LMdof_[LMdof_.size()-1] = LMId;
  return true;
}

/*----------------------------------------------------------------------*
 |  add a segment id to my adjacency list                    mwgee 06/05|
 | (checking whether already present)                                   |
 | WRN: This is NOT a collective call!                                  |
 *----------------------------------------------------------------------*/
bool MOERTEL::Node::AddSegment(int sid)
{
  if (seg_.size())
  {
    // search whether sid already exists in seg_
    for (int i=0; i<(int)seg_.size(); ++i)
      if (sid==seg_[i]) 
        return true;
    
    // resize seg_
    seg_.resize(seg_.size()+1);
    
    // add new sid to seg_
    seg_[seg_.size()-1] = sid;
    return true;
  }
  else
  {
    seg_.resize(seg_.size()+1);
    seg_[seg_.size()-1] = sid;
  } 
  return true;
}

/*----------------------------------------------------------------------*
 |                                                           mwgee 07/05|
 | construct ptrs to redundant segments from my segment id list         |
 *----------------------------------------------------------------------*/
bool MOERTEL::Node::GetPtrstoSegments(MOERTEL::Interface& interface)
{ 
  if (!interface.IsComplete()) return false;
  if (!interface.lComm()) return true;
  if (!seg_.size()) return false;
  
  // vector segptr_ might already exist, build new
  segptr_.resize(seg_.size());
  
  for (int i=0; i<(int)seg_.size(); ++i)
  {
    segptr_[i] = interface.GetSegmentView(seg_[i]).get();
    if (!segptr_[i])
    {
      cout << "***ERR*** MOERTEL::Node::GetPtrstoSegments:\n"
           << "***ERR*** Interface " << interface.Id() << ": GetSegmentView failed\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);
    }
  }
  return true;
}

/*----------------------------------------------------------------------*
 |  build nodal normal                                       mwgee 07/05|
 *----------------------------------------------------------------------*/
bool MOERTEL::Node::BuildAveragedNormal()
{ 
  // get segments adjacent to me
  int nseg = Nseg();
  int* sid = SegmentIds();
  
  for (int i=0; i<3; ++i) n_[i] = 0.0;
  double weight = 0.0;

#if 0
  cout << "Building normal for node\n" << *this;
#endif
  
  for (int i=0; i<nseg; ++i)
  {
    MOERTEL::Segment* seg = segptr_[i]; 

#if 0
    cout << "Now averaging from Segment\n" << *seg;
#endif    

    if (!seg)
    {
      cout << "***ERR*** MOERTEL::Node::BuildAveragedNormal:\n"
           << "***ERR*** Node " << Id() << ": Segment " << sid[i] << " not found -> fatal\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);
    }
    double* n   = seg->BuildNormalAtNode(Id());
    double  wgt = seg->Area(); 
    
    // add weighted normal to n_
    for (int i=0; i<3; ++i) n_[i] += wgt*n[i];
    
    // add weight to total weight    
    weight += wgt;
    
    delete [] n; n = NULL;
  } // for (int i=0; i<nseg; ++i)

  double length = sqrt(n_[0]*n_[0]+n_[1]*n_[1]+n_[2]*n_[2]);
  for (int i=0; i<3; ++i) n_[i] /= length;

#if 0
  cout << "Node " << Id() << ":"
       << " normal is " << setw(15) << n_[0] 
       << "   "<< setw(15) << n_[1] << "   " << setw(15) << n_[2] << endl;
#endif

  return true;
}

/*----------------------------------------------------------------------*
 |  set a projected node                                     mwgee 07/05|
 *----------------------------------------------------------------------*/
bool MOERTEL::Node::SetProjectedNode(MOERTEL::ProjectedNode* pnode)
{ 
  pnode_.resize(pnode_.size()+1);
  pnode_[pnode_.size()-1] = rcp(pnode);
  return true;
}

/*----------------------------------------------------------------------*
 |  get projected nodes                                      mwgee 07/05|
 *----------------------------------------------------------------------*/
RefCountPtr<MOERTEL::ProjectedNode>* MOERTEL::Node::GetProjectedNode(int& length)
{ 
  length = pnode_.size();
  if (length)
    return &(pnode_[0]);
  else
    return NULL;
}

/*----------------------------------------------------------------------*
 |  get projected node                                       mwgee 07/05|
 *----------------------------------------------------------------------*/
RefCountPtr<MOERTEL::ProjectedNode> MOERTEL::Node::GetProjectedNode()
{ 
  int length = pnode_.size();
  if (length)
    return pnode_[0];
  else
    return null;
}

/*----------------------------------------------------------------------*
 |  add a value to the Drow_ map                             mwgee 11/05|
 *----------------------------------------------------------------------*/
void MOERTEL::Node::AddDValue(double val, int col)
{ 
  if (Drow_ == null)
    Drow_ = rcp(new map<int,double>());
    
  map<int,double>* Dmap = Drow_.get();
  
  (*Dmap)[col] += val;

  return;
}

/*----------------------------------------------------------------------*
 |  add a value to the Drow_ map                             mwgee 11/05|
 *----------------------------------------------------------------------*/
void MOERTEL::Node::AddMValue(double val, int col)
{ 
  if (Mrow_ == null)
    Mrow_ = rcp(new map<int,double>());
    
  map<int,double>* Mmap = Mrow_.get();
  
  (*Mmap)[col] += val;

  return;
}

/*----------------------------------------------------------------------*
 |  add a value to the Drow_ map                             mwgee 02/06|
 *----------------------------------------------------------------------*/
void MOERTEL::Node::AddMmodValue(int row, double val, int col)
{ 
  if (Mmodrow_ == null)
    Mmodrow_ = rcp(new vector<map<int,double> >(Ndof()));

  if ((int)Mmodrow_->size() <= row)
    Mmodrow_->resize(row+1);
    
  map<int,double>& Mmap = (*Mmodrow_)[row];
  
  Mmap[col] += val;

  return;
}
